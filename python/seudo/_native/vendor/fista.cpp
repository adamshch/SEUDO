// The FISTA method implementation.

#include <stdio.h>
#include <math.h>
#include "fista.hpp"

namespace Fista {

// ------------------- ScaledGradient ------------------------

ScaledGradient::~ScaledGradient()
{ }

void ScaledGradient::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{ }

// ------------------- Limiter -------------------------------

Limiter::~Limiter()
{ }

const std::shared_ptr<Limiter> NoLimiter;

// ------------------- RangeLimiter --------------------------

bool RangeLimiter::limit(int idx, double &v)
{
	// idx doesn't matter, the limit is the same for all of them
	if (v < low_) {
		v = low_;
		return true;
	}
	if (v > high_) {
		v = high_;
		return true;
	}
	return false;
}

// ------------------- Run -----------------------------------


Run::Run(const std::shared_ptr<ScaledGradient> gradient,
	const std::shared_ptr<Limiter> limiter,
	int xsize, double diffEps)
	: gradient_(gradient), limiter_(limiter),
	diff_(xsize), x_(xsize), diffEps2_(diffEps * diffEps),
	grad_(xsize)
{ }

Run::Run(const std::shared_ptr<ScaledGradient> gradient,
	const std::shared_ptr<Limiter> limiter,
	Vector &x, double diffEps)
	: gradient_(gradient), limiter_(limiter),
	diff_(x.size()), x_(x), diffEps2_(diffEps * diffEps),
	grad_(x.size())
{ }

void Run::reset(double value)
{
	step_ = 0;
	t_ = 1.;
	x_.assign(x_.size(), value);
	diff_.assign(x_.size(), 0.);
}

void Run::reset(const Vector &x)
{
	step_ = 0;
	t_ = 1.;
	x_ = x;
	diff_.assign(x_.size(), 0.);
}

bool Run::oneStep()
{
	if (restart_every_ > 0 && step_ % restart_every_ == 0) {
		// "Hit the wall" and restart from a stop.
		t_ = 1.;
	}

	// step 4.2
	double t_next = (1. + sqrt(1. + t_ * t_ * 4.)) / 2.;
	// momentum coefficient
	double nu = (t_ - 1.) / t_next;
	// The classic algorithm doesn't update t on step 0.  It doesn't seem to
	// matte rmuch one way or another, in theory updating it should be better
	// because it gives a faster acceleration of momentum, but in reality it's
	// a toss-up by a small amount. So just keep compatibility to allow easier
	// direct comparisons.
	if (step_ != 0) {
		t_ = t_next;
	}

	if (fastBrakeNu_ && fastBrake_ && nu != 0.) {
		nu = 1.;
	}

	// step 4.3 - the momentum step
	for (int i = 0; i < x_.size(); i++) {
		x_[i] += nu * diff_[i];
		// printf("    momentum step %d, nu=%f, x[%d]=%f df[%d]=%f\n", step_, nu, i, x_[i], i, diff_[i]);
	}

	// TODO: it should be relatively easy to detect the L being too small and
	// adjust it. Basically, it would manifest by the step 0 going in one
	// direction and then the step 1 going in the opposite direction in many
	// dimensions, by a relatively close amount. Then we can increase L,
	// reset the step count to 0 and t to 1, and try again. For non-quadratic
	// functions, this could happen not only on steps 0 and 1. And we can also
	// experiment with increasing L only in the directions that have overshot
	// like this.

	Vector lastgrad;
	if (fastBrake_) {
		// Preserve the last gradient for comparison.
		lastgrad.swap(grad_);
		grad_.resize(lastgrad.size());
	}

	// step 4.1
	gradient_->computeGradient(x_, grad_, step_);

	double maxdiff = 0.; // max squared direction
	double maxgrad = 0.;
	double diffnorm = 0.;
	double xnorm = 0.;

	for (int i = 0; i < x_.size(); i++) {
		double grad_i = grad_[i];
		double xi = x_[i] -= grad_i;
		xnorm += xi * xi;

		double df = nu * diff_[i] - grad_i;
		if (limiter_ && limiter_->limit(i, x_[i])) {
			// If got limited, stop the inertia
			df = 0.;
			grad_i = 0.;
		}

		if (fastBrake_ && lastgrad[i] * grad_i < 0) {
			// See the description of fastBrake_ for explanation.
			diff_[i] = 0.;
		} else {
			diff_[i] = df;
		}

		double diff2 = df * df;
		diffnorm += diff2;
		if (diff2 > maxdiff) {
			maxdiff = diff2;
		}
		double grad2 = grad_i * grad_i;
		if (grad2 > maxgrad) {
			maxgrad = grad2;
		}

		// printf("  step %d, x[%d]=%f df[%d]=%f, grad[%d]=%f\n", step_, i, x_[i], i, df, i, grad_i);
	}

	// printf("  step %d, diffnorm^2=%f, gradnorm^2=%f, eps^2=%f\n", step_, diffnorm, gradnorm, diffEps2_);
	++step_;

	bool stop;
	switch (stopping_) {
	case StopEpsEveryDimension:
		// Both the gradient and the total step have to be under the
		// stopping condition, or it could spuriously hit the case
		// where the inertia overshoots and the gradient returns it
		// back to almost the previous point, so just looking at the
		// total step is not enough.
		stop = (maxdiff <= diffEps2_ && maxgrad <= diffEps2_);
		break;
	case StopEpsNorm2:
		// The norm2 of the last step has to be under the stopping
		// condition. This considers being a long way off by one dimension
		// equivalent to being a little way off by many dimensions.
		stop = (diffnorm <= diffEps2_);
		// printf("diffnorm=%f maxdiff=%f maxgrad=%f\n", sqrt(diffnorm), sqrt(maxdiff), sqrt(maxgrad));
		break;
	case StopEpsNorm2Rel:
		// The norm2 of the last step has to be under the stopping
		// condition. This considers being a long way off by one dimension
		// equivalent to being a little way off by many dimensions.
		// And also has worse precision at large X.
		// This mode is the same as TFOCS mode 1.
		stop = (diffnorm <= diffEps2_ * fmax(1., xnorm));
		// printf("diffnorm=%f limt=%f maxdiff=%f maxgrad=%f\n", sqrt(diffnorm), sqrt(diffEps2_ * fmax(1., xnorm)), sqrt(maxdiff), sqrt(maxgrad));
		break;
	default:
		stop = false;
		break;
	}

	return stop;
}

int Run::repeatedSteps(int n)
{
	// The member step_ also contains the number of steps but it counts
	// them from the start, while this method can be called for the incremental
	// steps.
	for (int i = 0; i < n; i++) {
		bool stop = oneStep();
		if (stop) {
			// stopped after the step, so add 1
			return i + 1;
		}
	}
	return n + 1;
}

}; // namespace Fista

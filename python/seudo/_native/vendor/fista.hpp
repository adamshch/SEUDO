#ifndef __FISTA_HPP__
#define __FISTA_HPP__

#include <vector>
#include <memory>
#include <limits>

#include "fista_types.hpp"
#include "fista_v1.hpp"

// The FISTA method implementation, based on
// https://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf
//
// This is the second version of FISTA implementation, the first version can be
// found declared in fista_v1.hpp. The first version is kept mainly as a
// historic artifact: it matches the original paper closer, so the code there
// is easier to match up with the paper, and the development of the ideas that
// led to v2 can be traced through the v1 code. But it's a pain to use.  The
// second version benefits from the use of proper software engineering
// principles that make it easier to use, and also from a mathematical
// transformation that streamlines the computation. It also got a much more
// solid computation of the stopping condition that became easy to implement
// after the redesign, and easily composable boundaries for X.  These niceties
// add some overhead, somewhere around 1% for smallish realistic examples (and
// less for larger examples), but at the same time they enable a better
// handling of the boundaries that tends to reduce the number of needed steps
// by 10% or so, so it's a win overall. Also, need be, the niceties can be
// changed from OOP and virtual functions to templates and traits to get rid of
// the overhead and save that last 1%.
namespace Fista {

// A subclass of this class contains the implementation of gradient for a
// particular function for minimization, i.e. for function function p_L
// (p.189):
// 
// p_L(y) = argmin(x){ g(x) - L/2 || x - (y - (1/L)grad_f(y)) ||^2 }
//
// It knows how to compute (1/L)grad_f(y), or more exactly, (1/L)grad_fg(y),
// where "fg" is a name of function: fg(x) = f(x) + g(x).
// For details of why this makes sense, see the discussion in fista_gradient.cpp
// on PosMatScaledGradient (this is the "streamlined math" mentioned above).
class ScaledGradient
{
public:
	virtual ~ScaledGradient();

	// Compute gradient of (f(x) + g(x)), denoted in other comments as grad_fg(x),
	// downscaled by the value of L: (1/L)grad_fg(x).
	// @param x - input values
	// @param grad - produced gradient, scaled by 1/L
	// @param - step of the computation, the subclass may want to initialize L on step 0.
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step) = 0;

	// Call this from computeGradient() to trigger the computation of L if needed.
	// In case if the computation of L uses the current gradient as-is, call this
	// function after the raw (unscaled by L) gradient has been computed and pass it
	// as an argument, so that it won't have to be computed twice.
	inline void tryInitializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
	{
		if (step == 0 || reinit_L_every_ > 0 && step % reinit_L_every_ == 0) {
			initializeL(x, grad, step);
		}
	}

	// Gives the subclass a chance to compute L on the step 0, 
	// and then if reinit_L_every_ > 0 then every so many steps.
	// 
	// By default does nothing.
	//
	// @param x - the initial point of the step
	// @param grad - the computed gradient at this point (this is a minor optimization,
	//		since both the normal course of computation and the estimation of L need it,
	//		it gets computed up front and then used for both). This is the raw gradient,
	//      not scaled by L.
	// @param step - the current step, starting from 0, in case if the estimation
	//      depends on it
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

	// On step 0 and then if reinit_L_every_ > 0 then every so many steps
	// initializeL() will be called to get a chance to get the new estimate
	// of L. This can be useful for functions of a higher degree than quadratic.
	void setReinitEvery(int n)
	{
		reinit_L_every_ = n;
	}

public:
	// On step 0 and then if reinit_L_every_ > 0 then every so many steps
	// initializeL() will be called to get a chance to get the new estimate
	// of L. Since reinitialization is not typical, this value defaults to none
	// of it, and if it is needed, change this value after construction.
	int reinit_L_every_ = 0;
};

// This is an interface to the code that knows how to limit the values
// to an allowed range.
class Limiter
{
public:
	virtual ~Limiter();

	// Limit the vector of values to an allowed range.
	//
	// @param idx - index of value in X 
	// @param v - value in X to limit, in-place
	// @return true if the value got limited
	virtual bool limit(int idx, double &v) = 0;
};

// Limit all the values to the same range, low <= x_i <= high.
class RangeLimiter : public Limiter
{
public:
	RangeLimiter(double low, double high)
		: low_(low), high_(high)
	{ }

	// from Limiter
	virtual bool limit(int idx, double &v);

public:
	double low_;
	double high_;
};

// Limit the values to >= 0 ("positive" but also 0).
class PosLimiter : public RangeLimiter
{
public:
	PosLimiter()
		: RangeLimiter(0., std::numeric_limits<double>::infinity())
	{ }
};

// A convenient constant.
extern const std::shared_ptr<Limiter> NoLimiter;

// One run of the FISTA algorithm.
class Run
{
public:
	// @param gradient - encapsulation of computation of (1/L)grad_fg(x), and
	//   also of the target values
	// @param limiter - (may be nullptr) limiter on the range of values in X
	// @param xsize - size of the vector X (this size is often referred in the comments
	//   as K), X gets initialized to all 0s
	// @param diffEps - "difference epsilon", the computation gets stopped when
	//   the absolute value of each dimension is less than this, both in the
	//   complete step and in the gradient part of it.  Looking at just one of
	//   them tends to stop too early, so both are used.
	Run(const std::shared_ptr<ScaledGradient> gradient,
		const std::shared_ptr<Limiter> limiter,
		int xsize, double diffEps);

	// @param gradient - encapsulation of computation of (1/L)grad_fg(x), and
	//   also of the target values
	// @param limiter - (may be nullptr) limiter on the range of values in X
	// @param x - the initial vector x
	// @param diffEps - "difference epsilon", the computation gets stopped when
	//   the absolute value of each dimension is less than this, both in the
	//   complete step and in the gradient part of it.  Looking at just one of
	//   them tends to stop too early, so both are used.
	Run(const std::shared_ptr<ScaledGradient> gradient,
		const std::shared_ptr<Limiter> limiter,
		Vector &x, double diffEps);

	// Reset the internal state and initial X for a new computation.
	// Does not need to be done before the first computation.
	// @param value - all the elements of X get set to it
	void reset(double value = 0.);
	// @param x - the initial value of vector X
	void reset(const Vector &x);

	// Run one step of FISTA algorithm.
	// The steps on p. 193 become more convenient when rearranged, placing 4.1 after 4.3:
	//   (4.2) t_{k+1} = ...
	//   (4.3) y_{k+1} = ...
	//   (4.1) x_{k+1} = ...
	// So this function produces the new t_{k+1} and x_{k+1} from t_k, x_{k-1}, and x_k.
	// x_{k-1} is represented implicitly, as (x_k - diff).
	//
	// @return - true if the stopping condition of diffEps gets satisfied
	bool oneStep();

	// Make repeated steps, until either the maximum number is reached or the
	// stopping condition of diffEps gets satisfied (the steps become too small).
	// @param n - maximum number of steps to take.
	// @return - the number of step taken; if the stopping condition was not 
	//   satisfied, will return (n+1).
	int repeatedSteps(int n);

	// If > 0, resets the parameter t to 1 after each this many steps, "hitting
	// the wall" and restarting the inertia from a stop. This can be used if
	// the algorithm is "circling the drain" too much. Which happens generally
	// when the value of L gets estimated way too high. For the gradient
	// computations of quadratic functions declared in fista_gradient.hpp,
	// especially those with names starting with Multi-, this should not be an
	// issue.
	//
	// If ever used, the reasonable values of n are somewhere above 1000,
	// anything too small just messes up the momentum too much.
	void setRestartEvery(int n)
	{
		restart_every_ = n;
	}

public:
	// Class that computes the gradient.
	std::shared_ptr<ScaledGradient> gradient_;
	// Class that knows how to limit the values (or nullptr).
	std::shared_ptr<Limiter> limiter_;
	
	// The difference from the last step. Will be initialized to 0s.
	Vector diff_;
	// The current position in the computation.
	Vector x_;
	// Square of diffEps, more convenient for computation.
	double diffEps2_;
	// The coefficient that determines the scale of the inertia.
	double t_ = 1.;
	// The current step, starts from 0 and gets increased after every step.
	int step_ = 0;

	// These fields are not adjustable in the constructor but can be changed
	// from default values afterwards.

	// If > 0, resets the parameter t to 1 after each this many steps,
	// "hitting the wall" and restarting the inertia from a stop. This can be used if
	// the algorithm is "circling the drain" too much.
	int restart_every_ = 0;

	// Choose the mode of stopping on eps.
	enum Stopping {
		// When abs every dimension of both ISTA step and FISTA step is less then eps.
		StopEpsEveryDimension,
		// When norm2 of last FISTA step is less than eps.
		StopEpsNorm2,
		// When norm2 of last FISTA step is less than eps*max(norm2(x), 1).
		// This mode is the same as TFOCS mode 1.
		StopEpsNorm2Rel,
	};
	Stopping stopping_ = StopEpsEveryDimension;

	// The "fast brake" mode. Puts the brakes on faster after an overshoot in
	// each dimension of X (treating each dimension independently).
	// The idea here is that in normal FISTA when the momentum substep overshoots
	// the minimum in some dimension, the gradient substep puts us back into
	// the close vicinity of the optimum. But a large part of momentum is preserved,
	// and on the next step the momentum substep takes us away from the optimum again.
	// Then the gradient sustep puts us back again, and this time the sum of both
	// substeps is almost 0, so the momentum stops in 2 steps. The "fast brake"
	// mode instead sets the momentum to 0 whenever it sees the change in the
	// sign of gradient, thus stopping it in 1 step, twice faster. The result
	// tends to be about twice fewer steps.
	bool fastBrake_ = true;

	// Fast brake controls the FISTA Nu parameter. When it and the fast brake
	// mode are both enabled, Nu does not gradually diminish. This can be done
	// because the goal of Nu is also to add braking, and the fast braking can
	// replace it.
	bool fastBrakeNu_ = false;

	// Temporary values used in computation (and also can be inspected for
	// entertainment value).
	
	// Gradient at the current point.
	Vector grad_;
};
	
}; // namespace Fista

#endif // __FISTA_HPP__

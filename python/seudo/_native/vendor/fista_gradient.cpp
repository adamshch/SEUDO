// The FISTA method implementation.

#include <stdio.h>
#include <math.h>
#include "fista.hpp"
#include "fista_gradient.hpp"
#include "par_fista.hpp"

namespace Fista {

// ------------------- PosMatScaledGradient ------------------

PosMatScaledGradient::PosMatScaledGradient(
	const Fista::Vector &k,
	const Fista::Vector &b,
	const Fista::Vector &lambda)
	: k_(k), b_(b), lambda_(lambda)
{ }

void PosMatScaledGradient::computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step)
{
	// initializeL() here doesn't use the gradient, so it's OK to call before
	// the gradient is computed.
	tryInitializeL(x, grad, step);

	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] = inv_L_ * lambda_[gidx];
	}

	// The values in parentheses in each "column" of terms for dimensions
	// of the gradient are the same, so they can be computed once per
	// "column" and then reused for this term in all the dimensions:
	//   grad_f(x) = [
	//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
	//     df/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1),
	//   ]
	// Each "column" of terms corresponds to a row in the matrix k_,
	// so as we go through them sequentially, we can just keep increasing
	// the current position, and as we increase it after the last column of
	// the current row, we automatically get to the first column of the next row.
	size_t pos = 0;
	for (size_t i = 0; i < b_.size(); i++) {
		double v = 0.;
		for (size_t j = 0; j < x.size(); j++) {
			v +=  k_[pos] * x[j];
			++pos;
		}
		v -= b_[i];

		// Apply to all dimensions, with 1/L.
		for (size_t gidx = 0; gidx < grad.size(); gidx++) {
			grad[gidx] += inv_L_ * 2. * v * k_[i*x.size() + gidx];
		}
	}
}

// Same logic as in v1's PosMatSquareMinimizer::computeLMaxGradient().
//
// This computation produces a very good result (further improved 
// by the "multi-L" approach as can be found in
// MultiPosMatScaledGradient::initializeL()). It starts with the question:
// why do we bother at all with L? The answer is that when we make a
// single gradient step (in basic ISTA) we won't overshoot the minimum.
// If we do overshoot it, and especially if we overshoot it by much,
// we can get into a situation when we're getting farther and farther
// from the minimum.
//
// Since the step by each dimension x_i is proportional to partial
// derivative df/dx_i in the gradient, we need to pick such a proportion
// that this step won't overshoot the minimum, and (1/L) is this proportion.
// Since we compute the step by each dimension separately, instead of
// using the norm2 of the gradient in computing L, we can consider
// each dimension separately, and then pick the largest L out of all
// dimensions. This saves us all the trouble with the squares and sums.
//
// We can approximate for any function:
//
// [df/dx_m](y) = [df/dx_m](x) + sum{i}([d^2f/dx_m dx_i](x) * (y_i-x_i))
//
// And since for the quadratic functions the second derivative is 
// constant, this would even be the exact expression. (If we ever care
// about the higher powers, we can use the fact that we're always trying
// to move towards the minimum, so we need to consider the upper boundary
// only in the region more or less between the initial point and the
// minimum, building in some safety margin for the possible detours).
//
// Then
//   L_m = ( [df/dx_m](y) - [df/dx_m](x) ) / (y_m - x_m)
//       = sum{i}([d^2f/dx_m dx_i](x) * (y_i - x_i)) / (y_m - x_m)
//
// Next, the (y_i - x_i) on top and (y_m - x_m) are begging to be eliminated
// but they're obviously not the same, so we need to explain it somehow.
// A simple handwavy way is to say that we'll choose to move by 1 in every
// dimension and they all will be the same. But in reality we don't move
// by 1 in all dimensions. Though we don't move by by arbitrary amounts
// either, we move by the amount proportional to the gradient.
// So we can say that
//   (y_i - x_i) / (y_m - x_m) = [df/dx_i](x) / [df/dx_m](x)
// and since we're looking for the maximum L_m, we should also be
// looking for the maximum of this ratio.
//
//     dfg/dx_i = 2 * sum{0<=j<N}( k_ji * (sum{0<=p<K}(k_jp * x_p) - b_j) )
//     dfg/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=p<K}(k_jp * x_p) - b_j) )
// 
// And there is gets complicated. But here is a handwavy argument:
// If there is no cross-dependency between x_i and x_m then
//   [d^2f/dx_m dx_i] = 0
// and the ratio (y_i - x_i) / (y_m - x_m) won't matter. If there is a strong
// cross-dependency between them, then since second derivatives are symmetric
//   [d^2f/dx_m dx_i] = [d^2f/dx_i dx_m]
// their gradients would generally be close, so we can reasonably assume
//   (y_i - // x_i) / (y_m - x_m) = 1
// It's a leap of faith, I have no strict proof of it, but it seems to work in
// practice.
//
// Then assuming that this works, we get (y_m - x_m) eliminated from
// both sides of the fraction, and
//   L_m = sum{0<=i<K}( [d^2f/dx_m dx_i](x) )
//
// All we need to do is find the second derivatives. Given
//   df/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) )
//
// we can produce
//   [d^2f/dx_m dx_i] = 2 * sum{0<=j<N}( k_jm * k_ji )
//
// And so
//   L_m = sum{0<=i<K}( 2 * sum{0<=j<N}( k_jm * k_ji ) )
//       = 2 * sum{0<=i<K}( sum{0<=j<N}( k_jm * k_ji ) )
//       = 2 * sum{0<=i<K}( sum{0<=j<N}( k_jm * k_ji ) )
//       = 2 * sum{0<=j<N}( k_jm * sum{0<=i<K}(k_ji) )
//
// And this nicely meshes with the logic of the pessimistic estimation:
//
// Looking at expression 
//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
//
// or in generalized form
//     df/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) )
//
// we can first eliminate b_j because they are constant and get:
//     df/dx_m = 2 * sum{0<=j<N}( k_jm * sum{0<=i<K}(k_ji * x_i)  ) 
//             = 2 * sum{0<=j<N}( sum{0<=i<K}(k_jm * k_ji * x_i)  ) 
//
// Here the maximum of multiplying any two values in k equals the maximum of square of
// any single value: max(k_ij*k_ab) = max(k_ij^2).
// And there are N of these values in each parenthesised term, and K terms.
// So for every dimension of the gradient computed on x that has every dimension x_i <= 1,
// grad_f(x)_i <= 2*K*N*max(k_ij^2). For K dimensions, the norm2 of gradient difference will
// be at most 2*K*N*max(k_ij^2)*sqrt(K).
// If the length of each dimension in x is 1, the norm2 of x, ||x|| = sqrt(K). So
//   L = 2*K*N*max(k_ij^2)*sqrt(K) / sqrt(K) = 2*K*N*max(k_ij^2)
// which would be the high estimate.
//
// When all k_ij are the same, as assumed in the worst case of the pessimistic
// estimation, both formulas will produce the same result. But if the dependency
// matrix k is sparse, the maximum L_m will be much lower than the pessimistic
// estimation.
//
// The L computed here works much better than the pessimistic
// estimation. It creates a much faster acceleration, and even more
// importantly, deceleration on overshooting, and reduces the
// required number of steps by an order of magnitude.
void PosMatScaledGradient::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	const size_t xsz = lambda_.size();

	// This is similar to "columns" in the gradient computation by drawing.
	Fista::Vector columns(b_.size(), 0.);

	int pos = 0;
	for (int i = 0; i < b_.size(); i++) {
		double v = 0.;
		for (int j = 0; j < lambda_.size(); j++) {
			v += k_[pos++];
		}
		columns[i] = v;
	}

	// these are L for each dimension of the gradient
	double l = 0.;
	for (int j = 0; j < xsz; j++) {
		double v = 0.;
		for (int i = 0; i < b_.size(); i++) {
			v += k_[i * xsz + j] * columns[i];
		}
		v *= 2.;
		// printf("  Lpart[%d] = %f\n", j, v);
		if (v > l)
			l = v;
	}
	// printf("better L = %f\n", l);

	if (l < 0.0000001)
		l = 0.0000001;
	inv_L_ = 1./l;
}

// ------------------- BaseDrawScaledGradient ----------------

BaseDrawScaledGradient::BaseDrawScaledGradient(
	int wd,
	int ht,
	std::shared_ptr<Drawing> drawing,
	const Fista::Vector &b,
	const Fista::Vector &lambda,
	int nthreads)
	: Drawable(wd, ht), drawing_(drawing), neg_b_(b.size()), lambda_(lambda),
	drawColumns_(b.size()), drawGrad_(nullptr)
{
	initialize(b, nthreads, drawing->isParallel());
}

BaseDrawScaledGradient::BaseDrawScaledGradient(
	int wd,
	int ht,
	std::shared_ptr<BaseDrawScaledGradient::DrawingToMe> drawingPass1,
	std::shared_ptr<BaseDrawScaledGradient::DrawingToMe> drawingPass2,
	const Fista::Vector &b,
	const Fista::Vector &lambda,
	int nthreads)
	: Drawable(wd, ht), drawingPass1_(drawingPass1), drawingPass2_(drawingPass2),
	neg_b_(b.size()), lambda_(lambda),
	drawColumns_(b.size()), drawGrad_(nullptr)
{
	initialize(b, nthreads, drawingPass1->isParallel());
}

void BaseDrawScaledGradient::initialize(
	const Fista::Vector &b,
	int nthreads,
	bool parallel)
{
	for (size_t i = 0; i < b.size(); i++) {
		neg_b_[i] = -b[i];
	}

	if (nthreads < 1) {
		nthreads = 1;
	}
	if (!parallel) {
		nthreads = 1;
	}

	// for pass 1
	int nthreads1 = nthreads > ht_? ht_ : nthreads;
	// for pass 2
	int nthreads2 = nthreads > lambda_.size()? lambda_.size() : nthreads;
	nthreads = std::max(nthreads1, nthreads2);

	// Fill the limits on pass 1, partitioned by destination, by whole rows of pixels.
	{
		limits_pass_1_.resize(nthreads1);
		int last = 0; // end of last interval
		int left = ht_; // how many rows are left to split
		// This bunches up the slightly shorter (by 1 unit) intervals up front,
		// but this bunching doesn't matter here.
		for (int i = 0; i < nthreads1; i++) {
			limits_pass_1_[i].start_input_ = 0;
			limits_pass_1_[i].end_input_ = lambda_.size();
			limits_pass_1_[i].start_dest_y_ = last;

			int step = left / (nthreads1 - i);
			last += step;
			left -= step;

			limits_pass_1_[i].end_dest_y_ = last;
			// printf("   pass1 %d of %d: [%d, %d)\n", i, nthreads1, limits_pass_1_[i].start_dest_y_, limits_pass_1_[i].end_dest_y_);
		}
	}

	// Fill the limits on pass 2, partitioned by input.
	{
		limits_pass_2_.resize(nthreads2);
		int last = 0; // end of last interval
		int left = lambda_.size(); // how many inputs are left to split
		// This bunches up the slightly shorter (by 1 unit) intervals up front,
		// but this bunching doesn't matter here.
		for (int i = 0; i < nthreads2; i++) {
			limits_pass_2_[i].start_dest_y_ = 0;
			limits_pass_2_[i].end_dest_y_ = ht_;
			limits_pass_2_[i].start_input_ = last;

			int step = left / (nthreads2 - i);
			last += step;
			left -= step;

			limits_pass_2_[i].end_input_ = last;
			// printf("   pass2 %d of %d: [%d, %d)\n", i, nthreads2, limits_pass_2_[i].start_input_, limits_pass_2_[i].end_input_);
		}
	}

	// Start the threads (one less than requested).
	for (int i = 1; i < nthreads; i++) {
		threads_.emplace_back(std::make_shared<GradientThread>(this));
	}
}

void BaseDrawScaledGradient::computeUnscaledGradient(bool zeroed, const Fista::Vector &x, Fista::Vector &grad)
{
	if (zeroed) {
		grad.assign(x.size(), 0.);
		drawColumns_.assign(neg_b_.size(), 0.);
	} else {
		// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
		grad = lambda_;
		
		// Initialize the columns to -b.
		drawColumns_ = neg_b_;
	}

	drawGrad_ = &grad;

	// Do the pass 1, collect the column data.
	drawPass1_ = true;
	for (int i = 1; i < limits_pass_1_.size(); i++) {
		threads_[i - 1]->compute(&x, &limits_pass_1_[i]);
	}
	// Do the fist partition in this thread.
	selfDraw(x, limits_pass_1_[0]);
	for (int i = 1; i < limits_pass_1_.size(); i++) {
		threads_[i - 1]->wait();
	}
	
	// Do the pass 2, add up the data by gradient dimensions.
	drawPass1_ = false;
	for (int i = 1; i < limits_pass_2_.size(); i++) {
		threads_[i - 1]->compute(&x, &limits_pass_2_[i]);
	}
	// Do the fist partition in this thread.
	selfDraw(x, limits_pass_2_[0]);
	for (int i = 1; i < limits_pass_2_.size(); i++) {
		threads_[i - 1]->wait();
	}

	drawGrad_ = nullptr;
}

void BaseDrawScaledGradient::drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
{
	// The way the drawing is done, to avoid collecting a huge matrix, it's repeated
	// twice, with the intermediate data collected in pass 1 and used in pass 2. The
	// easiest way to understand it is to look how the gradient is computed with
	// a matrix. This is the same computation but rearranged to write each pair
	// or (from_idx, to_idx) at most once per pass.
	int to_idx = to_y * wd_ + to_x;
	if (drawPass1_) {
		// The "columns" are per destination index.
		drawColumns_[to_idx] += k * x[from_idx];
	} else {
		// Now we can collect data from one column in one gradient dimension.
		(*drawGrad_)[from_idx] += 2. * k * drawColumns_[to_idx];
	}
}

void BaseDrawScaledGradient::selfDraw(const Fista::Vector &input, Drawing::Limits &limits)
{
	if (drawing_) {
		drawing_->draw(input, *this, limits);
	} else {
		if (drawPass1_) {
			drawingPass1_->draw(input, *this, limits);
		} else {
			drawingPass2_->draw(input, *this, limits);
		}
	}
}

// ------------------- PosDrawScaledGradient -----------------

void PosDrawScaledGradient::computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step)
{
	computeUnscaledGradient(/*zeroed*/ false, x, grad);

	// initializeL() here gets the unscaled gradient.
	tryInitializeL(x, grad, step);

	// Scale by 1/L.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] *= inv_L_;
	}
}

void PosDrawScaledGradient::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	if (step == 0) {
		// A copy of PosMatScaledGradient::computeLMaxGradient()
		// done via drawing.
		// The computation is the same as computing the gradient but with b=0,
		// lambda=0, and x=[all 1s]
		Fista::Vector lpart;
		Fista::Vector input_1(x.size(), 1.);
		computeUnscaledGradient(/*zeroed*/ true, input_1, lpart);

		double l = 0.;
		for (int j = 0; j < x.size(); j++) {
			double v = lpart[j];
			// printf("  Lpart[%d] = %f\n", j, v);
			if (v > l)
				l = v;
		}
		// printf("L = %f\n", l);
		inv_L_ = 1./l;
	} else {
		// Get here if asked to recompute on the following steps.
		// Like TFOCS, measure the actual ratio
		//   L = ||grad_f(x) - grad_f(y)|| / ||x - y||
		// using the current and the previous point.
		double norm_x = 0., norm_grad = 0.;
		for (int j = 0; j < x.size(); j++) {
			double v = x[j] - last_x_[j];
			norm_x += v*v;
		}
		for (int j = 0; j < grad.size(); j++) {
			double v = grad[j] - last_grad_[j];
			norm_grad += v*v;
		}
		if (norm_grad != 0. && norm_x != 0.) {
			double inv_l = sqrt(norm_x) / sqrt(norm_grad);
			
			// To leave a margin of safety, don't allow to reduce L
			// by a factor of more than max_drop per step;
			// This is 1/beta from TFOCS.
			static const double max_drop = 2.;

			if (inv_l < inv_L_) {
				// printf("   L grew %f -> %f\n", 1./inv_L_, 1./inv_l);
				inv_L_ = inv_l;
			} else if (inv_l < inv_L_ * max_drop) {
				// printf("   L dropped %f -> %f\n", 1./inv_L_, 1./inv_l);
				inv_L_ = inv_l;
			} else {
				// printf("   L dropped %f -> %f, limiting to %f\n", 1./inv_L_, 1./inv_l, 1./(max_drop*inv_L_));
				inv_L_ *= max_drop;
			}
		}
	}

	// remember for the next call
	last_x_ = x;
	last_grad_ = grad;
}

// ------------------- MultiPosMatScaledGradient -------------

MultiPosMatScaledGradient::MultiPosMatScaledGradient(
	const Fista::Vector &k,
	const Fista::Vector &b,
	const Fista::Vector &lambda)
	: k_(k), b_(b), lambda_(lambda), inv_L_(lambda.size())
{ }

void MultiPosMatScaledGradient::computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step)
{
	// initializeL() here doesn't use the gradient, so it's OK to call before
	// the gradient is computed.
	tryInitializeL(x, grad, step);

	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	grad = lambda_;

	// The values in parentheses in each "column" of terms for dimensions
	// of the gradient are the same, so they can be computed once per
	// "column" and then reused for this term in all the dimensions:
	//   grad_f(x) = [
	//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
	//     df/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1),
	//   ]
	// Each "column" of terms corresponds to a row in the matrix k_,
	// so as we go through them sequentially, we can just keep increasing
	// the current position, and as we increase it after the last column of
	// the current row, we automatically get to the first column of the next row.
	size_t pos = 0;
	for (size_t i = 0; i < b_.size(); i++) {
		double v = 0.;
		for (size_t j = 0; j < x.size(); j++) {
			v +=  k_[pos] * x[j];
			++pos;
		}
		v -= b_[i];

		// Apply to all dimensions.
		for (size_t gidx = 0; gidx < grad.size(); gidx++) {
			grad[gidx] += 2. * v * k_[i*x.size() + gidx];
		}
	}

	// Scale by 1/L.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] *= inv_L_[gidx];
	}
}

// Same logic as in v1's MultiPosMatSquareMinimizer::initializeL().
//
// This is a further development of the idea from
// PosMatScaledGradient::initializeL(): if we're considering each
// dimension x_i separately, why do we care to have a common L?
// We're only interested in the movement in each dimension not overshooting
// the minimum, so we can pick a separate L for each dimension.
// Yeah, it won't exactly be a "gradient descent", since the gradient
// will get skewed by these different Ls, but who cares if it works.
//
// The tests in t_mmsq_fista.cpp explore the examples of this.
// Such as, if we have the function
//   f(x) = 1 * x_0 + 4 * x_1
// then using the common maximum L gets one dimension right on point
// on the first step, while the first step gets there slowly. But with
// the separate L_0 and L_1, it gets both dimensions on point on the
// first step. Things go less great when all 4 elements of a 2*2 matrix
// get populated, since the dimensions get mixed. Theoretically speaking,
// it's possible to untangle the gradients in the general case and go
// to the minimum in one step by solving a system of linear equations. 
// But for all I can tell, solving this system will require more work
// than doing the descent. Also, there would be a complication with
// hitting a boundary of the domain for some dimension, that the gradient
// descent solves much better.
//
// So, depending on how varied the various dimensions of L are, this
// approach can give a massive win with the basic ISTA. When the
// "Fast" part gets added to make FISTA, the effect becomes smaller,
// since FISTA is good at accelerating the descent with the momentum.
// Still, using separate dimensions of L usually gives an improvement
// of 10-30% over using the common maximum L.
void MultiPosMatScaledGradient::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	const size_t xsz = lambda_.size();

	// This is similar to "columns" in the gradient computation by drawing.
	Fista::Vector columns(b_.size(), 0.);

	int pos = 0;
	for (int i = 0; i < b_.size(); i++) {
		double v = 0.;
		for (int j = 0; j < lambda_.size(); j++) {
			v += k_[pos++];
		}
		columns[i] = v;
	}

	// these are L for each dimension of the gradient
	for (int j = 0; j < xsz; j++) {
		double v = 0.;
		for (int i = 0; i < b_.size(); i++) {
			v += k_[i * xsz + j] * columns[i];
		}
		v *= 2.;
		// printf("  Lpart[%d] = %f\n", j, v);
		inv_L_[j] = 1. / v;
	}
}

// ------------------- MultiPosDrawScaledGradient ------------

void MultiPosDrawScaledGradient::computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step)
{
	computeUnscaledGradient(/*zeroed*/ false, x, grad);

	// initializeL() gets the unscaled gradient (and doesn't use it anyway)
	tryInitializeL(x, grad, step);

	// Scale by 1/L.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] *= inv_L_[gidx];
	}
}

// Same logic as in MultiPosMatScaledGradient::initializeL().
void MultiPosDrawScaledGradient::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	// The computation is the same as computing the gradient but with b=0,
	// lambda=0, and x=[all 1s]
	Fista::Vector input_1(x.size(), 1.);
	computeUnscaledGradient(/*zeroed*/ true, input_1, inv_L_);

	double l = 0.;
	for (int j = 0; j < x.size(); j++) {
		double v = inv_L_[j];
		// printf("  Lpart[%d] = %f\n", j, v);
		if (v == 0.) {
			inv_L_[j] = 1.;
		} else {
			inv_L_[j] = 1. / v;
		}
	}
}

}; // namespace Fista

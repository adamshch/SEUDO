// The FISTA method implementation.

#include <stdio.h>
#include <math.h>
#include "fista.hpp"
#include "fista_v1_minimizer.hpp"

namespace Fista {

// ------------------- GradientMinimizer ---------------------

bool GradientMinimizer::compute(Vector &y, Vector &x, int step)
{
	// Temporarily put the gradient into x.
	computeGradient(y, x);

	if (step == 0 || reinit_L_every_ > 0 && step % reinit_L_every_ == 0) {
		// Auto-compute L and 1/L.
		initializeL(y, x, step);
	}

	// This is really norm2(grad_fg(x))^2.
	double gradnorm = 0.;
	// Do the step by dimensions.
	for (size_t i = 0; i < x.size(); i++) {
		double gd = x[i]; // gradient dimension
		double xval = y[i] - inv_L_ * gd;
		if (xval <= 0.) {
			xval = 0.;
			if (positive_ && gd > 0.)  {
				// If we've hit the boundary of x>=0 and the gradient points
				// further downwards, we can't go there, so consider this
				// equivalent to gradient dimension becoming 0.
				gd = 0.;
			}
		}
		if (upperBounded_ && xval > upperBound_) {
			xval = upperBound_;
			gd = 0.;
		}
		x[i] = xval;
		gradnorm += gd*gd;
		verbose_ && printf("   p_L[%zd] (%f -> %f); gd=%f\n", i, y[i], xval, gd);
	}
	verbose_ && printf("   norm=%f\n", sqrt(gradnorm));

	// True if length (norm) of the gradient is less than epsilon.
	return gradnorm <= eps2_;
}

void GradientMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{ }

// ------------------- GradientSquareMinimizer ---------------

void GradientSquareMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
#if 0
	// Get one more point, in the direction against the gradient.
	int sz = x.size();
	double norm2 = 0.; // ||x - y||
	Fista::Vector y(sz);
	for (int i = 0; i < sz; i++) {
		double df = grad[i];
		y[i] = x[i] - df;
		norm2 += df * df;
	}
	norm2 = sqrt(norm2);

	Fista::Vector ygrad(sz);
	computeGradient(y, ygrad);

	double norm2grad = 0.; // ||grad_f(x) - grad_f(y)||
	for (int i = 0; i < sz; i++) {
		double df = grad[i] - ygrad[i];
		norm2grad += df * df;
	}
	norm2grad = sqrt(norm2grad);

	double l = (norm2 == 0? 1. : norm2grad / norm2);
	// printf("L = %f  norm2 = %f, norm2grad = %f\n", l, norm2, norm2grad);
	if (l < 0.00001)
		l = 0.00001;
	inv_L_ = 1. / l;
#endif
}

// ------------------- PosMatSquareMinimizer -----------------

PosMatSquareMinimizer::PosMatSquareMinimizer(
	const Fista::Vector &k,
	const Fista::Vector &b,
	const Fista::Vector &lambda,
	double eps)
	: GradientSquareMinimizer(eps, /*positive*/ true, /*reinit_L_every*/ 0),
	k_(k), b_(b), lambda_(lambda)
{
	// See the explanation of estimation for L above.
	double l;
	// l = computeLOptimistic();
	// printf("optimistic L = %f\n", l);
	// l = computeLHessian();
	// l = computeLPessimistic();
	// printf("pessimistic L = %f\n", l);
	l = computeLMaxGradient();
	// printf("better L = %f\n", l);

	// Adjust in parent after L becomes known.
	inv_L_ = 1. / l;
}

// This computation is WRONG and occasionally underestimates L with terrible
// consequences. It's kept for historic reasons, as an example of what not
// to do.
//
// The logic for it is:
// We can start be looking at L^2, which will be maximized in the same point as
// L.  Since grad_f(x) is linear, the max L will be the same for any point y
// and we can pick y=0 without loss of generality, and the gradient will scale
// linearly with x, so we can pick x = [1, 1, 1, ...] because I think all the
// dimensions of x being equal will be the worst case for the ratio of squares
// (at least with k >= 0, that we can assume for SEUDO). So essentially we can
// compute the gradient at points 0 and [1, 1, 1, ...], subtract one of them
// from another, and square the dimensions. Add them up, divide by the number
// of dimensions (since that would be the sum of all the squared values in x),
// that will give the estimation of L^2, and then take the square root to get
// L. And since grad_g(x) is constant for x >= 0, we can even reuse the
// gradient computation function that computes grad_fg(x).
double PosMatSquareMinimizer::computeLOptimistic()
{
	// See the explanation of estimation for L above.
	Fista::Vector grad_at_0(lambda_.size());
	Fista::Vector grad_at_1(lambda_.size());
	const Fista::Vector input_0(lambda_.size(), 0.);
	const Fista::Vector input_1(lambda_.size(), 1.);

	computeGradient(input_0, grad_at_0);
	computeGradient(input_1, grad_at_1);

	double l = 0.;
	for (size_t i = 0; i < grad_at_0.size(); i++) {
		double v = grad_at_1[i] - grad_at_0[i];
		l += v * v;
	}
	return sqrt(l / grad_at_0.size());
}

// This is a safe computation but tends to oversetimate L, leading to slow
// convergence.
// Looking at expression 
//
//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
//
// or in generalized form
//
//     df/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) )
//
// we can first eliminate b_j because they are constant and get:
//
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
double PosMatSquareMinimizer::computeLPessimistic()
{
	// See the explanation of estimation for L above.
	double l = 0;
	for (size_t i = 0; i < k_.size(); i++) {
		double v = k_[i] * k_[i];
		if (v > l)
			l = v;
	}
	return l * 2. * k_.size();
}

// This computation is WRONG and occasionally underestimates L with terrible
// consequences. It's kept for historic reasons, as an example of what not
// to do.
//
// Its premise is that
// suppose we move the point x a little, how would the expression
// ( ||grad_f(x) - grad_f(y)|| / ||x - y|| ) change? The shorter was the
// previous vector (x - y), the greater effect will be cause by this change.
// So to see the greatest effect, we should assume x = y. And the vector
// (grad_f(x) - grad_f(y)) can be computed from the second gradient (Hessian)
// of f. If we move x in exactly one dimension, (grad_f(x) - grad_f(y)) will
// be equal to one "slice" of the Hessian. And moving x in exacly one dimension
// would also give the greatest possible change of the gradient. Since the
// second gradient of this f() is a constant, moving x by any distance will
// give a proportional effect, so we can chose to move it by 1.  So if we try
// to move x in every dimension from 0 to 1 and pick the one that produces the
// change of gradient with the highest norm2, this will be the upper bound
// on the change of gradient's norm2, and will give L.
double PosMatSquareMinimizer::computeLHessian()
{
	// See the explanation of estimation for L above.
	double lsq = 0.;

	const size_t xsz = lambda_.size();
	// One "slice" of 2nd gradient (Hessian), by a single dimension of x
	Fista::Vector slice(xsz);

	for (size_t k = 0; k < xsz; k++) { // Slice for each x_k
		slice.assign(xsz, 0.);

		// The structure of this loop is the same as for the gradient computation,
		// so see the explanation there.
		for (size_t i = 0; i < b_.size(); i++) {
			// only x_k is equal to 1 and not 0, so the loop that computes
			// the "column term" for the gradient can be skipped
			double v = k_[i*xsz + k];

			// Apply to all dimensions of the slice. Skip the coefficient 2 for now.
			for (size_t j = 0; j < xsz; j++) {
				slice[j] += v * k_[i*xsz + j];
			}
		}

		// Now compute norm2 squared of the slice.
		double norm2sq = 0.;
		for (size_t j = 0; j < xsz; j++) {
			// The coefficient 2 got moved here.
			norm2sq += (2. * 2.) * (slice[j] * slice[j]);
		}
		
		if (norm2sq > lsq)
			lsq = norm2sq;
	}
	return sqrt(lsq);
}

// This computation produces a very good result (further improved in
// MultiPosMatSquareMinimizer::initializeL()). It starts with the question:
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
// And this nicely meshes with the logic in computeLPessimistic(),
// where all k_ij are the same, both formulas will produce the same
// result.
//
// The L computed here works much better than the pessimistic
// estimation. It creates a much faster acceleration, and even more
// importantly, deceleration on overshooting, and reduces the
// required number of steps by an order of magnitude.
double PosMatSquareMinimizer::computeLMaxGradient()
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

	return l;
}

void PosMatSquareMinimizer::computeGradient(const Fista::Vector &x, Fista::Vector &grad)
{
	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] = lambda_[gidx];
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

		// Apply to all dimensions.
		for (size_t gidx = 0; gidx < grad.size(); gidx++) {
			grad[gidx] += 2. * v * k_[i*x.size() + gidx];
		}
	}
}

// ------------------- PosDrawSquareMinimizer ----------------

void PosDrawSquareMinimizer::computeGradient(const Fista::Vector &x, Fista::Vector &grad)
{
	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] = lambda_[gidx];
	}
	
	// Initialize the columns to -b.
	for (size_t i = 0; i < b_.size(); i++) {
		drawColumns_[i] = -b_[i];
	}

	drawGrad_ = &grad;

	// Do the pass 1, collect the column data.
	drawPass1_ = true;
	drawing_->drawSimple(x, *this);
	
	// Do the pass 2, add up the data by gradient dimensions.
	drawPass1_ = false;
	drawing_->drawSimple(x, *this);

	drawGrad_ = nullptr;
}

void PosDrawSquareMinimizer::drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
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

#if 0
// this implementation produces a too optimiztic result
void PosDrawSquareMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	// A copy of PosMatSquareMinimizer::computeLOptimistic()

	// See the explanation of estimation for L above.
	Fista::Vector grad_at_0(lambda_.size());
	Fista::Vector grad_at_1(lambda_.size());
	const Fista::Vector input_0(lambda_.size(), 0.);
	const Fista::Vector input_1(lambda_.size(), 1.);

	computeGradient(input_0, grad_at_0);
	computeGradient(input_1, grad_at_1);

	double l = 0.;
	for (size_t i = 0; i < grad_at_0.size(); i++) {
		double v = grad_at_1[i] - grad_at_0[i];
		l += v * v;
	}
	l = sqrt(l / grad_at_0.size());

	inv_L_ = 1./l;
}
#endif

void PosDrawSquareMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	// A copy of PosMatSquareMinimizer::computeLMaxGradient()
	// done via drawing.
	// The computation is the same as computing the gradient but with b=0
	// and x=[all 1s]

	Fista::Vector lpart(x.size(), 0.);
	drawGrad_ = &lpart;

	drawColumns_.assign(b_.size(), 0);
	Fista::Vector input_1(x.size(), 1.);

	// Pass 1.
	drawPass1_ = true;
	drawing_->drawSimple(input_1, *this);

	// Pass 2.
	drawPass1_ = false;
	drawing_->drawSimple(input_1, *this);

	drawGrad_ = nullptr;

	double l = 0.;
	for (int j = 0; j < x.size(); j++) {
		double v = lpart[j];
		// printf("  Lpart[%d] = %f\n", j, v);
		if (v > l)
			l = v;
	}
	// printf("L = %f\n", l);
	inv_L_ = 1./l;
}

// ------------------- MultiGradientMinimizer ----------------

bool MultiGradientMinimizer::compute(Vector &y, Vector &x, int step)
{
	// Temporarily put the gradient into x.
	computeGradient(y, x);

	if (step == 0 || reinit_L_every_ > 0 && step % reinit_L_every_ == 0) {
		inv_L_.resize(x.size());
		// Auto-compute L and 1/L.
		initializeL(y, x, step);
	}

	// This is really norm2(grad_fg(x))^2.
	double gradnorm = 0.;
	// Do the step by dimensions.
	for (size_t i = 0; i < x.size(); i++) {
		double gd = x[i]; // gradient dimension
		double xval = y[i] - inv_L_[i] * gd;
		if (xval <= 0.) {
			xval = 0.;
			if (positive_ && gd > 0.)  {
				// If we've hit the boundary of x>=0 and the gradient points
				// further downwards, we can't go there, so consider this
				// equivalent to gradient dimension becoming 0.
				gd = 0.;
			}
		}
		if (upperBounded_ && xval > upperBound_) {
			xval = upperBound_;
			gd = 0.;
		}
		x[i] = xval;
		gradnorm += gd*gd;
		verbose_ && printf("   p_L[%zd] (%f -> %f); gd=%f\n", i, y[i], xval, gd);
	}
	verbose_ && printf("   norm=%f\n", sqrt(gradnorm));

	// True if length (norm) of the gradient is less than epsilon.
	return gradnorm <= eps2_;
}

// ------------------- MultiPosMatSquareMinimizer ------------

MultiPosMatSquareMinimizer::MultiPosMatSquareMinimizer(
	const Fista::Vector &k,
	const Fista::Vector &b,
	const Fista::Vector &lambda,
	double eps)
	: MultiGradientMinimizer(eps, /*positive*/ true, /*reinit_L_every*/ 0),
	k_(k), b_(b), lambda_(lambda)
{ }

void MultiPosMatSquareMinimizer::computeGradient(const Fista::Vector &x, Fista::Vector &grad)
{
	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] = lambda_[gidx];
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

		// Apply to all dimensions.
		for (size_t gidx = 0; gidx < grad.size(); gidx++) {
			grad[gidx] += 2. * v * k_[i*x.size() + gidx];
		}
	}
}

// This is a further development of the idea from
// PosMatSquareMinimizer::computeLMaxGradient(): if we're considering each
// dimension x_i separately, why do we care to have a common L?
// We're only interested in the movement in each dimension not overshooting
// the minimum, so we can pick a separate L for each dimension.
// Yeah, it won't exactly be a "gradient descent", since the gradient
// will get skewed by these different Ls, but who cares if it works.
//
// The tests in t_mmsq_fista_v1.cpp explore the examples of this.
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
void MultiPosMatSquareMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
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
		if (v == 0.) {
			inv_L_[j] = 1.;
		} else {
			inv_L_[j] = 1. / v;
		}
		// printf("  Lpart[%d] = %f\n", j, v);
	}
}

// ------------------- MultiPosDrawSquareMinimizer -----------

// Same as in PosDrawSquareMinimizer
void MultiPosDrawSquareMinimizer::computeGradient(const Fista::Vector &x, Fista::Vector &grad)
{
	// Start with g(x). Since we assume x >= 0, grad_g(x) is constant.
	for (size_t gidx = 0; gidx < grad.size(); gidx++) {
		grad[gidx] = lambda_[gidx];
	}
	
	// Initialize the columns to -b.
	for (size_t i = 0; i < b_.size(); i++) {
		drawColumns_[i] = -b_[i];
	}

	drawGrad_ = &grad;

	// Do the pass 1, collect the column data.
	drawPass1_ = true;
	drawing_->drawSimple(x, *this);
	
	// Do the pass 2, add up the data by gradient dimensions.
	drawPass1_ = false;
	drawing_->drawSimple(x, *this);

	drawGrad_ = nullptr;
}

// Same as in PosDrawSquareMinimizer
void MultiPosDrawSquareMinimizer::drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
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

// Slightly different than in PosDrawSquareMinimizer, saves the whole vector of 1/L.
// See MultiPosMatSquareMinimizer::initializeL() for the general explamation.
void MultiPosDrawSquareMinimizer::initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step)
{
	// A copy of PosMatSquareMinimizer::computeLMaxGradient()
	// done via drawing.
	// The computation is the same as computing the gradient but with b=0
	// and x=[all 1s]

	// Temporarily store L in inv_L_, then invert.
	inv_L_.assign(inv_L_.size(), 0.);
	drawGrad_ = &inv_L_;

	drawColumns_.assign(b_.size(), 0);
	Fista::Vector input_1(x.size(), 1.);

	// Pass 1.
	drawPass1_ = true;
	drawing_->drawSimple(input_1, *this);

	// Pass 2.
	drawPass1_ = false;
	drawing_->drawSimple(input_1, *this);

	drawGrad_ = nullptr;

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

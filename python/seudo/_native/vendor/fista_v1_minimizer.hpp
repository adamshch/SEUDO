#ifndef __FISTA_V1_MINIMIZER_HPP__
#define __FISTA_V1_MINIMIZER_HPP__

#include <memory>
#include "fista.hpp"

// Collection of minimizer subclasses that can be used with FISTA.

namespace Fista {

// A specialization or LipschitzMinimizer that computes p_L(y) from the
// gradients of the original functions. For the more concrete examples that
// were used to derive it, see the comments to SquareMinimizer2 in t_fista.cpp
// and PosMatSquareMinimizer below. The general logic goes as
// follows:
//
// We're looking (overall) to minimize the sum (f(x) + g(x)), where
//   g(x) = la * abs(x)
// ("la" stands for "lambda", a vector of coefficients for LASSO), and f(x) is
// some smooth function. [Note that the other comments that describe x by
// its elements like x_0 and x_1, use the notation |x_0| for absolute values,
// but it would be a wrong notation for the whole vector, so here it's abs(x)
// instead].
//
// Then the gradient g(x) is:
//   grad_g(x) = [
//     dg/dx_i = {
//       x_i >= 0: la_i,
//       x_i < 0: -la_i,
//     } = (+-la_i),
//   ]
//
// Then if we denote M = (y - (1/L)grad_f(y)),
// p_L(M) = argmin(x){ sum( la_i*abs(x_i) ) - L/2 * sum( (x_i - M_i)^2 ) }
//        [expand squares]
//        = argmin(x){ sum( la_i*abs(x_i) ) + L/2 * sum( - x_i^2 + 2*M_i*x_i - M_i^2 ) }
//        [get rid of constants that don't affect argmin]
//        = argmin(x){ sum( la_i*abs(x_i) ) + L/2 * sum( - x_i^2 + 2*M_i*x_i ) }
// then to find the minimizer, we find the zero for the gradient:
//   grad_p(x) = [
//     dp/dx_i = (+-la_i) + L/2 * ( -2*x_i + 2*M_i ) = -L*x_i + L*M_i + (+-la_i),
//   ]
// for grad_p(x) = 0, we have the point
//   x_i = M_i - (+-la_i)/L
//       = y_i - (1/L)grad_f(y_i) - (1/L)grad_g(y_i)
//       = y_i - (1/L)grad_fg(y_i)
// where grad_fg(x) = grad_f(x) + grad_g(x)
//
// The same value happens to work well in one more way. The "classic" ISTA
// does the step by -(1/L)grad_f(y) and then starts looking in this vicinity
// to do the LASSO part and do another step of size (+-la_i) in every dimension
// to minimize further. How do we find, in which direction (+ or -) to go
// for every dimension? A simple-minded approach would be to try every point,
// but this creates a combinatory explosion. Instead we can look at each dimension
// separately, find the gradient grad_fg(x) in it, and step down the gradient.
// But that's very much the same thing as replacing the "classic" ISTA step
// with -(1/L)grad_fg(y). There is a difference of computing the gradient at
// the point y instead of x but since the original step is small, the result
// is very close. Since this replaces two consecutive steps in potentially
// different directions with one step, the result is not quite as good, and
// the search gets slightly longer, but this allows to compute the gradient
// only once per step, so each step gets almost twice cheaper.
class GradientMinimizer : public LipschitzMinimizer
{
public:
	// @param inv_L - the value of 1/L ("inverse L"), if it's known. If
	//   unknown at construction time, pass 0., and update the field later.
	// @param eps - epsilon for detecting the stopping point, the
	//   stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param positive - expect and enforce x >= 0, since this is a common
	//   class of problems.
	// @param reinit_L_every - if > 0, initializeL() will called every
	//   so many steps (in addition to the call on step 0).
	//
	// This class also supports an upper bound for x, but since it's less typical,
	// it doesn't get enabled from the constructor. Instead set the fields
	// directly.
	GradientMinimizer(
		double inv_L,
		double eps,
		bool positive,
		int reinit_L_every = -1)
		: inv_L_(inv_L), eps2_(eps * eps), positive_(positive),
		reinit_L_every_(reinit_L_every), upperBounded_(false), upperBound_(0.)
	{ }

	// From LipschitzMinimizer, the common implementation based on the custom
	// gradient computation.
	virtual bool compute(Vector &y, Vector &x, int step);

	// Compute gradient of (f(x) + g(x)), denoted above as grad_fg(x).
	// @param x - input values
	// @param grad - produced gradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad) = 0;

	// Gives the subclass a chance to compute L on the step 0, 
	// and then if reinit_L_every_ > 0 then every so many steps.
	// 
	// By default does nothing.
	//
	// @param x - the initial point of the step
	// @param grad - the computed gradient at this point (this is a minor optimization,
	//		since both the normal course of computation and the estimation of L need it,
	//		it gets computed up front and then used for both).
	// @param step - the current step, starting from 0, in case if the estimation
	//      depends on it
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	// Value of 1/L.
	double inv_L_;
	// eps^2 for detecting the stopping point
	double eps2_;
	// Flag: assume always x >= 0.
	bool positive_;
	// On step 0 and then if reinit_L_every_ > 0 then every so many steps
	// initializeL() will be called to get a chance to get the new estimate
	// of L.
	int reinit_L_every_;
	// Definition of the upper bound. 
	bool upperBounded_;
	double upperBound_;
};

// A specialization of GradientMinimizer with computation of L for the
// square functions (where the function is convex and the gradient smoothly
// trends to 0 when approaching minimum, never growing).
class GradientSquareMinimizer : public GradientMinimizer
{
public:
	// @param eps - epsilon for detecting the stopping point, the
	//   stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param positive - expect and enforce x >= 0, since this is a common
	//   class of problems.
	// @param reinit_L_every - if > 0, initializeL() will called every
	//   so many steps (in addition to the call on step 0). Should not really
	//   be needed for a square function, because its L should be constant.
	//
	// This class also supports an upper bound for x, but since it's less typical,
	// it doesn't get enabled from the constructor. Instead set the fields
	// directly.
	GradientSquareMinimizer(
		double eps,
		bool positive,
		int reinit_L_every = -1)
		: GradientMinimizer(/*inv_L*/ 0.,  eps, positive, reinit_L_every)
	{ }

	// from GradientMinimizer
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);
};

// Minimizer for
//   f(x) = norm2(k*x - b)^2
//   g(x) = la*x
// where k is a matrix, x, b, and la are vectors, and x >= 0. Since x >= 0,
// we can assume |x| = x, and drop the absolute value computation from g(x).
//
// Start with visualising the simple case with vectors of size 2 and matrix 2*2,
// which then becomes obvious to generalize.
//
// The simple case functions:
//   f(x) = (k_00*x_0 + k_01*x_1 - b0)^2 + (k_10*x_0 + k_11*x_1 - b_1)^2
//   g(x) = la_0*x_0 + la_1*x_1
//
// Then following the lead from SquareMinimizer2 in t_fista.cpp :
//   grad_f(x) = [
//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
//     df/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1),
//   ]
//   grad_g(x) = [
//     dg/dx_0 = la_0,
//     dg/dx_1 = la_1,
//   ]
//   grad_fg(x) = grad_f(x) + grad_g(x) = [
//     dfg/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1) + la_0,
//     dfg/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1) + la_1,
//   ]
// (note that the expressions inside each "column" of parentheses are the same!).
//
// Then if we denote M = (y - (1/L)grad_f(y)),
// p_L(M) = argmin(x){ la_0*|x_0| + la_1*|x_1| - L/2 * ( (x_0 - M_0)^2 + (x_1 - M_1)^2 ) }
//        = argmin(x){ la_0*|x_0| + la_1*|x_1| + L/2 * ( - x_0^2 + 2*M_0*x_0 - M_0^2 - x_1^2 + 2*M_1*x_1 - M_1^2 ) }
//        = argmin(x){ la_0*|x_0| + la_1*|x_1| + L/2 * ( - x_0^2 + 2*M_0*x_0 - x_1^2 + 2*M_1*x_1 ) }
// then to find the minimizer
// grad_p(x) = [
//   dp/dx_0 = la_0 + L/2 * ( -2*x_0 + 2*M_0 ) = -L*x_0 + L*M_0 + la_0,
//   dp/dx_1 = la_1 + L/2 * ( -2*x_1 + 2*M_1 ) = -L*x_1 + L*M_1 + la_1,
// ]
// for grad_p(x) = 0, we have the point
//   x_0 = M_0 - la_0/L
//       = y_0 - (1/L)grad_f(y_0) - (1/L)grad_g(y_0)
//       = y_0 - (1/L)grad_fg(y_0)
//   x_1 = M_1 - la_1/L
//       = y_1 - (1/L)grad_f(y_1) - (1/L)grad_g(y_1)
//       = y_1 - (1/L)grad_fg(y_1)
//
// The approach described here for a matrix of 2*2 gets easily generalized to a
// matrix of K*N (where K is the size of x, and N is the size of b):
//
//   grad_fg(x) = [
//     {0<=m<N}
//     dfg/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) ) + la_m
//   ]
//
// Evaluating
//   L = max( ||grad_f(x) - grad_f(y)|| / ||x - y|| )
// I can see multiple approaches.
//
// The safest estimation for L is implemented in computeLPessimistic(), looking
// at expression 
//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
// we can eliminate b_i because they are constant and generalize for matrix dimension
// K*N (the expression assumes that as in C++ the fastest-changing index of k is the last one):
//     df/dx_i = 2*k_0i*(k_00*x_0 + k_01*x_1 + .. + k_0K*x_K) + ... + 2*k_Ni*(k_N0*x_0 + k_N1*x_1 + ... + k_NK*x_K),
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
// See t_msq_fista.cpp for an example of usage.
class PosMatSquareMinimizer : public GradientSquareMinimizer
{
public:
	// The sizes of all vectors must match as described.
	// @param k - a flattened matrix K*N, structured as N "rows"
	//		of K entries one after another:
	//            [ ... K ... ]  ...  [ ... K ... ] 
	//      that defines f = sum{i,j}
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	PosMatSquareMinimizer(
		const Fista::Vector &k,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps);

	// Various approaches to compute L. See the description of the
	// logic in the .cpp file. computeLOptimistic() and computeLHessian()
	// are wrong but kept for historic reasons. computeLPessimistic()
	// is correct but overestimates L and causes a slow convergence.
	// computeLMaxGradient() is the good one and the one currently 
	// hardcoded in the implementation.
	// See also MultiGradientMinimizer::initializeL().
	double computeLOptimistic();
	double computeLPessimistic();
	double computeLHessian();
	double computeLMaxGradient();

	// from GradientMinimizer
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad);

protected:
	Fista::Vector k_;
	Fista::Vector b_;
	Fista::Vector lambda_;
};

// Like PosMatSquareMinimizer but instead of computing gradient by
// matrix, computes it from code that "draws" the transformation.
// Since the most obvious usage is to do it with the images, they
// are called pixels in the comments, but they can be anything.
// This allows to avoid building the large matrices.
//
// See BlurDrawMinimizer in t_dsq_fista.cpp for an example of usage.
class PosDrawSquareMinimizer : public GradientMinimizer, Drawable
{
public:
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawing - specification of the dependency matrix via
	//      a "drawing". Its draw() method will be called twice per
	//      gradient computation.
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param highest_k - the highest possible value of coefficient k_ij that
	//      determines a dependency between the two pixels.
	//      (gets really used only if computeLOptimistic == false)
	// @param computeLOptimistic - estimate L using the logic of
	//      computeLOptimistic, otherwise of computeLPessimistic
	PosDrawSquareMinimizer(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps)
		: GradientMinimizer(/*inv_L*/ 1., eps, /*positive*/ true),
		Drawable(wd, ht),
		drawing_(drawing), b_(b), lambda_(lambda),
		drawColumns_(b.size()), drawGrad_(nullptr)
	{ }

	// Same, but also with an upper bound. The upper bound can also
	// be added manually afterwards but more convenient here.
	//
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawing - specification of the dependency matrix via
	//      a "drawing". Its draw() method will be called twice per
	//      gradient computation.
	// @param upper_bound - the upper bound for X
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param highest_k - the highest possible value of coefficient k_ij that
	//      determines a dependency between the two pixels.
	//      (gets really used only if computeLOptimistic == false)
	// @param computeLOptimistic - estimate L using the logic of
	//      computeLOptimistic, otherwise of computeLPessimistic
	PosDrawSquareMinimizer(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		double upper_bound,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps)
		: GradientMinimizer(/*inv_L*/ 1., eps, /*positive*/ true),
		Drawable(wd, ht),
		drawing_(drawing), b_(b), lambda_(lambda),
		drawColumns_(b.size()), drawGrad_(nullptr)
	{
		upperBounded_ = true;
		upperBound_ = upper_bound;
	}

	// from GradientMinimizer
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad);
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

	// From Drawable.
	virtual void drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);

protected:
	Fista::Vector b_;
	Fista::Vector lambda_;
	std::shared_ptr<Drawing> drawing_;

	// State of the current drawing.

	// This is pass 1 of drawing, otherwise pass 2.
	bool drawPass1_;
	// Collected values in parentheses in each "column" of terms for gradient. These
	// "column" values are the same for all "rows", so they get written in pass 1 and
	// read in pass 2.
	//   grad_f(x) = [
	//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
	//     df/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1),
	//                       ^^^^^^^^ "column" ^^^^^^^            ^^^^^^ "column" ^^^^^^^^^
	//   ]
	Fista::Vector drawColumns_;
	Fista::Vector *drawGrad_; // Gradient to compute by drawing.
};

// A version of gradient minimizer that computes a separate value of 1/L
// for each dimension of the gradient (or equivalently, for each dimension
// of X).
class MultiGradientMinimizer : public LipschitzMinimizer
{
public:
	// @param eps - epsilon for detecting the stopping point, the
	//   stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param positive - expect and enforce x >= 0, since this is a common
	//   class of problems.
	// @param reinit_L_every - if > 0, initializeL() will called every
	//   so many steps (in addition to the call on step 0).
	//
	// This class also supports an upper bound for x, but since it's less typical,
	// it doesn't get enabled from the constructor. Instead set the fields
	// directly.
	MultiGradientMinimizer(
		double eps,
		bool positive,
		int reinit_L_every = -1)
		: eps2_(eps * eps), positive_(positive),
		reinit_L_every_(reinit_L_every), upperBounded_(false), upperBound_(0.)
	{ }

	// From LipschitzMinimizer, the common implementation based on the custom
	// gradient computation.
	virtual bool compute(Vector &y, Vector &x, int step);

	// Compute gradient of (f(x) + g(x)), denoted above as grad_fg(x).
	// @param x - input values
	// @param grad - produced gradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad) = 0;

	// Gives the subclass a chance to compute L on the step 0, 
	// and then if reinit_L_every_ > 0 then every so many steps.
	// 
	// @param x - the initial point of the step
	// @param grad - the computed gradient at this point (this is a minor optimization,
	//		since both the normal course of computation and the estimation of L need it,
	//		it gets computed up front and then used for both).
	// @param step - the current step, starting from 0, in case if the estimation
	//      depends on it
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step) = 0;

protected:
	// Values of 1/L per dimension. The size will be set before calling initializeL().
	Fista::Vector inv_L_;
	// eps^2 for detecting the stopping point
	double eps2_;
	// Flag: assume always x >= 0.
	bool positive_;
	// On step 0 and then if reinit_L_every_ > 0 then every so many steps
	// initializeL() will be called to get a chance to get the new estimate
	// of L.
	int reinit_L_every_;
	// Definition of the upper bound. 
	bool upperBounded_;
	double upperBound_;
};

// Similar to PosMatSquareMinimizer but uses individual values
// of 1/L per dimension.
//  
// Optimizes the square function
//   f(x) = sum{0<=i<N}(( sum{0<=j<K}(k_ij*x_j) - b_i)^2 )
// The function is specified as a matrix of values k_ij.
// Here N is the number of "output variables" that are fitted to the
// expected result b_i, and K is the number of "input variables" x_j.
class MultiPosMatSquareMinimizer : public MultiGradientMinimizer
{
public:
	// See the meaning of arguments in the comment to the class.
	// The sizes of all vectors must match as described.
	// @param k - a flattened matrix K*N, structured as N "rows"
	//		of K entries one after another:
	//            [ ... K ... ]  ...  [ ... K ... ] 
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	MultiPosMatSquareMinimizer(
		const Fista::Vector &k,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps);

	double computeL();

	// from MultiGradientMinimizer,
	// same as in PosMatSquareMinimizer
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad);

	// from MultiGradientMinimizer,
	// Similar to PosMatSquareMinimizer::computeLMaxGradient().
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	Fista::Vector k_;
	Fista::Vector b_;
	Fista::Vector lambda_;
};

// Similar to PosDrawSquareMinimizer but uses individual values
// of 1/L per dimension.
class MultiPosDrawSquareMinimizer : public MultiGradientMinimizer, Drawable
{
public:
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawing - specification of the dependency matrix via
	//      a "drawing". Its draw() method will be called twice per
	//      gradient computation.
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param highest_k - the highest possible value of coefficient k_ij that
	//      determines a dependency between the two pixels.
	//      (gets really used only if computeLOptimistic == false)
	// @param computeLOptimistic - estimate L using the logic of
	//      computeLOptimistic, otherwise of computeLPessimistic
	MultiPosDrawSquareMinimizer(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps)
		: MultiGradientMinimizer(eps, /*positive*/ true),
		Drawable(wd, ht),
		drawing_(drawing), b_(b), lambda_(lambda),
		drawColumns_(b.size()), drawGrad_(nullptr)
	{ }

	// Same, but also with an upper bound. The upper bound can also
	// be added manually afterwards but more convenient here.
	//
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawing - specification of the dependency matrix via
	//      a "drawing". Its draw() method will be called twice per
	//      gradient computation.
	// @param upper_bound - the upper bound for X
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param eps - epsilon for detecting the stopping point, the
	//      stop will be triggered when norm2(grad_fg(x)) <= eps
	// @param highest_k - the highest possible value of coefficient k_ij that
	//      determines a dependency between the two pixels.
	//      (gets really used only if computeLOptimistic == false)
	// @param computeLOptimistic - estimate L using the logic of
	//      computeLOptimistic, otherwise of computeLPessimistic
	MultiPosDrawSquareMinimizer(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		double upper_bound,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		double eps)
		: MultiGradientMinimizer(eps, /*positive*/ true),
		Drawable(wd, ht),
		drawing_(drawing), b_(b), lambda_(lambda),
		drawColumns_(b.size()), drawGrad_(nullptr)
	{
		upperBounded_ = true;
		upperBound_ = upper_bound;
	}

	// from GradientMinimizer
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad);
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

	// From Drawable.
	virtual void drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);

protected:
	Fista::Vector b_;
	Fista::Vector lambda_;
	std::shared_ptr<Drawing> drawing_;

	// State of the current drawing.

	// This is pass 1 of drawing, otherwise pass 2.
	bool drawPass1_;
	// Collected values in parentheses in each "column" of terms for gradient. These
	// "column" values are the same for all "rows", so they get written in pass 1 and
	// read in pass 2.
	//   grad_f(x) = [
	//     df/dx_0 = 2*k_00*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_10*(k_10*x_0 + k_11*x_1 - b_1),
	//     df/dx_1 = 2*k_01*(k_00*x_0 + k_01*x_1 - b_0) + 2*k_11*(k_10*x_0 + k_11*x_1 - b_1),
	//                       ^^^^^^^^ "column" ^^^^^^^            ^^^^^^ "column" ^^^^^^^^^
	//   ]
	Fista::Vector drawColumns_;
	Fista::Vector *drawGrad_; // Gradient to compute by drawing.
};

}; // namespace Fista

#endif // __FISTA_V1_MINIMIZER_HPP__

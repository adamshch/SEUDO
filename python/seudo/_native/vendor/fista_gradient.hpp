#ifndef __FISTA_MINIMIZER_HPP__
#define __FISTA_MINIMIZER_HPP__

#include <memory>
#include "fista.hpp"

// Collection of minimizer subclasses that can be used with FISTA.

namespace Fista {

// defined in par_fista.hpp
class GradientThread;

// A specialization of ScaledGradient with computation of L for the
// square functions (where the function is convex and the gradient smoothly
// trends to 0 when approaching minimum, never growing), like PosMatSquareMinimizer of v1.
//
// This is used to minimize
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
// Then following the transformation similar to one described in SquareMinimizer2
// in t_fista_v1.cpp :
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
// So what this class needs to compute is (1/L)grad_fg().
//
// The approach described here for a matrix of 2*2 gets easily generalized to a
// matrix of K*N (where K is the size of x, and N is the size of b):
//
//   grad_fg(x) = [
//     {0<=m<N}
//     dfg/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) ) + la_m
//   ]
//
// See the discussion of L in the implementations of initializeL().
class PosMatScaledGradient : public ScaledGradient
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
	PosMatScaledGradient(
		const Fista::Vector &k,
		const Fista::Vector &b,
		const Fista::Vector &lambda);

	// from ScaledGradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step);

	// from ScaledGradient
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	// Value of 1/L.
	double inv_L_ = 0.;
	// parameters from constructor
	Fista::Vector k_;
	Fista::Vector b_;
	Fista::Vector lambda_;
};

// The common base class for quadratic functions that compute their gradients
// by "drawing" the transformation, expressing a sparse dependency matrix.
// This allows to avoid building the large matrices.
// Since the most obvious usage is to do it with the images, the elements of
// the target function are called pixels in the comments, but they can be
// anything.
class BaseDrawScaledGradient : public ScaledGradient, public Drawable, public SelfDrawable
{
public:
	typedef CustomDrawing<BaseDrawScaledGradient> DrawingToMe;

	// The basic constuctor with the general drawing.
	//
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawing - specification of the dependency matrix via
	//      a "drawing". Its draw() method will be called twice per
	//      gradient computation.
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	BaseDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int nthreads);

	// Constructor that uses the drawings templatized to passes 1 and 2,
	// so that they avoid the overhead of calling a virtual function on
	// each pixel.
	//
	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawingPass1 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 1 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass1>
	// @param drawingPass2 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 2 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass2>
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	BaseDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<DrawingToMe> drawingPass1,
		std::shared_ptr<DrawingToMe> drawingPass2,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int nthreads);

	// Compute the gradient as-is, not scaled by L.
	// @param zeroed - true if this is a computation that uses zeroes instead
	//    of b and lambda (to compute L), false for a normal gradient computation
	// @param x - input values
	// @param grad - produced gradient, unscaled
	void computeUnscaledGradient(
		bool zeroed,
		const Fista::Vector &x, Fista::Vector &grad);

	// From Drawable. Used with plain Drawing.
	virtual void drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);

	// From SelfDrawable. Knows how to select either plain or templatized drawings.
	virtual void selfDraw(const Fista::Vector &input, Drawing::Limits &limits);

	// Mainly for testing.
	typedef std::vector<Drawing::Limits> LimitsVector;
	const LimitsVector &limitsPass1() const
	{
		return limits_pass_1_;
	}
	const LimitsVector &limitsPass2() const
	{
		return limits_pass_2_;
	}

	// Templatized drawing.
	friend void 
	drawPixelPass1(BaseDrawScaledGradient &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);
	friend void 
	drawPixelPass2(BaseDrawScaledGradient &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);

protected:
	// Common initialization for all constructors.
	void initialize(
		const Fista::Vector &b,
		int nthreads,
		bool parallel);

protected:
	// Use either the general drawing or custom drawing by passes.
	std::shared_ptr<Drawing> drawing_;
	std::shared_ptr<DrawingToMe> drawingPass1_;
	std::shared_ptr<DrawingToMe> drawingPass2_;

	// negated b (pre-negation makes the computation easier)
	Fista::Vector neg_b_;
	Fista::Vector lambda_;

	// Thread information.

	// Note that potentially the sizes of vectors for pass 1 and pass 2
	// may be different, if the thread count requested is higher than
	// can be partitioned by either dimension.
	// Pass 1 partitions by destination.
	LimitsVector limits_pass_1_;
	// Pass 2 partitions by input.
	LimitsVector limits_pass_2_;
	// Thread size is the max of pass 1 and pass 2 minus 1 (since the caller's
	// thread is used as one of the threads).
	std::vector<std::shared_ptr<GradientThread>> threads_;

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
	//  or in generalized form:
	//      df/dx_m = 2 * sum{0<=j<N}( k_jm * (sum{0<=i<K}(k_ji * x_i) - b_j) ) + la_m
	//  the "columns" are (sum{0<=i<K}(k_ji * x_i) - b_j)
	Fista::Vector drawColumns_;
	Fista::Vector *drawGrad_; // Gradient to compute by drawing.
};

// Drawing of pass 1 for templatized implementation
inline void 
drawPixelPass1(BaseDrawScaledGradient &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
{
	int to_idx = to_y * dest.wd_ + to_x;
	dest.drawColumns_[to_idx] += k * x[from_idx];
}

// Drawing of pass 2 for templatized implementation
inline void 
drawPixelPass2(BaseDrawScaledGradient &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
{
	int to_idx = to_y * dest.wd_ + to_x;
	(*dest.drawGrad_)[from_idx] += 2. * k * dest.drawColumns_[to_idx];
}

// Like PosMatScaledGradient but instead of computing gradient by
// matrix, computes it from code that "draws" the transformation.
// Since the most obvious usage is to do it with the images, they
// are called pixels in the comments, but they can be anything.
// This allows to avoid building the large matrices.
class PosDrawScaledGradient : public BaseDrawScaledGradient
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
	// @param reinit_L_every - this class supports recomputing L dynamically
	//      from the empirical values at the recent points, like TFOCS does it. Set
	//      to 1 to do recomputation on every step like TFOCS. Setting it to a
	//      larger value to average over multiple steps is not a good idea,
	//      since then it ends up iterating over multiple steps with a value of L
	//      that is too low and takes longer to converge.
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	PosDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int reinit_L_every = 0,
		int nthreads = 1)
		: BaseDrawScaledGradient(wd, ht, drawing, b, lambda, nthreads)
	{ 
		reinit_L_every_ = reinit_L_every;
	}

	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawingPass1 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 1 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass1>
	// @param drawingPass2 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 2 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass2>
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param reinit_L_every - this class supports recomputing L dynamically
	//      from the empirical values at the recent points, like TFOCS does it. Set
	//      to 1 to do recomputation on every step like TFOCS. Setting it to a
	//      larger value to average over multiple steps is not a good idea,
	//      since then it ends up iterating over multiple steps with a value of L
	//      that is too low and takes longer to converge.
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	PosDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<DrawingToMe> drawingPass1,
		std::shared_ptr<DrawingToMe> drawingPass2,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int reinit_L_every = 0,
		int nthreads = 1)
		: BaseDrawScaledGradient(wd, ht, drawingPass1, drawingPass2, b, lambda, nthreads)
	{ 
		reinit_L_every_ = reinit_L_every;
	}

	// from ScaledGradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step);
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	// Value of 1/L.
	double inv_L_ = 0.;

	// Values from the last computation of L.
	Fista::Vector last_x_;
	Fista::Vector last_grad_;
};

// Similar to PosMatScaledGradient but uses individual values
// of 1/L per dimension.
class MultiPosMatScaledGradient : public ScaledGradient
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
	MultiPosMatScaledGradient(
		const Fista::Vector &k,
		const Fista::Vector &b,
		const Fista::Vector &lambda);

	// from ScaledGradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step);

	// from ScaledGradient
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	// Values of 1/L per dimension. The size will be set before calling initializeL().
	Fista::Vector inv_L_;
	// parameters from constructor
	Fista::Vector k_;
	Fista::Vector b_;
	Fista::Vector lambda_;
};

// Similar to PosDrawScaledGradient but uses individual values
// of 1/L per dimension.
class MultiPosDrawScaledGradient : public BaseDrawScaledGradient
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
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	MultiPosDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<Drawing> drawing,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int nthreads = 1)
		: BaseDrawScaledGradient(wd, ht, drawing, b, lambda, nthreads)
	{ }

	// @param wd - width of matrix that is simulated by drawing
	// @param ht - height of matrix that is simulated by drawing
	// @param drawingPass1 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 1 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass1>
	// @param drawingPass2 - specification of the dependency matrix via
	//      a "drawing", customized to to the pass 2 of drawing, such as
	//      SomeCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass2>
	// @param b - target values, of size N
	// @param lambda - LASSO coefficients, of size K, which matches
	// 		the size of vector X
	// @param nthreads - number of threads for parallel computation
	//       (it may be reduced if way too large, or always 1 for non-parallel
	//       drawings)
	MultiPosDrawScaledGradient(
		int wd,
		int ht,
		std::shared_ptr<DrawingToMe> drawingPass1,
		std::shared_ptr<DrawingToMe> drawingPass2,
		const Fista::Vector &b,
		const Fista::Vector &lambda,
		int nthreads = 1)
		: BaseDrawScaledGradient(wd, ht, drawingPass1, drawingPass2, b, lambda, nthreads)
	{ }

	// from ScaledGradient
	virtual void computeGradient(const Fista::Vector &x, Fista::Vector &grad, int step);
	virtual void initializeL(const Fista::Vector &x, const Fista::Vector &grad, int step);

protected:
	// Values of 1/L per dimension. The size will be set before calling initializeL().
	Fista::Vector inv_L_;
};

}; // namespace Fista

#endif // __FISTA_MINIMIZER_HPP__

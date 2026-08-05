#ifndef __FISTA_V1_HPP__
#define __FISTA_V1_HPP__

#include <vector>
#include <memory>

#include "fista_types.hpp"

// The FISTA method implementation, the first version patterned closely on
// https://people.rennes.inria.fr/Cedric.Herzet/Cedric.Herzet/Sparse_Seminar/Entrees/2012/11/12_A_Fast_Iterative_Shrinkage-Thresholding_Algorithmfor_Linear_Inverse_Problems_(A._Beck,_M._Teboulle)_files/Breck_2009.pdf
namespace Fista {

// Vector is defined in fista_types.hpp.

// This class implements the function p_L, computing the value from p.189:
// 
// p_L(y) = argmin(x){ g(x) - L/2 || x - (y - (1/L)grad_f(y)) ||^2 }
//
// The subclasses embed the values of L (Lipschitz constant of the gradient
// of function f()), of the function g() and of gradient of function f().
class LipschitzMinimizer
{
public:
	LipschitzMinimizer()
		: verbose_(false)
	{ }

	virtual ~LipschitzMinimizer();

	// Compute x = p_L(y)
	// @param y - the input, should generally match the argument size expected by
	//   the function embedded into this instance. If the size doesn't match, the
	//   implementation can either silently resize it (that's why this argument is not
	//   const), or throw an exception (which is probably a better option), or crash.
	// @param x - place to return the computed minimizer, its size will be equal to y
	// @param step - step number, gets used for the initial and periodic auto-tuning;
	//	 starts from 0
	// @return - the indication that the algorithm has achieved a close enough
	//   point and may stop (if returns true). How exactly it does that is up
	//   to the implementation.  Since the size of the step is based on the
	//   gradient, that is a convenient value that can be used in the
	//   estimation of closeness (don't forget to include g(x) into the
	//   gradient, or it may never get close enough).  An implementation may
	//   also always return false, then the algorithm will stop by the count of
	//   steps.
	virtual bool compute(Vector &y, Vector &x, int step) = 0;

public:
	// Set to true to get the logging of the details from teh subclasses
	// (provided that they check this field and do the logging).
	bool verbose_;

};

// Run one step of FISTA algorithm.
// The steps on p. 193 become more convenient when rearranged, placing 4.1 after 4.3:
//   (4.2) t_{k+1} = ...
//   (4.3) y_{k+1} = ...
//   (4.1) x_{k+1} = ...
// So this function produces the new t_{k+1} and x_{k+1} from t_k, x_{k-1}, and x_k.
//
// @param pl - encapsulation of the specific function to optimize
// @param t - on input t_k, on output t_{k+1}; the "braking coefficient" that
//   determines, how much of the inertia passes through from the previous step.
//   It gets gradually decreased to reduce the "inertia" and "put the brakes on"
//   as the algorithm gets closer to the minimum, to reduce the overshooting.
// @param x_prev - on input x_{k-1}, on output x_k
// @param x - on input x_x, on output x_{k+1}
// @param y - place to temporarily store y_{k+1} (which is also left there on output).
//   This could easily be a local variable in this function, but since there are usually
//   many steps run in a sequence, this allows to avoid the overhead of allocating and
//   freeing a Vector for every step. y gets automatically resized to the correct size,
//   so it can be initially an empty vector.
// @param step - step number, gets used for the initial and periodic auto-tuning;
//	 starts from 0
// @param nrestart - if > 0, resets the parameter t to 1 after each this many steps,
//   "hitting the wall" and restarting the inertia from a stop. This can be used if
//   the algorithm is "circling the drain" too much.
// @return - the Lipschitz minimizer has indicated that a close enough precision has
//   been achieved and the algorithm may stop.
bool oneStep(LipschitzMinimizer &pl, double &t, Vector &x_prev, Vector &x, Vector &y,
	int step, int nrestart);

// Run a fixed number of FISTA steps.
// @param pl - the encapsulation of function to minimize.
// @param x - on input contains the initial vector x_0, on output contains the
//   result x_n.
// @param n - number of steps to take.
// @param nrestart - if > 0, resets the parameter t to 1 after each this many steps,
//   "hitting the wall" and restarting the inertia from a stop. This can be used if
//   the algorithm is "circling the drain" too much.
// @return - true if stopped by the Lipschitz minimizer decising that the precision
//   is close enough, false if ran out of steps.
bool repeatedSteps(LipschitzMinimizer &pl, Vector &x, int n, int nrestart = 0);

}; // namespace Fista

#endif // __FISTA_V1_HPP__


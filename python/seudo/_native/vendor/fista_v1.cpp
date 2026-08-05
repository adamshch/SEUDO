// The FISTA method implementation, v1.

#include <stdio.h>
#include <math.h>
#include "fista_v1.hpp"

namespace Fista {

// ====================== FISTA as such v1 ==============================

LipschitzMinimizer::~LipschitzMinimizer()
{ }

bool oneStep(LipschitzMinimizer &pl, double &t, Vector &x_prev, Vector &x, Vector &y,
	int step, int nrestart)
{
	if (nrestart > 0 && step % nrestart == 0) {
		// "Hit the wall" and restart from a stop.
		t = 1.;
	}

	// step 4.2
	double t_next = (1. + sqrt(1. + t * t * 4.)) / 2.;
	// momentum coefficient
	double nu = (t - 1.) / t_next;
	// t is not used any more, so can assign the new value now
	t = t_next;

	// step 4.3 - the momentum step
	y.resize(x.size()); // just in case
	for (size_t i = 0; i < x.size(); i++) {
		y[i] = x[i] + nu * (x[i] - x_prev[i]);
	}

	// step 4.1
	x.swap(x_prev); // x_k gets moved to x_{k-1}
	y.resize(x.size());
	return pl.compute(y, x, step);
}

bool repeatedSteps(LipschitzMinimizer &pl, Vector &x, int n, int nrestart)
{
	if (n <= 0)
		return false;

	double t = 1.;
	Vector x_next(x.size());
	pl.compute(x, x_next, 0); // computes x_1 from x_0

	bool stopped = false;
	Vector y;
	// the loop goes from x_2 and up
	for (int i = 1; i < n; i++) {
		stopped = oneStep(pl, t, x, x_next, y, i, nrestart);
		if (stopped)
			break;
	}

	x.swap(x_next); // return the latest result
	return stopped;
}


}; // namespace Fista

// Support for parallel computation in FISTA, by computing the gradients
// in parallel.

#include <stdio.h>
#include "par_fista.hpp"

namespace Fista {

// ------------------- GradientThread ------------------------

GradientThread::GradientThread(SelfDrawable *dest)
	: dest_(dest), stop_(false), ev_ready_(true)
{
	start();
}

GradientThread::~GradientThread()
{
	// printf("Stopping GradientThread %p\n", this); fflush(stdout);
	stop_ = true;
	ev_run_.signal();
	join();
	// printf("Joined GradientThread %p\n", this); fflush(stdout);
}

void *
GradientThread::execute()
{
	// printf("Starting GradientThread %p\n", this); fflush(stdout);
	for (;;) {
		ev_run_.wait();
		if (stop_)
			break;

		dest_->selfDraw(*input_, *limits_);
		// printf("Computed GradientThread %p\n", this); fflush(stdout);
		ev_ready_.signal();
	}

	// printf("Exiting GradientThread %p\n", this); fflush(stdout);
	return nullptr;
}

}; // namespace Fista

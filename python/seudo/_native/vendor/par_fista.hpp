#ifndef __PAR_FISTA_HPP__
#define __PAR_FISTA_HPP__

// Support for parallel computation in FISTA, by computing the gradients
// in parallel.

#include "ptwrap.h"
#include "ptwrap2.h"
#include "fista_types.hpp"

namespace Fista {

class GradientThread : public pw::pwthread
{
public:
	// Starts the thread right away.
	// @param dest - destination for drawing (the object that owns this thread)
	GradientThread(SelfDrawable *dest);

	// Stops and joins the thread.
	virtual ~GradientThread();

	// From pw::pwthread, the main logic. Sets ready when started.
	virtual void *execute();

	// Wait for completion of computation.
	void wait()
	{
		ev_ready_.wait();
	}

	// Do a computation for given limits.
	void compute(const Vector *input, Drawing::Limits *limits)
	{
		input_ = input;
		limits_ = limits;
		ev_ready_.reset();
		ev_run_.signal();
	}

public:
	// Context for computing gradients.

	SelfDrawable *dest_; // Destination for drawing.
	bool stop_ = false;  // Don't compute anything, just stop the thread.

	// Input to draw.
	const Vector *input_ = nullptr;
	// Limits for drawing.
	Drawing::Limits *limits_ = nullptr;

	// Control for execution.

	// Signaled by thread when the thread is ready to do more work.
	// The control must reset it before initiating a new run.
	pw::event2 ev_ready_;

	// The control signals to initiate a new run or stop.
	pw::autoevent2 ev_run_;
};

}; // namespace Fista

#endif // __PAR_FISTA_HPP__

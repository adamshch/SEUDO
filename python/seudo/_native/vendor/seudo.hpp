#ifndef __SEUDO_HPP__
#define __SEUDO_HPP__

#include <vector>
#include <string>
#include <memory>
#include "fista.hpp"
#include "fista_types.hpp"
#include "fista_gradient.hpp"

namespace Fista {

// Arguments and methods for SEUDO minimizer.
class Seudo
{
public:
	// Most of values are to be filled after construction.
	// @param wd - width of target image
	// @param ht - height of target image
	Seudo(int wd, int ht)
		: image_(wd, ht)
	{ }

	// Do the FISTA computation on the arguments, updating the output values.
	void compute();

	// Set all weights to the same value. The ROIs must be already initialized.
	void setWeights(double v)
	{
		weights_.assign(rois_.size() + image_.img_.size(), v);
	}

	// Set all lambda to the same value. The ROIs must be already initialized.
	void setLambda(double v)
	{
		lambda_.assign(rois_.size() + image_.img_.size(), v);
	}

	// Set the number of threads for computation.
	void setNumThreads(int n)
	{
		nthreads_ = n;
	}

	// Draw the representation of the current weights into an image
	// (it will be of the same size as image_).
	Fista::DrawableImg draw() const;

	// Compute the norm2 distance between the current weights and the image to match.
	double distance() const;

public:
	// This part is really internal, bur used by the SeudoDrawing,
	// so it's easuer to make public.

	// Shrink the vector to exclude the blobs that don't fit onto the
	// spacing grid. Resizes the destination as needed.
	void shrinkToBlobSpacing(const Vector& src, Vector &dst);
	// Put the weights from the shrunk version back into the
	// full one, zeroing those that were not on the grid.
	void unshrinkToBlobSpacing(const Vector& src, Vector &dst);

	// Round the shrunk sizes up, because even one row of
	// pixels fittimg in means that it should be included.
	int getShrunkWd(int blockSpacing) const
	{
		return (image_.wd_ + blockSpacing - 1) / blockSpacing;
	}
	int getShrunkHt(int blockSpacing) const
	{
		return (image_.ht_ + blockSpacing - 1) / blockSpacing;
	}

public:
	// Image to match.
	Fista::DrawableImg image_;

	// Blob used for detecting the spurious lightings.
	std::shared_ptr<Fista::Sprite> blob_;

	typedef std::vector<std::shared_ptr<Fista::Sprite>> SpriteVector;
	// "Regions of Interest", sprites of the neurons.
	SpriteVector rois_;

	// The initial value of X that gets replaced with computed value.
	// The weights for neurons go first, followed by the weights for
	// blobs centered at each pixel of the image.
	Fista::Vector weights_;

	// Lambda values for each weight.
	Fista::Vector lambda_;

	// During computation, place blobs on a grid of this many pixels increase
	// between the points, both vertically and horizontally. The weights for
	// blobs in between will be filled with 0s.
	// 
	// This allows to speed up the computation by convolving fewer blobs.
	// Increasing this argument reduces the run time at least quadratically by
	// reducing the number of blobs quadratically. But this is also likely to
	// reduce the needed number of steps, so the time is reduced a little
	// faster than quadratically.
	//
	// The precision of the whole computed solution decreases a little (for a
	// reasonably small spacing, within 1/2 of the blob size or so) but it the
	// ROIs seem to be recognized just as well, it's the blobs that get
	// slightly less precise.
	int blobSpacing_ = 1;

	// Precision for optimization. Uses the StopEpsNorm2Rel condition, same as
	// TFOCS mode 1. The combination of this condition with auto-adjusted L
	// causes spurious early stops on overshooting. It would be better to use
	// a higher eps and more stable computation. But that's the baseline
	// to which we get compared, so match it closely.
	double eps_ = 0.001;

	// Maximum number of steps to take.
	int maxSteps_ = 1000;

	// If true, uses MultiPosDrawScaledGradient, if false, uses PosDrawScaledGradient
	// with L recomputed on every step.
	bool multiGrad_ = false;

	// Log debugging information.
	bool verbose_ = false;

	// Use the version of the code that builds a templatized custom drawing
	// function. This is the faster way, so there is no reason not to use
	// it, other than to test the comparative performance.
	bool useCustomDrawing_ = true;

	// Use fast braking in FISTA.
	bool fastBrake_ = true;
	bool fastBrakeNu_ = false;

	// Number of threads to use in computation.
	int nthreads_ = 1;

	// Stop mode. StopEpsNorm2Rel is same as in TFOCS.
	Fista::Run::Stopping stopping_ = Fista::Run::StopEpsNorm2Rel;

	// Returned count of steps taken until the precision is achieved, or
	// (maxSteps_+1) if it never got achieved and the computation stopped by
	// the step limit. If the error message is present, this value can be any.
	int stepsTaken_ = -1;

	// If not empty, returns an error message on invalid arguments.
	std::string error_;

	// If not empty, contains the verbose log.
	std::string log_;
};

// Knows how to draw a Seudo object. Templatized to inline the pixel drawing
// function and avoid doing a virtual call for every pixel. This about doubles
// the performance compared to a plain Drawing.
//
// TODO: This whole templatization is rather convoluted at the moment
// and can use some rethinking to straighten it, but it works.
template <class SubDrawable,
	void drawPixelFn(SubDrawable &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)>
class SeudoCustomDrawing : public Fista::CustomDrawing<SubDrawable>
{
public:
	typedef Fista::CustomDrawing<SubDrawable> Parent;

	// @param seudo - the parent to draw
	// @param blobSpacing - spacing, to which the weights are adapted.
	//   The parent shrinks the weights as instructed by its spacing at
	//   the start of the compute() and then restores back, while in draw()
	//   it leaves them as-is, so this is the way to tell, what the parent
	//   has done and expects.
	SeudoCustomDrawing(const Seudo &seudo, int blobSpacing)
		: Fista::CustomDrawing<SubDrawable>(/*parallel*/ true),
		seudo_(seudo), blobSpacing_(blobSpacing)
	{ }

	// from CustomDrawing
	virtual void drawCustom(const Fista::Vector &input, SubDrawable &dest, Drawing::Limits &limits) {
		// First go the neuron ROIs.
		int end_roi = seudo_.rois_.size();
		if (limits.end_input_ < end_roi)
			end_roi = limits.end_input_;

		int i;
		for (i = limits.start_input_; i < end_roi; i++) {
			Sprite::drawWith<SubDrawable, drawPixelFn>(
				*seudo_.rois_[i], input, /*from_idx*/i, dest, /*at_x*/ 0, /*at_y*/ 0,
				limits.start_dest_y_, limits.end_dest_y_, 1.);
		}

		// Then go the blobs. Continue from the last i.
		int shrunk_wd = seudo_.getShrunkWd(blobSpacing_);
		for (; i < limits.end_input_; i++) {
			int didx = i - end_roi;
			int at_x = blobSpacing_ * (didx % shrunk_wd);
			int at_y = blobSpacing_ * (didx / shrunk_wd);

			Sprite::drawWith<SubDrawable, drawPixelFn>(
				*seudo_.blob_, input, /*from_idx*/i, dest, at_x, at_y,
				limits.start_dest_y_, limits.end_dest_y_, 1.);
		}
	}

protected:
	const Seudo &seudo_;
	int blobSpacing_;
};

// Also keep the plainest version of the drawing for performance comparison.
// Even "plaining" the custom version is noticeably faster than this.
#define FISTA_SEUDO_PLAINEST 1
#if FISTA_SEUDO_PLAINEST // {
	// Knows how to draw a Seudo object.
	class SeudoDrawing : public Fista::Drawing
	{
	public:
		// @param seudo - the parent to draw
		// @param blobSpacing - spacing, to which the weights are adapted.
		//   The parent shrinks the weights as instructed by its spacing at
		//   the start of the compute() and then restores back, while in draw()
		//   it leaves them as-is, so this is the way to tell, what the parent
		//   has done and expects.
		SeudoDrawing(const Seudo &seudo, int blobSpacing)
			: Fista::Drawing(/*parallel*/ true),
			seudo_(seudo), blobSpacing_(blobSpacing)
		{ }

		// from Drawing
		virtual void draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits);

	protected:
		const Seudo &seudo_;
		int blobSpacing_;
	};
#else // } {
	// A plain version used for performance comparison.
	typedef SeudoCustomDrawing<Fista::Drawable, Fista::drawPixelBasic> SeudoDrawing;
#endif // }

}; // namespace Fista

#endif // __SEUDO_HPP__


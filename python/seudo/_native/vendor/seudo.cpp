#include <math.h>
#include "fista.hpp"
#include "fista_types.hpp"
#include "fista_gradient.hpp"
#include "seudo.hpp"
#include "draw.hpp"
#include "strprintf.hpp"

namespace Fista {

void Seudo::compute()
{
	// Validate the inputs.
	error_.clear();
	if (image_.img_.size() == 0) {
		error_ = "Empty image_.";
		return;
	}
	if (image_.img_.size() != image_.wd_ * image_.ht_) {
		error_ = "Image vector size != width * height.";
		return;
	}
	if (!blob_) {
		error_ = "Blob sprite is NULL.";
		return;
	}
	if (blob_->img_.size() != blob_->wd_ * blob_->ht_) {
		error_ = "Blob vector size != width * height.";
		return;
	}
	for (auto it = rois_.begin(); it != rois_.end(); ++it) {
		std::shared_ptr<Sprite> &s = *it;
		// Not very good diagnostics but all that is easy to do without formatting.
		if (!s) {
			error_ = "A ROIs sprite is NULL.";
			return;
		}
		if (s->img_.size() != s->wd_ * s->ht_) {
			error_ = "A ROIs sprite vector size != width * height.";
			return;
		}
	}
	if (weights_.size() != rois_.size() + image_.img_.size()) {
		error_ = "Weights size is incorrect.";
		return;
	}
	if (lambda_.size() != weights_.size()) {
		error_ = "Lambda size is incorrect.";
		return;
	}

	if (verbose_) {
		log_ += "Image to match:\n";
		log_ += formatImg(image_, 999.);
	}

	// Optimize the sprites
	blob_->cropWhitespace();
	if (verbose_) {
		log_ += strprintf("Blob after crop: wd=%d ht=%d origin_x=%d origin_y=%d pixels=%d\n",
			blob_->wd_, blob_->ht_, blob_->origin_x_, blob_->origin_y_, (int)blob_->img_.size());
		log_ += formatImg(*blob_, 999.);
	}
	for (int i = 0; i < rois_.size(); i++) {
		std::shared_ptr<Sprite> &s = rois_[i];
		s->cropWhitespace();
		
		if (verbose_) {
			log_ += strprintf("ROI %d after crop: wd=%d ht=%d origin_x=%d origin_y=%d pixels=%d\n",
				i, s->wd_, s->ht_, s->origin_x_, s->origin_y_, (int)s->img_.size());
			log_ += formatImg(*s, 999.);
		}
	}

	// Shrink the weights and lambda.
	Fista::Vector saved_weights;
	saved_weights.swap(weights_);
	shrinkToBlobSpacing(saved_weights, weights_);

	Fista::Vector saved_lambda;
	saved_lambda.swap(lambda_);
	shrinkToBlobSpacing(saved_lambda, lambda_);

	// Wipe out the weights for the blobs not on the used grid.
	for (int i = 0; i < weights_.size(); i++) {
		int at_x = i % image_.wd_;
		int at_y = i / image_.wd_;
		if (at_x % blobSpacing_ != 0 || at_y % blobSpacing_ != 0)
			weights_[i] = 0;
	}

	std::shared_ptr<Fista::ScaledGradient> gradient;
	if (useCustomDrawing_) {
		// Drawing with a custom per-pass code.
		auto drawingPass1 = std::make_shared< SeudoCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass1> >(*this, blobSpacing_);
		auto drawingPass2 = std::make_shared< SeudoCustomDrawing<Fista::BaseDrawScaledGradient, &Fista::drawPixelPass2> >(*this, blobSpacing_);

		if (multiGrad_) {
			gradient = std::make_shared<Fista::MultiPosDrawScaledGradient>(
				image_.wd_, image_.ht_, drawingPass1, drawingPass2, /*b*/ image_.img_, lambda_, nthreads_);
		} else {
			gradient = std::make_shared<Fista::PosDrawScaledGradient>(
				image_.wd_, image_.ht_, drawingPass1, drawingPass2, /*b*/ image_.img_, lambda_, 1, nthreads_);
		}
	} else {
		// Basic drawing.
		auto drawing = std::make_shared<SeudoDrawing>(*this, blobSpacing_);
		if (multiGrad_) {
			gradient = std::make_shared<Fista::MultiPosDrawScaledGradient>(
				image_.wd_, image_.ht_, drawing, /*b*/ image_.img_, lambda_, nthreads_);
		} else {
			gradient = std::make_shared<Fista::PosDrawScaledGradient>(
				image_.wd_, image_.ht_, drawing, /*b*/ image_.img_, lambda_, 1, nthreads_);
		}
	}

	Fista::Run run(gradient, std::make_shared<Fista::PosLimiter>(), weights_, /*diffEps*/ eps_);
	run.stopping_ = stopping_;
	run.fastBrake_ = fastBrake_;
	run.fastBrakeNu_ = fastBrakeNu_;

	stepsTaken_ = run.repeatedSteps(maxSteps_);
	weights_.swap(run.x_); // return the new weights

	// Restore the lambda to an unshrunk form.
	saved_lambda.swap(lambda_); // Lambda is a constant, to no unshrinking.
	saved_weights.swap(weights_);
	unshrinkToBlobSpacing(saved_weights, weights_);

	if (verbose_) {
		Fista::DrawableImg img(image_.wd_, image_.ht_);
		log_ += strprintf("Stopped after %d steps\n", stepsTaken_);
		double dist = distance();
		log_ += strprintf("Distance: %f, avg %f\n", dist, dist / sqrt(image_.wd_ * image_.ht_));

		for (int i = 0; i < rois_.size(); i++) {
			log_ += strprintf("%d: %f  ", i, weights_[i]);
			rois_[i]->draw(weights_, i, img, 0, 0);
		}
		log_ += strprintf("\n");

		log_ += strprintf("Image of recognized ROIs:\n");
		log_ += formatImg(img, 999.);

		log_ += strprintf("Recognized blobs:\n");
		img.img_.assign(weights_.begin() + rois_.size(), weights_.end()),
		log_ += formatImg(img, 999.);
	}
}

void Seudo::shrinkToBlobSpacing(const Vector& src, Vector &dst)
{
	int nrois = rois_.size();
	int shrunk_wd = getShrunkWd(blobSpacing_);
	int shrunk_ht = getShrunkHt(blobSpacing_);

	dst.resize(nrois + shrunk_wd * shrunk_ht);

	// The ROIs are copie dunchanged.
	for (int i = 0; i < nrois; i++) {
		dst[i] = src[i];
	}
	for (int y = 0; y < shrunk_ht; y++) {
		for (int x = 0; x < shrunk_wd; x++) {
			dst[nrois + y*shrunk_wd + x] = src[nrois + y*blobSpacing_*image_.wd_ + x*blobSpacing_];
		}
	}
}

void Seudo::unshrinkToBlobSpacing(const Vector& src, Vector &dst)
{
	int nrois = rois_.size();
	int shrunk_wd = getShrunkWd(blobSpacing_);
	int shrunk_ht = getShrunkHt(blobSpacing_);

	// The ROIs are copie dunchanged.
	for (int i = 0; i < nrois; i++) {
		dst[i] = src[i];
	}
	for (int y = 0; y < image_.ht_; y++) {
		for (int x = 0; x < image_.wd_; x++) {
			if (y % blobSpacing_ == 0 && x % blobSpacing_ == 0) {
				dst[nrois + y*image_.wd_ + x] = src[nrois + (y / blobSpacing_)*shrunk_wd + (x / blobSpacing_)];
			} else {
				dst[nrois + y*image_.wd_ + x] = 0.;
			}
		}
	}
}

Fista::DrawableImg Seudo::draw() const
{
	Fista::DrawableImg target(image_.wd_, image_.ht_);
	auto drawing = std::make_shared<SeudoDrawing>(*this, /*blobSpacing*/ 1);
	drawing->drawSimple(weights_, target);
	return target;
}

double Seudo::distance() const
{
	Fista::DrawableImg current = draw();
	return Fista::distance(image_.img_, current.img_);
}

#if FISTA_SEUDO_PLAINEST // {
void SeudoDrawing::draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits)
{
	// First go the neuron ROIs.
	int end_roi = seudo_.rois_.size();
	if (limits.end_input_ < end_roi)
		end_roi = limits.end_input_;

	int i;
	for (i = limits.start_input_; i < end_roi; i++) {
		seudo_.rois_[i]->draw(input, /*from_idx*/i, dest, /*at_x*/ 0, /*at_y*/ 0,
			limits.start_dest_y_, limits.end_dest_y_);
	}

	// Then go the blobs. Continue from the last i.
	int shrunk_wd = seudo_.getShrunkWd(blobSpacing_);
	for (; i < limits.end_input_; i++) {
		int didx = i - end_roi;
		int at_x = blobSpacing_ * (didx % shrunk_wd);
		int at_y = blobSpacing_ * (didx / shrunk_wd);

		seudo_.blob_->draw(input, /*from_idx*/i, dest, at_x, at_y,
			limits.start_dest_y_, limits.end_dest_y_);
	}
}
#endif // }

}; // namespace Fista

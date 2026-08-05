// FISTA helper types.

#include <math.h>
#include "fista_types.hpp"

namespace Fista {

// ------------------- Vector --------------------------------

double distance(const Vector &v1, const Vector &v2)
{
	if (v1.size() != v2.size())
		return -1.;

	double dist = 0.;
	for (size_t i = 0; i < v1.size(); i++) {
		double v = v1[i] - v2[i];
		dist += v * v;
	}

	return sqrt(dist);
}

// ====================== Drawing =======================================

// ------------------- Drawable ------------------------------

Drawable::Drawable(int wd, int ht)
	: wd_(wd), ht_(ht)
{ }

Drawable::~Drawable()
{ }

// ------------------- DrawableImg ---------------------------

DrawableImg::DrawableImg(int wd, int ht)
	: Drawable(wd, ht), img_(wd*ht, 0.)
{ }

void DrawableImg::drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
{
	img_[to_y * wd_ + to_x] += k * x[from_idx];
}

void DrawableImg::clear()
{
	img_.assign(wd_*ht_, 0.);
}

// ------------------- Drawing -------------------------------

Drawing::~Drawing()
{ }

void Drawing::drawSimple(const Fista::Vector &input, Fista::Drawable &dest)
{
	Limits lim;
	lim.start_input_ = 0;
	lim.end_input_ = input.size();
	lim.start_dest_y_ = 0;
	lim.end_dest_y_ = dest.ht_;

	draw(input, dest, lim);
}

// ------------------- DrawPipeline --------------------------

void DrawPipeline::draw(const Fista::Vector &x, Fista::Drawable &dest, Limits &limits)
{
	if (filters_.empty()) {
		head_->draw(x, dest, limits);
	} else {
		filters_.back()->setDest(&dest);
		head_->draw(x, *filters_.front(), limits);
	}
}

// ------------------- Sprite --------------------------------

void Sprite::draw(const Fista::Vector &input, int from_idx, Drawable &dest,
	int at_x, int at_y, int start_y, int end_y, double k) const
{
	drawWith<Drawable, drawPixelBasic>(*this, input, from_idx, dest,
		at_x, at_y, start_y, end_y, k);
}

void Sprite::cropWhitespace()
{
	// Start by finding the bounding box
	int min_x = wd_;
	int max_x = -1;
	int min_y = ht_;
	int max_y = -1;

	for (int y = 0; y < ht_; y++) {
		for (int x = 0; x < wd_; x++) {
			if (img_[y*wd_ + x] != 0.) {
				if (x < min_x)
					min_x = x;
				if (x > max_x)
					max_x = x;
				if (y < min_y)
					min_y = y;
				if (y > max_y)
					max_y = y;
			}
		}
	}

	// Move up from max used value to the boundary.
	max_x++;
	max_y++;

	// printf("min_x=%d max_x=%d min_y=%d max_y=%d\n", min_x, max_x, min_y, max_y);
	if (min_x == 0 && max_x == wd_ && min_y == 0 && max_y == ht_)
		return; // nothing to do

	int newwd = max_x - min_x;
	if (min_x >= max_x || min_y >= max_y) {
		// The sprite became empty.
		min_x = max_x = min_y = max_y = 0;
		newwd = 0;
	} else {
		// Move the img_ values in place. 
		for (int y = min_y; y < max_y; y++) {
			for (int x = min_x; x < max_x; x++) {
				// printf("Write at (%d, %d) = %d from (%d, %d) = %d, value %f\n",
					// (x - min_x), (y - min_y), (y - min_y)*newwd + (x - min_x),
					// x, y, y*wd_ + x, img_[y*wd_ + x]);
				img_[(y - min_y)*newwd + (x - min_x)] = img_[y*wd_ + x];
			}
		}
	}
	wd_ = newwd;
	ht_ = max_y - min_y;
	img_.resize(wd_*ht_);
	origin_x_ -= min_x;
	origin_y_ -= min_y;
}

}; // namespace Fista

#ifndef __FISTA_TYPES_HPP__
#define __FISTA_TYPES_HPP__

#include <vector>
#include <memory>

// Types used in the FISTA implementation

namespace Fista {

// ---------------------- Helper types ----------------------------------

// A vector of values.
typedef std::vector<double> Vector;

// Find the norm2 distance between two vectors. The sizes must be the same,
// returns -1 if the sizes are different.
double distance(const Vector &v1, const Vector &v2);

// ---------------------- Drawing ---------------------------------------
// Support for specifying the sparse dependency matrices by "drawing" the
// target vector from the values of the source vector. This get used to
// compute the gradient on the fly, without storing the whole dependency
// matrix, storing only the source and gradient. Or it can be used to draw
// the actual output image.

// Representation of the target, that might be either a straightforward image
// or a gradient computation (such as PosMatSquareMinimizer).
class Drawable {
public:
	// The width and height don't have to mean anything, the actual image might
	// be unidimensional or multi-dimensional. But since the typical images are
	// 2-dimensional, keeping the width and height here is convenient, and any
	// other number of dimensions can be mapped into them. The width has an
	// additional meaning: when the drawing is parallelized by rows, it is
	// split by whole rows, i.e. in the width units.
	//
	// @param wd - width of image
	// @param ht - height of image
	Drawable(int wd, int ht);

	virtual ~Drawable();

	// "Draw" a "pixel" to the destination by adding a scaled value of the source
	// "pixel" to the destination one.
	//
	// @param x - input values
	// @param from_idx - index of value in x that gets propagated
	// @param to_x - "X coordinate" of the drawable cell where the propagated value gets added
	// @param to_y - "Y coordinate" of the drawable cell where the propagated value gets added
	//   (internally the typical way is to compute the single (to_idx = to_y * wd_ + to_x).
	// @param k - coefficient to multiply the value from x by, before adding
	virtual void drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k) = 0;

public:
	// width of image
	int wd_ = 0;
	// height of image
	int ht_ = 0;
};

// A common way to store a plain image instead of generating the gradient.
class DrawableImg : public Drawable
{
public:
	// @param wd - width of image
	// @param ht - height of image
	DrawableImg(int wd, int ht);

	// from Drawable
	virtual void drawPixel(const Fista::Vector &x, int from_idx, int to_x, int to_y, double k);

	// Clear the image back to 0s, for a new drawing.
	void clear();

	// Image of size (wd_*ht_)
	Fista::Vector img_;
};

// An encapsulation of drawing logic that can draw to some Drawable.
class Drawing
{
public:
	// Limits to draw within. Honored only by the parallel-capable subclasses.
	// The input is limited with individual granularity, the destination
	// is limited in units of rows, AKA Y-coordinates: each row contains dest.wd_ pixels.
	struct Limits
	{
		// draw only the portion originating starting from this input
		int start_input_;
		// draw only the portion ending before this input, i.e. [start_input, end_input)
		int end_input_;
		// draw only the portion of the destination starting from this row;
		// must be >= 0
		int start_dest_y_;
		// draw only the portion of the destination ending before this row
		// i.e. [start_dest_y, end_dest_y); must be <= dest.ht_
		int end_dest_y_;
	};

	// @param parallel - the subclass supports the parallel processing.
	Drawing(bool parallel = false)
		: parallel_(parallel)
	{ }

	virtual ~Drawing();

	// The subclass must define it to do the "drawing", by calling drawPixel()
	// for every dependency between two "pixels".
	// @param input - the X source vector
	// @param dest - drawable destination
	// @param limits - the limits for drawing; they are intended for drawing
	//    in parallel by multiple threads, so the subclasses that are not
	//    parallel-capable, are not expected to honor them
	virtual void draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits) = 0;

	// The simple version that auto-sets the limits for the whole input
	// and whole destination.
	void drawSimple(const Fista::Vector &input, Fista::Drawable &dest);

	// Returns true if subclass supports the parallel processing.
	bool isParallel() const
	{
		return parallel_;
	}

public:
	// The subclass supports the parallel processing.
	bool parallel_;
};

// A custom version of Drawing optimized for a specific target.
// 
// TODO: This whole templatization is rather convoluted at the moment
// and can use some rethinking to straighten it, but it works.
//
// @param SubDrawable - a subclass of Drawable that will provide the
//   custom functions to be inlined into specializations of CustomDrawing,
//   to improve the performance.
template <class SubDrawable>
class CustomDrawing : public Drawing
{
public:
	// @param parallel - the subclass supports the parallel processing,
	//    typically if we bother to go to a CustomDrawing, should implement
	//    parllelism too
	CustomDrawing(bool parallel)
		: Drawing(parallel)
	{ }

	// From Drawing, a legacy version that will crash on wrong drawable.
	virtual void draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits)
	{
		drawCustom(input, dynamic_cast<SubDrawable &>(dest), limits);
	}

	virtual void drawCustom(const Fista::Vector &input, SubDrawable &dest, Limits &limits) = 0;
};

// A helper that allows to convert a custom drawing to one that uses
// a common Drawable.
//
// A custom drawing is commonly defined as:
//
// template <class SubDrawable,
//     void drawPixelFn(SubDrawable &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)>
// class SomeCustomDrawing : public Fista::CustomDrawing<SubDrawable>
// {
// public:
//     ...
// 
//     // from CustomDrawing
//     virtual void drawCustom(const Fista::Vector &input, SubDrawable &dest, Drawing::Limits &limits) {
//         ...
//         drawPixelFn(dest, input, from_idx, dest_x, dest_y, k); // for each pixel
//         ...
//     }
// };
//
// It can be converted to a common drawing as
//
// typedef SomeCustomDrawing<Fista::Drawable, Fista::drawPixelBasic> SomeDrawing;
//
// This converted drawing is usually still faster than a completely plain
// drawing, because more of its internal structure gets inlined and optimized.
inline void 
drawPixelBasic(Fista::Drawable &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)
{
	dest.drawPixel(x, from_idx, to_x, to_y, k);
}

// This is an encapsulation of both drawing and drawable that knows how
// to draw into itself. A gradient computation class is an example of such
// a self-drawable: it holds in itself the information about the specific 
// function as a Drawing, and it draws into itself as a Drawable to compute
// the gradient. This encapsulation is used to parallelize the drawing process,
// by calling the self-drawing from multiple threads with differnt limits.
class SelfDrawable
{
public:
	virtual void selfDraw(const Fista::Vector &input, Drawing::Limits &limits) = 0;
};

// A "filter" that allows to compose the transformations in drawing,
// such as draw a text and blur it.
class DrawFilter : public Drawable
{
public:
#if 0
	// The dimendions will be copied from destination in setDest().
	DrawFilter()
		: Drawable(0, 0)
	{ }
#endif

	// This constructor allows to set the image size early, before
	// setting the destination. Setting the destination will override
	// the dimenstions, but normally they should be the same.
	//
	// @param wd - width of image
	// @param ht - height of image
	DrawFilter(int wd, int ht)
		: Drawable(wd, ht)
	{ }

	// Set the destination where drawPixel() will translate the
	// drawing.
	// @param dest - destination that the following calls of drawPixel()
	//   will use for the nested calls of drawPixel(). The caller must
	//   make sure that the destination won't get destroyed while it's used.
	void setDest(Fista::Drawable *dest)
	{
		dest_ = dest;
		wd_ = dest->wd_;
		ht_ = dest->ht_;
	}

protected:
	// Consider valid only for the duration of a call of drawPixel().
	Fista::Drawable *dest_;
};

// A drawing that consists of a pipeline of filters.
class DrawPipeline : public Drawing
{
public:
	// Set the head of the pipeline.
	void setHead(std::shared_ptr<Fista::Drawing> head)
	{
		head_ = head;
	}

	// Add a filter to the tail of the pipeline.
	void addFilter(std::shared_ptr<DrawFilter> tail)
	{
		if (!filters_.empty()) {
			filters_.back()->setDest(tail.get());
		}
		filters_.emplace_back(tail);
	}

	// from Drawing
	virtual void draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits);

protected:
	// Head of the pipeline is a drawing.
	std::shared_ptr<Fista::Drawing> head_;
	// The rest of a pipeline are filters, going in sequence.
	std::vector<std::shared_ptr<DrawFilter>> filters_;
};

// A "rectangular blob" containing the "opaqueness values" per pixel,
// that can be drawn as a representation of a single source value.
// It inherits from DrawableImg to make the sprite creation easier.
class Sprite : public DrawableImg
{
public:
	// Creates a sprite image of all-0s.
	// The arguments have the same meaning as the fields.
	Sprite(int wd, int ht, int origin_x = 0, int origin_y = 0)
		: DrawableImg(wd, ht), origin_x_(origin_x), origin_y_(origin_y)
	{ }

	// "Draw" the sprite to the destination by adding the "sprite pixels" scaled
	// by the source value to the destination "pixels".
	//
	// This is really a convenience wrapper over specialization
	//   drawWith<Drawable, drawPixelBasic>(*this, ...)
	// But using the specialization directly is faster. It couls also be inlined
	// right here to avoid the penalty, but it kept as-is for comparison purposes
	// (and if you're not using the customized template, you probably don't care
	// much about a small performance penalty).
	//
	// @param input - input values, don't care here, passed through to Drawable dest
	// @param from_idx - index of value in x that gets propagated, don't care here,
	//   passed through to the Drawable dest
	// @param dest - the destination drawable
	// @param at_x - the X coordinate in the destination where the origin of this
	//   sprite will be placed
	// @param at_y - the Y coordinate in the destination where the origin of this
	//   sprite will be placed
	// @param start_y - draw only the portion of the destination starting from
	//   the row start_y; must be >= 0
	// @param end_y - draw only the portion of the destination ending before
	//   the row end_y, i.e. [start_y, end_y); must be <= dest.ht_
	// @param k - the values of "pixels" from the sprite get multiplied by this
	//   coefficient before passing them to the Drawable dest (which in turn will then
	//   also multiply them by the value x[from_idx]). This allows to adjust the
	//   "opacity" of the same single sprite in different uses.
	void draw(const Fista::Vector &input, int from_idx, Drawable &dest,
		int at_x, int at_y, int start_y, int end_y, double k = 1.) const;

	// "Draw" the sprite to the destination by adding the "sprite pixels" scaled
	// by the source value to the destination "pixels".
	//
	// @param input - input values, don't care here, passed through to Drawable dest
	// @param from_idx - index of value in x that gets propagated, don't care here,
	//   passed through to the Drawable dest
	// @param dest - the destination drawable
	// @param at_x - the X coordinate in the destination where the origin of this
	//   sprite will be placed
	// @param at_y - the Y coordinate in the destination where the origin of this
	//   sprite will be placed
	// @param k - the values of "pixels" from the sprite get multiplied by this
	//   coefficient before passing them to the Drawable dest (which in turn will then
	//   also multiply them by the value x[from_idx]). This allows to adjust the
	//   "opacity" of the same single sprite in different uses.
	void draw(const Fista::Vector &input, int from_idx, Drawable &dest,
		int at_x, int at_y, double k = 1.) const
	{
		draw(input, from_idx, dest, at_x, at_y, 0, dest.ht_, k);
	}

	// Same drawing, done as a template to hardcode the pixel-drawing function.
	template <class SubDrawable,
		void drawPixelFn(SubDrawable &dest, const Fista::Vector &x, int from_idx, int to_x, int to_y, double k)>
	static void drawWith(const Sprite &sprite, const Fista::Vector &input, int from_idx, SubDrawable &dest,
		int at_x, int at_y, int start_y, int end_y, double k = 1.)
	{
		int y = 0;
		int dest_y = at_y - sprite.origin_y_;
		if (dest_y < start_y) {
			y = start_y - dest_y;
			dest_y = start_y;
		}

		int limit_y = (end_y - dest_y) + y;
		if (limit_y > sprite.ht_) {
			limit_y = sprite.ht_;
		}
		if (y >= limit_y) {
			return; // sprite and destination don't intersect by y at all
		}

		int start_x = 0;
		int start_dest_x = at_x - sprite.origin_x_;
		if (start_dest_x < 0) {
			start_x = - start_dest_x;
			start_dest_x = 0;
		}

		int limit_x = (dest.wd_ - start_dest_x) + start_x;
		if (limit_x > sprite.wd_) {
			limit_x = sprite.wd_;
		}
		if (start_x >= limit_x) {
			return; // sprite and destination don't intersect by x at all
		}

		int offset_x = start_dest_x - start_x;

		for (; y < limit_y; y++, dest_y++) {
			int rowpos = y * sprite.wd_;
			for (int x = start_x; x < limit_x; x++) {
				double v = sprite.img_[rowpos + x];
				if (v != 0.) {
					drawPixelFn(dest, input, from_idx, /*dest_x*/ x + offset_x, dest_y, k * v);
				}
			}
		}
	}

	// Find if there is whitespace around this image and crop it, reducing the
	// width and height (if possible), resizing and repopulating img_ accordingly,
	// and keeping the origin logically the same by adjusting it for the cropping.
	// This makes the subsequent drawing more efficient, since it has to cover
	// less area.
	void cropWhitespace();

public:
	// The origin point relatively to the upper left corner of the image in
	// img_. A positive value typically can be used to represent an origin
	// at a center of the image, so for example an image of size (3 * 3)
	// with origin (1, 1) will have its origin at the center. A negative
	// value typically can be used to represent a cropped image, where the
	// original origin was at the upper left corner but after we've cropped
	// out the zero parts, that origin is to the left and up of the new
	// smaller image.
	int origin_x_;
	int origin_y_;
};

}; // namespace Fista

#endif // __FISTA_TYPES_HPP__

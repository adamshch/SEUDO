#ifndef __DRAW_HPP__
#define __DRAW_HPP__

#include <string>
#include "fista_types.hpp"

// Helper functions to draw text and such in a matrix, to create interesting
// images for blurring and such.

// Drawing that blurs the original image.
class BlurDrawing : public Fista::Drawing
{
public:
	// @param wd - width of the image
	// @param ht - height of the image
	// @param bradius - blur radius (i.e. 0 will average each pixel with itself
	//      and cause no blur, 1 will include up to 9 pixels, going +-1 by width
	//      and heihgt from the target pixel, unless it's at the edge of the image,
	//      and so on)
	BlurDrawing(int wd, int ht, int bradius)
		: wd_(wd), ht_(ht), bradius_(bradius)
	{ }

	// from Drawing
	virtual void draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits);

protected:
	// Image width.
	int wd_;
	// Image height.
	int ht_;
	// Blur radius.
	int bradius_;
};

// This blurring works slightly differently from BlurDrawing, by using the
// same coefficients at the corners of the image.
class BlurFilter : public Fista::DrawFilter
{
public:
	// @param wd - width of the image
	// @param ht - height of the image
	// @param bradius - blur radius (i.e. 0 will average each pixel with itself
	//      and cause no blur, 1 will include up to 9 pixels, going +-1 by width
	//      and heihgt from the target pixel, unless it's at the edge of the image,
	//      and so on)
	BlurFilter(int wd, int ht, int bradius)
		: DrawFilter(wd, ht), bradius_(bradius),
		reduction_(1./((bradius+1) * (bradius+1)))
	{ }

	// from Drawable
	virtual void drawPixel(const Fista::Vector &input, int from_idx, int to_x, int to_y, double k);

protected:
	// Blur radius.
	int bradius_;
	// When the blurred pixel gets distirbuted over many, this is the
	// reduction in "brightness" on each destination pipxel.
	double reduction_;
};

// Draw a text. The characters that are not present in the font as bitmaps
// are left as spaces. This is a convenience and historic function for
// simple drawing. It draws by writing a fixed value into each pixel if
// that pixel contained a lesser value.
//
// @param wd - width of the canvas to draw on
// @param ht - height of the canvas to draw on
// @param canvas - canvas to draw on (must have the size of wd*ht)
// @param charwd - width of each character in font bitmap
// @param charht - height of each character in font bitmap
// @param bitmaps - font bitmaps of characters; each entry in bitmap contains:
//    * 1 byte of the character code;
//    * (charwd*charht) bytes contatining the bitmap of the character row by row,
//      in each pixel '.' means an empty pixel, anything elase a drawable pixel;
//    * \0 at the end;
//    the bitmap array ends with a NULL
// @param color - the "color" to draw with, the pixels that are already "darker"
//    will be left unchanged
// @param text - string of text to draw
// @param startx - starting X position of top left corner (in pixels)
// @param startx - starting Y position of top left corner (in pixels)
void drawTextSimple(int wd, int ht, Fista::Vector &canvas,
	int charwd, int charht, const char **bitmaps,
	double color, const char *text, int startx = 0, int starty = 0);

// Draw a character into a Drawable. The drawing gets done by adding the
// character with coefficient/"weight"/"color" k to the existing contents of
// the canvas (this is different from drawTextSimple() that just writes a fixed
// value into each pixel).
//
// @param wd - width of the canvas to draw on
// @param ht - height of the canvas to draw on
// @param canvas - canvas to draw on (must have the size of wd*ht)
// @param charwd - width of each character in font bitmap
// @param charht - height of each character in font bitmap
// @param bitmaps - font bitmaps of characters; each entry in bitmap contains:
//    * 1 byte of the character code;
//    * (charwd*charht) bytes contatining the bitmap of the character row by row,
//      in each pixel '.' means an empty pixel, anything elase a drawable pixel;
//    * \0 at the end;
//    the bitmap array ends with a NULL
// @param x - X of the upper-left corner of the character on canvas
// @param y - Y of the upper-left corner of the character on canvas
// @param c - code of the character to draw
// @param input - gets passed through to Drawable, the source that controls
//    the drawing
// @param from_idx - gets passed through to Drawable, index in input that controls
//    the drawing of this character and the weight with which it's added
// @param k - gets passed through to Drawable, a coefficient by which the weight
//    from from_idx gets multiplied before drawing (1 is the typical case)
void drawChar(int wd, int ht, Fista::Drawable &canvas,
	int charwd, int charht, const char **bitmaps,
	int x, int y, char c,
	const Fista::Vector &input, int from_idx, double k = 1.);

// Draw every character in the font overlapping over the same position.
// Most of the arguments, except n_chars, are the same as in drawChar().
//
// @param wd - width of the canvas to draw on
// @param ht - height of the canvas to draw on
// @param canvas - canvas to draw on (must have the size of wd*ht)
// @param charwd - width of each character in font bitmap
// @param charht - height of each character in font bitmap
// @param bitmaps - font bitmaps of characters; each entry in bitmap contains:
//    * 1 byte of the character code;
//    * (charwd*charht) bytes contatining the bitmap of the character row by row,
//      in each pixel '.' means an empty pixel, anything elase a drawable pixel;
//    * \0 at the end;
//    the bitmap array ends with a NULL
// @param x - X of the upper-left corner of the character on canvas
// @param y - Y of the upper-left corner of the character on canvas
// @param n_chars - expected number of characters in the font, the function will
//    stop at the earliest of n_chars or end of font.
// @param input - gets passed through to Drawable, the source that controls
//    the drawing
// @param from_idx - gets passed through to Drawable, index in input that controls
//    the drawing and the weight with which it's added for the first character
//    in the font, the following elements (up to n_chars) determine the next
//    characters
// @param k - gets passed through to Drawable, a coefficient by which the weight
//    from from_idx gets multiplied before drawing (1 is the typical case)
void drawEveryChar(int wd, int ht, Fista::Drawable &canvas,
	int charwd, int charht, const char **bitmaps,
	int x, int y, int n_chars,
	const Fista::Vector &input, int from_idx, double k = 1.);

// Get the number of bitmaps in the font.
int getFontSize(const char **bitmaps);

// Print an image in text-bitmap for debugging.
// Scale allows to rescale the fractional values to where they look better
// as integers.
int
printImg(int wd, int ht, const Fista::Vector &img, double scale = 1.);

inline int printImg(const Fista::DrawableImg &img, double scale=1.)
{
	return printImg(img.wd_, img.ht_, img.img_, scale);
}

// Printing to a string.
std::string formatImg(int wd, int ht, const Fista::Vector &img, double scale = 1.);
inline std::string formatImg(const Fista::DrawableImg &img, double scale=1.)
{
	return formatImg(img.wd_, img.ht_, img.img_, scale);
}

// A limited font of 8x8 bitmaps.
extern const char *bitmaps8x8[];

#endif // __DRAW_HPP__

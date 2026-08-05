#include <ctype.h>
#include <stdio.h>
#include "draw.hpp"
#include "strprintf.hpp"

// -------------- BlurDrawing --------------------------------

void BlurDrawing::draw(const Fista::Vector &input, Fista::Drawable &dest, Limits &limits)
{
	for (int i = 0; i < wd_; i++) {
		for (int j = 0; j < ht_; j++) {
			int x0 = i - bradius_;
			if (x0 < 0) // crop on edges
				x0 = 0;
			int x1 = i + bradius_ + 1;
			if (x1 > wd_) // crop on edges
				x1 = wd_;

			int y0 = j - bradius_;
			if (y0 < 0) // crop on edges
				y0 = 0;
			int y1 = j + bradius_ + 1;
			if (y1 > ht_) // crop on edges
				y1 = ht_;

			double avg = 1. / ((x1 - x0) * (y1 - y0)); // averaging by num of pixels in a blurring region
			// double avg = 1. / ((bradius + 1) * (bradius + 1));
			
			for (int x = x0; x < x1; x++) {
				for (int y = y0; y < y1; y++) {
					dest.drawPixel(input, /*from_idx*/ y*wd_ + x, /*to_x*/ i, /*to_y*/ j, /*k*/ avg);
				}
			}
		}
	}
}

// -------------- BlurFilter ---------------------------------

void BlurFilter::drawPixel(const Fista::Vector &input, int from_idx, int to_x, int to_y, double k)
{
	int x0 = to_x - bradius_;
	if (x0 < 0) // crop on edges
		x0 = 0;
	int x1 = to_x + bradius_ + 1;
	if (x1 > wd_) // crop on edges
		x1 = wd_;

	int y0 = to_y - bradius_;
	if (y0 < 0) // crop on edges
		y0 = 0;
	int y1 = to_y + bradius_ + 1;
	if (y1 > ht_) // crop on edges
		y1 = ht_;

	for (int xi = x0; xi < x1; xi++) {
		for (int yi = y0; yi < y1; yi++) {
			// from_idx represents the ultimate source, so it passes through
			dest_->drawPixel(input, from_idx, xi, yi, k * reduction_);
		}
	}
}

// -------------- font drawing -------------------------------

void drawTextSimple(int wd, int ht, Fista::Vector &canvas,
	int charwd, int charht, const char **bitmaps,
	double color, const char *text, int startx, int starty)
{
	int x = startx;
	int y = starty;

	for (; *text != 0; text++) {
		char c = *text;
		if (c == '\n') {
			x = 0;
			y += charht;
			continue;
		}
		for (const char **mapptr = bitmaps; *mapptr != NULL; mapptr++) {
			const char *cmap = *mapptr;
			if (*cmap != c) {
				continue;
			}
			++cmap;
			for (int row = 0; row < charht; row++) {
				if (y + row >= ht)
					break;
				for (int col = 0; col < charwd; col++) {
					if (*cmap == 0)
						break;  // shouldn't happen
					bool draw = (*cmap++ != '.');
					if (!draw || x + col >= wd) {
						continue; // not break, to skip over pixels to the end of the row
					}
					if (canvas[(y + row)*wd + (x+col)] < color) {
						canvas[(y + row)*wd + (x+col)] = color;
					}
				}
			}
			break;
		}
		x += charwd;
		if (x >= wd) {
			x = 0;
			y += charht;
		}
		if (y >= ht)
			break;
	}
}

void drawChar(int wd, int ht, Fista::Drawable &canvas,
	int charwd, int charht, const char **bitmaps,
	int x, int y, char c,
	const Fista::Vector &input, int from_idx, double k)
{
	for (const char **mapptr = bitmaps; *mapptr != NULL; mapptr++) {
		const char *cmap = *mapptr;
		if (*cmap != c) {
			continue;
		}
		++cmap;
		for (int row = 0; row < charht; row++) {
			if (y + row >= ht)
				break;
			for (int col = 0; col < charwd; col++) {
				if (*cmap == 0)
					break;  // shouldn't happen
				bool draw = (*cmap++ != '.');
				if (!draw || x + col >= wd) {
					continue; // not break, to skip over pixels to the end of the row
				}
				canvas.drawPixel(input, from_idx, x+col, y + row, k);
			}
		}
		break;
	}
}

void drawEveryChar(int wd, int ht, Fista::Drawable &canvas,
	int charwd, int charht, const char **bitmaps,
	int x, int y, int n_chars,
	const Fista::Vector &input, int from_idx, double k)
{
	for (const char **mapptr = bitmaps; *mapptr != NULL && n_chars > 0; mapptr++, n_chars--, from_idx++) {
		const char *cmap = *mapptr;
		++cmap; // skip over the character's ASCII code
		for (int row = 0; row < charht; row++) {
			if (y + row >= ht)
				break;
			for (int col = 0; col < charwd; col++) {
				if (*cmap == 0)
					break;  // shouldn't happen
				bool draw = (*cmap++ != '.');
				if (!draw || x + col >= wd) {
					continue; // not break, to skip over pixels to the end of the row
				}
				canvas.drawPixel(input, from_idx, x+col, y + row, k);
			}
		}
	}
}

int getFontSize(const char **bitmaps)
{
	int n = 0;
	for (const char **mapptr = bitmaps; *mapptr != NULL; mapptr++) {
		n++;
	}
	return n;
}

int
printImg(int wd, int ht, const Fista::Vector &img, double scale)
{
	for (int i = 0; i < ht; i++) {
		printf("   ");
		for (int j = 0; j < wd; j++) {
			printf("%3.0f ", img[i*wd + j] * scale);
		}
		printf("\n");
	}
	printf("--\n");
	return 0;
}

std::string formatImg(int wd, int ht, const Fista::Vector &img, double scale)
{
	std::string result;

	for (int i = 0; i < ht; i++) {
		result += "   ";
		for (int j = 0; j < wd; j++) {
			result += strprintf("%3.0f ", img[i*wd + j] * scale);
		}
		result += "\n";
	}
	result += "--\n";
	return result;
}

// -------------- fonts --------------------------------------

const char *bitmaps8x8[] = {
	"H" 
	".*....*."
	".*....*."
	".*....*."
	".******."
	".*....*."
	".*....*."
	".*....*."
	"........",

	"h" 
	".*......"
	".*......"
	".*......"
	".*****.."
	".*....*."
	".*....*."
	".*....*."
	"........",

	"E" 
	".******."
	".*......"
	".*......"
	".*****.."
	".*......"
	".*......"
	".******."
	"........",

	"e" 
	"........"
	"........"
	"..****.."
	".*....*."
	".******."
	".*......"
	"..****.."
	"........",

	"L" 
	".*......"
	".*......"
	".*......"
	".*......"
	".*......"
	".*......"
	".******."
	"........",

	"l" 
	".*......"
	"..*....."
	"..*....."
	"..*....."
	"..*....."
	"..*..*.."
	"...**..."
	"........",

	"O" 
	"..****.."
	".*....*."
	".*....*."
	".*....*."
	".*....*."
	".*....*."
	"..****.."
	"........",

	"o" 
	"........"
	"........"
	"..****.."
	".*....*."
	".*....*."
	".*....*."
	"..****.."
	"........",

	/*
	// Space can be copied for a blank start of a new character.
	// But it should not be included into the font because it draws
	// nothing.
	" " 
	"........"
	"........"
	"........"
	"........"
	"........"
	"........"
	"........"
	"........",
	*/

	NULL
};


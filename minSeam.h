/*
 * Copyright, 2013, Aeron Buchanan
 *
 * This file is part of TexSynth, a digital inpainting resource.
 *
 * TexSynth is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TexSynth is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TexSynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <vector>

#include "CImg.h"

using cimg_library::CImg;

namespace TexSynth
{

template<class T>
std::vector<uint> findMinVerticalSeam(CImg<T> const & _img)
{
	CImg<double> costs(_img.width(), _img.height(), 1, 1, 0.f);
	CImg<uint> refs(_img.width(), _img.height(), 1, 1, 0);

	// top line
	for ( int x = 0, y = 0; x < _img.width(); ++x )
	{
		costs(x, y) = _img(x, y);
	}

	struct MinHelper
	{
		CImg<double> const * scores;
		MinHelper(CImg<double> const * _sc) : scores(_sc) {}
		uint xp(uint x, uint y) { return (*scores)(x + 1, y - 1) < (*scores)(x, y - 1) ? x + 1 : x; }
		uint xm(uint x, uint y) { return (*scores)(x - 1, y - 1) < (*scores)(x, y - 1) ? x - 1 : x; }
		uint x(uint x, uint y)
		{
			uint bestx = xm(x, y);
			bestx = xp(bestx, y);
			return bestx;
		}
	} best(&costs);

	// rest
	for ( int y = 1; y < _img.height(); ++y )
	{
		int x = 0;
		// left-most pixel
		uint bestx = best.xp(x, y);
		costs(x, y) = _img(x, y) + costs(bestx, y - 1);
		refs(x, y) = bestx;

		// middle pixels
		for ( ++x; x < _img.width() - 1; ++x )
		{
			bestx = best.x(x, y);
			costs(x, y) = _img(x, y) + costs(bestx, y - 1);
			refs(x, y) = bestx;
		}

		// right-most pixel
		bestx = best.xm(x, y);
		costs(x, y) = _img(x, y) + costs(bestx, y - 1);
		refs(x, y) = bestx;
	}

	// backtrack
	std::vector<uint> seam(costs.height());
	int y = costs.height() - 1;
	uint & bestx = seam[y];
	for ( int x = 1; x < costs.width(); ++x )
		if ( costs(x, y) < costs(bestx, y) ) bestx = x;

	for ( --y; y >= 0; --y ) seam[y] = refs(seam[y + 1], y + 1);

	// DEBUG
	for ( int y = 0; y < _img.height(); ++y )
	{
		printf("%u: ", seam[y]);

		for ( int x = 0; x < _img.width(); ++x )
			printf("%6.1f%s ", _img(x,y), uint(x) == seam[y] ? "*" : " ");

		printf("\n");
	}

	printf("   ");
	for ( int x = 0; x < _img.width(); ++x )
		printf(" %6.1f ", costs(x, costs.height() - 1));

	printf("\n");
	// DEBUG END

	return seam;
}



} // TexSynth

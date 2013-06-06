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
std::vector<uint> findMinVerticalSeam(CImg<T> const & _I)
{
	CImg<double> costs(_I.width(), _I.height(), 1, 1, 0.f);
	CImg<uint> refs(_I.width(), _I.height(), 1, 1, 0);

	// top line
	for ( uint x = 0, y = 0; x < _I.width(); ++x )
	{
		costs(x, y) = _I(x, y);
	}

	struct MinHelper
	{
		CImg<T> I;
		MinHelper(CImg<T> const & _I) : I(_I) {}
		uint xp(uint x, uint y) { return I(x + 1, y - 1) < I(x, y - 1) ? x + 1 : x; }
		uint xm(uint x, uint y) { return I(x - 1, y - 1) < I(x, y - 1) ? x - 1 : x; }
		uint x(uint x, uint y)
		{
			uint bestx = xm(I, x, y);
			bestx = xp(I, bestx, y);
			return bestx;
		}
	} best(_I);

	// rest
	for ( uint y = 1; y < _I.height(); ++y )
	{
		uint x = 0;
		// left-most pixel
		uint bestx = best.xp(x, y);
		costs(x, y) = _I(x, y) + costs(bestx, y - 1);
		refs(x, y) = bestx;

		// middle pixels
		for ( ++x; x < _I.width() - 1; ++x )
		{
			bestx = best.x(x, y);
			costs(x, y) = _I(x, y) + costs(bestx, y - 1);
			refs(x, y) = bestx;
		}

		// right-most pixel
		bestx = best.xm(x, y);
		costs(x, y) = _I(x, y) + costs(bestx, y - 1);
		refs(x, y) = bestx;
	}

	// backtrack
	std::vector<uint> seam(costs.height());
	int y = costs.height() - 1;
	uint & bestx = seam[y];
	for ( uint x = 1; x < costs.width(); ++x )
		if ( costs(x, y) < costs(bestx, y) ) bestx = x;

	for ( --y; y >= 0; --y ) seam[y] = refs(seam[y + 1], y + 1);

	// DEBUG
	for ( uint y = 0; y < _I.height(); ++y )
	{
		for ( uint x = 0; x < _I.width(); ++x )
			printf("%s%6d ", x == seam[y] ? "*" : " ", _I(x,y));

		printf("\n");
	}
	for ( uint x = 0; x < _I.width(); ++x )
		printf(" %6d ", costs(x, costs.height() - 1));
	// DEBUG END

	return seam;
}



} // TexSynth

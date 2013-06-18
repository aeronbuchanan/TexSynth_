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
#include <limits>

#include "CImg.h"
using cimg_library::CImg;

#include "table.h"
#include "patch.h"


namespace TexSynth
{

namespace Seams
{

enum Direction { Downwards = 1, Rightwards, Upwards, Leftwards };

template<Direction D>
struct Coord : public TexSynth::Coord
{
	typedef TexSynth::Coord::Type Type;

	Coord() : TexSynth::Coord(0, 0) {}
	Coord(Type _x, Type _y) : TexSynth::Coord(_x, _y) {}

	Coord & advNeg() { --orthComp(); return *this; }
	Coord & advPos() { ++orthComp(); return *this; }
	Coord & advPrev() { dirComp() -= dirStep(); return *this; }
	Coord & advNext() { dirComp() += dirStep(); return *this; }

	Coord & setOrth(Type _i) { orthComp() = _i; return *this; }
	Coord & setDir(Type _i) { dirComp() = _i; return *this; }

	Type & dirComp();
	Type & orthComp();

	int dirStep() const;
};

// Specializations
template<> TexSynth::Coord::Type & Coord<Downwards>::orthComp() { return x();  }
template<> TexSynth::Coord::Type & Coord<Downwards>::dirComp() { return y(); }
template<> int Coord<Downwards>::dirStep() const { return +1; }

template<> TexSynth::Coord::Type & Coord<Upwards>::orthComp() { return x();  }
template<> TexSynth::Coord::Type & Coord<Upwards>::dirComp() { return y(); }
template<> int Coord<Upwards>::dirStep() const { return -1; }

template<> TexSynth::Coord::Type & Coord<Rightwards>::orthComp() { return y(); }
template<> TexSynth::Coord::Type & Coord<Rightwards>::dirComp() { return x(); }
template<> int Coord<Rightwards>::dirStep() const { return +1; }

template<> TexSynth::Coord::Type & Coord<Leftwards>::orthComp() { return y(); }
template<> TexSynth::Coord::Type & Coord<Leftwards>::dirComp() { return x(); }
template<> int Coord<Leftwards>::dirStep() const { return -1; }

template<Direction D>
class Helper
{
public:
	Table<double> const * scores;

	explicit Helper(Table<double> const * _sc) : scores(_sc) {}

	bool coordIsValid(Coord<D> const & c) const { return c.x() >= 0 && c.y() >= 0 && c.x() < scores->width() && c.y() < scores->height(); }

	size_t dirSize() const { return Coord<D>(scores->width(), scores->height()).dirComp(); }
	size_t orthSize() const { return Coord<D>(scores->width(), scores->height()).orthComp(); }

	Coord<D> minCoord() const { return Coord<D>(0, 0); }
	Coord<D> maxCoord() const { return Coord<D>(scores->width() - 1, scores->height() - 1); }

	// relative to Direction
	//   TL | TR   TL   BL   BL ^ BR   BL   TL
	//      |F     ------>F     |R     <------R
	//   BL V BR   TR   BR   TL | TR   BR   TR
	// T at start; B at end in 'dir' Direction
	// L at zero; R at max in 'orth' Direction
	// F = Forward; R = Reverse
	Coord<D> topLeftCoord() const { Coord<D> c(maxCoord()); return c.dirStep() > 0 ? c.setOrth(0).setDir(0) : c.setOrth(0); }
	Coord<D> topRightCoord() const { Coord<D> c(maxCoord()); return c.dirStep() > 0 ? c.setDir(0)            : c;            }
	Coord<D> bottomLeftCoord() const { Coord<D> c(maxCoord()); return c.dirStep() > 0 ? c.setOrth(0) : c.setOrth(0).setDir(0); }
	Coord<D> bottomRightCoord() const { Coord<D> c(maxCoord()); return c.dirStep() > 0 ? c            : c.setOrth(0);           }

	double lookup(Coord<D> const & c) const { return (*scores)(c.x(), c.y()); }

	Coord<D> prevNeg(Coord<D> const & c) const { return Coord<D>(c).advPrev().advNeg(); }
	Coord<D> prevPos(Coord<D> const & c) const { return Coord<D>(c).advPrev().advPos(); }
	Coord<D> prevStr(Coord<D> const & c) const { return Coord<D>(c).advPrev(); }

	Coord<D> bestPrevNeg(Coord<D> const & c) const { return lookup(prevNeg(c)) < lookup(prevStr(c)) ? prevNeg(c) : prevStr(c); }
	Coord<D> bestPrevPos(Coord<D> const & c) const { return lookup(prevPos(c)) < lookup(prevStr(c)) ? prevPos(c) : prevStr(c); }

	Coord<D> bestPrevBoth(Coord<D> const & c) { return lookup(prevNeg(c)) < lookup(bestPrevPos(c)) ? prevNeg(c) : bestPrevPos(c); }
};

}

template<Seams::Direction D>
struct Seam
{
	std::vector<uint> indices;

	Seam(uint _s) : indices(_s) {}

	// NOTE: WARNING: no bounds checking
	uint & operator[](uint _i) { return indices[_i]; }
	uint operator[](uint _i) const { return indices[_i]; }

	std::vector<uint>::iterator begin() { return indices.begin(); }
	std::vector<uint>::iterator end() { return indices.end(); }

	// TODO: make all this use ImagePatches instead
	template<uint N>
	void modifyPatch(CImg<float> & _patch) const
	{
		CImg<float> ones(_patch.width(), _patch.height(), 1, _patch.spectrum(), 255);
		this->modifyPatch(_patch, ones);
	}

	// mask >= 1 marks areas that can be modified
	// i.e. where mask(x,y) is < 1.f, result is always patch(x,y)
	void modifyPatch(CImg<float> & _patch, CImg<float> const & _mask) const
	{
		Table<double> ref(_patch.width(), _patch.height());
		Seams::Helper<D> helper(&ref);

		if ( helper.dirSize() != indices.size() ) printf("DIMENSION MISMATCH: SILENT FAIL!\n");

		for ( Seams::Coord<D> base = helper.topLeftCoord(); helper.coordIsValid(base); base.advNext() )
			for ( Seams::Coord<D> c(base); c.orthComp() < helper.orthSize(); c.advPos() )
			{
				if ( _mask(c.x(), c.y(), 0) >= 1 )
				{
					float v = getVal(c, indices[c.dirComp()]);
					for ( int k = 0; k < _patch.spectrum(); ++k )
						_patch(c.x(), c.y(), k) = float(_patch(c.x(), c.y(), k)) * v;
				}
			}
	}

	float getVal(Seams::Coord<D> _c, uint _i) const
	{
		int d = _c.orthComp() > _i ? _c.orthComp() - _i : _i - _c.orthComp();
		float ret = 0;

		if ( d == 0 ) ret = 0.5;
		else if ( d == 1 ) ret = 0.2;

		if ( (D == Seams::Downwards || D == Seams::Leftwards) == (_c.orthComp() > indices[_c.dirComp()]) ) ret = 1.f - ret;
		return ret;
	}
};

//! Returns index of pixel on seam from either left edge (for Up and Down) or top edge (for Left or Right).
template<Seams::Direction D>
Seam<D> findMinSeam(Table<float> const & _costs)
{
	Table<double> cumulative(_costs.width(), _costs.height(), 0.f);
	Table<uint> refs(_costs.width(), _costs.height(), 0);

	Seams::Helper<D> helper(&cumulative);

	// first line
	for ( Seams::Coord<D> c = helper.topLeftCoord(); helper.coordIsValid(c); c.advPos() )
		cumulative(c.x(), c.y()) = _costs(c.x(), c.y());

	// rest
	Seams::Coord<D> base = helper.topLeftCoord();
	for ( base.advNext(); helper.coordIsValid(base); base.advNext() )
	{
		Seams::Coord<D> c(base);

		// first pixel
		Seams::Coord<D> m = helper.bestPrevPos(c);
		cumulative(c.x(), c.y()) = _costs(c.x(), c.y()) + cumulative(m.x(), m.y());
		refs(c.x(), c.y()) = m.orthComp() - c.orthComp();

		// middle pixels
		for ( c.advPos(); c.orthComp() < helper.orthSize() - 1; c.advPos() )
		{
			m = helper.bestPrevBoth(c);
			cumulative(c.x(), c.y()) = _costs(c.x(), c.y()) + cumulative(m.x(), m.y());
			refs(c.x(), c.y()) = m.orthComp() - c.orthComp();
		}

		// last pixel
		m = helper.bestPrevNeg(c);
		cumulative(c.x(), c.y()) = _costs(c.x(), c.y()) + cumulative(m.x(), m.y());
		refs(c.x(), c.y()) = m.orthComp() - c.orthComp();
	}

	// backtrack
	std::vector<Seams::Coord<D> > seamCoords(helper.dirSize());
	Seams::Coord<D> ptr = helper.bottomLeftCoord();
	seamCoords[ptr.dirComp()] = ptr;
	for ( ; helper.coordIsValid(ptr); ptr.advPos() )
		if ( cumulative(ptr.x(), ptr.y()) < cumulative(seamCoords[ptr.dirComp()].x(), seamCoords[ptr.dirComp()].y()) ) seamCoords[ptr.dirComp()] = ptr;

	ptr = seamCoords[ptr.dirComp()];
	Seams::Coord<D> pxy = ptr;
	for ( pxy.advPrev(); helper.coordIsValid(pxy); pxy.advPrev() )
	{
		ptr.orthComp() += refs(ptr.x(), ptr.y());
		ptr.advPrev();
		seamCoords[pxy.dirComp()] = ptr;
	}

	// copy to output
	Seam<D> seam(seamCoords.size());
	auto sit = seam.begin();
	auto cit = seamCoords.begin();
	for ( ; cit != seamCoords.end(); ++cit, ++sit ) *sit = (*cit).orthComp();

#if ( 0 )
	// DEBUG
	for ( int j = 0; j < _costs.height(); ++j )
	{
		printf("%i: ", j);

		for ( int i = 0; i < _costs.width(); ++i )
			printf("%6.1f  ", _costs(i, j));

		printf("\n");
	}
	printf("\n");

	for ( Seams::Coord<D> p = helper.topLeftCoord(); helper.coordIsValid(p); p.newOrth(0).advNext() )
	{
		printf("%u: ", seam[p.dirComp()]);

		for ( ; helper.coordIsValid(p); p.advPos() )
			printf("%6.1f%s ", _costs(p.x(), p.y()), p.orthComp() == seam[p.dirComp()] ? "*" : " ");

		printf("\n");
	}
	printf("\n");

	printf("   ");
	for ( Seams::Coord<D> p = helper.bottomLeftCoord(); helper.coordIsValid(p); p.advPos() )
		printf("%7.1f%s ", costs(p.x(), p.y()), p.orthComp() == seam[p.dirComp()] ? "*" : " ");

	printf("\n\n");
	// DEBUG END
#endif

	return seam;
}


} // namespace TexSynth

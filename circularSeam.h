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
#include <algorithm>
#include <utility>

#include "CImg.h"
using cimg_library::CImg;

#include "table.h"
#include "patch.h"


namespace TexSynth
{

namespace CircSeams
{

enum Direction { Downwards = 1, Rightwards, Upwards, Leftwards };
std::ostream & operator<<(std::ostream & _stream, Direction _dir)
{
	switch ( _dir )
	{
	case ( Downwards ):_stream << "Downwards"; break;
	case ( Upwards ):_stream << "Upwards"; break;
	case ( Rightwards ):_stream << "Rightwards"; break;
	case ( Leftwards ):_stream << "Leftwards"; break;
	default:_stream << "Unknown"; break;
	}
	return _stream;
}

bool isVertical(Direction d) { return d == Downwards || d == Upwards; }
bool isHorizontal(Direction d) { return !isVertical(d); }
bool isForward(Direction d) { return d == Downwards || d == Rightwards; }
bool isReverse(Direction d) { return !isForward(d); }

// CircularSeams are traversed in an anticlockwise direction
// (with the inside on the left)
//   SO | SI   SI   EI   EI ^ EO   EO   SO
//      |F     ------>F     |R     <------R
//   EO V EI   SO   EO   SI | SO   EI   SI
// S at start; E at end in 'dir' Direction
// I on left; O on right in 'orth' Direction

struct DirCoord : public Coord
{
	explicit DirCoord(Direction _d) : Coord(0, 0), direction(_d) {}
	DirCoord(Coord::Type _x, Coord::Type _y, Direction _d) : Coord(_x, _y), direction(_d) {}

	DirCoord & advIn() { orthComp() += orthStep(); return *this; }
	DirCoord & advOut() { orthComp() -= orthStep(); return *this; }
	DirCoord & advBack() { dirComp() -= dirStep(); return *this; }
	DirCoord & advForward() { dirComp() += dirStep(); return *this; }

	DirCoord & setOrth(Coord::Type _i) { orthComp() = _i; return *this; }
	DirCoord & setDir(Coord::Type _i) { dirComp() = _i; return *this; }

	Coord::Type & dirComp() { return isVertical(direction) ? y() : x(); }
	Coord::Type dirComp() const { return isVertical(direction) ? y() : x(); }
	Coord::Type & orthComp() { return isVertical(direction) ? x() : y(); }
	Coord::Type orthComp() const { return isVertical(direction) ? x() : y(); }

	int dirStep() const { return isForward(direction) ? +1 : -1; } //!< Outer to Inner
	int orthStep() const { return isVertical(direction) ? dirStep() : -dirStep(); } //!< Start to End

	DirCoord & setToInner(Coord::Type _max) { orthComp() = orthStep() < 0 ? 0 : _max; return *this; }
	DirCoord & setToOuter(Coord::Type _max) { orthComp() = orthStep() > 0 ? 0 : _max; return *this; }
	DirCoord & setToStart(Coord::Type _max) { dirComp() = dirStep() > 0 ? 0 : _max; return *this; }
	DirCoord & setToEnd  (Coord::Type _max) { dirComp() = dirStep() < 0 ? 0 : _max; return *this; }

	bool moreInnerThan(DirCoord const & _that) const { return this->orthStep() > 0 ? this->orthComp() > _that.orthComp() : this->orthComp() < _that.orthComp(); }
	bool moreOuterThan(DirCoord const & _that) const { return this->orthStep() > 0 ? this->orthComp() < _that.orthComp() : this->orthComp() > _that.orthComp(); }
	bool moreForwardThan(DirCoord const & _that) const { return this->dirStep() > 0 ? this->dirComp() > _that.dirComp() : this->dirComp() < _that.dirComp(); }
	bool moreBackwardThan(DirCoord const & _that) const { return this->dirStep() > 0 ? this->dirComp() < _that.dirComp() : this->dirComp() > _that.dirComp(); }

	Direction direction; // Could be templated for efficiency (see Seam::Coord)
};

std::ostream & operator<<(std::ostream & _stream, DirCoord const & _coord)
{
	_stream << "[" << _coord.direction << ": " << _coord.x() << ", " << _coord.y() << "]";
	return _stream;
}

class Helper
{
public:
	Helper(Table<double> const * _sc) : m_scores(_sc) {}

	bool coordIsValid(DirCoord const & c) const { return c.x() >= 0 && c.y() >= 0 && c.x() < m_scores->width() && c.y() < m_scores->height(); }

	double lookup(DirCoord const & c) const { return (*m_scores)(c.x(), c.y()); }

	Coord::Type dirSize(Direction _d) const { return DirCoord(m_scores->width(), m_scores->height(), _d).dirComp(); }
	Coord::Type orthSize(Direction _d) const { return DirCoord(m_scores->width(), m_scores->height(), _d).orthComp(); }

	DirCoord minCoord(Direction _d) const { return DirCoord(0, 0, _d); }
	DirCoord maxCoord(Direction _d) const { return DirCoord(m_scores->width() - 1, m_scores->height() - 1, _d); }
	Coord::Type maxOrth(Direction _d) const { return maxCoord(_d).orthComp(); }
	Coord::Type maxDir(Direction _d) const { return maxCoord(_d).dirComp(); }

	DirCoord startOutsideCoord(Direction _d) const { return DirCoord(_d).setToStart(maxDir(_d)).setToOuter(maxOrth(_d)); }
	DirCoord startInsideCoord(Direction _d) const { return DirCoord(_d).setToStart(maxDir(_d)).setToInner(maxOrth(_d)); }
	DirCoord endOuterCoord(Direction _d) const { return DirCoord(_d).setToEnd(maxDir(_d)).setToOuter(maxOrth(_d)); }
	DirCoord endInnerCoord(Direction _d) const { return DirCoord(_d).setToEnd(maxDir(_d)).setToInner(maxOrth(_d)); }

	Coord::Type distFromStart(DirCoord const & _c) const
	{
		Coord::Type s = startOutsideCoord(_c.direction).dirComp();
		Coord::Type c = _c.dirComp();
		return c > s ? c - s : s - c;
	}

	Coord::Type distFromOutside(DirCoord const & _c) const
	{
		Coord::Type o = startOutsideCoord(_c.direction).orthComp();
		Coord::Type c = _c.orthComp();
		return c > o ? c - o : o - c;
	}

	bool isWithinItsQuadrant(DirCoord const & _c) const
	{
		return distFromStart(_c) < dirSize(_c.direction) - distFromOutside(_c) && distFromStart(_c) > distFromOutside(_c) && coordIsValid(_c);
	}

	DirCoord prevIn(DirCoord const & c) const { return DirCoord(c).advBack().advIn(); }
	DirCoord prevOut(DirCoord const & c) const { return DirCoord(c).advBack().advOut(); }
	DirCoord prevStr(DirCoord const & c) const { return DirCoord(c).advBack(); }

	DirCoord bestPrevInStr(DirCoord const & c) const { return lookup(prevIn(c)) < lookup(prevStr(c)) ? prevIn(c) : prevStr(c); }
	DirCoord bestPrevOutStr(DirCoord const & c) const { return lookup(prevOut(c)) < lookup(prevStr(c)) ? prevOut(c) : prevStr(c); }

	DirCoord bestPrevBoth(DirCoord const & c) { return lookup(prevOut(c)) < lookup(bestPrevInStr(c)) ? prevOut(c) : bestPrevInStr(c); }

private:
	Table<double> const * m_scores;
};

} // end namespace CircSeams

struct CircSeam
{
	std::vector<uint> indices;

	CircSeam(uint _s) : indices(_s) {}

	// NOTE: WARNING: no bounds checking
	uint & operator[](uint _i) { return indices[_i]; }
	uint operator[](uint _i) const { return indices[_i]; }

	std::vector<uint>::iterator begin() { return indices.begin(); }
	std::vector<uint>::iterator end() { return indices.end(); }

/*
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
		CircSeams::Helper helper(&ref);

		if ( helper.dirSize(CircSeams::Downwards) != indices.size() ) printf("DIMENSION MISMATCH: SILENT FAIL!\n");

		for ( CircSeams::DirCoord base = helper.startOutsideCoord(CircSeams::Downwards); helper.coordIsValid(base); base.advForward() )
			for ( CircSeams::DirCoord c(base); c.orthComp() < helper.orthSize(CircSeams::Downwards); c.advIn() )
			{
				if ( _mask(c.x(), c.y(), 0) >= 1 )
				{
					float v = getVal(c, indices[c.dirComp()]);
					for ( int k = 0; k < _patch.spectrum(); ++k )
						_patch(c.x(), c.y(), k) = float(_patch(c.x(), c.y(), k)) * v;
				}
			}
	}

	float getVal(CircSeams::DirCoord _c, uint _i) const
	{
		int d = _c.orthComp() > _i ? _c.orthComp() - _i : _i - _c.orthComp();
		float ret = 0;

		if ( d == 0 ) ret = 0.5;
		else if ( d == 1 ) ret = 0.2;

		if ( (_c.direction == CircSeams::Downwards || _c.direction == CircSeams::Leftwards) == (_c.orthComp() > indices[_c.dirComp()]) ) ret = 1.f - ret;
		return ret;
	}
*/
};

//! Minimum cost seam all the way around, not using any entries where _mask < 255.f
CircSeam findMinCircSeam(Table<float> const & _costs, Table<float> & _mask)
{
	// increase size of the working table by 1 in every direction
	// to allow for prohibited areas going to the edge
	Table<float> costs(_costs.width() + 2, _costs.height() + 2, pow(2.f, 17.f)); // 256^2 * 2
	Table<double> cumulative(costs.width(), costs.height(), 0.f);
	Table<CircSeams::DirCoord> refs(costs.width(), costs.height(), CircSeams::DirCoord(CircSeams::Downwards));

	// TODO: should be some sort of global variable
	float epsilon = 0.1f;

	// fill costs
	for ( uint j = 0; j < _costs.height(); ++j )
		for ( uint i = 0; i < _costs.width(); ++i )
			costs(i + 1, j + 1) = _mask(i, j) < 255.f - epsilon ? std::numeric_limits<float>::infinity() : _costs(i, j);

	std::cout << "\n---\n";
	costs.printn();
	std::cout << "---" << std::endl;

	CircSeams::Helper helper(&cumulative);

	// order of quadrants
	std::vector<CircSeams::Direction> dirs = { CircSeams::Downwards, CircSeams::Rightwards, CircSeams::Upwards, CircSeams::Leftwards };

	// initialize cumulatives
	for ( CircSeams::DirCoord i = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(i); i.advBack().advIn() )
		cumulative(i.x(), i.y()) = costs(i.x(), i.y()) == std::numeric_limits<float>::infinity() ? std::numeric_limits<float>::infinity() : 0;

	// dynamic programming
	for ( CircSeams::Direction d : dirs )
	{
		for ( CircSeams::DirCoord base = helper.startOutsideCoord(d); !base.moreForwardThan(helper.endOuterCoord(d)); base.advForward() )
		{
			for ( CircSeams::DirCoord curr = base; helper.isWithinItsQuadrant(curr); curr.advIn() )
			{
				CircSeams::DirCoord b(d);
				// TODO: tidy up
				b = helper.distFromOutside(curr) > 0 ? helper.bestPrevBoth(curr) : helper.bestPrevInStr(curr);
				if ( !helper.isWithinItsQuadrant(b) )
					b = helper.distFromOutside(curr) > 0 ? helper.bestPrevOutStr(curr) : helper.prevStr(curr);
				cumulative(curr.x(), curr.y()) = costs(curr.x(), curr.y()) + cumulative(b.x(), b.y());
				refs(curr.x(), curr.y()) = b;
			}
		}
	}

	// backtracking
	typedef std::pair<double, CircSeams::DirCoord> DDC;
	struct compDDC { bool operator()(DDC const & a, DDC const & b) { return a.first > b.first; } };
	std::vector<DDC> heap;

	for ( CircSeams::DirCoord e = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(e); e.advBack().advIn() )
	{
		double v = cumulative(e.x(), e.y());
		if ( v < std::numeric_limits<double>::infinity() )
		{
			heap.push_back(std::make_pair(v, e));
			std::push_heap(heap.begin(), heap.end(), compDDC());
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compDDC());
	// heap in reverse order
	for ( auto seamEndPtr = heap.rbegin(); seamEndPtr != heap.rend(); ++seamEndPtr )
	{
		std::vector<CircSeams::DirCoord> seam;
		bool leftLastQuadrant = false;
		CircSeams::DirCoord curr = (*seamEndPtr).second;
		CircSeams::DirCoord test(curr.x(), curr.y(), dirs.back());
		while ( leftLastQuadrant != helper.isWithinItsQuadrant(test) )
		{
			std::cout << "Adding " << curr << std::endl;

			seam.push_back(curr);
			curr = refs(curr.x(), curr.y());
			test.x() = curr.x();
			test.y() = curr.y();
			if ( !leftLastQuadrant && !helper.isWithinItsQuadrant(test) )
				leftLastQuadrant = true;

			std::cout << "Check: test = " << test <<
						 "; left = " << (leftLastQuadrant ? "true" : "false") <<
						 "; within = " << (helper.isWithinItsQuadrant(test) ? "true" : "false") <<
						std::endl;
		}

		std::cout << "***" << (*seamEndPtr).first << ": " << std::endl;
		for ( auto itr = seam.rbegin(); itr != seam.rend(); ++itr )
			std::cout << *itr << std::endl;

	}

	return CircSeam(0);
}


} // namespace TexSynth

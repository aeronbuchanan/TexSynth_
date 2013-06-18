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

typedef std::vector<DirCoord> CircSeam;

template<typename T>
class Helper
{
public:
	Helper(Table<T> const * _sc) : m_scores(_sc) {}

	bool coordIsValid(DirCoord const & c) const { return c.x() >= 0 && c.y() >= 0 && c.x() < m_scores->width() && c.y() < m_scores->height(); }

	T lookup(DirCoord const & c) const { return (*m_scores)(c.x(), c.y()); }

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

	DirCoord bestInDirectionOf(DirCoord const & src, DirCoord const & tgt)
	{
		DirCoord s1 = src;
		DirCoord s2 = src;
		DirCoord s3 = src;
		if ( src.moreForwardThan(tgt) )
		{
			s1.advBack();
			s3 = s1;
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
				s3.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
				s3.advOut();
			}
			else
			{
				s2 = s1;
			}
		}
		else if ( src.moreBackwardThan(tgt) )
		{
			s1.advForward();
			s3 = s1;
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
				s3.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
				s3.advOut();
			}
			else
			{
				s2 = s1;
			}
		}
		else
		{
			if ( src.moreOuterThan(tgt) )
			{
				s2.advIn();
			}
			else if ( src.moreInnerThan(tgt) )
			{
				s2.advOut();
			}
			else
			{
				// src IS tgt
			}
			s1 = s2;
			s3 = s2;
		}

		T d1 = lookup(s1);
		T d2 = lookup(s2);
		T d3 = lookup(s3);

		return d1 > d2 ? (d1 > d3 ? s1 : (d3 > d2 ? s3 : s2)) : (d2 > d3 ? s2 : (d3 > d1 ? s3 : s1));
	}

	void display(CircSeam const & cs)
	{
		for ( uint j = 0; j < m_scores->height(); ++j )
		{
			for ( uint i = 0; i < m_scores->width(); ++i )
			{
				bool found = false;
				for ( auto itr : cs )
					found |= ( itr.x() == i && itr.y() == j );
				printf("%10.0f%s ", (*m_scores)(i,j), found ? "*" : " ");
			}
			printf("\n");
		}
	}

private:
	Table<T> const * m_scores;
};

//! Minimum cost seam all the way around, not using any entries where _mask < 255.f
Table<float> findMin(Table<float> const & _costs, Table<float> const & _mask)
{
	// Note: this is an approximation
	// The loop closure uses heuristics to get an answer rather than anything that can garrantee optimality

	// increase size of the working table by 1 in every direction
	// to allow for prohibited areas going to the edge
	Table<float> costs(_costs.width() + 2, _costs.height() + 2, pow(2.f, 17.f)); // 256^2 * 2
	Table<double> cumulative(costs.width(), costs.height(), 0.f);
	Table<DirCoord> refs(costs.width(), costs.height(), DirCoord(Downwards));

	// TODO: should be some sort of global variable
	float epsilon = 0.1f;

	// fill costs
	for ( uint j = 0; j < _costs.height(); ++j )
		for ( uint i = 0; i < _costs.width(); ++i )
			costs(i + 1, j + 1) = _mask(i, j) < 255.f - epsilon ? std::numeric_limits<float>::infinity() : _costs(i, j);

	std::cout << "\n---\n";
	costs.printn();
	std::cout << "---" << std::endl;

	Helper<double> helper(&cumulative);

	// order of quadrants
	std::vector<Direction> dirs = { Downwards, Rightwards, Upwards, Leftwards };

	// initialize cumulatives
	for ( DirCoord i = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(i); i.advBack().advIn() )
		cumulative(i.x(), i.y()) = costs(i.x(), i.y()) == std::numeric_limits<float>::infinity() ? std::numeric_limits<float>::infinity() : 0;

	double bestCost = 0;
	CircSeam bestSeam;

	// dynamic programming
	for ( Direction d : dirs )
	{
		for ( DirCoord base = helper.startOutsideCoord(d); !base.moreForwardThan(helper.endOuterCoord(d)); base.advForward() )
		{
			// create a default 'around the outside' seam
			bestSeam.push_back(base);
			bestCost += costs(base.x(), base.y());

			for ( DirCoord curr = base; helper.isWithinItsQuadrant(curr); curr.advIn() )
			{
				DirCoord b(d);
				// TODO: tidy up
				b = helper.distFromOutside(curr) > 0 ? helper.bestPrevBoth(curr) : helper.bestPrevInStr(curr);
				// Must not go into central pixels
				if ( !helper.isWithinItsQuadrant(b) )
					b = helper.distFromOutside(curr) > 0 ? helper.bestPrevOutStr(curr) : helper.prevStr(curr);
				cumulative(curr.x(), curr.y()) = costs(curr.x(), curr.y()) + cumulative(b.x(), b.y());
				refs(curr.x(), curr.y()) = b;
			}
		}
	}

	// backtracking
	typedef std::pair<double, DirCoord> DDC;
	struct compDDC { bool operator()(DDC const & a, DDC const & b) { return a.first > b.first; } };
	std::vector<DDC> heap;

	for ( DirCoord e = helper.endOuterCoord(dirs.back()); helper.isWithinItsQuadrant(e); e.advBack().advIn() )
	{
		double v = cumulative(e.x(), e.y());
		if ( v < std::numeric_limits<double>::infinity() )
		{
			heap.push_back(std::make_pair(v, e));
			std::push_heap(heap.begin(), heap.end(), compDDC());
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compDDC());

	Helper<float> costsHelper(&costs); // DEBUG

	// heap in reverse order
	for ( auto seamEndPtr = heap.rbegin(); seamEndPtr != heap.rend(); ++seamEndPtr )
	{
		printf("---\n");

		CircSeam thisSeam;
		bool leftLastQuadrant = false;
		double thisCost = (*seamEndPtr).first;
		DirCoord curr = (*seamEndPtr).second;
		DirCoord test(curr.x(), curr.y(), dirs.back());
		while ( leftLastQuadrant != helper.isWithinItsQuadrant(test) )
		{
			//std::cout << "Adding " << curr << std::endl;

			thisSeam.push_back(curr);
			curr = refs(curr.x(), curr.y());
			test.x() = curr.x();
			test.y() = curr.y();
			if ( !leftLastQuadrant && !helper.isWithinItsQuadrant(test) )
				leftLastQuadrant = true;

			//std::cout << "Check: test = " << test << "; left = " << (leftLastQuadrant ? "true" : "false") << "; within = " << (helper.isWithinItsQuadrant(test) ? "true" : "false") << std::endl;
		}

		//std::cout << "***" << (*seamEndPtr).first << ": " << std::endl;
		//for ( auto itr = thisSeam.rbegin(); itr != thisSeam.rend(); ++itr )
		//	std::cout << *itr << std::endl;
		//std::cout << "***" << std::endl;

		struct L1
		{
			static int dist(DirCoord const & s, DirCoord const & e)
			{
				return std::max( s.x() > e.x() ? s.x() - e.x() : e.x() - s.x(), s.y() > e.y() ? s.y() - e.y() : e.y() - s.y() );
			}
		};

		// cut corners on quadrant transitions
		// (because of isWithinItsQuadrant(b) during dynamic programming)
		for ( int i = 0; i < thisSeam.size(); ++i )
		{
			int j = i - 2;
			if ( j < 0 )
				j += thisSeam.size();

			if ( thisSeam[i].direction != thisSeam[j].direction )
			{
				int k = i - 1;
				if ( k < 0 )
					k += thisSeam.size();

				if ( L1::dist(thisSeam[i], thisSeam[j]) == 1 )
				{
					auto ptr = thisSeam.begin() + k;

					std::cout << "REMOVING " << k << ":" << *ptr << std::endl;

					thisCost -= costs((*ptr).x(), (*ptr).y());
					thisSeam.erase(ptr);
					i = std::max(i - 2, 0); // correct
				}
			}
		}

		std::cout << "###" << (*seamEndPtr).first << ": " << std::endl;
		//for ( auto itr = thisSeam.rbegin(); itr != thisSeam.rend(); ++itr )
		//	std::cout << *itr << std::endl;
		costsHelper.display(thisSeam);
		std::cout << "###" << std::endl;

		if (  thisCost < bestCost )
		{
			if ( L1::dist(thisSeam.front(), thisSeam.back()) == 1 )
			{
				// done
				bestCost = thisCost;
				bestSeam = thisSeam;
			}
			else
			{
				printf("This needs work...\n");



				// try connecting from back (greedy approach)
				while ( L1::dist(thisSeam.front(), thisSeam.back()) > 1 )
				{
					DirCoord s = thisSeam.front();
					DirCoord e = thisSeam.back();

					std::cout << "Advancing " << e << " towards " << s << std::endl;

					e = costsHelper.bestInDirectionOf(e, s);
					thisSeam.push_back( e );

					std::cout << " ==> " << e << std::endl;

					thisCost += costs(e.x(), e.y());

					std::cout << "dist = " << L1::dist(thisSeam.front(), thisSeam.back()) << " => cost = " << thisCost << std::endl;
				}

				if ( thisCost < bestCost )
				{
					bestCost = thisCost;
					bestSeam = thisSeam;
				}

				// try other things?
			}

		}

	}

	std::cout << "~~~" << bestCost << ": " << std::endl;
	//for ( auto itr = bestSeam.rbegin(); itr != bestSeam.rend(); ++itr )
	//	std::cout << *itr << std::endl;
	costsHelper.display(bestSeam);

	return Table<float>(0, 0);
}

} // end namespace CircSeams

} // namespace TexSynth

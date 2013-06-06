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

#include <iostream>
#include "CImg.h"
using cimg_library::CImg;

typedef unsigned char uchar;

namespace TexSynth
{

//! Another fixed vector class
template<class T, uint N>
class VecN
{

#define FOREACH_(I,N) for( uint (I) = 0; (I) < (N); ++(I) )

public:
	VecN() { setAll(T()); }
	explicit VecN(T const _v) { setAll(_v); }
	explicit VecN(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } }
	template<class Y, uint M>
	VecN(VecN<Y, M> const & _v) { setAll(0); FOREACH_(i, std::min(N,M)) { m_v[i] = _v[i]; } }

	VecN<T, N> & operator=(T const _v) { return setAll(_v); }
	VecN<T, N> & operator=(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }
	template<class Y>
	VecN<T, N> & operator=(VecN<Y,N> const & _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }

	// TODO: WARNING: silently returns last element if index out-of-bounds
	T & operator[](uint _i) { return m_v[std::min(_i, N - 1)]; }
	T operator[](uint _i) const { return m_v[std::min(_i, N - 1)]; }

	template<class Y>
	VecN<T, N> & operator+=(VecN<Y, N> const & _v) { FOREACH_(i,N) { m_v[i] += _v[i]; } return *this; }
	VecN<T, N> & operator+=(T _v) { return *this += VecN<T, N>(_v); }

	template<class Y>
	VecN<T, N> & operator-=(VecN<Y, N> const & _v) { return *this += (-_v); }
	VecN<T, N> & operator-=(T _v) { return *this -= VecN<T, N>(_v); }

	VecN<T, N> & operator*=(T _v) { FOREACH_(i,N) { m_v[i] *= _v; } return *this; }
	VecN<T, N> & operator/=(T _v) { FOREACH_(i,N) { m_v[i] /= _v; } return *this; }

	template<class Y>
	VecN<T, N> operator+(VecN<Y,N> const & _v) const { return VecN<Y,N>(*this) += _v; }
	VecN<T, N> operator+(T _v) const { return VecN<T, N>(*this) += _v; }

	template<class Y>
	VecN<T, N> operator-(VecN<Y, N> const & _v) const { return VecN<Y,N>(*this) -= _v; }
	VecN<T, N> operator-(T _v) const { return VecN<T, N>(*this) -= _v; }

	VecN<T, N> operator*(T _v) const { return VecN<T, N>(*this) *= _v; }
	VecN<T, N> operator/(T _v) const { return VecN<T, N>(*this) /= _v; }

	VecN<T, N> operator-() const { return VecN<T, N>(*this) *= -1; }

	VecN<T, N> & operator++() { return *this += 1; }
	VecN<T, N> & operator--() { return *this -= 1; }

	VecN<T, N> operator++(int) { VecN<T, N> pre(*this); ++(*this); return pre; }
	VecN<T, N> operator--(int) { VecN<T, N> pre(*this); --(*this); return pre; }

	template<class Y>
	T dot(VecN<Y,N> const & _v) const { T k = T(); FOREACH_(i,N) { k += m_v[i] * _v[i]; } return k; }

	//! Element-wise multiplication (dot-multiply of matlab)
	template<class Y>
	VecN<T, N> maskWith(VecN<Y,N> const & _m) { VecN<T, N> v(*this); FOREACH_(i,N) { v[i] *= _m[i]; } return v; }

	T magnitudeSqrd() const { return this->dot(*this); }
	T magnitude() const { return sqrt(magnitudeSqrd()); }
	template<class Y>
	T maskedMagSqrd(VecN<Y, N> const & _m) const { return this->maskWith(_m).dot(*this); }

	VecN<T, N> & makeAbsolute() { FOREACH_(i,N) { m_v[i] = std::abs(m_v[i]); } return *this; }
	VecN<T, N> asAbsolute() const { return VecN<T, N>(*this).makeAbsolute(); }

	void print() const { std::cout << "("; FOREACH_(i, (N - 1)) { std::cout << m_v[i] << ", "; } std::cout << m_v[N - 1] << ")"; }
	void printn() const { print(); std::cout << std::endl; }

	static inline uint size() { return N; }

private:
	VecN<T, N> & setAll(T _v) { FOREACH_(i,N) { m_v[i] = _v; } return *this; }

	T m_v[N];

#undef FOREACH_
};

template<uint N, uint M = 3, typename T = float>
class ImagePatch
{
	typedef VecN<T, N * N * M> VEC;

	// ordering (inner -> outer) is (c,i,j)
	// optimizes conversion to gray, etc
	// bad for conversion to/from CImg (I think)
	static inline uint kCoord(uint _i, uint _j, uint _c) { return (((_j * N) + _i) * M) + _c; } //!< convert (i,j,c) coord to linear index
	static inline uint iCoord(uint _k) { return ( _k / M) % N; } //!< convert linear index to horizontal patch coord
	static inline uint jCoord(uint _k) { return ( _k / M) / N; } //!< convert linear index to vertical patch coord
	static inline uint cCoord(uint _k) { return _k % M; } //!< convert linear index to color channel
	static inline uint xCoord(uint _x, uint, uint _k) { return _x + iCoord(_k); }
	static inline uint yCoord(uint, uint _y, uint _k) { return _y + jCoord(_k); }

public:
	ImagePatch() {}
	explicit ImagePatch(T _v) : m_v(_v) {}
	explicit ImagePatch(VEC _v) : m_v(_v) {}
	template<class U>
	ImagePatch(CImg<U> const & _img, uint _x, uint _y) { this->extractFrom(_img, _x, _y); }

	T & operator[](uint _k) { return m_v[_k]; }
	T operator[](uint _k) const { return m_v[_k]; }

	T & operator()(uint _i, uint _j, uint _c) { return m_v[kCoord(_i, _j, _c)]; }
	T operator()(uint _i, uint _j, uint _c) const { return m_v[kCoord(_i, _j, _c)]; }

	VEC asVecN() const { return m_v; }
	VecN<T, N * N> channelAsVecN(uint _c) const
	{
		VecN<T, N * N> v;
		for ( uint k = 0, i = 0; k < m_v.size(); ++k )
			if ( cCoord(k) == _c )
				v[i++] = m_v[k];
		return v;
	}

	template<class U>
	bool extractFrom(CImg<U> const & _img, uint _x, uint _y)
	{
		bool badCoords = false;

		if ( _x > _img.width() - N - 1 ) { printf(" [extract: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N - 1 ) { printf(" [extract: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < m_v.size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					m_v[k] = _img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k));
		}

		return badCoords;
	}

	template<class U>
	bool insertInto(CImg<U> & _img, uint _x, uint _y) const
	{
		bool badCoords = false;

		if ( _x > _img.width() - N - 1 ) { printf(" [insert: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N - 1 ) { printf(" [insert: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < m_v.size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					_img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k)) = m_v[k];
		}

		return badCoords;
	}

	ImagePatch<N, M, T> & convertToGray()
	{
		if ( M > 1 )
		{
			// TODO: seems to be coming out too dark...
			for ( uint j = 0; j < N; ++j )
			{
				for ( uint i = 0; i < N; ++i )
				{
					T v = T();
					for ( uint c = 0; c < M; ++c )
						v += m_v[kCoord(i, j, c)] / M;
					for ( uint c = 0; c < M; ++c )
						m_v[kCoord(i, j, c)] = v;
				}
			}
		}
		return *this;
	}

	ImagePatch<N, M, T> & grayVersion() { return ImagePatch<N, M, T>(*this).convertToGray(); }

	// TODO: refactor following for inter-type compatibility

	ImagePatch<N, M, T> & maskWith(ImagePatch<N, M, T> _m) { m_v.maskWith(_m.m_v); return *this; }
	//ImagePatch<N, M, T> & maskWith(ImagePatch<N / M, 1, U> _m) { m_v.maskWith(ImagePatch<N, M, T>(_m)); return *this; }

	static T sqrdError(ImagePatch<N, M, T> const & a, ImagePatch<N, M, T> const & b) { return diff(a,b).m_v.magnitudeSqrd(); }

	static T maskedSqrdError(ImagePatch<N, M, T> const & a, ImagePatch<N, M, T> const & b, ImagePatch<N, M, T> const & m) { return diff(a,b).m_v.maskWith(m.m_v).magnitudeSqrd(); }

	static ImagePatch<N, M, T> diff(ImagePatch<N, M, T> const & a, ImagePatch<N, M, T> const & b) { return ImagePatch<N, M, T>(a.m_v - b.m_v); }
	static ImagePatch<N, M, T> diffSqrd(ImagePatch<N, M, T> const & a, ImagePatch<N, M, T> const & b) { return ImagePatch<N, M, T>( (a.m_v - b.m_v).maskWith(a.m_v - b.m_v) ); }

	static ImagePatch<N, M, T> abs(ImagePatch<N, M, T> const & a) { return ImagePatch<N, M, T>(a.m_v.asAbsolute()); }

#if 0
	bool test()
	{
		bool good = true;
		for ( uint j = 0; j < N; ++j )
			for ( uint i = 0; i < N; ++i )
				for ( uint c = 0; c < M; ++c )
				{
					uint k = kCoord(i,j,c);
					uint ic = iCoord(k);
					uint jc = jCoord(k);
					uint cc = cCoord(k);
					if ( i != ic || j != jc || c != cc )
					{
						good = false;
						printf("ERROR: (%u,%u,%u) -> [%u] -> (%u,%u,%u)\n", i, j, c, k, iCoord(k), jCoord(k), cCoord(k));
					}
				}
		return good;
	}
#endif

private:
	VEC m_v;

};



} // end namespace TexSynth

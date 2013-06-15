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
#include <string>

// TODO: make TexSynth image-library agnostic
#include "CImg.h"
using cimg_library::CImg;

typedef unsigned char uchar;

namespace TexSynth
{

//! Useful debug typename identifier
template<typename T> std::string identifyType(T const &) { return "unknown"; }
template<> std::string identifyType(int const &) { return "int"; }
template<> std::string identifyType(uint const &) { return "uint"; }
template<> std::string identifyType(char const &) { return "char"; }
template<> std::string identifyType(uchar const &) { return "uchar"; }
template<> std::string identifyType(float const &) { return "float"; }
template<> std::string identifyType(double const &) { return "double"; }


//! Another fixed vector class
template<class T, uint N>
class VecN
{

#define FOREACH_(I,N) for( uint (I) = 0; (I) < (N); ++(I) )

public:
	VecN() { setAll(T()); }
	explicit VecN(T const & _v) { setAll(_v); }
	explicit VecN(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } }

	template<class Y, uint M>
	VecN(VecN<Y, M> const & _v) { setAll(0); FOREACH_(i, std::min(N,M)) { m_v[i] = _v[i]; } }

	// TODO: accept {} list constructors
	VecN(T const & _v0, T const & _v1)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
	}
	VecN(T const & _v0, T const & _v1, T const & _v2)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
		getElement<2>() = _v2;
	}
	VecN(T const & _v0, T const & _v1, T const & _v2, T const & _v3)
	{
		getElement<0>() = _v0;
		getElement<1>() = _v1;
		getElement<2>() = _v2;
		getElement<3>() = _v3;
	}

	VecN<T, N> & operator=(T const _v) { return setAll(_v); }
	VecN<T, N> & operator=(T const * _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }
	template<class Y>
	VecN<T, N> & operator=(VecN<Y,N> const & _v) { FOREACH_(i,N) { m_v[i] = _v[i]; } return *this; }

	// TODO: WARNING: silently returns last element if index out-of-bounds
	T & operator[](uint _i) { return m_v[std::min(_i, N - 1)]; }
	T operator[](uint _i) const { return m_v[std::min(_i, N - 1)]; }

private:
	template<uint I>
	T & getElement() { static_assert(I < N, "Index out-of-bounds."); return m_v[I]; }
	template<uint I>
	T getElement() const { static_assert(I < N, "Index out-of-bounds."); return m_v[I]; }

public:
	T & x() { return getElement<0>(); }
	T x() const { return getElement<0>(); }
	T & y() { return getElement<1>(); }
	T y() const { return getElement<1>(); }
	T & z() { return getElement<2>(); }
	T z() const { return getElement<2>(); }
	T & w() { return getElement<3>(); }
	T w() const { return getElement<3>(); }

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
	VecN<T, N> maskedWith(VecN<Y, N> const & _m) const { VecN<T, N> v(*this); FOREACH_(i,N) { v[i] *= _m[i]; } return v; }
	//! In-place element-wise multiplication (dot-multiply of matlab)
	template<class Y>
	VecN<T, N> & maskWith(VecN<Y, N> const & _m) { FOREACH_(i,N) { m_v[i] *= _m[i]; } return *this; }

	T magnitudeSqrd() const { return this->dot(*this); }
	T magnitude() const { return sqrt(magnitudeSqrd()); }
	template<class Y>
	T maskedMagSqrd(VecN<Y, N> const & _m) const { return this->maskedWith(_m).dot(*this); }

	VecN<T, N> & makeAbsolute() { FOREACH_(i,N) { m_v[i] = std::abs(m_v[i]); } return *this; }
	VecN<T, N> asAbsolute() const { return VecN<T, N>(*this).makeAbsolute(); }

	void print() const { std::cout << "("; FOREACH_(i, (N - 1)) { std::cout << m_v[i] << ", "; } std::cout << m_v[N - 1] << ")"; }
	void printn() const { print(); std::cout << std::endl; }

	static inline uint size() { return N; }

protected:
	VecN<T, N> & setAll(T _v) { FOREACH_(i,N) { m_v[i] = _v; } return *this; }

	T m_v[N];

#undef FOREACH_
};

template<uint N, uint M = 3>
class ImagePatch
{	
protected:
	typedef VecN<float, N * N * M> ThisVecN;

	// ordering (inner -> outer) is (c,i,j)
	// optimizes conversion to gray, etc
	// bad for conversion to/from CImg (I think)
	static inline uint kCoord(uint _i, uint _j, uint _c) { return (((_j * N) + _i) * M) + _c; } //!< convert (i,j,c) coord to linear index
	static inline uint iCoord(uint _k) { return ( _k / M) % N; } //!< convert linear index to horizontal patch coord
	static inline uint jCoord(uint _k) { return ( _k / M) / N; } //!< convert linear index to vertical patch coord
	static inline uint cCoord(uint _k) { return _k % M; } //!< convert linear index to color channel
	static inline uint xCoord(uint _x, uint, uint _k) { return _x + iCoord(_k); }
	static inline uint yCoord(uint, uint _y, uint _k) { return _y + jCoord(_k); }

private:
	ThisVecN m_vec;

public:
	ImagePatch() {}
	explicit ImagePatch(float _v) : m_vec(_v) {}
	explicit ImagePatch(ThisVecN _v) : m_vec(_v) {}
	ImagePatch(CImg<float> const & _img, uint _x, uint _y) { this->extractFrom(_img, _x, _y); }

	// TODO: use size_t everywhere else too
	size_t size() const { return N * N * M; }
	size_t width() const { return N; }
	size_t height() const { return N; }
	size_t depth() const { return M; }

	float & operator[](uint _k) { return m_vec[_k]; }
	float operator[](uint _k) const { return m_vec[_k]; }

	float & operator()(uint _i, uint _j, uint _c) { return (*this)[kCoord(_i, _j, _c)]; }
	float operator()(uint _i, uint _j, uint _c) const { return (*this)[kCoord(_i, _j, _c)]; }

	bool extractFrom(CImg<float> const & _img, uint _x, uint _y)
	{
		bool badCoords = false;

		if ( _x > _img.width() - N ) { printf(" [extract: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N ) { printf(" [extract: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < ImagePatch::size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					(*this)[k] = _img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k));
		}

		return badCoords;
	}

	bool insertInto(CImg<float> & _img, uint _x, uint _y) const
	{
		bool badCoords = false;

		if ( _x > _img.width() - N ) { printf(" [insert: bad x] "); badCoords = true; }
		if ( _y > _img.height() - N ) { printf(" [insert: bad y] "); badCoords = true; }

		if ( !badCoords )
		{
			for ( uint k = 0; k < this->size(); ++k )
				if ( static_cast<int>(cCoord(k)) < _img.spectrum() )
					_img(xCoord(_x, _y, k), yCoord(_x, _y, k), 0, cCoord(k)) = (*this)[k];
		}

		return badCoords;
	}

	ImagePatch<N, M> & convertToGray()
	{
		if ( M > 1 )
		{
			// TODO: seems to be coming out too dark...
			for ( uint j = 0; j < N; ++j )
			{
				for ( uint i = 0; i < N; ++i )
				{
					float v = 0;
					for ( uint c = 0; c < M; ++c )
						v += (*this)[kCoord(i, j, c)] / 3.f;
					for ( uint c = 0; c < M; ++c )
						(*this)[kCoord(i, j, c)] = v / 3.f;
				}
			}
		}
		return *this;
	}

	ImagePatch<N, M> grayVersion() const { return ImagePatch<N, M>(*this).convertToGray(); }

	void print() const { m_vec.print(); }
	void printn() const { m_vec.printn(); }
	void save(std::string const & _name) const { this->save(_name.c_str()); }
	void save(char * const _name) const { CImg<float> t(N, N, 1, M); insertInto(t, 0, 0); t.save(_name); }

	// TODO: refactor following for inter-type compatibility

	ImagePatch<N, M> blendedWith(ImagePatch<N, M> const & _q, ImagePatch<N, M> const & _m) const
	{
		ImagePatch<N, M> p(*this);
		for ( uint i = 0; i < this->size(); ++i )
		{
			float a = _m[i] / 255.f;
			p[i] = (1.f - a) * _q[i] +  a * p[i];
		}
		return p;
	}

	ImagePatch<N, M> & maskWith(ImagePatch<N, M> _m) { m_vec.maskWith(_m.m_vec); return *this; }

	static float sqrdError(ImagePatch<N, M> const & _a, ImagePatch<N, M> const & _b) { return diff(_a, _b).m_v.magnitudeSqrd(); }

	static float maskedSqrdError(ImagePatch<N, M> const & _a, ImagePatch<N, M> const & _b, ImagePatch<N, M> const & _m) { return diff(_a, _b).maskWith(_m).m_vec.magnitudeSqrd(); }

	static ImagePatch<N, M> diff(ImagePatch<N, M> const & a, ImagePatch<N, M> const & b) { return ImagePatch<N, M>(a.m_vec - b.m_vec); }
	static ImagePatch<N, M> diffSqrd(ImagePatch<N, M> const & a, ImagePatch<N, M> const & b) { return ImagePatch<N, M>( (a.m_vec - b.m_vec).maskWith(a.m_vec - b.m_vec) ); }
	static ImagePatch<N, M> abs(ImagePatch<N, M> const & a) { return ImagePatch<N, M>(a.m_vec.asAbsolute()); }

/*
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
*/

};

typedef VecN<uint, 2> Coord;

} // end namespace TexSynth

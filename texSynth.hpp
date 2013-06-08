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

#include <limits>
#include <algorithm>

template<uint N, typename T>
template<class U>
uint TexSynth::TexSynther<N, T>::selectPatchFor(cimg_library::CImg<U> & _I,cimg_library:: CImg<U> & _M, uint _x, uint _y) const
{
	// printf("selectPatchFor\n");

	if ( m_patches.size() == 0 ) return 0; // TODO: throw, etc.

	typedef std::pair<double, uint> DUIP;
	struct compareDUIP { bool operator()(DUIP const & a, DUIP const & b) { return a.first > b.first; } };

	Patch p(_I, _x, _y);
	Patch m(_M, _x, _y);

	std::vector<DUIP> heap;
	double threshold = std::numeric_limits<double>::max();

	for ( uint i = 0; i < m_patches.size(); ++i )
	{
		if ( heap.size() > 0) threshold = heap.front().first * 1.3; // (1.0 * m.magnitudeSqrd());
		double d = Patch::maskedSqrdError(p, m_patches[i], m);
		if ( d < threshold )
		{
			heap.push_back(std::make_pair(d, i));
			std::push_heap(heap.begin(), heap.end(), compareDUIP());
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compareDUIP());

	int i = heap.size(); // heap is in reverse order
	uint count = 0;
	uint CANDIDATE_MAX = 100;

	while (	++count < CANDIDATE_MAX && i > 0 &&	heap[--i].first < threshold	);

	--count;
	if ( count == 0 ) return 0; // TODO: throw, etc.

	uint s = heap.size() - 1 - (rand() % count); // selection

	return heap[s].second;
}


template<uint N, typename T>
template<class U>
void TexSynth::TexSynther<N, T>::extendTextureIn(cimg_library::CImg<U> & _I, cimg_library:: CImg<U> & _M) const
{
	//printf("extendTextureIn\n");

	Patch ones(1);

	uint shift = m_patchWidth * 0.6;

	printf("I(%dx%d) and %d\n", _I.width(), _I.height(), m_patchWidth);

	for ( uint y = 0; y < _I.height() - m_patchWidth; y += shift )
		for ( uint x = 0; x < _I.width() - m_patchWidth; x += shift )
		{
			Patch p = m_patches[selectPatchFor(_I, _M, x, y)];

			// TODO: put the following somewhere more general
			Patch ii(_I, x, y);
			Patch mm(_M, x, y);
			Patch dd = Patch::diffSqrd(ii, p);
			for ( uint i = 0; i < dd.size(); ++i ) dd[i] = mm[i] != 1.f ? std::numeric_limits<T>::max() : dd[i];
			CImg<float> d(N, N, 1, 3);
			dd.insertInto(d, 0, 0);
			// TODO: do all four directions
			auto vSeam = findMinSeam<Seams::Downwards>(d);
			auto hSeam = findMinSeam<Seams::Leftwards>(d);
			d.fill(1.f);
			vSeam.modifyPatch(d);
			hSeam.modifyPatch(d);
			// Merge patches
			for ( uint j = 0; j < N; ++j )
				for ( uint i = 0; i < N; ++i )
					for ( uint c = 0; c < 3; ++c )
						p(i, j, c) = d(i, j) * p(i, j, c) + (1.f - d(i, j)) * ii(i, j, c);

			p.insertInto(_I, x, y);

			ones.insertInto(_M, x, y);
		}
}


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
		double d = Patch::maskedSqrdError(p, m_patches[i], m);
		if ( d < threshold )
		{
			heap.push_back(std::make_pair(d, i));
			std::push_heap(heap.begin(), heap.end(), compareDUIP());
			if ( d != 0 )
			{
				// set threshold as multiple of second highest, non zero error
				uint j = 1;
				while ( j < heap.size() && heap[j].first == 0 )
					j++;
				if ( j < heap.size() )
					threshold = heap[j].first * 1.1;
			}
		}
	}

	std::sort_heap(heap.begin(), heap.end(), compareDUIP());

	int i = heap.size(); // heap is in reverse order
	uint count = 0;
	uint CANDIDATE_MAX = heap.size();

	while (	count < CANDIDATE_MAX && i > 0 && heap[--i].first <= threshold	)
		++count;

	// printf("DEBUG stats: %lu in total; ", m_patches.size());

	if ( count == 0 ) return printf("\n"); // TODO: throw, etc.

	uint s = heap.size() - 1 - (rand() % count); // selection

	// printf("%lu in selection; %d under threshold (%f); final selection = %i\n", heap.size(), count, threshold, heap[s].second);

	return heap[s].second;
}

template<uint N, typename T>
template<class U>
void TexSynth::TexSynther<N, T>::extendTextureIn(cimg_library::CImg<U> & _img, cimg_library:: CImg<U> & _msk) const
{
	//printf("extendTextureIn\n");

	Patch ones(255.f);

	CImg<U> debug(_img);

	uint shift = m_patchWidth * 0.6;

	printf("I(%dx%d) and %d\n", _img.width(), _img.height(), m_patchWidth);

	for ( uint y = 0; y < _img.height() - m_patchWidth; y += shift )
		for ( uint x = 0; x < _img.width() - m_patchWidth; x += shift )
		{
			Patch p = m_patches[selectPatchFor(_img, _msk, x, y)];

			// DEBUG
			p.insertInto(debug, x, y);

			// TODO: put the following somewhere more general
			Patch ii(_img, x, y);
			Patch mm(_msk, x, y);
			Patch dd = Patch::diffSqrd(p, ii);
			for ( uint i = 0; i < dd.size(); ++i ) dd[i] = mm[i] != 255.f ? std::numeric_limits<T>::max() : dd[i];
			CImg<float> d(N, N, 1, 3);
			dd.insertInto(d, 0, 0);

			// TODO: do all four directions
			CImg<float> m(N, N, 1, 3, 255);
			if ( x > 0 )
			{
				auto vSeam = findMinSeam<Seams::Downwards>(d);
				vSeam.modifyPatch(m);
			}
			if ( y > 0 )
			{
				auto hSeam = findMinSeam<Seams::Leftwards>(d);
				hSeam.modifyPatch(m);
			}

			//char name[100];
			//sprintf(name, "patch_mask_%d.png", (y / shift) * (_img.width() / shift) + (x / shift));
			//m.save(name);

			// Merge patches
			mm.extractFrom(m, 0, 0);
			p = p.blendedWith(ii, mm);
			p.insertInto(_img, x, y);

			ones.insertInto(_msk, x, y);
		}

	debug.save("debug.png");
}


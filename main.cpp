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

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <limits>

#include "CImg.h"
using cimg_library::CImg;

#include "patch.h"
#include "texSynth.h"
#include "minSeam.h"
using namespace TexSynth;

static uint const g_N = 22;

typedef TexSynther<g_N>::Patch Patch;

int main()
{
	std::string filestem = "mesh";
	std::string inputFilename = filestem + ".png";
	std::string outputFilename = filestem + "_out.png";

	CImg<float> image(inputFilename.c_str());

	/*
	Patch pp(image, 0, 0);
	Patch qq(image, 100, 100);
	Patch dd = abs(diff(pp,qq));
	CImg<float> II(g_N, g_N, 1, 1);
	dd.insertInto(II, 0, 0);
	findMinVerticalSeam()
	*/

#if 1

	Patch tmp;
	if ( tmp.test() ) printf("Patch: all good.\n");

	CImg<float> imgOut(image);
	imgOut.fill(0);

	for ( int y = 0; y < image.height(); y += g_N )
	{
		for ( int x = 0; x < image.width(); x += g_N )
		{
			Patch pp(image, x, y);
			pp.insertInto(imgOut, x, y);
		}
	}

	imgOut.save("test.png");

#else

	std::vector<CImg<float> > imagePyramid;

	CImg<float> imgInter(image);
	while ( std::min(imgInter.height(), imgInter.width()) > 50 )
	{
		printf("Adding I(%dx%d)\n", imgInter.width(), imgInter.height());
		imagePyramid.push_back(imgInter);
		imgInter.resize_halfXY();
	}

	imgInter = imagePyramid.back();
	CImg<float> imgOut(imgInter.width(), imgInter.height(), 1, 3, 0);

	int i = 0;
	for ( auto imgIt = imagePyramid.rbegin(); imgIt != imagePyramid.rend(); ++imgIt, ++i )
	{
		TexSynther<g_N>::PatchLibrary patches;

		imgInter = *imgIt;

		printf("Extracting patches for level %u/%d (%dx%d)... ", i + 1, imagePyramid.size(), imgInter.width(), imgInter.height());

		for ( uint y = 0; y < imgInter.height() - g_N; ++y )
			for ( uint x = 0; x < imgInter.width() - g_N; ++x )
				patches.push_back(Patch(imgInter, x, y));

		printf("done.\n");

		float alpha = 0.8; // weight of actual image in error blending with 'seed'
		CImg<float> M(imgInter.width(), imgInter.height(), 1, 3, 1.f - alpha);

		printf("Creating texture for level %u/%d (%d)...\n", i + 1, imagePyramid.size(), g_N);

		TexSynther<g_N> ts;
		ts.addPatches(patches);
		ts.extendTextureIn(imgOut, M);

		printf("done.\n");

		char tmpFilename[50];
		sprintf(tmpFilename, "%s_%u.png", filestem.c_str(), i);
		imgOut.save(tmpFilename);

		if ( i < imagePyramid.size() - 1 ) imgOut.resize_doubleXY();
	}

	imgOut.save(outputFilename.c_str());

#endif

	if ( 1 )
	{
		cimg_library::CImgDisplay display0(image, "Original");
		cimg_library::CImgDisplay display1(imgOut, "Reconstructed");

		while ( !display0.is_closed() )
		{
			display0.wait(1000);
		}
	}


	return 0;
}

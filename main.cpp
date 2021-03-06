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

#include "common.h"
#include "patch.h"
#include "texSynth.h"
#include "circularSeam.h"

using namespace TexSynth;

int main(int argc, char * * argv)
{
	srand(0);

	for ( int i = 0; i < argc; ++i )
		printf("%u/%i: '%s'\n", i + 1, argc, argv[i] );

	CImg<uchar> testUC("leaf.png"); printf("DEBUG: leaf<uchar>(0,0) = <%i, %i, %i>\n", testUC(0,0,0,0), testUC(0,0,0,1), testUC(0,0,0,2));
	CImg<float> testFL("leaf.png"); printf("DEBUG: leaf<float>(0,0) = <%f, %f, %f>\n", testFL(0,0,0,0), testFL(0,0,0,1), testFL(0,0,0,2));

	std::string filestem = argc > 1 ? argv[1] : "grass";
	std::string inputFilename = filestem + ".png";
	std::string outputFilename = filestem + "_out.png";

	CImg<uchar> image(inputFilename.c_str());

	static uint const N = 30;
	typedef TexSynther<N>::Patch Patch;

/*
	Table<float> tab(N, N);
	uint i = 0;
	uint j = 0;

	i = 0; tab(i++,   j) =   4; tab(i++, j) = 225; tab(i++, j) = 121; tab(i++, j) =   9; tab(i++, j) =   9; tab(i++, j) =   1; tab(i++, j) =  16; tab(i++, j) =   1; tab(i++, j) =   4; tab(i++, j) =   4; tab(i++, j) =   1; tab(i++, j) =  16; tab(i++, j) =   1; tab(i++, j) =  64; tab(i++, j) =   1; tab(i++, j) =   0; tab(i++, j) =  64; tab(i++, j) =  49; tab(i++, j) =  49; tab(i++, j) =  64; tab(i++, j) = 361; tab(i++, j) = 289; tab(i++, j) =   1; tab(i++, j) = 361; tab(i++, j) = 529; tab(i++, j) = 400; tab(i++, j) = 289; tab(i++, j) = 576; tab(i++, j) = 625; tab(i++, j) = 625;
	i = 0; tab(i++, ++j) = 121; tab(i++, j) =  81; tab(i++, j) =  72.25; tab(i++, j) =  36; tab(i++, j) =   4; tab(i++, j) =   4; tab(i++, j) =   4; tab(i++, j) =   0; tab(i++, j) =   1; tab(i++, j) =   1; tab(i++, j) =   9; tab(i++, j) =  49; tab(i++, j) =  49; tab(i++, j) =   9; tab(i++, j) =   4; tab(i++, j) =   1; tab(i++, j) =  64; tab(i++, j) =  36; tab(i++, j) =  25; tab(i++, j) =   1; tab(i++, j) =   4; tab(i++, j) = 169; tab(i++, j) = 1024; tab(i++, j) = 841; tab(i++, j) = 676; tab(i++, j) = 361; tab(i++, j) = 196; tab(i++, j) = 256; tab(i++, j) = 144; tab(i++, j) =  25;

	Table<float> mtab(tab);

	for ( uint i = 0; i < tab.size(); ++i )
		if ( tab[i] == -1 )
			mtab[i] = 0;
		else
			mtab[i] = 256;

	std::cout<<"."<<std::endl;

	CircSeams::findMin(tab, mtab);
*/


#if 0

	Patch pp(image, 0, 0);
	Patch qq(image, 100, 100);
	Patch dd = Patch::diffSqrd(pp, qq);
	CImg<float> d(N, N, 1, 3);
	dd.insertInto(d, 0, 0);
	printf("#### DOWNWARDS ####\n");
	auto down = findMinSeam<Seams::Downwards>(d);
	printf("#### UPWARDS ####\n");
	auto up = findMinSeam<Seams::Upwards>(d);
	printf("#### RIGHTWARDS ####\n");
	auto right = findMinSeam<Seams::Rightwards>(d);
	printf("#### LEFTWARDS ####\n");
	auto left = findMinSeam<Seams::Leftwards>(d);

	CImg<uchar> mask(N, N, 1, 1);

	mask.fill(255);
	down.modifyPatch(mask);
	mask.save("down.png");

	mask.fill(255);
	up.modifyPatch(mask);
	mask.save("up.png");

	mask.fill(255);
	right.modifyPatch(mask);
	mask.save("right.png");

	mask.fill(255);
	left.modifyPatch(mask);
	mask.save("left.png");

#elif 0

	Patch mm(255.f);
	Patch pp(image, 0, 0);
	uint qqx = 4; // rand() % image.width() - N;
	uint qqy = 0; // rand() % image.width() - N;
	Patch qq(image, qqx, qqy);

	CImg<float> tmp(N, N, 1, 3);
	pp.insertInto(tmp, 0, 0);
	tmp.save("pp1.png");
	qq.insertInto(tmp, 0, 0);
	tmp.save("qq.png");
	mm.insertInto(tmp, 0, 0);
	tmp.save("mm1.png");

	CImg<float> img(1.5 * N, N, 1, 3, 0.f);
	CImg<float> msk(1.5 * N, N, 1, 3, 0.f);
	pp.insertInto(img, 0, 0);
	mm.insertInto(msk, 0, 0);

	img.save("img.png");
	msk.save("msk.png");

	pp.extractFrom(img, N / 2.f - 1, 0);
	mm.extractFrom(msk, N / 2.f - 1, 0);

	pp.insertInto(tmp, 0, 0);
	tmp.save("pp2.png");
	mm.insertInto(tmp, 0, 0);
	tmp.save("mm2.png");

	Patch dd = Patch::diffSqrd(qq, pp);

	dd.insertInto(tmp, 0, 0);
	tmp.save("dd1.png");

	for ( uint i = 0; i < dd.size(); ++i ) dd[i] = mm[i] < 255.f ? std::numeric_limits<float>::max() : dd[i];

	dd.insertInto(tmp, 0, 0);
	tmp.save("dd2.png");

	CImg<float> d(N, N, 1, 3);
	dd.insertInto(d, 0, 0);
	// TODO: do all four directions
	auto vSeam = findMinSeam<Seams::Downwards>(d);
	d.fill(255.f);
	vSeam.modifyPatch(d);
	d.save("mask.png");

	dd.extractFrom(d, 0, 0);

	// Merge patches
	qq = qq.blendedWith(pp, dd);
	qq.insertInto(img, N / 2.f - 1, 0);

	img.save("seam_test.png");

#elif 1
	// *****

	int newWidth = N * (image.width() / N) * 3;
	int newHeight = N * (image.height() / N) * 3;
	CImg<float> imgOut(newWidth, newHeight, 1, 3, 0.f);

	float alpha = 1.f;
	CImg<float> M(imgOut.width(), imgOut.height(), 1, 3, 1.f - alpha);

	TexSynther<N>::PatchLibrary patches;

	printf("Extracting patches... ");

	for ( uint y = 0; y <= image.height() - N; ++y )
		for ( uint x = 0; x <= image.width() - N; ++x )
			patches.push_back(Patch(image, x, y));

	printf("done.\n");

	printf("Creating texture... ");

	TexSynther<N> ts;
	ts.addPatches(patches);
	ts.extendTextureIn(imgOut, M);

	printf("done.\n");

	imgOut.save(outputFilename.c_str());

#elif 0

	CImg<float> imgOut(image.width(), image.height(), 1, 3, 0);
	std::vector<CImg<float> > imagePyramid;

	uint minDim = 50; // std::min(image.height(), image.width());

	CImg<float> imgInter(image);
	while ( std::min(imgInter.height(), imgInter.width()) > minDim )
	{
		printf("Adding I(%dx%d)\n", imgInter.width(), imgInter.height());
		imagePyramid.push_back(imgInter);
		imgInter.resize_halfXY();
	}

	imgInter = imagePyramid.back();

	int i = 0;
	for ( auto imgIt = imagePyramid.rbegin(); imgIt != imagePyramid.rend(); ++imgIt, ++i )
	{
		TexSynther<N>::PatchLibrary patches;

		imgInter = *imgIt;

		printf("Extracting patches for level %u/%d (%dx%d)... ", i + 1, imagePyramid.size(), imgInter.width(), imgInter.height());

		for ( uint y = 0; y <= imgInter.height() - N; ++y )
			for ( uint x = 0; x <= imgInter.width() - N; ++x )
				patches.push_back(Patch(imgInter, x, y));

		printf("done.\n");

		float alpha = 1.f; // weight of actual image in error blending with 'seed'
		CImg<float> M(imgInter.width(), imgInter.height(), 1, 3, 1.f - alpha);

		printf("Creating texture for level %u/%d (%d)...\n", i + 1, imagePyramid.size(), N);

		TexSynther<N> ts;
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

#if DISPLAY
	if ( 0 )
	{
		cimg_library::CImgDisplay display0(image, "Original");
		//cimg_library::CImgDisplay display1(imgOut, "Reconstructed");

		while ( !display0.is_closed() )
		{
			display0.wait(1000);
		}
	}
#endif

	return 0;
}


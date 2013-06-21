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

	static uint const N = 15;
	typedef TexSynther<N>::Patch Patch;

	/*
	Table<float> tab(N, N);
	uint i = 0;
	uint j = 0;

		  tab(i++,   j) =    197; tab(i++, j) = 3291; tab(i++, j) = 1328; tab(i++, j) =  362; tab(i++, j) = 1726; tab(i++, j) =  694; tab(i++, j) =  883; tab(i++, j) =   47; tab(i++, j) =   20; tab(i++, j) = 1372; tab(i++, j) =  250; tab(i++, j) =  108; tab(i++, j) =   62; tab(i++, j) = 2818; tab(i++, j) =   81;
   i = 0; tab(i++, ++j) =   1210; tab(i++, j) =   77; tab(i++, j) =  190; tab(i++, j) =   27; tab(i++, j) =  347; tab(i++, j) =   28; tab(i++, j) =  428; tab(i++, j) =   79; tab(i++, j) =  694; tab(i++, j) = 1205; tab(i++, j) =  268; tab(i++, j) =   24; tab(i++, j) =  182; tab(i++, j) = 2038; tab(i++, j) =  373;
   i = 0; tab(i++, ++j) =    972; tab(i++, j) =  121; tab(i++, j) =  141; tab(i++, j) =  867; tab(i++, j) =  662; tab(i++, j) =   41; tab(i++, j) =  713; tab(i++, j) =  170; tab(i++, j) =    3; tab(i++, j) =  15; tab(i++, j) =  11; tab(i++, j) =  35; tab(i++, j) = 255; tab(i++, j) = 1382; tab(i++, j) =  105;
   i = 0; tab(i++, ++j) =    754; tab(i++, j) =  846; tab(i++, j) =  578; tab(i++, j) =   48; tab(i++, j) =    6; tab(i++, j) =    4; tab(i++, j) = 842; tab(i++, j) =  130; tab(i++, j) = 408; tab(i++, j) =   20; tab(i++, j) =   47; tab(i++, j) =   21; tab(i++, j) =  200; tab(i++, j) =  162; tab(i++, j) =  45;
   i = 0; tab(i++, ++j) =    415; tab(i++, j) =  267; tab(i++, j) =  246; tab(i++, j) =  581; tab(i++, j) =   67; tab(i++, j) =  17; tab(i++, j) =  415; tab(i++, j) = 695; tab(i++, j) =  737; tab(i++, j) =  543; tab(i++, j) =   33; tab(i++, j) =  786; tab(i++, j) =  103; tab(i++, j) =  492; tab(i++, j) =   18;
   i = 0; tab(i++, ++j) =    442; tab(i++, j) =  281; tab(i++, j) =   29; tab(i++, j) =   36; tab(i++, j) =  108; tab(i++, j) = 287; tab(i++, j) =  289; tab(i++, j) = 1553; tab(i++, j) = 1259; tab(i++, j) =  490; tab(i++, j) =  107; tab(i++, j) =  190; tab(i++, j) =  528; tab(i++, j) =  440; tab(i++, j) = 1409;
   i = 0; tab(i++, ++j) =    430; tab(i++, j) =   13; tab(i++, j) =  462; tab(i++, j) =  443; tab(i++, j) =  499; tab(i++, j) =   41; tab(i++, j) = -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    362; tab(i++, j) =  217; tab(i++, j) = 1449; tab(i++, j) =  983; tab(i++, j) =  489; tab(i++, j) =  329; tab(i++, j) = -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    489; tab(i++, j) =   36; tab(i++, j) = 1321; tab(i++, j) =  144; tab(i++, j) = 1081; tab(i++, j) =  175; tab(i++, j) = -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    789; tab(i++, j) =  681; tab(i++, j) =   47; tab(i++, j) =  498; tab(i++, j) =  221; tab(i++, j) =   51; tab(i++, j) = -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    731; tab(i++, j) = 1480; tab(i++, j) =   11; tab(i++, j) =  176; tab(i++, j) =   18; tab(i++, j) = 199; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    207; tab(i++, j) =  520; tab(i++, j) =   50; tab(i++, j) =   52; tab(i++, j) = 233; tab(i++, j) =  626; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    492; tab(i++, j) =  351; tab(i++, j) =   91; tab(i++, j) =   65; tab(i++, j) = 745; tab(i++, j) = 1000; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =    350; tab(i++, j) = 1161; tab(i++, j) = 1125; tab(i++, j) =  206; tab(i++, j) =  646; tab(i++, j) = 2588; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;
   i = 0; tab(i++, ++j) =     69; tab(i++, j) =  109; tab(i++, j) =  678; tab(i++, j) =  205; tab(i++, j) =  310; tab(i++, j) = 1960; tab(i++, j) = -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1; tab(i++, j) =  -1;


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

	int newWidth = 60; // N * (image.width() / N) * 3;
	int newHeight = 60; // N * (image.height() / N) * 3;
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


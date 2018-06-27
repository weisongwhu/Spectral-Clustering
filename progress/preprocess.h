//#include "stdafx.h"
#include "cv.h"
#include "highgui.h"
#include "highgui.hpp"
#include "core.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace cv;
using namespace std;

#ifndef _PREPRO
#define _PREPRO

#define HISTC 256 

void diffusion(Mat& srcimg, Mat& dstimg, int nite, int K, float dt=0.2, float sigmma=0);//niteµü´ú´ÎÊý£¬KÎªãÐÖµ
void calGrad(Mat& src,Mat& dN, Mat&  dS, Mat& dW, Mat& dE);

void sigmoid(Mat& src);

void Homomorphic(Mat& src);

void getLabelImgContours(Mat& origImg, Mat& labelImg);

void mor(IplImage* src);
// class preprocess
// {
// public:
// 	preprocess(void);
// 	
// 	
// 	~preprocess(void);
// }

#endif
// #pragma once
// 
// #include "cv.h"
// #include "highgui.h"
// #include <algorithm>
// #include <vector>
// #include <stdio.h>
// #include <ctype.h>
// using namespace cv;
// // #ifndef _PREPRO
// // #define _PREPRO
// // 
// // #define HISTC 256 
// 
// 
// 
// /*#endif*/
// class preprocess
// {
// public:
// 	preprocess(void);
// 	void Homomorphic(Mat& src);
// 	~preprocess(void);
// };
// 

#pragma once
#include <opencv2/core/core.hpp> 
#include<opencv2\opencv.hpp>
using namespace cv;

int pixelJudge(Mat &img, int ix, int jy, int* erosionModel, int modelSize = 3);
int ErosionDeal(Mat &img, int* erosionModel, Mat &imgDelt, int modelSize = 3);
int iterativeErosion(Mat& img, int* erosionModel, Mat& imgDelt, double RatioToStop = 0.1);
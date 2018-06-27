
#include "erosionProcess.h"
int ik8[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
int jk8[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
int pixelJudge(Mat &img, int ix, int jy, int* erosionModel, int modelSize)
{
	int x = ix - int(modelSize / 2);
	int y = jy - int(modelSize / 2);
	for (int i = 0; i < modelSize; i++)
	{
		for (int j = 0; j < modelSize; j++)
		{
			if ((x + i >= 0) && (y + j >= 0) && (x + i < img.rows) && (y + j <= img.cols))
			{
				if (erosionModel[i*modelSize + j] && img.at<uchar>(x + i, y + j)!=255)
					return 0;
			}
			else
				return 0;
		}
	}
	return 1;
}
int ErosionDeal(Mat &img, int* erosionModel, Mat &imgDelt, int modelSize)
{
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			if (pixelJudge(img, i, j, erosionModel))
			{
				imgDelt.at<uchar>(i, j) = 255;
			}
		}
	}
	return 0;
}
int iterativeErosion(Mat& img, int* erosionModel, Mat& imgDelt, double RatioToStop)
{
	int pixelSum = 0;
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			if (img.at<uchar>(i, j) == 255)
				pixelSum++;
		}
	}
	int ErosionFlag = 0;
	int pixelCount = pixelSum;
	Mat img0 = imgDelt.clone();
	Mat img1 = img.clone();
	Mat img2;
	//while (pixelCount>pixelSum*RatioToStop&&pixelCount>20)
	while (pixelCount!=0)
	{
		ErosionFlag = 1;
		imgDelt.release();
		imgDelt = img0.clone();
		ErosionDeal(img1, erosionModel, imgDelt);
		pixelCount = 0;
		for (int i = 0; i < imgDelt.rows; i++)
		{
			for (int j = 0; j < imgDelt.cols; j++)
			{
				if (imgDelt.at<uchar>(i, j) == 255)
					pixelCount++;
			}
		}
		img2.release();
		img2 = img1.clone();
		img1.release();
		img1 = imgDelt.clone();

	}
	imgDelt.release(); imgDelt = img2.clone();
	img1.release();
	img2.release();
	//if (ErosionFlag == 0) { imgDelt.release(); imgDelt = img1.clone(); }

	return 0;
}

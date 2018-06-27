#include "stdio.h"

#include "iostream"
#include "fstream"
#include<vector>
#include "cv.h"
#include "highgui.h"
#include "highgui.hpp"

#include "core.hpp"
#include "SLIC.h"
#include "extractTexture.h"
#include "autoClaster.h"
#include "pixelFunction.h"
#include "preprocess.h"
#define M_ALPHA 25.5
#define M_BETA 5
int iknB[8] = { -1, -1, -1, 0, 0, 1, 1, 1 };
int jknB[8] = { -1, 0, 1, -1, 1, -1, 0, 1 };
using namespace std;
using namespace cv;
int caculateContourNum1(Mat &img1, int di, int dj, int minX, int minY, int maxX, int maxY);
int caculateConPixelNum(vector<Point> &contourToCaculate, int sizeX, int sizeY);
int extureDestiContour(vector<int>& supixelInfmation, Mat& img, int K, int ROINum);
int findAllContour(vector<Point> &finalContourline, vector<vector<int>>& finalContourFlag, vector<vector<int>> &contourflag, vector<int>& boundryFlag, int di, int dj, int k);
int updatePQ(vector<int> &weight, vector<vector<int>> &positionQueue, int pXY, int od);
int showImg(vector<int>& supixelInfmation, Mat& img);
int showImg(vector<int>& supixelInfmation, Mat& img, vector<vector<int>>& neighFlag, vector<int>& superpixelName);
int clusterBasedSuperpixel(string pic_name, string sPathName, int ROI_pixel_num, int kn, string si, Mat &imgs);
int superpixelNum(vector<int> superpixelInfmation, Mat &imgs, int superpixelSize);
int main()
{

	string sPathName = "E:\\人工轮廓\\";
	string sOriginalPath = "E:\\图片\\原图\\";
	for (int i = 3; i <= 42; i++)
	{
		if (i==9||i == 10 || i == 16 || i == 19 || i == 35 || i == 39) { continue; }
		char buf0[10];
		sprintf_s(buf0, "%d", i);
		string si = buf0;
		string pic_name = sPathName + si + ".bmp";
		//string pic_name = path + string(ch[i]) + "_2.png";
		Mat    img = imread(pic_name, CV_LOAD_IMAGE_UNCHANGED);
		Mat imgs(img.rows, img.cols, CV_8UC1);
		int ROI_Pixel_num = caculateRoiNum(img, si, sPathName, imgs);
		fstream fst((sPathName + si + ".txt").c_str(), fstream::out);
		fst << ROI_Pixel_num << endl;
		fst.close();
		string origin_pic_name = sOriginalPath + si + ".jpg";
		//string origin_pic_name = path + string(ch[i]) + "_1.png";
		double pixel_num;
		if (ROI_Pixel_num > 10000) { pixel_num = 500; }
		else { pixel_num = 100; }

		for (int kn =10; kn < (double)ROI_Pixel_num / pixel_num; kn = kn + 2)
		{
			clusterBasedSuperpixel(origin_pic_name, sPathName, ROI_Pixel_num, kn, si, imgs);

		}
	}

	return 0;
}
int clusterBasedSuperpixel(string pic_name, string sPathName, int ROI_pixel_num, int kn, string si, Mat &imgs)
{
	char buf[10];
	sprintf_s(buf, "%d", kn);
	string skn = buf;
	Mat   grayimg = imread(pic_name, CV_LOAD_IMAGE_GRAYSCALE);
	Mat grayimg1 = grayimg.clone();
	//diffusion(grayimg, grayimg, 10, 10);
	//Homomorphic(grayimg);

	vector<int> supixelInfmation(grayimg.rows*grayimg.cols, 0);
	int supNum = kn* grayimg.rows*grayimg.cols / ROI_pixel_num;

	SLIC slic;
	int compactness = 2;
	slic.img = grayimg.clone();
	int* labels = new int[grayimg.rows*grayimg.cols];
	slic.PerformSuperPixel(supNum,compactness , labels);
	//mainImrg(grayimg, 0.1, imrgNum, 0.001, supixelInfmation);
	Mat suimg = grayimg.clone();

	vector<int> superpixelName;
	vector<vector<int>> coordinateOfSuperpixel;
	for (int i = 0; i < supixelInfmation.size(); i++)
	{
		supixelInfmation[i] = labels[i];
		updatePQ(superpixelName, coordinateOfSuperpixel, i, supixelInfmation[i]);
	}
	delete[]labels;
	int superpixelSize = superpixelName.size();
	int numOfsuperpixel = superpixelNum(supixelInfmation, imgs, superpixelSize);
	fstream fst((sPathName + si + ".txt").c_str(), fstream::app);
	fst <<skn<<": "<<numOfsuperpixel << endl;
	fst.close();
	vector<vector<double>> wOfNJW(superpixelName.size(), vector<double>(superpixelName.size(), 0));
	vector<vector<int>> neighFlag(superpixelName.size(), vector<int>(superpixelName.size(), 0));
	showImg(supixelInfmation, suimg,neighFlag, superpixelName);
	imwrite(sPathName + si + "_superpixel_" + skn + ".bmp", suimg);
	getAffinityMatrix(grayimg, superpixelName, coordinateOfSuperpixel, supixelInfmation, M_ALPHA, M_BETA, wOfNJW,neighFlag);
	vector<vector<int>> ClasterInfmation;
	mainAClaster(wOfNJW, ClasterInfmation);
	vector<int> supixelInfmation1(grayimg.rows*grayimg.cols, 0);
	for (int i = 0; i < ClasterInfmation.size(); i++)
	{
		for (int j = 0; j < ClasterInfmation[i].size(); j++)
		{
			int clusterFlag = ClasterInfmation[i][j];
			for (int k = 0; k < coordinateOfSuperpixel[clusterFlag].size(); k++)
			{
				supixelInfmation1[coordinateOfSuperpixel[clusterFlag][k]] = i;
			}
		}
	}
	if (ClasterInfmation.size() == 0)return 0;
	coordinateOfSuperpixel.clear();
	superpixelName.clear();
	extureDestiContour(supixelInfmation1, grayimg1, ClasterInfmation.size(), ROI_pixel_num);
	imwrite(sPathName + si + "_final_cluster_" + skn + ".bmp", grayimg1);
	showImg(supixelInfmation1, grayimg1);
	imwrite(sPathName + si + "_cluster_" + skn + ".bmp", grayimg1);
	supixelInfmation.clear();
	return 0;
}
int extureDestiContour(vector<int>& supixelInfmation, Mat& img, int K, int ROINum)
{
	int sizeX = img.rows;
	int sizeY = img.cols;
	vector<vector<int>> clusterFlag(sizeX, vector<int>(sizeY, 0));
	vector<vector<int>> contourflag(sizeX, vector <int>(sizeY, 1000));
	for (int i = 0; i < supixelInfmation.size(); i++)
	{
		int nX = i / sizeY;
		int nY = i %sizeY;
		clusterFlag[nX][nY] = supixelInfmation[i];
	}
	//	vector<vector<Point>> contourPixel(K, vector<Point>(0,Point(0,0)));
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				int ni = i + iknB[k];
				int nj = j + jknB[k];
				if (ni >= 0 && ni<sizeX&&nj >= 0 && nj<sizeY)
				{
					if (clusterFlag[ni][nj] != clusterFlag[i][j])
					{
						contourflag[i][j] = clusterFlag[i][j];
					}
				}
			}
		}
	}
	clusterFlag.clear();

	vector<vector<int>> finalContourFlag(sizeX, vector<int>(sizeY, 0));
	vector<vector<Point>> finalContour;
	vector<int> boundryFlag;
	int finalContourCnt = 0;
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{

			if (contourflag[i][j] >= 0 && contourflag[i][j] < K&&finalContourFlag[i][j] == 0)
			{
				vector<Point> line1;
				boundryFlag.push_back(0);
				findAllContour(line1, finalContourFlag, contourflag, boundryFlag, i, j, finalContourCnt);

				finalContour.push_back(line1);
				finalContourCnt++;
				line1.clear();
			}
		}
	}
	vector<int> pixelNumInContour(finalContourCnt, 0);
	for (int i = 0; i < finalContourCnt; i++)
	{
		if (boundryFlag[i] == 0)
		{
			pixelNumInContour[i] = caculateConPixelNum(finalContour[i], sizeX, sizeY);
		}

	}
	if (finalContour.size() == 0)
	{
		return 0;
	}
	finalContourFlag.clear();
	int destinationFlag = 0;
	int minContourLength = 100000;
	for (int i = 0; i < finalContourCnt; i++)
	{
		int NumChaZhi = abs(pixelNumInContour[i] - ROINum);
		if (boundryFlag[i] == 0 && minContourLength>NumChaZhi)
		{
			minContourLength = NumChaZhi;
			destinationFlag = i;
		}
	}


	for (int i = 0; i < finalContour[destinationFlag].size(); i++)
	{
		img.at<uchar>(finalContour[destinationFlag][i].x, finalContour[destinationFlag][i].y) = 255;
	}
	pixelNumInContour.clear();
	boundryFlag.clear();
	finalContour.clear();
	contourflag.clear();
	return 0;
}
int findAllContour(vector<Point> &finalContourline, vector<vector<int>>& finalContourFlag, vector<vector<int>> &contourflag, vector<int>& boundryFlag, int di, int dj, int k)
{
	int sizeX = contourflag.size();
	int sizeY = contourflag[0].size();
	/*
	if (di == 0 || di == (sizeX - 1) || dj == 0||dj==(sizeY-1))
	{
	/////边界轮廓
	boundryFlag[k] = 1;
	}*/
	finalContourline.push_back(Point(di, dj));
	finalContourFlag[di][dj] = 1;
	for (int i = 0; i < 8; i++)
	{
		int ni = iknB[i] + di;
		int nj = jknB[i] + dj;
		if (ni >= 0 && ni<sizeX&&nj >= 0 && nj<sizeY)
		{
			if (contourflag[di][dj] == contourflag[ni][nj] && finalContourFlag[ni][nj] == 0)
			{
				findAllContour(finalContourline, finalContourFlag, contourflag, boundryFlag, ni, nj, k);
			}
		}
		if (ni<0 || nj<0 || ni >= sizeX || nj >= sizeY)
		{
			boundryFlag[k] = 1;
		}
	}
	return 0;
}

int showImg(vector<int>& supixelInfmation, Mat& img, vector<vector<int>>& neighFlag, vector<int>& superpixelName)
{
	int maxsupixelName = 0;
	for (int i = 0; i < superpixelName.size(); i++)
	{
		if (maxsupixelName<superpixelName[i])
		{
			maxsupixelName = superpixelName[i];
		}
	}
	vector<int> superpixelNameFlag(maxsupixelName + 1, 0);
	for (int i = 0; i < superpixelName.size(); i++)
	{
		superpixelNameFlag[superpixelName[i]] = i;
	}
	int sizeX = img.rows;
	int sizeY = img.cols;
	vector<vector<int>> clusterFlag(sizeX, vector<int>(sizeY, 0));
	for (int i = 0; i < supixelInfmation.size(); i++)
	{
		int nX = i / sizeY;
		int nY = i %sizeY;
		clusterFlag[nX][nY] = supixelInfmation[i];
	}
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				int ik = i + iknB[k];
				int jk = j + jknB[k];
				if (ik >= 0 && ik < sizeX&&jk >= 0 && jk < sizeY)
				{
					if (clusterFlag[ik][jk] != clusterFlag[i][j])
					{
						img.at<uchar>(i, j) = 255;
						neighFlag[superpixelNameFlag[clusterFlag[i][j]]][superpixelNameFlag[clusterFlag[ik][jk]]] = -1;

					}
				}
			}
			/*
			if ((i - 1) >= 0 && (j - 1) >= 0)
			{
			if (clusterFlag[i - 1][j - 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; wOfNJW[clusterFlag[i - 1][j - 1]][clusterFlag[i][j]] = -1; }
			if (clusterFlag[i][j - 1] != clusterFlag[i][j]) {img.at<uchar>(i, j) = 255; wOfNJW[clusterFlag[i][j - 1]][clusterFlag[i][j]] = -1;}
			if (clusterFlag[i - 1][j ] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i-1)>=0&&(j+1)<sizeY)
			{
			if (clusterFlag[i - 1][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			if (clusterFlag[i][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i+1)<sizeX&&(j-1)>=0)
			{
			if (clusterFlag[i + 1][j - 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			if (clusterFlag[i + 1][j ] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i + 1) <= sizeX && (j + 1) < sizeY)
			{
			if (clusterFlag[i + 1][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}*/
		}
	}
	return 0;
}
int showImg(vector<int>& supixelInfmation, Mat& img)
{
	int sizeX = img.rows;
	int sizeY = img.cols;
	vector<vector<int>> clusterFlag(sizeX, vector<int>(sizeY, 0));
	for (int i = 0; i < supixelInfmation.size(); i++)
	{
		int nX = i / sizeY;
		int nY = i %sizeY;
		clusterFlag[nX][nY] = supixelInfmation[i];
	}
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			if ((i - 1) >= 0 && (j - 1) >= 0)
			{
				if (clusterFlag[i - 1][j - 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
				if (clusterFlag[i][j - 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
				if (clusterFlag[i - 1][j] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i - 1) >= 0 && (j + 1) < sizeY)
			{
				if (clusterFlag[i - 1][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
				if (clusterFlag[i][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i + 1) < sizeX && (j - 1) >= 0)
			{
				if (clusterFlag[i + 1][j - 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
				if (clusterFlag[i + 1][j] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
			else if ((i + 1) <= sizeX && (j + 1) < sizeY)
			{
				if (clusterFlag[i + 1][j + 1] != clusterFlag[i][j]) { img.at<uchar>(i, j) = 255; }
			}
		}
	}
	return 0;
}
int superpixelNum(vector<int> superpixelInfmation, Mat &imgs,int superpixelSize)
{

	int sizeX = imgs.rows;
	int sizeY = imgs.cols;
	double NumOfEachSuperpixel = (double)sizeX*sizeY / superpixelSize;
	vector<int> NumStore;
	for (int i = 0; i < sizeX;i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			if (imgs.at<uchar>(i, j) == CONTOUR_SIGN)
			{
				int ijp = i*sizeY + j;
				int k = 0;
				while (k < NumStore.size() && NumStore[k] != superpixelInfmation[ijp]) { k++; }
				if (k == NumStore.size())
				{
					NumStore.push_back(superpixelInfmation[ijp]);
				}
			}
		}
	}
	vector<vector<int>> pixelStore(NumStore.size(),vector<int>());
//	vector<int>  pixel_num(NumStore.size(), 0);
	vector<int> superInf(sizeX*sizeY, 0);
	int numOfsuperpixel = 0;
	for (int k = 0; k < NumStore.size();k++)
	{
		for (int i = 0; i < sizeX; i++)
			for (int j = 0; j < sizeY; j++)
			{
				if (superpixelInfmation[i*sizeY+j]==NumStore[k]&&imgs.at<uchar>(i,j)==CONTOUR_SIGN)
				{
					pixelStore[k].push_back(i*sizeY+j);
//					pixel_num[k]++;
				}
			}
		if (pixelStore[k].size()>=(NumOfEachSuperpixel/4))
		{
			numOfsuperpixel++;
			for(int m = 0; m < pixelStore[k].size(); m++)
			{
				superInf[pixelStore[k][m]] = k;
			}
			
		}
	}
	Mat imgToshow(sizeX, sizeY, CV_8UC1);
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			if (imgs.at<uchar>(i, j) == 255)
			{
				imgToshow.at<uchar>(i, j) = 125;
			}
			else
			{
				imgToshow.at<uchar>(i, j) = 0;
			}
		}
	}
	showImg(superInf, imgToshow);
	imwrite("1.bmp", imgToshow);
	return numOfsuperpixel;
}
int updatePQ(vector<int> &weight, vector<vector<int>> &positionQueue, int pXY, int od)
{
	//////////////////////////////////////////////////////////////////////////
	/////////////////////////////该函数可以跟据某个矩阵中某个值出现的个数od，按
	///////////////////////小到大的顺序排列个数值od，并且能够更新相应的数值对应的坐标PXY
	///////////////////////根据权值按D值大小更新坐标队列
	int x = 0;
	while (x < weight.size() && od != weight[x]) { x++; }
	if (x < weight.size())
	{
		positionQueue[x].push_back(pXY);
	}
	else
	{
		//////按大小找出od在weight中的位置并插入
		if (weight.size() == 0)
		{
			weight.push_back(od);
			vector<int> line1;
			line1.push_back(pXY);
			positionQueue.push_back(line1);
			line1.clear();
		}
		else
		{
			x = 0;
			while (x < weight.size() && od > weight[x]) { x++; }
			weight.insert(weight.begin() + x, od);
			vector<int> line1;
			line1.push_back(pXY);
			positionQueue.insert(positionQueue.begin() + x, line1);
			line1.clear();
		}
	}
	return 0;
}
int caculateConPixelNum(vector<Point> &contourToCaculate, int sizeX, int sizeY)
{
	Mat img1(sizeX, sizeY, CV_8UC1);
	for (int i = 0; i < sizeX; i++)
	{
		for (int j = 0; j < sizeY; j++)
		{
			img1.at<uchar>(i, j) = 0;

		}
	}
	int contourXmin = 10000, contourXmax = 0, contourYmin = 10000, contourYmax = 0;
	for (int i = 0; i < contourToCaculate.size(); i++)
	{
		img1.at<uchar>(contourToCaculate[i].x, contourToCaculate[i].y) = 255;
		if (contourToCaculate[i].x > contourXmax) { contourXmax = contourToCaculate[i].x; }
		if (contourToCaculate[i].x < contourXmin) { contourXmin = contourToCaculate[i].x; }
		if (contourToCaculate[i].y > contourYmax) { contourYmax = contourToCaculate[i].y; }
		if (contourToCaculate[i].y < contourYmin) { contourYmin = contourToCaculate[i].y; }
	}
	//	imwrite("1.bmp", img1);
	caculateContourNum1(img1, contourXmin - 1, contourYmin - 1, contourXmin - 1, contourYmin - 1, contourXmax + 1, contourYmax + 1);


	int ROI_pixel_num = 0;
	int contourLength = 0;
	for (int i = contourXmin - 1; i <= contourXmax + 1; i++)
	{
		for (int j = contourYmin - 1; j < contourYmax + 1; j++)
		{
			if (img1.at<uchar>(i, j) == 0)
			{
				ROI_pixel_num++;
			}
		}
	}
	return ROI_pixel_num;
}
int caculateContourNum1(Mat &img1, int di, int dj, int minX, int minY, int maxX, int maxY)
{
	int ik4[4] = { -1, 0, 0, 1 };
	int jk4[4] = { 0,-1,1,0 };
	img1.at<uchar>(di, dj) = 125;
	for (int i = minX; i <= maxX; i++)
	{
		for (int j = minY; j <= maxY; j++)
		{
			if (img1.at<uchar>(i, j) == 125)
			{
				for (int k = 0; k < 4; k++)
				{
					int neighdi = i + ik4[k];
					int neighdj = j + jk4[k];
					if (neighdi <= maxX&&neighdi >= minX && neighdj <= maxY&&neighdj >= minY)
					{
						if ((img1.at<uchar>(neighdi, neighdj) != 125 && img1.at<uchar>(neighdi, neighdj) != 255))
						{
							img1.at<uchar>(neighdi, neighdj) = 125;
						}
					}
				}
			}
		}
	}
	img1.at<uchar>(maxX, maxY) = 125;
	for (int i = maxX; i >= minX; i--)
	{
		for (int j = maxY; j >= minY; j--)
		{
			if (img1.at<uchar>(i, j) == 125)
			{
				for (int k = 0; k < 4; k++)
				{
					int neighdi = i + ik4[k];
					int neighdj = j + jk4[k];
					if (neighdi <= maxX&&neighdi >= minX && neighdj <= maxY&&neighdj >= minY)
					{
						if ((img1.at<uchar>(neighdi, neighdj) != 125 && img1.at<uchar>(neighdi, neighdj) != 255))
						{
							img1.at<uchar>(neighdi, neighdj) = 125;
						}
					}
				}
			}
		}
	}
	//	imwrite("2.bmp", img1);
	//	if (img1.at<uchar>(di,dj) ==0)
	//	{
	/*
	img1.at<uchar>(di, dj) = 125;
	for (int i = 0; i < 4; i++)
	{
	int neighdi = di + ik4[i];
	int neighdj = dj + jk4[i];
	if (neighdi <= maxX&&neighdi >= minX && neighdj <= maxY&&neighdj >= minY)
	{
	if(img1.at<uchar>(neighdi,neighdj)!=125&& img1.at<uchar>(neighdi, neighdj) != 255)
	caculateContourNum1(img1, neighdi, neighdj, minX, minY, maxX, maxY);
	}
	}*/
	//}
	return 0;
}
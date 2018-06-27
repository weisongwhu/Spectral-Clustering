
#include "pixelFunction.h"

#include<vector>

int caculateRoiNum(Mat& img,string &si,string &sPathName,Mat& imgs)
{
	vector<Point> sd;
	int ROIcontourLength=0;
//	fstream fst((sPathName + "c.txt").c_str(), ofstream::out);
	for (int i = 0; i < img.rows; i++)
	{
		for (int j = 0; j < img.cols; j++)
		{
			//			fst << img.at<Vec3b>(i, j)<<' ';

			if (img.at<Vec3b>(i, j)[1]!=img.at<Vec3b>(i,j)[2]) //== Vec3b(3, 3, 247))
			if (img.at<Vec3b>(i, j)== Vec3b(3, 3, 247))
			{
				sd.push_back(Point(i, j));
				imgs.at<uchar>(i,j) = 255;
				ROIcontourLength++;
			}
			else
			{
				imgs.at<uchar>(i, j) = 0;
			}
		}
		//		fst << endl;
	}
	//	fst.close();
	int contourXmin = 10000, contourXmax = 0, contourYmin = 10000, contourYmax = 0;
	for (int i = 0; i < sd.size(); i++)
	{
		if (sd[i].x > contourXmax) { contourXmax = sd[i].x; }
		if (sd[i].x < contourXmin) { contourXmin = sd[i].x; }
		if (sd[i].y > contourYmax) { contourYmax = sd[i].y; }
		if (sd[i].y < contourYmin) { contourYmin = sd[i].y; }
	}
	caculateContourNum(imgs, contourXmin - 1, contourYmin - 1, contourXmin - 1, contourYmin - 1, contourXmax + 1, contourYmax + 1);
/*
	namedWindow("Picture ROI Contours", CV_WINDOW_AUTOSIZE);
	imshow("Picture ROI Contours", imgs);*/
	
	int ROI_pixel_num = 0;

	for (int i = contourXmin - 1; i <=contourXmax + 1; i++)
	{
		for (int j = contourYmin - 1; j < contourYmax + 1; j++)
		{
			if (imgs.at<uchar>(i,j)==0)
			{
				ROI_pixel_num++;
				imgs.at <uchar>(i,j)= CONTOUR_SIGN;
			}
		}
	}
	sd.clear();
	imwrite(sPathName + si + "_Number" + ".bmp", imgs);
	return ROI_pixel_num;
}
char *allocCat(const char * s1, const char * s2)
{
	char *str = (char*)malloc(strlen(s1) + strlen(s2) + 1);
	if (str == NULL)
	{
		return NULL;
	}
	strcpy(str, s1);
	strcpy(str + strlen(s1), s2);
	*(str + strlen(s1) + strlen(s2)) = 0;
	return str;
}
int caculateContourNum(Mat& img1, int di, int dj, int minX, int minY, int maxX, int maxY)
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
/*
	int ik[4] = { -1, 0, 0, 1 };
	int jk[4] = { 0,-1,1,0 };
	if (img.at<uchar>(di,dj)!=125&&img.at<uchar>(di,dj)!=255)
	{
		img.at<uchar>(di, dj) = 125;
		for (int i = 0; i < 4;i++)
		{
			int neighdi = di + ik[i];
			int neighdj = dj + jk[i];
			if (neighdi <=maxX&&neighdi >= minX && neighdj <=maxY&&neighdj >= minY)
			{
				caculateContourNum(img, neighdi, neighdj,minX,minY,maxX,maxY);
			}
		}
	}*/
	return 0;
}

//////////////////////////////////////////////////////////////////////
#include "extractTexture.h"


int getAffinityMatrix(Mat& img,vector<int>& supixelName,vector<vector<int>>& coordinateOfSupixel,vector<int>& supixelInfmation,
	double alpha_NJW,double beta_NJW,vector<vector<double>>& wOfNJW, vector<vector<int>>& neighFlag)
{
/*
	vector<int> supixelName;
	vector<vector<int>> coordinateOfSupixel;
	for (int i = 0; i < supixelInfmation.size(); i++)
	{
		updatePQ(supixelName, coordinateOfSupixel, i, supixelInfmation[i]);
	}*/
	caculateW1(supixelName, coordinateOfSupixel, wOfNJW, img);
	vector<vector<double>> superpixelVector(supixelName.size(),vector<double>(16,0));
	getSuperpixelVector(img, supixelInfmation, coordinateOfSupixel, superpixelVector, supixelName);
	caculateW(wOfNJW, superpixelVector, alpha_NJW, beta_NJW,neighFlag);
	return 0;
}



int getSuperpixelName(vector<int>& supixelInfmation, vector<int>& supixelName)
{
	for (int i=0;i<supixelInfmation.size();i++)
	{
		int j = 0;
		while (j < supixelName.size() && supixelName[j] != supixelInfmation[i]) { j++; }
		if (j>=supixelInfmation.size())
		{
			supixelName.push_back(supixelInfmation[i]);
		}
	}
	return 0;
}


int caculateW1(vector<int>& supixelName,vector<vector<int>>& coordinateOfSupixel,vector<vector<double>>& wOfNJW,Mat& img)
{
	int sizeX = img.rows;
	int sizeY = img.cols;
	vector<double> averageGray;
	for (int i=0;i<supixelName.size();i++)
	{   
		int sumOfGray=0;
		int countOfI = 0;
		for (int j=0;j<coordinateOfSupixel[i].size();j++)
		{
			int pixelX = coordinateOfSupixel[i][j] / sizeY;
			int pixelY = coordinateOfSupixel[i][j] % sizeY;
			sumOfGray += img.at<uchar>(pixelX, pixelY);
			countOfI++;
		}
		double aveGrayI = (double)sumOfGray / countOfI;
		averageGray.push_back(aveGrayI);
	}
	
	for (int i=0;i<supixelName.size();i++)
	{
		for (int j=0;j<supixelName.size();j++)
		{
			wOfNJW[i][j]=abs(averageGray[i] - averageGray[j]);
		}
	}
	averageGray.clear();
	return 0;
}



int getCoorMatrix(Mat &img,double *coorMatrix,vector<int>& supixelInfmation,vector<int>& pixelOfFlag,int directionFlag,int distanceFlag)
{
//////////////////////////////////////////////////////////////////////////
// coorMatrix为一个大小为256*256的矩阵
//directionFlag        xita
//      0               0
//      1               45
//      2               90
//      3               135
//////////////////////////////////////////////////////////////////////////
	int sumCoorMatrix = 0;
	int sizeY = img.cols;
	int sizeX = img.rows;
	switch (directionFlag)
	{
	case 0: {
		for (int i=0;i<pixelOfFlag.size();i++)
		{

				int pixelX = pixelOfFlag[i] / sizeY;
				int pixelY = pixelOfFlag[i]%sizeY;
				int neighX = pixelX;
				int neighY = pixelY + distanceFlag;
				int neighXY = neighX*sizeY + neighY;
				if (neighX >= 0 && neighX < sizeX&&neighY >= 0 && neighY < sizeY&&supixelInfmation[pixelOfFlag[i]] == supixelInfmation[neighXY])
				{
					coorMatrix[img.at<uchar>(pixelX, pixelY) * 256 + img.at<uchar>(neighX, neighY)] += 1;
					coorMatrix[img.at<uchar>(pixelX, pixelY) + 256 * img.at<uchar>(neighX, neighY)] += 1;
					sumCoorMatrix++;
				}
		}
		break;
	}
	case 1: {
		for (int i = 0; i < pixelOfFlag.size(); i++)
		{

			int pixelX = pixelOfFlag[i] / sizeY;
			int pixelY = pixelOfFlag[i] % sizeY;
			int neighX = pixelX+distanceFlag;
			int neighY = pixelY + distanceFlag;
			int neighXY = neighX*sizeY + neighY;
			if (neighX >= 0 && neighX < sizeX&&neighY >= 0 && neighY < sizeY&&supixelInfmation[pixelOfFlag[i]] == supixelInfmation[neighXY])
			{
				coorMatrix[img.at<uchar>(pixelX, pixelY) * 256 + img.at<uchar>(neighX, neighY)] += 1;
				coorMatrix[img.at<uchar>(pixelX, pixelY) + 256 * img.at<uchar>(neighX, neighY)] += 1;
				sumCoorMatrix++;
			}
		}
		break;
	}
	case 2: {
		for (int i = 0; i < pixelOfFlag.size(); i++)
		{

			int pixelX = pixelOfFlag[i] / sizeY;
			int pixelY = pixelOfFlag[i] % sizeY;
			int neighX = pixelX + distanceFlag;
			int neighY = pixelY;
			int neighXY = neighX*sizeY + neighY;
			if (neighX >= 0 && neighX < sizeX&&neighY >= 0 && neighY < sizeY&&supixelInfmation[pixelOfFlag[i]] == supixelInfmation[neighXY])
			{
				coorMatrix[img.at<uchar>(pixelX, pixelY) * 256 + img.at<uchar>(neighX, neighY)] += 1;
				coorMatrix[img.at<uchar>(pixelX, pixelY) + 256 * img.at<uchar>(neighX, neighY)] += 1;
				sumCoorMatrix++;
			}
		}
		break;
	}
	case 3: {
		for (int i = 0; i < pixelOfFlag.size(); i++)
		{

			int pixelX = pixelOfFlag[i] / sizeY;
			int pixelY = pixelOfFlag[i] % sizeY;
			int neighX = pixelX + distanceFlag;
			int neighY = pixelY-distanceFlag;
			int neighXY = neighX*sizeY + neighY;
			if (neighX >= 0 && neighX < sizeX&&neighY >= 0 && neighY < sizeY&&supixelInfmation[pixelOfFlag[i]] == supixelInfmation[neighXY])
			{
				coorMatrix[img.at<uchar>(pixelX, pixelY) * 256 + img.at<uchar>(neighX, neighY)] += 1;
				coorMatrix[img.at<uchar>(pixelX, pixelY) + 256 * img.at<uchar>(neighX, neighY)] += 1;
				sumCoorMatrix++;
			}
		}
		break;
	}
	default:
		break;
	}
	if (sumCoorMatrix != 0)
	{
		for (int i = 0; i < 256; i++)
		{
			for (int j = 0; j < 256; j++)
			{
				coorMatrix[i * 256 + j] = coorMatrix[i * 256 + j] / (double)sumCoorMatrix;
			}
		}
	}
	return 0;
}

int getSuperpixelVector(Mat& img,vector<int>& supixelInfmation,vector<vector<int>> coordinateOfSupixel,vector<vector<double>>& superpixelVector,vector<int> & supixelName)
{
	for (int i=0;i<supixelName.size();i++)
	{
		//vector<double> line1(16, 0);
		for (int k = 0; k < 4; k++)
		{
			double ASM_sum = 0;
			double ENT_sum = 0;
			double COR_sum = 0;
			double IDM_sum = 0;
			int countNM = 0;
			for (int d = 1; d < 4; d++)
			{
				double ASM = 0;
				double ENT = 0;
				double COR = 0;
				double IDM = 0;
				double coorMatrix[256 * 256] = { 0 };
/*
				for (int i = 0; i < 256; i++)
				{
					for (int j = 0; j < 256; j++)
					{
						coorMatrix[i * 256 + j] = 0;
					}
				}*/
				getCoorMatrix(img, coorMatrix, supixelInfmation,coordinateOfSupixel[i], k,d);
				double averageM = 0;
				double averageN= 0;
				double varianceM = 0;
				double varianceN = 0;
				for (int m=0;m<256;m++)
				{
					//sumM:一行的和
					//sumN:一列的和
					double aveSumM = 0;
					double aveSumN = 0;
					for (int n=0;n<256;n++)
					{
						aveSumM += coorMatrix[m * 256 + n];
						aveSumN += coorMatrix[m + n * 256];
						
					}
					averageM += m*aveSumM;
					averageN += m*aveSumN;

				}
				for (int m = 0; m < 256; m++)
				{
					//sumM:一行的和
					//sumN:一列的和
					double aveSumM = 0;
					double aveSumN = 0;
					for (int n = 0; n < 256; n++)
					{
						aveSumM += coorMatrix[m * 256 + n];
						aveSumN += coorMatrix[m + n * 256];

					}
					varianceM += (m-averageM)*(m-averageM)*aveSumM;
					varianceN += (m-averageN)*(m-averageN)*aveSumN;

				}

				for (int m=0;m<256;m++)
				{
					for (int n = 0; n < 256; n++)
					{
						int positionXY = m * 256 + n;
						ASM += coorMatrix[positionXY] * coorMatrix[positionXY];
						if (coorMatrix[positionXY] > 0)
						{
							ENT += (-coorMatrix[positionXY])*log2(coorMatrix[positionXY]);
						}
						if (varianceN != 0 && varianceM != 0)
						{
							COR += (m - averageM)*(n - averageN)*coorMatrix[positionXY] / sqrt(varianceN*varianceM);
						}
						IDM += coorMatrix[positionXY] / (1 + abs(m - n));
						//IDM = IDM + (m - n)*(m - n)*coorMatrix[positionXY];

					}
				}
				if (varianceN != 0 && varianceM != 0)
				{
					countNM++;
				}
				ASM_sum += ASM;
				ENT_sum += ENT;
				COR_sum += COR;
				IDM_sum += IDM;
				//delete coorMatrix;

			}
			superpixelVector[i][k*4] = (ASM_sum / 3);
			superpixelVector[i][k*4 + 1] = (ENT_sum / 3);
			if (countNM!=0)
			{
				superpixelVector[i][k * 4 + 2] = (COR_sum / countNM);
			}
			superpixelVector[i][k*4 + 3] = (IDM_sum / 3);
		}
/*
		superpixelVector.push_back(line1);
		line1.clear();*/
	}
	return 0;
}

int caculateW(vector<vector<double>>& wOfNJW,vector<vector<double>>& superpixelVector,double alpha_NJW,double beta_NJW, vector<vector<int>>& neighFlag)
{
/*
	if (altoParaFlag == 1)
	{
		vector < vector<double>> vDx(wOfNJW.size(), vector<double>(wOfNJW.size(), 0));
		double sumIJ = 0;
		double sumW = 0;
		for (int i = 0; i < wOfNJW.size(); i++)
		{
			for (int j = 0; j < wOfNJW[i].size(); j++)
			{
				double Dx = 0;
				for (int k = 0; k < 16; k++)
				{
					double fvi = superpixelVector[i][k];
					double fvj = superpixelVector[j][k];
					if ((fvj + fvj) != 0)
					{
						double fvni = fvi * 2 / (fvi + fvj);
						double fvnj = fvj * 2 / (fvi + fvj);
						Dx += (fvni - fvnj)*(fvni - fvnj) / (fvni + fvnj);
					}
					//Dx += ((superpixelVector[i][k] - superpixelVector[j][k]) / (superpixelVector[i][k] + superpixelVector[j][k]))*(superpixelVector[i][k] - superpixelVector[j][k]);
				}
				Dx = Dx*0.5;
				vDx[i][j] = Dx;
				sumIJ += Dx;
				sumW += wOfNJW[i][j];
			}
		}
		double alphaSquare = sumIJ / (wOfNJW.size()*wOfNJW.size());
		double betaSquare = sumW / (wOfNJW.size()*wOfNJW.size());
		for (int i = 0; i < wOfNJW.size(); i++)
		{
			for (int j = 0; j < wOfNJW.size(); j++)
			{

				double DeDx = wOfNJW[i][j]*wOfNJW[i][j] / (alphaSquare)+vDx[i][j] / (betaSquare);
				wOfNJW[i][j] = exp(-DeDx);
			}
		}
		vDx.clear();
	}
	else
	{*/
		for (int i = 0; i < wOfNJW.size(); i++)
		{
			for (int j = 0; j < wOfNJW[i].size(); j++)
			{
				double Dx = 0;
				for (int k = 0; k < 16; k++)
				{
					double fvi = superpixelVector[i][k];
					double fvj = superpixelVector[j][k];
					if ((fvj + fvj) != 0)
					{
						double fvni = fvi * 2 / (fvi + fvj);
						double fvnj = fvj * 2 / (fvi + fvj);
						Dx += (fvni - fvnj)*(fvni - fvnj) / (fvni + fvnj);
					}
					//Dx += ((superpixelVector[i][k] - superpixelVector[j][k])*(superpixelVector[i][k] - superpixelVector[j][k]) / (superpixelVector[i][k] + superpixelVector[j][k]));
				}
				Dx = Dx*0.5;
				double DeDx= wOfNJW[i][j]*wOfNJW[i][j] / (alpha_NJW*alpha_NJW)+Dx / (beta_NJW*beta_NJW);
				if (neighFlag[i][j] == -1)
				{

						wOfNJW[i][j] = exp(-DeDx);

				}
				else
				{
					if (i != j)
					{
						wOfNJW[i][j] = 0;
					}
					else
					{
						wOfNJW[i][j] = 1;
					}
				}
			}
		}
/*	}*/
	return 0;
}

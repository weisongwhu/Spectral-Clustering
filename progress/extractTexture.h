#pragma once
#include <vector>
using namespace std;
#include "cv.h"
#include"highgui.h"
#include <highgui.hpp>
using namespace cv;
int getAffinityMatrix(Mat& img,vector<int>& supixelName,vector<vector<int>>& coordinateOfSupixel,vector<int>& supixelInfmation,
	double alpha_NJW,double bete_NJW,vector<vector<double>>& wOfNJW, vector<vector<int>>& neighFlag);


int getSuperpixelName(vector<int>& supixelInfmation,vector<int>& supixelName);

int caculateW1(vector<int>& supixelName,vector<vector<int>>& coordinateOfSupixel,vector<vector<double>>& wOfNJW,Mat& img);

int getCoorMatrix(Mat &img,double *coorMatrix,vector<int>& supixelInfmation,vector<int>& pixelOfFlag,int directionFlag,int distanceFlag);

int getSuperpixelVector(Mat& img,vector<int>& supixelInfmation,vector<vector<int>> coordinateOfSupixel,vector<vector<double>>& superpixelVector,vector<int> & supixelName);


int caculateW(vector<vector<double>>& wOfNJW,vector<vector<double>>& superpixelVector,double alpha_NJW,double beta_NJW, vector<vector<int>>& neighFlag);

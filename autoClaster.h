#pragma once

#include<string>
#include <vector>
#include "cv.h"
#include "highgui.h"
#include "highgui.hpp"
using namespace std;
using namespace cv;

int mainAClaster(vector<vector<double>>& wOfNJW, vector<vector<int>>& ClasterInfmation);


int caculateX(vector<vector<double>>& wOfNJW, vector<vector<double>>& xOfNJW, vector<vector<double>>& dOfNJW);

int matrixMtpl(vector<vector<double>>& MtplA, vector<vector<double>>& MtplB, vector<vector<double>>& MtplResult);

int evrotMain(vector<vector<double>>& vdVector, vector<vector<int>>& cluster_array, int method_rot,double& J);

double evqual(vector<vector<double>>& vdVector, int dim, int ndata);

int rotate_givens(vector < vector<double>>& vdVector, vector < vector<double>> &vdVectorRot, 
	vector<double>& vtheta, int* ik, int* jk, int angle_num, int dim);

double evqualitygrad(vector<vector<double>>& vdVector, vector<double>& vtheta,
	int* ik, int* jk, int angle_num, int angle_index, int dim, int ndata);

int computeUab(vector<double>& vtheta, vector<vector<double>>& Uab, int a, int b, int* ik, int* jk, int dim);
int computeV(vector<double>& vtheta, vector<vector<double>>& matrixV, int k, int* ik, int* jk, int dim);


int computeA(vector<vector<double>>& vdVector, vector < vector<double>>& matrixA, vector<vector<double>>& matrixU1, vector<vector<double>>& matrixV, vector<vector<double>>& matrixU2);

int autoClaster(vector<vector<double>>& vdVector, vector<vector<int>>& clasterd_array, int dim, int ndata);

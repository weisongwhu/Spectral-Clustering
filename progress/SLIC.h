#pragma once
#include <vector>
#include "cv.h"
#include "highgui.h"
#include <highgui.hpp>
#include <opencv.hpp>
#include <iostream>
//#include "Kmeans.h"
using std::vector;
using namespace cv;
class SLIC
{
public:
	SLIC(void);
	
	~SLIC(void);
	Mat img;
	void PerformSuperPixel(const int& ,const int& ,int*);
	int s;
	int reallabel(int* );
private:
	//============================================================================
	// Get CIELAB image's each l,a,b channel and store in vectors
	//============================================================================
	void GetData(
		double*&					lvec);
	//============================================================================
	// Detect color edges, to help PerturbSeeds()
	//============================================================================
	void GetEdgeMap(
		const double*				lvec,
		vector<double>&				edges);
	//============================================================================
	// Pick seeds for superpixels when number of superpixels is input.
	//============================================================================
	void InitialKSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&	    edges);
	//============================================================================
	// Move the seeds to low gradient positions to avoid putting seeds at region boundaries.
	//============================================================================
	void PerturbSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);
	//============================================================================
	// The main SLIC algorithm for generating superpixels
	//============================================================================
	void PerformIteration(
		vector<double>&				kseedsl,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		int*						klabels,
		const int&					STEP,
		const vector<double>&		edgemag,
		const double&				m);
	//============================================================================
	// Post-processing of SLIC segmentation, to avoid stray labels.
	//============================================================================
	void EnforceLabelConnectivity(
		const int*					labels,
		const int&					width,
		const int&					height,
		int*						nlabels,//input labels that need to be corrected to remove stray labels
		int&						numlabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user

	//============================================================================
	// Find next superpixel label; helper for EnforceLabelConnectivity()
	//============================================================================
	void FindNext(
		const int*					labels,
		int*						nlabels,
		const int&					height,
		const int&					width,
		const int&					h,
		const int&					w,
		const int&					lab,
		int*						xvec,
		int*						yvec,
		int&						count);

	//===========================================================================
	// Adjust labels for my personal use 
	//===========================================================================
	void AdjustLabels(int* labels);
	
private:

	int										m_width;
	int										m_height;
	int										m_depth;
	int										m_NumPixel;
	double*									m_lvec;
};


//===========================================================================
/// DrawContourAroundSegments
///
/// Draw contour and get the boundary pixel, use a flag to indicate whether 
/// to get boundary pixel or not
//===========================================================================

void DrawContoursAroundSegments(
	Mat &				Img,
	const int*				labels,const int & b,const int & g,const int & r);
void Label2IMG(const int* labels, Mat& img);
void F_Gray2Color(Mat& gray_mat, Mat& color_mat);


//	std::vector<int> pImgcx;
//	std::vector<int> pImgcy;

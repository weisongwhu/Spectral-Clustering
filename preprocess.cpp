
#include "preprocess.h"

void diffusion (Mat& srcimg, Mat& dstimg, int nite, int K, float dt, float sigmma)
{
	CvSize size = srcimg.size();//Size of a rectangle or an image
	Mat src1(size,CV_8UC1);//OpenCV C++ n-dimensional dense array class

	Mat dN(size,CV_32FC1);
	Mat dS(size,CV_32FC1);
	Mat dW(size,CV_32FC1);
	Mat dE(size,CV_32FC1);

	Mat cN(size,CV_32FC1);
	Mat cS(size,CV_32FC1);
	Mat cW(size,CV_32FC1);
	Mat cE(size,CV_32FC1);


	Mat src(size,CV_32FC1);
	Mat srctmp(size,CV_32FC1);

	if (srcimg.channels() > 1);
		//cvtColor(srcimg,src1,CV_BGR2GRAY);//Converts an image from one color space to another.
	else
		srcimg.copyTo(src1);
	src1.convertTo(src,CV_32FC1);//convertTo:Converts an array to another datatype with optional scaling

	for(int i=0; i<nite; i++)
	{
/*
		if(sigmma > 0)//.h里面sigmma=0
		{
			src.copyTo(srctmp);		 
			GaussianBlur(src,src,cvSize(3,3),sigmma,sigmma);//GaussianBlur:Smoothes an image using a Gaussian ?lter.
		}*/
		calGrad(src,dN,dS,dW,dE);

		double s = 1.0/double(K);
		dN.convertTo(dN,CV_32FC1,s);
		dS.convertTo(dS,CV_32FC1,s);
		dW.convertTo(dW,CV_32FC1,s);
		dE.convertTo(dE,CV_32FC1,s);

		pow(dN,2,dN);
		pow(dS,2,dS);
		pow(dW,2,dW);
		pow(dE,2,dE);

		add(dN,Mat::ones(size,CV_32FC1),dN);
		add(dS,Mat::ones(size,CV_32FC1),dS);
		add(dW,Mat::ones(size,CV_32FC1),dW);
		add(dE,Mat::ones(size,CV_32FC1),dE);
		divide(1.0,dN,cN);
		divide(1.0,dS,cS);
		divide(1.0,dW,cW);
		divide(1.0,dE,cE);


		if(sigmma > 0)
		{
			//  cvSetData(src,srctmp->imageData,size.width);
			calGrad(srctmp,dN,dS,dW,dE);
			srctmp.copyTo(src);
		}
		else
			calGrad(src,dN,dS,dW,dE);

		multiply(dN,cN,dN);
		multiply(dE,cE,dE);
		multiply(dW,cW,dW);
		multiply(dS,cS,dS);

		// Add the all 4 matrices
		add(dN,dS,dN);
		add(dN,dW,dN);
		add(dN,dE,dN);

		dN.convertTo(dN,CV_32FC1,dt);
		add(src,dN,src);

	}

	src.convertTo(src1,CV_8UC1);

	src1.copyTo(dstimg);

}

void calGrad(Mat &src, Mat& dN, Mat& dS, Mat& dW, Mat& dE)
{
	CvSize size = src.size();

	//dN 上边界
	for(int i=0; i< size.width; i++)
		dN.at<float>(0,i) = 0;
	//dN
	for(int y=1; y<size.height; y++)
		for(int x=0;x< size.width; x++)
			dN.at<float>(y,x) = src.at<float>(y-1,x) - src.at<float>(y,x);
	//dS下边界
	for(int i=0;i< size.width; i++)
		dS.at<float>(size.height-1,i) = 0;
	//dS
	for(int y=0; y<size.height-1; y++)
		for(int x=0; x<size.width; x++)
			dS.at<float>(y,x) =  src.at<float>(y+1,x) -  src.at<float>(y,x);
	//dW左边界
	for(int i=0; i<size.height; i++)
		dW.at<float>(i,0) = 0;
	//dW
	for(int y=0; y<size.height; y++)
		for(int x=1; x<size.width; x++)
			dW.at<float>(y,x) =  src.at<float>(y,x-1) -  src.at<float>(y,x);
	//dE右边界
	for(int i=0; i<size.height; i++)
		dE.at<float>(i,size.width-1) = 0;
	//dE
	for(int y=0; y<size.height; y++)
		for(int x=0; x<size.width-1; x++)
			dE.at<float>(y,x) = src.at<float>(y,x+1) -  src.at<float>(y,x);
}

void sigmoid(Mat& src)
{
	int alpha = 1;
	CvSize size = src.size();
	Mat sig(size,CV_32FC1);
	Mat out1(size,CV_32FC1);
	Mat out2(size,CV_32FC1);
	double max = 1;
	double min = 0;

	double *maxVal = &max;
	double *minVal = &min;
	minMaxLoc(src,minVal,maxVal);
	double diff = *maxVal - *minVal;
	double a = 1.0/diff;
	double b =0 -*minVal/diff;
	src.convertTo(sig,CV_32FC1,a,b);

	multiply(sig,Mat::ones(size,CV_32FC1),sig,-2);
	exp(sig,sig);
	add(sig,Mat::ones(size,CV_32FC1),out1);
	subtract(Mat::ones(size,CV_32FC1),sig,out2);
	divide(out2,out1,sig);

	minMaxLoc(sig,minVal,maxVal);
	diff = *maxVal - *minVal;
	a = 1.0/diff;
	b =0 -*minVal/diff;
	sig.convertTo(sig,CV_32FC1,a,b);
	sig.convertTo(src,CV_8UC1,255);

	namedWindow("SigmoidImage",CV_WINDOW_AUTOSIZE);
	imshow("SigmoidImage",src);
	waitKey(0);
}

void Homomorphic(Mat& src)
{
	double cutoff=0.45;
	double order = 3.0;
	double boost =2.0;
	CvSize sz =src.size() ;
	int row = sz.height;
	int col = sz.width;
	//Get source image's spectrum
	Mat real(sz,CV_64FC1);
	Mat imagine = Mat::zeros(sz,CV_64FC1);
	Mat complex(sz,CV_64FC2); 
	Mat Hu(sz,CV_64FC2);
	src.convertTo(real,CV_64FC1,1.0,0.01);
	log(real,real);

	vector<Mat> m;
	m.push_back(real);
	m.push_back(imagine);
	merge(m,complex);//merge:Composes a multi-channel array from several single-channel arrays.
	Mat F;     //spectrum matrix
	dft(complex,F,DFT_COMPLEX_OUTPUT);

	//Make high-boost filter
	// 	
	for (int i=0;i<row;i++)
		for(int j=0;j<col;j++)
			real.at<double>(i,j)= (-0.5+(double)i/(row-1))*(-0.5+(double)i/(row-1))+\
			(-0.5+(double)j/(col-1))*(-0.5+(double)j/(col-1));
	sqrt(real,real);	 // Here real is radius
	real.convertTo(real,CV_64FC1,1.0/cutoff);
	pow(real,2.0*order,real);
	add(real,Mat::ones(sz,CV_64FC1),real);
	divide(1.0,real,real); //real for low pass filter

	Mat H(sz,CV_64FC1);
	for (int x=0;x<row;x++)
		for (int y=0;y<col;y++)
			H.at<double>(x,y) =real.at<double>((x+row/2)%row,(y+col/2)%col);

	scaleAdd(H,1.0/boost-1,Mat::ones(sz,CV_64FC1),H);
	Mat sv;
	H.convertTo(sv,CV_8UC1,255);
//	imwrite("HighBoostFilter.jpg",sv);
	split(F,m);
	multiply(m[0],H,m[0]);
	multiply(m[1],H,m[1]);
	merge(m,F);

	dft(F,real,DFT_INVERSE+DFT_REAL_OUTPUT);


	real.convertTo(real,CV_64FC1,1.0/(row*col));
	exp(real,real);
	//Saturate image,remove the low||high grayscale pixels
	vector<double> sortedVec;
	for (int i=0;i<row;i++)
		for (int j=0;j<col;j++)
			sortedVec.push_back(real.at<double>(i,j));
	sort(sortedVec.begin(),sortedVec.end());
	double minVal = sortedVec[row*col*10/100];
	double maxVal = sortedVec[row*col*98/100];
	for (int i=0;i<row;i++)
		for (int j=0;j<col;j++)
		{
			double tmp = real.at<double>(i,j);
			if(tmp<minVal)
				real.at<double>(i,j) = minVal;
			else if(tmp>maxVal)
				real.at<double>(i,j) = maxVal;
		}

		minMaxLoc(real,&minVal,&maxVal);
		real.convertTo(src,CV_8UC1,255.0/(maxVal-minVal),-255.0*minVal/(maxVal-minVal));
/*
		namedWindow("HomoImg",CV_WINDOW_AUTOSIZE);
		imshow("HomoImg",src);
		imwrite("D:\\twoproject\\slicncut2\\HomoImg.bmp",src);*/
	//	waitKey(0);
}

void mor(IplImage* src)
{
	cvThreshold(src,src,128,255,CV_THRESH_BINARY);
	cvErode(src,src,NULL,2);
	cvDilate(src,src,NULL,2);
	cvNamedWindow("src_mor");
	cvShowImage("src_mor",src);

}


// #include "StdAfx.h"
// #include "preprocess.h"
// 
// 
// preprocess::preprocess(void)
// {
// }
// 
// 
// preprocess::~preprocess(void)
// {
// }
// 
// void preprocess::Homomorphic(Mat& src)
// {
// 	double cutoff=0.45;
// 	double order = 3.0;
// 	double boost =2.0;
// 	CvSize sz =src.size() ;
// 	int row = sz.height;
// 	int col = sz.width;
// 	//Get source image's spectrum
// 	Mat real(sz,CV_64FC1);
// 	Mat imagine = Mat::zeros(sz,CV_64FC1);
// 	Mat complex(sz,CV_64FC2); 
// 	Mat Hu(sz,CV_64FC2);
// 	src.convertTo(real,CV_64FC1,1.0,0.01);
// 	log(real,real);
// 
// 	vector<Mat> m;
// 	m.push_back(real);
// 	m.push_back(imagine);
// 	merge(m,complex);
// 	Mat F;     //spectrum matrix
// 	dft(complex,F,DFT_COMPLEX_OUTPUT);
// 
// 	//Make high-boost filter
// 	// 	
// 	for (int i=0;i<row;i++)
// 		for(int j=0;j<col;j++)
// 			real.at<double>(i,j)= (-0.5+(double)i/(row-1))*(-0.5+(double)i/(row-1))+\
// 			(-0.5+(double)j/(col-1))*(-0.5+(double)j/(col-1));
// 	sqrt(real,real);	 // Here real is radius
// 	real.convertTo(real,CV_64FC1,1.0/cutoff);
// 	pow(real,2.0*order,real);
// 	add(real,Mat::ones(sz,CV_64FC1),real);
// 	divide(1.0,real,real); //real for low pass filter
// 
// 	Mat H(sz,CV_64FC1);
// 	for (int x=0;x<row;x++)
// 		for (int y=0;y<col;y++)
// 			H.at<double>(x,y) =real.at<double>((x+row/2)%row,(y+col/2)%col);
// 
// 	scaleAdd(H,1.0/boost-1,Mat::ones(sz,CV_64FC1),H);
// 	Mat sv;
// 	H.convertTo(sv,CV_8UC1,255);
// 	imwrite("HighBoostFilter.jpg",sv);
// 	split(F,m);
// 	multiply(m[0],H,m[0]);
// 	multiply(m[1],H,m[1]);
// 	merge(m,F);
// 
// 	dft(F,real,DFT_INVERSE+DFT_REAL_OUTPUT);
// 
// 
// 	real.convertTo(real,CV_64FC1,1.0/(row*col));
// 	exp(real,real);
// 	//Saturate image,remove the low||high grayscale pixels
// 	vector<double> sortedVec;
// 	for (int i=0;i<row;i++)
// 		for (int j=0;j<col;j++)
// 			sortedVec.push_back(real.at<double>(i,j));
// 	sort(sortedVec.begin(),sortedVec.end());
// 	double minVal = sortedVec[row*col*10/100];
// 	double maxVal = sortedVec[row*col*98/100];
// 	for (int i=0;i<row;i++)
// 		for (int j=0;j<col;j++)
// 		{
// 			double tmp = real.at<double>(i,j);
// 			if(tmp<minVal)
// 				real.at<double>(i,j) = minVal;
// 			else if(tmp>maxVal)
// 				real.at<double>(i,j) = maxVal;
// 		}
// 
// 		minMaxLoc(real,&minVal,&maxVal);
// 		real.convertTo(src,CV_8UC1,255.0/(maxVal-minVal),-255.0*minVal/(maxVal-minVal));
// 		namedWindow("HomoImg",CV_WINDOW_AUTOSIZE);
// 		imshow("HomoImg",src);
// 		imwrite("homoimg.jpg",src);
// //		waitKey(0);
// }

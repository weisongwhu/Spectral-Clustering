
#include "SLIC.h"
#include "stdio.h"
#include "minmax.h"
#include <fstream>
#include <core.hpp>


const int dx4[4] = {-1,  0,  1,  0};
const int dy4[4] = { 0, -1,  0,  1};
const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

SLIC::SLIC(void) 
{
	m_lvec = NULL;
	
}


SLIC::~SLIC(void)
{
	if(m_lvec!=NULL)
		delete[] m_lvec;

	
}

void SLIC::PerformSuperPixel(const int& superpixelNumbers ,const int& compactness ,int* klabels)
{
	m_height = img.rows;
	m_width  = img.cols;
	m_NumPixel=m_height*m_width;

	int K = superpixelNumbers;  //the number of the superpixels
	int m = compactness;   //the spatial weight
	bool perturbseeds = true;
	//if(0 == klabels) klabels = new int[m_NumPixel];
	for( int s = 0; s < m_NumPixel; s++ ) klabels[s] = -1;

	vector<double> kseedsl(0);

	vector<double> kseedsx(0);
	vector<double> kseedsy(0);

	vector<double> edgemag(0);

	GetData(m_lvec);
	GetEdgeMap(m_lvec,edgemag);
	InitialKSeeds(kseedsl,  kseedsx, kseedsy, K, perturbseeds, edgemag);
	int STEP = (int)sqrt(double(m_NumPixel)/double(K)) + 2.0;//adding a small value in the even the STEP size is too small.
	PerformIteration(kseedsl, kseedsx, kseedsy, klabels, STEP, edgemag, m);
	
	int numlabels = kseedsl.size();
	int* nlabels = new int[m_NumPixel];
	EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, K);
	{for(int i = 0; i < m_NumPixel; i++ ) klabels[i] = nlabels[i];}
	if(nlabels) delete [] nlabels;
	//AdjustLabels(klabels);	 // For a test, the labels are continous from 0~N-1
}
//====================================================================
// Get three channel and store them in vector 
//====================================================================
void SLIC::GetData(
	double*&					lvec)
{
	lvec = new double[m_NumPixel];
	int i=0;
	for (int row=0;row<m_height;row++)
	{
		uchar* data= img.ptr<uchar>(row);
		for(int col = 0;col<m_width;col++)
			lvec[i++] = *data++;
	}
}

//====================================================================
// Get the edge map of the lab image , according to this formate: 
// G(x;y) =|I(x+1;y)-I(x-1;y)|2 + |I(x;y+1)-I(x;y-1)|2
//====================================================================
void SLIC::GetEdgeMap(
	const double*				lvec,
	vector<double>&				edges)
{
	edges.resize(m_NumPixel,0);
	for( int j = 1; j < m_height-1; j++ )
	{
		for( int k = 1; k < m_width-1; k++ )
		{
			int i = j*m_width+k;

			double dx = (lvec[i-1]-lvec[i+1])*(lvec[i-1]-lvec[i+1]) ;

			double dy = (lvec[i-m_width-1]-lvec[i+m_width-1])*(lvec[i-m_width-1]-lvec[i+m_width-1]) ;
			//edges[i] = (sqrt(dx) + sqrt(dy));
			edges[i] = (dx + dy);
		}
	}
}

void SLIC::InitialKSeeds(
	vector<double>&				kseedsl,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const int&					K,
	const bool&					perturbseeds,
	const vector<double>&	        edgemag)
{
	double step = sqrt(double(m_NumPixel)/double(K));
	//int T = step;
	int xoff = (int)step/2;
	int yoff = (int)step/2;

	int n(0);
	for( int y = 0; y < m_height; y++ )
	{
		int Y = y*step + yoff;
		if( Y > m_height-1 ) break;

		for( int x = 0; x < m_width; x++ )
		{
			int X = x*step + xoff;
			if(X > m_width-1) break;

			int i = Y*m_width + X;

			//_ASSERT(n < K);

			//kseedsl[n] = m_lvec[i];
			//kseedsa[n] = m_avec[i];
			//kseedsb[n] = m_bvec[i];
			//kseedsx[n] = X;
			//kseedsy[n] = Y;
			kseedsl.push_back(m_lvec[i]);
			kseedsx.push_back(X);
			kseedsy.push_back(Y);
			n++;
		}
	}

	if(perturbseeds)
	{
		PerturbSeeds(kseedsl,  kseedsx, kseedsy, edgemag);
	}
}
void SLIC::PerturbSeeds(
	vector<double>&				kseedsl,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	const vector<double>&	edges)
{
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

	int numseeds = kseedsl.size();

	for( int n = 0; n < numseeds; n++ )
	{
		int ox = kseedsx[n];//original x
		int oy = kseedsy[n];//original y
		int oind = oy*m_width + ox;

		int storeind = oind;

		for( int i = 0; i < 8; i++ )
		{
			int nx = ox+dx8[i];//new x
			int ny = oy+dy8[i];//new y

			if( nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
			{
				int nind = ny*m_width + nx;
				if( edges[nind] < edges[storeind])
				{
					storeind = nind;
				}
			}
		}
		if(storeind != oind)
		{
			kseedsx[n] = storeind%m_width;
			kseedsy[n] = storeind/m_width;
			kseedsl[n] = m_lvec[storeind];
		}
	}
}


//===========================================================================
///	PerformSuperpixelSLIC
///
///	Performs k mean segmentation. It is fast because it looks locally, not
/// over the entire image.
//===========================================================================
void SLIC::PerformIteration(
	vector<double>&				kseedsl,
	vector<double>&				kseedsx,
	vector<double>&				kseedsy,
	int*						klabels,
	const int&					STEP,
	const vector<double>&		edgemag,
	const double&				M)
{
	int sz = m_width*m_height;
	const int numk = kseedsl.size();
	//----------------
	int offset = STEP;
	//----------------

	vector<double> clustersize(numk, 0);
	vector<double> inv(numk, 0);//to store 1/clustersize[k] values

	vector<double> sigmal(numk, 0);
	vector<double> sigmax(numk, 0);
	vector<double> sigmay(numk, 0);
	vector<double> distvec(sz, DBL_MAX);

	double invwt = 1.0/((STEP/M)*(STEP/M));

	int x1, y1, x2, y2;
	double l;
	double dist;
	double distxy;
	for( int itr = 0; itr < 10; itr++ )
	{
		distvec.assign(sz, DBL_MAX);
		for( int n = 0; n < numk; n++ )
		{
			y1 = max(0,	kseedsy[n]-offset);
			y2 = min(m_height,	kseedsy[n]+offset);
			x1 = max(0,			kseedsx[n]-offset);
			x2 = min(m_width,	kseedsx[n]+offset);


			for( int y = y1; y < y2; y++ )
			{
				for( int x = x1; x < x2; x++ )
				{
					int i = y*m_width + x;
					l = m_lvec[i];

					dist =	(l - kseedsl[n])*(l - kseedsl[n]) ;

					distxy =(x - kseedsx[n])*(x - kseedsx[n]) +
							(y - kseedsy[n])*(y - kseedsy[n]);

					//------------------------------------------------------------------------
					dist += distxy*invwt;//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
					//------------------------------------------------------------------------
					if( dist < distvec[i] )
					{
						distvec[i] = dist;
						klabels[i]  = n;
					}
				}
			}
		}
		//-----------------------------------------------------------------
		// Recalculate the centroid and store in the seed values
		//-----------------------------------------------------------------
		//instead of reassigning memory on each iteration, just reset.

		sigmal.assign(numk, 0);
		sigmax.assign(numk, 0);
		sigmay.assign(numk, 0);
		clustersize.assign(numk, 0);
		//------------------------------------
		//edgesum.assign(numk, 0);
		//------------------------------------

		{int ind(0);
		for( int r = 0; r < m_height; r++ )
		{
			for( int c = 0; c < m_width; c++ )
			{
				sigmal[klabels[ind]] += m_lvec[ind];
				sigmax[klabels[ind]] += c;
				sigmay[klabels[ind]] += r;
				//------------------------------------
				//edgesum[klabels[ind]] += edgemag[ind];
				//------------------------------------
				clustersize[klabels[ind]] += 1.0;
				ind++;
			}
		}}

		{for( int k = 0; k < numk; k++ )
		{
			if( clustersize[k] <= 0 ) clustersize[k] = 1;
			inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
		}}

		{for( int k = 0; k < numk; k++ )
		{
			kseedsl[k] = sigmal[k]*inv[k];
			kseedsx[k] = sigmax[k]*inv[k];
			kseedsy[k] = sigmay[k]*inv[k];
			//------------------------------------
			//edgesum[k] *= inv[k];
			//------------------------------------
		}}
	}
}
//===========================================================================
///	FindNext
///
///	Helper function for EnforceLabelConnectivity. Called recursively.
//===========================================================================
void SLIC::FindNext(
	const int*					labels,
	int*						nlabels,
	const int&					height,
	const int&					width,
	const int&					h,
	const int&					w,
	const int&					lab,
	int*						xvec,
	int*						yvec,
	int&						count)
{
	int oldlab = labels[h*width+w];
	for( int i = 0; i < 4; i++ )
	{
		int y = h+dy4[i];int x = w+dx4[i];
		if((y < height && y >= 0) && (x < width && x >= 0) )
		{
			int ind = y*width+x;
			if(nlabels[ind] < 0 && labels[ind] == oldlab )
			{
				xvec[count] = x;
				yvec[count] = y;
				count++;
				nlabels[ind] = lab;
				FindNext(labels, nlabels, height, width, y, x, lab, xvec, yvec, count);
			}
		}
	}
}

//===========================================================================
///	EnforceLabelConnectivity
///
///	Some superpixels may be unconnected, Relabel them. Recursive algorithm
/// used here, can crash if stack overflows. This will only happen if the
/// superpixels are very large, otherwise safe.
///		STEPS:
///		1. finding an adjacent label for each new component at the start
///		2. if a certain component is too small, assigning the previously found
///		    adjacent label to this component, and not incrementing the label.
//===========================================================================
void SLIC::EnforceLabelConnectivity(
	const int*					labels,
	const int&					width,
	const int&					height,
	int*						nlabels,
	int&						numlabels,
	const int&					K)
{
	int sz = width*height;		
	{for( int i = 0; i < sz; i++ ) nlabels[i] = -1;}

	const int SUPSZ = sz/K;
	//------------------
	// labeling
	//------------------
	int lab(0);
	int i(0);
	int adjlabel(0);//adjacent label
	int* xvec = new int[sz];//worst case size
	int* yvec = new int[sz];//worst case size
	{for( int h = 0; h < height; h++ )
	{
		for( int w = 0; w < width; w++ )
		{
			if(nlabels[i] < 0)
			{
				nlabels[i] = lab;
				//-------------------------------------------------------
				// Quickly find an adjacent label for use later if needed
				//-------------------------------------------------------
				{for( int n = 0; n < 4; n++ )
				{
					int x = w + dx4[n];
					int y = h + dy4[n];
					if( (x >= 0 && x < width) && (y >= 0 && y < height) )
					{
						int nindex = y*width + x;
						if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
					}
				}}
				xvec[0] = w; yvec[0] = h;
				int count(1);
				FindNext(labels, nlabels, height, width, h, w, lab, xvec, yvec, count);
				//-------------------------------------------------------
				// If segment size is less then a limit, assign an
				// adjacent label found before, and decrement label count.
				//-------------------------------------------------------
				if(count <= (SUPSZ >> 2))
				{
					for( int c = 0; c < count; c++ )
					{
						int ind = yvec[c]*width+xvec[c];
						nlabels[ind] = adjlabel;
					}
					lab--;
				}
				lab++;
			}
			i++;
		}
	}}
	//------------------
	numlabels = lab;
	//------------------
	if(xvec) delete [] xvec;
	if(yvec) delete [] yvec;
}




void SLIC::AdjustLabels(int* labels)
{
	bool index[2000] = {};
	//std::cout<<index[1]<<std::endl;
	for(int i=0;i<m_NumPixel;i++)
	{
		if (index[labels[i]]==false)
		{
			index[labels[i]] = true;
			std::cout<<labels[i]<<std::endl;
		}
	}
}

void F_Gray2Color(Mat& gray_mat, Mat& color_mat)
{
// 	if(color_mat)
// 		cvZero(color_mat);

//	int stype = CV_MAT_TYPE(gray_mat->type), dtype = CV_MAT_TYPE(color_mat->type);
	if (gray_mat.channels()>1)
	{
		cvtColor(gray_mat,gray_mat,CV_BGR2GRAY);
	}
	int rows = gray_mat.rows, cols = gray_mat.cols;
//将灰度分布延生到0~255的区间上
	int maxVal = 0;
	for(int i=0;i<rows;i++)
		for(int j=0;j<cols;j++)
		{
			if (gray_mat.at<uchar>(i,j)>maxVal)
			{
				maxVal = gray_mat.at<uchar>(i,j);
			}
		}
	int ratio = 1;
	if (maxVal)
	{
		ratio = 255/maxVal;
	}
		gray_mat = gray_mat*ratio;
	// 判断输入的灰度图和输出的伪彩色图是否大小相同、格式是否符合要求   
//	if (CV_ARE_SIZES_EQ(gray_mat, color_mat) && stype == CV_8UC1 && dtype == CV_8UC3)
//	{
// 		CvMat* red = cvCreateMat(rows, cols, CV_8U);      // 红色分量
// 		CvMat* green = cvCreateMat(rows, cols, CV_8U);   // 绿色分量
// 		CvMat* blue = cvCreateMat(rows, cols, CV_8U);    // 蓝色分量
// 		CvMat* mask = cvCreateMat(rows, cols, CV_8U);   
//		cvSubRS(gray_mat, cvScalar(255), blue);          // blue = 255 - gray
// 		cvCopy(gray_mat, red);                           // red = gray
// 		cvCopy(gray_mat, green);                         // green = gray , if gray < 128
// 		cvCmpS(green, 128, mask, CV_CMP_GE );            //
// 		cvSubRS(green, cvScalar(255), green, mask);      // green = 255 - gray , if gray >= 128
// 		cvConvertScale(green, green, 2.0, 0.0);          // green = 2 * green
// 
// 		// 将蓝绿红三色融合为一幅伪彩色图
// 		cvMerge(blue, green, red, NULL, color_mat);

// 		cvReleaseMat( &red );
// 		cvReleaseMat( &green );
// 		cvReleaseMat( &blue );
// 		cvReleaseMat( &mask );
// 	}
		Mat mask(rows,cols,CV_8UC1);
		vector<Mat> rgb(3);
		Mat& blue = rgb[0];
		Mat& green = rgb[1];
		Mat& red = rgb[2];

		blue = 255-gray_mat;
		red = gray_mat;
		green = gray_mat;
		compare(gray_mat,Mat::ones(rows,cols,CV_8UC1)*128,mask,CV_CMP_GE);
		subtract(Mat::ones(rows,cols,CV_8UC1)*255,green,green,mask);
		green = green*2;

		merge(rgb,color_mat);
}

void DrawContoursAroundSegments(
	Mat&					Img,
	const int*				labels,const int & b,const int & g,const int & r)
{
	if(Img.channels()==1)
		cvtColor(Img,Img,CV_GRAY2BGR);
	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};//八邻域
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};//八邻域
	int width = Img.cols;
	int height =Img.rows;
	int sz = width*height;
	vector<bool> isBoundary(sz,false);
	int j ,k;
	int mainindex(0);
	for(  j = 0; j < height; j++ )
		for( k = 0; k < width; k++ )
		{
			int np(0);
			int cha = 1;
			for( int i = 0; i < 8; i++ )
			{
				int x = k + dx8[i];
				int y = j + dy8[i];

				if( (x >= 0 && x < width) && (y >= 0 && y < height) )
				{
					int index = y*width + x;
					if( false == isBoundary[index] )//comment this to obtain internal contours
						if( labels[mainindex] != labels[index] ) 
							np++;
				}
			}
			if (np >1)//change to 2 or 3 for thinner lines
			{
				Img.at<Vec3b>(j,k) =Vec3b (b,g,r);
				isBoundary[mainindex] = true;
			}
			mainindex++;
		}
}
void Label2IMG(const int* labels, Mat& img)
{
	for (int i=0;i<img.rows;i++)
		for (int j=0;j<img.cols;j++)
		{
			img.at<uchar>(i,j) = labels[j+i*img.cols];
			Scalar intensity = img.at<uchar>(i,j);  
			std::cout<<intensity.val[0]<<std::endl;//输出每个像素点的label
			std::cout<<intensity.val[1]<<std::endl;//0
			std::cout<<intensity.val[2]<<std::endl;//0
		}

}
int SLIC::reallabel(int* labels)
{

	s=0;
	bool index[2000] = {false};
	for(int i=0;i<img.cols*img.rows;i++)
	{
		if (index[labels[i]]==false)
		{
			index[labels[i]] = true;
			s=s+1;
//			std::cout<<labels[i]<<std::endl;

		}
	}
//	std::cout<<"lastlabel"<<" "<<s<<std::endl;
	return 0;
}
	
	//for(int m=0;m<pImgMatrix.size();m++)
	//{
	//	std::cout<<pImgMatrix[m]<<std::endl;
	//}


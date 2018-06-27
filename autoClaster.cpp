
#include<iostream>
#include<fstream>
#include"engine.h"
#include "autoClaster.h"


int mainAClaster(vector<vector<double>>& wOfNJW, vector<vector<int>>& ClasterInfmation)
{
	vector<vector<double>> xOfNJW;
	vector<vector<double>> dOfNJW(wOfNJW.size(), vector<double>(wOfNJW.size(), 0));
	caculateX(wOfNJW, xOfNJW,dOfNJW);
	vector<vector<double>> vdVector;
	vector<vector<vector<int>>> claster_array;
	vector<double> clusterJ;
	vdVector.push_back(xOfNJW[0]);
	int miniter;
	if (wOfNJW.size()>20)
	{
		miniter = 20;
	}
	else
	{
		miniter = wOfNJW.size();
	}
/*
	ofstream wof1("E:\\宋伟\\IMRGandAutoNum\\ImrgAutoNum\\ImrgAutoNum\\x.txt");
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		for (int j = 0; j < wOfNJW[0].size(); j++)
		{
			wof1 << xOfNJW[i][j] << ' ';
		}
		wof1 << std::endl;
	}
	wof1.close();*/
	for (int i=1;i<miniter;i++)
	{
		vdVector.push_back(xOfNJW[i]);
/*
		int Xflag = 0;//标志，X一行不全为0的超像素的数目
		for (int n = 0; n < vdVector[0].size(); n++)
		{
			int j = 0;
			while (j < vdVector.size() && vdVector[j][n] == 0) { j++; }
			if (j < vdVector.size()) { Xflag++; }
		}
		if (Xflag == vdVector[0].size())
		{*/
			vector<vector<int>> cluster_array_i;
			double J = 0;
			evrotMain(vdVector, cluster_array_i, 1, J);
			claster_array.push_back(cluster_array_i);
			clusterJ.push_back(J);
			cluster_array_i.clear();
/*		}*/
	}
//	int maxi = 4;
	double maxJ = 0;
	for (int i=0;i<clusterJ.size();i++)
	{
		if (maxJ<clusterJ[i])
		{
			maxJ = clusterJ[i];
	//		maxi = i+2;
		}
	}
/*
	fstream fsd("G:\\MATLAB\\MATLAB2015b\\交叉编译\\MNCUT\\D.txt", fstream::out);
	fstream fsx("G:\\MATLAB\\MATLAB2015b\\交叉编译\\MNCUT\\V.txt", fstream::out);
	for (int i = 0; i < maxi; i++)
	{
		for (int j = 0; j < xOfNJW[0].size(); j++)
		{
			fsx << xOfNJW[i][j] << ' ';
		}
		fsx << endl;
	}
	for (int i = 0; i < xOfNJW.size(); i++)
	{
		for (int j = 0; j < xOfNJW[0].size(); j++)
		{
			fsd << dOfNJW[i][j] << ' ';
		}
		fsd << endl;
	}
	fsd.close();
	fsx.close();
//////////////////////////////////////////////////////////////////////////
	Engine *eig1;
	while (!(eig1 = engOpen(NULL)))
	{
		std::cout << endl << "MATLAB引擎启动失败！" << endl;
	}
	engEvalString(eig1, "MNcutK");
	if (engClose(eig1))
	{
		std::cout << "close failure" << endl;
	}
//	delete[] eig1;
	fstream  ifXk("G:\\MATLAB\\MATLAB2015b\\交叉编译\\MNCUT\\Xk.txt", fstream::in);
	for (int i = 0; i < maxi; i++)
	{
		vector<int> line1;
		for (int j = 0; j < xOfNJW[0].size(); j++)
		{
			int Xkij = 0;
			ifXk >> Xkij;
			if (Xkij==1)
			{
				line1.push_back(j);
			}
		}
		ClasterInfmation.push_back(line1);
		line1.clear();
	}
	ifXk.close();*/
//////////////////////////////////////////////////////////////////////////
	for (int i=0;i<clusterJ.size();i++)
	{
		if (abs(maxJ - clusterJ[i]) < 0.001)
		{
			ClasterInfmation.clear();
			ClasterInfmation = claster_array[i];
		}
	}
	vdVector.clear();
	claster_array.clear();
	xOfNJW.clear();
	dOfNJW.clear();
	return 0;
}

int caculateX(vector<vector<double>>& wOfNJW, vector<vector<double>>& xOfNJW, vector<vector<double>>& dOfNJW)
{

/*
	ofstream wof1("E:\\宋伟\\IMRGandAutoNum\\ImrgAutoNum\\ImrgAutoNum\\w.txt");
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		for (int j = 0; j < wOfNJW[0].size(); j++)
		{
			wof1 << wOfNJW[i][j] << ' ';
		}
		wof1 << std::endl;
	}
	wof1.close();
	Engine *ep;
	while (!(ep = engOpen(NULL)))
	{
		//fprintf(stderr, "\n matlab引擎启动失败!\n");
		std::cout << endl << "MATLAB引擎启动失败！" << endl;
	}
	engEvalString(ep, "clculateXfromW");
	if (engClose(ep))
	{
		std::cout << "close failure" << endl;
	}
	std::fstream V_in("E:\\宋伟\\IMRGandAutoNum\\ImrgAutoNum\\ImrgAutoNum\\V.txt", std::ofstream::in);
	double vtd = 0;
	xOfNJW.clear();
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		vector<double> line1;
		for (int j = 0; j < wOfNJW[0].size(); j++)
		{
			V_in >> vtd;
			line1.push_back(vtd);
		}
		xOfNJW.push_back(line1);
		line1.clear();
	}*/

///////////////////////////////////////////////////	

	vector<vector<double>> lOfNJW(wOfNJW.size(), vector<double>(wOfNJW.size(), 0));
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		double sumI = 0;
		for (int j = 0; j < wOfNJW.size(); j++)
		{
			sumI += wOfNJW[i][j];
		}
		dOfNJW[i][i] = 1/sqrt(sumI);     //直接计算D矩阵的-1/2次方
	}
	vector<vector<double>> dl(wOfNJW.size(), vector<double>(wOfNJW.size(), 0));
	matrixMtpl(dOfNJW, wOfNJW, dl);
	matrixMtpl(dl, dOfNJW, lOfNJW);
	//Mat mLOfNJW(wOfNJW.size(), wOfNJW[0].size(), CV_64FC1);
	CvMat* mLOfNJW = cvCreateMat(wOfNJW.size(), wOfNJW[0].size(), CV_64FC1);
	for (int i=0;i<wOfNJW.size();i++)
	{
		for (int j = 0; j < wOfNJW[0].size(); j++)
		{
	//		mLOfNJW.at<double>(i, j) = lOfNJW[i][j];
			cvmSet(mLOfNJW, i, j, lOfNJW[i][j]);
		}
	}
//	Mat eigenValue(wOfNJW.size(), 1, CV_64FC1);
//	Mat eigenVector(wOfNJW.size(), wOfNJW.size(), CV_64FC1);
	CvMat* eigenValue = cvCreateMat(wOfNJW.size(), 1, CV_64FC1);
	CvMat* eigenVector = cvCreateMat(wOfNJW.size(), wOfNJW.size(), CV_64FC1);
	cvEigenVV(mLOfNJW, eigenVector, eigenValue);
	xOfNJW.clear();
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		
		vector<double> line2;
		for (int j = 0; j < wOfNJW.size(); j++)
		{
			line2.push_back(cvmGet(eigenVector, i, j));
		}
		xOfNJW.push_back(line2);
		line2.clear();
	}


	vector<double> vdValue(wOfNJW.size(), 0);
	vector<vector<double>> vdVector(wOfNJW.size(), vector<double>(wOfNJW.size(), 0));
	for (int i=0;i<wOfNJW.size();i++)
	{
		for (int j = 0; j < wOfNJW.size(); j++)
		{
		//	vdVector[i][j] = eigenVector.at<double>(i, j);
			vdVector[i][j] = cvmGet(eigenVector, i, j);
		}
		//vdValue[i] = eigenValue.at<double>(i);
		vdValue[i] = cvmGet(eigenValue, i,0);
	}
	lOfNJW.clear();
	dl.clear();

/*
	mLOfNJW.release();
	eigenValue.release();
	eigenVector.release();*/
	/*vector<int> flagSort(wOfNJW.size(), 0);
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		flagSort[i] = i;
	}
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		for (int j=0;j<wOfNJW.size()-1-i;j++)
		{
			if (vdValue[j] < vdValue[j + 1])
			{
				double dchange = vdValue[j];
				vdValue[j] = vdValue[j + 1];
				vdValue[j + 1] = dchange;
				int ichange = flagSort[j];
				flagSort[j] = flagSort[j + 1];
				flagSort[j + 1] = ichange;
			}
		}
	}
	for (int i = 0; i < wOfNJW.size(); i++)
	{
		xOfNJW.push_back(vdVector[flagSort[i]]);
	}*/
	vdVector.clear();
	vdValue.clear();
	return 0;
}


int  matrixMtpl(vector<vector<double>>& MtplA, vector<vector<double>>& MtplB, vector<vector<double>>& MtplResult) 
{
	vector<vector<double>> MtplA1 = MtplA;
	vector<vector<double>> MtplB1 = MtplB;
	MtplResult.clear();
	for (int i=0;i<MtplA1.size();i++)
	{
		vector<double> line1;
		for (int n = 0; n < MtplB1[0].size(); n++)
		{
			double sumMtpl = 0;
			for (int j = 0; j < MtplA1[i].size(); j++)
			{
					sumMtpl += MtplA1[i][j] * MtplB1[j][n];
			}
			line1.push_back(sumMtpl);
		}
		MtplResult.push_back(line1);
		line1.clear();
	}
	MtplA1.clear();
	MtplB1.clear();
	return 0;
}


int evrotMain(vector<vector<double>>& vdVector, vector<vector<int>>& cluster_array, int method_rot,double& J)
{
	int dim = vdVector.size();
	int ndata = vdVector[0].size();
	int angle_num = (int)(dim*(dim - 1) / 2);
	/*build index mapping*/
	int *ik = new int[angle_num];
	int *jk = new int[angle_num];
	int cntk = 0;
	for (int i = 0; i < dim - 1; i++)
	{
		for (int j = i + 1; j <= dim - 1; j++)
		{
			ik[cntk] = i;
			jk[cntk] = j;
			cntk++;
		}
	}

	int max_iter = 200;
	double dJ, J_new, J_old1, J_old2, J_up, J_down;
	double alpha;
	int iterJ=0;
	vector<double> vtheta(angle_num, 0);
	vector<double> vtheta_new(angle_num, 0);
	/*initial quality*/
	J = evqual(vdVector, dim, ndata);
	J_old1 = J;
	J_old2 = J;
	while (iterJ<max_iter)
	{
		iterJ++;
		for (int d=0;d<angle_num;d++)
		{
			if (method_rot==0)
			{
			/*用数值导数梯度下降*/
				alpha = 0.1;
				/*move up*/
				vtheta_new[d] = vtheta[d] + alpha;
				vector<vector<double>> vdVectorRot;
				rotate_givens(vdVector, vdVectorRot, vtheta_new, ik, jk, angle_num, dim);
				J_up = evqual(vdVectorRot, dim, ndata);
				vdVectorRot.clear();
				/*move down*/
				vtheta_new[d] = vtheta[d] - alpha;
				rotate_givens(vdVector, vdVectorRot, vtheta_new, ik, jk, angle_num, dim);
				J_down = evqual(vdVectorRot, dim, ndata);
				vdVectorRot.clear();
				if (J_up > J || J_down > J)
				{
					if (J_up > J_down)
					{
						vtheta[d] = vtheta[d] + alpha;
						vtheta_new[d] = vtheta[d];
						J = J_up;
					}
					else {
						vtheta[d] = vtheta[d] - alpha;
						vtheta_new[d] = vtheta[d];
						J = J_down;
					}
				}
			}
			else/*矩阵导数梯度下降*/
			{
				alpha = 1.0;
				dJ = evqualitygrad(vdVector, vtheta, ik, jk, angle_num, d, dim, ndata);
				vtheta_new[d] = vtheta[d] - alpha*dJ;
				vector<vector<double>> vdVectorRot;
				rotate_givens(vdVector, vdVectorRot, vtheta_new, ik, jk, angle_num, dim);
				J_new = evqual(vdVectorRot, dim, ndata);
				if (J_new>J)
				{
					vtheta[d] = vtheta_new[d];
					J = J_new;
				}
				else
				{
					vtheta_new[d] = vtheta[d];
				}
				vdVectorRot.clear();
		//////////////////////////////////////////////////////////////////////////
			}
		}
		if (iterJ>2)
		{
			if (J - J_old2 < 0.001) { break; }
		}
		J_old2 = J_old1;
		J_old1 = J;
	}
	vector<vector<double>> vdVectorRot;
	rotate_givens(vdVector, vdVectorRot, vtheta_new, ik, jk, angle_num, dim);
	autoClaster(vdVectorRot, cluster_array, dim, ndata);
	vdVector.clear();
	vdVector = vdVectorRot;
	delete[] ik;
	delete[] jk;
	vtheta.clear();
	vtheta_new.clear();
	vdVectorRot.clear();
	return 0;
}

double evqual(vector<vector<double>>& vdVector, int dim, int ndata)
{
/////////根据矩阵dim*ndata,求J值
	vector<double> max_values(ndata,0);
	vector<int> max_index(ndata, 0);
	for (int i=0;i<ndata;i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (max_values[i]*max_values[i]<=vdVector[j][i]*vdVector[j][i])
			{
				max_values[i] = vdVector[j][i];
				max_index[i] = j;
			}
		}
	}
	double J = 0;
	for (int i=0;i<ndata;i++)
	{
		for (int j = 0; j < dim; j++)
		{
			double vM = vdVector[j][i] / max_values[i];
			J += vM*vM;

		}
	}

	J =1.0-(J / ndata - 1.0) / dim;
	max_values.clear();
	max_index.clear();
	return J;

}

int rotate_givens(vector < vector<double>>& vdVector, vector < vector<double>> &vdVectorRot, 
	vector<double>& vtheta, int* ik, int* jk, int angle_num, int dim)
{
////做吉文斯旋转，Y:dim*ndata;
	vector<vector<double>> Uab(dim, vector<double>(dim, 0));
	computeUab(vtheta, Uab, 0, angle_num - 1, ik, jk, dim);
	vector<vector<double>> vdVectorZhuanZhi(vdVector[0].size(), vector<double>(vdVector.size(), 0));
	for (int i = 0; i < vdVector[0].size(); i++)
	{
		for (int j=0;j<vdVector.size();j++)
		{
			vdVectorZhuanZhi[i][j] = vdVector[j][i];
		}
	}
	vector<vector<double>> vdVectorRotZhuanZhi;
	matrixMtpl(vdVectorZhuanZhi, Uab, vdVectorRotZhuanZhi);
	vdVectorRot.clear();
	for (int i = 0; i < vdVector.size(); i++)
	{
		vector<double> line1;
		for (int j = 0; j < vdVector[0].size(); j++)
		{
			line1.push_back(vdVectorRotZhuanZhi[j][i]);
		}
		vdVectorRot.push_back(line1);
	}

	Uab.clear();
	vdVectorZhuanZhi.clear();
	vdVectorRotZhuanZhi.clear();
	return 0;
}

int computeUab(vector<double>& vtheta, vector<vector<double>>& Uab, int a, int b, int* ik, int* jk, int dim)
{
//////根据角度theta，求吉文斯旋转的旋转矩阵
	for (int i=0;i<dim;i++)
	{
		Uab[i][i] = 1.0;
	}
	if (b < a) { return 0; }
	for (int i = a; i <= b; i++)
	{
		double tt = vtheta[i];
		for (int j=0;j<dim;j++)
		{
/*
			double u_ik = Uab[ik[i]][j] * cos(tt) - Uab[jk[i]][j] * sin(tt);
			Uab[jk[i]][j] = Uab[ik[i]][j] * sin(tt) + Uab[jk[i]][j] * cos(tt);
			Uab[ik[i]][j] = u_ik;*/
			double u_ik = Uab[j][ik[i]] * cos(tt) - Uab[j][jk[i]]*sin(tt);
			Uab[j][jk[i]] = Uab[j][ik[i]] * sin(tt) + Uab[j][jk[i]] * cos(tt);
			Uab[j][ik[i]] = u_ik;
		}
	}
	////////////////////////[cos(t) -sin(t);sin(t) cos(t)]*[A;B]
	return 0;
}


double evqualitygrad(vector<vector<double>>& vdVector, vector<double>& vtheta,
	int* ik, int* jk, int angle_num, int angle_index, int dim, int ndata)
{ 
//////////////dJ/dtheta
	vector<vector<double>> matrixV(dim, vector<double>(dim,0));
	computeV(vtheta, matrixV, angle_index, ik, jk, dim);
	vector<vector<double>> matrixU1(dim, vector<double>(dim, 0));
	vector<vector<double>> matrixU2(dim, vector<double>(dim, 0));
	computeUab(vtheta, matrixU1, 0, angle_index - 1, ik, jk, dim);
	computeUab(vtheta, matrixU2, angle_index + 1, angle_num - 1, ik, jk, dim);
	vector<vector<double>> matrixA(dim, vector<double>(ndata, 0));
	computeA(vdVector, matrixA, matrixU1, matrixV, matrixU2);
	matrixU1.clear();
	matrixU2.clear();
	matrixV.clear();
	vector<vector<double>> matrixY(dim, vector<double>(ndata,0));
	rotate_givens(vdVector, matrixY, vtheta, ik, jk, angle_num, dim);
////////////////////////////////////////////////////////////////////////////////
	vector<double> max_value(ndata, 0);
	vector<int>    max_index(ndata, 0);
	for (int i = 0; i < ndata; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (max_value[i] * max_value[i] <= (matrixY[j][i] * matrixY[j][i]))
			{
				max_value[i] = matrixY[j][i];
				max_index[i] = j;
			}
		}
	}
///////found max of each column;
	double dJ = 0;
	for (int i = 0; i <dim; i++)
	{
		for (int j=0;j<ndata;j++)
		{
			double tmp1 = (matrixA[i][j] / max_value[j])*(matrixY[i][j]/ max_value[j]);
			double tmp2 = (matrixA[max_index[j]][j]/ max_value[j]) * (matrixY[i][j] / max_value[j] )* (matrixY[i][j] /  max_value[j]);
			dJ += tmp1 - tmp2;

		}
	}
	dJ = 2 * dJ / ndata / dim;
	max_index.clear();
	max_value.clear();
	matrixY.clear();
	matrixA.clear();

///////待续；



	return dJ;
}

int computeV(vector<double>& vtheta, vector<vector<double>>& matrixV, int k, int* ik, int* jk, int dim)
{
	matrixV[ik[k]][ik[k]] = -sin(vtheta[k]);
	matrixV[ik[k]][jk[k]] = cos(vtheta[k]);
	matrixV[jk[k]][ik[k]] = -cos(vtheta[k]);
	matrixV[jk[k]][jk[k]] = -sin(vtheta[k]);
	return 0;
}

int computeA(vector<vector<double>>& vdVector,vector<vector<double>>& matrixA, vector<vector<double>>& matrixU1, vector<vector<double>>& matrixV, vector<vector<double>>& matrixU2)
/////////////A=X*U*V*U,A:ndata*dim
{
	int dim = vdVector.size();
	int ndata = vdVector[0].size();
	vector<vector<double>> matrixA1(ndata, vector<double>(dim, 0));
	vector<vector<double>> matrixUVU(dim, vector<double>(dim, 0));
	vector<vector<double>> vdVectorZhuanZhi(vdVector[0].size(),vector<double>(vdVector.size(),0));
	for (int i=0;i<vdVector.size();i++)
	{
		for (int j=0;j<vdVector[0].size();j++)
		{
			vdVectorZhuanZhi[j][i] = vdVector[i][j];
		}
	}
	matrixMtpl(matrixU1, matrixV, matrixUVU);
	matrixMtpl(matrixUVU, matrixU2, matrixUVU);
	matrixMtpl(vdVectorZhuanZhi, matrixUVU, matrixA1);
/*
	matrixMtpl(vdVectorZhuanZhi, matrixU1, matrixA1);
	matrixMtpl(matrixA1, matrixV, matrixA1);
	matrixMtpl(matrixA1, matrixU2, matrixA1);*/
	vdVectorZhuanZhi.clear();
/////A=X*U1*V*U2;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < ndata; j++)
		{
			matrixA[i][j] =matrixA1[j][i];
		}
	}
	matrixUVU.clear();
	matrixA1.clear();
	return 0;


}

int autoClaster(vector<vector<double>>& vdVector, vector<vector<int>>& clasterd_array, int dim, int ndata)
{
	clasterd_array.clear();
	vector<double> max_value(ndata, 0);
	vector<int> max_index(ndata, 0);
//	vector<int> claster_count(dim, 0);
////找出每一行中最大的数目，并统计属于每一个dim的像素的个数
	
	for (int i=0;i<ndata;i++)
	{
		for (int j = 0; j < dim; j++)
		{
//			if (j == 0) { max_index[i] = -1; }
			if (max_value[i] <= vdVector[j][i] * vdVector[j][i])
			{
/*
				if (max_index[i] >= 0)
				{
					claster_count[max_index[i]]--;
				}
//				claster_count[j]++;*/
				max_value[i] = vdVector[j][i] * vdVector[j][i];
				max_index[i] = j;
			}
		}
	}
//////////////////////////////////////////////////////////////////////////
	for (int i = 0; i < dim; i++)
	{
		vector<int> line1;
		for (int j = 0; j < ndata; j++)
		{
			if (max_index[j] == i)
			{
				line1.push_back(j);
			}
		}
		clasterd_array.push_back(line1);
		line1.clear();
	}
	max_index.clear();
	max_value.clear();
	return 0;
}





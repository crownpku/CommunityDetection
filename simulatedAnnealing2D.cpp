#include "simulatedAnnealing2D.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
using namespace std;

simulatedAnnealing::simulatedAnnealing(void)
{
}


simulatedAnnealing::~simulatedAnnealing(void)
{
}

void simulatedAnnealing::evolve(const vector<vector<double > >& idist, vector<int >& itagx, vector<int >& itagy )		
	// input the adjacency matrix, out the relabeled tag x and y
{
	static const double mDescentRatio = DescentRatio;
	static const int mMaxStep = MaxStep;
	vector<vector<double > > dist = idist;
	//idist is the adjacent matrix, described by vector<vector<int>>
	//So now the job is too input the edge list and their weight into this idist
	mNodeNum = NumberofNode;
	mRange = Range;
	vector<int  > mtagx(mNodeNum, 0);
	vector<int  > mtagy(mNodeNum, 0);
	for(int i = 0; i < mRange; i++)
	{
		for (int j = 0; j < mRange; j++){
		mtagx[i*mRange+j] = i;
		mtagy[i*mRange+j] = j;
	}
	}
	this->mTemp = this->findEnergy(dist, mtagx, mtagy);

	for(int i = 0; i < mMaxStep; i++)
	{
		cout << i << endl;
		for(int k = 0; k < mNodeNum; k++)
		{
			int a = this->pickIndex();
			int b = this->pickIndex();//pick two nodes
			if(a == b)
			{
				continue;
			}
			double eb = this->totalEnergy(a, b, mtagx, mtagy, dist);
			vector<int> tmtagx = mtagx;
			vector<int> tmtagy = mtagy;
			swap(tmtagx[a], tmtagx[b]);
			swap(tmtagy[a], tmtagy[b]);
			double ea = this->totalEnergy(a, b, tmtagx, tmtagy, dist);
			if(ea < eb)
			{
				mtagx = tmtagx;
				mtagy = tmtagy;
			}
			else
			{
				double t = (double)(rand() % 32767) / (double)32767;
				double threshold = exp(-(ea - eb) / mTemp);
				if(t < threshold)
				{
					mtagx = tmtagx;
					mtagy = tmtagy;
				}
			}
		}
		mTemp *= mDescentRatio;
		double s = this->findEnergy(dist, mtagx, mtagy);
	}
	itagx = mtagx;
	itagy = mtagy;
}

double simulatedAnnealing::findEnergy(const vector<vector<double> >& dist, const vector<int >& itagx, const vector<int >& itagy)
{
	double s = 0.;
	for(int k = 0; k < mNodeNum; k++)
	{
		for(int l = 0; l < mNodeNum; l++)
		{
			s = dist[k][l] * ( double( (itagx[k]-itagx[l])*(itagx[k]-itagx[l]) + (itagy[l]-itagy[k])*(itagy[l]-itagy[k]) ) ) ;
		}
	}
	return s;
}



double simulatedAnnealing::hamilton(int ix, int jx, int iy, int jy, double w)
{
	double d =  ( double( (ix-jx)*(ix-jx) + (iy-jy)*(iy-jy)) );	
	double r = d * w;
	return r;
}

double simulatedAnnealing::totalEnergy(int i, int j, const vector< int> &x, const vector< int> &y, const vector<vector<double> >& w)
{
	double s = 0;
	for(int k = 0; k < mNodeNum; k++)
	{
		if(k != i && k != j)
		{
			s += hamilton(x[k], x[i], y[k], y[i], w[k][i]);
			s += hamilton(x[k], x[j], y[k], y[j], w[k][j]);
		}
	}
	return s;
}

int simulatedAnnealing::pickIndex(void)
{
	int n = (int)((double)(rand() % 32767) / (double)32767 * mNodeNum);
	return n;

}

int main()
{
	int nNode=NumberofNode;
	int nEdge=NumberofEdge;
	
	int tempNode1;
	int tempNode2;
  double tempWeight;
	int i,j;

  string InputFile="footballTS.pairs";
	string OutputFile="footballTS_2DStructualProfile_22_5000.txt";
	ifstream fin;
  ofstream fout;
	fin.open(InputFile.c_str());
	if (!fin.is_open()) { cout<<"Input File Open Failed"<<endl; exit(0);}	
	fout.open(OutputFile.c_str());
	if (!fout.is_open()) { cout<<"OutputOut File Open Failed"<<endl; exit(0);}

  simulatedAnnealing crown;
	vector<int> itagx;
	vector<int> itagy;
	vector<vector<double> > matrix; //Adjacent Matrix
	vector<vector<int> > layout;
	for(i=0; i<nNode; i++) matrix.push_back(vector<double>(nNode));
	for(i=0; i<Range; i++) layout.push_back(vector<int> (Range));
	
	for(i=0; i<nNode; i++){
		for(j=0; j<nNode; j++){
			matrix[i][j]=0;
		}
	}
	for(i=0; i<Range; i++){
		for (j=0; j<Range; j++){
			layout[i][j]=0;
			}
		}
	for(i=0; i<nEdge; i++){
		fin>>tempNode1;
		//cout<<tempNode1<<" ";
		fin>>tempNode2;
	//	cout<<tempNode2<<" ";
		//cout<<endl;
		//fin>>tempWeight;
		//tempNode1=tempNode1-1;
		//tempNode2=tempNode2-1;
		matrix[tempNode1][tempNode2]=1;
		matrix[tempNode2][tempNode1]=1;
	}
	
	crown.evolve(matrix, itagx, itagy);
	
	for(i=0; i<RealNumberofNode; i++) layout[itagx[i]][itagy[i]]=i;
  for (i=0; i<Range; i++){
  	for (j=0; j<Range; j++){
  		fout<<layout[i][j]<<" ";
  	}
  	fout<<endl;
  }
	fin.close();
	fout.close();
	return 0;
}

#pragma once
#include <vector>
#define NumberofNode 484//Usually 4*RealNumberofNode
#define Range 22
#define NumberofEdge 613
#define RealNumberofNode 114
#define DescentRatio 0.99
#define MaxStep 5000

using namespace std;


class simulatedAnnealing
{
public:
	simulatedAnnealing(void);
	virtual ~simulatedAnnealing(void);


public:

	void evolve(const vector<vector<double > >& idist, vector<int>& itagx, vector<int>& itagy);
	//idist is the adjacent matrix, itagx[i] and itagy[i] is the position of the node i



private:
	vector<vector<double> > mDist;
	int mNodeNum;
	int mRange;//mNodeNum should be mRange*mRange
	double mTemp;


private:
	double hamilton(int ix, int jx, int iy, int jy, double w);
	double totalEnergy(int i, int j, const vector< int> &x, const vector< int> &y, const vector<vector<double> >& w);
	int pickIndex(void);
	double findEnergy(const vector<vector<double> >& dist, const vector<int >& itagx, const vector<int >& itagy);
};


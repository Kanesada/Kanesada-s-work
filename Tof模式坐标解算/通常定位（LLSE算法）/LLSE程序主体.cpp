#include <armadillo>
#include <math.h>
using namespace arma;
using namespace std;

double distance_2(rowvec a, rowvec b)     // Calculates the distance between two coords
{
	rowvec temp_vec = a - b;
	mat temp_mat = temp_vec * temp_vec.st();
	return temp_mat(0);
}

/*
This multilateration algorithm is called Linear Least squre Estimation(LLSE), programed by Kai.
*If anchor coordinates have 3 dimension, we will calculate 3D results. 4 or more anchors are requested.
Otherwise we will calculate 2D results. 3 or more anchors are requested.
*/
colvec LLSE(mat &anchors, rowvec &tof, int clesolver3D)       
{
	mat A;
	mat S;
	colvec b;
	rowvec doa;
	rowvec origin;
	colvec result;
	int i;
	#define CMperS				((double)299702547)				/**< Speed of Light m/s */
	
	if(clesolver3D)
	{
		S = anchors;
	}
	else
	{
		S = anchors.cols(0,1);
	}

	
	//int numAnchors = anchors.n_rows;
	doa = tof*CMperS;
	A.zeros(numAnchors-1, S.n_cols);
	b.zeros(numAnchors-1, 1);
	origin.zeros(1, S.n_cols);
	
	for(i=0;i<numAnchors-1;i++)
	{
		A.row(i) = 2*(S.row(i+1) - S.row(0));
		b(i) = pow(doa(i+1), 2) - pow(doa(0), 2) 
		       + distance_2(S.row(0), origin) - distance_2(S.row(i+1), origin);
	}
	
	result = pinv(A)*b;
	
	return -1*result;
}
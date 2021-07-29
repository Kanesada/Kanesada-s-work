#include <iostream>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>

using namespace arma;
using namespace std;

int       numAnchors;
rowvec    tagposition;
rowvec    tof;
mat       anchors;
#define CMperS				((double)299702547)				/**< Speed of Light m/s */


double gaussrand(float E, float V)    //Generate the gaussian noise
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;
             
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);
         
    phase = 1 - phase;
	X = X * V + E;     
 
    return X;
}


double distance_2(rowvec a, rowvec b)     // Calculates the distance between two coords
{
	rowvec temp_vec = a - b;
	mat temp_mat = temp_vec * temp_vec.st();
	return temp_mat(0);
}

/*
Generate test data.
Including anchors coordinates and tof measurements with gaussian noise.
*/

void Generate_data()             
{
	int i;
	/* If the anchors are deployed correctly.
	anchors<<0<<0<<0<<endr
	       <<20<<0<<20<<endr
		   <<0<<20<<20<<endr
		   <<20<<20<<0<<endr; */
	
/* If the anchors are deployed nearly in a line. */	
	anchors<<0<<-0.01<<0<<endr
	       <<20<<0.01<<20<<endr
	       <<0<<0.1<<20<<endr
	       <<20<<0.1<<0<<endr;
		   
	numAnchors = (int)anchors.n_rows;
	
	tagposition<<5<<10<<10;

	tof.zeros(1, numAnchors);
	for(i=0;i<numAnchors;i++)
	{
		tof(i) = (sqrt(distance_2(tagposition, anchors.row(i))) + gaussrand(0,0.1))
					/ CMperS;
	}
	cout<<"The anchors' coordinates are: "<<endl<<anchors<<endl;
	cout<<"The tagposition is : "<<endl<<tagposition;
	cout<<"The tof measurements with 10cm noise are: "<<tof<<endl;
}

/*
Calculate the 2 nomal of A matrix's Cols.
If the 2 nomal is less than threshold, return 0.
*/
int Detect_Line(mat A, double threshold)
{	int i;
	for(i=0;i<A.n_cols;i++)
	{
		colvec temp_col = A.col(i);
		mat temp_mat = temp_col.st() * temp_col;
		double x = temp_mat(0);
		if(x<threshold) return 0;
	}
	return 1;
}




/*
This multilateration algorithm is called Linear Least squre Estimation(LLSE), programed by Kai.
If anchor coordinates have 3 dimension, we will calculate 3D results. 4 or more anchors are requested.
Otherwise we will calculate 2D results. 3 or more anchors are requested.
Use Line_detect to check if the anhors are deployed in a line. 
If colvecs' 2 nomals in A are less than threshold, means the anchors are deployed in line.
LLSE can not cope with this case.
In any circumstance, we should avoid deploying anchors in a line.
*/
colvec LLSE(mat &anchors, rowvec &tof, int clesolver3D, int Line_detect, double threshold)       
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
	
	if(Line_detect && (Detect_Line(A, threshold)==0))   
	{                 
		//如果检测到选中基站全部被布置在一条直线上
		//则跳转使用三角质心算法或两基站定位算法
		//但需要依靠先验知识例如地图信息
		//确定定位点在基站连线的哪边
		//实际使用中应当尽量避免这种情况
		cout<<"The anchors are deployed in a line!"<<endl<<
		"Jump to another Algorithm!"<<endl
		<<"The LLSE result will not be the best solution!"<<endl<<endl;	
	}
	
	if (clesolver3D)
	{
		result = pinv(A)*b;
	}
	else
	{
		result = pinv(A)*b;
		result.set_size(3, 1);
		result.at(2, 0) = 0; 
	}
	
	return -1*result;
}

int main()
{
	colvec    result1;
	Generate_data();
	result1 = LLSE(anchors, tof,1, 1, 0.1);
	cout<<"The result of LLSE is : "<<result1<<endl;
	return 0;
}
	
	
	
	
	







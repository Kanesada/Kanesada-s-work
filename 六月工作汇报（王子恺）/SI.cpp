#include <iostream>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>
#include <string>
using namespace arma;
using namespace std;

#define CMperS				((double)299702547)				/**< Speed of Light m/s */
//#define DWT_TIME_UNITS      ((double) 1.0/499.2e6/128.0)	/**< 15.65e-12 s */

int       numAnchors;									/**< number of anchor in the set to be used in the multilateration*/ 
mat       S, R, P_Id, W, I;
rowvec    refAnchorCoords;
rowvec    tagposition;
rowvec	  tdoa;
colvec ddoa, delta;
mat result1, B;



//Generate the gaussian noise with 0 mean mu and 1 square error theta
//E is the mean,  V is the square error
//X = X * V + E;

//#include <stdlib.h>                            
//#include <math.h>
double gaussrand(float E, float V)    //adjusted
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
	X = X * V + E;     //adjusted
 
    return X;
}



void init(void)
{	
	double n0, n1, n2;      //noise in ddoa
	rowvec d0, d1, d2 , d3 , d4;     //the distance between tag and anchors
	
	numAnchors = 4;
	tagposition<<1.2<<1.2<<1;      //Set the tag position
	cout<<"The tag position is : "<<tagposition<<endl;

	refAnchorCoords<<0<<0<<1;  //Set the master anchor to the origin
	cout<<"The master anchor's coordinates are: "<<endl<<refAnchorCoords<<endl;
	//Set the slave anchors' coordinates
	S<<2<<0<<0<<endr
	 <<2<<2<<1<<endr
	 <<0<<2<<0<<endr;
	 cout<<"The slave anchor coordinates matrix S is: "<<endl<<S<<endl;
	 
	//Calculate the R matrix
	R = sqrt(sum(square(S), 1));
	
	//Calculate the distance between tag and anchors
	d0 = square(tagposition - refAnchorCoords);
	d1 = square(tagposition - S.row(0));
	d2 = square(tagposition - S.row(1));
	d3 = square(tagposition - S.row(2)); 
	d0(0) = sqrt(accu(d0));
	d1(0) = sqrt(accu(d1));
	d2(0) = sqrt(accu(d2));
	d3(0) = sqrt(accu(d3));
	
	ddoa<<d1(0)-d0(0)<<d2(0)-d0(0)<<d3(0)-d0(0);          //Calculate the true ddoa value
	cout<<"The true ddoa is :"<<endl<<ddoa;
	
	n0 = gaussrand(0,0.01);     //Generate the ddoa measurement with gaussian noise
	n1 = gaussrand(0,0.01);
	n2 = gaussrand(0,0.01);
	ddoa(0) += n0;            
	ddoa(1) += n1;
	ddoa(2) += n2;	
	cout<<"ddoa measurement is: "<<endl<<ddoa<<endl;
	
	tdoa = (ddoa.st()) / CMperS;     //Generate the tdoa measurement
	cout<<"tdoa measurement is: "<<endl<<tdoa<<endl;



}








rowvec multilaterate(rowvec tdoa, /*mat anchor_coords,*/ int cleSolver3D)
{
	
	/*colvec ddoa, delta;
	double a;
	mat result1, B;*/
	double a;
	
	
	ddoa = tdoa.st() * CMperS; //convert time to distance (s to m)
	//ddoa = ddoa.rows(1, numAnchors - 1);
	delta = square(R) - square(ddoa);
	I.eye(numAnchors - 1, numAnchors - 1);

	W = I;																/*W is a weight matrix, which is surely be symmetric positive definite. 
																		For now , W is an identity matrix. It can be a true weight matrix
																		if a NLOS/MultiPath detection is implemented.*/
	
	a = as_scalar(ddoa.st() * ddoa);
	mat P_Id = I - ((ddoa * ddoa.st() )/ a);
	B = S.st() * P_Id * W * P_Id ;

	
	if (cleSolver3D)
	{
		result1 = refAnchorCoords.st() + 0.5 * solve(B * S,	B * delta);       //transpotion may be needed.																							    
	}
	else
	{
		result1 = refAnchorCoords.st() + 0.5 * solve(B * S,	B * delta);     //transpotion may be needed.
		result1(2,0) = 0;
	}
	
	result1 = result1.st();
	ofstream coordLog("Coordinates.csv", ios::app);           //Record the coordinate
	//ofstream fs("Coordinates.csv");
	result1.raw_print(coordLog);
	coordLog.close();
	
	ofstream ddoaLog("DDOA.csv", ios::app);           //Record the DDOA
	//ofstream ddoaLog("DDOA.csv");
	(ddoa.st()).raw_print(ddoaLog);
	ddoaLog.close();

	return result1;
}



int main(int argc, char** argv)
{	rowvec result;
	init();
	result = multilaterate(tdoa,0);
	cout<<"The result is "<<endl<<result;
	return 0;	
}
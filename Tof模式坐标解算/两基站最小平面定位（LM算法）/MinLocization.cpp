#include <iostream>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>
#include <string>
#include <stdio.h>
#include <levmar.h>

#ifndef LM_DBL_PREC
#error Example program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


/* the following macros concern the initialization of a random number generator for adding noise */
#undef REPEATABLE_RANDOM
#define DBL_RAND_MAX (double)(RAND_MAX)

#ifdef _MSC_VER // MSVC
#include <process.h>
#define GETPID  _getpid
#elif defined(__GNUC__) // GCC
#include <sys/types.h>
#include <unistd.h>
#define GETPID  getpid
#else
#warning Do not know the name of the function returning the process id for your OS/compiler combination
#define GETPID  0
#endif /* _MSC_VER */

#ifdef REPEATABLE_RANDOM
#define INIT_RANDOM(seed) srandom(seed)
#else
#define INIT_RANDOM(seed) srandom((int)GETPID()) // seed unused
#endif



using namespace arma;
using namespace std;

#define CMperS				((double)299702547)				/**< Speed of Light m/s */
//#define DWT_TIME_UNITS      ((double) 1.0/499.2e6/128.0)	/**< 15.65e-12 s */

int       solver3D = 1;
int       numAnchors;									/**< number of anchor in the set to be used in the multilateration*/ 
mat       S;
rowvec    refAnchorCoords;
rowvec    last_result;
rowvec    tagposition;
colvec    doa;
double    set_height;



/*Generate the gaussian noise with 0 mean mu and 1 square error theta
E is the mean,  V is the square error
X = X * V + E;
*/
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


/*
Calculates the distance between two coords
*/
double Distance(rowvec a, rowvec b)     
{
	rowvec temp_vec = a - b;
	mat temp_mat = temp_vec * temp_vec.st();
	return sqrt(temp_mat(0));
}



void init(void)
{	
	mat temp;
	register int i;
		
	tagposition<<2<<2<<1.6;      //Set the tag position
	cout<<"The tag position is : "<<tagposition<<endl;

	//Set the anchors' coordinates
	S<<0<<0<<1<<endr
	<<2<<0<<0<<endr
	<<2<<2<<1<<endr
	<<0<<2<<0<<endr;
	
	
		/*S<<0<<0<<1<<endr
		<<2<<0<<0<<endr;*/
	
	cout<<"The anchor coordinates matrix S is: "<<endl<<S<<endl;
	 
	numAnchors = S.n_rows;
	doa.zeros(numAnchors, 1);      
	

	
	for(i=0;i<numAnchors;i++)
	{
		doa(i) = Distance(tagposition, S.row(i));
	}
	cout<<"The true doa is :"<<endl<<doa;
	
	
	for(i=0;i<numAnchors;i++)               //Generate the tdoa measurement
	{
		doa(i) += gaussrand(0,0.1); 	
	}	
	cout<<"doa measurement is: "<<endl<<doa<<endl;	   
	set_height = tagposition(2) + 0.0535;
}


void doafunc(double *p, double *x, int m, int n, void *data)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
{
	register int i;
	double d0;
	if(!solver3D)
	{
		for(i=0; i<n; i++)
		{
			x[i]=sqrt( pow(p[0] - S(i,0), 2) 
				+ pow(p[1] - S(i,1), 2) + pow(set_height - S(i,2), 2));
		}
	}
	else
	{
		for(i=0; i<n; i++)
		{
			x[i]=sqrt( pow(p[0] - S(i,0), 2) 
				+ pow(p[1] - S(i,1), 2) + pow(p[2] - S(i,2), 2));
		}
}
}

void jacdoafunc(double *p, double *jac, int m, int n, void *data)
{   if (solver3D)
	{
		register int i, j; 
		/* fill Jacobian row by row */
		for(i=j=0; i<n; i++)
		{
			jac[j++] = 2*(p[0] - S(i,0));
			jac[j++] = 2*(p[1] - S(i,1));
			jac[j++] = 2*(p[2] - S(i,2));
		}
	}
	else
	{
		register int i, j; 
		/* fill Jacobian row by row */
		for(i=j=0; i<n; i++)
		{
			jac[j++] = 2*(p[0] - S(i,0));
			jac[j++] = 2*(p[1] - S(i,1));
		}
	}
}


/*
'last' is the last_result
'raw_result' is the result before determination, maybe symmetrical to the real result.
'S' is the anchoors' coordinates matrix.  
This function is used to determin whether the result need to symmetrized.
*/
/*int Determin(rowvec raw_result, mat S)   
{
	double k;
	double y;
	if((S(0,0) == S(1,0)))                 // The anchor line parallel to Y axis.  
	{
		if((S(0,0) - raw_result(0)) > 0) return -1;
		else return 1;
	}
	else if(S(0,1) == S(1,1))         // The anchor line parallel to X axis.
	{
		if((S(0,1) - raw_result(1)) > 0) return -1;
		else return 1;
	}
	else                             // Normal Cases
	{
		k = (S(0,1) - S(1,1)) / (S(0,0) - S(1,0));  //Gradient of the line made by two anchors.
		y = k*(raw_result(0) - S(1,0)) + S(1,1);
		if(raw_result(1) > y)
		{
			return 1;	     //That means the row result is higer than the line made by two anchors
		}
		else
		{
			return -1;        //That means the row result is lower than the line made by two anchors
		}
	}
}*/
	


	
/*
'raw_result' is the result before symmetrized.
'S' is the anchoors' coordinates matrix.
This function is used to symmetrize a 2D point about a line.  
*/
rowvec Trans_Side(mat S, rowvec raw_result)
{
	rowvec trans_result = raw_result;
	
	if((S(0,0) == S(1,0)))                 // The anchor line parallel to Y axis.            
	{    
		trans_result(0) = 2*S(0,0) - raw_result(0);
	}
	else if(S(0,1) == S(1,1))         // The anchor line parallel to X axis.
	{
		trans_result(1) = 2*S(0,1) - raw_result(1);
	}
	else                             // Normal Cases
	{
		double A = S(0,1) - S(1,1);
		double B = S(1,0) - S(0,0);
		double C = S(0,0)*S(1,1) - S(0,1)*S(1,0);

		trans_result(0) = raw_result(0) - 2 * A * ((A * raw_result(0) 
						+ B * raw_result(1) + C) / (A * A + B * B));
		trans_result(1) = raw_result(1) - 2 * B * ((A * raw_result(0) 
						+ B * raw_result(1) + C) / (A * A + B * B));
	}

	/*double d = 0.5 * Distance(raw_result, trans_result);
	
	if(d <= range) 
	{
		trans_result = raw_result;
	}
	*/
	cout<<"Trans result is : "<<trans_result<<endl;
	
	return trans_result;
}


rowvec solve_doa()
{
	if(!solver3D)
	{		
		const int n=numAnchors, m=2; // n measurements, 2 parameters:x,y
		double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
		register int i;
		int ret;
		int det_last, det_this;   //The determination result of last result and this result
		rowvec result;
		result.zeros(1,3);
		float range = 0.1;           //This is the distance range that ignore the symmetrical phenomenon.


		/* initial parameters estimate: (1.0, 0.0, 0.0) */
		p[0]=last_result(0); p[1]=last_result(1);


		for(i=0;i<n;i++)             //Transport ddoa to x 
		{
			x[i] = doa(i);
		}

		/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
		opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
		opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

		/* invoke the optimization function */
		ret=dlevmar_der(doafunc, jacdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // with analytic Jacobian
		//ret=dlevmar_dif(doafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // without Jacobian
		printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
		printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], set_height);
		if(ret== -1) result = last_result;
		else 
		{
			result(0) = p[0];
			result(1) = p[1];
			result(2) = set_height;
		}
		
		/*det_last = Determin(last_result, S);
		det_this = Determin(result, S);
		if( (det_last*det_this) < 0)      // Detect wheter the result need to symmetrize.     
		{
			result = Trans_Side(S, result, range);    // Symmetrize the result about the anchor line.
		}*/
		
		/*result(0) = 0;
		result(1) = 0;
		result(2) = set_height;
		cout<<endl<<"Set the result to (0,0), simulate the wrong result."<<endl;*/
		if(S.n_rows<3)
		{
			rowvec trans_result = Trans_Side(S, result);
			double d_tag_line = 0.5 * Distance(result, trans_result);    // d_tag_line is the distance between tag and anchor line. 
			cout<<"The distance between result and anchor line is: "<<d_tag_line<<endl;
			if(d_tag_line > range)
			{
				double d1 = Distance(result, last_result);
				double d2 = Distance(trans_result, last_result);
				cout<<"The distance between raw_result and last_result is: "<<d1<<endl;
				cout<<"The distance between trans_result and last_result is: "<<d2<<endl;
				if(d1>d2)
				{			
					result = trans_result;
					cout<<endl<<"Choose the symmetrized point!"<<endl;
				}
				else cout<<endl<<"No need to symmetrize ~"<<endl;
			}
		}
		
		last_result = result;
		return result;	
	}
	else
	{
		const int n=numAnchors, m=3; // n measurements, 3 parameters:x,y,z
		double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
		register int i;
		int ret;
		int det_last, det_this;   //The determination result of last result and this result
		rowvec result;
		result.zeros(1,3);

		/* initial parameters estimate: (1.0, 0.0, 0.0) */
		p[0]=last_result(0); p[1]=last_result(1); p[2]=last_result(2);

		for(i=0;i<n;i++)             //Transport ddoa to x 
		{
			x[i] = doa(i);
		}

		/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
		opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
		opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

		/* invoke the optimization function */
		//ret=dlevmar_der(doafunc, jacdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // with analytic Jacobian
		ret=dlevmar_dif(doafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // without Jacobian
		printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
		printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);
		if(ret== -1) result = last_result;
		else 
		{
			result(0) = p[0];
			result(1) = p[1];
			result(2) = p[2];
		}
		
		last_result = result;
		return result;	
		
	}
}


int main(int argc, char** argv)
{	
	rowvec result;
	init();
	last_result<<1.6<<1.6<<1;
	cout<<"The Last result is: "<<last_result<<endl;
	result = solve_doa();
	cout<<"The final result is "<<endl<<result;
	return 0;	
}

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

int cleSolver3D = 0;
int       numAnchors;									/**< number of anchor in the set to be used in the multilateration*/ 
mat       S, R, P_Id, W, I;
rowvec    refAnchorCoords;
rowvec    last_result;
rowvec    tagposition;
rowvec	  tdoa;
colvec ddoa, delta;
mat result1, B;
double set_height;



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
	//rowvec d0, d1, d2 , d3 , d4;     //the distance between tag and anchors
	rowvec d;
	mat temp;
	register int i;
	
	tagposition<<1.2<<1.2<<1;      //Set the tag position
	cout<<"The tag position is : "<<tagposition<<endl;

	refAnchorCoords<<0<<0<<1;  //Set the master anchor to the origin
	cout<<"The master anchor's coordinates are: "<<endl<<refAnchorCoords<<endl;
	//Set the slave anchors' coordinates
	S<<2<<0<<0<<endr
	 <<2<<2<<1<<endr
	 <<0<<2<<0<<endr;
	cout<<"The slave anchor coordinates matrix S is: "<<endl<<S<<endl;
	 
	numAnchors = S.n_rows + 1;
	d.zeros(1,numAnchors);
	ddoa.zeros(numAnchors - 1, 1);
	
	//Calculate the R matrix
	R = sqrt(sum(square(S), 1));
	

	temp = sqrt((tagposition - refAnchorCoords) * (tagposition - refAnchorCoords).st());
	d(0) = temp(0);
	
	for(i=1;i<numAnchors;i++)
	{
		temp = sqrt((tagposition - S.row(i-1)) * (tagposition - S.row(i-1)).st());
		d(i) = temp(0);
	}
	
	for(i=0;i<numAnchors-1;i++)
	{
		ddoa(i) = d(i+1)-d(0);	              //Calculate the true ddoa value
	}
	cout<<"The true ddoa is :"<<endl<<ddoa;
	
	for(i=0;i<numAnchors-1;i++)               //Generate the tdoa measurement
	{
		ddoa(i) += gaussrand(0,0.01); 	
	}	
	cout<<"ddoa measurement is: "<<endl<<ddoa<<endl;	
	tdoa = (ddoa.st()) / CMperS;     
	cout<<"tdoa measurement is: "<<endl<<tdoa<<endl;
	set_height = tagposition(2) + gaussrand(0,0.01);
}


void tdoafunc(double *p, double *x, int m, int n, void *data)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
{
	if (cleSolver3D)
	{
		register int i;
		double d0;
		d0 = sqrt(pow(p[0] - refAnchorCoords[0], 2) + pow(p[1] - refAnchorCoords[1], 2) 
					+ pow(p[2] - refAnchorCoords[2], 2));    
					
		for(i=0; i<n; i++)
		{
			x[i]=sqrt( pow(p[0] - S(i,0), 2) + pow(p[1] - S(i,1), 2) + pow(p[2] - S(i,2), 2)) - d0;
		}
	}
	else
	{
		register int i;
		double d0;
		d0 = sqrt(pow(p[0] - refAnchorCoords[0], 2) + pow(p[1] - refAnchorCoords[1], 2) 
					+ pow(set_height - refAnchorCoords[2], 2));  
		
		for(i=0; i<n; i++)
		{
			x[i]=sqrt( pow(p[0] - S(i,0), 2) + pow(p[1] - S(i,1), 2) + pow(set_height - S(i,2), 2)) - d0;
		}
		
	}
}



void jactdoafunc(double *p, double *jac, int m, int n, void *data)
{   if (cleSolver3D)
	{
		register int i, j; 
		register double temp1, temp2;
  
		/* fill Jacobian row by row */
		for(i=j=0; i<n; i++)
		{
			temp1 = sqrt( pow(p[0] - S(i,0),2) + pow(p[1] - S(i,1),2) + pow(p[2] - S(i,2),2) );
			temp2 = sqrt( pow(p[0] - refAnchorCoords(0),2) + pow(p[1] - refAnchorCoords(1),2) + pow(p[2] - refAnchorCoords(2),2));  
			jac[j++] = (p[0] - S(i,0))/temp1 - (p[0] - refAnchorCoords(0))/temp2;
			jac[j++] = (p[1] - S(i,1))/temp1 - (p[1] - refAnchorCoords(1))/temp2;
			jac[j++] = (p[2] - S(i,2))/temp1 - (p[2] - refAnchorCoords(2))/temp2;
		}
	}
	else
	{
		register int i, j; 
		register double temp1, temp2;
		
		/* fill Jacobian row by row */
		for(i=j=0; i<n; i++)
		{
			temp1 = sqrt( pow(p[0] - S(i,0),2) + pow(p[1] - S(i,1),2) + pow(set_height - S(i,2),2) );
			temp2 = sqrt( pow(p[0] - refAnchorCoords(0),2) + pow(p[1] - refAnchorCoords(1),2) + pow(p[2] - refAnchorCoords(2),2));  
			jac[j++] = (p[0] - S(i,0))/temp1 - (p[0] - refAnchorCoords(0))/temp2;
			jac[j++] = (p[1] - S(i,1))/temp1 - (p[1] - refAnchorCoords(1))/temp2;

		}
	}
}


rowvec solve_tdoa()
{
	if (cleSolver3D)
	{
		const int n=numAnchors-1, m=3; // n measurements, 3 parameters:x,y,z
		double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
		register int i;
		int ret;
		rowvec result;
		result.zeros(1,3);

		/* initial parameters estimate: (1.0, 0.0, 0.0) */
		p[0]=last_result(0); p[1]=last_result(1); p[2]=last_result(2);


		for(i=0;i<n;i++)             //Transport ddoa to x 
		{
			x[i] = ddoa(i);
		}

		/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
		opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
		opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

		/* invoke the optimization function */
		ret=dlevmar_der(tdoafunc, jactdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // with analytic Jacobian
		//ret=dlevmar_dif(tdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // without Jacobian
		printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
		printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], p[2]);

		result(0) = p[0];
		result(1) = p[1];
		result(2) = p[2];
		
	    return result;
	}
	else
	{
		const int n=numAnchors-1, m=2; // n measurements, 3 parameters:x,y
		double p[m], x[n], opts[LM_OPTS_SZ], info[LM_INFO_SZ];
		register int i;
		int ret;
		rowvec result;
		result.zeros(1,3);


		/* initial parameters estimate: (1.0, 0.0, 0.0) */
		p[0]=last_result(0); p[1]=last_result(1);


		for(i=0;i<n;i++)             //Transport ddoa to x 
		{
			x[i] = ddoa(i);
		}

		/* optimization control parameters; passing to levmar NULL instead of opts reverts to defaults */
		opts[0]=LM_INIT_MU; opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
		opts[4]=LM_DIFF_DELTA; // relevant only if the finite difference Jacobian version is used 

		/* invoke the optimization function */
		ret=dlevmar_der(tdoafunc, jactdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // with analytic Jacobian
		//ret=dlevmar_dif(tdoafunc, p, x, m, n, 100, opts, info, NULL, NULL, NULL); // without Jacobian
		printf("Levenberg-Marquardt returned in %g iter, reason %g, sumsq %g [%g]\n", info[5], info[6], info[1], info[0]);
		printf("Best fit parameters: %.7g %.7g %.7g\n", p[0], p[1], set_height);

		result(0) = p[0];
		result(1) = p[1];
		result(2) = set_height;
		
	    return result;

	}
}


rowvec multilaterate(rowvec tdoa /*mat anchor_coords,*/)
{
	double a;
	ddoa = tdoa.st() * CMperS; //convert time to distance (s to m)
	cout<<"DDOA measurement is : "<<ddoa<<endl;
	delta = square(R) - square(ddoa);
	I.eye(numAnchors - 1, numAnchors - 1);
	W = I;																
	a = as_scalar(ddoa.st() * ddoa);
	mat P_Id = I - ((ddoa * ddoa.st() )/ a);
	B = S.st() * P_Id * W * P_Id ;

	if (cleSolver3D)
	{
		result1 = refAnchorCoords.st() + 0.5 * solve(B * S,	B * delta);       																						    
	}
	else
	{
		result1 = refAnchorCoords.st() + 0.5 * solve(B * S,	B * delta);     
		result1(2,0) = 0;
	}
	
	result1 = result1.st();
	last_result = result1;         /////Transport the SI result

	return result1;
}



int main(int argc, char** argv)
{	
	rowvec result_SI;
	rowvec result;
	init();
	//result_SI = multilaterate(tdoa);
	//cout<<"The SI result is "<<endl<<result_SI;
	result = solve_tdoa();
	cout<<"The final result is "<<endl<<result;
	return 0;	
}
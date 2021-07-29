/*! ------------------------------------------------------------------------------------------------------------------
* @file multilateration.cpp
* @brief DecaWave RTLS multilateration algorithm classes and functions
*
* @attention
* Copyright 2008, 2014-2015 (c) DecaWave Ltd, Dublin, Ireland.
*
* All rights reserved.
*
*/
#include "multilateration.h"

/**
* @brief MultiLateration calculateSwTSw function
*        Calculate the pseudoinverse matrix of S. S is the relative coordinate matrix between anchors and reference anchor.
*/
int MultiLateration::calculateSwTSw(void)
{
	mat x = S.st()*S;  //x is the product of S * the transpose of S matrix

	double det = 0;    //determinant

	if (x.n_rows == 3) 
	{
		det =
			(x(0, 0)*x(1, 1) - x(0, 1)*x(1, 0))*x(2, 2) +
			(x(0, 2)*x(1, 0) - x(0, 0)*x(1, 2))*x(2, 1) +
			(x(0, 1)*x(1, 2) - x(0, 2)*x(1, 1))*x(2, 0);
	}
	else
	{
		det = x(0, 0)*x(1, 1) - x(0, 1)*x(1, 0);
	}

	if (det == 0)      //need to exit as cannot divide by zero 
	{
		return 0;
	}

	mat r(x.n_rows, x.n_rows);      // inverse of x calculated below      //inv(S.st()*S)
	
	if (x.n_rows == 3)              // If 3d solver is used
	{
		r(0, 0) = (x(1, 1)*x(2, 2) - x(1, 2)*x(2, 1)) / det;
		r(0, 1) = -(x(0, 1)*x(2, 2) - x(0, 2)*x(2, 1)) / det;
		r(0, 2) = (x(0, 1)*x(1, 2) - x(0, 2)*x(1, 1)) / det;

		r(1, 0) = -(x(1, 0)*x(2, 2) - x(1, 2)*x(2, 0)) / det;
		r(1, 1) = (x(0, 0)*x(2, 2) - x(0, 2)*x(2, 0)) / det;
		r(1, 2) = -(x(0, 0)*x(1, 2) - x(0, 2)*x(1, 0)) / det;

		r(2, 0) = (x(1, 0)*x(2, 1) - x(1, 1)*x(2, 0)) / det;
		r(2, 1) = -(x(0, 0)*x(2, 1) - x(0, 1)*x(2, 0)) / det;
		r(2, 2) = (x(0, 0)*x(1, 1) - x(0, 1)*x(1, 0)) / det;

	}
	else if (x.n_rows == 2)         // If 2d solver is used
	{
		r(0, 0) = x(1, 1) / det;
		r(0, 1) = -x(0, 1) / det;
		r(1, 0) = -x(1, 0) / det;
		r(1, 1) = x(0, 0) / det;
	}
	else                            // If other cases happen, there is errors. We will not multilaterate.
	{
		return 0;
	}

	R = sqrt(sum(square(S), 1));		
	Sw = r * S.st();               // Sw = inv(S.st()*S) * S.st();  Sw is the Pseudoinverse matrix of S. S is the relative coordinate matrix between anchors and reference anchor.
	SwTSw = Sw.st() * Sw;          // SwTSw is the product of the transpose of Sw matrix * Sw 

	return 1;
}

/**
* @brief MultiLateration setAnchors function
*        Sets the anchor's coordinates (of the anchors used in the multilateration function)
* @param anchor_coords calculate anchor coordinates relative to point 1 (reference anchor)
*/
int MultiLateration::setAnchors(mat &anchor_coords, int cleSolver3D)
{
	numAnchors = anchor_coords.n_rows;

	if (cleSolver3D)
	{
		S.set_size(numAnchors - 1, 3);
		refAnchorCoords = anchor_coords.row(0);  
	}
	else
	{
		S.set_size(numAnchors - 1, 2);
		refAnchorCoords = anchor_coords(0, 0, size(1, 2));
	}


	for (int r = 1; r < numAnchors; r++)
	{
		for (uword c = 0; c < S.n_cols; c++)
		{
			S(r-1,c) = anchor_coords(r,c) - anchor_coords(0,c);
		}
	}

	return calculateSwTSw();
}

/**
* @brief MultiLateration multilaterate function
*        This is the spherical-intersection (SX) method described in Julius O. Smith and Jonathan S. Abel, 
*        "Closed-Form Least-Squares Source Location Estimation from Range-Difference Measurements", IEEE Trans.
*        Acoustics, Speech and Signal Processing, Vol. ASSP-35, pp. 1661-1669, Dec. 1987.
*
*		 Verifyed to DecaWave Matlab model v.1.1 29/05/2015
*
* @brief 2d solver 
*        If 2D anchor coordinates are passed into the algorithm then the result will be 2D. Since the solver 
*        does not need to consider reducing the error in the Z dimension it will often have improved performance 
*        in the X and Y dimensions. The implementation of 2d solver come from the Matlab function in revision 1624.
* @param tdoa is the tdoa array, which contains reference anchor (tdoa of 0) in 1st position and then the rest of
*             the tdoas (w.r.t. reference anchor) from the set of the anchors to be used in this mulitlateration 
* @param anchor_coords are the coordinates of the anchors used in the tdoa array, relative to reference anchor.
* @param cleSolver3D is the flag to specify if 3D or 2D solver is used in the system
*
*/
rowvec MultiLateration::multilaterate(rowvec tdoa, mat anchor_coords, int cleSolver3D)
{
	colvec ddoa, delta;
	vec a, b, c, delta2rsd1, delta2rsd2;
	double rs1;
	mat result1;

	ddoa = (tdoa.st() - tdoa(0)) * CMperS; //convert time to distance (s to m)
	ddoa = ddoa.rows(1, numAnchors - 1);
	delta = square(R) - square(ddoa);

	a = ddoa.st()*SwTSw*ddoa; a(0) = 4 - 4 * a(0);
	b = ddoa.st() * SwTSw * delta; b(0) = 4 * b(0);
	c = delta.st() * SwTSw * delta; c(0) = -1 * c(0);

	double t = b(0)*b(0) - 4 * a(0)*c(0);

	
	if (t < 0)                     //if < 0 cannot mulitlaterate - no answer
	{
		result1.set_size(1, 1);
		result1(0, 0) = 0;
	}
	else
	{
		rs1 = (-b(0) + sqrt(t)) / (2 * a(0));

		if (rs1 < 0)			//if < 0 should not mulitlaterate according Matlab model v.1.1
		{
			result1.set_size(1, 1);
			result1(0, 0) = 0;
		}
		else
		{

			delta2rsd1 = delta - ddoa * 2.0 * rs1;

			if (cleSolver3D)
			{
				result1 = (refAnchorCoords.st() + ((Sw*delta2rsd1) * 0.5)).st();  //all calculation is relative to p1, so we need to add with p1 coordinate
			}
			else
			{
				result1 = (refAnchorCoords.st() + ((Sw*delta2rsd1) * 0.5)).st();  //all calculation is relative to p1, so we need to add with p1 coordinate
				result1.set_size(1, 3);
				result1.at(0, 2) = 0;                                             // 2D solver set Height to be 0
			}
		}
	}
	return result1.row(0);
}

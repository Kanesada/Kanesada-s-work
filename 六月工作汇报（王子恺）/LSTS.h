/*! ------------------------------------------------------------------------------------------------------------------
* @file	deca_cskalman.h
* @brief Kalman filter for wireless clock synchronisation
*
* @attention
*
* Copyright 2013 (c) DecaWave Ltd, Dublin, Ireland.
*
* All rights reserved.
*
* @author
*/

#ifndef __DECA_CSKALMAN_H
#define __DECA_CSKALMAN_H

#include <iostream>
#include <math.h>

#define TIME_UNITS              (1.0/499.2e6/128.0)


#define CLOCK_PERIOD            (0x10000000000LL) //period is 40 bits in DW1000

#define CLOCK_PERIOD_SEC        (CLOCK_PERIOD*TIME_UNITS)
#define CLOCK_HALF_PERIOD_SEC   (CLOCK_PERIOD_SEC/2.0)

#define KAL_START_OREJ			(30) //how many runs before outlier rejection can start

//Notation follows Welch & Bishop
class CS_LSTS
{
public:
    // Constructor
    LSTS_init(void);

	void LSTS_update(double timeRxs, double timeTxs);
    double syncTOA_FULL(double blinkTime);

    void log(double timeRxs, double timeTxs, int seqnum);

	int  id;

private:

    double prevTxTime;
	double prevRxTime;
	double aij, bij;
	int l;
	double rou_l;
	float rou_a, rou_b;
	float u;
	double tao_i0, taoj_0;
	double tao_i, tao_j;
	int Tx_period_num, Rx_period_num;
	double timeTxs_FULL, timeRxs_FULL;
	double weight;        //The fusion weight of two method
	double drifting_fused   //The fused drifting
	int temp 
	int counter;			//< Counts the number of runs of the process function from start
	double Tx_dt, Rx_dt;						     

    bool init;
	bool iterate;




};

#endif //__DECA_CSKALMAN_H

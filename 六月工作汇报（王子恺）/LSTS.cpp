/*! ------------------------------------------------------------------------------------------------------------------
* @file	deca_cskalman.cpp
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
#include "deca_cskalman.h"
#include "dsPrint.h"

#define MAX_OUTLIERS  (9)	  //if we get this many outliers in a row then reset the rejection algorithm	
#define CHISQUAREDMAX (150.0) //maximum for anomalous input detection

static const double ONE_OVER_ROOT_TWO = 0.70710678118654752440084436210485;


/**
* @brief CS_Kalman Constructor - initialises the LSTS filter variables
*
*/
CS_LSTS::LSTS_init(void)
{
    
    prevTxTime = 0.0;
	prevRxTime = 0.0;
    init = false;
	iterate = false;
    counter = 0;


    aij = 1;   
	bij = 0;    
	l = 1;    
	rou_l = 0;   
	rou_a = rou_b = 0.5;    
	u = 0.5;    
	tao_i0 = tao_j0 = 0;
	tao_i = tao_j = 0;
	temp = 0;
	Tx_period_num = Rx_period_num = 0;
	timeTxs_FULL = timeRxs_FULL = 0;
	weight = 0.5;
	drifting_fused =1;
	Tx_dt = Rx_dt = 0;

	
}

/**
* @brief CS_Kalman updateFilter function runs the Kalman filter
*
* @param timeRxs is the slave anchor's time in seconds of CCP frame reception
* @param timeTxs is, the time in seconds, that is sum of the master anchor's CCP frame 
*        transmission time and the constant ToF between the master and the slave
*/
void CS_LSTS::LSTS_update(double timeRxs, double timeTxs)           
{
    if (counter < 75)
    {
        counter++;
    }
	else
	{
		iterate = true;
	}
	


    if (init)  // Check if filter has been initialised
    {

        double rou_tmp;
        double _tao_i, _tao_j;       //the transfer time and receive time of the last sample
		

		//Clock wrapping detection
        Tx_dt = fmod(timeTxs, CLOCK_PERIOD_SEC) - prevTxTime;
		Rx_dt = fmod(timeRxs, CLOCK_PERIOD_SEC) - prevRxTime;
   

        if (Tx_dt < 0)                           //calculate the time without the clock wrapping
        {
			timeTxs_FULL += CLOCK_PERIOD_SEC;
			Tx_dt += CLOCK_PERIOD_SEC; 
			Tx_period_num = floor(timeTxs / CLOCK_PERIOD_SEC);          //calculate the number of clock period
        }
		
	    if (Rx_dt < 0)
        {
			timeRxs_FULL += CLOCK_PERIOD_SEC;
			Rx_dt += CLOCK_PERIOD_SEC; 
			Rx_period_num = timeRxs / CLOCK_PERIOD_SEC;        //same as above

        }
		
		tao_j = timeRxs_FULL;
		tao_i = timeTxs_FULL;


		for(i = 1; i<=l; i++)              //calculate rou_l
		{
			temp += i*i;
		}
		rou_l = (l*l)/(temp);
			
		aij = (1 - rou_l)*aij + rou_l*((tao_j - tao_j0)/(tao_i - tao_i0));
		bij = (1 - rou_l)*bij + rou_l*((tao_j-tao_j0) - aij*(tao_i - tao_i0));     

	        

			//check for DW1000 clock wrapping
            prevTxTime = timeTxs;
			prevRxTime = timeRxs;
        }
    }
    else
    {                         //  Initialize
        //check for DW1000 clock wrapping
        prevTxTime = fmod(timeTxs, CLOCK_PERIOD_SEC);
		prevRxTime = fmod(timeRxs, CLOCK_PERIOD_SEC);
		
		tao_j0 = timeRxs;           //The origin of LSTS time
		tao_i0 = timeTxs;
        aij = 1.0;
		bij = timeRxs - timeTxs;
        init = true;

		/*check for DW1000 clock wrapping
        prevTxTime = fmod(timeTxs, CLOCK_PERIOD_SEC);  */
    }

} // end  process()


/**
* @brief CS_LSTS syncTOA_FULL function converts the slave anchor's blink ToA to master anchor's time base in Full time 
*
* @param blinkTime is the slave anchor's time of receiving of a blink frame in seconds
*            
* @return the blink time in seconds converted to master anchor's timebase
*/
double CS_LSTS::syncTOA_FULL(double blinkTime)
{
	double timeSinceLast_FULL, syncTime_FULL;

	timeSinceLast = blinkTime + CLOCK_PERIOD_SEC*Rx_period_num;    ///////   the time since last sync     (unless there is a wrapping between the blink with the last sync)

	//if blink came after CCP rx time or the clock has wrapped
	
	if (prevRxTime <= blinkTime)
	{
		//if the clock has wrapped   (Blink first,wrapped,RX next)   
		
		if ((blinkTime - prevRxTime) > CLOCK_HALF_PERIOD_SEC)
		{
			timeSinceLast_FULL -= CLOCK_PERIOD_SEC;
		}

	}
	else //blinkTime later than prevRxTime (last CCP Rx time)     (RX first , wrapped, Blink next)
	{
		//has the clock wrapped
	
		if ((prevRxTime - blinkTime) > CLOCK_HALF_PERIOD_SEC)
		{
			timeSinceLast_FULL += CLOCK_PERIOD_SEC;
		} 
	}

	syncTime_FULL = timeSinceLast_FULL/aij - bij;             //////The high precision division need to be manipulated!!!!!!


	return syncTime_FULL;
}// end  syncTOA_FULL()



double CS_LSTS::LSTS_Fusion( )
{

	double Rx_dt_cali;    //calibrated delta Rx by kalman filter
	double Rx_dt_aij;     //calibrated delta Rx by LSTS

	
	Rx_dt_cali = Rx_dt/x1;
	Rx_dt_aij = Rx_dt/aij;
	weight = pow((Tx_dt - Rx_dt_cali),2) / (pow((Tx_dt - Rx_dt_cali),2) + pow((Tx_dt - Rx_dt_aij),2));
	drifting_fused = weight*aij + (1 - weight)*x1;
	
}
	


/**
* @brief CS_Kalman log function -- log data for performance checking
*
* - called immediately after call to updateFilter function if we are logging
*   data to check that the wireless sync is operating correctly
*
* @param timeRxs is used in the call to the updateFilter function
* @param timeTxs is used in the call to the updateFilter function
* @param seqnum is frame sequence number (so can see if there are any missing CCP frames)
*/
void CS_LSTS::log(double timeRxs, double timeTxs, int seqnum)
{

	DS_PRINT_INFO("CS_Kalman : "
		"%d %3d %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e %+15.15e\n",
		id, seqnum, timeTxs, timeRxs, dt, measuredError, x0, x1, aij, bij, drifting_fused, weight);    //add aij , bij , fused drifting and weight into the log
			
} // end log()


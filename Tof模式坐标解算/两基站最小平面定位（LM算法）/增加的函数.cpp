/*
'last' is the last_result
'raw_result' is the result before determination, maybe symmetrical to the real result.
'S' is the anchoors' coordinates matrix.  
This function is used to determin whether the result need to symmetrized.
*/
int Determin(rowvec raw_result, mat S)   
{
	double k;
	double y;
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
	
/*
'raw_result' is the result before symmetrized.
'S' is the anchoors' coordinates matrix.
This function is used to symmetrize a 2D point about a line.  
*/
rowvec Trans_Side(mat S, rowvec raw_result)
{
	rowvec trans_result = raw_result;
	
	if(S(0,0) == S(1,0))                 // The anchor line parallel to Y axis              
	{    
		trans_result(0) = 2*S(0,0) - raw_result(0);
	}
	else if(S(0,1) == S(1,1))         // The anchor line parallel to X axis
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
	return trans_result;
}



//solve_doa中增加如下语句：
	
	det_last = Determin(last_result, S);
	det_this = Determin(result, S);
	if( (det_last*det_this) < 0)      // Detect wheter the result need to symmetrize.       
	{
		result = Trans_Side(S, result);    // Symmetrize the result about the anchor line.
	}
	result(2) = set_height;



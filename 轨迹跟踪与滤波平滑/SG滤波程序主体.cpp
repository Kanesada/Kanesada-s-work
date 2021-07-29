#include <gram_savitzky_golay/gram_savitzky_golay.h>
#include <boost/circular_buffer.hpp>
#include <iostream>
using namespace std;

// Window size is 2*m+1
const size_t m = 4;
// Initial Point Smoothing (ie evaluate polynomial at first point in the window)
// Points are defined in range [-m;m]
const size_t t = 0;
// Derivation order? 0: no derivation, 1: first derivative, 2: second derivative...
const int d = 0;
boost::circular_buffer<double> window_x(2*m+1);
boost::circular_buffer<double> window_y(2*m+1);
boost::circular_buffer<double> window_z(2*m+1);

// Polynomial Order
const size_t x_order = 2;
const size_t y_order = 2;
const size_t z_order = 1;



double SG_Filter_x(double new_data, const size_t n)   // n is the x_order
{
double result;

gram_sg::SavitzkyGolayFilter xfilter(m, t, n, d);
window_x.push_back(new_data);
if(window_x.size() < 2*m+1)
{
	result = new_data;
}
else
{
	result = xfilter.filter(window_x);
	window_x[t+m] = result;
}
return result;
}



double SG_Filter_y(double new_data, const size_t n)   // n is the y_order
{
double result;

gram_sg::SavitzkyGolayFilter yfilter(m, t, n, d);
window_y.push_back(new_data);
if(window_y.size() < 2*m+1)
{
	result = new_data;
}
else
{
	result = yfilter.filter(window_y);
	window_y[t+m] = result;
}
return result;
}


double SG_Filter_z(double new_data, const size_t n)   // n is the z_order
{
double result;

gram_sg::SavitzkyGolayFilter zfilter(m, t, n, d);
window_z.push_back(new_data);
if(window_z.size() < 2*m+1)
{
	result = new_data;
}
else
{
	result = zfilter.filter(window_y);
	window_z[t+m] = result;
}
return result;
}



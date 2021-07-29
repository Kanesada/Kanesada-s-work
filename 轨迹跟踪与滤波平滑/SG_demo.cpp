#include <gram_savitzky_golay/gram_savitzky_golay.h>
#include <boost/circular_buffer.hpp>
#include <iostream>
using namespace std;

// Window size is 2*m+1
const size_t m = 4;
boost::circular_buffer<double> window_x(2*m+1);
boost::circular_buffer<double> window_y(2*m+1);
vector<double> filter_result_x(40);
vector<double> filter_result_y(40);


double SG_Filter_x(double new_data)   //Need to add interface   *Tag, data, etc
{

// Polynomial Order
const size_t n = 2;
// Initial Point Smoothing (ie evaluate polynomial at first point in the window)
// Points are defined in range [-m;m]
const size_t t = 0;
// Derivation order? 0: no derivation, 1: first derivative, 2: second derivative...
const int d = 0;
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
cout<<result<<",";
return result;
}


double SG_Filter_y(double new_data)   //Need to add interface   *Tag, data, etc
{

// Polynomial Order
const size_t n = 2;
// Initial Point Smoothing (ie evaluate polynomial at first point in the window)
// Points are defined in range [-m;m]
const size_t t = 0;
// Derivation order? 0: no derivation, 1: first derivative, 2: second derivative...
const int d = 0;
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
cout<<result<<endl;
return result;
}



int main()
{
	double x_data[40] = 
{
1.697771842466192505e+00,
2.223740481762199739e+00,
2.403425401559721841e+00,
2.787475730377451999e+00,
3.009450508686439374e+00,
3.451415658637193129e+00,
3.862031054434258870e+00,
3.983171806715033014e+00,
4.289112082477998911e+00,
4.561781656020608366e+00,
4.860957980941171819e+00,
5.227454038999738373e+00,
5.466881422913535182e+00,
5.977869921646986384e+00,
6.117724868325330689e+00,
6.521097558909278114e+00,
6.549495608084993314e+00,
6.978565071218401528e+00,
7.183502940061363695e+00,
7.374183301563155268e+00,
7.879342235277730921e+00,
7.873512052587728682e+00,
7.953779565865314538e+00,
7.956660338901272667e+00,
7.894094977838506111e+00,
7.923980120591702203e+00,
7.881824044211978908e+00,
7.873720110887153290e+00,
7.930169653976525623e+00,
7.665105037118113529e+00,
7.787514373413790381e+00,
7.703199095321857648e+00,
7.831747140081703584e+00,
7.829443864547868337e+00,
7.909154527722242811e+00,
7.606170024952064956e+00,
7.667285382246745939e+00,
7.830800431732551914e+00,
7.791383362616896235e+00,
7.775773584617216017e+00
};

double y_data[40] = 
{
1.703732110365783425e+00,
1.757778148701133913e+00,
1.993294734904166532e+00,
2.099602564773610958e+00,
1.900299667634785861e+00,
1.951131326225185258e+00,
1.968274025932890359e+00,
1.797911243423779171e+00,
1.886319914173889778e+00,
1.696076508385804926e+00,
2.027235619906603947e+00,
1.880952466827050307e+00,
1.898553407656218495e+00,
1.802368802476498111e+00,
1.852880143217625841e+00,
1.755143390196274744e+00,
1.938470939652203962e+00,
1.871383807150470258e+00,
1.788288441525292871e+00,
1.850621198748258589e+00,
1.945972600057446211e+00,
2.140125990932341349e+00,
2.425255998173669880e+00,
2.668196724699992473e+00,
2.938408857909811989e+00,
3.143047699661112304e+00,
3.524417954537535369e+00,
3.908168863899904810e+00,
4.163193264888079170e+00,
4.559410495834526955e+00,
4.912806842795315099e+00,
5.161897784777379705e+00,
5.443066153432718046e+00,
5.655803356078370037e+00,
6.026599112005740722e+00,
6.254994893669071132e+00,
6.753293675797245932e+00,
6.718212557896033665e+00,
7.340794283150849253e+00,
7.358940302119123977e+00	
};

	register int i;
	for(i=0;i<40;i++)
	{
		filter_result_x.push_back(SG_Filter_x(x_data[i]));
		filter_result_y.push_back(SG_Filter_y(y_data[i]));
	}
	//cout<<filter_result<<endl;
	return 0;
}



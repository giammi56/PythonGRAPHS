#include "ADC_helper_functions.h"
#include "math.h"
#include <iostream>

#define PI	3.14159265359
#define TWOPI	(2.0*PI)

int linear_gaus::fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max,  double &x_0, double &sigma, double &y_0){
	if(data_max-data_min<=0 && data_max>=nbrPoints-1)return 0;
	
	int array_size = data_max - data_min - 1;

	double *x = new double[array_size];
	double *y = new double[array_size];
	double *w = new double[array_size];

	int j=0;
	for(int i = data_min; i < data_max-1; i++){
		x[j]=i+1;
		y[j]=log(data[i]/data[i+2]);
		w[j]=1./(1/data[i+2]+1/data[i]);

		//std::cout << "\n" << x[j] << " , " << y[j];
		j++;
	}

	double m = 0.;
	double c = 0.;

	linear_weighted(x,y,w,array_size,m,c);

	sigma	= sqrt(2/m);
	x_0		= -c/m;

	get_y0(data, nbrPoints, data_min, data_max, m, c, x_0, sigma, y_0);

	delete[] x;
	delete[] y;
	delete[] w;
	
	if(m>0)
		return 1;
	return 0;
}
int linear_gaus::fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max,  double &x_0, double &sigma){
	if(data_max-data_min<=0 && data_max>=nbrPoints-1)return 0;

	int array_size = data_max - data_min - 1;

	double *x = new double[array_size];
	double *y = new double[array_size];
	double *w = new double[array_size];

	int j=0;
	for(int i = data_min; i < data_max-1; i++){
		x[j]=i+1;
		y[j]=log(data[i]/data[i+2]);
		w[j]=1./(1./data[i+2]+1/data[i]);

		//std::cout << "\n" << x[j] << " , " << y[j];
		j++;
	}

	double m = 0.;
	double c = 0.;

	linear_weighted(x,y,w,array_size,m,c);

	delete[] x;
	delete[] y;
	delete[] w;

	sigma	= sqrt(2/m);
	x_0		= -c/m;

	if(m>0)
		return 1;
	return 0;
}
double linear_gaus::gaus_fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max){
	if (data_max - data_min <= 0 && data_max >= nbrPoints - 1)return 0;

	int array_size = data_max - data_min - 1;

	double *x = new double[array_size];
	double *y = new double[array_size];
	double *w = new double[array_size];

	int j = 0;
	for (int i = data_min; i < data_max - 1; i++) {
		x[j] = i + 1;
		y[j] = log(data[i] / data[i + 2]);
		w[j] = 1. / (1. / data[i + 2] + 1 / data[i]);
		j++;
	}

	double m = 0.;
	double c = 0.;

	linear_weighted(x, y, w, array_size, m, c);

	delete[] x;
	delete[] y;
	delete[] w;

	if (m<0)
		return -1;
	
	return -c / m;
}
void linear_gaus::linear_weighted(double *x, double *y, double *w, __int32 nbrPoints,double &m, double &c){
	double SumW=0., SumWX=0., SumWY=0., SumWXsq=0., SumWXY=0.;
	for(int i=0; i<nbrPoints; i++){
		SumW    += w[i];
		SumWX   += w[i] * x[i];
		SumWY   += w[i] * y[i];
		SumWXsq += w[i] * x[i] * x[i];
		SumWXY  += w[i] * x[i] * y[i];
	}

	//-- normal equations (Numerical Methods with C++ Programming - Nita H. Shah):	--//
	//-- SumW  * m + SumWX   * c = SumWY											--//
	//-- SumWX * m + SumWXsq * c = SumWXY											--//

	m = (SumWXY-SumWX*SumWY/SumW)/(SumWXsq-SumWX*SumWX/SumW); 
	c = (SumWY-SumWX*m)/SumW;

	return;
}
void linear_gaus::get_y0(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max, double &m, double &c, double &x_0, double &sigma, double &y_0){
	double SumX=0.,SumXsq=0.; double sigma_sq=0.;
	int N = data_max-data_min-2;

	//-- P. R. Bevington, Data reduction and error analysis for the	--//
	//-- physical science (McGraw-Hill, New York, 1969) Ch. 6.		--//

	for(int i = data_min; i < data_max-1; i++){
		SumX		+= i+1.;
		SumXsq		+= (i+1.)*(i+1.);
		sigma_sq	+= (data[i+1]-c-m*(i+1.))*(data[i+1]-c-m*(i+1.));
	}

	double delta = N*SumXsq - (SumX*SumX);
	sigma_sq = 1./(N-2.)*sigma_sq;

	double SumWplus=0., SumW=0., W=0.;
	double sigma4 = 1./(sigma*sigma*sigma*sigma);
	double inv_sigma2 = 1./sigma/sigma;
	for(int i = data_min; i < data_max-1; i++) {
		W = 1./(1./data[i+1]+((i+1.)-x_0)*((i+1.)-x_0)*sigma4*(delta*delta+((i+1)-x_0)*((i+1.)-x_0)*inv_sigma2*sigma_sq));                                 
		SumW	 += W;
		SumWplus += W*log(data[i+1]+((i+1.)-x_0)*((i+1.)-x_0)*0.5*inv_sigma2);
	}
	
	y_0 = exp(SumWplus/SumW);	
	return;
}

double parabola_spline::parabola_fit(double y1, double y2, double y3, double x_pos)
{
	double a = 0.5*(y1+y3) - y2;
	double b = 2*y2 - 1.5*y1 - 0.5*y3;
	return (a*x_pos + b)*x_pos + y1;
}
double parabola_spline::get_y_at_x(int * y, __int32 number_of_points, double x)
{
	if (x < 0. || x > number_of_points) return 1.e200;
	if (number_of_points < 3) return 1.e200;

	if (x <= 1.) {
		return parabola_fit(y[0], y[1], y[2], x);
	}
	if (x >= number_of_points-2) {
		return parabola_fit(y[number_of_points-3], y[number_of_points-2], y[number_of_points-1], x-number_of_points+3);
	}
	__int32 index = __int32(x);
	double dx = x-index;
	double yi = y[index];
	double yip1 = y[index+1];
	double y1 = parabola_fit(y[index-1], yi, yip1, dx+1.);
	double y2 = parabola_fit(yi, yip1, y[index+2], dx);

	return y1*(1.-dx) + y2*dx;
}
double parabola_spline::get_y_at_x(double * y, __int32 number_of_points, double x)
{
	if (x < 0. || x > number_of_points) return 1.e200;
	if (number_of_points < 3) return 1.e200;

	if (x <= 1.) {
		return parabola_fit(y[0], y[1], y[2], x);
	}
	if (x >= number_of_points-2) {
		return parabola_fit(y[number_of_points-3], y[number_of_points-2], y[number_of_points-1], x-number_of_points+3);
	}
	__int32 index = __int32(x);
	double dx = x-index;
	double yi = y[index];
	double yip1 = y[index+1];
	double y1 = parabola_fit(y[index-1], yi, yip1, dx+1.);
	double y2 = parabola_fit(yi, yip1, y[index+2], dx);

	return y1*(1.-dx) + y2*dx;
}
double parabola_spline::get_parabola_vertex(double y1, double y2, double y3) {
	double a = 0.5*(y1 + y3) - y2;
	double b = 2 * y2 - 1.5*y1 - 0.5*y3;
	return y1 - (b*b*0.25) / a;
}
double parabola_spline::get_parabola_vertex(double y1, double y2, double y3, double y4, double dx) {
	return get_parabola_vertex(y1,y2,y3)*(1.-dx) + get_parabola_vertex(y2,y3,y4)*dx;
}

void bezier_class::bezier3_make_curve_parameters(double *x, double *y, double weight)
{
	if (weight == 0.) weight = 1.;

	xA = x[0];
	xB = x[2];
	yA = y[0];
	yB = y[2];

	bez_xa2 = x[0]+x[2]-2.*x[1]*weight;
	bez_xa1 = 2.*(x[1]*weight-x[0]);
	bez_xa0 = x[0];

	bez_ya2 = y[0]+y[2]-2.*y[1]*weight;
	bez_ya1 = 2.*(y[1]*weight-y[0]);
	bez_ya0 = y[0];

	w2 = 2.-2.*weight;
	w1 = 2.*weight-2.;
	w0 = 1.;
}
void bezier_class::bezier4_make_curve_parameters(double *x, double *y, double weight)
{
	if (weight == 0.) weight = 1.;

	xA = x[0];
	xB = x[3];
	yA = y[0];
	yB = y[3];

	bez_xa3 = x[3]-x[0] +3.*weight*(x[1]-x[2]);
	bez_xa2 = 3.*(x[0]+x[2]*weight) - 6.*x[1]*weight;
	bez_xa1 = 3.*(x[1]*weight - x[0]);
	bez_xa0 = x[0];

	bez_ya3 = y[3]-y[0]  +3.*weight*(y[1]-y[2]);
	bez_ya2 = 3.*(y[0]+y[2]*weight) - 6.*y[1]*weight;
	bez_ya1 = 3.*(y[1]*weight - y[0]);
	bez_ya0 = y[0];

	w3 = 0.;
	w2 = 3.-3.*weight;
	w1 = 3.*weight-3.;
	w0 = 1;
}
void bezier_class::calculate_bezier3_point(double t, double &x, double &y)
{
	double tt =  t*t;
	double weight = tt*w2 + t*w1 + w0;

	x = (bez_xa2*tt + bez_xa1*t + bez_xa0) / weight;
	y = (bez_ya2*tt + bez_ya1*t + bez_ya0) / weight;
}
void bezier_class::calculate_bezier4_point(double t, double &x, double &y)
{
	double tt =  t*t;
	double ttt = tt*t;
	double weight = ttt*w3 + tt*w2 + t*w1 + w0;

	x = (bez_xa3*ttt + bez_xa2*tt + bez_xa1*t + bez_xa0)/weight;
	y = (bez_ya3*ttt + bez_ya2*tt + bez_ya1*t + bez_ya0)/weight;
}
double bezier_class::get_x_at_y_bezier3(double y_thresh)
{
	double y0 = yA;
	double y1 = yB;
	double x0 = xA;
	double x1 = xB;

	double t0 = 0.;
	double t1 = 1.;

	double tn, yn, x_new;

	if (y0 > y1) {
		if (y_thresh > y0 || y_thresh < y1) {
			last_command_was_successful = false;
			return 0.;
		}
		while (true) {
			tn = (t1+t0)*0.5;
			calculate_bezier3_point(tn, x_new, yn);
			if (yn > y_thresh) {
				t0 = tn; y0 = yn; x0 = x_new;
			} else {
				t1 = tn; y1 = yn; x1 = x_new;
			}
			if (x1 - x0 < 80.) break;
		}
	} else {
		if (y_thresh < y0 || y_thresh > y1) {
			last_command_was_successful = false;
			return 0.;
		}
		while (true) {
			tn = (t1+t0)*0.5;
			calculate_bezier3_point(tn, x_new, yn);
			if (yn < y_thresh) {
				t0 = tn; y0 = yn; x0 = x_new;
			} else {
				t1 = tn; y1 = yn; x1 = x_new;
			}
			if (x1 - x0 < 80.) break;
		}
	}

	last_command_was_successful = true;
	return x0+(x1-x0)*(y_thresh-y0)/(y1-y0);
}
double bezier_class::get_x_at_y_bezier4(double y_thresh)
{
	double y0 = yA;
	double y1 = yB;
	double x0 = xA;
	double x1 = xB;

	double t0 = 0.;
	double t1 = 1.;

	double tn, yn, x_new;


	if (y0 > y1) {
		if (y_thresh > y0 || y_thresh < y1) {
			last_command_was_successful = false;
			return 0.;
		}
		while (true) {
			tn = (t1+t0)*0.5;
			calculate_bezier4_point(tn, x_new, yn);
			if (yn > y_thresh) {
				t0 = tn; y0 = yn; x0 = x_new;
			} else {
				t1 = tn; y1 = yn; x1 = x_new;
			}
			if (x1 - x0 < 80.) break;
		}
	} else {
		if (y_thresh < y0 || y_thresh > y1) {
			last_command_was_successful = false;
			return 0.;
		}
		while (true) {
			tn = (t1+t0)*0.5;
			calculate_bezier4_point(tn, x_new, yn);
			if (yn < y_thresh) {
				t0 = tn; y0 = yn; x0 = x_new;
			} else {
				t1 = tn; y1 = yn; x1 = x_new;
			}
			if (x1 - x0 < 80.) break;
		}
	}

	last_command_was_successful = true;
	return x0+(x1-x0)*(y_thresh-y0)/(y1-y0);
}

void Trilinear::resetErrorMatrix() {
	for (int x = 0; x< x_bins; x++) {
		for (int y = 0; y < y_bins; y++) {
			for (int z = 0; z < z_bins; z++) {
				m_matrix[x][y][z].error = 0;
			}
		}
	}
}
void Trilinear::errorMatrix() {
	for (int x = 0; x< x_bins; x++) {
		for (int y = 0; y < y_bins; y++) {
			for (int z = 0; z < z_bins; z++) {
				if (m_matrix[x][y][z].nbr)
					m_matrix[x][y][z].error = sqrt(m_matrix[x][y][z].error / m_matrix[x][y][z].nbr);
			}
		}
	}
	return;
}
void Trilinear::resetNbrMatrix() {
	for (int x = 0; x< x_bins; x++) {
		for (int y = 0; y < y_bins; y++) {
			for (int z = 0; z < z_bins; z++) {
				m_matrix[x][y][z].nbr = 0;
			}
		}
	}
}
int Trilinear::fillErrorMatrix(double x, double y, double z, double error) {
	if (!m_matrix)return 0;

	int x_bin = int((x - x_min) / x_range*x_bins + 0.5);
	int y_bin = int((y - y_min) / y_range*y_bins + 0.5);
	int z_bin = int((z - z_min) / z_range*z_bins + 0.5);

	if (x_bin >= 0 && x_bin < x_bins && y_bin >= 0 && y_bin < y_bins && z_bin >= 0 && z_bin < z_bins) {
		m_matrix[x_bin][y_bin][z_bin].error += (m_matrix[x_bin][y_bin][z_bin].val - error) * (m_matrix[x_bin][y_bin][z_bin].val - error);
		m_matrix[x_bin][y_bin][z_bin].nbr++;
		return 0;
	}
	return 0;
}

Trilinear::Trilinear() {
	x_bins = 0;
	x_min = 0.;
	x_max = 0.;
	x_delta = 0.;
	x_range = 0.;

	y_bins = 0;
	y_min = 0.;
	y_max = 0.;
	y_delta = 0.;
	y_range = 0.;

	z_bins = 0;
	z_min = 0.;
	z_max = 0.;
	z_delta = 0.;
	z_range = 0.;

	matrix_volume = 0.;

	m_matrix = nullptr;

	yz_cnt_array = nullptr;
}
Trilinear::~Trilinear() {
	// BACON
	//for (int x = 0; x < x_bins; x++) {
	//	for (int y = 0; y < y_bins; y++) {
	//		if (!m_matrix[x][y]) {
	//			delete m_matrix[x][y];
	//		}
	//	}
	//}
}
void Trilinear::createMatrix(int x_bins_in, double x_min_in, double x_max_in, int y_bins_in, double y_min_in, double y_max_in, int z_bins_in, double z_min_in, double z_max_in) {
	x_bins = x_bins_in;
	x_min = x_min_in;
	x_max = x_max_in;
	x_delta = (x_max - x_min) / x_bins;
	x_range = x_max - x_min;

	y_bins = y_bins_in;
	y_min = y_min_in;
	y_max = y_max_in;
	y_delta = (y_max - y_min) / y_bins;
	y_range = y_max - y_min;

	z_bins = z_bins_in;
	z_min = z_min_in;
	z_max = z_max_in;
	z_delta = (z_max - z_min) / z_bins;
	z_range = z_max - z_min;

	matrix_volume = x_delta * y_delta * z_delta;

	m_matrix = new Entry**[x_bins];
	for (int i = 0; i< x_bins; i++) {
		m_matrix[i] = new Entry*[y_bins];
		for (int j = 0; j < y_bins; j++) {
			m_matrix[i][j] = nullptr; // new Entry[z_bins];
		}
	}

	//for (int x = 0; x< x_bins; x++) {
	//	for (int y = 0; y < y_bins; y++) {
	//		for (int z = 0; z < z_bins; z++) {
	//			m_matrix[x][y][z].error = 0;
	//			m_matrix[x][y][z].nbr = 0;
	//			m_matrix[x][y][z].val = 0.;
	//		}
	//	}
	//}


	yz_cnt_array = new int*[y_bins];
	for (int i = 0; i < y_bins; i++) {
		yz_cnt_array[i] = new int[z_bins];
		for (int j = 0; j < z_bins; j++) {
			yz_cnt_array[i][j] = 0;
		}
	}
	return;
}
void Trilinear::resetMatrix() {
	for (int x = 0; x< x_bins; x++) {
		for (int y = 0; y < y_bins; y++) {
			for (int z = 0; z < z_bins; z++) {
				m_matrix[x][y][z].error = 0;
				m_matrix[x][y][z].nbr = 0;
				m_matrix[x][y][z].val = 0.;
			}
		}
	}

	x_bins = 0;
	x_min = 0.;
	x_max = 0.;
	x_delta = 0.;
	y_bins = 0;
	y_min = 0.;
	y_max = 0.;
	y_delta = 0.;
	z_bins = 0;
	z_min = 0.;
	z_max = 0.;
	z_delta = 0.;
}
int Trilinear::fillMatrix(double x, double y, double z, double val) {
	if (!m_matrix)return 0;

	int x_bin = int((x - x_min) / x_range*x_bins - 1 + 0.5);
	int y_bin = int((y - y_min) / y_range*y_bins - 1 + 0.5);
	int z_bin = int((z - z_min) / z_range*z_bins - 1 + 0.5);

	if (x_bin < 0 || x_bin >= x_bins || y_bin < 0 || y_bin >= y_bins || z_bin < 0 || z_bin >= z_bins)
		return 0;

	if (!m_matrix[x_bin][y_bin]) {
		m_matrix[x_bin][y_bin] = new Entry[z_bins];
		for (int z = 0; z < z_bins; z++) {
			m_matrix[x_bin][y_bin][z].error = 0;
			m_matrix[x_bin][y_bin][z].nbr = 0;
			m_matrix[x_bin][y_bin][z].val = 0.;
		}
	}

	if (x_bin >= 0 && x_bin < x_bins && y_bin >= 0 && y_bin < y_bins && z_bin >= 0 && z_bin < z_bins) {
		m_matrix[x_bin][y_bin][z_bin].val += val;
		m_matrix[x_bin][y_bin][z_bin].nbr++;

		yz_cnt_array[y_bin][z_bin]++;
		return 0;
	}
	return 0;
}
void Trilinear::averageMatrix() {
	for (int x = 0; x< x_bins; x++) {
		for (int y = 0; y < y_bins; y++) {
			if (!m_matrix[x][y]) {
				for (int z = 0; z < z_bins; z++) {
					m_matrix[x][y][z].val = m_matrix[x][y][z].val / m_matrix[x][y][z].nbr;
				}
			}
		}
	}
	return;
}
double Trilinear::interpolate3D(double x, double y, double z) {
	if (!m_matrix)return 0;

	int x_bin = int((x - x_min) / x_range*x_bins + 0.5);
	int y_bin = int((y - y_min) / y_range*y_bins + 0.5);
	int z_bin = int((z - z_min) / z_range*z_bins + 0.5);

	if (x_bin >= 0 && x_bin < x_bins - 1 && y_bin >= 0 && y_bin < y_bins - 1 && z_bin >= 0 && z_bin < z_bins - 1) {
		if (!m_matrix[x_bin][y_bin] || !m_matrix[x_bin][y_bin][z_bin].val)
			return 0;

		double x_r = x - (x_min + x_bin * x_delta);
		double x_i = x_delta - x_r;
		double y_r = y - (y_min + y_bin * y_delta);
		double y_i = y_delta - y_r;
		double z_r = z - (z_min + z_bin * z_delta);
		double z_i = z_delta - z_r;

		//double error = (m_matrix[x_bin + 1][y_bin + 1][z_bin + 1].error * x_r * y_r * z_r
		//	+ m_matrix[x_bin][y_bin + 1][z_bin + 1].error * x_i * y_r * z_r
		//	+ m_matrix[x_bin][y_bin + 1][z_bin].error * x_i * y_r * z_i
		//	+ m_matrix[x_bin][y_bin][z_bin + 1].error * x_i * y_i * z_r
		//	+ m_matrix[x_bin + 1][y_bin + 1][z_bin].error * x_r * y_r * z_i
		//	+ m_matrix[x_bin][y_bin][z_bin].error * x_i * y_i * z_i
		//	+ m_matrix[x_bin + 1][y_bin][z_bin].error * x_r * y_i * z_i
		//	+ m_matrix[x_bin + 1][y_bin][z_bin + 1].error * x_r * y_i * z_r) / matrix_volume;

		//if (m_matrix[x_bin + 1][y_bin + 1][z_bin + 1].error > 0.4 ||
		//	m_matrix[x_bin][y_bin + 1][z_bin + 1].error > 0.4 ||
		//	m_matrix[x_bin][y_bin + 1][z_bin].error > 0.4 ||
		//	m_matrix[x_bin][y_bin][z_bin + 1].error > 0.4 ||
		//	m_matrix[x_bin + 1][y_bin + 1][z_bin].error > 0.4 ||
		//	m_matrix[x_bin][y_bin][z_bin].error > 0.4 ||
		//	m_matrix[x_bin + 1][y_bin][z_bin].error > 0.4 ||
		//	m_matrix[x_bin + 1][y_bin][z_bin + 1].error > 0.4)return -1.;

		//if (error>0.4)return -1.;

		//if a relevant neighbour cell is emty fall into emergency mode:
		if (!m_matrix[x_bin][y_bin + 1] || !m_matrix[x_bin + 1][y_bin] || !m_matrix[x_bin + 1][y_bin + 1]) {
			return m_matrix[x_bin][y_bin][z_bin].val + (m_matrix[x_bin][y_bin][z_bin + 1].val - m_matrix[x_bin][y_bin][z_bin].val)*(z - z_bin) / z_delta;
		}

		return  /*double v_111	=*/  (m_matrix[x_bin + 1][y_bin + 1][z_bin + 1].val * x_r * y_r * z_r	/*;*/
			/*double v_011	=*/ + m_matrix[x_bin][y_bin + 1][z_bin + 1].val * x_i * y_r * z_r		/*;+/
																									/*double v_001	=*/ + m_matrix[x_bin][y_bin][z_bin + 1].val * x_i * y_i * z_r			/*;+/
																									/*double v_110	=*/ + m_matrix[x_bin + 1][y_bin + 1][z_bin].val * x_r * y_r * z_i		/*;+/
																									/*double v_000	=*/ + m_matrix[x_bin][y_bin][z_bin].val * x_i * y_i * z_i				/*;+/
																									/*double v_100	=*/ + m_matrix[x_bin + 1][y_bin][z_bin].val * x_r * y_i * z_i			/*;+/
																									/*double v_101	=*/ + m_matrix[x_bin + 1][y_bin][z_bin + 1].val * x_r * y_i * z_r) / matrix_volume;
	}
	return 0.;
}
double Trilinear::GetVal(int x, int y, int z) {
	if (x >= 0 && x < this->x_bins && y >= 0 && y < this->y_bins && z >= 0 && z < this->z_bins) {
		if (m_matrix[x][y])
			return m_matrix[x][y][z].val;
	}
	return 0.;
}
void Trilinear::SetVal(int x, int y, int z, double val) {
	if (x >= 0 && x < this->x_bins && y >= 0 && y < this->y_bins && z >= 0 && z < this->z_bins) {
		if (m_matrix[x][y])
			m_matrix[x][y][z].val = val;
	}
	return;
}

Bilinear::Bilinear() {
	x_bins = 0;
	x_min = 0.;
	x_max = 0.;
	x_delta = 0.;
	x_range = 0.;

	y_bins = 0;
	y_min = 0.;
	y_max = 0.;
	y_delta = 0.;
	y_range = 0.;

	matrix_volume = 0.;

	m_matrix = nullptr;
}
Bilinear::~Bilinear() {
	if (m_matrix) {
		for (int x = 0; x < x_bins; x++) {
			delete m_matrix[x];
			m_matrix[x] = 0;
		}
		delete m_matrix;
		m_matrix = 0;
	}
}
void Bilinear::createMatrix(int x_bins_in, double x_min_in, double x_max_in, int y_bins_in, double y_min_in, double y_max_in) {
	x_bins = x_bins_in;
	x_min = x_min_in;
	x_max = x_max_in;
	x_delta = (x_max - x_min) / x_bins;
	x_range = x_max - x_min;

	y_bins = y_bins_in;
	y_min = y_min_in;
	y_max = y_max_in;
	y_delta = (y_max - y_min) / y_bins;
	y_range = y_max - y_min;

	matrix_volume = x_delta * y_delta;

	m_matrix = new Entry*[x_bins];
	for (int i = 0; i< x_bins; i++) {
		m_matrix[i] = new Entry[y_bins];
		for (int j = 0; j < y_bins; j++) {
			m_matrix[i][j].val = 0;
			m_matrix[i][j].error = 0;
			m_matrix[i][j].nbr = 0;
		}
	}

	return;
}
int Bilinear::fillMatrix(double x, double y, double val) {
	if (!m_matrix)return 0;

	int x_bin = int((x - x_min) / x_range*x_bins - 1 + 0.5);
	int y_bin = int((y - y_min) / y_range*y_bins - 1 + 0.5);

	if (x_bin < 0 || x_bin >= x_bins || y_bin < 0 || y_bin >= y_bins)
		return 0;

	//if (!m_matrix[x_bin][y_bin]) {
	//	m_matrix[x_bin][y_bin] = new Entry[z_bins];
	//	for (int z = 0; z < z_bins; z++) {
	//		m_matrix[x_bin][y_bin][z].error = 0;
	//		m_matrix[x_bin][y_bin][z].nbr = 0;
	//		m_matrix[x_bin][y_bin][z].val = 0.;
	//	}
	//}

	if (x_bin >= 0 && x_bin < x_bins && y_bin >= 0 && y_bin < y_bins) {
		m_matrix[x_bin][y_bin].val += val;
		m_matrix[x_bin][y_bin].nbr++;
		return 0;
	}
	return 0;
}
double Bilinear::GetVal(int x, int y) {
	if (x >= 0 && x < this->x_bins && y >= 0 && y < this->y_bins) {
		if (m_matrix[x])
			return m_matrix[x][y].val;
	}
	return 0.;
}
void Bilinear::SetVal(int x, int y, double val) {
	if (x >= 0 && x < this->x_bins && y >= 0 && y < this->y_bins) {
		if (m_matrix[x])
			m_matrix[x][y].val = val;
	}
	return;
}
void Bilinear::averageMatrix() {
	for (int x = 0; x< x_bins; x++) {
		if (m_matrix[x]) {
			for (int y = 0; y < y_bins; y++) {
				m_matrix[x][y].val = m_matrix[x][y].val / m_matrix[x][y].nbr;
			}
		}
	}
	return;
}
double Bilinear::interpolate(double x, double y) {
	if (!m_matrix)return 0;

	int x_bin = int((x - x_min) / x_range*x_bins + 0.5);
	int y_bin = int((y - y_min) / y_range*y_bins + 0.5);

	if (x_bin >= 0 && x_bin < x_bins - 1 && y_bin >= 0 && y_bin < y_bins - 1) {
		if (!m_matrix[x_bin])
			return 0;

		double x_r = x - (x_min + x_bin * x_delta);		// x1------x
		double x_i = x_delta - x_r;						//         x-----------x2
		double y_r = y - (y_min + y_bin * y_delta);		// y1---------y
		double y_i = y_delta - y_r;						//            y--------y2

		double value = 0., volume = 0.;

		if (m_matrix[x_bin][y_bin].nbr) {
			value += m_matrix[x_bin][y_bin].val * x_i * y_i;
			volume += x_i * y_i;
		}
		if (m_matrix[x_bin + 1][y_bin].nbr) {
			value += m_matrix[x_bin + 1][y_bin].val * x_r * y_i;
			volume += x_r * y_i;
		}
		if (m_matrix[x_bin + 1][y_bin + 1].nbr) {
			value += m_matrix[x_bin + 1][y_bin + 1].val * x_r * y_r;
			volume += x_r * y_r;
		}
		if (m_matrix[x_bin][y_bin + 1].nbr) {
			value += m_matrix[x_bin][y_bin + 1].val * x_i * y_r;
			volume += x_i * y_r;
		}

		if (volume > 0.)
			return value / volume;


		//return ( m_matrix[x_bin  ][y_bin  ].val * x_i * y_i +
		//		 m_matrix[x_bin+1][y_bin  ].val * x_r * y_i +
		//		 m_matrix[x_bin+1][y_bin+1].val * x_r * y_r +
		//		 m_matrix[x_bin  ][y_bin+1].val * x_i * y_r) / matrix_volume;
	}
	return 0.;
}

int FFT_class::getNbrPointsFFT(int nbrPoints) {
	return (int)pow(2.0, ceil(log((double)nbrPoints) / log(2.0)));
}
void FFT_class::fft(double *data_in, int nbrPoints, double *&data_fft_out, int nbrPointsFFT)
{
	int n, mmax, m, j, istep, i, NFFT;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	/* calculate NFFT as the next higher power of 2 >= nbrPoints */
	if (nbrPointsFFT > nbrPoints)
		NFFT = getNbrPointsFFT(nbrPointsFFT);
	else
		NFFT = getNbrPointsFFT(nbrPoints);

	/* allocate memory for NFFT complex numbers (note the +1) */
	data_fft_out = new double[2 * NFFT + 1];

	/* Storing x(n) in a complex array to make it work with four1.
	This is needed even though x(n) is purely real in this case. */
	for (i = 0; i<nbrPoints; i++)
	{
		data_fft_out[2 * i + 1] = data_in[i];
		data_fft_out[2 * i + 2] = 0.0;
	}
	/* pad the remainder of the array with zeros (0 + 0 j) */
	for (i = nbrPoints; i<NFFT; i++)
	{
		data_fft_out[2 * i + 1] = 0.0;
		data_fft_out[2 * i + 2] = 0.0;
	}

	// Do the FFT..
	n = NFFT << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			tempr = data_fft_out[j];     data_fft_out[j] = data_fft_out[i];     data_fft_out[i] = tempr;
			tempr = data_fft_out[j + 1]; data_fft_out[j + 1] = data_fft_out[i + 1]; data_fft_out[i + 1] = tempr;
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2 * mmax;
		theta = TWOPI / mmax;
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data_fft_out[j] - wi*data_fft_out[j + 1];
				tempi = wr*data_fft_out[j + 1] + wi*data_fft_out[j];
				data_fft_out[j] = data_fft_out[i] - tempr;
				data_fft_out[j + 1] = data_fft_out[i + 1] - tempi;
				data_fft_out[i] += tempr;
				data_fft_out[i + 1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}

	data_fft_out[0] = NFFT;
}
void FFT_class::ifft(double *data_fft_in, int nbrPoints, double *&data_out, double *&real_out, double *&imag_out, int nbrPointsFFT) {
	int n, mmax, m, j, istep, i, NFFT;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	/* calculate NFFT as the next higher power of 2 >= nbrPoints */
	// zero padding
	double *data_fft_tmp = 0;
	if (nbrPointsFFT > nbrPoints) {
		NFFT = getNbrPointsFFT(nbrPointsFFT);

		(data_fft_tmp) = new double[2 * NFFT + 1];

		for (i = 0; i < int(data_fft_in[0] / 2); i++) {
			data_fft_tmp[2 * i + 1] = data_fft_in[2 * i + 1];
			data_fft_tmp[2 * i + 2] = data_fft_in[2 * i + 2];
		}
		for (i = int(data_fft_in[0] / 2); i < int(NFFT - data_fft_in[0] / 2); i++) {
			data_fft_tmp[2 * i + 1] = 0;
			data_fft_tmp[2 * i + 2] = 0;
		}
		for (i = int(NFFT - data_fft_in[0] / 2); i < NFFT; i++) {
			data_fft_tmp[2 * i + 1] = data_fft_in[2 * i + 1];
			data_fft_tmp[2 * i + 2] = data_fft_in[2 * i + 2];
		}
	}
	else {
		NFFT = getNbrPointsFFT(nbrPoints);
		data_fft_tmp = new double[2 * NFFT + 1];
		data_fft_tmp[0] = data_fft_in[0];
		for (i = 0; i < NFFT; i++)
		{
			data_fft_tmp[2 * i + 1] = data_fft_in[2 * i + 1];
			data_fft_tmp[2 * i + 2] = data_fft_in[2 * i + 2];
		}
	}
	
	// Do the inverse FFT..
	n = NFFT << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			tempr = data_fft_tmp[j];     data_fft_tmp[j] = data_fft_tmp[i];     data_fft_tmp[i] = tempr;
			tempr = data_fft_tmp[j + 1]; data_fft_tmp[j + 1] = data_fft_tmp[i + 1]; data_fft_tmp[i + 1] = tempr;
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2 * mmax;
		theta = -TWOPI / mmax;
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr*data_fft_tmp[j] - wi*data_fft_tmp[j + 1];
				tempi = wr*data_fft_tmp[j + 1] + wi*data_fft_tmp[j];
				data_fft_tmp[j] = data_fft_tmp[i] - tempr;
				data_fft_tmp[j + 1] = data_fft_tmp[i + 1] - tempi;
				data_fft_tmp[i] += tempr;
				data_fft_tmp[i + 1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}

	data_out = new double[nbrPoints];
	real_out = new double[nbrPoints];
	imag_out = new double[nbrPoints];

	for (i = 0; i < nbrPoints; i++) {
		data_out[i] = /*data_fft_tmp[2 * i + 1]*/sqrt(data_fft_tmp[2 * i + 1]* data_fft_tmp[2 * i + 1] + data_fft_tmp[2 * i + 2]* data_fft_tmp[2 * i + 2]) / NFFT;
		real_out[i] = data_fft_tmp[2 * i + 1];
		imag_out[i] = data_fft_tmp[2 * i + 2];
		//data_fft_tmp[2 * i + 2] /= NFFT;
	}

	delete[] data_fft_tmp;
}
void FFT_class::fft_shift(double *data_fft, int nbrPoints = 0, int nbrPointsFFT) {
	if (nbrPoints > 0 && nbrPoints > nbrPointsFFT)
		nbrPointsFFT = getNbrPointsFFT(nbrPoints);

	double real_tmp = 0., imag_tmp = 0.;
	for (int i = 0; i < nbrPointsFFT / 2; i++) {
		real_tmp = data_fft[2 * i + 1];
		imag_tmp = data_fft[2 * i + 2];

		data_fft[2 * i + 1] = data_fft[2 * i + 1 + (nbrPointsFFT)];
		data_fft[2 * i + 2] = data_fft[2 * i + 2 + (nbrPointsFFT)];

		data_fft[2 * i + 1 + (nbrPointsFFT)] = real_tmp;
		data_fft[2 * i + 2 + (nbrPointsFFT)] = imag_tmp;
	}
}
void FFT_class::fft_shift_undo(double *data_fft_shifted, int nbrPoints = 0, int nbrPointsFFT) {
	if (nbrPoints > 0 && nbrPoints > nbrPointsFFT)
		nbrPointsFFT = getNbrPointsFFT(nbrPoints);

	double real_tmp = 0., imag_tmp = 0.;
	for (int i = 0; i < nbrPointsFFT / 2; i++) {
		real_tmp = data_fft_shifted[2 * i + 1];
		imag_tmp = data_fft_shifted[2 * i + 2];

		data_fft_shifted[2 * i + 1] = data_fft_shifted[2 * i + 1 + (nbrPointsFFT / 2)];
		data_fft_shifted[2 * i + 2] = data_fft_shifted[2 * i + 2 + (nbrPointsFFT / 2)];

		data_fft_shifted[2 * i + 1 + (nbrPointsFFT / 2)] = real_tmp;
		data_fft_shifted[2 * i + 2 + (nbrPointsFFT / 2)] = imag_tmp;
	}
}
void FFT_class::getAmplitude(double *data_fft_in, int nbrPoints, double *&data_Aftt_out, int nbrPointsFFT) {
	//if (nbrPoints > 0 && nbrPointsFFT != 0 && nbrPoints > nbrPointsFFT)
	//	nbrPointsFFT = getNbrPointsFFT(nbrPoints);

	nbrPointsFFT = int(data_fft_in[0] + 0.1);	//getNbrPointsFFT(nbrPoints);
	data_Aftt_out = new double[nbrPointsFFT];

	for (int i = 0; i < nbrPointsFFT; i++) {
		data_Aftt_out[i] = sqrt(data_fft_in[2 * i + 1] * data_fft_in[2 * i + 1] + data_fft_in[2 * i + 2] * data_fft_in[2 * i + 2]);
	}
}
void FFT_class::fft_divide(double *fft_A, double *fft_B) {
	int PointsA = int(fft_A[0] + 0.1);	//getNbrPointsFFT(nbrPoints);
	int PointsB = int(fft_B[0] + 0.1);

	if (PointsA == PointsB) {
		double Are, Aim, Bre, Bim, nenner;

		for (int i = 0; i < PointsA; i++) {
			Are = fft_A[2 * i + 1];
			Aim = fft_A[2 * i + 2];
			Bre = fft_B[2 * i + 1];
			Bim = fft_B[2 * i + 2];

			nenner = Bre*Bre + Bim*Bim;
			if (nenner != 0) {
				fft_A[2 * i + 1] = (Are*Bre + Aim*Bim) / nenner;
				fft_A[2 * i + 2] = (Aim*Bre - Are*Bim) / nenner;
			}
		}
	}
}

void lma2root::ReadLMAHeader() {
	//flag that will show wether something has changed
	bool changed = false;

	//--read all the values in the header for the Event--//
	short  tempS;
	long   tempL;
	double tempD;
//	bool   tempB;

	//the headersize in bytes//
	long headersize = 0;
	*this >> headersize;
	std::cout << "Headersize: " << headersize << std::endl;
	//nbr of Channels in Instrument, check wether it changed is below//
	short nbrChannels = 0;
	*this >> nbrChannels;
	hi.NbrChannels = nbrChannels;
	std::cout << "nbrChannels: " << nbrChannels << std::endl;
	//the nbr of bytes (in file is the nbr of bits)
	*this >> tempS; tempS /= 8;
	//	changed = (tempS != fNbrBytes) || changed;
	hi.NbrBytes = tempS;
	std::cout << "Nbr of bytes in file: " << tempS << std::endl;
	*this >> tempD;
	//	changed = (TMath::Abs(tempD - fSampInter) > 1.e-4) || changed;
	hi.SampInter = tempD;
	//the nbr of samples in the window
	*this >> tempL;
	//	changed = (tempL != fNbrSamples) || changed;
	hi.NbrSamples = tempL;
	//the delaytime	
	*this >> tempD;
	//	changed = (TMath::Abs(tempD - fDelayTime) > 1.e-4) || changed;
	hi.DelayTime = tempD;
	//the trigger Channel
	*this >> tempS;
	//	changed = (tempS != fTrigChan) || changed;
	hi.TrigChan = tempS;
	std::cout << "Trigger channel: " << tempS << std::endl;
	//the trigger Level
	*this >> tempD;
	//	changed = (TMath::Abs(tempD - fTrigLevel) > 1.e-4) || changed;
	hi.TrigLevel = tempD;
	std::cout << "Trigger level: " << tempD << std::endl;
	//the trigger Slope
	*this >> tempS;
	//	changed = (tempS != fTrigSlope) || changed;
	hi.TrigSlope = tempS;
	//the used Channels Bitmask
	*this >> tempL;
	//	changed = (tempL != fUsedChans) || changed;
	hi.UsedChans = tempL;
	std::cout << "Used channels: " << tempL << std::endl;
	//the Channel Combination Bitmask
	*this >> tempL;
	//	changed = (tempL != fChanCombUsedChans) || changed;
	hi.ChanCombUsedChans = tempL;
	//the Nbr of Converters per Channel
	*this >> tempS;
	//	changed = (tempS != fNbrConPerCh) || changed;
	hi.NbrConPerCh = tempS;

	hi.ChI.clear();
	for (int i = 0; i < hi.NbrChannels; i++)
	{
		hi.ChI.push_back(ChannelInfo(i));
		ReadChannelInfo(i);
	}
}
void lma2root::ReadChannelInfo(int chNbr) {
	short  tempS;
	long   tempL;
	double tempD;
	//the fullscale in mV
	*this >> tempS;
	std::cout << "FullScale in mV: " << tempS << std::endl;
	//changed = (tempS != fFullscale) || changed;
	hi.ChI[chNbr].FullScale = tempS;
	//the offset in mV
	*this >> tempS;
	//changed = (tempS != fOffset) || changed;
	hi.ChI[chNbr].Offset = tempS;
	//the vertical Gain (conversion factor to convert the bits to mV
	*this >> tempD;
	//changed = (TMath::Abs(tempD - fGain) > 1.e-4) || changed;
	hi.ChI[chNbr].Gain = tempD;
	//the Baseline for zero substraction
	*this >> tempS;
	//changed = (tempS != fBaseline) || changed;
	hi.ChI[chNbr].Baseline = tempS;
	//the Noise Level for zero substraction
	*this >> tempS;
	//changed = (tempS != fNoise) || changed;
	hi.ChI[chNbr].Noise = tempS;
	//the stepsize of the zero substraction
	*this >> tempL;
	//changed = (tempL != fStsi) || changed;
	hi.ChI[chNbr].Stsi = tempL;
	//the BackSize of the zero substraction
	*this >> tempL;
	//changed = (tempL != fBs) || changed;
	hi.ChI[chNbr].Bs = tempL;
}
bool lma2root::newFile(const char* NewFileName) {

	//--first check wether there is a file already open--//
	if (file.is_open())
	{
		//if this Archive is writing we need to flush before we close the file
		if (isWriting)
		{
			file.flush();
			file.close();
		}
		//otherwise we just close the file here
		else
		{
			file.close();
		}
	}

	file.open(NewFileName, std::ios::in | std::ios::binary);
	//if the file has not been opened give an error message and return false
	if (!file.is_open())
	{
		std::cerr << "something went wrong opening the file \"" << NewFileName << std::endl;
		return false;
	}

	//--get the filesize--//
	file.seekg(0, std::ios::end);
	filesize = __int64(file.tellg());
	file.seekg(0, std::ios::beg);

	//if the file was opened fine return true
	return true;
}
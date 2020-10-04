#pragma once

#ifndef __ADC_helper_functions_h_
#define __ADC_helper_functions_h_

#ifdef FADC_SIGNAL_ANALYSIS_EXPORTS
#define FADC_SIGNAL_ANALYSIS_API __declspec(dllexport)
#else
#define FADC_SIGNAL_ANALYSIS_API __declspec(dllimport)
#endif

#include <complex>
#include <valarray>
#include <math.h>
#include <vector>
#include <fstream>


class bezier_class
{
public:
	bezier_class() {};
	~bezier_class() {};
	
private:
	double bez_xa0, bez_xa1, bez_xa2, bez_xa3, bez_xa4, bez_xa5;
	double bez_ya0, bez_ya1, bez_ya2, bez_ya3, bez_ya4, bez_ya5;
	double w0,w1,w2,w3,w4,w5,w6;

	double xA,xB,yA,yB;

public:
	bool last_command_was_successful;

	void bezier3_make_curve_parameters(double *x, double *y, double weight);
	void bezier4_make_curve_parameters(double *x, double *y, double weight);

	void calculate_bezier3_point(double t, double &x, double &y);
	void calculate_bezier4_point(double t, double &x, double &y);

	double get_x_at_y_bezier3(double y);
	double get_x_at_y_bezier4(double y);
};

class parabola_spline
{
public:
	parabola_spline(void) {};
	~parabola_spline(void) {};

	double parabola_fit(double y1, double y2, double y3, double x_pos);
	
	double get_parabola_vertex(double y1, double y2, double y3);
	double get_parabola_vertex(double y1, double y2, double y3, double y4, double dx);

	double get_y_at_x(int * array, __int32 number_of_points, double x);
	double get_y_at_x(double * array, __int32 number_of_points, double x);
};

class linear_gaus
{
	//-- Implementation of Nuclear Instruments and Methods 125 (1975) 289-291	--//
	//-- Fitting of Gaussian to Peaks by Non-Iterative methode					--//
	//-- T. Mukoyama															--//

public:
	int fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max,  double &x_0, double &sigma);
	int fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max, double &x_0, double &sigma, double &y_0);
	double gaus_fit(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max);

protected:
	void linear_weighted(double *x, double *y, double *w, __int32 nbrPoints, double &m, double &c);
	void get_y0(double *data, __int32 nbrPoints, __int32 data_min, __int32 data_max, double &m, double &c, double &x_0, double &sigma, double &y_0);
};

struct Entry {
	double	val;
	int		nbr;
	double	error;
};

class Trilinear
{
public:

	Trilinear();
	~Trilinear();

	void createMatrix(int x_bins, double x_min, double x_max, int y_bins, double y_min, double y_max, int z_bins, double z_min, double zmax);
	int fillMatrix(double x, double y, double z, double val);
	int fillErrorMatrix(double x, double y, double z, double error);

	void averageMatrix();
	void errorMatrix();
	void resetMatrix();
	void resetErrorMatrix();
	void resetNbrMatrix();

	double	GetVal(int x, int y, int z);
	void	SetVal(int x, int y, int z, double val);

	double interpolate3D(double x, double y, double z);

	//private:
	int x_bins;
	double x_min, x_max, x_delta, x_range;
	int y_bins;
	double y_min, y_max, y_delta, y_range;
	int z_bins;
	double z_min, z_max, z_delta, z_range;

	double matrix_volume;

	Entry ***m_matrix;
	int **yz_cnt_array;
};

class Bilinear
{
public:

	Bilinear();
	~Bilinear();

	void createMatrix(int x_bins, double x_min, double x_max, int y_bins, double y_min, double y_max);
	int fillMatrix(double x, double y, double val);
	//	int fillErrorMatrix(double x, double y, double error);

	void averageMatrix();
	//	void errorMatrix();
	//	void resetMatrix();
	//	void resetErrorMatrix();
	//	void resetNbrMatrix();

	double	GetVal(int x, int y);
	void	SetVal(int x, int y, double val);

	double interpolate(double x, double y);

	//private:
	int x_bins;
	double x_min, x_max, x_delta, x_range;
	int y_bins;
	double y_min, y_max, y_delta, y_range;

	double matrix_volume;
	Entry **m_matrix;
};

typedef std::valarray<std::complex<double>> CArray;

class FFT_class {
public:
	void fft(double *data_in, int nbrPoints, double *&data_fft_out, int nbrPointsFFT = 0);
	void ifft(double *data_fft_in, int nbrPoints, double *&data_out, double *&real_out, double *&imag_out, int nbrPointsFFT = 0);

	void fft_shift(double *data_fft, int nbrPoints, int nbrPointsFFT = 0);

	void fft_shift_undo(double *data_fft_shifted, int nbrPoints, int nbrPointsFFT = 0);

	void getAmplitude(double *data_fft_in, int nbrPoints, double *&data_Aftt_out, int nbrPointsFFT = 0);
	int getNbrPointsFFT(int nbrPoints);

	void fft_divide(double *fft_A, double *fft_B);
};

struct ChannelInfo
{
	ChannelInfo(int nbr) { ChNbr = nbr; };

	int ChNbr;
	short FullScale;
	short Offset;
	double Gain;
	short Baseline;
	short Noise;
	long Stsi;
	long Bs;
};

struct HeaderInfo
{
	short					 NbrChannels;
	short					 NbrBytes;								//! Nbr of bytes of the adc values (either 1 or 2)
	double					 SampInter;								//! the time between two consecutive points (in ns)
	long					 NbrSamples;							//! Nbr of Points (multiplied by the fSampInter it will give the timewindow in ns)
	double					 DelayTime;								//! the delay of the trigger with respect to the window
	short					 TrigChan;								//! the fTriggering Channel
	double					 TrigLevel;								//! the trigger Level from the Offset
	short					 TrigSlope;								//! which Slope was used by the fTrigger
	long					 UsedChans;								//! Bitmask discribing which Channels have been recorded 
	long					 ChanCombUsedChans;						//! Bitmask discribing which Converters per Channel have been used
	short					 NbrConPerCh;							//! tells how many converts per channel have been used
	std::vector<ChannelInfo> ChI;
};

class lma2root
{
public:
	lma2root(const char* FileName) {
		fileOpen = newFile(FileName);
		ReadLMAHeader();
	};
	~lma2root() { if (file.is_open()) file.close(); };

	lma2root & operator>>(char &c) { file.read(&c, sizeof(char)); return *this; }
	lma2root & operator>>(short &s) { file.read((char*)&s, sizeof(short)); return *this; }
	lma2root & operator>>(long &l) { file.read((char*)&l, sizeof(long)); return *this; }
	lma2root & operator>>(double &d) { file.read((char*)&d, sizeof(double)); return *this; }

	bool newFile(const char* NewFileName);

	void ReadLMAHeader();
	void readArray(void * data, long sizeOfData) { file.read((char*)data, sizeOfData); }
	void ReadChannelInfo(int chNbr);

	HeaderInfo hi;

	bool fileOpen;

private:
	std::fstream file;
	__int64 filesize;
	bool isWriting;
};

#endif
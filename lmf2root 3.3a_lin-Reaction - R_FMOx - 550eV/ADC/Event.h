#ifndef __Event_h_
#define __Event_h_

#include "../Src/OS_Version.h"
#include "MyPunkt.h"
#include "ADC_helper_functions.h"


#include <conio.h>			//TBB: alle drei zum tastendruck abfangen 	
#include <TSystem.h>
#include <TCanvas.h>

#include <iostream>

#include <TH1.h>
#include <TColor.h>
#include <TStyle.h>
#include <TGClient.h>
#include <TArrow.h>



#define MAX_NBR_PULSES		30
#define MAX_NBR_PEAKS		30
#define MAX_NBR_SAMPLES		500
#define	MAX_DATALENGHT		1000

//#include "Ueberstruct.h"

#pragma warning(disable : 4800)



struct Peak
{
	double	TimefromFit;					//Time extracted from fit routine
	double	Time;							//the time of the peaks, calculated from either cfd or com
	double	Cfd;							//the bin within the pulse array calculated from cfd (ss)
	double	Com;							//the bin within the pulse array calculated form com (ss)

	double	Time_ns;
	double	Time_ss;
	bool	is_TDC_signal;

	bool	Polarity;						//the polarity of the peak 0 neg, 1 pos
	double	Slope;							//the slope of the leading edge
	double	Slope_fe;						//the slope of the falling edge

	long	Maxpos;							//the position where the maximum of peak is
	double	Maximum;						//the height in bits
	double	Height;							//the height in mV
	double	HeightfromFit;					//the height when you use the substraction cfd

	double	Height_mV;
	double	Height_ss;

	double	Fwhm;							//the fwhm of the peak
	double	Width;							//the width at the bottom of the peak
	double	PosHalfLeft;					//the pos where the left edge crosses the half of the height
	double	PosHalfRight;					//the pos where the right edge crosses the half of the height

	double	Integral;						//the integral of the peak

	int	Startpos;						//the start postion of the peak
	int	Stoppos;						//the stop position of the peak

	double	scalefactor;					//Scalefactor obtained for substraction of meanpulse
	int		calcbinoffset;					//Caluclated bin offset
	double	runtime;						//if it is a anode signal, store the runtime for refined subtraction, initialized with -100.

	bool	is_valid;

	bool	dpa_is_potential_double;

	double	dpa_max1_pos;
	double	dpa_max1_hight;
	double	dpa_min_pos;
	double	dpa_min_hight;
	double	dpa_max2_pos;
	double	dpa_max2_hight;
	double	dpa_fwhm1_pos;
	double	dpa_fwhm1_hight;
	double	dpa_fwhm1_slope;
	double	dpa_fwhm2_pos;
	double	dpa_fwhm2_hight;
	double	dpa_fwhm2_slope;

	int		left_curved_cnt;				// points left curved
};

class Pulse
{
public:
	Pulse();
	~Pulse();

	Pulse*				Pulse_cpy();
	void				Pulse_init(Pulse *orig);

	// Pulse
	int					Pu_from_Ch;
	long				DataLength;						
	long long			timestamp;
	double				*Waveform;
	bool				is_TDC_signal;
	double				GetWaveformAt(double time_ss);		//FUNCTION

	// Peaks
	Peak				*Peaks;
	Peak				*Peaks_backup;
	int					NbrPeaks;
	int					NbrPeaks_backup;

	void				Peak_delete(int pk);			//FUNCTION			
	void				Peak_insert(int pk, Peak &newPeak);
	void				resetPeaks();					//FUNCTION
	
	void				Peak_backup();


	// MixedEvents
 	int					used_in_DPA;
	double				mixed_dTime_ss;
	double				mixed_Time1_ss;
	double				mixed_Time2_ss;
	int					mixed_Event;
	
	bool				is_Meanpulse;
	int					meanpulse_cnt;

	void				resetPulse();					//FUNCTION



private:
	double				parabola_fit(double y1, double y2, double y3, double x_pos);
	void				Peak_cpy(int pk_from, int pk_to);
	void				Peak_cpy(Peak *pk_from, Peak* pk_to);
};

struct DPAinfos {
	TH2D *hist_ADCpuls1_error;	// error puls 1 in mixed double pulse without double puls analysis
	TH2D *hist_ADCpuls2_error;	// error puls 2 in mixed double pulse without double puls analysis

	TH2D *hist_DPApuls1_error;	// error puls 1 in mixed double pulse with double puls analysis
	TH2D *hist_DPApuls2_error;	// error puls 2 in mixed double pulse with double puls analysis

	double ADC_risky;
};

struct Channel
{
	Pulse *Pulses;
	int NbrPulses;
};

class Event
{
public:
	Event(/*Ueberstruct *Ueber, ADCinfos *adc_info*/);
	~Event();

	long id;

	Ueberstruct *Ueber_ptr; 
	parabola_spline *parabola;
	bezier_class *bezier;

	Channel *channels[NUM_CHANNELS];
	
	struct ADCinfos *adc_set;

//	void showEvent();
//	void showPulse(int ch, int pu);
//	void showPulse(Pulse *curPulse);
//	void showPulseFFT(Pulse *curPulse, Pulse *meanPulse);

	TCanvas *EventCanvas;



	void prepair_pulse(Pulse *curPulse, int NbrPointsFFT);
	void shift_pulse(Pulse *curPulse);
	void compress_pulse(Pulse *curPulse, double compress_rate);

	void resetEvent();
	
private:
	//void startstop	(Pulse *pulse, Peak *newPeak, double threshold_mV = 0);
	//void maximum	(Pulse *pulse, Peak *newPeak);
	//void fwhm		(Pulse *pulse, Peak *newPeak);
	//void CoM		(Pulse *pulse, Peak *newPeak);

	//bool FindPeaksStartStop(Pulse *pulse);
	//bool FindPeaksMaxima(Pulse *pulse);
};



#endif
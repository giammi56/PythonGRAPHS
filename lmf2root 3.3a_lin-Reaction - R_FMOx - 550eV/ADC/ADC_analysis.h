#ifndef __ADC_analysis_h_
#define __ADC_analysis_h_

#define ADCVERSION 1.0

#include "../Src/OS_Version.h"
#include "../Src/Ueberstruct.h"
//#include "ADC_helper_functions.h"
#include "Event.h"
#include "ADC_settings.h"
#include "H1D.h"

#pragma warning(disable : 4800)

class Pulse;
class Histo;
class TCanvas;
class Trilinear;
class Bilinear;
class H1D;
class H2D;

class ADC_analysis
{
public:
	ADC_analysis(double threshold_mV, int delay_ns, int use_dpa, double fraction, int polarity, int method, double edge_level, double	gain, double samplerate, double fullscale);
	~ADC_analysis();

	Event		*m_event;

	int			initialize_DPA_analysis(Pulse *pu);	
	int			analyze_Pulse(Pulse *pu);

	ADCsettings	*m_adc_settings;

private:
	int			adc_for_channel;	// adc analysis for channel
	bool		dpa_initialized;	// dpa analysis is initialized

	// ***************************************************************************************
	// ADC Analysis
	void		FindPeaks(Pulse *pulse, int start = 0, int stop = 0);
	void		FindPeaksNoDPA(Pulse *pulse, int start = 0, int stop = 0);
	void		FindPeaksCOM(Pulse *pulse);
	void		FindPeaksLE(Pulse *pulse);
	void		FindPeaksDLE(Pulse *pulse);
	bool		FindPeaksStartStop(Pulse *pulse);
	bool		FindPeaksMaxima(Pulse *pulse);
	void		FindPeaksCFD(Pulse *pulse, long delay_ns_ext, double walk_ext, double threshold_mV_ext, double fraction_ext, int start, int stop);
	void		startstop(Pulse *pulse, Peak *newPeak, double threshold_mV = 0);
	void		maximum(Pulse *pulse, Peak *newPeak);
	void		fwhm(Pulse *pulse, Peak *newPeak);
	void		CoM(Pulse *pulse, Peak *newPeak);

	// ***************************************************************************************
	// Double Pulse Analysis
	void		FindDoublePulse(Pulse *PulseIn);
	void		FindMaximaDPA(Pulse *PulseIn);
	int			DPA_initialisation_status;	//-2: collect ADC infos, -1: accumulate meanpulse, 0: collect ADC and DPA error and evaluate, 1: is initialized

	// DPA - ADC infos
	int			DPA_init_nbrPulses_ADCinfos;
	int			DPA_init_curPulses_ADCinfos;
	void		collect_ADCinfos(Pulse *pu);		
	void		evaluate_ADCinfos();				

	// DPA - Meanpulse
	int			DPA_init_nbrPulses_Meanpulse;
	int			DPA_init_curPulses_Meanpulse;
	bool		GetMeanpulse(Pulse* MeanPulse, Pulse *pu, int pk);
	void		accumulate_Meanpulse(Pulse *pu);
	void		evaluate_Meanpulses();	
	Trilinear	*m_adv_meanpulse;
	Trilinear	*m_adv_meanpulse_falling;
	Pulse*		m_meanpulse;
	Pulse		*MeanPulse1, *MeanPulse2, *subtractedPulse1, *subtractedPulse2;

	// DPA - DPA infos
	int			DPA_init_nbrPulses_DPAerror;
	int			DPA_init_curPulses_DPAerror;
	void		collect_ADCandDPA_PulseError(Pulse *mixedPulse);
	void		GetMixedPulse(Pulse *pu1, Pulse *pu2, Pulse *mixedPulse, int dSS);
	void		evaluate_ADCandDPA_PulseError();
//	void		evaluateADCpulses_error();				
	Pulse		*mixedPulse, *lastPulse;

	// ***************************************************************************************
	// helper functions
	bool		substractPulse(Pulse* substractedPulse, Pulse* original, Pulse* substract);
	void		GetRangeHisto(H1D* hist, double *range_low, double *range_high, double percentage = 95, double savety_per = 3);
	bezier_class *bezier;
	parabola_spline *parabola;
};

#endif
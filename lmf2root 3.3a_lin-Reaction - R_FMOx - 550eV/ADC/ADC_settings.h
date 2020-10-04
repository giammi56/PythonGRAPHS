#pragma once

class H1D;
class H2D;

class ADCsettings {
public:
	ADCsettings();
	ADCsettings(double threshold_mV, int delay_ns, int use_dpa, double fraction, int polarity, int method, double edge_level, double	gain, double samplerate, double fullscale);

	double		get_gain() { return gain; }
	double		get_samplerate() { return samplerate; }
	double		get_fullscale() { return fullscale; }

	double		get_threshold() { return threshold; }
	double		get_delay() { return delay; }
	int			get_polarity() { return polarity; }


	bool		adc_initialized;		//if adc channl is used

	double		threshold;				//in mV
	int			delay;					//in ns
	int			use_dpa;				//use double pulse analysis
	double		fraction;				//Fraction Ratio
	int			polarity;				//Polarity
	int			method;					//0: Cfd, 1: Com, 2: LE, 3: FE, 4: Max
	double		edge_level;				//ratio or edge mV

	double		gain;
	double		samplerate;				// 
	double		fullscale;

	H1D			*height_hist;
	H1D			*fwhm_hist;
	H1D			*slope_hist;
	H1D			*slope_fe_hist;
	H1D			*integral_hist;
	double		height_min, height_max;
	double		fwhm_min, fwhm_max;
	double		slope_min, slope_max;
	double		slope_fe_min, slope_fe_max;
	double		integral_min, integral_max;

	H2D			*hist_ADCpuls1_error;	// error puls 1 in mixed double pulse without double puls analysis
	H2D			*hist_ADCpuls2_error;	// error puls 2 in mixed double pulse without double puls analysis
	H2D			*hist_DPApuls1_error;	// error puls 1 in mixed double pulse with double puls analysis
	H2D			*hist_DPApuls2_error;	// error puls 2 in mixed double pulse with double puls analysis
	H1D			*hist_ADCvsDPA_puls1;
	H1D			*hist_ADCvsDPA_puls2;

	H1D			*hist_integral_dpa;		// integral of double pulses
	double		integral_dpa_min, integral_dpa_max;

	H1D			*hist_fwhm_dpa;
	H1D			*hist_height_dpa;

	double		Cfd_risky;				//below peak distance cfd algorithm gets risky (ss)
	double		Cfd_risky_ss;
	double		DPA_is_favorable_below_ss;

										// (D)ouble (P)ulse (A)nalysis
	bool		dpa_initialized;

};
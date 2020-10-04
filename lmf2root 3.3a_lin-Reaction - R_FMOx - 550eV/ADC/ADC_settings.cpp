#include "ADC_settings.h"

ADCsettings::ADCsettings() {
	this->threshold = 0;
	this->delay		= 0;
	this->use_dpa	= 0;
	this->fraction	= 0;
	this->polarity	= 0;
	this->method	= 0;
	this->edge_level= 0;
	this->gain		= 0;
	this->samplerate= 0;
	this->fullscale = 0;

	this->height_hist = nullptr;
	this->fwhm_hist = nullptr;
	this->slope_hist = nullptr;
	this->slope_fe_hist = nullptr;
	this->integral_hist = nullptr;
	this->height_min = 0;
	this->height_max = 0;
	this->fwhm_min = 0;
	this->fwhm_max = 0;
	this->slope_min = 0;
	this->slope_max = 0;
	this->slope_fe_min = 0;
	this->slope_fe_max = 0;
	this->integral_min = 0;
	this->integral_max = 0;

	this->hist_ADCpuls1_error = nullptr;	// error puls 1 in mixed double pulse without double puls analysis
	this->hist_ADCpuls2_error = nullptr;	// error puls 2 in mixed double pulse without double puls analysis
	this->hist_DPApuls1_error = nullptr;	// error puls 1 in mixed double pulse with double puls analysis
	this->hist_DPApuls2_error = nullptr;	// error puls 2 in mixed double pulse with double puls analysis
	this->hist_ADCvsDPA_puls1 = nullptr;
	this->hist_ADCvsDPA_puls2 = nullptr;

	this->hist_integral_dpa = 0;			// integral of double pulses
	this->integral_dpa_min = 0;
	this->integral_dpa_max = 0;

	this->hist_fwhm_dpa = nullptr;
	this->hist_height_dpa = nullptr;

	this->Cfd_risky = 0;					//below peak distance cfd algorithm gets risky (ss)
	this->Cfd_risky_ss = 0;
	this->DPA_is_favorable_below_ss = 0;

	this->dpa_initialized = false;
}

ADCsettings::ADCsettings(double threshold_mV, int delay_ns, int use_dpa, double fraction, int polarity, int method, double edge_level, double	gain, double samplerate, double fullscale){
	this->threshold = threshold_mV;
	this->delay = delay_ns;
	this->use_dpa = use_dpa;
	this->fraction = fraction;
	this->polarity = polarity;
	this->method = method;
	this->edge_level = edge_level;
	this->gain = gain;
	this->samplerate = samplerate;
	this->fullscale = fullscale;


	this->height_hist = nullptr;
	this->fwhm_hist = nullptr;
	this->slope_hist = nullptr;
	this->slope_fe_hist = nullptr;
	this->integral_hist = nullptr;
	this->height_min = 0;
	this->height_max = 0;
	this->fwhm_min = 0;
	this->fwhm_max = 0;
	this->slope_min = 0;
	this->slope_max = 0;
	this->slope_fe_min = 0;
	this->slope_fe_max = 0;
	this->integral_min = 0;
	this->integral_max = 0;

	this->hist_ADCpuls1_error = nullptr;	// error puls 1 in mixed double pulse without double puls analysis
	this->hist_ADCpuls2_error = nullptr;	// error puls 2 in mixed double pulse without double puls analysis
	this->hist_DPApuls1_error = nullptr;	// error puls 1 in mixed double pulse with double puls analysis
	this->hist_DPApuls2_error = nullptr;	// error puls 2 in mixed double pulse with double puls analysis
	this->hist_ADCvsDPA_puls1 = nullptr;
	this->hist_ADCvsDPA_puls2 = nullptr;

	this->hist_integral_dpa = nullptr;			// integral of double pulses
	this->integral_dpa_min = 0;
	this->integral_dpa_max = 0;

	
	this->hist_fwhm_dpa = nullptr;
	this->hist_height_dpa = nullptr;

	this->Cfd_risky = 0;					//below peak distance cfd algorithm gets risky (ss)
	this->Cfd_risky_ss = 0;
	this->DPA_is_favorable_below_ss = 0;

	this->dpa_initialized = false;
}
#include "..\Src\OS_Version.h"

#include "../Src/Ueberstruct.h"
#include "ADC_analysis.h"
//#include "DPA_analysis.h"
#include "ADC_helper_functions.h"
//#include <afx.h>
#include <Windows.h>
#include <iostream>



using namespace std;
TCanvas* Canv1(const char *name, int hor, int ver)
{
	TCanvas * MyC = new TCanvas(name,"MyCanvas",-1,1,hor*500,ver*500);
	MyC->Divide(hor,ver,static_cast<float>(0.0001),static_cast<float>(0.0001));
	MyC->ToggleEventStatus();
	for (int i=1;i<=hor*ver;++i)
	{
		MyC->cd(i);
		gPad->SetBottomMargin(static_cast<float>(0.15));
		gPad->SetRightMargin(static_cast<float>(0.15));
		gPad->SetLeftMargin(static_cast<float>(0.15));
		gPad->SetTopMargin(static_cast<float>(0.15));
	}
	return MyC;
}
void linearRegression(const int nbrPoints, const double x[], const double y[], double &m, double &c);

ADC_analysis::ADC_analysis(double threshold_mV, int delay_ns, int use_dpa, double fraction, int polarity, int method, double edge_level, double	gain, double samplerate, double fullscale)
{
	m_adc_settings = new ADCsettings(threshold_mV, delay_ns, use_dpa, fraction, polarity, method, edge_level, gain, samplerate, fullscale);
	lastPulse					= nullptr;
	mixedPulse					= new Pulse();
	m_adv_meanpulse				= nullptr;
	m_meanpulse					= nullptr;
	m_adv_meanpulse_falling		= nullptr;
	DPA_initialisation_status	= -2;

	dpa_initialized				= false;

	DPA_init_nbrPulses_ADCinfos		= 20000;
	DPA_init_curPulses_ADCinfos		= 0;
	DPA_init_nbrPulses_Meanpulse	= 100000;
	DPA_init_curPulses_Meanpulse	= 0;
	DPA_init_nbrPulses_DPAerror		= 10000;
	DPA_init_curPulses_DPAerror		= 0;

	MeanPulse1					= nullptr;
	MeanPulse2					= nullptr;
	subtractedPulse1			= nullptr; 
	subtractedPulse2			= nullptr;

	if (m_adc_settings->use_dpa) {
		MeanPulse1 = new Pulse();
		MeanPulse2 = new Pulse();
		subtractedPulse1 = new Pulse();
		subtractedPulse2 = new Pulse();
	}

	bezier = new bezier_class();
	parabola = new parabola_spline();
}
ADC_analysis::~ADC_analysis()
{
	//	delete parabola;
	//	delete bezier;

	if (m_adc_settings)
		delete m_adc_settings;
	m_adc_settings = 0;

	if (lastPulse)
		delete lastPulse;
	lastPulse = 0;
	//	delete lma;
	//	if(!canv)
	//		delete canv;

	if (!MeanPulse1)
		delete MeanPulse1;
	MeanPulse1 = 0;
	if (!MeanPulse2)
		delete MeanPulse2;
	MeanPulse2 = 0;
	if (!subtractedPulse1)
		delete subtractedPulse1;
	subtractedPulse1 = 0;
	if (!subtractedPulse2)
		delete subtractedPulse2;
	subtractedPulse2 = 0;

	delete bezier;
	bezier = 0;
	delete parabola;
	parabola = 0;
}

int ADC_analysis::analyze_Pulse(Pulse *pu) {
	this->FindPeaks(pu);
	return 1;
}
int ADC_analysis::initialize_DPA_analysis(Pulse *pu) {
	// 1: keep feeding data, 0: initialisation finished 
	
	if (1 == this->DPA_initialisation_status)
		return 0;
	
	if (!this->m_adc_settings)
		int something_is_wrong = 0;

	this->FindPeaks(pu);

	adc_for_channel = pu->Pu_from_Ch;								// todo: müsste eigentlich nur 1mal init werden..

	// collect ADC infos
	if (-2 == this->DPA_initialisation_status) {	
		if (DPA_init_curPulses_ADCinfos < DPA_init_nbrPulses_ADCinfos) {
			collect_ADCinfos(pu);
			DPA_init_curPulses_ADCinfos++;
			return 1;
		}
		evaluate_ADCinfos();
		this->DPA_initialisation_status = -1;
		std::cout << "\n1/3 - ADC infos " << pu->Pu_from_Ch << " initialized";
	}

	// accumulate meanpulse
	if (-1 == this->DPA_initialisation_status) {
		if (DPA_init_curPulses_Meanpulse < DPA_init_nbrPulses_Meanpulse) {
			accumulate_Meanpulse(pu);
			DPA_init_curPulses_Meanpulse++;
			return 1;
		}
		evaluate_Meanpulses();
		this->DPA_initialisation_status = 0;
		std::cout << "\n2/3 - Meanpulse " << pu->Pu_from_Ch << " initialized";
	}

	// calc DPA error
	if (!lastPulse) {
		lastPulse = pu->Pulse_cpy();
		return 1;
	}

	if (0 == this->DPA_initialisation_status) {
		if (DPA_init_curPulses_DPAerror < DPA_init_nbrPulses_DPAerror) {
			for (int dSamplesteps = 0; dSamplesteps < m_adc_settings->fwhm_max * 2.5; dSamplesteps++) {	// 2.5*fwhm = pulses are not overlapping
				GetMixedPulse(lastPulse, pu, mixedPulse, dSamplesteps);
				collect_ADCandDPA_PulseError(mixedPulse);
			}

			DPA_init_curPulses_DPAerror++;
			return 1;
		}
		evaluate_ADCandDPA_PulseError();
		this->DPA_initialisation_status = 1;
		std::cout << "\n3/3 - DPA  initialized";
	}

	dpa_initialized = true;
	this->DPA_initialisation_status = 1;
	return 0;
}

// #########################################################PRIVATE###############################################################
void ADC_analysis::FindPeaks(Pulse *pulse, int start, int stop) {
	//--Methodes--//
	// 0: Constand fraction discrimination
	// 1: Center of mass
	// 2: Leading Edge 
	// 3: Dynamic Leading Edge
	// 4: Falling Edge
	// 5: 
	//if (m_adc_settings->use_dpa && dpa_initialized) {
	//	FindDoublePulse(pulse);
	//	return;
	//}

	if (pulse->is_TDC_signal) {
		if (!pulse->Peaks) {
			pulse->Peaks = new Peak[MAX_NBR_PEAKS];
			pulse->NbrPeaks = 0;
		}
		
		pulse->NbrPeaks = 1;
		pulse->Peaks[0].is_valid = false;
		pulse->Peaks[0].Time_ss = 0;
		pulse->Peaks[0].Time_ns = 0;
		pulse->Peaks[0].Time = 0;
		pulse->Peaks[0].is_TDC_signal = true;
		return;
	}

	if (m_adc_settings->method == 0) {
		FindPeaksCFD(pulse, 0, 0, 0, 0, start, stop);
		if(dpa_initialized)
			FindDoublePulse(pulse);
		return;
	}
	else if (m_adc_settings->method == 1) {
		FindPeaksCOM(pulse);
		if (dpa_initialized)
			FindDoublePulse(pulse);
		return;
	}
	else if (m_adc_settings->method == 2) {
		FindPeaksLE(pulse);
		if (dpa_initialized)
			FindDoublePulse(pulse);
		return;
	}
	else if (m_adc_settings->method == 3) {
		FindPeaksDLE(pulse);
		if (dpa_initialized)
			FindDoublePulse(pulse);
		return;
	}
	else if (m_adc_settings->method == 4) {
		FindPeaksMaxima(pulse);
		return;
	}
	std::cout << "\n\n<X> ERROR in ADC_analysis::FindPeaks: pulse (" << pulse->Pu_from_Ch << "), method (" << m_adc_settings->method << ") unknown\n\n";
	return;
}
void ADC_analysis::FindPeaksNoDPA(Pulse *pulse, int start, int stop) {
	if (m_adc_settings->method == 0) {
		FindPeaksCFD(pulse, 0, 0, 0, 0, start, stop);
		return;
	}
	else if (m_adc_settings->method == 1) {
		FindPeaksCOM(pulse);
		return;
	}
	else if (m_adc_settings->method == 2) {
		FindPeaksLE(pulse);
		return;
	}
	else if (m_adc_settings->method == 3) {
		FindPeaksDLE(pulse);
		return;
	}
	else if (m_adc_settings->method == 4) {
		FindPeaksMaxima(pulse);
		return;
	}
	std::cout << "\n\n<X> ERROR in ADC_analysis::FindPeaks: pulse (" << pulse->Pu_from_Ch << "), method (" << m_adc_settings->method << ") unknown\n\n";
	return;
}
void ADC_analysis::FindPeaksCOM(Pulse *pulse) {
	//--go through the puls--//
	for (int i = 0; i < pulse->DataLength; i++) {
		//--above threshhold--//
		if (pulse->Waveform[i] > m_adc_settings->threshold * m_adc_settings->gain) {
			if (!pulse->Peaks) {
				pulse->Peaks = new Peak[MAX_NBR_PEAKS];
				pulse->NbrPeaks = 0;
			}
			if (pulse->NbrPeaks < MAX_NBR_PEAKS) {
				pulse->Peaks[pulse->NbrPeaks].is_valid = false;
				pulse->Peaks[pulse->NbrPeaks].is_TDC_signal = false;
				pulse->Peaks[pulse->NbrPeaks].TimefromFit = -1.;
				pulse->Peaks[pulse->NbrPeaks].Startpos = i;
			}

			while (i < pulse->DataLength && pulse->Waveform[i] > m_adc_settings->threshold * m_adc_settings->gain) {
				pulse->Peaks[pulse->NbrPeaks].Stoppos = i;
				i++;
			}

			if (pulse->Peaks[pulse->NbrPeaks].Stoppos - pulse->Peaks[pulse->NbrPeaks].Startpos > 2) {
				maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				CoM(pulse, &pulse->Peaks[pulse->NbrPeaks]);

				pulse->Peaks[pulse->NbrPeaks].Time_ns = ((pulse->Peaks[pulse->NbrPeaks].Com * m_adc_settings->samplerate) + pulse->timestamp) *0.001;
				pulse->Peaks[pulse->NbrPeaks].Time_ss = pulse->Peaks[pulse->NbrPeaks].Com;
				pulse->NbrPeaks++;
				//--too many peaks detected.. something is wrong--//
				if (pulse->NbrPeaks >= MAX_NBR_PEAKS) {
					std::cout << "\n\n<x> ERROR in ADC_analysis::FindPeaksCOM: pulse (" << pulse->Pu_from_Ch << ") too many Peaks detected.. aborted!\n\n";
					return;
				}
			}
		}
	}
	return;
}
void ADC_analysis::FindPeaksLE(Pulse *pulse) {
	//--go through the puls--//
	for (int i = 2; i < pulse->DataLength - 1; i++) {
		if (pulse->Waveform[i] > m_adc_settings->edge_level * m_adc_settings->gain) {
			if (pulse->Waveform[i + 1] > pulse->Waveform[i]) {

				if (!pulse->Peaks) {
					pulse->Peaks = new Peak[MAX_NBR_PEAKS];
					pulse->NbrPeaks = 0;
				}
				if (pulse->NbrPeaks < MAX_NBR_PEAKS) {
					pulse->Peaks[pulse->NbrPeaks].Startpos = i;
				}

				double m = 0, c = 0;
				double x[4];
				double y[4];
				x[0] = i - 2;	y[0] = abs(pulse->Waveform[i - 2]);
				x[1] = i - 1;	y[1] = abs(pulse->Waveform[i - 1]);
				x[2] = i;		y[2] = abs(pulse->Waveform[i]);
				x[3] = i + 1;	y[3] = abs(pulse->Waveform[i + 1]);
				linearRegression(4, x, y, m, c);

				pulse->Peaks[pulse->NbrPeaks].Time_ss = (m_adc_settings->edge_level * m_adc_settings->gain - c) / m;
				// strange peak / pulse
				if (pulse->Peaks[pulse->NbrPeaks].Time_ss < 0 || pulse->Peaks[pulse->NbrPeaks].Time_ss > pulse->DataLength) {
					pulse->Peaks[pulse->NbrPeaks].Time_ss = 0;
					pulse->Peaks[pulse->NbrPeaks].is_valid = false;
					pulse->Peaks[pulse->NbrPeaks].is_TDC_signal = false;
					return;
				}


				i++;	//catch up
				while (i < pulse->DataLength - 1 && pulse->Waveform[i] > m_adc_settings->edge_level * m_adc_settings->gain &&
					pulse->Waveform[i + 1] > m_adc_settings->edge_level * m_adc_settings->gain) {
					i++;
				}
				pulse->Peaks[pulse->NbrPeaks].Stoppos = i;

				startstop(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);

				pulse->NbrPeaks++;
			}
		}
	}

}
void ADC_analysis::FindPeaksDLE(Pulse *pulse) {
	double hight = 0;
	int max_pos = 0;
	double threshhold = m_adc_settings->threshold * m_adc_settings->gain;

	//--go through pulse--//
	for (int i = 0; i < pulse->DataLength - 4; i++) {
		if (pulse->Waveform[i] > hight) {
			hight = pulse->Waveform[i];
			max_pos = i;
			if (pulse->Waveform[i + 1] <= hight && pulse->Waveform[i + 2] < hight && hight > m_adc_settings->threshold * m_adc_settings->gain) {
				if (!pulse->Peaks) {
					pulse->Peaks = new Peak[MAX_NBR_PEAKS];
					memset(pulse->Peaks, 0, MAX_NBR_PEAKS*sizeof(pulse->Peaks));
					pulse->NbrPeaks = 0;
				}
				if (pulse->NbrPeaks < MAX_NBR_PEAKS) {
					pulse->Peaks[pulse->NbrPeaks].Maxpos = long(double(max_pos) + 0.5 * (pulse->Waveform[max_pos + 1] - pulse->Waveform[max_pos - 1]) / (2 * pulse->Waveform[max_pos] - pulse->Waveform[max_pos - 1] - pulse->Waveform[max_pos + 1]));
					pulse->Peaks[pulse->NbrPeaks].Maximum = parabola->get_y_at_x(pulse->Waveform, pulse->DataLength, pulse->Peaks[pulse->NbrPeaks].Maxpos);

					//--maxima found.. go left until under edge limit--// 
					int j = max_pos;
					while (j > 0) {
						if (pulse->Waveform[j] < m_adc_settings->edge_level * pulse->Peaks[pulse->NbrPeaks].Maximum) {
							double x[3], y[3];
							x[0] = j - 1;  y[0] = pulse->Waveform[j - 1];
							x[1] = j;	   y[1] = pulse->Waveform[j];
							x[2] = j + 1;  y[2] = pulse->Waveform[j + 1];

							bezier->bezier3_make_curve_parameters(x, y, 1.);
							pulse->Peaks[pulse->NbrPeaks].Time_ss = bezier->get_x_at_y_bezier3(m_adc_settings->edge_level * pulse->Peaks[pulse->NbrPeaks].Maximum);

							if (pulse->Peaks[pulse->NbrPeaks].Time_ss < 0)
								int stop = 4;

							break;
						}
						j--;
					}

					startstop(pulse, &pulse->Peaks[pulse->NbrPeaks]);
					maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
					fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);
					CoM(pulse, &pulse->Peaks[pulse->NbrPeaks]);

					i = max_pos + 4;

					if (pulse->Peaks[pulse->NbrPeaks].Time_ss < 0) {
						pulse->Peaks[pulse->NbrPeaks].Time_ss = 0;
						pulse->Peaks[pulse->NbrPeaks].is_valid = false;
						pulse->Peaks[pulse->NbrPeaks].is_TDC_signal = false;
					}

					pulse->NbrPeaks++;
					continue;
				}
			}
		}
	}
}
bool ADC_analysis::FindPeaksStartStop(Pulse *pulse) {
	for (int i = 0; i < pulse->DataLength; i++) {
		if (pulse->Waveform[i] > m_adc_settings->threshold * m_adc_settings->gain) {
			if (!pulse->Peaks) {
				pulse->Peaks = new Peak[MAX_NBR_PEAKS];
				pulse->NbrPeaks = 0;
			}
			if (pulse->NbrPeaks < MAX_NBR_PEAKS) {
				pulse->Peaks[pulse->NbrPeaks].Startpos = i;
				pulse->NbrPeaks++;
			}
			while (i < pulse->DataLength && pulse->Waveform[i] > m_adc_settings->threshold * m_adc_settings->gain) {
				i++;
				pulse->Peaks[pulse->NbrPeaks].Stoppos = i;
				break;
			}
			break;
		}
	}
	if (pulse->Peaks[pulse->NbrPeaks].Stoppos > MAX_DATALENGHT) {
		pulse->Peaks[pulse->NbrPeaks].Stoppos = MAX_DATALENGHT;
		std::cout << "\n some errors..";
	}

	return true;
}
bool ADC_analysis::FindPeaksMaxima(Pulse *pulse) {
	//--go through the pulse--//
	int height = 0;
	int max_height = 0;
	int max_pos = 0;

	if (!pulse->Peaks) {
		pulse->Peaks = new Peak[MAX_NBR_PEAKS];
		pulse->NbrPeaks = 0;
	}

	for (int i = 2; i < pulse->DataLength - 1; i++) {

		if (pulse->Waveform[i] > this->m_adc_settings->threshold*this->m_adc_settings->gain) {	// above threshold

			if (pulse->Waveform[i] > max_height && pulse->Waveform[i-1] < pulse->Waveform[i]) {	// maxima candidate
				max_height = int(pulse->Waveform[i]);
				max_pos = i;
			}

			if (pulse->Waveform[i] < max_height * 0.95 && max_pos != 0) {	// check if point is 5 % lower than max + bigger threshold
				pulse->Peaks[pulse->NbrPeaks].Time_ss = double(max_pos) + 0.5 * (abs(pulse->Waveform[max_pos + 1]) - abs(pulse->Waveform[max_pos - 1])) / (2 * abs(pulse->Waveform[max_pos]) - abs(pulse->Waveform[max_pos - 1]) - abs(pulse->Waveform[max_pos + 1]));
				
				startstop(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				CoM(pulse, &pulse->Peaks[pulse->NbrPeaks]);

				pulse->NbrPeaks++;

				max_height = int(pulse->Waveform[i]);
				max_pos = 0;
			}
		}
	}



	//for (int i = PulseIn->Peaks[pk].Startpos; i < PulseIn->Peaks[pk].Stoppos; i++) {
	//	if (abs(Data[i]) > hight) {
	//		hight = abs(Data[i]);
	//		pos = i;
	//		if (abs(Data[i + 1]) <= hight && abs(Data[i + 2]) < hight && hight > this->m_adc_settings->threshold*this->m_adc_settings->gain) {
	//			PulseIn->Peaks[pk].dpa_max1_pos = double(pos) + 0.5 * (abs(Data[pos + 1]) - abs(Data[pos - 1])) / (2 * abs(Data[pos]) - abs(Data[pos - 1]) - abs(Data[pos + 1]));
	//			PulseIn->Peaks[pk].dpa_max1_hight = parabola->get_y_at_x(Data, pLength, PulseIn->Peaks[pk].dpa_max1_pos);
	//			break;
	//		}
	//	}
	//}



	return true;
}
void ADC_analysis::FindPeaksCFD(Pulse *pulse, long delay_ns_ext, double walk_ext, double threshold_mV_ext, double fraction_ext, int start, int stop) {
	//--get the right cfd settings--//
	long delay = (long)(m_adc_settings->delay * 1000 / m_adc_settings->samplerate);	//ns -> sampleinterval units
	if (delay_ns_ext > 0)
		delay = (long)(delay_ns_ext * 1000 / m_adc_settings->samplerate);

	double walk = 0; //m_adc_settings->walk			*	m_adc_settings->gain;										//mV -> ADC Bytes
	if (walk_ext > 0)
		walk = walk_ext;

	double threshold = m_adc_settings->threshold		*	m_adc_settings->gain;								//mV -> ADC Bytes
	if (threshold_mV_ext > 0)
		threshold = threshold_mV_ext *	m_adc_settings->gain;

	double fraction = m_adc_settings->fraction;
	if (fraction_ext > 0)
		fraction = fraction_ext;

	const double *Data = &pulse->Waveform[0];
	long pLength = pulse->DataLength;

	if (pLength > MAX_DATALENGHT || pLength < 3)	// 3 is arb.. maybe posible..
		return;

	//--only look for peaks in certain area--//
	if (start > 0 && start < stop && stop < pLength) {
		pLength = stop;
	}
	else {
		start = 0;
		stop = 0;
	}


	//--go through the puls--//
	for (int i = delay + 1 + start; i<pLength - 2; ++i)
	{
		const double fx = pulse->Waveform[i];			//the original Point at i
		const double fxd = pulse->Waveform[i - delay];	//the delayed Point	at i 
		const double fsx = -fx*fraction + fxd;			//the calculated CFPoint at i

		const double fx_1 = pulse->Waveform[i + 1]/* - static_cast<double>(baseline)*/;		//original Point at i+1
		const double fxd_1 = pulse->Waveform[i + 1 - delay]/* - static_cast<double>(baseline)*/;	//delayed Point at i+1
		const double fsx_1 = -fx_1*fraction + fxd_1;							//calculated CFPoint at i+1

																				//check wether the criteria for a Peak are fullfilled
		if (((fsx - walk) * (fsx_1 - walk)) <= 0) {	//one point above one below the walk
			if (TMath::Abs(fx) > threshold)				//original point above the threshold
			{
				//--it could be that the first criteria is 0 because	--//
				//--one of the Constant Fraction Signal Points or both	--//
				//--are exactly where the walk is						--//
				if (TMath::Abs(fsx - fsx_1) < 1e-8)	//both points are on the walk
				{
					//--go to next loop until at least one is over or under the walk--//
					continue;
				}
				else if ((fsx - walk) == 0)		//only first is on walk
				{
					//--Only the fist is on the walk, this is what we want--//
					//--so:do nothing--//
				}
				else if ((fsx_1 - walk) == 0)		//only second is on walk
				{
					//--we want that the first point will be on the walk,--//
					//--so in the next loop this point will be the first.--//
					continue;
				}
				//does the peak have the right polarity?//
				//if two pulses are close together then the cfsignal goes through the walk//
				//three times, where only two crossings are good. So we need to check for//
				//the one where it is not good//
				if (fsx > fsx_1)		//neg polarity
					if (pulse->Waveform[i] > 0.)	//but pos Puls .. skip
						continue;
				if (fsx < fsx_1)		//pos polarity
					if (pulse->Waveform[i] < 0.)	//but neg Puls .. skip
						continue;


				//--later we need two more points, create them here--//
				const double fx_m1 = pulse->Waveform[i - 1];		//the original Point at i-1
				const double fxd_m1 = pulse->Waveform[i - 1 - delay];	//the delayed Point	at i-1 
				const double fsx_m1 = -fx_m1*fraction + fxd_m1;		//the calculated CFPoint at i-1

				const double fx_2 = pulse->Waveform[i + 2];			//original Point at i+2
				const double fxd_2 = pulse->Waveform[i + 2 - delay];	//delayed Point at i+2
				const double fsx_2 = -fx_2*fraction + fxd_2;		//calculated CFPoint at i+2


																	// Achim / Stefan
																	// Wir wollen, dass fsx_m1 auf der gleichen Seite von Null ist wie fsx,
																	// und fsx_2 genauso über oder unter Null ist wie fsx_1.
																	// Also: (-)	   *	(-)   sei >  0
				if (((fsx_m1 - walk) *  (fsx - walk)) <= 0)	continue;	//Sonst: Schleifendurchlauf dieses Peaks abbrechen
																		// Und:  (+)	   *	(+)   sei >  0
				if (((fsx_1 - walk) *(fsx_2 - walk)) <= 0)	continue;	//Sonst: Schleifendurchlauf dieses Peaks abbrechen

																		//--find x_0a (Nulldurchgang_a) with a linear interpolation between the two points--//							
				const double m_a = fsx_1 - fsx;						//Steigung bestimmen	//(fsx_1 - fsx)/(i+1 - i);
				if (m_a == 0.) continue;								//per continue wird die loop-expr-Klausel ausgeführt, dieser Schleifendurchlauf wird abgebrochen und der nächste Schleifendurchlauf gestartet
				const double xLin_a = i + (walk - fsx) / m_a;			//Nullstelle berechnen	//PSF fx = (x - i)*m + cfs[i]
																		//--find x_0b (Nulldurchgang_b) with a linear interpolation between the two other points--//
				const double m_b = (fsx_2 - fsx_m1) / 3.;				//Steigung bestimmen	//(fsx_2 - fsx_m1)/(i+2 -(i-1));
				if (m_b == 0.) continue;
				const double xLin_b = i - 1 + (walk - fsx_m1) / m_b;	//Nullstelle berechnen	//PSF fx_m1 = (x - (i-1))*m + cfs[i]
																		//arithmetisches Mittel zwischen den beiden Nullstellen ergibt noch bessere/genauere Nullstelle (da hier die Unsicherheit, die auf jedem Messpunkt liegt, reduziert wird)
				const double pos = (xLin_a + xLin_b) / 2.;
				const double m = (m_a + m_b) / 2.;

				if (fsx < fsx_1 && fsx_m1 < fsx_2)	//check polarity first, only pulses with correct polarity are added as true pulses
				{
					if (pulse->NbrPeaks < MAX_NBR_PEAKS) {
						pulse->Peaks[pulse->NbrPeaks].is_valid = true;
						pulse->Peaks[pulse->NbrPeaks].is_TDC_signal = false;
						pulse->Peaks[pulse->NbrPeaks].TimefromFit = -1.;
						pulse->Peaks[pulse->NbrPeaks].Cfd = pos;				// TDC-Zeit des Pulses in [sampling-steps]
						pulse->Peaks[pulse->NbrPeaks].Time_ss = pos;
						//pulse->Peaks[pulse->NbrPeaks].Slope = m;
						pulse->Peaks[pulse->NbrPeaks].runtime = -100.;
						pulse->Peaks[pulse->NbrPeaks].Time = ((pulse->Peaks[pulse->NbrPeaks].Cfd * m_adc_settings->samplerate) + pulse->timestamp) *0.001;

						startstop(pulse, &pulse->Peaks[pulse->NbrPeaks]);
						maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
						fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);
						CoM(pulse, &pulse->Peaks[pulse->NbrPeaks]);

						pulse->NbrPeaks++;
					}
					else {
						std::cout << "\nToo many pulses detected.. something is wrong!";
					}
				}
			}
		}
	}

	//--if no peak was found try emergency cfd mode (more complex)--//
	//	the normal cfd mode has problems with fluxtuations..
	//    #
	//    # #    <-- normal cfd fails.. 
	//   #####
	//	#######
	if (pulse->NbrPeaks == 0) {
		int cfd_array_length = pLength - delay - 1;
		if (cfd_array_length < 1) return; // no peak found ACHIM TODO: check this here!
		double *cfd_array = new double[cfd_array_length];

		// calculate full cfd array
		int j = 0;
		for (int i = delay + 1 + start; i < pLength; ++i) {
			cfd_array[j] = -pulse->Waveform[i] * fraction + pulse->Waveform[i - delay];
			j++;
		}

		// find min and max
		int min_pos = 0, max_pos = 0;
		double min = cfd_array[0], max = cfd_array[0];
		for (int i = 1; i < cfd_array_length; i++) {
			if (min > cfd_array[i]) {
				min = cfd_array[i];
				min_pos = i;
			}
			if (max < cfd_array[i]) {
				max = cfd_array[i];
				max_pos = i;
			}
		}

		// find first point
		int left = 0, right = 0;
		for (int i = max_pos; i > min_pos; i--) {
			if (cfd_array[i] < walk) {
				right = i + 1;
				break;
			}
		}
		for (int i = min_pos; i < max_pos; i++) {
			if (cfd_array[i] > walk) {
				left = i - 1;
				break;
			}
		}

		if (left == right)
			right++;

		if ((right - left + 1) != 0) {
			const double m_a = (cfd_array[right + 1] - cfd_array[left]) / (right - left + 1);	//Steigung bestimmen	//(fsx_1 - fsx)/(i+1 - i);
			if (m_a == 0.) return;																//per continue wird die loop-expr-Klausel ausgeführt, dieser Schleifendurchlauf wird abgebrochen und der nächste Schleifendurchlauf gestartet
			const double xLin_a = (walk - cfd_array[left]) / m_a;								//Nullstelle berechnen	//PSF fx = (x - i)*m + cfs[i]
																								//--find x_0b (Nulldurchgang_b) with a linear interpolation between the two other points--//
			const double m_b = (cfd_array[right] - cfd_array[left - 1]) / (right - left + 1);	//Steigung bestimmen	//(fsx_2 - fsx_m1)/(i+2 -(i-1));
			if (m_b == 0.) return;
			const double xLin_b = (walk - cfd_array[left - 1]) / m_b - 1;						//Nullstelle berechnen	//PSF fx_m1 = (x - (i-1))*m + cfs[i]
																								//arithmetisches Mittel zwischen den beiden Nullstellen ergibt noch bessere/genauere Nullstelle (da hier die Unsicherheit, die auf jedem Messpunkt liegt, reduziert wird)
			const double pos = 1 + left + delay + start + (xLin_a + xLin_b) / 2.;

			if (pos > delay && pos < pLength - delay) {
				pulse->Peaks[pulse->NbrPeaks].is_valid = true;
				pulse->Peaks[pulse->NbrPeaks].is_TDC_signal = false;
				pulse->Peaks[pulse->NbrPeaks].TimefromFit = -1.;
				pulse->Peaks[pulse->NbrPeaks].Cfd = pos;				// TDC-Zeit des Pulses in [sampling-steps]
				pulse->Peaks[pulse->NbrPeaks].Time_ss = pos;
				//pulse->Peaks[pulse->NbrPeaks].Slope = m;
				pulse->Peaks[pulse->NbrPeaks].runtime = -100.;
				pulse->Peaks[pulse->NbrPeaks].Time = ((pulse->Peaks[pulse->NbrPeaks].Cfd * m_adc_settings->samplerate) + pulse->timestamp) *0.001;

				startstop(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				maximum(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				fwhm(pulse, &pulse->Peaks[pulse->NbrPeaks]);
				CoM(pulse, &pulse->Peaks[pulse->NbrPeaks]);

				pulse->NbrPeaks++;

				//this->showPulse(pulse);
			}
		}
		delete[] cfd_array;
	}
	return;
}
void ADC_analysis::startstop(Pulse *pulse, Peak *newPeak, double threshold_mV)
{
	//--this function will find the start and the stop of the peak--//
	const double *Data = &pulse->Waveform[0];
	const long pulslength = pulse->DataLength;

	long threshold = long(m_adc_settings->threshold * m_adc_settings->gain);
	if (threshold_mV != 0)
		threshold = long(threshold_mV * m_adc_settings->gain);

	const long center = long(newPeak->Time_ss);

	//go left from center until either i == 0, or the datapoint is inside the noise
	//or we go from the previous one (i+1) to the actual one (i) through the baseline
	int i = 0;
	for (i = center; i>0; i--)
		if ((abs(Data[i]) < threshold) || !((Data[i] >= 0) ^ (Data[i - 1] < 0)))
			break;
	int start = i;

	//go right form center until either i < pulslength, or the datapoint is inside the noise
	//or we go from the previous one (i-1) to the actual one (i) through the baseline
	for (i = center; i< pulslength; i++)
		if ((abs(Data[i]) < threshold) || !((Data[i] >= 0) ^ (Data[i - 1] < 0)))
			break;
	int stop = i;

	//found center is stuck below threshold, try to escape by touching forward in both directions
	if (start == stop) {
		int l = start, r = stop, t = 1;
		while (l > 0 && r < pulslength) {
			l--;
			r++;
			if (l > 0 && Data[l] > threshold) {
				for (i = l; i>0; i--)
					if ((abs(Data[i]) < threshold) || !((Data[i] >= 0) ^ (Data[i - 1] < 0)))
						break;
				start = i;
				break;
			}
			else if (r < pulslength && Data[r] > threshold) {
				for (i = r; i< pulslength; i++)
					if ((abs(Data[i]) < threshold) || !((Data[i] >= 0) ^ (Data[i - 1] < 0)))
						break;
				stop = i;
				break;
			}
		}
	}

	newPeak->Startpos = start;
	newPeak->Stoppos = stop;
	if (stop >= MAX_DATALENGHT) {
		newPeak->Stoppos = MAX_DATALENGHT;
		std::cout << "\n some errors..";
	}
}
void ADC_analysis::maximum(Pulse *pulse, Peak *newPeak)
{
	//--this function will find the maximum of the peak and its position--//
	const double Gain = m_adc_settings->gain;

	const double *Data = &pulse->Waveform[0];
	const long pLength = pulse->DataLength;
	const int start = newPeak->Startpos;
	const int stop = newPeak->Stoppos;

	double maximum_ = 0;
	int maxpos = 0;

	for (int i = start; i <= stop; ++i)
	{
		if (abs(Data[i]) > maximum_)
		{
			maximum_ = abs(Data[i]/*-baseline*/);
			maxpos = i;
		}
	}

	if (maxpos > 1 && maxpos + 1 != pulse->DataLength) {
		newPeak->Maxpos = long(double(maxpos) + 0.5 * (Data[maxpos + 1] - Data[maxpos - 1]) / (2 * Data[maxpos] - Data[maxpos - 1] - Data[maxpos + 1]));

		if (newPeak->Maxpos < 0) {
			//something is strange with the pulse..
			newPeak->is_valid = false;
			newPeak->Maxpos = 0;
		}

		newPeak->Maximum = parabola->get_y_at_x(pulse->Waveform, pLength, newPeak->Maxpos);
		newPeak->Height = (double)maximum_ / m_adc_settings->gain;
		newPeak->Height_mV = (double)maximum_ / m_adc_settings->gain;
		newPeak->Height_ss = (double)maximum_;

		if (newPeak->Maxpos < 0)
			int stop = 5;

		return;
	}
	else {
		newPeak->Maxpos = maxpos;
		newPeak->Maximum = maximum_;
		newPeak->Height = (double)maximum_ / m_adc_settings->gain;

		if (newPeak->Maxpos < 0) {
			newPeak->is_valid = false;
			newPeak->Maxpos = 0;
		}
		return;
	}
}
void ADC_analysis::fwhm(Pulse *pulse, Peak *newPeak)
{
	const double *Data = &pulse->Waveform[0];
	const long pLength = pulse->DataLength;

	//--get peak fwhm--//
	long fwhm_l = 0;
	long fwhm_r = 0;
	const double HalfMax = 0.5* newPeak->Maximum;

	////--go from middle to left until 0.5*height find first point that is above 0.5 height--//
	for (int i = newPeak->Maxpos; i >= 0; --i)
	{
		if (abs(Data[i]) < HalfMax)
		{
			fwhm_l = i + 1;
			break;
		}
	}

	//--go from middle to right until 0.5*height (find last point that is still above 0.5 Height--//
	for (int i = newPeak->Maxpos; i<pLength; ++i)
	{
		if (abs(Data[i]) < HalfMax)
		{
			fwhm_r = i - 1;
			break;
		}
	}

	//--if we found a right side and a left side, then--//
	//--compute the fwhm with a linear interpolation--//
	//--between the points that are left and right from--//
	//--where the fwhm is, else return here--//
	if (!fwhm_r || !fwhm_l)
		return;

	double lx[4];
	double ly[4];
	lx[0] = fwhm_l - 2;	ly[0] = abs(Data[fwhm_l - 2]);
	lx[1] = fwhm_l - 1;	ly[1] = abs(Data[fwhm_l - 1]);
	lx[2] = fwhm_l - 0;	ly[2] = abs(Data[fwhm_l - 0]);
	lx[3] = fwhm_l + 1;	ly[3] = abs(Data[fwhm_l + 1]);

	double rx[4];
	double ry[4];
	rx[0] = fwhm_r - 1;	ry[0] = abs(Data[fwhm_r - 1]);
	rx[1] = fwhm_r - 0;	ry[1] = abs(Data[fwhm_r - 0]);
	rx[2] = fwhm_r + 1;	ry[2] = abs(Data[fwhm_r + 1]);
	rx[3] = fwhm_r + 2;	ry[3] = abs(Data[fwhm_r + 2]);

	double mLeft, cLeft, mRight, cRight;
	linearRegression(4, lx, ly, mLeft, cLeft);
	linearRegression(4, rx, ry, mRight, cRight);

	//y = m*x+c => x = (y-c)/m;
	const double fwhm_L = (HalfMax - cLeft) / mLeft;
	const double fwhm_R = (HalfMax - cRight) / mRight;

	const double fwhm = fwhm_R - fwhm_L;
	//--set all found parameters--//
	newPeak->Slope = mLeft;
	newPeak->Slope_fe = mRight;
	newPeak->Fwhm = fwhm;
	newPeak->Width = newPeak->Stoppos - newPeak->Startpos;
	newPeak->PosHalfLeft = fwhm_L;
	newPeak->PosHalfRight = fwhm_R;
}
void ADC_analysis::CoM(Pulse *pulse, Peak *newPeak)
{
	const double threshold = m_adc_settings->threshold * m_adc_settings->gain;	//mV -> ADC Bytes


	const double *Data = &pulse->Waveform[0];
	const long pLength = pulse->DataLength;

	//--this function goes through the puls from start to stop and finds the center of mass--//
	double integral = 0;
	double wichtung = 0;
	const int start = newPeak->Startpos;
	const int stop = newPeak->Stoppos;

	for (int i = start; i <= stop; ++i)
	{
		integral += (abs(-Data[i]) - threshold);			//calc integral
		wichtung += ((abs(-Data[i]) - threshold)*i);		//calc weight
	}
	newPeak->Integral = integral;
	newPeak->Com = wichtung / integral;
}

// ***************************************************************************************
void ADC_analysis::FindDoublePulse(Pulse *PulseIn) {
	static int dpa_id = 0;								// dpa id for debugging (finding specific pulses)


	int NbrPeaks_before_DPA = PulseIn->NbrPeaks;
	double TimePeaks_before_DPA[10];
	if (NbrPeaks_before_DPA < 10) {
		for (int i = 0; i < NbrPeaks_before_DPA; i++) {
			TimePeaks_before_DPA[i] = PulseIn->Peaks[i].Time_ss;
		}
	}
	



	// check if pulse is potential DPA candidate
	FindMaximaDPA(PulseIn);

	//if (newEvent->channels[ch]->Pulses[pu].NbrPeaks == 1 && newEvent->channels[ch]->Pulses[pu].Peaks[0].dpa_is_potential_double && DPA_is_initialized) {
	//	show = false;
	//}




	//double *TimePeaks_before_DPA = new double[NbrPeaks_before_DPA];
	//for (int i = 0; i < NbrPeaks_before_DPA; i++) {
	//	TimePeaks_before_DPA[i] = PulseIn->Peaks[i].Time_ss;
	//}


	for (int pk = 0; pk < PulseIn->NbrPeaks; pk++) {
		if (!PulseIn->Peaks[pk].dpa_is_potential_double)
			continue;

		if (!this->GetMeanpulse(MeanPulse1, PulseIn, pk))
			continue;

		this->substractPulse(subtractedPulse1, PulseIn, MeanPulse1);
		subtractedPulse1->resetPeaks();
		this->FindPeaksNoDPA(subtractedPulse1, PulseIn->Peaks[pk].Startpos, PulseIn->Peaks[pk].Stoppos);

		if (subtractedPulse1->NbrPeaks >= 1) {
			subtractedPulse1->Peaks[0].dpa_fwhm1_slope = subtractedPulse1->Peaks[0].Slope;
			subtractedPulse1->Peaks[0].dpa_max1_hight = subtractedPulse1->Peaks[0].Height_ss;	// TO DO
			subtractedPulse1->Peaks[0].dpa_max1_pos = subtractedPulse1->Peaks[0].Maxpos;		// 1. wenn peaks zu nahe sind.. einen löschen und los


			if (!this->GetMeanpulse(MeanPulse2, subtractedPulse1, 0))
				continue;
			// 2. checken ob peaks richtig untersucht werden..
			this->substractPulse(subtractedPulse2, PulseIn, MeanPulse2);

			// find peak only in range of peak..					
			subtractedPulse2->resetPeaks();
			this->FindPeaksNoDPA(subtractedPulse2, PulseIn->Peaks[pk].Startpos, PulseIn->Peaks[pk].Stoppos);

			int something_is_wrong = false;

			if (subtractedPulse2->NbrPeaks >= 1) {
				// Two Peaks found..
				if (subtractedPulse1->NbrPeaks <= 2 && subtractedPulse2->NbrPeaks <= 3) {
					if (!PulseIn->Peaks_backup)
						PulseIn->Peak_backup();

					// delete old peak
					PulseIn->Peak_delete(pk);

					// insert new peaks
					if (subtractedPulse1->NbrPeaks == 1)
						PulseIn->Peak_insert(pk, subtractedPulse2->Peaks[0]);
					else if (subtractedPulse1->NbrPeaks == 2)
						PulseIn->Peak_insert(pk, subtractedPulse2->Peaks[1]);

					PulseIn->Peak_insert(pk + 1, subtractedPulse1->Peaks[0]);
					PulseIn->used_in_DPA = true;
					// one peak added.. 
					pk += 2;
				}
				else {
					std::cout << "\nError substractedPulse 2";
					//show = true;
				}
			}
		}
		else {
			std::cout << "\nError substractedPulse 2";
			//show = true;
		}


		if (true && dpa_initialized && subtractedPulse1 && subtractedPulse2) {
			if(NbrPeaks_before_DPA < PulseIn->NbrPeaks && PulseIn->NbrPeaks == 3){
			//if ((subtractedPulse1->NbrPeaks >= 1 || subtractedPulse2->NbrPeaks >= 1) && NbrPeaks_before_DPA < PulseIn->NbrPeaks && dpa_initialized) {
				// >>> debug - debug - debug - debug - debug
				if (subtractedPulse1->NbrPeaks >= 1 && subtractedPulse2->NbrPeaks >= 1) {
					TCanvas	*MixedCanvas = Canv1("EventCanvas", 2, 2);

					//std::cout << "\n\n-----------------------------------------------\nADC Debug Infos: ";
					//std::cout << "\nChannel: " << ch;
					//std::cout << "\nPulse:   " << pu;
					//std::cout << "\nNbrPeaks:" << newEvent->channels[ch]->Pulses[pu].NbrPeaks;
					//std::cout << "\nRUNNING COUNTER: " << running_counter < "\n\n";

					TH1D *Hist1, *Hist2, *Hist3;
					Hist1 = new TH1D("tmp1", "black: original, red: meanpulse, blue: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					Hist2 = new TH1D("tmp2", "black: original, red: meanpulse, blue: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					Hist3 = new TH1D("tmp3", "black: original, red: meanpulse, blue: result", PulseIn->DataLength, 1, PulseIn->DataLength);

					//TH2D *Hist2D_1;
					//char title[1000];
					//sprintf(title, "dE1: %f, dE2: %f, dT: %f", newEvent->channels[ch]->Pulses[pu].mixed_Time1_ss - subtractedPulse1->Peaks[0].Time_ss, newEvent->channels[ch]->Pulses[pu].mixed_Time2_ss - subtractedPulse1->Peaks[0].Time_ss, newEvent->channels[ch]->Pulses[pu].mixed_dTime_ss);
					//Hist2D_1 = new TH2D("tmp", title, 100, -10, 10, 40, 0, 100);

					for (int i = 0; i < PulseIn->DataLength; i++) {
						Hist1->SetBinContent(i, PulseIn->Waveform[i]);
						if (MeanPulse1)
							Hist2->SetBinContent(i, MeanPulse1->Waveform[i]);
						if (subtractedPulse1)
							Hist3->SetBinContent(i, subtractedPulse1->Waveform[i]);
					}

					//Hist2D_1->SetBinContent(1, 1, 2);
					//Hist2D_1->Fill(newEvent->channels[ch]->Pulses[pu].mixed_Time1_ss - subtractedPulse1->Peaks[0].Time_ss, newEvent->channels[ch]->Pulses[pu].mixed_dTime_ss, 1);
					//Hist2D_1->Fill(newEvent->channels[ch]->Pulses[pu].mixed_Time2_ss - subtractedPulse2->Peaks[0].Time_ss, newEvent->channels[ch]->Pulses[pu].mixed_dTime_ss, 2);

					// PAD 1
					MixedCanvas->cd(1);		//change to the right pad
					Hist1->Draw("");
					Hist2->SetLineColor(4);
					Hist2->Draw("same");
					Hist3->SetLineColor(2);
					Hist3->Draw("same");

					TArrow *arrow_1 = nullptr;
					if (subtractedPulse1) {
						arrow_1 = new TArrow[subtractedPulse1->NbrPeaks];
						for (int i = 0; i < subtractedPulse1->NbrPeaks; i++) {
							arrow_1[i].SetX1(subtractedPulse1->Peaks[i].Time_ss);
							arrow_1[i].SetY1(0);
							arrow_1[i].SetX2(subtractedPulse1->Peaks[i].Time_ss);
							arrow_1[i].SetY2(10000000000);
							arrow_1[i].SetOption("<|");
							arrow_1[i].SetArrowSize(float(0.01));
							arrow_1[i].Draw();
						}
					}

					// ----
					//std::cout << "\nPeaks Found without DPA (" << NbrPeaks_before_DPA << "):";
					TArrow *cfd_arrows = nullptr;
					cfd_arrows = new TArrow[PulseIn->NbrPeaks_backup];
					for (int i = 0; i < PulseIn->NbrPeaks_backup; i++) {
						//std::cout << "\nPeak[" << i << "]: " << PulseIn->Peaks_backup[i].Time_ss;

						cfd_arrows[i].SetX1(PulseIn->Peaks_backup[i].Time_ss);
						cfd_arrows[i].SetY1(0);
						cfd_arrows[i].SetX2(PulseIn->Peaks_backup[i].Time_ss);
						cfd_arrows[i].SetY2(-10000000000);
						cfd_arrows[i].SetOption("<|");
						cfd_arrows[i].SetArrowSize(float(0.01));
						cfd_arrows[i].SetLineColor(2);
						cfd_arrows[i].Draw();

						if (PulseIn->Peaks[i].Time_ss < 0)
							int stoping = 0;
					}

					std::cout << "\n\n******************************\n";
					std::cout << "Peaks found: (Channel: " << PulseIn->Pu_from_Ch << ", ID: " << dpa_id << ")";
					std::cout << "\nADC (" << NbrPeaks_before_DPA << ")";
					for (int i = 0; i < NbrPeaks_before_DPA; i++) {
						std::cout << "\t" << TimePeaks_before_DPA[i];
					}
					std::cout << "\nDPA (" << PulseIn->NbrPeaks << ")";
					for (int i = 0; i < PulseIn->NbrPeaks; i++) {
						std::cout << "\t" << PulseIn->Peaks[i].Time_ss;
					}
					if (!dpa_initialized) {
						std::cout << "\nMIX (2)\t" << PulseIn->mixed_Time1_ss << "\t" << PulseIn->mixed_Time2_ss;
						if(2 == PulseIn->NbrPeaks)
							std::cout << "\nErr (2)\t" << PulseIn->Peaks[0].Time_ss - PulseIn->mixed_Time1_ss << "\t" << PulseIn->Peaks[1].Time_ss - PulseIn->mixed_Time2_ss;
						if (1 == PulseIn->NbrPeaks)
							std::cout << "\nErr (1)\t" << PulseIn->Peaks[0].Time_ss - PulseIn->mixed_Time1_ss;

					}
					std::cout << "\n******************************\n";

					TArrow *final_arrows = nullptr;
					final_arrows = new TArrow[PulseIn->NbrPeaks];
					for (int i = 0; i < PulseIn->NbrPeaks; i++) {		
						final_arrows[i].SetX1(PulseIn->Peaks[i].Time_ss);
						final_arrows[i].SetY1(0);
						final_arrows[i].SetX2(PulseIn->Peaks[i].Time_ss);
						final_arrows[i].SetY2(10000000000);
						final_arrows[i].SetOption("<|");
						final_arrows[i].SetArrowSize(float(0.01));
						final_arrows[i].SetLineColor(3);
						final_arrows[i].Draw();

						if (PulseIn->Peaks[i].Time_ss < 0)
							int stoping = 0;
					}

					if (PulseIn->Peaks[pk].dpa_max1_pos > 0) {
						TArrow *arrow_max1 = new TArrow(PulseIn->Peaks[pk].dpa_max1_pos, 0., PulseIn->Peaks[pk].dpa_max1_pos, 1000000000, float(0.01), "<|");
						arrow_max1->Draw();
					}

					// PAD 2
					MixedCanvas->cd(2);
					TH1D *hist1_2, *hist2_2, *hist3_2;
					hist1_2 = new TH1D("tmp4", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist2_2 = new TH1D("tmp5", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist3_2 = new TH1D("tmp6", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					for (int i = 0; i < PulseIn->DataLength; i++) {
						hist1_2->SetBinContent(i, PulseIn->Waveform[i]);
						if (MeanPulse2)
							hist2_2->SetBinContent(i, MeanPulse2->Waveform[i]);
						if (subtractedPulse2)
							hist3_2->SetBinContent(i, subtractedPulse2->Waveform[i]);
					}
					hist1_2->Draw("");
					hist2_2->SetLineColor(4);	// 1: black, 2:red, 3: green, 4: blue, 6: magenta
					hist2_2->Draw("same");
					hist3_2->SetLineColor(2);
					hist3_2->Draw("same");

					TArrow *arrow_2 = nullptr;
					if (subtractedPulse2) {
						arrow_2 = new TArrow[subtractedPulse2->NbrPeaks];
						for (int i = 0; i < subtractedPulse2->NbrPeaks; i++) {
							arrow_2[i].SetX1(subtractedPulse2->Peaks[i].Time_ss);
							arrow_2[i].SetY1(0);
							arrow_2[i].SetX2(subtractedPulse2->Peaks[i].Time_ss);
							arrow_2[i].SetY2(10000000000);
							arrow_2[i].SetOption("<|");
							arrow_2[i].SetArrowSize(float(0.01));
							arrow_2[i].Draw();
						}
					}

					//pad 2
					//MixedCanvas->cd(3);
					//Hist2D_1->Draw("colz");


					// PAD 3
					MixedCanvas->cd(3);
					TH1D *hist_orig, *hist_sum, *hist_diff;
					hist_orig	= new TH1D("tmp7", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist_sum	= new TH1D("tmp8", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist_diff	= new TH1D("tmp9", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					for (int i = 0; i < PulseIn->DataLength; i++) {
						hist_orig->SetBinContent(i, PulseIn->Waveform[i]);
						if (subtractedPulse2) {
							hist_sum->SetBinContent(i, subtractedPulse1->Waveform[i] + subtractedPulse2->Waveform[i]);
							hist_diff->SetBinContent(i, PulseIn->Waveform[i] - subtractedPulse1->Waveform[i] - subtractedPulse2->Waveform[i]);
						}
					}
					hist_orig->Draw("");
					hist_sum->SetLineColor(4);	// 1: black, 2:red, 3: green, 4: blue, 6: magenta
					hist_sum->Draw("same");
					hist_diff->SetLineColor(2);
					hist_diff->Draw("same");


					// PAD 4
					MixedCanvas->cd(4);
					TH1D *hist_original, *hist_mean1, *hist_mean2;
					hist_original = new TH1D("tmp10", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist_mean1 = new TH1D("tmp11", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					hist_mean2 = new TH1D("tmp12", "black: original, blue: meanpulse, red: result", PulseIn->DataLength, 1, PulseIn->DataLength);
					for (int i = 0; i < PulseIn->DataLength; i++) {
						hist_original->SetBinContent(i, PulseIn->Waveform[i]);
						hist_mean1->SetBinContent(i, MeanPulse1->Waveform[i]);
						hist_mean2->SetBinContent(i, MeanPulse2->Waveform[i]);
					}
					hist_original->Draw("");
					hist_mean1->SetLineColor(4);	// 1: black, 2:red, 3: green, 4: blue, 6: magenta
					hist_mean1->Draw("same");
					hist_mean1->SetLineWidth(2);
					hist_mean2->SetLineColor(2);
					hist_mean2->Draw("same");
					hist_mean2->SetLineWidth(4);




					gPad->Update();

					while (!_kbhit())
					{
						gSystem->Sleep(50);
						gSystem->ProcessEvents();
					}

					char ch_;
					while (_kbhit()) ch_ = _getch();

					if (arrow_1)
						delete[] arrow_1;
					if (arrow_2)
						delete[] arrow_2;
					if (cfd_arrows)
						delete[] cfd_arrows;

					delete Hist1;			Hist1=0;
					delete Hist2;			Hist2=0;
					delete Hist3;			Hist3=0;

					delete hist1_2;			hist1_2=0;
					delete hist2_2;			hist2_2=0;
					delete hist3_2;			hist3_2=0;

					delete hist_orig;		hist_orig=0;
					delete hist_sum;		hist_sum=0;
					delete hist_diff;		hist_diff=0;

					delete hist_original;	hist_original=0;
					delete hist_mean1;		hist_mean1=0;
					delete hist_mean2;		hist_mean2=0;

					delete MixedCanvas;		MixedCanvas=0;
				}
			}
			// debug - debug - debug - debug - debug <<<
		}
	}

	dpa_id++;
}
void ADC_analysis::FindMaximaDPA(Pulse *PulseIn) {
	// if peaks are to close delete second peak and extend range of first peak
	if (dpa_initialized) {
		for (int pk = 0; pk < PulseIn->NbrPeaks - 1; pk++) {
			if (PulseIn->Peaks[pk + 1].Time_ss - PulseIn->Peaks[pk].Time_ss <= m_adc_settings->DPA_is_favorable_below_ss) {
				PulseIn->Peak_backup();
				PulseIn->Peaks[pk].Stoppos = PulseIn->Peaks[pk + 1].Stoppos;
				PulseIn->Peaks[pk].Integral += PulseIn->Peaks[pk + 1].Integral;
				PulseIn->Peak_delete(pk + 1);
			}
		}

		//--check if stop and start of neighboring peaks overlap
		for (int pk = 0; pk < PulseIn->NbrPeaks - 1; pk++) {
			if (PulseIn->Peaks[pk].Stoppos >= PulseIn->Peaks[pk + 1].Startpos) {
				if (!PulseIn->Peaks_backup)
					PulseIn->Peak_backup();
				PulseIn->Peaks[pk].Stoppos = PulseIn->Peaks[pk + 1].Stoppos;
				PulseIn->Peaks[pk].Integral += PulseIn->Peaks[pk + 1].Integral;
				PulseIn->Peak_delete(pk + 1);
			}
		}
	}
	else { //!DPA_is_initialized
		PulseIn->Peak_backup();
		if (PulseIn->NbrPeaks >= 2) {
			PulseIn->Peaks[0].Stoppos = PulseIn->DataLength;
			PulseIn->Peaks[0].Integral += PulseIn->Peaks[1].Integral;
			PulseIn->Peak_delete(1);
		}
	}

	double *Data;
	long pLength = 0;
	double hight = 0;
	int pos = 0;

	for (int pk = 0; pk < PulseIn->NbrPeaks; pk++) {
		// check if integral fulfils integral criteria for double pulse
		if (PulseIn->Peaks[pk].Integral <= this->m_adc_settings->integral_dpa_min)
			continue;

		if (PulseIn->Peaks[pk].Stoppos >= MAX_DATALENGHT) {
			//PulseIn->Peaks[pk].Stoppos = MAX_DATALENGHT;
			break;
		}

		PulseIn->Peaks[pk].dpa_max1_pos = 0;
		PulseIn->Peaks[pk].dpa_max1_hight = 0;
		PulseIn->Peaks[pk].dpa_max2_pos = 0;
		PulseIn->Peaks[pk].dpa_max2_hight = 0;
		PulseIn->Peaks[pk].dpa_is_potential_double = false;
		PulseIn->Peaks[pk].dpa_min_pos = 0;
		PulseIn->Peaks[pk].dpa_min_hight = 0;

		Data = &PulseIn->Waveform[0];
		pLength = PulseIn->DataLength;
		hight = 0;
		pos = 0;

		// find first maxima
		for (int i = PulseIn->Peaks[pk].Startpos; i < PulseIn->Peaks[pk].Stoppos; i++) {
			if (abs(Data[i]) > hight) {
				hight = abs(Data[i]);
				pos = i;
				if (abs(Data[i + 1]) <= hight && abs(Data[i + 2]) < hight && hight > this->m_adc_settings->threshold*this->m_adc_settings->gain) {
					PulseIn->Peaks[pk].dpa_max1_pos = double(pos) + 0.5 * (abs(Data[pos + 1]) - abs(Data[pos - 1])) / (2 * abs(Data[pos]) - abs(Data[pos - 1]) - abs(Data[pos + 1]));
					PulseIn->Peaks[pk].dpa_max1_hight = parabola->get_y_at_x(Data, pLength, PulseIn->Peaks[pk].dpa_max1_pos);
					break;
				}
			}
		}

		// find first min if max1 is found                         
		if (PulseIn->Peaks[pk].dpa_max1_pos > 0) {
			for (int i = pos; i < PulseIn->Peaks[pk].Stoppos; i++) {
				if (abs(Data[i]) < hight) {
					hight = abs(Data[i]);
					pos = i;
					if (abs(Data[i + 1]) >= hight && abs(Data[i + 2]) >= hight) {
						PulseIn->Peaks[pk].dpa_min_pos = double(pos) + 0.5 * (abs(Data[pos + 1]) - abs(Data[pos - 1])) / (2 * abs(Data[pos]) - abs(Data[pos - 1]) - abs(Data[pos + 1]));
						PulseIn->Peaks[pk].dpa_min_hight = parabola->get_y_at_x(Data, pLength, PulseIn->Peaks[pk].dpa_min_pos);
						break;
					}
				}
			}
		}

		// find second max if min1 is found
		if (PulseIn->Peaks[pk].dpa_min_pos > 0) {
			for (int i = pos; i < PulseIn->Peaks[pk].Stoppos; ++i) {
				if (abs(Data[i]) > hight) {
					hight = abs(Data[i]);
					pos = i;
					if (abs(Data[i + 1]) <= hight && abs(Data[i + 2]) <= hight && hight > this->m_adc_settings->threshold*this->m_adc_settings->gain) {
						PulseIn->Peaks[pk].dpa_max2_pos = double(pos) + 0.5 * (abs(Data[pos + 1]) - abs(Data[pos - 1])) / (2 * abs(Data[pos]) - abs(Data[pos - 1]) - abs(Data[pos + 1]));
						PulseIn->Peaks[pk].dpa_max2_hight = parabola->get_y_at_x(Data, pLength, PulseIn->Peaks[pk].dpa_max2_pos);
						PulseIn->Peaks[pk].dpa_is_potential_double = true;
						break;
					}
				}
			}
		}

		// bacon
		// problem.. if two peaks are close we change stop to second peak.. but what if one puls is double pulse.. 


		// find slope leading and falling edge
		if (PulseIn->Peaks[pk].dpa_is_potential_double) {
			if (PulseIn->Peaks[pk].dpa_max1_pos > 0 && PulseIn->Peaks[pk].dpa_max2_pos > 0) {

				// fwhm left (1) with max 1
				__int32 fwhm_left = 0;// int(curEvent->channels[ch]->Pulses[pu].Peaks[pk].dpa_max1_pos);
				double half_max_1 = PulseIn->Peaks[pk].dpa_max1_hight*0.5;
				for (fwhm_left = int(PulseIn->Peaks[pk].dpa_max1_pos); fwhm_left > PulseIn->Peaks[pk].Startpos; --fwhm_left) {
					if (Data[fwhm_left] <= half_max_1)
						break;
				}

				//pos half left
				double xbez[4];
				double ybez[4];
				xbez[0] = fwhm_left - 1;	ybez[0] = Data[fwhm_left - 1];
				xbez[1] = fwhm_left;		ybez[1] = Data[fwhm_left];
				xbez[2] = fwhm_left + 1;	ybez[2] = Data[fwhm_left + 1];
				xbez[3] = fwhm_left + 2;	ybez[3] = Data[fwhm_left + 2];
				bezier->bezier4_make_curve_parameters(xbez, ybez, 1.);
				PulseIn->Peaks[pk].dpa_fwhm1_pos = bezier->get_x_at_y_bezier4(half_max_1);

				//slope for mean pulse
				double mLeft, cLeft;
				double ylin[4];
				ylin[0] = Data[fwhm_left - 1];
				ylin[1] = Data[fwhm_left];
				ylin[2] = Data[fwhm_left + 1];
				ylin[3] = Data[fwhm_left + 2];
				linearRegression(4, xbez, ylin, mLeft, cLeft);
				PulseIn->Peaks[pk].dpa_fwhm1_slope = mLeft;

				//fwhm right (2) with max2
				int fwhm_right = 0;
				double half_max_2 = PulseIn->Peaks[pk].dpa_max2_hight*0.5;
				for (fwhm_right = int(PulseIn->Peaks[pk].dpa_max2_pos); fwhm_right < PulseIn->Peaks[pk].Stoppos; ++fwhm_right) {
					if (Data[fwhm_right] <= half_max_2)
						break;
				}
				double m_1 = abs(Data[fwhm_right]) - abs(Data[fwhm_right - 1]);
				double a_1 = abs(Data[fwhm_right - 1]) - m_1 * (fwhm_right - 1);
				double m_2 = (abs(Data[fwhm_right + 1]) - abs(Data[fwhm_right - 2])) / 3;
				double a_2 = abs(Data[fwhm_right - 2]) - m_2 * (fwhm_right - 2);
				PulseIn->Peaks[pk].dpa_fwhm2_pos = ((half_max_2 - a_1) / m_1 + (half_max_2 / 2 - a_2) / m_2) / 2;

				//slope falling edge for mean pulse
				double mRight, cRight;
				ylin[0] = Data[fwhm_right - 2];
				ylin[1] = Data[fwhm_right - 1];
				ylin[2] = Data[fwhm_right];
				ylin[3] = Data[fwhm_right + 1];
				linearRegression(4, xbez, ylin, mRight, cRight);
				PulseIn->Peaks[pk].dpa_fwhm2_slope = mRight;
			}
		}
	}
}

void ADC_analysis::collect_ADCinfos(Pulse *pu) {
	if (pu->NbrPeaks == 1) {
		if (!this->m_adc_settings->height_hist) {
				char hist_name[100];
				sprintf(hist_name, "ch%d_height", pu->Pu_from_Ch);
				this->m_adc_settings->height_hist = new H1D(hist_name, hist_name, 4000, -2000, 2000);
				sprintf(hist_name, "ch%d_fwhm", pu->Pu_from_Ch);
				this->m_adc_settings->fwhm_hist = new H1D(hist_name, hist_name, 250, 0, 50);
				sprintf(hist_name, "ch%d_slope", pu->Pu_from_Ch);
				this->m_adc_settings->slope_hist = new H1D(hist_name, hist_name, 1100, -1000, 10000);
				sprintf(hist_name, "ch%d_slope_fe", pu->Pu_from_Ch);
				this->m_adc_settings->slope_fe_hist = new H1D(hist_name, hist_name, 1100, -10000, 1000);
				sprintf(hist_name, "ch%d_integral", pu->Pu_from_Ch);
				this->m_adc_settings->integral_hist = new H1D(hist_name, hist_name, 1100, -1000, 1000000);
			}

			this->m_adc_settings->height_hist->Fill(	pu->Peaks[0].Height);
			this->m_adc_settings->fwhm_hist->Fill(		pu->Peaks[0].Fwhm);
			this->m_adc_settings->slope_hist->Fill(		pu->Peaks[0].Slope);
			this->m_adc_settings->slope_fe_hist->Fill(	pu->Peaks[0].Slope_fe);
			this->m_adc_settings->integral_hist->Fill(	pu->Peaks[0].Integral);		
		}
	}
void ADC_analysis::evaluate_ADCinfos() {
	if (m_adc_settings->height_hist) {
		this->GetRangeHisto(m_adc_settings->height_hist, &m_adc_settings->height_min, &m_adc_settings->height_max);
		this->GetRangeHisto(m_adc_settings->fwhm_hist, &m_adc_settings->fwhm_min, &m_adc_settings->fwhm_max);
		this->GetRangeHisto(m_adc_settings->slope_hist, &m_adc_settings->slope_min, &m_adc_settings->slope_max);
		this->GetRangeHisto(m_adc_settings->slope_fe_hist, &m_adc_settings->slope_fe_min, &m_adc_settings->slope_fe_max);
		this->GetRangeHisto(m_adc_settings->integral_hist, &m_adc_settings->integral_min, &m_adc_settings->integral_max);
	}
}

bool ADC_analysis::GetMeanpulse(Pulse* MeanPulse, Pulse *pu, int pk) {
	// find max of meanpulse (search around center of pulse)
	double mean_max_array[10];
	double max_height = 0;
	int max_pos = 0;
	int j = 0;
	int pulse_center = int(0.5*MEAN_PULSE_LENGTH);
	for (int i = pulse_center - 6; i < pulse_center + 4; i++) {	// find maximum around center
		mean_max_array[j] = m_adv_meanpulse->interpolate3D(pu->Peaks[pk].dpa_fwhm1_slope, pu->Peaks[pk].dpa_max1_hight / m_adc_settings->gain, i);
		if (max_height < mean_max_array[j]) {
			max_height = mean_max_array[j];
			max_pos = j;
		}
		j++;
	}

	// calculate exact position of maximum
	double mean_max_pos = double(max_pos) + 0.5 * (abs(mean_max_array[max_pos + 1]) - abs(mean_max_array[max_pos - 1])) / (2 * abs(mean_max_array[max_pos]) - abs(mean_max_array[max_pos - 1]) - abs(mean_max_array[max_pos + 1]));
	mean_max_pos += pulse_center - 6;

	//double mein_max2 = parabola->get_parabola_vertex(mean_max_array[max_pos - 1], mean_max_array[max_pos], mean_max_array[max_pos + 1], mean_max_array[max_pos + 2], mean_max_pos-double(max_pos));

	// create meanpulse
	MeanPulse->DataLength = pu->DataLength;
	MeanPulse->Pu_from_Ch = pu->Pu_from_Ch;

	if (0 == max_height) {
		MeanPulse->is_Meanpulse = false;
		return false;
	}

	MeanPulse->Waveform[0] = m_adv_meanpulse->interpolate3D(pu->Peaks[pk].dpa_fwhm1_slope, pu->Peaks[pk].dpa_max1_hight / m_adc_settings->gain, mean_max_pos - pu->Peaks[pk].dpa_max1_pos);
	double max_height_2 = MeanPulse->Waveform[0];
	int max_pos_2 = 0;
	for (int i = 1; i < pu->DataLength; i++) {
		MeanPulse->Waveform[i] = m_adv_meanpulse->interpolate3D(pu->Peaks[pk].dpa_fwhm1_slope, pu->Peaks[pk].dpa_max1_hight / m_adc_settings->gain, mean_max_pos + i - pu->Peaks[pk].dpa_max1_pos);
		if (max_height_2 < MeanPulse->Waveform[i]) {
			max_height_2 = MeanPulse->Waveform[i];
			max_pos_2 = i;
		}
	}

	double scale = pu->Peaks[pk].dpa_max1_hight / max_height_2;

	for (int i = 0; i < pu->DataLength; i++) {
		MeanPulse->Waveform[i] = MeanPulse->Waveform[i] * scale;
	}

	MeanPulse->is_Meanpulse = true;
	return true;
}
void ADC_analysis::accumulate_Meanpulse(Pulse* pu) {
	//prepair
	//leading edge
	if (!m_adv_meanpulse) {
		m_adv_meanpulse = new Trilinear();
		m_adv_meanpulse->createMatrix(100, m_adc_settings->slope_min, m_adc_settings->slope_max, 50, m_adc_settings->height_min, m_adc_settings->height_max, MEAN_PULSE_LENGTH, 0, MEAN_PULSE_LENGTH);
	}
	
	// normal meanpulse
	if (!m_meanpulse) {
		m_meanpulse = new Pulse();
		m_meanpulse->DataLength = 100;
		m_meanpulse->Pu_from_Ch = pu->Pu_from_Ch;
		m_meanpulse->is_Meanpulse = 1;
		m_meanpulse->meanpulse_cnt = 0;
		for (int i = 0; i < MAX_DATALENGHT; i++) {
			m_meanpulse->Waveform[i] = 0;
		}
	}
	
	
	//falling edge
	if (!m_adv_meanpulse_falling) {
		m_adv_meanpulse_falling = new Trilinear();
		//m_adv_meanpulse_falling[ch]->createMatrix(100, m_adc_info[ch].slope_fe_min, m_adc_info[ch].slope_fe_max, 50, m_adc_info[ch].height_min, m_adc_info[ch].height_max, MEAN_PULSE_LENGTH, 0, MEAN_PULSE_LENGTH);
		m_adv_meanpulse_falling->createMatrix(100, m_adc_settings->slope_min, m_adc_settings->slope_max, 50, m_adc_settings->height_min, m_adc_settings->height_max, MEAN_PULSE_LENGTH, 0, MEAN_PULSE_LENGTH);
	}
	
	//accumulate
	if (pu->NbrPeaks == 1) {
		int slope_bin = 0, height_bin = 0;
		slope_bin = int((pu->Peaks[0].Slope - m_adv_meanpulse->x_min) / (m_adv_meanpulse->x_range) * 10 + 0.5);
		height_bin = int((pu->Peaks[0].Height - m_adv_meanpulse->y_min) / (m_adv_meanpulse->y_range) * 10 + 0.5);
	
		int slope_fe_bin = 0;
		slope_fe_bin = int((pu->Peaks[0].Slope - m_adv_meanpulse_falling->x_min) / (m_adv_meanpulse_falling->x_range) * 10 + 0.5);
	
		//m_adv_meanpulse_hist[ch]->Fill(newEvent->channels[ch]->Pulses[0].Peaks[0].Slope, newEvent->channels[ch]->Pulses[0].Peaks[0].Height);
	
		//leading edge
		if (slope_bin >= 0 && slope_bin < m_adv_meanpulse->x_bins && height_bin >= 0 && height_bin < m_adv_meanpulse->y_bins) {
			int offset = -MEAN_PULSE_LENGTH / 2;
			for (int i = 0; i < m_adv_meanpulse->z_bins; i++) {
				double value = pu->GetWaveformAt(pu->Peaks[0].Time_ss + offset);
				if (value > -1e10) {
					m_adv_meanpulse->fillMatrix(pu->Peaks[0].Slope, pu->Peaks[0].Height, i, value / m_adc_settings->gain / pu->Peaks[0].Height);
					//m_adv_meanpulse_cnt[slope_bin][height_bin]++;
	
					// normal meanpulse
					m_meanpulse->Waveform[i] += value / m_adc_settings->gain / pu->Peaks[0].Height;
					m_meanpulse->meanpulse_cnt++;
				}
				offset++;
			}
		}
	
		//falling edge
		if (slope_fe_bin >= 0 && slope_fe_bin < m_adv_meanpulse->x_bins && height_bin >= 0 && height_bin < m_adv_meanpulse->y_bins) {
			int offset = -MEAN_PULSE_LENGTH / 2;
			for (int i = 0; i < m_adv_meanpulse_falling->z_bins; i++) {
				double value = pu->GetWaveformAt(pu->Peaks[0].Time_ss + offset);
				if (value > -1e10) {
					m_adv_meanpulse_falling->fillMatrix(pu->Peaks[0].Slope, pu->Peaks[0].Height, i, value / m_adc_settings->gain / pu->Peaks[0].Height);
					//m_adv_meanpulse_falling_cnt[slope_bin][height_bin]++;
				}
				offset++;
			}
		}
	}
}
void ADC_analysis::evaluate_Meanpulses() {
	if (m_adv_meanpulse) {
		for (int x = 0; x < m_adv_meanpulse->x_bins - 1; x++) {
			for (int y = 0; y < m_adv_meanpulse->y_bins - 1; y++) {
				if (m_adv_meanpulse->m_matrix[x][y]) {
					double max = 0;
					for (int z = 0; z < m_adv_meanpulse->z_bins; z++) {
						if (max < m_adv_meanpulse->m_matrix[x][y][z].val)
							max = m_adv_meanpulse->m_matrix[x][y][z].val;
					}
					for (int z = 0; z < m_adv_meanpulse->z_bins; z++) {
						m_adv_meanpulse->m_matrix[x][y][z].val = m_adv_meanpulse->m_matrix[x][y][z].val / max;
					}
				}
			}
		}
	}

	if (m_adv_meanpulse_falling) {
		for (int x = 0; x < m_adv_meanpulse_falling->x_bins - 1; x++) {
			for (int y = 0; y < m_adv_meanpulse_falling->y_bins - 1; y++) {
				if (m_adv_meanpulse_falling->m_matrix[x][y]) {
					double max = 0;
					for (int z = 0; z < m_adv_meanpulse_falling->z_bins; z++) {
						if (max < m_adv_meanpulse_falling->m_matrix[x][y][z].val)
							max = m_adv_meanpulse_falling->m_matrix[x][y][z].val;
					}
					for (int z = 0; z < m_adv_meanpulse_falling->z_bins; z++) {
						m_adv_meanpulse_falling->m_matrix[x][y][z].val = m_adv_meanpulse_falling->m_matrix[x][y][z].val / max;
					}
				}
			}
		}
	}

	if (m_meanpulse) {
		for (int i = 0; i < m_meanpulse->DataLength; i++) {
			m_meanpulse->Waveform[i] /= double(m_meanpulse->meanpulse_cnt);
		}
	}
}

void ADC_analysis::collect_ADCandDPA_PulseError(Pulse *mixedPulse) {
	if (!mixedPulse)
		return;

	if (!mixedPulse->mixed_Event) 
		return; 


	if (!m_adc_settings->hist_ADCpuls1_error || !m_adc_settings->hist_ADCpuls2_error) {
		char hist_name[100];
		sprintf(hist_name, "ch%d_ADCpuls1_error", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_ADCpuls1_error = new H2D(hist_name, hist_name, 1001, -10, 10, 100, 0.5, 100.5);
		sprintf(hist_name, "ch%d_ADCpuls2_error", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_ADCpuls2_error = new H2D(hist_name, hist_name, 1001, -10, 10, 100, 0.5, 100.5);
		sprintf(hist_name, "ch%d_DPApuls1_error", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_DPApuls1_error = new H2D(hist_name, hist_name, 1001, -10, 10, 100, 0.5, 100.5);
		sprintf(hist_name, "ch%d_DPApuls2_error", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_DPApuls2_error = new H2D(hist_name, hist_name, 1001, -10, 10, 100, 0.5, 100.5);

		sprintf(hist_name, "ch%d_Integal_DPA", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_integral_dpa = new H1D(hist_name, hist_name, 1100, -1000, 1000000);
		sprintf(hist_name, "ch%d_fwhm_DPA", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_fwhm_dpa = new H1D(hist_name, hist_name, 250, 0, 50);
		sprintf(hist_name, "ch%d_height_DPA", mixedPulse->Pu_from_Ch);
		m_adc_settings->hist_height_dpa = new H1D(hist_name, hist_name, 4000, -2000, 2000);
	}

	// ADC_part
	mixedPulse->resetPeaks();
	FindPeaks(mixedPulse);
	if (mixedPulse->NbrPeaks >= 1) {
		m_adc_settings->hist_ADCpuls1_error->Fill(mixedPulse->mixed_Time1_ss - mixedPulse->Peaks[0].Time_ss, mixedPulse->mixed_dTime_ss);
	}
	if (mixedPulse->NbrPeaks >= 2) {
		m_adc_settings->hist_ADCpuls2_error->Fill(mixedPulse->mixed_Time2_ss - mixedPulse->Peaks[1].Time_ss, mixedPulse->mixed_dTime_ss);
	}

	// if only one puls is found but we know its a double pulse
	if (mixedPulse->NbrPeaks == 1) {
		m_adc_settings->hist_integral_dpa->Fill(mixedPulse->Peaks[0].Integral);
		m_adc_settings->hist_fwhm_dpa->Fill(mixedPulse->Peaks[0].Fwhm);
		m_adc_settings->hist_height_dpa->Fill(mixedPulse->Peaks[0].Height_mV);
	}
	
	// DPA_part	
	FindDoublePulse(mixedPulse);
	if (m_adc_settings->use_dpa && mixedPulse->mixed_Event && mixedPulse->used_in_DPA) {
		if (mixedPulse->NbrPeaks >= 1) {
			m_adc_settings->hist_DPApuls1_error->Fill(mixedPulse->mixed_Time1_ss - mixedPulse->Peaks[0].Time_ss, mixedPulse->mixed_dTime_ss);
		}
		if (mixedPulse->NbrPeaks >= 2) {
			m_adc_settings->hist_DPApuls2_error->Fill(mixedPulse->mixed_Time2_ss - mixedPulse->Peaks[1].Time_ss, mixedPulse->mixed_dTime_ss);
		}
	}
}
void ADC_analysis::GetMixedPulse(Pulse *pu1, Pulse *pu2, Pulse *mixedPulse, int dSS) {
	if (1 == pu1->NbrPeaks && 1 == pu2->NbrPeaks) {
		mixedPulse->resetPulse();
		mixedPulse->DataLength = pu1->DataLength + pu2->DataLength;
		
		if (mixedPulse->DataLength > MAX_DATALENGHT) { return; }

		mixedPulse->timestamp = 0;
		mixedPulse->Pu_from_Ch = pu1->Pu_from_Ch;
		mixedPulse->mixed_Event = 0;

		for (int i = 0; i < mixedPulse->DataLength; i++) {
			mixedPulse->Waveform[i] = 0.;
		}

		double dTime_samplesteps = pu2->Peaks[0].Time_ss - pu1->Peaks[0].Time_ss + dSS;

		if (dSS >= 0 && dTime_samplesteps >= 0) {
			for (int i = 0; i < pu1->DataLength; i++){
				mixedPulse->Waveform[i] = pu1->Waveform[i];
			}
			for (int i = 0; i < pu2->DataLength && i + dSS < mixedPulse->DataLength; i++){
				mixedPulse->Waveform[i + dSS] += pu2->Waveform[i];
			}

			mixedPulse->mixed_Time1_ss = pu1->Peaks[0].Time_ss;			// position stays the same
			mixedPulse->mixed_Time2_ss = pu2->Peaks[0].Time_ss + dSS;	// position is shifted by dsamplesteps
			mixedPulse->mixed_Event = 1;
			mixedPulse->mixed_dTime_ss = pu2->Peaks[0].Time_ss - pu2->Peaks[0].Time_ss + dSS;
		}
	}
}void ADC_analysis::evaluate_ADCandDPA_PulseError() {
// function to determine  if it is better to use normal adc analysis or dpa analysis..
	if (m_adc_settings->hist_ADCpuls1_error &&
		m_adc_settings->hist_ADCpuls2_error &&
		m_adc_settings->hist_DPApuls1_error &&
		m_adc_settings->hist_DPApuls1_error) {

		if (!m_adc_settings->hist_ADCvsDPA_puls1 &&
			!m_adc_settings->hist_ADCvsDPA_puls2) {
			char hist_name[100];
			sprintf(hist_name, "ch%d_ADCvsDPA_puls1", mixedPulse->Pu_from_Ch);
			m_adc_settings->hist_ADCvsDPA_puls1 = new H1D(hist_name, hist_name, m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetNbins(), m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetXlow(), m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetXup());
			sprintf(hist_name, "ch%d_ADCvsDPA_puls2", mixedPulse->Pu_from_Ch);
			m_adc_settings->hist_ADCvsDPA_puls2 = new H1D(hist_name, hist_name, m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetNbins(), m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetXlow(), m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetXup());

			for (int y = 1; y <= m_adc_settings->hist_ADCpuls1_error->GetYaxis()->GetNbins(); y++) {
				double error_adc_puls1 = 0, error_dpa_puls1 = 0;
				int error_adc_n_puls1 = 0, error_dpa_n_puls1 = 0;;
				for (int x = 1; x <= m_adc_settings->hist_ADCpuls1_error->GetXaxis()->GetNbins(); x++) {
					error_adc_n_puls1 += int(m_adc_settings->hist_ADCpuls1_error->GetBinContent(x, y));
					error_adc_puls1 += m_adc_settings->hist_ADCpuls1_error->GetBinContent(x, y) * m_adc_settings->hist_ADCpuls1_error->GetXaxis()->GetBinCenter(x) * m_adc_settings->hist_ADCpuls1_error->GetXaxis()->GetBinCenter(x);

					error_dpa_n_puls1 += int(m_adc_settings->hist_DPApuls1_error->GetBinContent(x, y));
					error_dpa_puls1 += m_adc_settings->hist_DPApuls1_error->GetBinContent(x, y) * m_adc_settings->hist_DPApuls1_error->GetXaxis()->GetBinCenter(x) * m_adc_settings->hist_DPApuls1_error->GetXaxis()->GetBinCenter(x);
				}

				if (double(error_adc_n_puls1)*sqrt(error_adc_puls1) > 0 && double(error_dpa_n_puls1)*sqrt(error_dpa_puls1) > 0) {
					error_adc_puls1 = 1. / double(error_adc_n_puls1)*sqrt(error_adc_puls1);
					error_dpa_puls1 = 1. / double(error_dpa_n_puls1)*sqrt(error_dpa_puls1);
					m_adc_settings->hist_ADCvsDPA_puls1->SetBinContent(y, error_adc_puls1 - error_dpa_puls1);
				}

				double error_adc_puls2 = 0, error_dpa_puls2 = 0;
				int error_adc_n_puls2 = 0, error_dpa_n_puls2 = 0;;
				for (int x = 1; x <= m_adc_settings->hist_ADCpuls2_error->GetXaxis()->GetNbins(); x++) {
					error_adc_n_puls2 += int(m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y));
					error_adc_puls2 += m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y) * m_adc_settings->hist_ADCpuls2_error->GetXaxis()->GetBinCenter(x) * m_adc_settings->hist_ADCpuls2_error->GetXaxis()->GetBinCenter(x);

					error_dpa_n_puls2 += int(m_adc_settings->hist_DPApuls2_error->GetBinContent(x, y));
					error_dpa_puls2 += m_adc_settings->hist_DPApuls2_error->GetBinContent(x, y) * m_adc_settings->hist_DPApuls2_error->GetXaxis()->GetBinCenter(x) * m_adc_settings->hist_DPApuls2_error->GetXaxis()->GetBinCenter(x);
				}

				if (double(error_adc_n_puls2)*sqrt(error_adc_puls2) > 0 && double(error_dpa_n_puls2)*sqrt(error_dpa_puls2) > 0) {
					error_adc_puls2 = 1. / double(error_adc_n_puls2)*sqrt(error_adc_puls2);
					error_dpa_puls2 = 1. / double(error_dpa_n_puls2)*sqrt(error_dpa_puls2);
					m_adc_settings->hist_ADCvsDPA_puls2->SetBinContent(y, error_adc_puls2 - error_dpa_puls2);
				}
			}
		}
	}

	// find cfd risky settings..
	for (int i = m_adc_settings->hist_ADCvsDPA_puls2->GetXaxis()->GetNbins(); i > 1; i--) {
		if (m_adc_settings->hist_ADCvsDPA_puls2->GetBinContent(i - 1) > m_adc_settings->hist_ADCvsDPA_puls2->GetBinContent(i)) {
			m_adc_settings->DPA_is_favorable_below_ss = m_adc_settings->hist_ADCvsDPA_puls2->GetBinCenter(i);
			break;
		}
	}


	//GetRangeHisto(m_adc_info[ch].hist_integral_dpa, &m_adc_info[ch].integral_dpa_min, &m_adc_info[ch].integral_dpa_max);

}
/*
void ADC_analysis::evaluateADCpulses_error() {
	if (m_adc_settings->hist_ADCpuls2_error) {
		int bin_min = m_adc_settings->hist_ADCpuls2_error->GetXaxis()->FindBin(-0.125);	// 1 = 800 ps -> 0.2 = 160 ps
		int bin_max = m_adc_settings->hist_ADCpuls2_error->GetXaxis()->FindBin(+0.125);
		
		m_adc_settings->Cfd_risky = 0;

		for (int y = m_adc_settings->hist_ADCpuls2_error->GetYaxis()->GetNbins(); y > 1; y--) {
			int sum = 0;
			int sum_below_error = 0;
				
			for (int x = 1; x < m_adc_settings->hist_ADCpuls2_error->GetXaxis()->GetNbins(); x++) {
				sum += m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y);
			}
			for (int x = bin_min; x <= bin_max; x++) {
				sum_below_error += m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y);
			}
				
			if (sum > 0 && sum_below_error / sum < 0.98) {	// if not 98% in confidence interval, check also next row..
				int second_sum_below_error = 0;
				int second_sum = 0;
				for (int x = 1; x < m_adc_settings->hist_ADCpuls2_error->GetXaxis()->GetNbins(); x++) {
					second_sum += m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y-1);
				}
				for (int x = bin_min; x <= bin_max; x++) {
					second_sum_below_error += m_adc_settings->hist_ADCpuls2_error->GetBinContent(x, y-1);
				}
				if (second_sum > 0 && second_sum_below_error / second_sum < 0.98) {
					m_adc_settings->Cfd_risky = m_adc_settings->hist_ADCpuls2_error->GetYaxis()->GetBinCenter(y) + 0.5; // add 5 ss security
					break;
				}
			}		
		}
	}
}
*/

bool ADC_analysis::substractPulse(Pulse* substractedPulse, Pulse* original, Pulse* substract) {
	substractedPulse->DataLength = original->DataLength;
	substractedPulse->Pu_from_Ch = original->Pu_from_Ch;

	for (int i = 0; i < original->DataLength; i++) {
		substractedPulse->Waveform[i] = original->Waveform[i] - substract->Waveform[i];
	}
	

	return true;
}
void ADC_analysis::GetRangeHisto(H1D* hist, double *range_low, double *range_high, double percentage, double safety_per) {
	if (hist) {
		int xbins = hist->GetXaxis()->GetNbins();

		int sum = 0;
		for (int i = 1; i <= xbins; i++) {
			sum += int(hist->GetBinContent(i));
		}

		if (0 == sum) {
			OutputDebugString("\nError - GetHistoRange not possible for empty hist.");
			*range_low	= 0;
			*range_high = 0;
			return;
		}

		int sum_center = 0;
		int bin_center = 0;
		for (int i = 1; sum_center <= sum / 2; i++) {
			sum_center += int(hist->GetBinContent(i));
			bin_center = i;
		}

		double per = percentage / 100 / 2;

		// go left 
		double sum_temp = hist->GetBinContent(bin_center) / 2;
		for (int i = bin_center - 1; i > 0; i--) {
			sum_temp += hist->GetBinContent(i);
			if (sum_temp / sum >= per) {
				*range_low = hist->GetBinCenter(i);
				break;
			}
		}

		//go right
		sum_temp = hist->GetBinContent(bin_center) / 2;
		for (int i = bin_center + 1; i < xbins; i++) {
			sum_temp += hist->GetBinContent(i);
			if (sum_temp / sum >= per) {
				*range_high = hist->GetBinCenter(i);
				break;
			}
		}

		// give safty 
		*range_low = *range_low - (safety_per/100) * abs(*range_high - *range_low);
		*range_high = *range_high + (safety_per/100) * abs(*range_high - *range_low);
	}
}
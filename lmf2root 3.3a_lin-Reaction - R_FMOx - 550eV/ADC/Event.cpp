#include "../Src/Ueberstruct.h"
#include "Event.h"

void linearRegression(const int nbrPoints, const double x[], const double y[], double &m, double &c);
double findXForGivenY(const double * x, const double * coeff, const double Y, const double Start);
void createNewtonPolynomial(const double * x, const double * y, double * coeff);
double evalNewtonPolynomial(const double * x, const double * coeff, double X);

TCanvas* Canv(const char *name, int hor, int ver, int hor_size = 500, int ver_size = 200)
{
	TCanvas * MyC = new TCanvas(name,"MyCanvas",-1,1,hor*hor_size,ver*ver_size);
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

Event::Event(/*Ueberstruct *Ueber, ADCinfos *adc_info*/){
//	Ueber_ptr = Ueber;
	id = 0;
//	adc_set = adc_info;

	for(int i = 0; i < NUM_CHANNELS; i++){
		channels[i] = nullptr;
	}

	parabola = new parabola_spline;
	bezier = new bezier_class;
}
Event::~Event(){
	for (int i = 0; i < NUM_CHANNELS; i++) {
		if (channels[i]) {
			delete channels[i];
			channels[i] = 0;
		}
	}

	delete parabola; parabola = 0;
	delete bezier; bezier = 0;
}

void Event::resetEvent(){
	for(int ch = 0; ch < NUM_CHANNELS; ch++){
		if(channels[ch]){
			for(int pu = 0; pu < MAX_NBR_PULSES; pu++){
				//channels[ch]->Pulses[pu].NbrPeaks = 0;
				//channels[ch]->Pulses[pu].DataLength = 0;
				//memset(channels[ch]->Pulses[pu].Waveform, 0, sizeof(channels[ch]->Pulses[pu].Waveform));
				channels[ch]->Pulses[pu].resetPeaks();
				channels[ch]->Pulses[pu].resetPulse();
			}
			channels[ch]->NbrPulses = 0;
		}
	}
}

Pulse::Pulse(){
	NbrPeaks				= 0;
	DataLength				= 0;					
	timestamp				= 0;
	is_TDC_signal			= false;
	
	used_in_DPA				= 0;
	is_Meanpulse			= 0;

	//ZeroSurpressedWaveform	= new int[MAX_DATALENGHT];
	Waveform				= new double[MAX_DATALENGHT];
	//memset(Waveform, 0, MAX_DATALENGHT*sizeof(double));

	Peaks					= new Peak[MAX_NBR_PEAKS];  //nullptr;//new Peak[MAX_NBR_PEAKS];
	Peaks_backup			= nullptr;

	NbrPeaks_backup			= 0;
}

Pulse::~Pulse(){
	NbrPeaks				= 0;
	NbrPeaks_backup			= 0;
	DataLength				= 0;
	timestamp				= 0;

	used_in_DPA				= 0;
	is_Meanpulse			= 0;

	if (Waveform) delete[] Waveform; Waveform = 0;
	if (Peaks) delete[] Peaks; Peaks = 0;
	if (Peaks_backup) delete[] Peaks_backup; Peaks_backup = 0;
}

Pulse* Pulse::Pulse_cpy() {
	Pulse *copy = new Pulse();
	copy->Pu_from_Ch = this->Pu_from_Ch;
	copy->DataLength = this->DataLength;
	copy->timestamp = this->timestamp;

	copy->used_in_DPA = this->used_in_DPA;
	copy->is_Meanpulse = this->is_Meanpulse;

	//this->Waveform = new double[MAX_DATALENGHT];
	for (int i = 0; i < copy->DataLength; i++)
		copy->Waveform[i] = this->Waveform[i];

	//Peaks = new Peak[MAX_NBR_PEAKS];
	copy->NbrPeaks = this->NbrPeaks;
	for (int pk = 0; pk < NbrPeaks; pk++) {
		Peak_cpy(&this->Peaks[pk], &copy->Peaks[pk]);
	}

	copy->Peaks_backup = nullptr;

	copy->NbrPeaks_backup = this->NbrPeaks_backup;
	return copy;
}

void Pulse::Pulse_init(Pulse *orig) {
	this->DataLength = orig->DataLength;
	this->timestamp = orig->timestamp;

	this->used_in_DPA = orig->used_in_DPA;
	this->is_Meanpulse = orig->is_Meanpulse;

	//this->Waveform = new double[MAX_DATALENGHT];
	for (int i = 0; i < this->DataLength; i++)
		this->Waveform[i] = orig->Waveform[i];

	//Peaks = new Peak[MAX_NBR_PEAKS];
	this->NbrPeaks = orig->NbrPeaks;
	for (int pk = 0; pk < this->NbrPeaks; pk++) {
		Peak_cpy(&orig->Peaks[pk], &this->Peaks[pk]);
	}

	this->Peaks_backup = nullptr;

	this->NbrPeaks_backup = orig->NbrPeaks_backup;
}

void Pulse::resetPeaks() {
	//for (int i = 0; i < this->NbrPeaks; i++) {
	//	this->Peaks[i].TimefromFit = 0;
	//	this->Peaks[i].Time = 0;
	//	this->Peaks[i].Cfd = 0;
	//	this->Peaks[i].Com = 0;
	//	this->Peaks[i].Time_ns = 0;
	//	this->Peaks[i].Time_ss = 0;
	//	this->Peaks[i].Polarity = 0;
	//	this->Peaks[i].Slope = 0;
	//	this->Peaks[i].Slope_fe = 0;
	//	this->Peaks[i].Maxpos = 0;
	//	this->Peaks[i].Maximum = 0;
	//	this->Peaks[i].Height = 0;
	//	this->Peaks[i].HeightfromFit = 0;
	//	this->Peaks[i].Height_mV = 0;
	//	this->Peaks[i].Height_ss = 0;
	//	this->Peaks[i].Fwhm = 0;
	//	this->Peaks[i].Width = 0;
	//	this->Peaks[i].PosHalfLeft = 0;
	//	this->Peaks[i].PosHalfRight = 0;
	//	this->Peaks[i].Integral = 0;
	//	this->Peaks[i].Startpos = 0;
	//	this->Peaks[i].Stoppos = 0;
	//	this->Peaks[i].scalefactor = 0;
	//	this->Peaks[i].calcbinoffset = 0;
	//	this->Peaks[i].runtime = 0;
	//	this->Peaks[i].is_valid = 0;
	//  this->Peaks[i].is_TDC_signal = false;
	//	this->Peaks[i].dpa_is_potential_double = 0;
	//	this->Peaks[i].dpa_max1_pos = 0;
	//	this->Peaks[i].dpa_max1_hight = 0;
	//	this->Peaks[i].dpa_min_pos = 0;
	//	this->Peaks[i].dpa_min_hight = 0;
	//	this->Peaks[i].dpa_max2_pos = 0;
	//	this->Peaks[i].dpa_max2_hight = 0;
	//	this->Peaks[i].dpa_fwhm1_pos = 0;
	//	this->Peaks[i].dpa_fwhm1_hight = 0;
	//	this->Peaks[i].dpa_fwhm1_slope = 0;
	//	this->Peaks[i].dpa_fwhm2_pos = 0;
	//	this->Peaks[i].dpa_fwhm2_hight = 0;
	//	this->Peaks[i].dpa_fwhm2_slope = 0;
	//	this->Peaks[i].left_curved_cnt = 0;
	//}
	
	NbrPeaks = 0;
	NbrPeaks_backup = 0; // BACON no real reset for backup peaks..
}

void Pulse::resetPulse(){
	Pu_from_Ch		= 0;
	DataLength		= 0;
	timestamp		= 0;
	
	//memset(Waveform, 0, MAX_DATALENGHT*sizeof(Waveform));
	
	resetPeaks();	

	used_in_DPA		= 0;
	is_Meanpulse    = 0;
	
	mixed_dTime_ss	= 0;
	mixed_Time1_ss	= 0;
	mixed_Time2_ss	= 0;
	mixed_Event		= 0;
}

double Pulse::parabola_fit(double y1, double y2, double y3, double x_pos)
{
	double a = 0.5*(y1 + y3) - y2;
	double b = 2 * y2 - 1.5*y1 - 0.5*y3;
	return (a*x_pos + b)*x_pos + y1;
}

// smooth interpolation between discrete values relative to first peak
double Pulse::GetWaveformAt(double x){
	if (x < 0. || x > this->DataLength) return -1.e200;
	if (this->DataLength < 3) return -1.e200;

	if (x <= 1.) {
		return parabola_fit(this->Waveform[0], this->Waveform[1], this->Waveform[2], x);
	}

	if (x >= this->DataLength - 2) {
		return parabola_fit(this->Waveform[this->DataLength - 3], this->Waveform[this->DataLength - 2], this->Waveform[this->DataLength - 1], x - this->DataLength + 3);
	}

	__int32 index = __int32(x);
	double dx = x - index;
	double yi = this->Waveform[index];
	double yip1 = this->Waveform[index + 1];
	double y1 = parabola_fit(this->Waveform[index - 1], yi, yip1, dx + 1.);
	double y2 = parabola_fit(yi, yip1, this->Waveform[index + 2], dx);

	return y1*(1. - dx) + y2*dx;
} 

void Pulse::Peak_delete(int pk) {
	if (pk >= 0 && pk < NbrPeaks) {
		for (int i = pk; i < NbrPeaks - 1; i++) {
			Peak_cpy(i + 1, i);
		}
		NbrPeaks--;
	}
	else if (pk == NbrPeaks - 1) {
		NbrPeaks--;
	}
}

void Pulse::Peak_cpy(int pk_from, int pk_to) {
	if (pk_from >= 0 && pk_from < NbrPeaks && pk_to >= 0 && pk_to < NbrPeaks) {
		Peaks[pk_to].TimefromFit			= Peaks[pk_from].TimefromFit;
		Peaks[pk_to].Time					= Peaks[pk_from].Time;
		Peaks[pk_to].Time_ns				= Peaks[pk_from].Time_ns;
		Peaks[pk_to].Time_ss				= Peaks[pk_from].Time_ss;
		Peaks[pk_to].Cfd					= Peaks[pk_from].Cfd;
		Peaks[pk_to].Com					= Peaks[pk_from].Com;
		Peaks[pk_to].Polarity				= Peaks[pk_from].Polarity;
		Peaks[pk_to].Slope					= Peaks[pk_from].Slope;
		Peaks[pk_to].Maxpos					= Peaks[pk_from].Maxpos;
		Peaks[pk_to].Maximum				= Peaks[pk_from].Maximum;
		Peaks[pk_to].Height					= Peaks[pk_from].Height;
		Peaks[pk_to].HeightfromFit			= Peaks[pk_from].HeightfromFit;
		Peaks[pk_to].Width					= Peaks[pk_from].Width;
		Peaks[pk_to].Fwhm					= Peaks[pk_from].Fwhm;
		Peaks[pk_to].PosHalfLeft			= Peaks[pk_from].PosHalfLeft;
		Peaks[pk_to].PosHalfRight			= Peaks[pk_from].PosHalfRight;
		Peaks[pk_to].Integral				= Peaks[pk_from].Integral;
		Peaks[pk_to].Startpos				= Peaks[pk_from].Startpos;
		Peaks[pk_to].Stoppos				= Peaks[pk_from].Stoppos;
		Peaks[pk_to].scalefactor			= Peaks[pk_from].scalefactor;
		Peaks[pk_to].calcbinoffset			= Peaks[pk_from].calcbinoffset;
		Peaks[pk_to].runtime				= Peaks[pk_from].runtime;
		Peaks[pk_to].is_valid				= Peaks[pk_from].is_valid;
		Peaks[pk_to].is_TDC_signal			= Peaks[pk_from].is_TDC_signal;
		Peaks[pk_to].dpa_is_potential_double= Peaks[pk_from].dpa_is_potential_double;
		Peaks[pk_to].dpa_max1_pos			= Peaks[pk_from].dpa_max1_pos;
		Peaks[pk_to].dpa_max1_hight			= Peaks[pk_from].dpa_max1_hight;
		Peaks[pk_to].dpa_min_pos			= Peaks[pk_from].dpa_min_pos;
		Peaks[pk_to].dpa_min_hight			= Peaks[pk_from].dpa_min_hight;
		Peaks[pk_to].dpa_max2_pos			= Peaks[pk_from].dpa_max2_pos;
		Peaks[pk_to].dpa_max2_hight			= Peaks[pk_from].dpa_max2_hight;
		Peaks[pk_to].dpa_fwhm1_pos			= Peaks[pk_from].dpa_fwhm1_pos;
		Peaks[pk_to].dpa_fwhm1_hight		= Peaks[pk_from].dpa_fwhm1_hight;
		Peaks[pk_to].dpa_fwhm1_slope		= Peaks[pk_from].dpa_fwhm1_slope;
		Peaks[pk_to].dpa_fwhm2_pos			= Peaks[pk_from].dpa_fwhm2_pos;
		Peaks[pk_to].dpa_fwhm2_hight		= Peaks[pk_from].dpa_fwhm2_hight;
		Peaks[pk_to].dpa_fwhm2_slope		= Peaks[pk_from].dpa_fwhm2_slope;
		Peaks[pk_to].left_curved_cnt		= Peaks[pk_from].left_curved_cnt;
	}
}

void Pulse::Peak_cpy(Peak *pk_from, Peak *pk_to) {
	pk_to->TimefromFit = pk_from->TimefromFit;
	pk_to->Time = pk_from->Time;
	pk_to->Time_ns = pk_from->Time_ns;
	pk_to->Time_ss = pk_from->Time_ss;
	pk_to->Cfd = pk_from->Cfd;
	pk_to->Com = pk_from->Com;
	pk_to->Polarity = pk_from->Polarity;
	pk_to->Slope = pk_from->Slope;
	pk_to->Maxpos = pk_from->Maxpos;
	pk_to->Maximum = pk_from->Maximum;
	pk_to->Height = pk_from->Height;
	pk_to->HeightfromFit = pk_from->HeightfromFit;
	pk_to->Width = pk_from->Width;
	pk_to->Fwhm = pk_from->Fwhm;
	pk_to->PosHalfLeft = pk_from->PosHalfLeft;
	pk_to->PosHalfRight = pk_from->PosHalfRight;
	pk_to->Integral = pk_from->Integral;
	pk_to->Startpos = pk_from->Startpos;
	pk_to->Stoppos = pk_from->Stoppos;
	pk_to->scalefactor = pk_from->scalefactor;
	pk_to->calcbinoffset = pk_from->calcbinoffset;
	pk_to->runtime = pk_from->runtime;
	pk_to->is_valid = pk_from->is_valid;
	pk_to->dpa_is_potential_double = pk_from->dpa_is_potential_double;
	pk_to->dpa_max1_pos = pk_from->dpa_max1_pos;
	pk_to->dpa_max1_hight = pk_from->dpa_max1_hight;
	pk_to->dpa_min_pos = pk_from->dpa_min_pos;
	pk_to->dpa_min_hight = pk_from->dpa_min_hight;
	pk_to->dpa_max2_pos = pk_from->dpa_max2_pos;
	pk_to->dpa_max2_hight = pk_from->dpa_max2_hight;
	pk_to->dpa_fwhm1_pos = pk_from->dpa_fwhm1_pos;
	pk_to->dpa_fwhm1_hight = pk_from->dpa_fwhm1_hight;
	pk_to->dpa_fwhm1_slope = pk_from->dpa_fwhm1_slope;
	pk_to->dpa_fwhm2_pos = pk_from->dpa_fwhm2_pos;
	pk_to->dpa_fwhm2_hight = pk_from->dpa_fwhm2_hight;
	pk_to->dpa_fwhm2_slope = pk_from->dpa_fwhm2_slope;
	pk_to->left_curved_cnt = pk_from->left_curved_cnt;
}

// function will insert a new peak
void Pulse::Peak_insert(int pk, Peak &newPeak) {
	for (int i = this->NbrPeaks; i >= pk; i--) {
		Peak_cpy(i, i + 1);
	}
	this->NbrPeaks++;
	
	Peaks[pk].TimefromFit = newPeak.TimefromFit;
	Peaks[pk].Time = newPeak.Time;
	Peaks[pk].Time_ns = newPeak.Time_ns;
	Peaks[pk].Time_ss = newPeak.Time_ss;
	Peaks[pk].Cfd = newPeak.Cfd;
	Peaks[pk].Com = newPeak.Com;
	Peaks[pk].Polarity = newPeak.Polarity;
	Peaks[pk].Slope = newPeak.Slope;
	Peaks[pk].Maxpos = newPeak.Maxpos;
	Peaks[pk].Maximum = newPeak.Maximum;
	Peaks[pk].Height = newPeak.Height;
	Peaks[pk].HeightfromFit = newPeak.HeightfromFit;
	Peaks[pk].Fwhm = newPeak.Fwhm;
	Peaks[pk].Width = newPeak.Width;
	Peaks[pk].PosHalfLeft = newPeak.PosHalfLeft;
	Peaks[pk].PosHalfRight = newPeak.PosHalfRight;
	Peaks[pk].Startpos = newPeak.Startpos;
	Peaks[pk].Integral = newPeak.Integral;
	Peaks[pk].Stoppos = newPeak.Stoppos;
	Peaks[pk].scalefactor = newPeak.scalefactor;
	Peaks[pk].calcbinoffset = newPeak.calcbinoffset;
	Peaks[pk].runtime = newPeak.runtime;
	Peaks[pk].is_valid = newPeak.is_valid;
	Peaks[pk].is_TDC_signal = newPeak.is_TDC_signal;
	Peaks[pk].dpa_is_potential_double = newPeak.dpa_is_potential_double;
	Peaks[pk].dpa_max1_pos = newPeak.dpa_max1_pos;
	Peaks[pk].dpa_max1_hight = newPeak.dpa_max1_hight;
	Peaks[pk].dpa_min_pos = newPeak.dpa_min_pos;
	Peaks[pk].dpa_min_hight = newPeak.dpa_min_hight;
	Peaks[pk].dpa_max2_pos = newPeak.dpa_max2_pos;
	Peaks[pk].dpa_max2_hight = newPeak.dpa_max2_hight;
	Peaks[pk].dpa_fwhm1_pos = newPeak.dpa_fwhm1_pos;
	Peaks[pk].dpa_fwhm1_hight = newPeak.dpa_fwhm1_hight;
	Peaks[pk].dpa_fwhm1_slope = newPeak.dpa_fwhm1_slope;
	Peaks[pk].dpa_fwhm2_pos = newPeak.dpa_fwhm2_pos;
	Peaks[pk].dpa_fwhm2_hight = newPeak.dpa_fwhm2_hight;
	Peaks[pk].dpa_fwhm2_slope = newPeak.dpa_fwhm2_slope;
	Peaks[pk].left_curved_cnt = newPeak.left_curved_cnt;
}

// will copy all peaks into peaks_backup
void Pulse::Peak_backup() {
	this->NbrPeaks_backup = this->NbrPeaks;

	if(!Peaks_backup)
		Peaks_backup = new Peak[MAX_NBR_PEAKS];

	for (int pk = 0; pk < this->NbrPeaks_backup; pk++) {
		Peaks_backup[pk].TimefromFit = Peaks[pk].TimefromFit;
		Peaks_backup[pk].Time = Peaks[pk].Time;
		Peaks_backup[pk].Time_ns = Peaks[pk].Time_ns;
		Peaks_backup[pk].Time_ss = Peaks[pk].Time_ss;
		Peaks_backup[pk].Cfd = Peaks[pk].Cfd;
		Peaks_backup[pk].Com = Peaks[pk].Com;
		Peaks_backup[pk].Polarity = Peaks[pk].Polarity;
		Peaks_backup[pk].Slope = Peaks[pk].Slope;
		Peaks_backup[pk].Maxpos = Peaks[pk].Maxpos;
		Peaks_backup[pk].Maximum = Peaks[pk].Maximum;
		Peaks_backup[pk].Height = Peaks[pk].Height;
		Peaks_backup[pk].HeightfromFit = Peaks[pk].HeightfromFit;
		Peaks_backup[pk].Fwhm = Peaks[pk].Fwhm;
		Peaks_backup[pk].Width = Peaks[pk].Width;
		Peaks_backup[pk].PosHalfLeft = Peaks[pk].PosHalfLeft;
		Peaks_backup[pk].PosHalfRight = Peaks[pk].PosHalfRight;
		Peaks_backup[pk].Startpos = Peaks[pk].Startpos;
		Peaks_backup[pk].Integral = Peaks[pk].Integral;
		Peaks_backup[pk].Stoppos = Peaks[pk].Stoppos;
		Peaks_backup[pk].scalefactor = Peaks[pk].scalefactor;
		Peaks_backup[pk].calcbinoffset = Peaks[pk].calcbinoffset;
		Peaks_backup[pk].runtime = Peaks[pk].runtime;
		Peaks_backup[pk].is_valid = Peaks[pk].is_valid;
		Peaks_backup[pk].is_TDC_signal = Peaks[pk].is_TDC_signal;
		Peaks_backup[pk].dpa_is_potential_double = Peaks[pk].dpa_is_potential_double;
		Peaks_backup[pk].dpa_max1_pos = Peaks[pk].dpa_max1_pos;
		Peaks_backup[pk].dpa_max1_hight = Peaks[pk].dpa_max1_hight;
		Peaks_backup[pk].dpa_min_pos = Peaks[pk].dpa_min_pos;
		Peaks_backup[pk].dpa_min_hight = Peaks[pk].dpa_min_hight;
		Peaks_backup[pk].dpa_max2_pos = Peaks[pk].dpa_max2_pos;
		Peaks_backup[pk].dpa_max2_hight = Peaks[pk].dpa_max2_hight;
		Peaks_backup[pk].dpa_fwhm1_pos = Peaks[pk].dpa_fwhm1_pos;
		Peaks_backup[pk].dpa_fwhm1_hight = Peaks[pk].dpa_fwhm1_hight;
		Peaks_backup[pk].dpa_fwhm1_slope = Peaks[pk].dpa_fwhm1_slope;
		Peaks_backup[pk].dpa_fwhm2_pos = Peaks[pk].dpa_fwhm2_pos;
		Peaks_backup[pk].dpa_fwhm2_hight = Peaks[pk].dpa_fwhm2_hight;
		Peaks_backup[pk].dpa_fwhm2_slope = Peaks[pk].dpa_fwhm2_slope;
		Peaks_backup[pk].left_curved_cnt = Peaks[pk].left_curved_cnt;
	}

}









void Event::prepair_pulse(Pulse *curPulse, int NbrPointsFFT) {
	//int NbrPointsFFT = getNbrPointsFFT(curPulse->DataLength);

	// add zeroes
	for (int i = curPulse->DataLength - 1; i < NbrPointsFFT; i++) {
		curPulse->Waveform[i] = 0.;
	}
	curPulse->DataLength = NbrPointsFFT;


	// find center of mass
	double integral = 0;
	double wichtung = 0;
	for (int i = 0; i < curPulse->DataLength; ++i)
	{
		integral += curPulse->Waveform[i];	//calc integral
		wichtung += curPulse->Waveform[i] * i;		//calc weight
	}
	double COM = wichtung / integral;

	// shift data
	double *copy = new double[curPulse->DataLength];
	for (int i = 0; i < curPulse->DataLength; i++) {
		copy[i] = 0;
	}
	int shift = (int)(curPulse->DataLength / 2 - COM);
	if (shift >= 0) {
		for (int i = 0; i < curPulse->DataLength; i++) {
			copy[(i + shift) % (curPulse->DataLength)] = curPulse->Waveform[i];
		}
	}
	else if(shift < 0) {
		for (int i = 0; i < curPulse->DataLength; i++) {
			copy[i] = curPulse->Waveform[(i - shift) % (curPulse->DataLength-1)];
		}
	}
	// temp copy data back
	for (int i = 0; i < curPulse->DataLength; i++) {
		curPulse->Waveform[i] = copy[i];
	}

	delete[] copy;
}
void Event::shift_pulse(Pulse *curPulse) {
	double tmp = 0.;
	for (int i = 0; i < curPulse->DataLength / 2; i++) {
		tmp = curPulse->Waveform[i];
		
		curPulse->Waveform[i] = curPulse->Waveform[i + curPulse->DataLength / 2];
		curPulse->Waveform[i + curPulse->DataLength / 2] = tmp;
	}
}
void Event::compress_pulse(Pulse *curPulse, double compress_rate) {
	double *compressed = new double[curPulse->DataLength];

	for (int i = 0; i < curPulse->DataLength; i++) {
		compressed[i] = curPulse->GetWaveformAt(i*(1 / compress_rate));
	}

	for (int i = 0; i < curPulse->DataLength; i++) {
		if (compressed[i] >= 0.)
			curPulse->Waveform[i] = compressed[i];
		else
			curPulse->Waveform[i] = 0.;
	}

	delete[] compressed;
}
//
//void Event::showPulseFFT(Pulse *curPulse, Pulse *meanPulse) {
//	EventCanvas = Canv("EventCanvas", 2, 6);
//	EventCanvas->SetWindowSize(600, 1500);
//
//	const double *data;
//	const double *dataCFD;
//	TH1D* ChHisto;
//	TH1D* CFDHisto;
//	TH1D* SecHisto;
//
//	TArrow *arrow = 0;
//
//	long chlength = 300; //tbb unknown
//	int nbrpeaks = 0;
//
//	///////////////////////////////////
//	//Initialize Histos////////////////
//	///////////////////////////////////
//	char tmp[256];
//	char name[64];
//	sprintf(tmp, "channel%i", curPulse->Pu_from_Ch);
//	sprintf(name, "Channel %i", curPulse->Pu_from_Ch);
//	ChHisto = new TH1D(tmp, name, 1, 0, 1);
//	CFDHisto = new TH1D(tmp, name, 1, 0, 1);
//	SecHisto = new TH1D(tmp, name, 1, 0, 1);
//
//	//labeling the axis//
//	ChHisto->SetXTitle(tmp);
//	ChHisto->GetXaxis()->CenterTitle();
//	ChHisto->SetYTitle("U [mV]");
//	ChHisto->GetYaxis()->CenterTitle();
//	ChHisto->SetBit(TH1::kNoStats);
//
//	//--set the histogram to the right size--//
//	ChHisto->Reset();
//	ChHisto->SetBins(chlength, 0 - 0.5, chlength - 0.5);
//	ChHisto->SetMaximum(adc_set[curPulse->Pu_from_Ch].fullscale + 100);
//	ChHisto->SetMinimum(-100);
//
//	CFDHisto->Reset();
//	CFDHisto->SetBins(chlength, 0 - 0.5, chlength - 0.5);
//	CFDHisto->SetMaximum(adc_set[curPulse->Pu_from_Ch].fullscale + 100);
//	CFDHisto->SetMinimum(-100);
//
//	SecHisto->Reset();
//	SecHisto->SetBins(chlength, 0 - 0.5, chlength - 0.5);
//	SecHisto->SetMaximum(adc_set[curPulse->Pu_from_Ch].fullscale + 100);
//	SecHisto->SetMinimum(-100);
//
//	//////////////////////////////
//	//Fill the Histos/////////////
//	//////////////////////////////
//	for (size_t j = 0; j<chlength; ++j)
//	{
//		ChHisto->SetBinContent(j + 1, 0);
//	}
//
//	double *x1 = 0, *x2 = 0, *y1 = 0, *y2 = 0;
//
//	data = &curPulse->Waveform[0];
//
//	// CFD
//	const long delay = (long)adc_set[curPulse->Pu_from_Ch].delay * 1000 / adc_set[curPulse->Pu_from_Ch].samplerate;
//	int bin = 0;
//	int cfd_k = 0;
//	double height = 0, heightCFD = 0;
//
//	for (int k = 0; k < curPulse->DataLength; ++k)
//	{
//		bin = int(curPulse->timestamp / adc_set[curPulse->Pu_from_Ch].samplerate) + k;
//
//		height = data[k] / adc_set[curPulse->Pu_from_Ch].gain;
//
//		cfd_k = k -/*+*/ delay;
//		if (cfd_k >= 0 && cfd_k < curPulse->DataLength)
//			heightCFD = (1.)*data[cfd_k] / adc_set[curPulse->Pu_from_Ch].gain;
//		else
//			heightCFD = 0.;
//
//		ChHisto->SetBinContent(bin, height);
//		CFDHisto->SetBinContent(bin, height - heightCFD);
//	}
//
//	// second derivitiv
//	for (int k = 1; k < curPulse->DataLength; ++k) {
//		bin = int(curPulse->timestamp / adc_set[curPulse->Pu_from_Ch].samplerate) + k;
//
//		SecHisto->SetBinContent(bin, (data[k - 1] - 2 * data[k] + data[k + 1]) / adc_set[curPulse->Pu_from_Ch].gain);
//	}
//
//
//	// draw arrows
//	nbrpeaks = curPulse->NbrPeaks;
//	if (curPulse->NbrPeaks > 0) {
//		arrow = new TArrow[curPulse->NbrPeaks];
//		x1 = new double[curPulse->NbrPeaks];
//		x2 = new double[curPulse->NbrPeaks];
//		y1 = new double[curPulse->NbrPeaks];
//		y2 = new double[curPulse->NbrPeaks];
//
//		for (int peak = 0; peak < curPulse->NbrPeaks; peak++) {
//			arrow[peak].SetAngle(40);
//			arrow[peak].SetFillColor(kRed);
//			arrow[peak].SetLineColor(kRed);
//
//			x1[peak] = curPulse->Peaks[peak].Time_ss + curPulse->timestamp / adc_set[curPulse->Pu_from_Ch].samplerate;
//			x2[peak] = x1[peak];
//			y1[peak] = 0;
//			y2[peak] = 1000;
//		}
//	}
//
//	EventCanvas->cd(1);		//change to the right pad
//	ChHisto->Draw("");
//	ChHisto->SetLineColor(2);
//	CFDHisto->Draw("same");
//	SecHisto->SetLineColor(3);
//	SecHisto->Draw("same");
//
//	for (int peak = 0; peak < nbrpeaks; peak++)
//	{
//		if (arrow) arrow[peak].DrawArrow(x1[peak], y1[peak], x2[peak], y2[peak], 0.01, "<|");
//		std::cout << "\n" << curPulse->Peaks[peak].left_curved_cnt;
//	}
//
//	if (x1) delete[] x1;
//	if (x2) delete[] x2;
//	if (y1) delete[] y1;
//	if (y2) delete[] y2;
//	nbrpeaks = 0;
//
//
//
//	
//
//
//	
//	// ######################################################################################################################
//	// FFT stuff
//	FFT_class *fft = new FFT_class;
//
//	if (curPulse->DataLength > 128)
//		return;
//
//	// original pulse..
//	EventCanvas->cd(3);
//	TH1D *original = new TH1D("before", "before", curPulse->DataLength, 0, curPulse->DataLength);
//	for (int i = 0; i < curPulse->DataLength; i++) {
//		original->SetBinContent(i + 1, curPulse->Waveform[i]);
//	}
//	original->Draw("");
//
//
//	//ofstream data_txt;
//	//static int running_id = 0;
//	//char data_txt_name[100];
//	//sprintf(data_txt_name, "%d.txt", running_id);
//	//running_id++;
//
//	//data_txt.open(data_txt_name, std::ios::out);
//	//for (int i = 0; i < meanPulse->DataLength; i++) {
//	//	data_txt << meanPulse->Waveform[i] << " ";
//	//}
//	//data_txt << std::endl;
//	//for (int i = 0; i < curPulse->DataLength; i++) {
//	//	data_txt << curPulse->Waveform[i] << " ";
//	//}
//	//data_txt << std::endl;
//	//data_txt << curPulse->mixed_dTime_ss << " " << curPulse->mixed_Time1_ss << " " << curPulse->mixed_Time2_ss;
//
//	//data_txt.close();
//
//
//
//
//
//
//	// prepaired pulse..
//	EventCanvas->cd(4);
//	prepair_pulse(curPulse, fft->getNbrPointsFFT(curPulse->DataLength));
//	TH1D *prepaired = new TH1D("before", "before", curPulse->DataLength, 0, curPulse->DataLength);
//	for (int i = 0; i < curPulse->DataLength; i++) {
//		prepaired->SetBinContent(i + 1, curPulse->Waveform[i]);
//	}
//	prepaired->Draw("");  
//
//	TH1D *shifted = new TH1D("shifted", "shifted", curPulse->DataLength, 0, curPulse->DataLength);
//	shift_pulse(curPulse);
//	for (int i = 0; i < curPulse->DataLength; i++) {
//		shifted->SetBinContent(i + 1, curPulse->Waveform[i]);
//	}
//	shifted->Draw("same");
//	shifted->SetLineColor(2);
//
//	// FFT
//	EventCanvas->cd(5);
//	double *data_fft = 0;
//	fft->fft(curPulse->Waveform, curPulse->DataLength, data_fft);
//	double *data_Afft = 0;
//	fft->getAmplitude(data_fft, curPulse->DataLength, data_Afft);
//	int nbrPointsFFT = fft->getNbrPointsFFT(curPulse->DataLength);
//	TH1D *Afft = new TH1D("FFT", "FFT", nbrPointsFFT, 0, nbrPointsFFT);
//	for (int i = 0; i < nbrPointsFFT; i++) {
//		Afft->SetBinContent(i + 1, data_Afft[i]);
//	}
//
//	double *data_Afft_shifted = 0;
//	fft->fft_shift(data_fft, curPulse->DataLength);
//	fft->getAmplitude(data_fft, curPulse->DataLength, data_Afft_shifted);
//	TH1D *Afft_shifted = new TH1D("FFT_shifted", "FFT_shifted", nbrPointsFFT, 0, nbrPointsFFT);
//	for (int i = 0; i < nbrPointsFFT; i++) {
//		Afft_shifted->SetBinContent(i + 1, data_Afft_shifted[i]);
//	}
//
//	Afft->Draw("");
//	Afft_shifted->Draw("same");
//	Afft_shifted->SetLineColor(2);
//
//	// MeanPulse
//	EventCanvas->cd(7);
//	Pulse *myMean = new Pulse();
//	myMean->DataLength = meanPulse->DataLength;
//	for (int i = 0; i < myMean->DataLength; i++) {
//		myMean->Waveform[i] = meanPulse->Waveform[i];
//	}
//	for (int i = myMean->DataLength; i < MAX_DATALENGHT; i++) {
//		myMean->Waveform[i] = 0;
//	}
//	myMean->DataLength = curPulse->DataLength;
//	TH1D *myMean_hist = new TH1D("myMean", "myMean", myMean->DataLength, 0, myMean->DataLength);
//	for (int i = 0; i < myMean->DataLength; i++) {
//		myMean_hist->SetBinContent(i, myMean->Waveform[i]);
//	}
//
//	// Second Mean Pulse for testing
//	Pulse *secMean = new Pulse();
//	secMean->DataLength = meanPulse->DataLength;
//	for (int i = 0; i < secMean->DataLength; i++) {
//		secMean->Waveform[i] = meanPulse->Waveform[i];
//	}
//	for (int i = secMean->DataLength; i < MAX_DATALENGHT; i++) {
//		secMean->Waveform[i] = 0;
//	}
//	secMean->DataLength = curPulse->DataLength;
//	compress_pulse(secMean, 0.78);
//
//
//	
//	myMean_hist->Draw("");
//	// compress pulse
//	compress_pulse(myMean, 0.8);
//	TH1D *myMean_compressed = new TH1D("myMean", "myMean", myMean->DataLength, 0, myMean->DataLength);
//	for (int i = 0; i < myMean->DataLength; i++) {
//		myMean_compressed->SetBinContent(i, myMean->Waveform[i]);
//	}
//	myMean_compressed->Draw("same");
//	myMean_compressed->SetLineColor(2);
//
//
//	// prepair
//	EventCanvas->cd(8);
//	prepair_pulse(myMean, fft->getNbrPointsFFT(myMean->DataLength));
//	TH1D *mean_prepaired = new TH1D("mean_prepaired", "mean_prepaired", myMean->DataLength, 0, myMean->DataLength);
//	for (int i = 0; i < myMean->DataLength; i++) {
//		mean_prepaired->SetBinContent(i + 1, myMean->Waveform[i]);
//	}
//	mean_prepaired->Draw("");
//
//	TH1D *mean_shifted = new TH1D("mean_shifted", "mean_shifted", curPulse->DataLength, 0, curPulse->DataLength);
//	shift_pulse(myMean);
//	for (int i = 0; i < myMean->DataLength; i++) {
//		mean_shifted->SetBinContent(i + 1, myMean->Waveform[i]);
//	}
//	mean_shifted->Draw("same");
//	mean_shifted->SetLineColor(2);
//
//	// fft
//	EventCanvas->cd(9);
//	double *mean_fft = 0;
//	fft->fft(myMean->Waveform, myMean->DataLength, mean_fft);
//	double *mean_Afft = 0;
//	fft->getAmplitude(mean_fft, curPulse->DataLength, mean_Afft);
//	int mean_nbrPointsFFT = fft->getNbrPointsFFT(myMean->DataLength);
//	TH1D *mean_Afft_hist = new TH1D("meanFFT", "meanFFT", mean_nbrPointsFFT, 0, mean_nbrPointsFFT);
//	for (int i = 0; i < mean_nbrPointsFFT; i++) {
//		mean_Afft_hist->SetBinContent(i + 1, mean_Afft[i]);
//	}
//	double *mean_Afft_shifted = 0;
//	fft->fft_shift(mean_fft, myMean->DataLength);
//	fft->getAmplitude(mean_fft, curPulse->DataLength, mean_Afft_shifted);
//	TH1D *mean_Afft_shifted_hist = new TH1D("FFT_shifted", "FFT_shifted", mean_nbrPointsFFT, 0, mean_nbrPointsFFT);
//	for (int i = 0; i < mean_nbrPointsFFT; i++) {
//		mean_Afft_shifted_hist->SetBinContent(i + 1, mean_Afft_shifted[i]);
//	}
//		
//	mean_Afft_hist->Draw("");
//	mean_Afft_shifted_hist->Draw("same");
//	mean_Afft_shifted_hist->SetLineColor(2);
//
//
//	//EventCanvas->cd(11);
//	//double *data_recover = 0;
//	//fft->fft_shift(data_fft, curPulse->DataLength);
//	//fft->ifft(data_fft, curPulse->DataLength, data_recover);
//	//TH1D *recover_hist = new TH1D("recover", "recover", curPulse->DataLength, 0, curPulse->DataLength);
//	//for (int i = 0; i < curPulse->DataLength; i++) {
//	//	recover_hist->SetBinContent(i + 1, data_recover[i]);
//	//}
//	//recover_hist->Draw();
//
//	// division
//	EventCanvas->cd(12);
//	fft->fft_shift(data_fft, curPulse->DataLength);
//	fft->fft_divide(data_fft, mean_fft);
//	double *result = 0;
//	double *r_real = 0;
//	double *r_imag = 0;
//	TH1D *result_hist = new TH1D("iResult", "iResult", curPulse->DataLength, 0, curPulse->DataLength);
//	TH1D *result_real = new TH1D("iResult", "iResult", curPulse->DataLength, 0, curPulse->DataLength);
//	TH1D *result_imag = new TH1D("iResult", "iResult", curPulse->DataLength, 0, curPulse->DataLength);
//	fft->ifft(data_fft, curPulse->DataLength, result, r_real, r_imag);
//	for (int i = 0; i < curPulse->DataLength; i++) {
//		result_hist->SetBinContent(i + 1, result[i]);
//		result_real->SetBinContent(i + 1, r_real[i]);
//		result_imag->SetBinContent(i + 1, r_imag[i]);
//	}
//	result_hist->Draw();
//	result_real->Draw("same");
//	result_imag->Draw("same");
//
//
//	// secMean Test
//	//double *secMean_fft = 0;
//	//fft->fft(secMean->Waveform, secMean->DataLength, secMean_fft);
//	//fft->fft_shift(secMean_fft, secMean->DataLength);
//	//fft->fft_divide(mean_fft, secMean_fft);
//	//double *result_mean = 0;
//	//TH1D *result_mean_hist = new TH1D("meanResult", "meanResult", curPulse->DataLength, 0, curPulse->DataLength);
//	//fft->ifft(mean_fft, curPulse->DataLength, result_mean/*, 1000*/);
//	//EventCanvas->cd(2);
//	//for (int i = 0; i < curPulse->DataLength; i++) {
//	//	result_mean_hist->SetBinContent(i + 1, result_mean[i]);
//	//}
//	//result_mean_hist->Draw();
//
//
//
//
//
//
//
//
//
//	gPad->Update();
//
//	while (!_kbhit())
//	{
//		gSystem->Sleep(50);
//		gSystem->ProcessEvents();
//	}
//
//	char ch;
//	while (_kbhit()) ch = _getch();
//
//	delete[] arrow;
//	delete ChHisto;
//	delete CFDHisto;
//	delete SecHisto;
//
//	//delete result_hist;
//}
//





//Helper Functions
void createNewtonPolynomial(const double * x, const double * y, double * coeff)
{
	//**this function creates the coefficients for Newton interpolating Polynomials	**//
	//**Newton Polynomial are Created from n Points and have the form				**//
	//**p(x) = c0 + c1(x-x0) + c2(x-x0)(x-x1)+...+c_(n-1)(x-x0)(x-x1)...(x-x_(n-2))	**//
	//**given that you have n Points (x0,y0), (x1,y1), ..., (x_(n-1),y_(n-1))		**//

	double f_x0_x1 = (y[1]-y[0]) / (x[1]-x[0]);
	double f_x1_x2 = (y[2]-y[1]) / (x[2]-x[1]);
	double f_x2_x3 = (y[3]-y[2]) / (x[3]-x[2]);

	double f_x0_x1_x2 = (f_x1_x2 - f_x0_x1) / (x[2]-x[0]);
	double f_x1_x2_x3 = (f_x2_x3 - f_x1_x2) / (x[3]-x[1]);

	double f_x0_x1_x2_x3 = (f_x1_x2_x3 - f_x0_x1_x2) / (x[3]-x[0]);

	coeff[0] = y[0];
	coeff[1] = f_x0_x1;
	coeff[2] = f_x0_x1_x2;
	coeff[3] = f_x0_x1_x2_x3;
}
double evalNewtonPolynomial(const double * x, const double * coeff, double X)
{
	//**this function evaluates the Newton Polynomial that was created from n Points**//
	//** (x0,y0),..., (x(n-1),y(n-1)) with coefficients (c0,...,c(n-1))				**//
	//**using Horner's Rule															**//

	double returnValue = coeff[3];
	returnValue = returnValue * (X - x[2]) + coeff[2];
	returnValue = returnValue * (X - x[1]) + coeff[1];
	returnValue = returnValue * (X - x[0]) + coeff[0];

	return returnValue;
}
double findXForGivenY(const double * x, const double * coeff, const double Y, const double Start)
{
	//initialisiere die Grenzen//
	MyPunkt Low(x[1], evalNewtonPolynomial(x,coeff,x[1]));
	MyPunkt Up (x[2], evalNewtonPolynomial(x,coeff,x[2]));

	//initialisiere den iterierenden Punkt mit dem Startwert//
	MyPunkt p (Start, evalNewtonPolynomial(x,coeff,Start));

	//ist der Startpunkt schon der richtige Punkt//
	//liefere den dazugehörigen x-Wert zurück//
	if (p.y() == Y)
		return p.x();

	//finde heraus ob es ein positiver oder ein negativer Durchgang ist//
	bool Neg = (Low.y() > Up.y())?true:false;

	//der Startpunkt soll die richtige neue Grenze bilden//
	if (Neg)	//wenn es ein negativer Druchgang ist
	{
		if (p.y() > Y)		//ist der y-Wert grösser als der gewollte
			Low = p;		//bildet der Punkt die neue untere Grenze
		else if (p.y() < Y)	//ist der y-Wert ist kleiner als der gewollte
			Up = p;			//bildet der Punkt die neue obere Grenze
		else				//ist der Punkt genau getroffen
			return p.x();	//liefer den dazugehörigen x-Wert zurück
	}
	else		//wenn es ein positiver Druchgang ist
	{
		if (p.y() > Y)		//und der y-Wert grösser als der gewollte
			Up = p;			//bildet der Punkt die neue obere Grenze
		else if (p.y() < Y)	//und y-Wert ist kleiner als der gewollte
			Low = p;		//bildet der Punkt die neue untere Grenze
		else				//ist der Punkt genau getroffen
			return p.x();	//liefer den dazugehörigen x-Wert zurück
	}


	while((Up.x()-Low.x()) > 0.005) //iteriere solange bis der Abstand zwischen den x-Werten
		//kleiner als 0.005
	{
		//bilde das arithmetische Mittel zwischen beiden Grenzen//
		//das ist der neue x-Wert unseres Punktes//
		p.x() = 0.5 * (Up.x()+Low.x());
		//finde den dazugehörigen y-Wert//
		p.y() = evalNewtonPolynomial(x,coeff,p.x());

		if (Neg)	//wenn es ein negativer Druchgang ist
		{
			if (p.y() > Y)		//und der y-Wert grösser als der gewollte
				Low = p;		//bildet der Punkt die neue untere Grenze
			else if (p.y() < Y)	//und der y-Wert ist kleiner als der gewollte
				Up = p;			//bildet der Punkt die neue obere Grenze
			else				//ist der Punkt genau getroffen
				return p.x();		//liefer den dazugehörigen x-Wert zurück
		}
		else		//wenn es ein positiver Druchgang ist
		{
			if (p.y() > Y)		//und der y-Wert grösser als der gewollte
				Up = p;			//bildet der Punkt die neue obere Grenze
			else if (p.y() < Y)	//und y-Wert ist kleiner als der gewollte
				Low = p;		//bildet der Punkt die neue untere Grenze
			else				//ist der Punkt genau getroffen
				return p.x();		//liefer den dazugehörigen x-Wert zurück
		}
		//std::cout<<"("<<Low.x<<","<<Low.y<<")   ("<<p.x<<","<<p.y<<")   ("<<Up.x<<","<<Up.y<<") "<<Y<<std::endl;
	}
	return ((Up.x() + Low.x())*0.5);
}
void linearRegression(const int nbrPoints, const double x[], const double y[], double &m, double &c)
{
	//--this funktion does a linear regression of 4 points--//
	//--getting a line that follows the form: y(x) = m*x + c--//
	//--return the x value for a given y(x) (yb)--//
	//-- => x = (y(x)-c)/m--//
	double SumXsq=0.,SumX=0.,SumY=0.,SumXY=0.;
	for (int i=0;i<nbrPoints;++i)
	{
		SumX	+=  x[i];
		SumY	+=  y[i];
		SumXY	+= (x[i]*y[i]);
		SumXsq	+= (x[i]*x[i]);
	}

	double a1 = ((SumX*SumX) - (nbrPoints*SumXsq));

	m = ((SumX*SumY) - (nbrPoints*SumXY)) / a1;
	c = ((SumX*SumXY) - (SumY*SumXsq)) / a1;
}
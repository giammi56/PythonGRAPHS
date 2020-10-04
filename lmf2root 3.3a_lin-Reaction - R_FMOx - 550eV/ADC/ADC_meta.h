#pragma once

#ifndef __ADC_meta_h_
#define __ADC_meta_h_

#include "ADC_analysis.h"
#include "../Src/Ueberstruct.h"
#include "Event.h"

class ADC_meta {
public:
	Ueberstruct			*Ueber_ptr;

	// singelton instance
	static ADC_meta		*instance(Ueberstruct *Ueber_ptr, int max_nbr_channels, int max_nbr_pulses, int max_nbr_peaks);
	static void			deleteInstance() { delete inst; inst = 0;}
	
	ADC_analysis		**adc;

	bool				GetTDCandCNTArray(double *tdc_ns, __int32 cnt[]);	//lmf2root interface

	Event				*m_event;								// working event

protected:
	ADC_meta(Ueberstruct *Ueber_ptr, int max_nbr_channels, int max_nbr_pulses, int max_nbr_peaks);
	~ADC_meta();	

	void				plot_adc_infos();

	void				InitializeADCanalysis();
	void				InitializeDPAanalysis();
	
	bool				ReadNextEvent(Event *newEvent);
	bool				ReadNextfADC4event(Event *newEvent);	// read next fADC4 event (todo: compatible for fadc4 and fadc8 data)
	void				ReadNextfADC8event(Event *newEvent);
	void				ReadNextAGATevent(Event *newEvent);	// read next AGAT event 
	bool				ReadNextHDF5event(Event *newEvent);

	void				analyseEvent(Event *ev);

private:
	static	ADC_meta	*inst;

	lma2root			*lma;							// extention for AGAT files

	__int32				event_counter;						// counts every read event (ReadNextfADC4event)
	
	void				save_lmf_startpos();				// save the current lmf file beginnning
	void				goto_lmf_startpos();				// goto saved lmf startpos
	__int64				lmf_startpos;						// position of the lmf file beginning (initalized in constructor)

	int					max_nbr_channels;					// max number of channels 
	int					max_nbr_pulses;						// max number of pulses per channel
	int					max_nbr_peaks;						// max number of peaks per pulse

	TH1D*				conv_hist(H1D* hist);
	TH2D*				conv_hist(H2D* hist);

	bool				showSinglePulse(Pulse *pu);
	bool				showAllPulses(Event *ev);
	void				WriteRawPulseInfo(int nbr_events);
	void				WriteOscilliscopeMode(int nbr_events, int trigger_ch);
	void				WriteSignalInspector(int nbr_events);

	int					*used_channels_list;	// 0 not used - 1 adc used - 2 dpa used
	int					used_channels_nbr;

	void				InspectADC();
};

#endif
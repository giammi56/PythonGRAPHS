#include "ADC_meta.h"



TDirectory* getDir(TFile *rootfile, TString dirName);
TCanvas* Canv(const char *name, int hor, int ver, int hor_size = 500, int ver_size = 200);

//singelton
ADC_meta *ADC_meta::instance(Ueberstruct *Ueber_ptr, int max_nbr_channels, int max_nbr_pulses, int max_nbr_peaks) {
	if (!inst)
		inst = new ADC_meta(Ueber_ptr, max_nbr_channels, max_nbr_pulses, max_nbr_peaks);
	return inst;
}

// public
ADC_meta::ADC_meta(Ueberstruct *Ueber_ptr, int max_nbr_channels, int max_nbr_pulses, int max_nbr_peaks) {
	this->Ueber_ptr = Ueber_ptr;

	this->max_nbr_channels	= max_nbr_channels;
	this->max_nbr_pulses	= max_nbr_pulses;
	this->max_nbr_peaks		= max_nbr_peaks;

	m_event = new Event();

	lma = nullptr;											// AGAT input file.. no multiple files supported jet
	if (Ueber_ptr->LMF_input) {
	if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_AGAT)
		lma = new lma2root(Ueber_ptr->input_data_file_name);
	}

	adc = new ADC_analysis*[max_nbr_channels];
	for (int ch = 0; ch < max_nbr_channels; ch++) adc[ch] = nullptr;

	event_counter = 0;
	lmf_startpos = 0;
	if (Ueber_ptr->LMF_input) save_lmf_startpos();

	used_channels_list = new int[max_nbr_channels];
	memset(used_channels_list, 0, max_nbr_channels*sizeof(int));
	used_channels_nbr = 0;
	InitializeADCanalysis();								// initialize ADC analysis (channel settings, used channels..)

	if (1 == Ueber_ptr->parameter[810]) {
		InspectADC();
	}

	if (Ueber_ptr->parameter[900]) {						// use dpa analysis (config)	
		InitializeDPAanalysis();							// initialize DPA analysis
		plot_adc_infos();
	}

	if (Ueber_ptr->LMF_input) goto_lmf_startpos();
	if (Ueber_ptr->HDF5_input) Ueber_ptr->HDF5_input->seek_to_begin();
}

ADC_meta::~ADC_meta() {
	delete m_event;
	m_event = 0;

	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if(adc[ch])
			delete adc[ch];
		adc[ch] = nullptr;
	}

}

bool ADC_meta::GetTDCandCNTArray(double * tdc_ns, __int32 cnt[]) {
	__int32 max_channel = NUM_CHANNELS;
	__int32 max_hits = NUM_IONS;

	if (!ReadNextEvent(m_event))
		return false;

	analyseEvent(m_event);

	for (int ch = 0; ch < NUM_CHANNELS; ch++) {
		cnt[ch] = 0;
		if (!m_event->channels[ch])continue;
		for (int pu = 0; pu < m_event->channels[ch]->NbrPulses; pu++) {
			for (int pk = 0; pk < m_event->channels[ch]->Pulses[pu].NbrPeaks; pk++) {
				tdc_ns[ch*NUM_CHANNELS + cnt[ch]] = (m_event->channels[ch]->Pulses[pu].timestamp + m_event->channels[ch]->Pulses[pu].Peaks[pk].Time_ss * adc[ch]->m_adc_settings->get_samplerate())*0.001;
				cnt[ch]++;
			}
		}
	}
	return true;
}

void ADC_meta::analyseEvent(Event *ev) {
	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if (ev->channels[ch] && ev->channels[ch]->NbrPulses > 0 && adc[ch]) {
			for (int pu = 0; pu < ev->channels[ch]->NbrPulses; pu++) {
				adc[ch]->analyze_Pulse(&ev->channels[ch]->Pulses[pu]);
			}
		}
	}
}

void ADC_meta::InspectADC() {
	this->save_lmf_startpos();

	if (Ueber_ptr->parameter[811]) {				// show all pulses
		while (true) {
			ReadNextEvent(m_event);
			
			for (int ch = 0; ch < max_nbr_channels; ch++) {
				if (m_event->channels[ch] && m_event->channels[ch]->NbrPulses)
					adc[ch]->analyze_Pulse(&m_event->channels[ch]->Pulses[0]);
			}

			if (0 == this->showAllPulses(m_event))
				break;
		}
	}
	this->goto_lmf_startpos();

	if (Ueber_ptr->parameter[812] >= 0) {			// show pulse in channel
		int ch = int(Ueber_ptr->parameter[812]);
		while (true) {
			if (!used_channels_list[ch])
				break;

			ReadNextEvent(m_event);
			//ReadNextfADC4event(m_event);
			
			if (m_event->channels[ch] && 0 < m_event->channels[ch]->NbrPulses) {
				if (!m_event->channels[ch]->Pulses[0].is_TDC_signal) {
				this->adc[ch]->analyze_Pulse(&m_event->channels[ch]->Pulses[0]);
				if (m_event->channels[ch]->Pulses[0].NbrPeaks == 2) {
					if (0 == this->showSinglePulse(&m_event->channels[ch]->Pulses[0]))
						break;
				}
			}
		}
	}
	}
	this->goto_lmf_startpos();
	
	if (Ueber_ptr->parameter[820] && (Ueber_ptr->parameter[821] || Ueber_ptr->parameter[822] || Ueber_ptr->parameter[824])) {
		if (1 == Ueber_ptr->parameter[821]) {
			WriteRawPulseInfo(int(Ueber_ptr->parameter[820]));
			this->goto_lmf_startpos();
		}

		if (1 == Ueber_ptr->parameter[822]) {
			WriteOscilliscopeMode(int(Ueber_ptr->parameter[820]), int(Ueber_ptr->parameter[823]));
			this->goto_lmf_startpos();
		}

		if (1 == Ueber_ptr->parameter[824]) {
			WriteSignalInspector(int(Ueber_ptr->parameter[820]));
			this->goto_lmf_startpos();
		}
	}
}

bool ADC_meta::ReadNextEvent(Event *newEvent) {
	if (Ueber_ptr->LMF_input) {
		if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_FADC4) {
			if (!ReadNextfADC4event(newEvent))
				return false;
			if (Ueber_ptr->LMF_input->errorflag) return false;
			return true;
		}
		if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_FADC8) {
			ReadNextfADC8event(newEvent);
			if (Ueber_ptr->LMF_input->errorflag) return false;
			return true;
		}
		if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_AGAT) {
			ReadNextAGATevent(newEvent);
			if (Ueber_ptr->LMF_input->errorflag) return false;
			return true;
		}
	}
	if (Ueber_ptr->input_is_HDF5) {
		if (!ReadNextHDF5event(newEvent))
			return false;
		if (Ueber_ptr->HDF5_input->error_flag) return false;
		if (Ueber_ptr->HDF5_input->eof) return false;
		return true;
	}
	return false;	// no ReadNextEvent function found for DAQ_ID
}


bool ADC_meta::ReadNextfADC4event(Event *newEvent) {
	newEvent->resetEvent();

	__int32 * cnt = Ueber_ptr->cnt;
	if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_FADC4)
	{
		bool bEnd_of_group_detected = false;
		__int32 TriggerChannel = Ueber_ptr->LMF_input->fADC4.iTriggerChannel;
		__int16 i16buffer[i16_BUFFERSIZE];

		event_counter++;
		newEvent->id = event_counter;

		while (true) {
			if (Ueber_ptr->LMF_input->ReadNextfADC4packet(&Ueber_ptr->n_packet, bEnd_of_group_detected, i16buffer, i16_BUFFERSIZE)) {
				__int32 cardchannel = Ueber_ptr->n_packet.card * 6 + Ueber_ptr->n_packet.channel;
				if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_16_BIT_SIGNED || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_DATA || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
					if (0 == newEvent->channels[cardchannel]) {
						newEvent->channels[cardchannel] = new Channel;
						newEvent->channels[cardchannel]->Pulses = new Pulse[MAX_NBR_PULSES];
						newEvent->channels[cardchannel]->NbrPulses = 0;
						//adc[cardchannel] = new ADC_analysis(Ueber_ptr->parameter[1000 + cardchannel * 7], Ueber_ptr->parameter[1005 + cardchannel * 7], Ueber_ptr->parameter[1003 + cardchannel * 7], Ueber_ptr->parameter[1004 + cardchannel * 7], Ueber_ptr->parameter[1001 + cardchannel * 7], Ueber_ptr->parameter[1002 + cardchannel * 7], Ueber_ptr->parameter[1006 + cardchannel * 7], Ueber_ptr->LMF_input->fADC4.bits_per_mVolt[0], Ueber_ptr->LMF_input->fADC4.ndigo_parameters[0].sample_period, 500);
					}
				}

				if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_16_BIT_SIGNED) {
					if (newEvent->channels[cardchannel]->NbrPulses < MAX_NBR_PULSES) {
						double offset = Ueber_ptr->LMF_input->fADC4.GND_level[cardchannel];
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Pu_from_Ch = cardchannel;
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength = Ueber_ptr->n_packet.length * 4;		//convert to bins
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].timestamp = Ueber_ptr->n_packet.timestamp;

						for (int i = 0; i < newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength; i++)
						{
							newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Waveform[i] = double(adc[cardchannel]->m_adc_settings->get_polarity() * (i16buffer[i] - int(offset)));
						}

						newEvent->channels[cardchannel]->NbrPulses++;
					}
				}

				if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_DATA || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
					if (newEvent->channels[cardchannel]->NbrPulses < MAX_NBR_PULSES) {
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Pu_from_Ch = cardchannel;
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength = 0;		//convert to bins
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].timestamp = Ueber_ptr->n_packet.timestamp;
						newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].is_TDC_signal = true;
						newEvent->channels[cardchannel]->NbrPulses++;
					}
				}

				if (bEnd_of_group_detected) break;
			}
			else {
				return false;
			}
		}
	}
	return true;
}


bool ADC_meta::ReadNextHDF5event(Event *newEvent) {
	newEvent->resetEvent();

	if (!Ueber_ptr->input_is_HDF5) return false;
	if (!Ueber_ptr->HDF5_input) return false;


	__int32 * cnt = Ueber_ptr->cnt;

	bool bEnd_of_group_detected = false;
	__int32 TriggerChannel = Ueber_ptr->HDF5_input->TriggerChannel;
	if (TriggerChannel < -1) {
		Ueber_ptr->stop_reading_files = true;
		printf("HDF5: No trigger channel given.");
		return false;
	}
	__int16 i16buffer[i16_BUFFERSIZE];

	event_counter++;
	newEvent->id = event_counter;

	while (true) {
		bool packet_read;
		if (Ueber_ptr->HDF5_input->do_grouping) {
			packet_read = Ueber_ptr->HDF5_input->ReadNextHDF5packet_grouping(&Ueber_ptr->n_packet, bEnd_of_group_detected, i16buffer, i16_BUFFERSIZE);
		} else {
			packet_read = Ueber_ptr->HDF5_input->ReadNextHDF5packet_slicing(&Ueber_ptr->n_packet, bEnd_of_group_detected, i16buffer, i16_BUFFERSIZE);
		}
		if (packet_read) {
			if (bEnd_of_group_detected) break;
			__int32 cardchannel = Ueber_ptr->n_packet.card * 4 + Ueber_ptr->n_packet.channel;
			if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_16_BIT_SIGNED || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_DATA || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
				if (0 == newEvent->channels[cardchannel]) {
					newEvent->channels[cardchannel] = new Channel;
					newEvent->channels[cardchannel]->Pulses = new Pulse[MAX_NBR_PULSES];
					newEvent->channels[cardchannel]->NbrPulses = 0;
					//adc[cardchannel] = new ADC_analysis(Ueber_ptr->parameter[1000 + cardchannel * 7], Ueber_ptr->parameter[1005 + cardchannel * 7], Ueber_ptr->parameter[1003 + cardchannel * 7], Ueber_ptr->parameter[1004 + cardchannel * 7], Ueber_ptr->parameter[1001 + cardchannel * 7], Ueber_ptr->parameter[1002 + cardchannel * 7], Ueber_ptr->parameter[1006 + cardchannel * 7], Ueber_ptr->LMF_input->fADC4.bits_per_mVolt[0], Ueber_ptr->LMF_input->fADC4.ndigo_parameters[0].sample_period, 500);
				}
			}

			if (newEvent->channels[cardchannel]->NbrPulses < MAX_NBR_PULSES) {

			if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_16_BIT_SIGNED) {
					double offset = -1500; // GND offset
				// TODO:  use the first samples to measure the GND level?
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Pu_from_Ch = cardchannel;
					newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength = Ueber_ptr->n_packet.length;		//convert to bins
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].timestamp = Ueber_ptr->n_packet.timestamp;

				for (int i = 0; i < newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength; i++)
				{
					newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Waveform[i] = double(adc[cardchannel]->m_adc_settings->get_polarity() * (i16buffer[i] - int(offset)));
				}

				newEvent->channels[cardchannel]->NbrPulses++;
			}

			if (Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_DATA || Ueber_ptr->n_packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].Pu_from_Ch = cardchannel;
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].DataLength = 0;		//convert to bins
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].timestamp = Ueber_ptr->n_packet.timestamp;
				newEvent->channels[cardchannel]->Pulses[newEvent->channels[cardchannel]->NbrPulses].is_TDC_signal = true;
				newEvent->channels[cardchannel]->NbrPulses++;
			}

		}
		} else {
			return false;
		}
	}
	
	return true;
}



void ADC_meta::ReadNextfADC8event(Event *newEvent){}
void ADC_meta::ReadNextAGATevent(Event *newEvent) {
	newEvent->resetEvent();
	__int32 * cnt = Ueber_ptr->cnt;


	long EventIDtemp = 0;
	double Horpostemp = 0.;
	short NbrPulses = 0;

	//--read the event values from the archive--//
	*lma >> (long)EventIDtemp;				//event time stamp
	newEvent->id = EventIDtemp;
	*lma >> (double)Horpostemp;				//horpos from acqiris
											//ev->Horpos = Horpostemp;
	Horpostemp *= 1e9;


	//ReadNumberof Pulses in Channel
	for (int i = 0; i < lma->hi.NbrChannels; i++)
	{
		short tempNbrPulses;
		*lma >> tempNbrPulses;

		if (0 == newEvent->channels[i]) {
			newEvent->channels[i] = new Channel;
			newEvent->channels[i]->Pulses = new Pulse[MAX_NBR_PULSES];
			newEvent->channels[i]->NbrPulses = 0;
		}

		newEvent->channels[i]->NbrPulses = tempNbrPulses;

		for (int j = 0; j < newEvent->channels[i]->NbrPulses; j++)
		{
			newEvent->channels[i]->Pulses[j].Pu_from_Ch = i;

			long IdxFiPOrig;
			*lma >> IdxFiPOrig;
			long DataLength;
			*lma >> DataLength;
			newEvent->channels[i]->Pulses[j].DataLength = DataLength;
			newEvent->channels[i]->Pulses[j].timestamp = __int64((Horpostemp + IdxFiPOrig) * 1000);

			size_t DataSize = newEvent->channels[i]->Pulses[j].DataLength * lma->hi.NbrBytes;
			size_t NewDataSize = newEvent->channels[i]->Pulses[j].DataLength;// * sizeof(float);
			char *data = new char[newEvent->channels[i]->Pulses[j].DataLength];

			//read Data to array
			lma->readArray(data, DataSize);
			//copy Data to array
			for (int h = 0; h < newEvent->channels[i]->Pulses[j].DataLength; h++) {
				newEvent->channels[i]->Pulses[j].Waveform[h] = adc[i]->m_adc_settings->polarity * ((double)data[h] - lma->hi.ChI[i].Baseline);
			}
		}
	}
	//ev->NbrBytes = sizeof(float);

	return;
}



// initialization of double pulse analysis
void ADC_meta::InitializeADCanalysis() {
	int channel_dpa_init_needed[NUM_CHANNELS] = { 0 };

	if (Ueber_ptr->LMF_input) {
		if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_FADC4 || Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_FADC8) {
			for (int ch = 0; ch < max_nbr_channels; ch++) {
				if (!adc[ch]) {
					adc[ch] = new ADC_analysis(Ueber_ptr->parameter[1000 + ch * 7], int(Ueber_ptr->parameter[1005 + ch * 7]), int(Ueber_ptr->parameter[1003 + ch * 7]), Ueber_ptr->parameter[1004 + ch * 7], int(Ueber_ptr->parameter[1001 + ch * 7]), int(Ueber_ptr->parameter[1002 + ch * 7]), Ueber_ptr->parameter[1006 + ch * 7], Ueber_ptr->LMF_input->fADC4.bits_per_mVolt[0], Ueber_ptr->LMF_input->fADC4.ndigo_parameters[0].sample_period
						, 1000);  // full scale in mV
					if (adc[ch]->m_adc_settings->use_dpa != 0 && adc[ch]->m_adc_settings->use_dpa != 1) {
						std::cout << "\nADC Error - ADC_parameters.txt: 'Use DPA for channel' only 0 and 1 allowed!";
						return;
					}
				}
			}
		}
		else if (Ueber_ptr->LMF_input->DAQ_ID == DAQ_ID_AGAT) {
			for (int ch = 0; ch < lma->hi.NbrChannels /*max_nbr_channels*/ ; ch++) {
				if (!adc[ch]) {
					adc[ch] = new ADC_analysis(Ueber_ptr->parameter[1000 + ch * 7], int(Ueber_ptr->parameter[1005 + ch * 7]), int(Ueber_ptr->parameter[1003 + ch * 7]), Ueber_ptr->parameter[1004 + ch * 7], int(Ueber_ptr->parameter[1001 + ch * 7]), int(Ueber_ptr->parameter[1002 + ch * 7]), Ueber_ptr->parameter[1006 + ch * 7], 1 / lma->hi.ChI[ch].Gain, lma->hi.SampInter * 1e9 * 1000, lma->hi.ChI[ch].FullScale);
					if (adc[ch]->m_adc_settings->use_dpa != 0 && adc[ch]->m_adc_settings->use_dpa != 1) {
						std::cout << "\nADC Error - ADC_parameters.txt: 'Use DPA for channel' only 0 and 1 allowed!";
						return;
					}
				}
			}
		}
	}
	
	if (Ueber_ptr->HDF5_input) {
		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (!adc[ch]) {
				adc[ch] = new ADC_analysis(Ueber_ptr->parameter[1000 + ch * 7], int(Ueber_ptr->parameter[1005 + ch * 7]), int(Ueber_ptr->parameter[1003 + ch * 7]), Ueber_ptr->parameter[1004 + ch * 7], int(Ueber_ptr->parameter[1001 + ch * 7]), int(Ueber_ptr->parameter[1002 + ch * 7]), Ueber_ptr->parameter[1006 + ch * 7]
								, Ueber_ptr->HDF5_input->bits_per_mVolt
								, Ueber_ptr->HDF5_input->sample_period_ps
								, Ueber_ptr->HDF5_input->full_scale_mV); // full scale in mV
				if (adc[ch]->m_adc_settings->use_dpa != 0 && adc[ch]->m_adc_settings->use_dpa != 1) {
					std::cout << "\nADC Error - ADC_parameters.txt: 'Use DPA for channel' only 0 and 1 allowed!";
					return;
				}
			}
		}
	}

	__int64 start_pos = 0;
	if (Ueber_ptr->LMF_input) start_pos = Ueber_ptr->LMF_input->input_lmf->tell();

	// find out what channels are active and need dpa init..
	for (int i = 0; i < 10000; i++) {					// ACHIM: Check only first 10000 events to see if channels are used..
		if (!ReadNextEvent(m_event)) {
			event_counter = 0;
			if (Ueber_ptr->LMF_input) Ueber_ptr->LMF_input->input_lmf->seek(start_pos);
			if (Ueber_ptr->HDF5_input) Ueber_ptr->HDF5_input->seek_to_begin();
			return;
		}
		
		for (int ch = 0; ch < NUM_CHANNELS; ch++) {
			if (m_event->channels[ch]) {
				if(adc[ch]->m_adc_settings->use_dpa)
					used_channels_list[ch] = 2;			// used channel with dpa
				else
					used_channels_list[ch] = 1;			// used channel without dpa
			}
		}
	}

	event_counter = 0;
	if (Ueber_ptr->LMF_input) Ueber_ptr->LMF_input->input_lmf->seek(start_pos);
	if (Ueber_ptr->HDF5_input) Ueber_ptr->HDF5_input->seek_to_begin();

	for (int ch = 0; ch < max_nbr_channels; ch++) {		// number of used channels
		if (used_channels_list[ch] > 0)
			used_channels_nbr++;
	}
}



void ADC_meta::InitializeDPAanalysis() {
	int *keep_feeding_array = new int[max_nbr_channels];
	memset(keep_feeding_array, 0, max_nbr_channels*sizeof(int));
	int channels_to_feed = 0;
	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if (2 == used_channels_list[ch]) {
			keep_feeding_array[ch] = 1;
			channels_to_feed++;
		}
	}
	
	std::cout << "\nStart DPA Initialisation for " << channels_to_feed << " channels:";

	int keep_feeding = 0;
	while (event_counter < 10000000) {					// ACHIM: maximum number of events to initialize dpa analysis
		keep_feeding = 0; 
		ReadNextEvent(m_event);
		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (keep_feeding_array[ch]) {
				if (m_event->channels[ch] && 1 == m_event->channels[ch]->NbrPulses) {
					if (0 == adc[ch]->initialize_DPA_analysis(&m_event->channels[ch]->Pulses[0]))
						keep_feeding_array[ch] = 0;
				}
			}
		}
	
		for (int ch = 0; ch < max_nbr_channels; ch++) {	// check if nomore channel needs feeding
			keep_feeding += keep_feeding_array[ch];
		}
		if (!keep_feeding)								
			break;										// feeding fished
	}
}

// private
void ADC_meta::save_lmf_startpos() {
	if (this->Ueber_ptr->LMF_input->DAQ_ID != DAQ_ID_AGAT) {
		lmf_startpos = this->Ueber_ptr->LMF_input->input_lmf->tell();
	}
}
void ADC_meta::goto_lmf_startpos() {
	if (Ueber_ptr->HDF5_input) Ueber_ptr->HDF5_input->seek_to_begin();

	if (this->Ueber_ptr->LMF_input) {
	if (this->Ueber_ptr->LMF_input->DAQ_ID != DAQ_ID_AGAT) {
		event_counter = 0;
		m_event->id = 0;
		if (Ueber_ptr->LMF_input->errorflag == 18) {
			Ueber_ptr->LMF_input->errorflag = 0;
			Ueber_ptr->LMF_input->input_lmf->error = 0;
			Ueber_ptr->LMF_input->input_lmf->eof = false;
		}
		Ueber_ptr->LMF_input->int64_number_of_read_events = 0;
		Ueber_ptr->LMF_input->int64_LMF_EventCounter = 0;
		Ueber_ptr->LMF_input->input_lmf->seek(lmf_startpos);
		}
	}
}

void ADC_meta::plot_adc_infos() {
	TFile *output = new TFile("new6.root", "RECREATE");
	//TH1D *histos[100];

	for (int ch = 0; ch < NUM_CHANNELS; ch++) {
		char dir[100];
		sprintf(dir, "/channel_%d", ch);
		if (adc[ch]) {
			if (adc[ch]->m_adc_settings->fwhm_hist) {
				TH1D *hist = conv_hist(adc[ch]->m_adc_settings->fwhm_hist);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}
			if (adc[ch]->m_adc_settings->hist_ADCpuls1_error) {
				TH2D *hist = conv_hist(adc[ch]->m_adc_settings->hist_ADCpuls1_error);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}

			if (adc[ch]->m_adc_settings->hist_DPApuls1_error) {
				TH2D *hist = conv_hist(adc[ch]->m_adc_settings->hist_DPApuls1_error);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}

			if (adc[ch]->m_adc_settings->hist_ADCpuls2_error) {
				TH2D *hist = conv_hist(adc[ch]->m_adc_settings->hist_ADCpuls2_error);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}

			if (adc[ch]->m_adc_settings->hist_DPApuls2_error) {
				TH2D *hist = conv_hist(adc[ch]->m_adc_settings->hist_DPApuls2_error);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}
			
			if (adc[ch]->m_adc_settings->hist_ADCvsDPA_puls1) {
				TH1D *hist = conv_hist(adc[ch]->m_adc_settings->hist_ADCvsDPA_puls1);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}

			if (adc[ch]->m_adc_settings->hist_ADCvsDPA_puls2) {
				TH1D *hist = conv_hist(adc[ch]->m_adc_settings->hist_ADCvsDPA_puls2);
				hist->SetDirectory(getDir(output, dir));
				hist->Write();
				delete hist; hist = 0;
			}

		}
	}

	output->Close();
	delete output; output = 0;
}
TH1D* ADC_meta::conv_hist(H1D* hist) {
	TH1D* root_hist = new TH1D(hist->GetName(), hist->GetTitle(), hist->GetXaxis()->GetNbins(), hist->GetXaxis()->GetXlow(), hist->GetXaxis()->GetXup());

	for (int x = 0; x < hist->GetNbins() + 2; x++) {
		root_hist->SetBinContent(x, hist->GetBinContent(x));
	}

	return root_hist;
}
TH2D* ADC_meta::conv_hist(H2D* hist) {
	TH2D* root_hist = new TH2D(hist->GetName(), hist->GetTitle(), hist->GetXaxis()->GetNbins(), hist->GetXaxis()->GetXlow(), hist->GetXaxis()->GetXup(), hist->GetYaxis()->GetNbins(), hist->GetYaxis()->GetXlow(), hist->GetYaxis()->GetXup());

	for (int x = 0; x < hist->GetXaxis()->GetNbins() + 2; x++) {
		for (int y = 0; y < hist->GetYaxis()->GetNbins() + 2; y++) {
			root_hist->SetBinContent(x, y, hist->GetBinContent(x, y));
		}
	}

	return root_hist;
}

bool ADC_meta::showSinglePulse(Pulse *pu) {
	TCanvas *EventCanvas = Canv("EventCanvas", 1, 1, 600, 600);

	const double *data;
//	const double *dataCFD;
	TH1D* ChHisto;
	TH1D* CFDHisto;

	TArrow *arrow = 0;

	long chlength = 300; //tbb unknown
	int nbrpeaks = 0;

	///////////////////////////////////
	//Initialize Histos////////////////
	///////////////////////////////////
	char tmp[256];
	char name[64];
	sprintf(tmp, "channel%i", pu->Pu_from_Ch);
	sprintf(name, "Channel %i", pu->Pu_from_Ch);
	ChHisto = new TH1D(tmp, name, 1, 0, 1);
	CFDHisto = new TH1D(tmp, name, 1, 0, 1);

	//labeling the axis//
	ChHisto->SetXTitle(tmp);
	ChHisto->GetXaxis()->CenterTitle();
	ChHisto->SetYTitle("U [mV]");
	ChHisto->GetYaxis()->CenterTitle();
	ChHisto->SetBit(TH1::kNoStats);

	//--set the histogram to the right size--//
	ChHisto->Reset();
	ChHisto->SetBins(chlength, 0 - 0.5, chlength - 0.5);
	ChHisto->SetMaximum(adc[pu->Pu_from_Ch]->m_adc_settings->fullscale + 100);
	ChHisto->SetMinimum(-100);

	CFDHisto->Reset();
	CFDHisto->SetBins(chlength, 0 - 0.5, chlength - 0.5);
	CFDHisto->SetMaximum(adc[pu->Pu_from_Ch]->m_adc_settings->fullscale + 100);
	CFDHisto->SetMinimum(-100);

	//////////////////////////////
	//Fill the Histos/////////////
	//////////////////////////////
	for (int j = 0; j<chlength; ++j)
	{
		ChHisto->SetBinContent(j + 1, 0);
	}

	double *x1 = 0, *x2 = 0, *y1 = 0, *y2 = 0;

	data = &pu->Waveform[0];

	// CFD
	const long delay = (long)(adc[pu->Pu_from_Ch]->m_adc_settings->delay * 1000 / adc[pu->Pu_from_Ch]->m_adc_settings->samplerate);
	int bin = 0;
	int cfd_k = 0;
	double height = 0, heightCFD = 0;

	for (int k = 0; k < pu->DataLength; ++k)
	{
		bin = int(pu->timestamp / adc[pu->Pu_from_Ch]->m_adc_settings->samplerate) + k + 1;

		height = data[k] / adc[pu->Pu_from_Ch]->m_adc_settings->gain;

		cfd_k = k -/*+*/ delay;
		if (cfd_k >= 0 && cfd_k < pu->DataLength)
			heightCFD = (1.)*data[cfd_k] / adc[pu->Pu_from_Ch]->m_adc_settings->gain;
		else
			heightCFD = 0.;

		ChHisto->SetBinContent(bin, height);
		CFDHisto->SetBinContent(bin, height - heightCFD);
	}

	// draw arrows
	nbrpeaks = pu->NbrPeaks;
	if (pu->NbrPeaks > 0) {
		arrow = new TArrow[pu->NbrPeaks];
		x1 = new double[pu->NbrPeaks];
		x2 = new double[pu->NbrPeaks];
		y1 = new double[pu->NbrPeaks];
		y2 = new double[pu->NbrPeaks];

		for (int peak = 0; peak < pu->NbrPeaks; peak++) {
			arrow[peak].SetAngle(40);
			arrow[peak].SetFillColor(kRed);
			arrow[peak].SetLineColor(kRed);

			x1[peak] = pu->Peaks[peak].Time_ss + pu->timestamp / adc[pu->Pu_from_Ch]->m_adc_settings->samplerate;
			x2[peak] = x1[peak];
			y1[peak] = 0;
			y2[peak] = 1000;
		}
	}

	EventCanvas->cd(1);		//change to the right pad
	ChHisto->Draw("p0c");
	ChHisto->SetLineColor(2);
	CFDHisto->Draw("p0c,same");

	for (int peak = 0; peak < nbrpeaks; peak++)
	{
		if (arrow) arrow[peak].DrawArrow(x1[peak], y1[peak], x2[peak], y2[peak], float(0.01), "<|");
	}




	EventCanvas->cd(0);
	gPad->Update();

	while (!_kbhit())
	{
		gSystem->Sleep(50);
		gSystem->ProcessEvents();
	}

	char ch;
	while (_kbhit()) ch = _getch();

	delete ChHisto;		ChHisto=0; 
	delete CFDHisto;	CFDHisto=0;
	delete EventCanvas;	EventCanvas=0;
	if (x1) delete[] x1;
	if (x2) delete[] x2;
	if (y1) delete[] y1;
	if (y2) delete[] y2;
	if (arrow) delete[] arrow;

	if (ch == 'q' || ch == 'Q')
		return false;				// Quit was pressed..

	return true;
}

bool ADC_meta::showAllPulses(Event *ev) {
	int nbr_channels = 0;
	for (int i = 0; i<max_nbr_channels; i++) {
		if (ev->channels[i] != 0) {
			nbr_channels++;
		}
	}

	TCanvas *EventCanvas = Canv("EventCanvas", 1, nbr_channels);

	const double *data;
	const double *dataCFD;
	std::vector<TH1D*> ChHisto;
	std::vector<TH1D*> CFDHisto;

	TArrow *arrow = 0;

	long chlength = 300;
	int ch_pos = 0;

	int nbrpeaks = 0;

	for (int ch = 0; ch<NUM_CHANNELS; ch++) {
			if (ev->channels[ch] != 0) {
				///////////////////////////////////
				//Initialize Histos////////////////
				///////////////////////////////////
				char tmp[256];
				char name[64];
				char tmp2[256];
				char name2[64];
				sprintf(tmp, "channel%i", ch);
				sprintf(name, "Channel %i", ch);
				sprintf(tmp2, "channelCFD%i", ch);
				sprintf(name2, "ChannelCFD %i", ch);
				ChHisto.push_back(new TH1D(tmp, name, 1, 0, 1));
				CFDHisto.push_back(new TH1D(tmp2, name2, 1, 0, 1));

				//labeling the axis//
				ChHisto[ch_pos]->SetXTitle(tmp);
				ChHisto[ch_pos]->GetXaxis()->CenterTitle();
				ChHisto[ch_pos]->SetYTitle("U [mV]");
				ChHisto[ch_pos]->GetYaxis()->CenterTitle();
				ChHisto[ch_pos]->SetBit(TH1::kNoStats);

				//--set the histogram to the right size--//
				ChHisto[ch_pos]->Reset();
				ChHisto[ch_pos]->SetBins(chlength, 0 - 0.5, chlength - 0.5);
				ChHisto[ch_pos]->SetMaximum(adc[ch]->m_adc_settings->fullscale + 100);
				ChHisto[ch_pos]->SetMinimum(-100);

				CFDHisto[ch_pos]->Reset();
				CFDHisto[ch_pos]->SetBins(chlength, 0 - 0.5, chlength - 0.5);
				CFDHisto[ch_pos]->SetMaximum(adc[ch]->m_adc_settings->fullscale + 100);
				CFDHisto[ch_pos]->SetMinimum(-100);

				//////////////////////////////
				//Fill the Histos/////////////
				//////////////////////////////
				for (int j = 0; j<chlength; ++j)
				{
					ChHisto[ch_pos]->SetBinContent(j + 1, 0);
					CFDHisto[ch_pos]->SetBinContent(j + 1, 0);
				}

				double *x1 = 0, *x2 = 0, *y1 = 0, *y2 = 0;

				if (ev->channels[ch]->NbrPulses > 0)
				{
					for (int pu = 0; pu < ev->channels[ch]->NbrPulses; ++pu)
					{
						data = &ev->channels[ch]->Pulses[pu].Waveform[0];//&channels[ch]->Pulses[pu].CopyOfWaveform[0];
						dataCFD = &ev->channels[ch]->Pulses[pu].Waveform[0];

						// CFD
						const long delay = (long)(adc[ch]->m_adc_settings->delay * 1000 / adc[ch]->m_adc_settings->samplerate);
						int bin = 0;
						int cfd_k = 0;
						double height = 0, heightCFD = 0;

						for (int k = 0; k < ev->channels[ch]->Pulses[pu].DataLength; ++k)
						{
							bin = int(ev->channels[ch]->Pulses[0].timestamp / adc[ch]->m_adc_settings->samplerate) + k;

							height = data[k] / adc[ch]->m_adc_settings->gain;

							cfd_k = k -/*+*/ delay;
							if (cfd_k >= 0 && cfd_k < ev->channels[ch]->Pulses[0].DataLength)
								heightCFD = (1.)*data[cfd_k] / adc[ch]->m_adc_settings->gain;
							else
								heightCFD = 0.;

							ChHisto[ch_pos]->SetBinContent(bin, height);
							CFDHisto[ch_pos]->SetBinContent(bin, height - heightCFD);
						}

						// draw arrows
						nbrpeaks = ev->channels[ch]->Pulses[pu].NbrPeaks;
						if (ev->channels[ch]->Pulses[pu].NbrPeaks > 0) {
							arrow = new TArrow[ev->channels[ch]->Pulses[pu].NbrPeaks];
							x1 = new double[ev->channels[ch]->Pulses[pu].NbrPeaks];
							x2 = new double[ev->channels[ch]->Pulses[pu].NbrPeaks];
							y1 = new double[ev->channels[ch]->Pulses[pu].NbrPeaks];
							y2 = new double[ev->channels[ch]->Pulses[pu].NbrPeaks];

							for (int peak = 0; peak < ev->channels[ch]->Pulses[pu].NbrPeaks; peak++) {
								arrow[peak].SetAngle(40);
								arrow[peak].SetFillColor(kRed);
								arrow[peak].SetLineColor(kRed);

								x1[peak] = ev->channels[ch]->Pulses[pu].Peaks[peak].Cfd + ev->channels[ch]->Pulses[pu].timestamp / adc[ch]->m_adc_settings->samplerate;
								x2[peak] = x1[peak];
								y1[peak] = 0;
								y2[peak] = 1000;

							}
						}
					}
				}

				EventCanvas->cd(ch_pos + 1);		//change to the right pad
				ChHisto[ch_pos]->Draw("p0c");
				ChHisto[ch_pos]->SetLineColor(2);
				ChHisto[ch_pos]->Draw("p0c,same");
				CFDHisto[ch_pos]->SetLineColor(4);
				CFDHisto[ch_pos]->Draw("p0c,same");

				for (int peak = 0; peak < nbrpeaks; peak++)
				{
					if (arrow) arrow[peak].DrawArrow(x1[peak], y1[peak], x2[peak], y2[peak], float(0.01), "<|");
				}

				if (x1) delete[] x1;
				if (x2) delete[] x2;
				if (y1) delete[] y1;
				if (y2) delete[] y2;
				nbrpeaks = 0;


				EventCanvas->cd(0);
				gPad->Update();

				ch_pos++;
			}
		}

		while (!_kbhit())
		{
			gSystem->Sleep(50);
			gSystem->ProcessEvents();
		}

		char ch;
		while (_kbhit()) ch = _getch();

		if (ch == 'q' || ch == 'Q')
			return false;

		ChHisto.clear();
		CFDHisto.clear();

		return true;

}

void ADC_meta::WriteRawPulseInfo(int nbr_events) {
	TFile *RawPulseInfo = new TFile("RawPulseInfo.root", "RECREATE");
	
	TH1D **hist1D[3];
	for (int i = 0; i < 3; i++) {
		hist1D[i] = new TH1D*[max_nbr_channels];
		memset(hist1D[i], 0, max_nbr_channels*sizeof(TH1D*));
	}

	TH2D **hist2D[2];
	for (int i = 0; i < 2; i++) {
		hist2D[i] = new TH2D*[max_nbr_channels];
		memset(hist2D[i], 0, max_nbr_channels*sizeof(TH1D*));
	}

	char name[50];
	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if (used_channels_list[ch] > 0) {
			sprintf(name, "%d_Height", ch);
			hist1D[0][ch] = new TH1D(name, "", 1000, 0, 1000);				//height
			hist1D[0][ch]->SetXTitle("Time [Sample Steps]");
			hist1D[0][ch]->SetYTitle("Height [mV]");
			sprintf(name, "%d_FWHM", ch);
			hist1D[1][ch] = new TH1D(name, "", 250, 0, 50);					//FWHM
			sprintf(name, "%d_Slope", ch);
			hist1D[2][ch] = new TH1D(name, "", 200, 0, 20);					//Slope

			sprintf(name, "%d_FWHM_vs_Height", ch);
			hist2D[0][ch] = new TH2D(name, "", 250, 0, 50, 500, 0, 500);
			sprintf(name, "%d_Slope_vs_Height", ch);
			hist2D[1][ch] = new TH2D(name, "", 200, -1000, 10000, 500, 0, 500);
		}
	}

	for (int i = 0; i < nbr_events; i++) {
		ReadNextEvent(m_event);
		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (m_event->channels[ch] && m_event->channels[ch]->NbrPulses && used_channels_list[ch] > 0) {
				if (!m_event->channels[ch]->Pulses[0].is_TDC_signal) {
				adc[ch]->analyze_Pulse(&m_event->channels[ch]->Pulses[0]);

				if (m_event->channels[ch]->Pulses[0].NbrPeaks == 1) {
					hist1D[0][ch]->Fill(m_event->channels[ch]->Pulses[0].Peaks[0].Height);
					hist1D[1][ch]->Fill(m_event->channels[ch]->Pulses[0].Peaks[0].Fwhm);		
					hist1D[2][ch]->Fill(m_event->channels[ch]->Pulses[0].Peaks[0].Slope);

					hist2D[0][ch]->Fill(m_event->channels[ch]->Pulses[0].Peaks[0].Fwhm, m_event->channels[ch]->Pulses[0].Peaks[0].Height);
					hist2D[1][ch]->Fill(m_event->channels[ch]->Pulses[0].Peaks[0].Slope, m_event->channels[ch]->Pulses[0].Peaks[0].Height);
				}
				}
			}
		}
	}
	
	for (int ch = 0; ch < NUM_CHANNELS; ch++) {
		if (used_channels_list[ch] > 0) {
			hist1D[0][ch]->Write();
			delete hist1D[0][ch]; hist1D[0][ch] = 0;
			hist1D[1][ch]->Write();
			delete hist1D[1][ch]; hist1D[1][ch] = 0;
			hist1D[2][ch]->Write();
			delete hist1D[2][ch]; hist1D[2][ch] = 0;
			hist2D[0][ch]->Write();
			delete hist2D[0][ch]; hist2D[0][ch] = 0;
			hist2D[1][ch]->Write();
			delete hist2D[1][ch]; hist2D[1][ch] = 0;
		}
	}
	
	RawPulseInfo->Close();
	delete RawPulseInfo; RawPulseInfo = 0;
}

void ADC_meta::WriteOscilliscopeMode(int nbr_events, int trigger_ch) {
	TFile *OscilliscopeMode = new TFile("OscilliscopeMode.root", "RECREATE");
	
	TH1D **hist1D[1];
	for (int i = 0; i < 1; i++) {
		hist1D[i] = new TH1D*[max_nbr_channels];
		memset(hist1D[i], 0, max_nbr_channels*sizeof(TH1D*));
	}

	TH2D **hist2D[1];
	for (int i = 0; i < 1; i++) {
		hist2D[i] = new TH2D*[max_nbr_channels];
		memset(hist2D[i], 0, max_nbr_channels*sizeof(TH1D*));
	}

	char name[50];
	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if (used_channels_list[ch] > 0) {
			sprintf(name, "%d_Osz1D", ch);
			hist1D[0][ch] = new TH1D(name, "", 10000, -500, 500);
		
			sprintf(name, "%d_Osz2D", ch);
			hist2D[0][ch] = new TH2D(name, "", 600, -150, 150, 400, -50, 350);
		}
	}

	for (int events = 0; events < nbr_events; events++) {
		ReadNextEvent(m_event);
		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (m_event->channels[ch] && m_event->channels[ch]->NbrPulses && used_channels_list[ch] > 0) {
				if (!m_event->channels[ch]->Pulses[0].is_TDC_signal) {
				adc[ch]->analyze_Pulse(&m_event->channels[ch]->Pulses[0]);
			}
		}
		}

		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (m_event->channels[ch] && m_event->channels[trigger_ch]) {
				if (m_event->channels[ch]->NbrPulses == 1 && m_event->channels[trigger_ch]->NbrPulses > 0) {
					for (int i = 0; i < m_event->channels[ch]->Pulses[0].DataLength; i++) {
						if (m_event->channels[ch]->Pulses[0].NbrPeaks == 1) {
							hist1D[0][ch]->Fill((m_event->channels[ch]->Pulses[0].timestamp + i * adc[ch]->m_adc_settings->samplerate)* 0.001);

							hist2D[0][ch]->Fill((m_event->channels[ch]->Pulses[0].timestamp - m_event->channels[trigger_ch]->Pulses[0].Peaks[0].Cfd + i * adc[ch]->m_adc_settings->samplerate)* 0.001, m_event->channels[ch]->Pulses[0].Waveform[i] / adc[ch]->m_adc_settings->gain, m_event->channels[ch]->Pulses[0].Waveform[i] / adc[ch]->m_adc_settings->gain);
						}
					}
				}
			}
		}
	}

	for (int ch = 0; ch < NUM_CHANNELS; ch++) {
		if (used_channels_list[ch] > 0) {
			hist1D[0][ch]->Write();
			delete hist1D[0][ch]; hist1D[0][ch] = 0;

			hist2D[0][ch]->Write();
			delete hist2D[0][ch]; hist2D[0][ch] = 0;
		}
	}

	OscilliscopeMode->Close();
	delete OscilliscopeMode; OscilliscopeMode = 0;
}

void ADC_meta::WriteSignalInspector(int nbr_events) {
	TFile *SignalInspector = new TFile("SignalInspector.root", "RECREATE");

	TH2D **hist2D[3];
	for (int i = 0; i < 3; i++) {
		hist2D[i] = new TH2D*[max_nbr_channels];
		memset(hist2D[i], 0, max_nbr_channels*sizeof(TH1D*));
	}

	char name[50];
	for (int ch = 0; ch < max_nbr_channels; ch++) {
		if (used_channels_list[ch] > 0) {
			sprintf(name, "%d_signal", ch);
			hist2D[0][ch] = new TH2D(name, "", 200, -100, 100, 200, -50, 350);
			sprintf(name, "%d_norm", ch);
			hist2D[1][ch] = new TH2D(name, "", 2000, -100, 100, 1000, -.5, 1.5);
			sprintf(name, "%d_cfd", ch);
			hist2D[2][ch] = new TH2D(name, "", 2000, -100, 100, 1000, -200., 400.);	
		}
	}


	for (int i = 0; i < nbr_events; i++) {
		ReadNextEvent(m_event);

		for (int ch = 0; ch < max_nbr_channels; ch++) {
			if (m_event->channels[ch] && m_event->channels[ch]->NbrPulses && used_channels_list[ch] > 0) {
				if (!m_event->channels[ch]->Pulses[0].is_TDC_signal) {
				adc[ch]->analyze_Pulse(&m_event->channels[ch]->Pulses[0]);

				if (m_event->channels[ch]->Pulses[0].NbrPeaks == 1) {
					for (int pu = 0; pu < m_event->channels[ch]->NbrPulses; pu++) {
						if (1 == m_event->channels[ch]->Pulses[pu].NbrPeaks) {
							for (int i = 0; i < m_event->channels[ch]->Pulses[0].DataLength; i++) {
								hist2D[0][ch]->Fill(((i - m_event->channels[ch]->Pulses[pu].Peaks[0].Time_ss) * adc[ch]->m_adc_settings->samplerate) *0.001, m_event->channels[ch]->Pulses[0].Waveform[i] / adc[ch]->m_adc_settings->gain);
								hist2D[1][ch]->Fill(((i - m_event->channels[ch]->Pulses[pu].Peaks[0].Time_ss) * adc[ch]->m_adc_settings->samplerate) *0.001, m_event->channels[ch]->Pulses[0].Waveform[i] / adc[ch]->m_adc_settings->gain / m_event->channels[ch]->Pulses[pu].Peaks[0].Height);

								int cfd_time = int(i - adc[ch]->m_adc_settings->delay * 1000 / adc[ch]->m_adc_settings->samplerate);
								if (cfd_time >= 0 && cfd_time < m_event->channels[ch]->Pulses[pu].DataLength)
									hist2D[2][ch]->Fill(((i - m_event->channels[ch]->Pulses[pu].Peaks[0].Time_ss) * adc[ch]->m_adc_settings->samplerate) *0.001, (m_event->channels[ch]->Pulses[0].Waveform[i] - m_event->channels[ch]->Pulses[0].Waveform[cfd_time]) / adc[ch]->m_adc_settings->gain);
								}
							}
						}
					}
				}
			}
		}
	}

	for (int ch = 0; ch < NUM_CHANNELS; ch++) {
		if (used_channels_list[ch] > 0) {
			hist2D[0][ch]->Write();
			delete hist2D[0][ch]; hist2D[0][ch] = 0;

			hist2D[1][ch]->Write();
			delete hist2D[1][ch]; hist2D[1][ch] = 0;

			hist2D[2][ch]->Write();
			delete hist2D[2][ch]; hist2D[2][ch] = 0;
		}
	}

	SignalInspector->Close();
	delete SignalInspector;  SignalInspector = 0;
}


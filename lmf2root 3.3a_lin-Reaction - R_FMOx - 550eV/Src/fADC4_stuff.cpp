#include "OS_Version.h"

#include "resort64c.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
//#include "TMinuit.h"
#include "rootstuff.h"
#ifdef _DEBUG
	#include <assert.h>
#endif
#include "Histo.h"

#include <sys/stat.h> 

#include <math.h>

#include "console.h"

#ifdef LINUX
	#include <stdlib.h>
	#include <string.h>
	#include <termios.h>
#else
	#include <conio.h>
#endif

#ifndef dont_use_MFC
	#include "parallel_sort_class.h"
#endif

#pragma warning(disable : 4996)

#include "LMF_IO.h"
#include "Ueberstruct.h"

using namespace std;

/*
	#ifdef _DEBUG
	#define new DEBUG_NEW
	#undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
	#endif
*/

extern unsigned __int32 ui32buffer[ui32_BUFFERSIZE];
extern __int16 i16buffer[i16_BUFFERSIZE];

__int32 analysis(__int64 eventcounter, double  parameter[],  TTree * Data, Ueberstruct * Ueber);
__int32 sort_and_write_NTuple(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[], double timestamp);

double	get_time_s();
bool FileExists(const char * strFilename);

#ifndef dont_use_MFC
	char * GetApplicationPath();
	bool init_MFC();
#endif

#ifndef LINUX
	__int32 my_kbhit();
#else
	__int32 my_kbhit(void);
#endif




// xxx Robert:
///////////////////////////////////////////////////////////
bool fADC4_LMF_read_event_loop(Ueberstruct * Ueber, double tdc_ns[][NUM_IONS])
///////////////////////////////////////////////////////////
{
	double * parameter = Ueber->parameter;
	__int32 * cnt = Ueber->cnt;
	bool	write_sorted_file = ((parameter[49] > 0.5) ? true:false);
	bool	entmixen =			((parameter[50] > 0.5) ? true:false);
	__int32 merged_channel = __int32(parameter[51]+0.1);
	double	timelimit =			parameter[52];
	__int32 new_demerged_channel = __int32(parameter[53]+0.1);
	__int64 max_events_to_write_into_file = (unsigned __int32)(parameter[56]+0.1);

	double root_file_flush_time = 0;
	double timestamp = 0.;
	__int32 rate_counter_limit = 1000;
	__int32 rate_counter = 0;
	char root_output_file_name_with_extension[300]; root_output_file_name_with_extension[0] = 0;

	while(!Ueber->stop_reading_files) {
		++Ueber->fast_mode;
		if (Ueber->fast_mode >= __int32(parameter[11]+0.01)) Ueber->fast_mode = 0;
		if (parameter[11] < 0.9) Ueber->fast_mode = 1;


		if (Ueber->start_new_root_file) {  // start new output root file
			Ueber->start_new_root_file = false;
			if(Ueber->filenumber == 0) {
				sprintf(root_output_file_name_with_extension,"%s.root",Ueber->root_output_file_name_without_extension);
			} else {
				if (Ueber->outputRootFile_handle) {
					Ueber->outputRootFile_handle->Write();
					Ueber->Hist->Reset();
					Ueber->outputRootFile_handle->Close();
					Ueber->outputRootFile_handle = 0;
				}
				sprintf(root_output_file_name_with_extension,"%s_%03d.root",Ueber->root_output_file_name_without_extension,Ueber->filenumber);
			}

			if (Ueber->config_info) {delete Ueber->config_info; Ueber->config_info = 0;}

			Ueber->outputRootFile_handle =Ueber->rt->RecreateRootFile(root_output_file_name_with_extension,"");
			if (!Ueber->outputRootFile_handle->IsWritable()) {
				printf("Could not open output root file for writing.\n");
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"Could not open output root file for writing.\n"); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				return true;
			}
			if (Ueber->outputRootFile_handle) {
				Ueber->outputRootFile_handle->cd();
				Ueber->Hist->rootfile = Ueber->outputRootFile_handle;
			}

			if (Ueber->parser) {
				Ueber->config_info = new TText(0,0,Ueber->parser->all_read_bytes);
				Ueber->config_info->Write("config_info");
			}

			++Ueber->filenumber;
		}



		ndigo_packet packet;

		memset(cnt,0,sizeof(__int32)*NUM_CHANNELS);
		bool bEnd_of_group_detected = false;
		double timestamp_trigger_ns = 0;
		__int32 TriggerChannel = Ueber->LMF_input->fADC4.iTriggerChannel;
			
		while (true) {
			if (Ueber->LMF_input->ReadNextfADC4packet(&packet, bEnd_of_group_detected, i16buffer, i16_BUFFERSIZE)) {
				bool signal_is_valid = true;

				__int32 channel = packet.card * 5 + packet.channel;
				__int64 signal_timestamp_ps = 0;

				if (packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
					signal_timestamp_ps = packet.timestamp;
				}
			
				if (packet.type == CRONO_PACKET_TYPE_TDC_RISING) {
					signal_timestamp_ps = packet.timestamp;
				}
			
				if (CRONO_PACKET_TYPE_16_BIT_SIGNED) {
					signal_timestamp_ps = packet.timestamp;
					//for (unsigned __int32 i=0;i<packet.length*4;i++) printf("bin %i = %i\n",i,i16buffer[i]);
				}

				if (signal_is_valid) {
					if (cnt[channel] < NUM_IONS) {
						tdc_ns[channel][cnt[channel]] = signal_timestamp_ps *0.001;
						cnt[channel]++;
					}
				} else {
					if (cnt[channel] == 0) {
						if (timestamp_trigger_ns == __int64(0)) timestamp_trigger_ns = signal_timestamp_ps * 0.001;
					}
				}

				if (bEnd_of_group_detected) break;
			} else break;
		}

		if (Ueber->LMF_input->input_lmf->eof) break;
			
		if (cnt[TriggerChannel] > 0) {
			timestamp_trigger_ns = tdc_ns[TriggerChannel][0];
			for (__int32 ch = 0;ch<NUM_CHANNELS;ch++) {
				for (__int32 hit = 0;hit<cnt[ch];hit++) tdc_ns[ch][hit] -= timestamp_trigger_ns;
			}
		} else {
			memset(cnt,0,sizeof(__int32)*NUM_CHANNELS);
		}

		timestamp = timestamp_trigger_ns*1.e9; // convert to seconds
// end of fADC4-stuff




		// demix:
		if (entmixen) {
			double temp_tdc[NUM_IONS];
			__int32 n = cnt[merged_channel];
			for (__int32 i = 0;i < n;++i) temp_tdc[i] = tdc_ns[merged_channel][i];
			cnt[new_demerged_channel] = 0;
			cnt[merged_channel]   = 0;

			for (__int32 i = 0;i<n;++i) {
				double a = temp_tdc[i];
				if (a > timelimit) {
					tdc_ns[new_demerged_channel][cnt[new_demerged_channel]] = a;
					cnt[new_demerged_channel] = cnt[new_demerged_channel] + 1;
				} else {
					tdc_ns[merged_channel][cnt[merged_channel]] = a;
					cnt[merged_channel] = cnt[merged_channel] + 1;
				}
			}
		}

		// Anti-Moire :
		if (parameter[15] > 0.5) {
			for (__int32 i=0;i<NUM_CHANNELS;++i) {for (__int32 j=0;j<cnt[i];++j) 	tdc_ns[i][j] += (double(rand())/RAND_MAX-0.5)*Ueber->LMF_input->tdcresolution;	}
		}

		if (Ueber->stop_reading_files) return true;
		
		if ((++rate_counter) == rate_counter_limit) {
			rate_counter = 0;

			__int32 c = my_kbhit();
			if (c) 	{
					while (my_kbhit());
					bool do_break = false;
					if(c=='q') {printf("\nq was pressed -> skipping this file\npress Q to skip all files.\n"); do_break = true;}
					if(c=='Q') {printf("\nQ was pressed -> skipping all files.\n"); Ueber->stop_reading_files = true; do_break = true;}
					if (do_break) {
						printf("event %I64i of %I64i\n",Ueber->LMF_input->int64_number_of_read_events, Ueber->LMF_input->int64_Numberofevents);
						if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",Ueber->LMF_input->int64_number_of_read_events, Ueber->LMF_input->int64_Numberofevents); fflush(Ueber->logfile_handle);}
						if (Ueber->stop_reading_files) return true;
						break;
					}
			}

			double new_time = get_time_s();
			static double print_time = 0;
			if (new_time > print_time + 1.) {
				print_time = new_time;
				__int32 stars_in_progress_bar = 0;
				if (Ueber->LMF_input->int64_Numberofevents > 0) {
					stars_in_progress_bar = __int32(52. *Ueber->LMF_input->int64_number_of_read_events / Ueber->LMF_input->int64_Numberofevents);
				}
				printf("\r");
				for (__int32 i=0;i<stars_in_progress_bar;++i) printf("*");
				printf("| %i/s   ",__int32(rate_counter_limit/(new_time-Ueber->old_time)));
			}
			rate_counter_limit = __int32(0.25*rate_counter_limit / (new_time-Ueber->old_time));
			if (rate_counter_limit > 10000) rate_counter_limit = 10000;
			if (rate_counter_limit < 100) rate_counter_limit = 100;
			Ueber->old_time = new_time;
			
			if (new_time > root_file_flush_time + 30) {
				root_file_flush_time = new_time;
				//if (Ueber->outputRootFile_handle) {Ueber->outputRootFile_handle->Flush();}  // Flushing does not yield any advantage
			}
		}

		++Ueber->eventcounter;
		if (parameter[58] > 0 && Ueber->eventcounter > (__int64)(parameter[58])) return true;


		// xxx Robert: hier weiter. sortieren/Analyse...


	} // end while ...

	return true;
}












// xxx Robert:
///////////////////////////////////////////////////////////
bool extract_mean_pulse_fADC4(char * input_LMF_file_name, Ueberstruct * Ueber, double tdc_ns[][NUM_IONS])
///////////////////////////////////////////////////////////
{
	sprintf(Ueber->input_data_file_name, input_LMF_file_name);
	if (Ueber->LMF_output_file_name_without_extension[0] == 0) Ueber->output_LMF_file_index = 0;
	double * parameter = Ueber->parameter;
	__int32 * cnt = Ueber->cnt;
	if (Ueber->stop_reading_files) return true;

	__int64 max_events_to_write_into_file = (unsigned __int32)(parameter[56]+0.1);

	__int32 new_demerged_channel = __int32(parameter[53]+0.1);
	bool entmixen = ((parameter[50] > 0.5) ? true:false);
	__int32 merged_channel = __int32(parameter[51]+0.1);
	double timelimit = parameter[52];

	if (Ueber->LMF_input) delete Ueber->LMF_input;
	Ueber->LMF_input = new LMF_IO(NUM_CHANNELS,NUM_IONS);
	Ueber->fast_mode = parameter[11] > 0.5 ? 1:0;

	Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions = __int32(parameter[9] + 0.1);
	if (parameter[9] < -0.1) Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions--;

	if ( parameter[58] > 0 ) printf("\n\n%G events left (out of %G)\n",parameter[58] - Ueber->eventcounter,parameter[58]);

	printf("\n\nInput file is: %s\n",Ueber->input_data_file_name);
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n\nreading input file %s\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
	if(!FileExists(Ueber->input_data_file_name)) {
		printf("\n \n %s could not be opened!\n\n",Ueber->input_data_file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
		return false;
	}

	if (!Ueber->LMF_input->OpenInputLMF(Ueber->input_data_file_name)) {
		printf("\n \n %s could not be opened!\n\n",Ueber->input_data_file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
		Ueber->stop_reading_files = true;
		return false;
	}

	Ueber->number_of_channels	= ((__int32(parameter[54]+0.1) == 0) ? (Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2) : __int32(parameter[54]+0.1));
	__int32 number_of_hits		= ((__int32(parameter[55]+0.1) == 0) ? (Ueber->LMF_input->max_number_of_hits) : __int32(parameter[55]+0.1));

	double timestamp = 0.;

	Ueber->eventcounter = 0;
	while (Ueber->eventcounter < 10000) {
		++Ueber->eventcounter;

		ndigo_packet packet;

		memset(cnt,0,sizeof(__int32)*NUM_CHANNELS);
		bool bEnd_of_group_detected = false;
		double timestamp_trigger_ns = 0;
		__int32 TriggerChannel = Ueber->LMF_input->fADC4.iTriggerChannel;
			
		while (true) {
			if (Ueber->LMF_input->ReadNextfADC4packet(&packet, bEnd_of_group_detected, i16buffer, i16_BUFFERSIZE)) {
				bool signal_is_valid = true;

				__int32 channel = packet.card * 5 + packet.channel;
				__int64 signal_timestamp_ps = 0;

				if (packet.type == CRONO_PACKET_TYPE_TDC_FALLING) {
					signal_timestamp_ps = packet.timestamp;
				}
			
				if (packet.type == CRONO_PACKET_TYPE_TDC_RISING) {
					signal_timestamp_ps = packet.timestamp;
				}
			
				if (CRONO_PACKET_TYPE_16_BIT_SIGNED) {
					signal_timestamp_ps = packet.timestamp;
					//for (unsigned __int32 i=0;i<packet.length*4;i++) printf("bin %i = %i\n",i,i16buffer[i]);
				}

				if (signal_is_valid) {
					if (cnt[channel] < NUM_IONS) {
						tdc_ns[channel][cnt[channel]] = signal_timestamp_ps *0.001;
						cnt[channel]++;
					}
				} else {
					if (cnt[channel] == 0) {
						if (timestamp_trigger_ns == __int64(0)) timestamp_trigger_ns = signal_timestamp_ps * 0.001;
					}
				}

				if (bEnd_of_group_detected) break;
			} else break;
		}

		if (Ueber->LMF_input->input_lmf->eof) break;
			
		if (cnt[TriggerChannel] > 0) {
			timestamp_trigger_ns = tdc_ns[TriggerChannel][0];
			for (__int32 ch = 0;ch<NUM_CHANNELS;ch++) {
				for (__int32 hit = 0;hit<cnt[ch];hit++) tdc_ns[ch][hit] -= timestamp_trigger_ns;
			}
		} else {
			memset(cnt,0,sizeof(__int32)*NUM_CHANNELS);
		}

		timestamp = timestamp_trigger_ns*1.e9; // convert to seconds




		// xxx Robert: hier gehts weiter



	} // end of while loop


	if (Ueber->LMF_input) {
		Ueber->LMF_input->CloseInputLMF();
		delete Ueber->LMF_input;	Ueber->LMF_input = 0;
	}

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"finished creating mean pulse form\n"); fflush(Ueber->logfile_handle);}


	return true;
}




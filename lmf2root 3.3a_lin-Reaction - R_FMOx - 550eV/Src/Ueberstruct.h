#ifndef UEBERSTRUCT_ALREADY_INCLUDED
	#define UEBERSTRUCT_ALREADY_INCLUDED

	#define ui32_BUFFERSIZE 10000
	#define i16_BUFFERSIZE 10000

	#include "OS_Version.h"
	#include "config_file_parser.h"
	#include "CH_Interface.h"
    #include "../ipa/ipa.h"
	#include "Histo.h"
	#include "LMF_IO.h"
	#include "TText.h"
	#include "rootstuff.h"
	#include "DetStruct.h"
	#include "HDF5_IO.h"

	using namespace CH;

	struct Ueberstruct {
		void * adc_meta;
		CH_Int *CH;
		bool use_CH;
		int CH_status;
		ndigo_packet n_packet;

		sort_class * current_proj_Sorter;
		sort_class * current_elec_Sorter;
		sort_class * current_rec_Sorter;

		sort_class * proj_Sorter;
		sort_class * elec_Sorter;
		sort_class * rec_Sorter;

		Det_struct *proj;
		Det_struct *elec;
		Det_struct *rec;

		double * parameter;

		peak_tracker_class * peak_tracker1;
		peak_tracker_class * peak_tracker2;

			void * parallel_sorter;

		bool start_new_root_file;
		bool stop_reading_files;
		bool write_this_event_into_sorted_LMF;

		Histo * Hist;

		FILE * logfile_handle;
		bool ext_logging;

		__int64 config_version;
		__int32 output_LMF_file_index;
		__int64 EntriesInFile;
		__int64 eventswritten;
		__int64 eventcounter;
		__int32 filenumber;

		bool input_is_HDF5;

		__int32 number_of_channels;

		double old_time;

		config_file_parser * parser;

		TText * config_info;
		TFile * outputRootFile_handle;
		char applicationfilepath[300];
		char logfilepathname[300];
		char root_output_file_name_without_extension[300];
		char LMF_output_file_name_without_extension[300];
		LMF_IO * LMF_output;
		LMF_IO * LMF_input;

		HDF5_IO * HDF5_input;

		rootstuff * rt;
		__int32 * cnt;
		double * tdc_ns;
		__int32 fast_mode;

		char config_file_name[300];

		int ElecInNTuple;
		int RecInNTuple;
		int ProjInNTuple;
		
		int scan_channel;
		int scan_num_vals;
		int scan_factor;
		bool scan_in_data;
		int scan_par_nums[16];

		ipa_class * IPA;
		bool ipa_enable;
		__int32 ipa_channel;
		__int32 ipa_num_events;

		__int32 reverse_time_direction_in_event;

		char ** processed_data_file_names;
		char input_data_file_name[300];

		__int32 bot_mode;

		__int32 error;

		__int32 DEBUG_counter_xxx;
	};

#endif


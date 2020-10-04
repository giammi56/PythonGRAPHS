
#define LMF2ROOTVERSION (3.34)

#include "OS_Version.h"



#include "resort64c.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TSystem.h"
//#include "TMinuit.h"
#include "rootstuff.h"
#ifdef _DEBUG
	#include <assert.h>
#endif
#include "Histo.h"




#include <sys/stat.h> 

#include <math.h>
#include <iostream>

#include "direct.h"
#include "console.h"
#include <stdlib.h>



#ifdef LINUX
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
#include "HDF5_IO.h"
#include "Ueberstruct.h"
#include "../ADC/ADC_meta.h"


#include "../ipa/ipa.h"

//tbb
//#include "../ADC/Event.h"
//#include "../ADC/ADC_analysis.h"


using namespace std;

/*
	#ifdef _DEBUG
	#define new DEBUG_NEW
	#undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
	#endif
*/


ADC_meta *ADC_meta::inst = 0;



#define ui32_BUFFERSIZE 10000
unsigned __int32 ui32buffer[ui32_BUFFERSIZE];

#define i16_BUFFERSIZE 10000
__int16 i16buffer[i16_BUFFERSIZE];

void extract_tdc_data(unsigned __int32 data[],__int32 StartED,__int32 StopED);
__int32 analysis(__int64 eventcounter, double  parameter[],  TTree * Data, Ueberstruct * Ueber);
__int32 sort_and_write_NTuple(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[], double timestamp);
void print_help();

// fADC4 functions:
bool extract_mean_pulse_fADC4(char * input_LMF_file_name, Ueberstruct * Ueber, double tdc_ns[][NUM_IONS]);
bool fADC4_LMF_read_event_loop(Ueberstruct * Ueber, double tdc_ns[][NUM_IONS]);





#ifndef dont_use_MFC
	bool init_MFC();
#endif


#ifndef LINUX

	__int32 my_kbhit()
	{
		if (!_kbhit()) return 0;
		__int32 c = _getch();
		while (_kbhit()) c = _getch(); // empty keyboard buffer
		return c;
	}

#else

	__int32 my_kbhit(void)
	{
		struct termios term, oterm;
		__int32 fd = 0;
		__int32 c = 0;
		tcgetattr(fd, &oterm);
		memcpy(&term, &oterm, sizeof(term));
		term.c_lflag = term.c_lflag & (!ICANON);
		term.c_cc[VMIN] = 0;
		term.c_cc[VTIME] = 1;
		tcsetattr(fd, TCSANOW, &term);
		c = getchar();
		tcsetattr(fd, TCSANOW, &oterm);
		if (c == -1) return 0;
		return c;
	}

#endif



////////////////////////////////////////////////////
double	get_time_s()
////////////////////////////////////////////////////
{
	static unsigned __int64 start_time = 0;
	if (!start_time) start_time = GetTickCount64();
	return 0.001 * double(GetTickCount64()-start_time);
}








/////////////////////////////////////////////////////////////////////////////
bool FileExists(const char * strFilename)
/////////////////////////////////////////////////////////////////////////////
{
	#ifdef LINUX
		if ( access(strFilename,F_OK) )
			return(false);
		else if ( access(strFilename,R_OK) ) {
			return(false);
		} else	return(true);

	#else
		FILE * fi;
		fi = fopen(strFilename,"rb");
		if (!fi) return false;
		if (ferror(fi)) {
			clearerr(fi); fclose(fi);	fi = 0;
			return false;
		}

		if (fi) {
			fclose(fi);	fi = 0;
		}
		return true;
	#endif
}



void delete_file(const char* name)
{
	if (!FileExists(name)) return;
	char name_[500];
#ifdef LINUX
	sprintf(name_,"rm %s",name);
#else
	sprintf(name_,"del %s",name);
#endif
	system(name_);
}




/////////////////////////////////////////////////////////////////////////////
bool File_is_readable(const char * strFilename)
/////////////////////////////////////////////////////////////////////////////
{
	#ifdef LINUX
		return FileExists(strFilename);
	#else
		FILE * fi;
		fi = fopen(strFilename,"rb");
		if (!fi) return false;
		if (ferror(fi)) {
			clearerr(fi); fclose(fi);	fi = 0;
			return false;
		}

		__int8 byte_dummy = 0;
		if (fread(&byte_dummy,1,1,fi) != 1) {
			clearerr(fi); fclose(fi);	fi = 0;
			return false;
		}
			fclose(fi);	fi = 0;
		return true;
	#endif
}



///////////////////////////////////////////////////////////////////////////
void check_for_pause_file(Ueberstruct * Ueber)
///////////////////////////////////////////////////////////////////////////
{
	static double old_time = 0.;
	double new_time = get_time_s();
	if (new_time - old_time > 5) {
		old_time = new_time;
		char c = 0;
		delete_file("are_you_alive.txt");
		while (true) {
			if (FileExists("pause.txt")) {
				for (__int32 i=0;i<5*4;i++) { // wait for 5 seconds
					gSystem->Sleep(250);
					c = my_kbhit();
					if (c != 'q' && c != 'Q') {
						c = 0;
					} else {
						Ueber->stop_reading_files = true;
						break;
					}
				}
				if (!c) continue;
			}
			old_time = get_time_s();
			break;
		}
	}
}










/////////////////////////////////////////////////////////////////////////////
void write_correction_tables(Ueberstruct * Ueber, std::string det_name, Det_struct * det)
/////////////////////////////////////////////////////////////////////////////
{
	char fname[500];

	// now write sum correction table:
	fname[0] = 0;

#ifndef dont_use_MFC
	sprintf(fname,"%scorrection_table_%s.txt", Ueber->applicationfilepath, det_name.c_str());
#endif
	printf("writing sum correction tables for detector \"%s\"\ninto the file:\n\"%s\"\n",det_name.c_str(),fname);
	FILE * ofile = 0;
	ofile = fopen(fname,"wt");
	if (!ofile || !FileExists(fname)) {
		printf("\ncould not open file\n\"%s\"\nfor writing correction table%c\n",fname,7);
		gSystem->Sleep(2000);
		return;
	}

	std::string layer_name;
	profile_class * prof;

	if (det->sorter->sum_walk_calibrator) {
		for (__int32 layer = 0;layer < 3; ++layer) {
			if (layer == 0 && ((!det->sorter->use_u1) || (!det->sorter->use_u2))) continue;
			if (layer == 1 && ((!det->sorter->use_v1) || (!det->sorter->use_v2))) continue;
			if (layer == 2 && ((!det->sorter->use_w1) || (!det->sorter->use_w2))) continue;
			if (layer == 0) {layer_name = "u"; prof = det->sorter->sum_walk_calibrator->sumu_profile;}
			if (layer == 1) {layer_name = "v"; prof = det->sorter->sum_walk_calibrator->sumv_profile;}
			if (layer == 2) {layer_name = "w"; prof = det->sorter->sum_walk_calibrator->sumw_profile;}
			__int32 n = prof ? prof->number_of_columns : 0;
			fprintf(ofile,"\n//sum correction_points layer %s on %s detector:\n\n",layer_name.c_str(),det_name.c_str());
			for (__int32 i=0;i < n; ++i) {
				double x,y;
				det->sorter->sum_walk_calibrator->get_correction_point(x,y,i,layer);
				if (fabs(y) > 50) continue;
				fprintf(ofile,"set_point sum %s %s %lf \t%lf;\n",det_name.c_str(),layer_name.c_str(),x,y);
			}
		}
	}

	fprintf(ofile,"\n");

	if (det->sorter->scalefactors_calibrator && det->sorter->pos_walk_calibrator) {
		if (! (det->sorter->pos_walk_calibrator->error || det->sorter->scalefactors_calibrator->error) ) {
			if (det->sorter->pos_walk_calibrator->calibration_succeeded && det->sorter->scalefactors_calibrator->calibration_succeeded) {
				for (__int32 layer = 0;layer < 3; ++layer) {
					if (layer == 0) {layer_name = "u";}
					if (layer == 1) {layer_name = "v";}
					if (layer == 2) {layer_name = "w";}

					fprintf(ofile,"\n//pos correction_points layer %s on %s detector:\n\n",layer_name.c_str(),det_name.c_str());
					for (__int32 i=0;i < det->sorter->pos_walk_calibrator->number_of_columns; ++i) {
						double x,y;
						det->sorter->pos_walk_calibrator->get_correction_point(x,y,i,layer);
						if (fabs(y) > 50) continue;
						fprintf(ofile,"set_point pos %s %s %lf \t%lf;\n",det_name.c_str(),layer_name.c_str(),x,y);
					}
				}
			}
		}
	}

	if (ofile) {fclose(ofile); ofile = 0;}
}











/////////////////////////////////////////////////////////////////////////////
bool fill_parameters_into_sorters(sort_class * sorter,Det_struct * detector, __int32 parameter_offset)
/////////////////////////////////////////////////////////////////////////////
{
	#ifdef _DEBUG
		assert(detector);
		assert(detector->Ueberstruct_pointer);
	#endif

	if (!detector) return false;

	Ueberstruct * Ueber = (Ueberstruct*)detector->Ueberstruct_pointer;

	if (!Ueber) return false;

	__int32 poffs = 0;
	if (Ueber->config_version < 201401221558) poffs = -2;

	double * parameter = Ueber->parameter;

	if (!parameter) return false;

	detector->use_this_detector = false;

	if (parameter[parameter_offset] < 0.5) return true;

	detector->use_this_detector = true;
	detector->hex_offset_w			= 0.;
	detector->center_X				= parameter[parameter_offset+26+poffs];
	detector->center_Y				= parameter[parameter_offset+27+poffs];
	detector->number_of_reconstructed_hits = 0;
	detector->auto_calibration		= ((parameter[parameter_offset+28+poffs]>0.5)?1:0);
	detector->use_reconstruction		= ((parameter[parameter_offset+29+poffs]>0.5) ? true:false);
	if (detector->auto_calibration) detector->use_reconstruction = false;

	if (sorter) {
		sorter->use_HEX		= false;
		sorter->common_start_mode	= ((parameter[parameter_offset+1] < 0.5) ? true : false);
		sorter->Cu1			= __int32(parameter[parameter_offset+2]*1.0001);
		sorter->Cu2			= __int32(parameter[parameter_offset+3]*1.0001);
		sorter->Cv1			= __int32(parameter[parameter_offset+4]*1.0001);
		sorter->Cv2			= __int32(parameter[parameter_offset+5]*1.0001);
		sorter->Cmcp		= __int32(parameter[parameter_offset+8]*1.0001);

		bool channel_err = false;
		if (sorter->Cu1 >= 0 && sorter->Cu1 == sorter->Cu2) channel_err = true;
		if (sorter->Cu1 >= 0 && sorter->Cu1 == sorter->Cv1) channel_err = true;
		if (sorter->Cu1 >= 0 && sorter->Cu1 == sorter->Cv2) channel_err = true;
		if (sorter->Cu1 >= 0 && sorter->Cu1 == sorter->Cmcp) channel_err = true;
		if (sorter->Cu2 >= 0 && sorter->Cu2 == sorter->Cv1) channel_err = true;
		if (sorter->Cu2 >= 0 && sorter->Cu2 == sorter->Cv2) channel_err = true;
		if (sorter->Cu2 >= 0 && sorter->Cu2 == sorter->Cmcp) channel_err = true;
		if (sorter->Cv1 >= 0 && sorter->Cv1 == sorter->Cv2) channel_err = true;
		if (sorter->Cv1 >= 0 && sorter->Cv1 == sorter->Cmcp) channel_err = true;
		if (sorter->Cv2 >= 0 && sorter->Cv2 == sorter->Cmcp) channel_err = true;

		if (sorter->Cu1<0) {sorter->use_u1 = false; sorter->Cu1 *= -1;}
 		if (sorter->Cu2<0) {sorter->use_u2 = false; sorter->Cu2 *= -1;}
		if (sorter->Cv1<0) {sorter->use_v1 = false; sorter->Cv1 *= -1;}
		if (sorter->Cv2<0) {sorter->use_v2 = false; sorter->Cv2 *= -1;}
		if (sorter->Cmcp<0) {sorter->use_MCP = false; sorter->Cmcp *= -1;}

		sorter->use_MCP		= ((parameter[parameter_offset+9] > 0.5) ? true:false);
		sorter->fu			= parameter[parameter_offset+10];
		sorter->fv			= parameter[parameter_offset+11];
		sorter->uncorrected_time_sum_half_width_u	= parameter[parameter_offset+17];
		sorter->uncorrected_time_sum_half_width_v	= parameter[parameter_offset+18];
		if (sorter->uncorrected_time_sum_half_width_u < 0.001) {printf("\nError: parameter %i is too small\n",parameter_offset+17); return false;}
		if (sorter->uncorrected_time_sum_half_width_v < 0.001) {printf("\nError: parameter %i is too small\n",parameter_offset+18); return false;}
		sorter->dead_time_anode	= parameter[parameter_offset+20];
		sorter->dead_time_mcp	= parameter[parameter_offset+21];
		sorter->MCP_radius		= parameter[parameter_offset+22];
		sorter->runtime_u		= parameter[parameter_offset+23];
		sorter->runtime_v = sorter->runtime_w = 0.;
		if (Ueber->config_version >= 201401221558) {
			sorter->runtime_v		= parameter[parameter_offset+24];
			sorter->runtime_w		= parameter[parameter_offset+25];	
		}
		if (sorter->runtime_v == 0.) sorter->runtime_v = sorter->runtime_u;
		if (sorter->runtime_w == 0.) sorter->runtime_w = sorter->runtime_u;
		
		sorter->use_sum_correction		= ((parameter[parameter_offset+30+poffs]>0) ? true:false);
		sorter->use_pos_correction		= ((parameter[parameter_offset+31+poffs]>0) ? true:false);
		if (detector->auto_calibration) sorter->use_sum_correction = false;
		if (detector->auto_calibration) sorter->use_pos_correction = false;

		sorter->count = Ueber->cnt;
		sorter->tdc_array_row_length = NUM_IONS;
		sorter->tdc_pointer = Ueber->tdc_ns;

		if (parameter[parameter_offset] > 1.5) {// parameter for HEX_-Detector:
			sorter->use_HEX		= true;
			sorter->Cw1			= __int32(parameter[parameter_offset+6]*1.0001);
			sorter->Cw2			= __int32(parameter[parameter_offset+7]*1.0001);

			if (sorter->Cw1 >= 0) {
				if (sorter->Cw1 == sorter->Cu1) channel_err = true;
				if (sorter->Cw1 == sorter->Cu2) channel_err = true;
				if (sorter->Cw1 == sorter->Cv1) channel_err = true;
				if (sorter->Cw1 == sorter->Cv2) channel_err = true;
				if (sorter->Cw1 == sorter->Cw2) channel_err = true;
				if (sorter->Cw1 == sorter->Cmcp) channel_err = true;
			}
			if (sorter->Cw2 >= 0) {
				if (sorter->Cw2 == sorter->Cu1) channel_err = true;
				if (sorter->Cw2 == sorter->Cu2) channel_err = true;
				if (sorter->Cw2 == sorter->Cv1) channel_err = true;
				if (sorter->Cw2 == sorter->Cv2) channel_err = true;
				if (sorter->Cw2 == sorter->Cmcp) channel_err = true;
			}

			if (sorter->Cw1<0) {sorter->use_w1 = false; sorter->Cw1 *= -1;}
			if (sorter->Cw2<0) {sorter->use_w2 = false; sorter->Cw2 *= -1;}

			sorter->fw			= parameter[parameter_offset+12];
			sorter->uncorrected_time_sum_half_width_w	= parameter[parameter_offset+19];
			if (sorter->uncorrected_time_sum_half_width_w < 0.001) {printf("\nError: parameter %i is too small\n",parameter_offset+19); return false;}
		}
		if (channel_err) {
			printf("error: same channel numbers in one of the detectors\n");
			return false;
		}

		sorter->run_without_sorting_flag = !detector->use_reconstruction;
		if (sorter->run_without_sorting_flag) sorter->dont_overwrite_original_data = true;
	}

	// detector:
	detector->offset_timesum_U		= parameter[parameter_offset+13];
	detector->offset_timesum_V		= parameter[parameter_offset+14];
	detector->offset_timesum_W		= 0.;

	if (parameter[parameter_offset] > 1.5) {// parameter for HEX_-Detector:
		// detector for HEX:
		detector->offset_timesum_W		= parameter[parameter_offset+15];
		detector->hex_offset_w			= parameter[parameter_offset+16];
	}

	detector->use_sum_tracker_layer_u = parameter[parameter_offset+32+poffs] > 0.5 ? true : false;
	detector->use_sum_tracker_layer_v = parameter[parameter_offset+33+poffs] > 0.5 ? true : false;
	detector->use_sum_tracker_layer_w = parameter[parameter_offset+34+poffs] > 0.5 ? true : false;
	detector->sumu_watchdog->offset_to_user_in_first_run = parameter[parameter_offset+35+poffs];
	detector->sumv_watchdog->offset_to_user_in_first_run = parameter[parameter_offset+36+poffs];
	detector->sumw_watchdog->offset_to_user_in_first_run = parameter[parameter_offset+37+poffs];

	detector->sorter = sorter;

	return true;
}















////////////////////////////////////////////////////////////////////
void ProcessRootFile(char * name, Ueberstruct * Ueber)
////////////////////////////////////////////////////////////////////
{
	if (Ueber->stop_reading_files) return;
	double * parameter = Ueber->parameter;

	Ueber->fast_mode = parameter[11] > 0.5 ? 1:0;

	double root_file_flush_time = 0.;
	__int32 rate_counter_limit = 5000;
#ifdef _DEBUG
	rate_counter_limit = 200;
#endif
	__int32 rate_counter = 0;

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"reading %s ...\n",name); fflush(Ueber->logfile_handle);}

	if (parameter[58] > 0)	printf("\n\n%G events left (out of %G)\n",parameter[58] - Ueber->eventcounter,parameter[58]);

	printf("\n\nInput file is: %s\n",name);
	
	if(!File_is_readable(name)) {
		Red(true);
		printf("\n \n Error: %s does not exist!\n\n",name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s does not exist!\n\n",name); fflush(Ueber->logfile_handle);}
		White(false);
		return;
	}
	char root_output_file_name_with_extension[517];
	static TTree * RootTree = 0;
	if (RootTree) {delete RootTree; RootTree = 0;}
	RootTree = Ueber->rt->OpenRootFileGetTree(name,"Data");
	if(!RootTree){
		printf("\n could not find NTuple \"Data\" in this rootfile:\n%s\n",name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n could not find NTuple \"Data\" in this rootfile:\n%s\n",name); fflush(Ueber->logfile_handle);}
		return;
	}

	Ueber->EntriesInFile = 0;

	__int64 NumEvents = (__int64)(RootTree->GetEntries());

	printf("Number of events in file: %.0lf\n", double(NumEvents));
	printf("\n0%%           25%%          50%%          75%%        100%%\n");
	printf("|------------|------------|------------|-----------|\n");
	
	for(__int64 NumEventsSoFar=0;NumEventsSoFar<NumEvents;++NumEventsSoFar) {
		++Ueber->fast_mode;
		if (Ueber->fast_mode >= __int32(parameter[11]+0.01)) Ueber->fast_mode = 0;
		if (parameter[11] < 0.9) Ueber->fast_mode = 1;

		if(Ueber->start_new_root_file) {  // start new output root file
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

			Ueber->outputRootFile_handle = Ueber->rt->RecreateRootFile(root_output_file_name_with_extension,"");
			if (!Ueber->outputRootFile_handle->IsWritable()) {
				printf("Could not open output root file for writing.\n");
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"Could not open output root file for writing.\n"); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				return;
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
			printf ("Output file name is %s\n", root_output_file_name_with_extension);
		}

analysis(Ueber->eventcounter,Ueber->parameter,RootTree,Ueber);		// jumps into the analysis

		++Ueber->eventcounter;

		if (Ueber->parameter[58] > 0 && Ueber->eventcounter > (__int64)(Ueber->parameter[58])) {
			Ueber->stop_reading_files = true;
			return;
		}

		// Invoke IPA as number of requested events were read...
		if(Ueber->ipa_enable) {
			if(Ueber->IPA->events_in_memory() >= Ueber->ipa_num_events) {
//				Ueber->IPA->main_loop();
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",NumEventsSoFar, NumEvents); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				return;
			}
		}

		if ((++rate_counter) == rate_counter_limit) {
			rate_counter = 0;

			__int32 c = my_kbhit();
			if (c) {
				bool do_break = false;
				if (c=='q') {printf("\nq was pressed -> skipping this file\npress Q to skip all files.\n"); do_break = true;}
				if (c=='Q') {printf("\nQ was pressed -> skipping all files.\n"); Ueber->stop_reading_files = true; do_break = true;}
				if (c=='p') {
					printf("Pause. Hit any key to continue.");
					__int32 i=0;
					while (true) {
						i++;
						if (i > 4*4) {
							i=0;
							delete_file("are_you_alive.txt");
						}
						char c = my_kbhit();
						if (!c) {gSystem->Sleep(250); continue;}
						printf("\r");
						if (c=='q') {printf("\nq was pressed -> skipping this file\npress Q to skip all files.\n"); do_break = true;}
						if (c=='Q') {printf("\nQ was pressed -> skipping all files.\n"); Ueber->stop_reading_files = true; do_break = true;}
						break;
					}
				}
				// User invoke IPA
				if(Ueber->ipa_enable) {
					if(c=='i' ||c=='I') {
//						Ueber->IPA->main_loop();
						Ueber->stop_reading_files = true; do_break = true;
					}
				}

				if (do_break) {
					printf("event %I64i of %I64i\n",NumEventsSoFar, NumEvents);
					if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",NumEventsSoFar, NumEvents); fflush(Ueber->logfile_handle);}
					if (Ueber->stop_reading_files) return;
					break;
				}
			}

			double new_time = get_time_s();
			static double print_time = 0;
			if (new_time > print_time + 1.) {
				check_for_pause_file(Ueber);
				print_time = new_time;
				__int32 stars_in_progress_bar = NumEvents > 0 ? (__int32(52. * NumEventsSoFar / NumEvents)) : 0; 
				printf("\r");
				for (__int32 i=0;i<stars_in_progress_bar;++i) printf("*");
				printf("| %i/s   ",__int32(rate_counter_limit/(new_time-Ueber->old_time)));
			}
#ifndef _DEBUG
			rate_counter_limit = __int32(0.25*rate_counter_limit / (new_time-Ueber->old_time));
			if (rate_counter_limit > 10000) rate_counter_limit = 10000;
			if (rate_counter_limit < 100) rate_counter_limit = 100;
#endif
			Ueber->old_time = new_time;
			
			if (new_time > root_file_flush_time + 30) {
				root_file_flush_time = new_time;
				if (Ueber->outputRootFile_handle) {Ueber->outputRootFile_handle->Flush();}
			}
		}

	}

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"finished reading the root file\n"); fflush(Ueber->logfile_handle);}

/*
	if (Ueber->bot_mode) {
		char command[300]; command[0] = 0;
		sprintf(command,"del %s", name);
		system(command);
	}
*/
}









///////////////////////////////////////////////////////////
bool LMF_read_event_loop(Ueberstruct * Ueber)
///////////////////////////////////////////////////////////
{
	double * parameter = Ueber->parameter;
	double * tdc_ns				= Ueber->tdc_ns;
	__int32 * cnt = Ueber->cnt;
	Ueber->current_proj_Sorter	= Ueber->proj_Sorter;
	Ueber->current_elec_Sorter	= Ueber->elec_Sorter;
	Ueber->current_rec_Sorter	= Ueber->rec_Sorter;

	bool	write_sorted_file = ((parameter[49] > 0.5) ? true:false);
	bool	entmixen =			((parameter[50] > 0.5) ? true:false);
	__int32 merged_channel = __int32(parameter[51]+0.1);
	double	timelimit =			parameter[52];
	__int32 new_demerged_channel = __int32(parameter[53]+0.1);
	__int64 max_events_to_write_into_file = (unsigned __int32)(parameter[56]+0.1);

	double root_file_flush_time = 0.;
	double timestamp = 0.;
	__int32 rate_counter_limit = 1000;
#ifdef _DEBUG
	rate_counter_limit = 200;
#endif
	__int32 rate_counter = 0;
	__int32 sleep0_counter			= 0;
	char root_output_file_name_with_extension[300]; root_output_file_name_with_extension[0] = 0;


	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT) Ueber->stop_reading_files = true;  	// ADC - lmf2root loophole since no AGAT files are supported..
	}







	////////////////////////////////////////////////////

	while(!Ueber->stop_reading_files) {

		// if (Ueber->eventcounter > 10000) Ueber->stop_reading_files = true; // xxx


		if (Ueber->parallel_sorter) {
			parallel_sort_class * parallel = (parallel_sort_class *)Ueber->parallel_sorter;
	
			__int32 ID = parallel->next_thread_ID_to_be_sorted;

			Ueber->tdc_ns				= parallel->events[ID]->tdc_ns;
			Ueber->cnt					= parallel->events[ID]->cnt;
			tdc_ns						= Ueber->tdc_ns;
			cnt							= Ueber->cnt;
			Ueber->current_proj_Sorter	= parallel->events[ID]->sorter_proj;
			Ueber->current_rec_Sorter	= parallel->events[ID]->sorter_rec;
			Ueber->current_elec_Sorter	= parallel->events[ID]->sorter_elec;
		}

		sleep0_counter++;
		if (sleep0_counter > 1000) {
			sleep0_counter = 0;
			#ifndef dont_use_MFC
				Sleep(0); // gives other processes the chance to do get CPU access
			#endif
		}
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

		bool event_is_valid = true;



//// TDC - TDC - TDC - TDC - TDC - TDC - TDC - TDC - TDC - TDC - TDC
	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8HPRAW ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8HP ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8 ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_2TDC8 ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_RAW32BIT ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_SIMPLE)
			{
				if (!Ueber->LMF_input->ReadNextEvent()) {
					break;
				}
				Ueber->LMF_input->GetTDCDataArray(tdc_ns);
				Ueber->LMF_input->GetNumberOfHitsArray(cnt);
				__int32 max_chan;
				if (NUM_CHANNELS < Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2) {
					max_chan = NUM_CHANNELS;
				} else {
					max_chan = Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2;
				}
			
				for (__int32 i=0;i<max_chan;++i) {
					__int32 n = cnt[i];
					for (__int32 j=0;j<n;++j)  {
						if(i!=Ueber->scan_channel) 
							tdc_ns[i*NUM_CHANNELS +j] *= Ueber->LMF_input->tdcresolution;
					}
				}
		
				if (Ueber->LMF_input->common_mode == 1) {
					if (!Ueber->reverse_time_direction_in_event) {
						for (__int32 ch = 0;ch< max_chan;ch++) {
							if (cnt[ch] > 1) {
								if (tdc_ns[ch*NUM_CHANNELS+0] < tdc_ns[ch*NUM_CHANNELS + 1]) {Ueber->reverse_time_direction_in_event =  1; break;}
								if (tdc_ns[ch*NUM_CHANNELS+0] > tdc_ns[ch*NUM_CHANNELS + 1]) {Ueber->reverse_time_direction_in_event = -1; break;}
							}
						}
					} else {
						for (__int32 ch = 0;ch< max_chan;ch++) {
							__int32 n = cnt[ch] >> 1;
							for (__int32 i=0;i<n;i++) {
								double temp = tdc_ns[ch*NUM_CHANNELS + i];
								tdc_ns[ch*NUM_CHANNELS + i] = tdc_ns[ch*NUM_CHANNELS + cnt[ch]-1-i];
								tdc_ns[ch*NUM_CHANNELS + cnt[ch]-1-i] = temp;
							}
						}
					}
				}
			
				timestamp = Ueber->LMF_input->GetDoubleTimeStamp();
			}
		

			// ADC - ADC - ADC - ADC - ADC - ADC - ADC - ADC -ADC - ADC -ADC - ADC
			//OutputDebugString("\nLMF_read_event_loop1");

			ADC_meta * adc = (ADC_meta*)Ueber->adc_meta;
					
			if(adc && (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4 || Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8 || Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT)){
				if (!adc->GetTDCandCNTArray(tdc_ns,cnt))
					break;
			}

			//if (Ueber->parameter[850] > 0.5) {
			//	std::cout << "\nLMA file mode selected..";
			//}
			// ADC - ADC - ADC - ADC - ADC - ADC - ADC - ADC -ADC - ADC -ADC - ADC

			// start of fADC8-stuff
			if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8) {
				printf("please check the source code. This part is commented out.");
				Ueber->stop_reading_files = true;
			}

		} // end if (Ueber->LMF_input)
	
		if (Ueber->HDF5_input) {
			ADC_meta * adc = (ADC_meta*)Ueber->adc_meta;			
			if (!adc->GetTDCandCNTArray(tdc_ns,cnt))
				break;
			if (Ueber->HDF5_input->error_flag) Ueber->stop_reading_files;
			if (Ueber->HDF5_input->eof) break;
		}

		// demix:
		if (entmixen) {
			double temp_tdc[NUM_IONS];
			__int32 n = cnt[merged_channel];
			for (__int32 i = 0;i < n;++i) temp_tdc[i] = tdc_ns[merged_channel*NUM_CHANNELS+i];
			cnt[new_demerged_channel] = 0;
			cnt[merged_channel]   = 0;

			for (__int32 i = 0;i<n;++i) {
				double a = temp_tdc[i];
				if (a > timelimit) {
					tdc_ns[new_demerged_channel*NUM_CHANNELS + cnt[new_demerged_channel]] = a;
					cnt[new_demerged_channel] = cnt[new_demerged_channel] + 1;
				} else {
					tdc_ns[merged_channel*NUM_CHANNELS + cnt[merged_channel]] = a;
					cnt[merged_channel] = cnt[merged_channel] + 1;
				}
			}
		}


		__int32 twin_counter = 0;
		for (__int32 i=0;i<NUM_CHANNELS;++i) {
			if(i!=Ueber->scan_channel) {
				double old_value = tdc_ns[i*NUM_CHANNELS + 0];
				for (__int32 j=1;j<cnt[i];++j) {	
					if (tdc_ns[i*NUM_CHANNELS + j] == old_value) {
						twin_counter++;
						if (twin_counter > 2) {cnt[i]=0; break;}
					}
					old_value = tdc_ns[i*NUM_CHANNELS + j];
				}
			}
		}



		// Anti-Moire :
		if (parameter[15] > 0.5 && Ueber->LMF_input) {
			double resol = Ueber->LMF_input->tdcresolution;
			for (__int32 i=0;i<NUM_CHANNELS;++i) {
				for (__int32 j=0;j<cnt[i];++j) {	
					if(i!=Ueber->scan_channel) {
						tdc_ns[i*NUM_CHANNELS + j] += (double(rand())/RAND_MAX-0.5) * resol;
					}
				}
			}
		}

		if (Ueber->stop_reading_files) {
			return true;
		}
		
		if ((++rate_counter) == rate_counter_limit) {
			rate_counter = 0;

			__int32 c = my_kbhit();
			if (c) 	{
				bool do_break = false;
				if (c=='q') {printf("\nq was pressed -> skipping this file\npress Q to skip all files.\n"); do_break = true;}
				if (c=='Q') {printf("\nQ was pressed -> skipping all files.\n"); Ueber->stop_reading_files = true; do_break = true;}

				if (c=='p') {
					printf("Pause. Hit any key to continue.");
					__int32 i=0;
					while (true) {
						i++;
						if (i > 4*4) {
							i=0;
							delete_file("are_you_alive.txt");
						}
						char c = my_kbhit();
						if (!c) {gSystem->Sleep(250); continue;}
						printf("\r");
						if (c=='q') {printf("\nq was pressed -> skipping this file\npress Q to skip all files.\n"); do_break = true;}
						if (c=='Q') {printf("\nQ was pressed -> skipping all files.\n"); Ueber->stop_reading_files = true; do_break = true;}
						break;
					}
				}

				// User invoke IPA
				if(Ueber->ipa_enable) {
					if(c=='i' ||c=='I') {
//						Ueber->IPA->main_loop();
						Ueber->stop_reading_files = true; do_break = true;
					}
				}

				if (do_break) {
					if (Ueber->LMF_input) {
						printf("event %I64i of %I64i\n",Ueber->LMF_input->int64_number_of_read_events, Ueber->LMF_input->int64_Numberofevents);
						if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",Ueber->LMF_input->int64_number_of_read_events, Ueber->LMF_input->int64_Numberofevents); fflush(Ueber->logfile_handle);}
					}
					if (Ueber->HDF5_input) {
						printf("event %I64i of %I64i\n",Ueber->HDF5_input->int64_number_of_read_events, Ueber->HDF5_input->int64_Numberofevents);
						if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",Ueber->HDF5_input->int64_number_of_read_events, Ueber->HDF5_input->int64_Numberofevents); fflush(Ueber->logfile_handle);}
					}
					if (Ueber->stop_reading_files) return true;
					break;
				}
			}

			// Invoke IPA as number of requested events were read...
			if(Ueber->ipa_enable) {
				if(Ueber->IPA->events_in_memory() >= Ueber->ipa_num_events) {
//					Ueber->IPA->main_loop();
					Ueber->stop_reading_files = true; return true;
				}
			}

			double new_time = get_time_s();
			static double print_time = 0;
			if (new_time > print_time + 1.) {
				check_for_pause_file(Ueber);
				print_time = new_time;
				__int32 stars_in_progress_bar = 0;
				if (Ueber->LMF_input) {
					if (Ueber->LMF_input->int64_Numberofevents > 0) {
						stars_in_progress_bar = __int32(52. *Ueber->LMF_input->int64_number_of_read_events / Ueber->LMF_input->int64_Numberofevents);
					}
				}
				if (Ueber->HDF5_input) {
					if (Ueber->HDF5_input->int64_Numberofevents > 0) {
						stars_in_progress_bar = __int32(52. *Ueber->HDF5_input->int64_number_of_read_events / Ueber->HDF5_input->int64_Numberofevents);
					}
				}
				printf("\r");
				for (__int32 i=0;i<stars_in_progress_bar;++i) printf("*");
				printf("| %i/s   ",__int32(rate_counter_limit/(new_time - Ueber->old_time)));
			}
#ifndef _DEBUG
			rate_counter_limit = __int32(0.25*rate_counter_limit / (new_time - Ueber->old_time));
			if (rate_counter_limit > 10000) rate_counter_limit = 10000;
			if (rate_counter_limit < 100) rate_counter_limit = 100;
#endif
			Ueber->old_time = new_time;
			
			if (new_time > root_file_flush_time + 30) {
				root_file_flush_time = new_time;
				if (Ueber->outputRootFile_handle) {Ueber->outputRootFile_handle->Flush();}
			}
		}

if (event_is_valid) {
	sort_and_write_NTuple(Ueber->eventcounter,Ueber,parameter, timestamp); // process TDC Data in lmf-File
	++Ueber->eventcounter;
}

		if (parameter[58] > 0 && Ueber->eventcounter > (__int64)(parameter[58])) return true;

		__int32 temp_tdc_ch[NUM_CHANNELS*NUM_IONS];
		unsigned __int32 uicnt[NUM_CHANNELS];

		if (event_is_valid && write_sorted_file) { // write sorted LMF-file

//			Ueber->write_this_event_into_sorted_LMF = true; // XXX

			if(Ueber->write_this_event_into_sorted_LMF)
			{
				if (!Ueber->LMF_output) {
					if (!Ueber->LMF_input) {
						write_sorted_file = false;
					} else {
						Ueber->LMF_output = new LMF_IO(NUM_CHANNELS,NUM_IONS);
						Ueber->LMF_input->Clone(Ueber->LMF_output);
						if (Ueber->LMF_output) {
							char name[552];
							if (strlen(Ueber->LMF_output_file_name_without_extension) > 0) {
								sprintf(name,"%s_%03d.lmf",Ueber->LMF_output_file_name_without_extension,Ueber->output_LMF_file_index);
							} else {
								sprintf(name,"%s_%03d.lmf",Ueber->input_data_file_name, Ueber->output_LMF_file_index);
							}

							if (!Ueber->LMF_output->OpenOutputLMF(name)) {
								printf("Error: Could not open output LMF for writing.\n");
								if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"Error: Could not open output LMF for writing.\n"); fflush(Ueber->logfile_handle);}
								Ueber->stop_reading_files = true;
								break;
							}

							++Ueber->output_LMF_file_index;
						} else write_sorted_file = false;
					}
				}


				// fADC8:
				if (Ueber->LMF_output) {
					if (Ueber->LMF_output->DAQ_ID_output == DAQ_ID_FADC8) {
						// Ueber->LMF_input->input_lmf->seek(current_file_position);
						while (true) {
							fADC8_signal_info_struct signal_info;
							bool bEnd_of_group_detected;
							Ueber->LMF_input->ReadNextfADC8Signal(signal_info, bEnd_of_group_detected, ui32buffer, ui32_BUFFERSIZE);
						
							Ueber->LMF_output->output_lmf->write((char*)(&signal_info), sizeof(signal_info));
							if (signal_info.signal_type != 8) {
								Ueber->LMF_output->output_lmf->write((char*)ui32buffer, signal_info.signallength_including_header_in_32bit_words*4);
							}

							if (bEnd_of_group_detected) {
								*Ueber->LMF_output->output_lmf << Ueber->LMF_input->fADC8.GroupEndMarker;
								Ueber->LMF_output->uint64_number_of_written_events++;
								break;
							}
							__int32 int_dummy = 0;
							*Ueber->LMF_output->output_lmf << int_dummy;
						}
				
						if (Ueber->LMF_output->uint64_number_of_written_events >= (unsigned __int64)max_events_to_write_into_file && max_events_to_write_into_file > 0)
						{
							Ueber->LMF_output->CloseOutputLMF();
							Ueber->LMF_output->CloseInputLMF();
							delete Ueber->LMF_output;
							Ueber->LMF_output = 0;
						}
					}



					if ((Ueber->LMF_output->DAQ_ID_output != DAQ_ID_FADC8) && (Ueber->LMF_output->DAQ_ID_output != DAQ_ID_FADC4)) {
						if (parameter[100] > 0.5) {
							Ueber->current_proj_Sorter->shift_sums(-1,Ueber->proj->offset_timesum_U,Ueber->proj->offset_timesum_V,Ueber->proj->offset_timesum_W);
							Ueber->current_proj_Sorter->shift_position_origin(-1,Ueber->proj->center_X,Ueber->proj->center_Y);	// shift all layers so that the position picture is centered around X=zero,Y=zero
							if (Ueber->current_proj_Sorter->use_HEX) Ueber->current_proj_Sorter->shift_layer_w(-1,Ueber->proj->hex_offset_w);						// shift layer w so that all 2 center lines of the layers meet in one point
						}
						if (parameter[200] > 0.5) {
							Ueber->current_rec_Sorter->shift_sums(-1,Ueber->rec->offset_timesum_U,Ueber->rec->offset_timesum_V,Ueber->rec->offset_timesum_W);
							Ueber->current_rec_Sorter->shift_position_origin(-1,Ueber->rec->center_X,Ueber->rec->center_Y);	// shift all layers so that the position picture is centered around X=zero,Y=zero
							if (Ueber->current_rec_Sorter->use_HEX) Ueber->current_rec_Sorter->shift_layer_w(-1,Ueber->rec->hex_offset_w);						// shift layer w so that all 2 center lines of the layers meet in one point
						}
						if (parameter[300] > 0.5) {
							Ueber->current_elec_Sorter->shift_sums(-1,Ueber->elec->offset_timesum_U,Ueber->elec->offset_timesum_V,Ueber->elec->offset_timesum_W);
							Ueber->current_elec_Sorter->shift_position_origin(-1,Ueber->elec->center_X,Ueber->elec->center_Y);	// shift all layers so that the position picture is centered around X=zero,Y=zero
							if (Ueber->current_elec_Sorter->use_HEX) Ueber->current_elec_Sorter->shift_layer_w(-1,Ueber->elec->hex_offset_w);						// shift layer w so that all 2 center lines of the layers meet in one point
						}

						for (__int32 i=0;i<NUM_CHANNELS;++i) 
						{
							uicnt[i] = cnt[i];
							for (__int32 j=0;j<cnt[i];++j) 
							{
								double a = 0.;
								if (Ueber->LMF_input)  a = tdc_ns[i*NUM_IONS+j]/Ueber->LMF_input->tdcresolution;
								if (Ueber->HDF5_input) a = tdc_ns[i*NUM_IONS+j]/0.025; // don't know if this is right... (Achim)

								if (a >= 0.) temp_tdc_ch[i*NUM_IONS+j] = __int32(a+0.001);
								if (a < 0. ) temp_tdc_ch[i*NUM_IONS+j] = __int32(a-0.001);
							}
						}

						Ueber->LMF_output->WriteTDCData(timestamp,uicnt,temp_tdc_ch);
						if (Ueber->LMF_output->uint64_number_of_written_events >= (unsigned __int64)max_events_to_write_into_file && max_events_to_write_into_file > 0)
						{
							Ueber->LMF_output->CloseOutputLMF();
							Ueber->LMF_output->CloseInputLMF();
							delete Ueber->LMF_output;
							Ueber->LMF_output = 0;
						}
					}
				}
			}
		}
	} // end while ...

	return true;
}




















void ProcessLMAFile(char *input_LMA_file_name, Ueberstruct * Ueber, double tds_ns[][NUM_IONS])
{
	//lma2root *lma = new lma2root(input_LMA_file_name);
	//lma->ReadLMAHeader();
}















////////////////////////////////////////////////////////////////////
void ProcessInputDataFile(char * _input_data_file_name, Ueberstruct * Ueber)
////////////////////////////////////////////////////////////////////
{
	sprintf(Ueber->input_data_file_name, _input_data_file_name);
	Ueber->reverse_time_direction_in_event = 0;
	if (Ueber->LMF_output_file_name_without_extension[0] == 0) Ueber->output_LMF_file_index = 0;
	double * parameter = Ueber->parameter;
	__int32 * cnt = Ueber->cnt;
	if (Ueber->stop_reading_files) return;

	bool write_sorted_file = ((parameter[49] > 0.5) ? true:false);

	if (Ueber->LMF_input)  {delete Ueber->LMF_input;  Ueber->LMF_input  = 0;}
	if (Ueber->HDF5_input) {delete Ueber->HDF5_input; Ueber->HDF5_input = 0;}

	if (!Ueber->input_is_HDF5) {
		Ueber->LMF_input = new LMF_IO(NUM_CHANNELS,NUM_IONS);
	}

	if (Ueber->input_is_HDF5) {
		Ueber->HDF5_input = new HDF5_IO(
			int(parameter[830]+0.1),				// trigger_channel
			parameter[831],							// group_range_start
			parameter[832],							// group_range_end
			// parameter[833],							// trigger dead time
			int(parameter[834]+0.001),				// time slice width in ns
			int(parameter[835]+0.001),				// ADC binsize in ps
			parameter[836],							// full scale in mV
			(parameter[837] > 0.5) ? true:false);	// use_grouping = 1 = true, use time slicing = 0 = false
	}

	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4) {
			Ueber->LMF_input->CloseInputLMF();
			delete Ueber->LMF_input;	Ueber->LMF_input = 0;
			//extract_mean_pulse_fADC4(input_data_file_name, Ueber, tdc_ns); // xxx Robert
			Ueber->LMF_input = new LMF_IO(NUM_CHANNELS,NUM_IONS);
		}
	}
	Ueber->fast_mode = parameter[11] > 0.5 ? 1:0;

	if (Ueber->LMF_input) {
		Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions = __int32(parameter[9] + 0.1);
		if (parameter[9] < -0.1) Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions--;
	}

	for (__int32 det = 0; det < 3;++det) {
		Det_struct * detector = 0;
		sort_class * sorter = 0;
		__int32 parameter_offset = 0;
		if (det == 0) {detector = Ueber->proj; sorter = Ueber->proj_Sorter; parameter_offset = 100;}
		if (det == 1) {detector = Ueber->rec;  sorter = Ueber->rec_Sorter;  parameter_offset = 200;}
		if (det == 2) {detector = Ueber->elec; sorter = Ueber->elec_Sorter; parameter_offset = 300;}
		if (detector) {
			if (detector->sumu_watchdog) detector->sumu_watchdog->reset();
			if (detector->sumv_watchdog) detector->sumv_watchdog->reset();
			if (detector->sumw_watchdog) detector->sumw_watchdog->reset();
			if (!fill_parameters_into_sorters(sorter, detector, parameter_offset)) {
				Ueber->stop_reading_files = true;
				if (Ueber->LMF_input) {
					Ueber->LMF_input->CloseInputLMF();
					delete Ueber->LMF_input;	Ueber->LMF_input = 0;
				}
				if (Ueber->HDF5_input) {
					Ueber->HDF5_input->CloseInputHDF5();
					delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
				}
				return;
			}
		}
	}

	if ( parameter[58] > 0 ) printf("\n\n%G events left (out of %G)\n",parameter[58] - Ueber->eventcounter,parameter[58]);

	printf("\n\nInput file is: %s\n",Ueber->input_data_file_name);
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n\nreading input file %s\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
	if(!File_is_readable(Ueber->input_data_file_name)) {
		printf("\n \n %s could not be opened!\n\n",Ueber->input_data_file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
		if (Ueber->LMF_input) {
			Ueber->LMF_input->CloseInputLMF();
			delete Ueber->LMF_input;	Ueber->LMF_input = 0;
		}
		if (Ueber->HDF5_input) {
			Ueber->HDF5_input->CloseInputHDF5();
			delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
		}
		return;
	}

	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID != DAQ_ID_AGAT && !Ueber->LMF_input->OpenInputLMF(Ueber->input_data_file_name)) {	// AGAT files are handled separate..
			printf("\n \n %s could not be opened!\n\n",Ueber->input_data_file_name);
			if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
			Ueber->stop_reading_files = true;
			if (Ueber->LMF_input) {
				Ueber->LMF_input->CloseInputLMF();
				delete Ueber->LMF_input;	Ueber->LMF_input = 0;
			}
			return;
		}
	}

	if (Ueber->HDF5_input) {
		if (!Ueber->HDF5_input->OpenInputHDF5(Ueber->input_data_file_name)) {
			printf("\n \n %s could not be opened!\n\n",Ueber->input_data_file_name);
			if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_data_file_name); fflush(Ueber->logfile_handle);}
			Ueber->stop_reading_files = true;
			if (Ueber->HDF5_input) {
				Ueber->HDF5_input->CloseInputHDF5();
				delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
			}
			return;
		}
	}


	if (write_sorted_file) {
		char name[552];
		if (strlen(Ueber->LMF_output_file_name_without_extension) > 0) {
			sprintf(name,"%s_%03d.lmf",Ueber->LMF_output_file_name_without_extension,Ueber->output_LMF_file_index);
		} else {
			sprintf(name,"%s_%03d.lmf",Ueber->input_data_file_name,Ueber->output_LMF_file_index);
		}
		printf("\nwriting new LMF-file\n\"%s\"\n",name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\nwriting new LMF-file\n\"%s\"\n",name); fflush(Ueber->logfile_handle);}
	}

	__int32 number_of_hits = 0;
	if (Ueber->LMF_input) {
		Ueber->number_of_channels	= ((__int32(parameter[54]+0.1) == 0) ? (Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2) : __int32(parameter[54]+0.1));
		number_of_hits				= ((__int32(parameter[55]+0.1) == 0) ? (Ueber->LMF_input->max_number_of_hits) : __int32(parameter[55]+0.1));
	
		if (Ueber->proj_Sorter) Ueber->proj_Sorter->common_start_mode = Ueber->LMF_input->common_mode == 0 ? true : false;
		if (Ueber->rec_Sorter)  Ueber->rec_Sorter-> common_start_mode = Ueber->LMF_input->common_mode == 0 ? true : false;
		if (Ueber->elec_Sorter) Ueber->elec_Sorter->common_start_mode = Ueber->LMF_input->common_mode == 0 ? true : false;

		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8HP) {
			if (Ueber->proj_Sorter) Ueber->proj_Sorter->common_start_mode = true;
			if (Ueber->rec_Sorter)  Ueber->rec_Sorter-> common_start_mode = true;
			if (Ueber->elec_Sorter) Ueber->elec_Sorter->common_start_mode = true;
		}

		if (parameter[10] > 0.5) { // xxx
			if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4) {
				parameter[10] = 0.;
				Red();
				printf("parameter[10]: Multicore does not work with ADC-data.\n");
				White();
			}
		}
	}

	if (Ueber->HDF5_input) {
		Ueber->number_of_channels	= ((__int32(parameter[54]+0.1) == 0) ? Ueber->HDF5_input->number_of_channels : __int32(parameter[54]+0.1));
		number_of_hits				= ((__int32(parameter[55]+0.1) == 0) ? Ueber->HDF5_input->max_number_of_hits : __int32(parameter[55]+0.1));
	
		if (Ueber->proj_Sorter) Ueber->proj_Sorter->common_start_mode = true;
		if (Ueber->rec_Sorter)  Ueber->rec_Sorter-> common_start_mode = true;
		if (Ueber->elec_Sorter) Ueber->elec_Sorter->common_start_mode = true;

		if (parameter[10] > 0.5) { // xxx
			parameter[10] = 0.;
			Red();
			printf("parameter[10]: Multicore does not work with ADC-data.\n");
			White();
		}
	}


	for (__int32 detector = 0; detector < 3; ++detector) {
		__int32 para_offset;
		sort_class * sorter		= 0;
		Det_struct * det_struct = 0;
		if (detector == 0) {para_offset = 100; sorter = Ueber->proj_Sorter; det_struct = Ueber->proj;}
		if (detector == 1) {para_offset = 200; sorter = Ueber->rec_Sorter;  det_struct = Ueber->rec;}
		if (detector == 2) {para_offset = 300; sorter = Ueber->elec_Sorter; det_struct = Ueber->elec;}
		if (sorter) {
			if (Ueber->LMF_input) sorter->TDC_resolution_ns = Ueber->LMF_input->tdcresolution;
			if (Ueber->HDF5_input) {sorter->TDC_resolution_ns = 0.025;}
			if (!sorter->initialization_successful) {
				__int32 err = sorter->init_after_setting_parameters();
				if (err) {
					char error_message[200];
					sorter->get_error_text(err,200,error_message);
					if (detector == 0) printf("\nerror on detector \"proj\":\n%s\n",error_message);
					if (detector == 1) printf("\nerror on detector \"rec\":\n%s\n", error_message);
					if (detector == 2) printf("\nerror on detector \"elec\":\n%s\n",error_message);
					Ueber->stop_reading_files = true;
				}
			}
		}
	}

	///////////////////////////////////////////////////////
	//    multi-core support
	///////////////////////////////////////////////////////

		#ifndef dont_use_MFC
	if (parameter[10] > 0.5) { // use multi-core processing?

		parameter[10] = 0.; printf("multicore does not work at the moment\n"); Ueber->stop_reading_files = true;// xxx

		if (!Ueber->parallel_sorter) {
			if (atoi(getenv("NUMBER_OF_PROCESSORS")) > 1) {
				__int32 number_of_processors = atoi(getenv("NUMBER_OF_PROCESSORS"))-1;
				if ((parameter[10] > 1.5) && parameter[10] < number_of_processors) number_of_processors = __int32(parameter[10]+0.1);

				// number_of_processors = 1; // xxx remove this line (used only for debugging)

					static parallel_sort_class * parallel = new parallel_sort_class(number_of_processors,NUM_CHANNELS,NUM_IONS,Ueber->proj_Sorter,Ueber->proj->use_reconstruction,Ueber->rec_Sorter,Ueber->rec->use_reconstruction,Ueber->elec_Sorter,Ueber->elec->use_reconstruction);
					Ueber->parallel_sorter = (void*)parallel;

					if (Ueber->parallel_sorter == 0) {
						printf("could not initialize multi core support.\n");
						Ueber->stop_reading_files = true;
						if (Ueber->LMF_input) {
							Ueber->LMF_input->CloseInputLMF();
							delete Ueber->LMF_input;	Ueber->LMF_input = 0;
						}
						if (Ueber->HDF5_input) {
							Ueber->HDF5_input->CloseInputHDF5();
							delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
						}
						return;
					}
					parallel->start_sort_threads();
			} else parameter[10] = 0.;
		}
	}
#endif

	if (parameter[2] > 0.5) {
		if (Ueber->LMF_input) {
			printf("\nInput file header: LMF-Type :%i",Ueber->LMF_input->data_format_in_userheader);
			if (Ueber->LMF_input->data_format_in_userheader==2)  printf(" = 16bit unsigned integer ");
			if (Ueber->LMF_input->data_format_in_userheader==6) {
				printf(" = 24bit unsigned integer not supported");
				Ueber->stop_reading_files = true;
				if (Ueber->LMF_input) {
					Ueber->LMF_input->CloseInputLMF();
					delete Ueber->LMF_input;	Ueber->LMF_input = 0;
				}
				return;
			}
			if (Ueber->LMF_input->data_format_in_userheader==3)  printf(" = 32bit signed integer ");
			if (Ueber->LMF_input->data_format_in_userheader==-1) printf(" = 32bit signed integer variable event size");
			#ifndef LINUX
					//printf("\n            Coordinates     :%11ld ",Ueber->LMF_input->Numberofcoordinates);     
					printf("\n            channels        :%11ld ",Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2);
					printf("\n            hits per channel:%11ld ",Ueber->LMF_input->max_number_of_hits);     
					printf("\n            HeaderSize      :%11ld ",Ueber->LMF_input->Headersize);
					//printf("\n            UserHeaderSize  :%11ld ",Ueber->LMF_input->User_header_size);
					printf("\n            Events in file  :%11I64i ",Ueber->LMF_input->int64_Numberofevents);
			#endif
			printf("\n            LMVersionString : %s    ",Ueber->LMF_input->Versionstring.c_str());
		}
	}

	printf("\n0%%           25%%          50%%          75%%        100%%\n");
	printf("|------------|------------|------------|-----------|\n");


	ADC_meta * adc = 0;
	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4 || Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8 || Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT) { // is this a LMF with ADC-data?
			// ADC stuff
			if (fabs(parameter[819]-1234567) > 0.1) {
				Red();
				printf("error reading ADC-data: You must have the line 'execute \"ADC\\ADC_parameters.txt\"' in the config-file.\n");
				White();
				if (Ueber->LMF_input) {Ueber->LMF_input->CloseInputLMF(); delete Ueber->LMF_input;	Ueber->LMF_input = 0;}
				if (adc) {adc->deleteInstance(); adc = 0;}
				return;
			} else {
				if(!Ueber->adc_meta) Ueber->adc_meta = ADC_meta::instance(Ueber, NUM_CHANNELS, MAX_NBR_PULSES, MAX_NBR_PEAKS);//(void*)new ADC_analysis(Ueber);
				adc = (ADC_meta*)Ueber->adc_meta;
				if (parameter[950] > 0.5 && adc) {
				//	adc->inspectDPA(parameter[951]);
				}
			}
		}
	}

	if (Ueber->HDF5_input) {
		// ADC stuff	 
		if(!Ueber->adc_meta) Ueber->adc_meta = ADC_meta::instance(Ueber, NUM_CHANNELS, MAX_NBR_PULSES, MAX_NBR_PEAKS);//(void*)new ADC_analysis(Ueber);
		adc = (ADC_meta*)Ueber->adc_meta;
		if (parameter[950] > 0.5 && adc) {
		//	adc->inspectDPA(parameter[951]);
		}
	}


	OutputDebugString("\ngo into read_event_loop1");
	LMF_read_event_loop(Ueber);
	
	if (adc) {adc->deleteInstance(); adc = 0; Ueber->adc_meta = 0;}

	printf("\nworst sum misalignments in this file:\n");
	printf("proj: u=%lf \tv=%lf \tw=%lf\n",Ueber->proj->sumu_watchdog->maximum_total_correction_up_to_now,Ueber->proj->sumv_watchdog->maximum_total_correction_up_to_now,Ueber->proj->sumw_watchdog->maximum_total_correction_up_to_now);
	printf("rec:  u=%lf \tv=%lf \tw=%lf\n",Ueber->rec->sumu_watchdog->maximum_total_correction_up_to_now ,Ueber->rec->sumv_watchdog->maximum_total_correction_up_to_now, Ueber->rec->sumw_watchdog->maximum_total_correction_up_to_now);
	printf("elec: u=%lf \tv=%lf \tw=%lf\n\n",Ueber->elec->sumu_watchdog->maximum_total_correction_up_to_now,Ueber->elec->sumv_watchdog->maximum_total_correction_up_to_now,Ueber->elec->sumw_watchdog->maximum_total_correction_up_to_now);

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\nworst sum misalignments in this file:\n"); fflush(Ueber->logfile_handle);}
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"proj: u=%lf \tv=%lf \tw=%lf\n",Ueber->proj->sumu_watchdog->maximum_total_correction_up_to_now,Ueber->proj->sumv_watchdog->maximum_total_correction_up_to_now,Ueber->proj->sumw_watchdog->maximum_total_correction_up_to_now); fflush(Ueber->logfile_handle);}
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"rec:  u=%lf \tv=%lf \tw=%lf\n",Ueber->rec->sumu_watchdog->maximum_total_correction_up_to_now ,Ueber->rec->sumv_watchdog->maximum_total_correction_up_to_now, Ueber->rec->sumw_watchdog->maximum_total_correction_up_to_now); fflush(Ueber->logfile_handle);}
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"elec: u=%lf \tv=%lf \tw=%lf\n\n",Ueber->elec->sumu_watchdog->maximum_total_correction_up_to_now,Ueber->elec->sumv_watchdog->maximum_total_correction_up_to_now,Ueber->elec->sumw_watchdog->maximum_total_correction_up_to_now); fflush(Ueber->logfile_handle);}

	if (Ueber->LMF_input) {
		Ueber->LMF_input->CloseInputLMF();
		delete Ueber->LMF_input;	Ueber->LMF_input = 0;
	}
	if (Ueber->HDF5_input) {
		Ueber->HDF5_input->CloseInputHDF5();
		delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
	}

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"finished reading LMF-file.\n"); fflush(Ueber->logfile_handle);}

	
/*	if (Ueber->bot_mode) {
		char command[300]; command[0] = 0;
		sprintf(command,"del %s", input_data_file_name);
		system(command);
	}*/
}



































/////////////////////////////////////////////////////////////////////////////
void split_pathname(char * input, char * path, char * name)
/////////////////////////////////////////////////////////////////////////////
{
	if (!input) return;
	if (path) path[0] = 0;
	if (name) name[0] = 0;
	sprintf(path,input);
	while (true) {
		if (!strlen(path)) break;
		char c = path[strlen(path)-1];
		if (c != '\\' && c != '/') {
			path[strlen(path)-1] = 0;
			continue;
		}
		c = path[strlen(path)-1];
		if (c == '\\' || c == '/') path[strlen(path)-1] = 0;
		break;
	}
	if (name) {
		for (unsigned __int32 i=strlen(path);i<strlen(input);i++) {
			name[i-strlen(path)] = input[i];
			name[i+1-strlen(path)] = 0;
		}
	}
	if (path[strlen(path)-1] == '.') path[strlen(path)-1] = 0;
}






void get_name_without_path(char * in_string, char * out_string) {
	__int32 slash_index = -1;
	for (__int32 i=0;i<__int32(strlen(in_string));i++) {
		if (in_string[i] == '/' || in_string[i] == '\\') slash_index = i;
	}
	for (__int32 i=slash_index+1;i<=__int32(strlen(in_string));i++) { // <= is ok here
		out_string[i-(slash_index+1)] = in_string[i];
		out_string[i-(slash_index+1)+1] = 0;
	}
}








void get_path_without_name(char * in_string, char * out_string) {
	__int32 slash_index = -1;
	for (__int32 i=0;i<__int32(strlen(in_string));i++) {
		if (in_string[i] == '/' || in_string[i] == '\\') slash_index = i;
	}
	for (__int32 i=0;i<slash_index;i++) {
		out_string[i] = in_string[i];
	}
	if (slash_index < 0) slash_index = 0;
	out_string[slash_index] = 0;
}




bool get_next_new_data_file_name(char * new_name, char * main_name, Ueberstruct * Ueber)
{
#ifdef LINUX
	printf("wildcards are not yet allowed in the Linux version of LMF2root\n");
	return false;
#else

	if (!new_name) return false;
	if (new_name) new_name[0] = 0;
	if (!main_name) return false;
	if (!Ueber) return false;

	WIN32_FIND_DATA FindFileData;
	HANDLE hFind = FindFirstFile(main_name, &FindFileData);
	if (hFind == INVALID_HANDLE_VALUE) return false;

	bool a_new_file_was_found = true;
	DWORD last_status = 1;

	char new_short_name[300];
	char new_name_[250]; new_name[0]  = 0;
	char path_name[250]; path_name[0] = 0;
	get_path_without_name(main_name,path_name);

	while (true) {
		if (!last_status) {
			if (GetLastError() == ERROR_NO_MORE_FILES) break;
		}

		if (FindFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			last_status = FindNextFile(hFind, &FindFileData);
			continue;
		}

		sprintf(new_name_, FindFileData.cFileName);
		get_name_without_path(new_name_, new_short_name);

		a_new_file_was_found = true;
		for (__int32 i=0;i<5000;i++) {
			if (!Ueber->processed_data_file_names[i]) break;
			char temp_string[300];
			get_name_without_path(Ueber->processed_data_file_names[i], temp_string);
			if (!strcmp(temp_string,new_short_name)) {a_new_file_was_found = false; break;} // identical
		}
		if (a_new_file_was_found) {
			if (strstr(new_short_name,"emplate")) a_new_file_was_found = false;
		}
		if (a_new_file_was_found) {
			
			if (strlen(path_name) > 0) {
				#ifdef LINUX
					sprintf(new_name,"%s/%s",path_name,new_name_);
				#else
					sprintf(new_name,"%s\\%s",path_name,new_name_);
				#endif
			} else {
				#ifdef LINUX
						sprintf(new_name,"%s",new_name_);
				#else
						sprintf(new_name,"%s",new_name_);
				#endif
			}
			FILE * file = fopen(new_name,"a"); // check that the file is not still open
			if (!file) a_new_file_was_found = false;
			if (file) {
				if (ferror(file)) a_new_file_was_found = false;
				fclose(file);
			}
		}
		if (a_new_file_was_found) break;
		last_status = FindNextFile(hFind, &FindFileData);
	}

	FindClose(hFind);
	new_name[0] = 0;
	if (!a_new_file_was_found) return false;

	if (strlen(path_name) > 0) {
		#ifdef LINUX
			sprintf(new_name,"%s/%s",path_name,new_name_);
		#else
			sprintf(new_name,"%s\\%s",path_name,new_name_);
		#endif
	} else {
		#ifdef LINUX
				sprintf(new_name,"%s",new_name_);
		#else
				sprintf(new_name,"%s",new_name_);
		#endif
	}

	for (__int32 i=0;i<5001;i++) {
		if (i == 5000) {
			new_name[0] = 0;
			Ueber->stop_reading_files = true;
			return false;
		}
		if (!Ueber->processed_data_file_names[i]) {
			Ueber->processed_data_file_names[i] = new char[strlen(new_name_)+1];
			sprintf(Ueber->processed_data_file_names[i],new_name_);
			break;
		}
	}
	
	return true;
#endif
}






















void get_current_directory(char * buffer, __int32 max_len)
{
	#ifdef LINUX
		getwd(buffer);
		__int32 len = strlen(buffer);
		if (buffer[len-1] != '/') {
			buffer[len] = '/';
			buffer[len+1] = 0;
			return;
		}
		return;
	#else
		GetCurrentDirectory(max_len, buffer);
		__int32 len = strlen(buffer);
		if (buffer[len-1] != '\\') {
			buffer[len] = '\\';
			buffer[len+1] = 0;
			return;
		}
	#endif
}









/////////////////////////////////////////////////////////////////////////////
int lmf2root_main(__int32 argc, char* argv[])
/////////////////////////////////////////////////////////////////////////////
{
	double * parameter = 0;
	__int32	master_cnt[NUM_CHANNELS];
	double	master_tdc_ns[NUM_CHANNELS][NUM_IONS];

	#ifdef _DEBUG
		Red(true);
//		printf("\n\n\nLMF2ROOT: USING MEGA SLOW DEBUGGING MODE !!!");	_CrtSetDbgFlag(_CrtSetDbgFlag( _CRTDBG_REPORT_FLAG ) | _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF | _CRTDBG_CHECK_ALWAYS_DF);
		printf("\n***********************\n    SLOW DEBUG VERSION !\n***********************\n\n");
		White(false);
	#endif

	version_number_class sorter_version = sort_class::get_version_number();
	printf("LMF2ROOT   version %.2f (by Achim Czasch) using\n",LMF2ROOTVERSION); // XXX
	printf("resort64c  version %i.%i.%i.%i and\n",sorter_version.A,sorter_version.B,sorter_version.C,sorter_version.D);
	printf("LMF_IO     version %i\n",LMF_IO::GetVersionNumber());
	printf("ColAHelL   version %.2f (by T. Jahnke)\n",COLAHELLVERSION);
	printf("IPA        version %.2f (by T. Jahnke)\n",IPAVERSION);
	printf("ADC support version %.2f (by T. Bauer / R. Wallauer)\n\n",ADCVERSION);


	Red(true);
	printf("Don't panic! Everything will be fine.\n");
	White(false);

	FILE * ascii_file = 0;

	//////////////////////////////////////////
	// Achim's root stuff:
	// start the Root framework:
	char * root_argv[3];
	char argv1[50],argv2[50],argv3[50];
	__int32 root_argc = 3;
	sprintf(argv1,"root_App");	sprintf(argv2,"-l");  sprintf(argv3,"-n");
	root_argv[0] = argv1;	root_argv[1] = argv2;	root_argv[2] = argv3;
	TApplication theRootApp("root_App", &root_argc, root_argv);
	//------------------------------------------

	delete_file("are_you_alive.txt");
	delete_file("I_am_busy.txt");
	if (!FileExists("I_am_alive.txt")) {
		FILE * fo = fopen("I_am_alive.txt","wt");
		if (fo) {
			if (!ferror(fo)) {
				fprintf(fo,"I am active");
				fclose(fo);	fo = 0;
			}
		}
	}

/*
	std::string tags_ascii[24] = {"X","Y","TOF","PX","PY","PZ","P","PHI_DEG","PHIPOS","THETA_DEG","ENERGY",
								  "PSUMX","PSUMY","PSUMZ","PSUM","PSUM_PHI_DEG","PSUM_THETA_DEG",
								  "PRELX","PRELY","PRELZ","PREL","PREL_PHI_DEG","PREL_THETA_DEG","KER"};
*/
	Ueberstruct * Ueber = new Ueberstruct();
	Ueber->DEBUG_counter_xxx = 0;
	Ueber->config_file_name[0] = 0;
	Ueber->logfilepathname[0] = 0;
	Ueber->proj = Ueber->rec = Ueber->elec = 0;
	Ueber->rt = 0;
	Ueber->input_is_HDF5 = false;
	Ueber->config_version = 0;
	Ueber->Hist = 0;
	Ueber->current_proj_Sorter = Ueber->current_elec_Sorter = Ueber->current_rec_Sorter = 0;
	Ueber->proj_Sorter = Ueber->rec_Sorter = Ueber->elec_Sorter = 0;
	Ueber->parser = 0;
	Ueber->processed_data_file_names = 0;
	Ueber->input_data_file_name[0] = 0;
	Ueber->bot_mode = 0;
	Ueber->scan_channel = -1;
	Ueber->scan_in_data = false;
	Ueber->scan_num_vals = 0;
	Ueber->scan_factor = 1;
	Ueber->ext_logging = false;

	Ueber->adc_meta = 0;
	Ueber->CH = 0;
	Ueber->CH_status = -100;

	//split_pathname(argv[0],Ueber->applicationfilepath,0);
	get_current_directory(Ueber->applicationfilepath, 300);
	
	#ifndef dont_use_MFC
		printf("Running in folder:\n\"%s\"\n",Ueber->applicationfilepath);
	#endif

L100:

	if (Ueber->processed_data_file_names) {
		for (__int32 i=0;i<5000;i++) {
			if (Ueber->processed_data_file_names[i]) {delete[] Ueber->processed_data_file_names[i]; Ueber->processed_data_file_names[i]=0; continue;}
			break;
		}
	} else {
		Ueber->processed_data_file_names = new char * [5000];
		memset(Ueber->processed_data_file_names, 0, sizeof(char*)*5000);
	}

	if (!Ueber->bot_mode) {
		if (argc == 1) {
			printf("\nEnter \"lmf2root -help\" to display help screen.\n");
			printf("I will try to read \"config.txt\" now...\n");
			char filename[400];
			sprintf(filename,"%sconfig.txt",Ueber->applicationfilepath);
			sprintf(Ueber->config_file_name,filename);
			ascii_file = fopen(filename,"rt");
			if (!ascii_file || !File_is_readable(filename)) {Ueber->stop_reading_files = true; Ueber->bot_mode = false; goto L666;}
			fclose(ascii_file); ascii_file = 0;
			Ueber->parser = new config_file_parser(filename);
			if (!Ueber->parser) {Ueber->stop_reading_files = true; goto L666;}
			sprintf(Ueber->logfilepathname,"%s.log",filename);
		} else {
			if (argc>1) {
				if (!strcmp(argv[1],"help") || !strcmp(argv[1],"-help")) {
					print_help();
					goto L666;
				}
				if (!strcmp(argv[1],"bot") || !strcmp(argv[1],"-bot")) {
					printf("running in bot mode.\n");
					sprintf(Ueber->config_file_name,"bot_config.txt");
					Ueber->bot_mode = Ueber->bot_mode = true;
				}
				if (!Ueber->bot_mode) {
					ascii_file = fopen(argv[1],"rt");
					sprintf(Ueber->config_file_name,argv[1]);
					if (!ascii_file || !File_is_readable(argv[1])) {Ueber->stop_reading_files = true; Ueber->bot_mode = false; goto L666;}
					fclose(ascii_file); ascii_file = 0;
					Ueber->parser = new config_file_parser(argv[1]);
					if (!Ueber->parser) {Ueber->stop_reading_files = true; Ueber->bot_mode = false; goto L666;}
					sprintf(Ueber->logfilepathname,"%s.log",argv[1]);
				}
			} else {
				printf("wrong number of arguments.\n");
				gSystem->Sleep(2000);
				Ueber->stop_reading_files = true;
				Ueber->bot_mode = false;
				goto L666;
			}
		}
	}

	while (Ueber->bot_mode) {
		for (__int32 i=0;i<5*5;i++) gSystem->Sleep(200); // wait for 5 seconds
		delete_file("are_you_alive.txt");
		if (_kbhit()) {
			char c= _getch();
			while(_kbhit()) _getch();
			if (c == 'q' || c == 'Q') {
				Ueber->stop_reading_files = true;
				Ueber->bot_mode = false;
				goto L666;
			}
		}
		if (!FileExists("I_am_alive.txt")) {
			FILE * fo = fopen("I_am_alive.txt","wt");
			if (fo) {
				if (!ferror(fo)) {
					fprintf(fo,"I am active");
					fclose(fo);	fo = 0;
				}
			}
		}
		if (FileExists("bot_config.txt")) {
			Sleep(200);
			ascii_file = fopen("bot_config.txt","rt");
			if (!ascii_file) continue;
			if (ferror(ascii_file)) {
				fclose(ascii_file); ascii_file = 0;
				continue;
			}
			fclose(ascii_file); ascii_file = 0;
			if (Ueber->parser) {delete Ueber->parser; Ueber->parser = 0;}
			Ueber->parser = new config_file_parser("bot_config.txt");
			if (!Ueber->parser) {Ueber->stop_reading_files = true; Ueber->bot_mode = false; goto L666;}
			sprintf(Ueber->logfilepathname,"%s.log","bot_config.txt");
			break;
		}
	}

	
		Ueber->parallel_sorter = 0;

	Ueber->proj  =	new Det_struct();
	Ueber->rec  =	new Det_struct();
	Ueber->elec  =	new Det_struct();

	Ueber->peak_tracker1 = new peak_tracker_class(10000,-10.,10.,200); // number of events, left margin, right margin, number of bins
	Ueber->peak_tracker2 = new peak_tracker_class(10000,-10.,10.,200); // number of events, left margin, right margin, number of bins

	Ueber->old_time = get_time_s();

	Ueber->number_of_channels = 0;

	if (Ueber->rt) {delete Ueber->rt; Ueber->rt = 0;}
	Ueber->rt = new rootstuff();

	sprintf(Ueber->root_output_file_name_without_extension,"output");
	Ueber->LMF_output_file_name_without_extension[0] = 0;
	Ueber->stop_reading_files = false;
	Ueber->eventswritten = 0;

	Ueber->write_this_event_into_sorted_LMF = false;

	Ueber->applicationfilepath[0] = 0;

	Ueber->config_info = 0;
	Ueber->output_LMF_file_index = 0;
	Ueber->LMF_output = 0;
	Ueber->LMF_input = 0;
	Ueber->outputRootFile_handle = 0;
	Ueber->start_new_root_file = true;
	Ueber->error = 0;

	for (__int32 det = 0; det < 3;++det) {
		Det_struct * detector = 0;
		if (det == 0) detector = Ueber->proj;
		if (det == 1) detector = Ueber->rec;
		if (det == 2) detector = Ueber->elec;
		double sw = 10.;
		detector->sumu_watchdog = new peak_tracker_class(20000,-sw,sw,__int32(2.*sw/0.025+1.e-6)+1);
		detector->sumv_watchdog = new peak_tracker_class(20000,-sw,sw,__int32(2.*sw/0.025+1.e-6)+1);
		detector->sumw_watchdog = new peak_tracker_class(20000,-sw,sw,__int32(2.*sw/0.025+1.e-6)+1);
		detector->auto_calibration	= 0;
		detector->use_reconstruction	= false;
		detector->sorter = 0;
		detector->use_this_detector = false;
		detector->use_sum_tracker_layer_u = false;
		detector->use_sum_tracker_layer_v = false;
		detector->use_sum_tracker_layer_w = false;
		detector->Ueberstruct_pointer = (void*)Ueber;
	}


	double start_time = get_time_s();

	Ueber->fast_mode = 0;
	Ueber->Hist = new Histo(500000);

	parameter = new double[10000];
	Ueber->parameter = parameter;

	memset(parameter,0,10000*8);
	memset(master_cnt,0,NUM_CHANNELS*sizeof(int));
	memset(master_tdc_ns,0,NUM_CHANNELS*NUM_IONS*sizeof(double));

	Ueber->cnt = master_cnt;
	Ueber->tdc_ns = &master_tdc_ns[0][0];
	Ueber->eventcounter = 0;
	Ueber->filenumber = 0;

	Ueber->ipa_enable = false;
	Ueber->ipa_channel = 0;
	Ueber->ipa_num_events = 0;


//////////////////////////////////////////////////////////////////////////////////////


	Ueber->proj_Sorter = 0;
	Ueber->rec_Sorter  = 0;
	Ueber->elec_Sorter = 0;

	if (Ueber->logfilepathname[0]) {
		Ueber->logfile_handle = fopen(Ueber->logfilepathname,"wt");
		if (Ueber->logfile_handle) {
			if (ferror(Ueber->logfile_handle)) {
				fclose(Ueber->logfile_handle);
				Ueber->logfile_handle = 0;
			}
		}
	}

	while(!Ueber->stop_reading_files && !Ueber->parser->error_flag) {
		Ueber->parser->run();
		if (Ueber->parser->error_flag) {
			printf("\n%s\n",Ueber->parser->error_text.c_str());
			if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n%s\n",Ueber->parser->error_text.c_str()); fflush(Ueber->logfile_handle);}
			break;
		}

		if (Ueber->parser->external_command == "config_version") {
			Ueber->config_version = __int64(Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[0])+0.001);
			continue;
		}


		if (Ueber->parser->external_command == "readHDF5") {
			#ifdef _MSC_VER
				#if _MSC_VER < 1700
					Red();
					printf("This program was compiled using VS2010. So HDF5 is not supported.\nUse VS2015 or higher.\n");
					White();
				#endif
			#endif	
			if (Ueber->LMF_input) {
				Ueber->LMF_input->CloseInputLMF();
				delete Ueber->LMF_input;	Ueber->LMF_input = 0;
			}
			if (Ueber->HDF5_input) {
				Ueber->HDF5_input->CloseInputHDF5();
				delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
			}
			#ifdef LINUX
				if (Ueber->parser->number_of_external_arguments < 1) {
					printf("the multi readHDF5 command is not yet supported in Linux.\n");
					continue;
				}
			#endif

			if (fabs(parameter[819]-1234567) > 0.1) {
				Red();
				printf("error reading ADC-data: You must have the line 'execute \"ADC\\ADC_parameters.txt\"' in the config-file.\n");
				White();
			}

			if (Ueber->parser->number_of_external_arguments == 1) {
				char main_name[300]; main_name[0] = 0;
				if (Ueber->parser->number_of_external_arguments == 1) {
					Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
					sprintf(main_name, Ueber->parser->external_arguments[0].c_str());
				} else {
					sprintf(main_name, "*.h5*");
				}

				bool single_file_mode = true;
				for (unsigned __int32 i=0;i<strlen(main_name);i++) {
					if (main_name[i] == '*' || main_name[i] == '?') {single_file_mode = false; break;}
				}

				while (true) {
					char new_name[300]; new_name[0] = 0;
					bool new_file_found = false;
					if (!single_file_mode) {
						if (get_next_new_data_file_name(new_name, main_name, Ueber))
							new_file_found = true;
						else {
							White(true);
							printf("\n \n No further files of name '%s' found!\n\n",main_name);
							White(false);
						}
					} else {
						new_file_found = true;
						sprintf(new_name,main_name);
					}

					char c = my_kbhit();
					if (c) {
						if (c == 'q') break;
						if (c == 'Q') {Ueber->stop_reading_files = true; break;}
					}

					if (!new_file_found) {
						if (Ueber->parser->number_of_external_arguments == 1) break;
						for (__int32 i=0;i<40 && (!Ueber->stop_reading_files);i++) {gSystem->Sleep(100);} // wait for 4 seconds
						if (!Ueber->stop_reading_files) continue;
					}

					delete_file("I_am_busy.txt");
					FILE * fo = fopen("I_am_busy.txt","wt");
					if (fo) {
						if (!ferror(fo)) {
							fprintf(fo,"I am processing data");
							fclose(fo); fo = 0;
						}
					}
					Ueber->input_is_HDF5 = true;
					#ifdef LINUX
						ProcessInputDataFile((char *)(Ueber->parser->convert_backslash_to_slash(new_name), Ueber);
					#else
						ProcessInputDataFile((char *)(new_name), Ueber);
					#endif
					delete_file("I_am_busy.txt");

					if (single_file_mode) break;
					if(Ueber->stop_reading_files) break;
				}
				continue;
			}
			if (Ueber->parser->number_of_external_arguments > 1) {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
		}


		if (Ueber->parser->external_command == "readLMF") {
			if (Ueber->LMF_input) {
				Ueber->LMF_input->CloseInputLMF();
				delete Ueber->LMF_input;	Ueber->LMF_input = 0;
			}
			if (Ueber->HDF5_input) {
				Ueber->HDF5_input->CloseInputHDF5();
				delete Ueber->HDF5_input;	Ueber->HDF5_input = 0;
			}
			#ifdef LINUX
				if (Ueber->parser->number_of_external_arguments < 1) {
					printf("the multi readLMF command is not yet supported in Linux.\n");
					continue;
				}
			#endif

			if (Ueber->parser->number_of_external_arguments == 1) {
				char main_name[300]; main_name[0] = 0;
				if (Ueber->parser->number_of_external_arguments == 1) {
					Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
					sprintf(main_name, Ueber->parser->external_arguments[0].c_str());
				} else {
					sprintf(main_name, "*.lmf*");
				}

				bool single_file_mode = true;
				for (unsigned __int32 i=0;i<strlen(main_name);i++) {
					if (main_name[i] == '*' || main_name[i] == '?') {single_file_mode = false; break;}
				}

				while (true) {
					char new_name[300]; new_name[0] = 0;
					bool new_file_found = false;
					if (!single_file_mode) {
						if (get_next_new_data_file_name(new_name, main_name, Ueber))
							new_file_found = true;
						else {
							White(true);
							printf("\n \n No further files of name '%s' found!\n\n",main_name);
							White(false);
						}
					} else {
						new_file_found = true;
						sprintf(new_name,main_name);
					}

					char c = my_kbhit();
					if (c) {
						if (c == 'q') break;
						if (c == 'Q') {Ueber->stop_reading_files = true; break;}
					}

					if (!new_file_found) {
						if (Ueber->parser->number_of_external_arguments == 1) break;
						for (__int32 i=0;i<40 && (!Ueber->stop_reading_files);i++) {gSystem->Sleep(100);} // wait for 4 seconds
						if (!Ueber->stop_reading_files) continue;
					}

					delete_file("I_am_busy.txt");
					FILE * fo = fopen("I_am_busy.txt","wt");
					if (fo) {
						if (!ferror(fo)) {
							fprintf(fo,"I am processing data");
							fclose(fo); fo = 0;
						}
					}

					Ueber->input_is_HDF5 = false;
					#ifdef LINUX
						ProcessInputDataFile((char *)(Ueber->parser->convert_backslash_to_slash(new_name), Ueber);
					#else
						ProcessInputDataFile((char *)(new_name), Ueber);
					#endif
					delete_file("I_am_busy.txt");

					if (single_file_mode) break;
					if(Ueber->stop_reading_files) break;
				}
				continue;
			}
			if (Ueber->parser->number_of_external_arguments > 1) {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
		}


		if (Ueber->parser->external_command == "botreadlmf") continue;


		if (Ueber->parser->external_command == "readROOTfile") {
			if (Ueber->parser->number_of_external_arguments == 1) {
				char main_name[300]; main_name[0] = 0;
				if (Ueber->parser->number_of_external_arguments == 1) {
					Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
					sprintf(main_name, Ueber->parser->external_arguments[0].c_str());
				} else {
					sprintf(main_name, "*.root*");
				}
				
				bool single_file_mode = true;
				for (unsigned __int32 i=0;i<strlen(main_name);i++) {
					if (main_name[i] == '*' || main_name[i] == '?') {single_file_mode = false; break;}
				}

				while (true) {
					char new_name[300]; new_name[0] = 0;
					bool new_file_found = false;
					if (!single_file_mode) {
						if (get_next_new_data_file_name(new_name, main_name, Ueber)) {
							new_file_found = true;
						} else {
							White(true);
							printf("\n \n No further files of name '%s' found!\n\n",main_name);
							White(false);
						}
					} else {
						new_file_found = true;
						sprintf(new_name,main_name);
					}

					char c = my_kbhit();
					if (c) {
						if (c == 'q') break;
						if (c == 'Q') {Ueber->stop_reading_files = true; break;}
					}

					if (!new_file_found) {
						if (Ueber->parser->number_of_external_arguments == 1) break;
						for (__int32 i=0;i<40 && (!Ueber->stop_reading_files);i++) {gSystem->Sleep(100);} // wait for 4 seconds
						if (!Ueber->stop_reading_files) continue;
					}

					delete_file("I_am_busy.txt");
					FILE * fo = fopen("I_am_busy.txt","wt");
					if (fo) {
						if (!ferror(fo)) {
							fprintf(fo,"I am processing data");
							fclose(fo); fo = 0;
						}
					}
					#ifdef LINUX
						ProcessRootFile((char *)(Ueber->parser->convert_backslash_to_slash(new_name), Ueber);
					#else
						ProcessRootFile((char *)(new_name), Ueber);
					#endif
					delete_file("I_am_busy.txt");

					if (single_file_mode) break;
				}
				continue;
			}
			if (Ueber->parser->number_of_external_arguments != 1) {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
		}

		if (Ueber->parser->external_command == "set_root_output_file_name") {
			if (Ueber->parser->number_of_external_arguments == 1) {
				if (Ueber->outputRootFile_handle) {
					Ueber->outputRootFile_handle->Write();
					Ueber->Hist->Reset();
					Ueber->outputRootFile_handle->Close();
					Ueber->outputRootFile_handle = 0;
				}
				Ueber->start_new_root_file = true;
				Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
				#ifdef LINUX
					sprintf(Ueber->root_output_file_name_without_extension,"%s",Ueber->parser->convert_backslash_to_slash(Ueber->parser->external_arguments[0]).c_str());
				#else
					sprintf(Ueber->root_output_file_name_without_extension,"%s",Ueber->parser->external_arguments[0].c_str());

					if(Ueber->ext_logging) {

						char buffer[MAX_PATH];
						GetCurrentDirectory(MAX_PATH,buffer);
						char logdir[255];
						sprintf(logdir,"%s\\logs",buffer);
						CreateDirectory (logdir, NULL);
						sprintf(logdir,"%s\\logs\\%s.root",buffer,Ueber->root_output_file_name_without_extension);
						CreateDirectory (logdir, NULL);

						char tmpfilename[400];
                        sprintf(tmpfilename,"%s//config.txt",logdir);
                        CopyFile(Ueber->parser->file_name, tmpfilename, FALSE);
                        sprintf(tmpfilename,"%s//%s",logdir,"ColAHelL.cfg");
                        CopyFile("ColAHelL.cfg", tmpfilename, FALSE);
                        sprintf(tmpfilename,"%s//%s",logdir,"USER_Analysis.cpp");
                        CopyFile("ColAHelL\\UserSrc\\USER_Analysis.cpp", tmpfilename, FALSE);
					}
				#endif
				for (__int32 i=0;i<5000;i++) {
					if (!Ueber->processed_data_file_names[i]) {
						Ueber->processed_data_file_names[i] = new char[strlen(Ueber->root_output_file_name_without_extension)+1+5];
						sprintf(Ueber->processed_data_file_names[i],"%s.root",Ueber->root_output_file_name_without_extension);
						break;
					}
				}
				continue;
			} else {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
		}

		if (Ueber->parser->external_command == "set_LMF_output_file_name") {
			if (Ueber->parser->number_of_external_arguments == 1) {
				Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
				#ifdef LINUX
					sprintf(Ueber->LMF_output_file_name_without_extension,"%s",Ueber->parser->convert_backslash_to_slash(Ueber->parser->external_arguments[0]).c_str());
				#else
					sprintf(Ueber->LMF_output_file_name_without_extension,"%s",Ueber->parser->external_arguments[0].c_str());
				#endif
				for (__int32 i=0;i<5000;i++) {
					if (!Ueber->processed_data_file_names[i]) {
						Ueber->processed_data_file_names[i] = new char[strlen(Ueber->LMF_output_file_name_without_extension)+1+4];
						sprintf(Ueber->processed_data_file_names[i],"%s.lmf",Ueber->LMF_output_file_name_without_extension);
						break;
					}
				}
				continue;
			} else {
				if (Ueber->parser->number_of_external_arguments == 0) {
					Ueber->LMF_output_file_name_without_extension[0] = 0;
					continue;
				} else {
					printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
					if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
					Ueber->stop_reading_files = true;
					break;
				}
			}
		}

		if (Ueber->parser->external_command == "cd") {
			if (Ueber->parser->number_of_external_arguments == 0) {
				printf("\ncurrent working path is:\n%s\n",Ueber->applicationfilepath);
			}
			if (Ueber->parser->number_of_external_arguments == 1) {
				Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
				if (Ueber->parser->external_arguments[0].find(' ', 0)>=0 || Ueber->parser->external_arguments[0].find('\t', 0)>=0) {
					Ueber->parser->external_arguments[0] = "\"" + Ueber->parser->external_arguments[0] + "\"";
				}
				if (!_chdir(Ueber->parser->external_arguments[0].c_str())) {
					get_current_directory(Ueber->applicationfilepath, 300);
					continue;
				}
				continue;
			}
			printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
			if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
			Ueber->stop_reading_files = true;
			break;
		}

		if (Ueber->parser->external_command == "system") {
			if (Ueber->parser->number_of_external_arguments == 1) {
				Ueber->parser->strip_quotation_marks(Ueber->parser->external_arguments[0]);
				if (Ueber->logfile_handle) fflush(Ueber->logfile_handle);
				if (Ueber->outputRootFile_handle) Ueber->outputRootFile_handle->Flush();
				if (Ueber->LMF_output) {
					if (Ueber->LMF_output->output_lmf) Ueber->LMF_output->output_lmf->flush();
				}
				__int32 result;
				result = system(Ueber->parser->external_arguments[0].c_str());
				if (result != -1) {
					get_current_directory(Ueber->applicationfilepath, 300);
					continue;
				}
			}
			printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
			if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
			Ueber->stop_reading_files = true;
			break;
		}


		if (Ueber->parser->external_command == "beep") {
			if (Ueber->parser->number_of_external_arguments == 0) {
				printf("%c",7);
				continue;
			} else {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
		}
// Interactive parameter adjustment: ask T.Jahnke for support
//--------------------------------------------------------------------------------------


		// we want to use the interactive parameter adjustment feature
		if (Ueber->parser->external_command == "invoke_ipa") {
			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			
			if(Ueber->ipa_enable) { // ipa already invoked before!
				Red(true);
				printf("Only one IPA-instance can be used!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			else {
				Green(true);
				printf("\nInteractive parameter adjustment (IPA) requested:");
				White(false);

				if (Ueber->parser->number_of_external_arguments > 0) {

					
					Ueber->ipa_enable = true;
					Ueber->ipa_num_events = __int32(pargs[0]+0.001);
					if(Ueber->parser->number_of_external_arguments > 1) {
						Ueber->ipa_channel = __int32(pargs[1]+0.001);
						printf("\nWill collect %d events belonging to channel %d...\n",(int)Ueber->ipa_num_events, (int)Ueber->ipa_channel);
					} else {
						Ueber->ipa_channel = -1;
						printf("\nWill collect %d events belonging to any channel...\n",(int)Ueber->ipa_num_events);
					}

					printf("(Or press 'i' to start IPA at any time with number of events read so far..)\n");
					White(true);
					printf("Hint: the sequence ColAHelL uses to calculate the momenta is...\n");
					printf("  1. Shift, then stretch the detector image.\n");
					printf("  2. Rotate the detector image.\n");
					printf("  3. Apply jet offset (ions only).\n");
					printf("  4. Calculate the momenta.\n");
					printf("  5. Shift, then stretch the momenta.\n");
					printf("  6. Mirror the momenta.\n");
					Red(true);
					printf("\n Order of 1. + 2. might be reversed depending on user settings.\n");
					White(false);

					Ueber->IPA = new ipa_class();
					Ueber->IPA->init(Ueber->ipa_channel,Ueber->ipa_num_events, Ueber->rt);
					Ueber->IPA->Setup(Ueber->CH->GetColAHelL());

					continue;

				} else {
					Red(true);
					printf(" ..failed!\n");
					printf("Wrong number of arguments!\n");
					White(false);
					printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
					Ueber->stop_reading_files = true;
					break;
				}
			}
		} // end of  if (Ueber->parser->external_command == "invoke_interactive_parameter")

		// Add parameters to IPA
		if (Ueber->parser->external_command == "add_ipa_parameter") {
			double pargs[25];
			for(int i=1;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			

			if (Ueber->parser->number_of_external_arguments == 5) {

				char * name = new char[Ueber->parser->external_arguments[0].length()+1];
				strcpy(name,Ueber->parser->external_arguments[0].c_str());

				Green(true);
				printf("\nAdding IPA parameter: ");
				White(true);					
				printf("%s",name);
				White(false);
					

				Ueber->IPA->add_par(pargs[1],pargs[2],pargs[3],pargs[4],name);

				continue;

			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "add_ipa_parameter")

		// request IPA PIPICO tool
		if (Ueber->parser->external_command == "ipa_pipico_tool") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA PIPICO tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}

			double pargs[25];
			for(int i=1;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			

			if (Ueber->parser->number_of_external_arguments == 0) {

				Green(true);
				printf("\nIPA PIPICO tool requested using standard settings...");
				White(true);					

				Ueber->IPA->enable_PPCO(1000,500.0,10000.0,true);

				continue;
			}
			if (Ueber->parser->number_of_external_arguments == 1) {

				Green(true);
				printf("\nIPA PIPICO tool requested");
				White(true);					

				if(Ueber->parser->external_arguments[0]=="true") 
					Ueber->IPA->enable_PPCO(500,500.0,5000.0,true);
				else
					Ueber->IPA->enable_PPCO(500,500.0,5000.0,false);

				continue;
			}
			if (Ueber->parser->number_of_external_arguments == 4) {

				Green(true);
				printf("\nIPA PIPICO tool requested");
				White(true);					

				if(Ueber->parser->external_arguments[0]=="true") 
					Ueber->IPA->enable_PPCO((int)pargs[1],pargs[2],pargs[3],true);
				else
					Ueber->IPA->enable_PPCO((int)pargs[1],pargs[2],pargs[3],false);

				continue;
			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_pipico_tool")

		// request IPA LUT tool
		if (Ueber->parser->external_command == "ipa_LUT_tool_eV") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA LookUpTable tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}

			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			

			if (Ueber->parser->number_of_external_arguments > 0) {

				double energy = pargs[0];
				double mom=sqrt(2*energy/27.2);

				int electrons = 1;
				int bins_mom = 20;
				int bins_phi = 36;
				int bins_ctheta = 18;

				Green(true);
				printf("\nIPA LUT tool requested for electrons of energy = %.1f eV (p=%.1f au)...",energy, mom);
				White(true);					

				if (Ueber->parser->number_of_external_arguments > 3) {
					bins_mom = (int)pargs[1];
					bins_phi = (int)pargs[2];
					bins_ctheta = (int)pargs[3];					
				}

				if (Ueber->parser->number_of_external_arguments == 5) {
					electrons = (int)pargs[4];
				}

				Ueber->IPA->enable_LUT(mom, bins_mom, bins_phi, bins_ctheta, electrons);
				continue;			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_LUT_tool")
		// request IPA LUT tool
		if (Ueber->parser->external_command == "ipa_LUT_tool_au") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA LookUpTable tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}

			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			

			if (Ueber->parser->number_of_external_arguments > 0) {

				double mom=pargs[0];
				double energy = mom*mom/2.0*27.2;

				int electrons = 1;
				int bins_mom = 20;
				int bins_phi = 36;
				int bins_ctheta = 18;

				Green(true);
				printf("\nIPA LUT tool requested for electrons of momentum = %.1f au (E=%.1f eV)...",mom,energy);
				White(true);					

				if (Ueber->parser->number_of_external_arguments > 3) {
					bins_mom = (int)pargs[1];
					bins_phi = (int)pargs[2];
					bins_ctheta = (int)pargs[3];					
				}

				if (Ueber->parser->number_of_external_arguments == 5) {
					electrons = (int)pargs[4];
				}

				Ueber->IPA->enable_LUT(mom, bins_mom, bins_phi, bins_ctheta, electrons);
				continue;			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_LUT_tool_au")

		// request IPA wiggle tool
		if (Ueber->parser->external_command == "ipa_wiggle_tool") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA wiggle tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
/*			if(!Ueber->CTof) {
				Red(true);
				printf("IPA wiggle tool needs proper time-of-flight calculation parameters !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
*/
			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			
			// use standard parameters
			if (Ueber->parser->number_of_external_arguments == 0) {

				Green(true);
				printf("\nIPA wiggle tool requested ...");
				White(true);
				Ueber->IPA->enable_wiggletool(400, -10.0, 200.0);
				continue;
			}

			if (Ueber->parser->number_of_external_arguments == 3) {

				Green(true);
				printf("\nIPA wiggle tool requested ...");
				White(true);
				Ueber->IPA->enable_wiggletool(int(pargs[0]), pargs[1], pargs[2]);
				continue;
			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_wiggle_tool")

		// request IPA fish tool
		if (Ueber->parser->external_command == "ipa_fish_tool") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA fish tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}

			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			
			// use standard parameters
			if (Ueber->parser->number_of_external_arguments == 0) {

				Green(true);
				printf("\nIPA fish tool requested ...");
				White(true);
				Ueber->IPA->enable_fishtool(400, -10.0, 200.0);
				continue;
			}

			if (Ueber->parser->number_of_external_arguments == 3) {

				Green(true);
				printf("\nIPA fish tool requested ...");
				White(true);
				Ueber->IPA->enable_fishtool(int(pargs[0]), pargs[1], pargs[2]);
				continue;
			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_fish_tool")

		// request IPA ion tof tool
		if (Ueber->parser->external_command == "ipa_ion_tof_tool") {

			if(!Ueber->ipa_enable) {
				Red(true);
				printf("IPA ion tof tool request needs to occur after 'invoke_ipa' !\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}
			
			// use standard parameters
			if (Ueber->parser->number_of_external_arguments == 0) {

				Green(true);
				printf("\nIPA ion tof tool requested ...");
				White(true);
				Ueber->IPA->enable_iontoftool(1000, 10.0, 10000.0);
				continue;
			}

			if (Ueber->parser->number_of_external_arguments == 3) {

				Green(true);
				printf("\nIPA ion tof tool requested ...");
				White(true);
				Ueber->IPA->enable_iontoftool(int(pargs[0]), pargs[1], pargs[2]);
				continue;
			
			} else {
				Red(true);
				printf(" ..failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			
		} // end of  if (Ueber->parser->external_command == "ipa_ion_tof_tool")
// end of interactive paramter adjustment stuff

//--------------------------------------------------------------------------------------
// ColAHelL definitions: ask T.Jahnke for support
//--------------------------------------------------------------------------------------
		if (Ueber->parser->external_command == "use_ColAHelL" || Ueber->parser->external_command == "use_ColaHell" || Ueber->parser->external_command == "use_colahell" ) {
			
			White(true);
			printf("\nYou're using ColAHelL! Good job! ColAHelL is your friend!\n");
			White(false);

			if (Ueber->parser->number_of_external_arguments == 0) {

				Ueber->CH = new CH_Int();
				Ueber->use_CH = true;
				Ueber->CH_status = Ueber->CH->ReadConfig();
				if(Ueber->CH_status < -89) {
					Red(true);
					printf(" Shutting down LMF2Root!\n");
					White(false);
					Ueber->stop_reading_files = true;
					break;
				}
				continue;
			}
			if (Ueber->parser->number_of_external_arguments == 1) {

				Ueber->CH = new CH_Int();
				Ueber->use_CH = true;
				Green(true);
				printf("Reading '%s'...\n\n",Ueber->parser->external_arguments[0].c_str());
				White(true);
				Ueber->CH_status = Ueber->CH->ReadConfig(Ueber->parser->external_arguments[0].c_str());
				if(Ueber->CH_status<0) {
					Red(true);
					printf(" Shutting down LMF2Root!\n");
					White(false);
					Ueber->stop_reading_files = true;
					break;
				}
				continue;
			
			} else {
				Red(true);
				printf(" ..Oh no!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("Error in line %i in file %s.\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
		} // end if (Ueber->parser->external_command == "use_ColAHelL")

		if (Ueber->parser->external_command == "replace_ColAHelL_cfg" || Ueber->parser->external_command == "replace_ColaHell_cfg" || Ueber->parser->external_command == "replace_colahell_cfg" ) {

			if(!Ueber->use_CH) {
				Red(true);
				printf(" ..Oh no!\n");
				printf("You need to init ColAHelL before replacing its configuration!\n");
				White(false);
				printf("Error in line %i in file %s.\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			White(true);
			printf("\nReplacing ColAHelL configuration!\n");
			White(false);

			if (Ueber->parser->number_of_external_arguments == 1) {

				Green(true);
				printf("Reading '%s'...\n\n",Ueber->parser->external_arguments[0].c_str());
				White(true);
				Ueber->CH_status = Ueber->CH->ReplaceConfig(Ueber->parser->external_arguments[0].c_str());
				if(Ueber->CH_status<0) {
					Red(true);
					printf(" Shutting down LMF2Root!\n");
					White(false);
					Ueber->stop_reading_files = true;
					break;
				}
				continue;
			
			} else {
				Red(true);
				printf(" ..Oh no!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("Error in line %i in file %s.\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
		} // end if (Ueber->parser->external_command == "replace_ColAHelL_cfg")

		if (Ueber->parser->external_command == "scan_channel") {
			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}

			White(true);
			printf("TDC channel used for scanning values is: ");
			White(false);

			if (Ueber->parser->number_of_external_arguments > 0) {
				Ueber->scan_channel = (int)(pargs[0]+0.001);
				Green(true);
				printf("%i.\n",Ueber->scan_channel); 
				White(false);
			} else {
				Red(true);
				printf(" .. this failed!\n");
				printf("Wrong number of arguments!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			if (Ueber->parser->number_of_external_arguments > 1) {
				Ueber->scan_num_vals = (int)(pargs[1]+0.001);
				if(Ueber->scan_num_vals>16) {
					Red(true);
					printf("Only up to 16 scan parameters are supported!");
					White(false);
					Ueber->scan_num_vals=16;
				}
			}
			if (Ueber->parser->number_of_external_arguments > 2) {
				Ueber->scan_factor = (int)(pargs[2]+0.001);
			}

			printf("Number of scanned parameters: ");
			Green(true);
			printf("%i.\n",Ueber->scan_num_vals); 
			White(false);
			
			printf("The scanning conversion factor is: ");
			Green(true);
			printf("%i.\n",Ueber->scan_factor); 
			White(false);

			continue;
		} // end of  if (Ueber->parser->external_command == "scan_channel")

		if (Ueber->parser->external_command == "scan_parameter") {
			double pargs[25];
			for(int i=0;i<Ueber->parser->number_of_external_arguments;i++) {
				pargs[i] = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[i]);
			}

			White(true);
			printf("LMF parameters were used to store scanning values.\n");
			White(false);

			if (Ueber->parser->number_of_external_arguments > 0) {
				Ueber->scan_channel = 999;
			} else {
				Red(true);
				printf(" .. this failed!\n");
				printf("No parameter numbers provided!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			if (Ueber->parser->number_of_external_arguments > 1) {
				Ueber->scan_num_vals = (int)(pargs[0]+0.001);
				if(Ueber->scan_num_vals>16) {
					Red(true);
					printf("Only up to 16 scan parameters are supported!");
					White(false);
					Ueber->scan_num_vals=16;
				}
			}

			printf("Number of stored parameters: ");
			Green(true);
			printf("%i\n",Ueber->scan_num_vals); 
			White(false);
			
			if (Ueber->parser->number_of_external_arguments == Ueber->scan_num_vals + 1) {
				printf("Stored parameters are: ");
				Green(true);
				for(int i=0;i<Ueber->scan_num_vals;i++){
					Ueber->scan_par_nums[i] = (int)(pargs[1 + i]+0.001);
					printf("%i ",Ueber->scan_par_nums[i]);
				}
				printf("\n");
				White(false);
			} else {
				Red(true);
				printf(" .. this failed!\n");
				printf("Wrong numner of parameters provided!\n");
				White(false);
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				Ueber->stop_reading_files = true;
				break;
			}
			continue;
		} // end of  if (Ueber->parser->external_command == "scan_channel")

		if (Ueber->parser->external_command == "extensive_log") {
			
			White(true);
			printf("\nExtensive logging enabled!\n");
			White(false);
			Ueber->ext_logging = true;
			continue;
		} // end if (Ueber->parser->external_command == "extensive_log")

		if (Ueber->parser->external_command == "set_master_folder") {
			if(Ueber->use_CH) {
				if (Ueber->parser->number_of_external_arguments == 1) {

					Green(true);
					printf("\nSetting master folder to '%s'...\n",Ueber->parser->external_arguments[0].c_str());
					White(false);
					Ueber->CH->SetMasterFolder(Ueber->parser->external_arguments[0].c_str());
					continue;
			
				} else {
					Red(true);
					printf(" ..Oh no!\n");
					printf("Wrong number of arguments!\n");
					White(false);
					printf("Error in line %i in file %s.\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name);
					Ueber->stop_reading_files = true;
					break;
				}
			}
		} // end if (Ueber->parser->external_command == "extensive_log")
				////////////////////////////
		if (Ueber->parser->external_command == "parameter") {
			__int32 index = __int32(Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[0])+0.001);
			double value =  Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[1]);
			if (index < 10000 && Ueber->parser->number_of_external_arguments == 2) {
				parameter[index] = value;
			} else {
				printf("error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				if (index >= 10000) printf("(parameter index > 9999)\n");
				Ueber->stop_reading_files = true;
				break;
			}
			continue;
		}


		////////////////////////////
		if (Ueber->parser->external_command == "feed_parameters_into_sorters") {
			if (Ueber->parser->number_of_external_arguments != 0) {
				printf("wrong number of arguments in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
				Ueber->stop_reading_files = true;
				break;
			}
			if (parameter[100]>0.5) {
				#ifdef _DEBUG
					assert(Ueber->proj_Sorter == 0);
				#endif
				Ueber->proj->sorter = Ueber->proj_Sorter = new sort_class();
				if (!fill_parameters_into_sorters(Ueber->proj_Sorter,Ueber->proj,100)) {
					Ueber->stop_reading_files = true;
					break;
				}
			}
			if (parameter[200]>0.5) {
				#ifdef _DEBUG
					assert(Ueber->rec_Sorter == 0);
				#endif
				Ueber->rec->sorter = Ueber->rec_Sorter = new sort_class();
				if (!fill_parameters_into_sorters(Ueber->rec_Sorter ,Ueber->rec ,200)) {
					Ueber->stop_reading_files = true;
					break;
				}
			}
			if (parameter[300]>0.5) {
				#ifdef _DEBUG
					assert(Ueber->elec_Sorter == 0);
				#endif

				Ueber->elec->sorter = Ueber->elec_Sorter = new sort_class();
				if (!fill_parameters_into_sorters(Ueber->elec_Sorter,Ueber->elec,300)) {
					Ueber->stop_reading_files = true;
					break;
				}
			}

			__int32 i = (Ueber->elec->auto_calibration?1:0) +
					(Ueber->rec->auto_calibration ?1:0) +
					(Ueber->proj->auto_calibration?1:0);
			if (i>1) {
				printf("Don't calibrate more than one detector at the same time.\n");
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"Don't calibrate more than one detectors at the same time.\n"); fflush(Ueber->logfile_handle);}
				gSystem->Sleep(5000);
				Ueber->stop_reading_files = true;
				break;
			}

			continue;
		}



		
		////////////////////////////
		while (true) {
			if (Ueber->parser->external_command == "set_point") {
				if (Ueber->parser->number_of_external_arguments != 5) {
					printf("wrong number of arguments in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
					if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
					Ueber->stop_reading_files = true;
					break;
				}
				std::string sum_pos   = Ueber->parser->extract_first_string(Ueber->parser->external_arguments[0]);

				std::string sort_name = Ueber->parser->extract_first_string(Ueber->parser->external_arguments[1]);
				sort_class * sorter   = 0;
				if (sort_name == "proj" && Ueber->proj->use_this_detector) sorter = Ueber->proj_Sorter;
				if (sort_name == "rec" && Ueber->rec->use_this_detector)  sorter = Ueber->rec_Sorter;
				if (sort_name == "elec" && Ueber->elec->use_this_detector) sorter = Ueber->elec_Sorter;
				
				std::string layer = Ueber->parser->extract_first_string(Ueber->parser->external_arguments[2]);
				double pos_ns = Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[3]);
				double corr =   Ueber->parser->read_double_from_string(Ueber->parser->external_arguments[4]);

				if (fabs(corr) < 50.) {
					if (sum_pos == "sum") {
						if (sorter) {
							if (sorter->use_sum_correction) {
								if (layer == "u") sorter->signal_corrector->sum_corrector_U->set_point(pos_ns,corr);
								if (layer == "v") sorter->signal_corrector->sum_corrector_V->set_point(pos_ns,corr);
								if (sorter->use_HEX) {
									if (layer == "w") sorter->signal_corrector->sum_corrector_W->set_point(pos_ns,corr);
								}
							}
						}
					}
					if (sum_pos == "pos") {
						if (sorter) {
							if (sorter->use_pos_correction) {
								if (layer == "u") sorter->signal_corrector->pos_corrector_U->set_point(pos_ns,corr);
								if (layer == "v") sorter->signal_corrector->pos_corrector_V->set_point(pos_ns,corr);
								if (sorter->use_HEX) {
									if (layer == "w") sorter->signal_corrector->pos_corrector_W->set_point(pos_ns,corr);
								}
							}
						}
					}
				}
				{Ueber->parser->external_command = "known_command"; break;}
			}


			std::string &command = Ueber->parser->external_command;
			__int32 i1,i2,i3;
			i1 = command.find("reset_");
			i2 = command.find("_det_correction_table");
			if (i1 > -1 && i2 > -1 && i2 == i1 + 7) {Ueber->parser->external_command = "known_command"; break;}

			i1 = command.find("check_f");
			i2 = command.find("_on_");
			i3 = command.find("_det");
			if (i1 >-1 && i2 >-1 && i3>-1 && i2 == i1 + 8 && i3 == i2 + 5) {
				Ueber->parser->external_command = "known_command";
				break;
			}

			i1 = command.find("check_w_offset_on_");
			if (i1  > -1) {Ueber->parser->external_command = "known_command"; break;}			

			if (Ueber->parser->external_command != "known_command") break;

			std::string new_command = "";
			__int32 step_counter = 0;
			i1 = command.find(",set_correction_point,sum,");
			if (i1 > -1) {new_command = "set_point"; Ueber->parser->external_arguments[0] = "sum"; step_counter++;}
			i1 = command.find(",set_correction_point,pos,");
			if (i1 > -1) {new_command = "set_point"; Ueber->parser->external_arguments[0] = "pos"; step_counter++;}
			if (step_counter == 1) {
				i1 = command.find(",i_det,");
				if (i1  > -1) {Ueber->parser->external_arguments[1] = "rec ";   step_counter++;}
				i1 = command.find(",e_det,");
				if (i1 > -1) {Ueber->parser->external_arguments[1] = "elec ";  step_counter++;}
				i1 = command.find(",p_det,");
				if (i1 > -1) {Ueber->parser->external_arguments[1] = "proj ";  step_counter++;}
				i1 =command.find(",p2_det,");
				if (i1 >-1) {Ueber->parser->external_arguments[1] = "proj2 "; step_counter++;}
			}
			std::string values = "";
			if (step_counter == 2) {
				i1 = command.find(",u,");
				if (i1 > -1) {Ueber->parser->external_arguments[2] = "u";  values = command.substr(i1+3,999); step_counter++;}
				i1 = command.find(",v,");
				if (i1 > -1) {Ueber->parser->external_arguments[2] = "v";  values = command.substr(i1+3,999); step_counter++;}
				i1 = command.find(",w,");
				if (i1 > -1) {Ueber->parser->external_arguments[2] = "w";  values = command.substr(i1+3,999); step_counter++;}
			}

			if (step_counter == 3) {
				i1 = values.find(",");
				if (i1 > -1) {
					step_counter++;
					Ueber->parser->external_arguments[3] = values.substr(0,i1);
					Ueber->parser->external_arguments[4] = values.substr(i1+1,999);
				}
			}
			if (step_counter == 4) {
				Ueber->parser->external_command = new_command;
				Ueber->parser->number_of_external_arguments = 5;
				continue;
			}

		}

		if (Ueber->parser->external_command == "known_command") continue;

		if (Ueber->stop_reading_files) break;

		if (Ueber->parser->number_of_open_files == 0 || Ueber->parser->external_command == "end" || Ueber->parser->external_command == "exit" || Ueber->parser->external_command == "quit") {
			Ueber->stop_reading_files = true;
			break;
		}

		printf("unknown command:\n%s\nin line %i in file %s\n",Ueber->parser->external_command.c_str(), Ueber->parser->GetCurrentLineNumber()+1,Ueber->parser->file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"error in line %i in file %s\n",Ueber->parser->GetCurrentLineNumber(),Ueber->parser->file_name); fflush(Ueber->logfile_handle);}
		Ueber->stop_reading_files = true;
		break;
	}

	double stop_time = get_time_s();
	printf("\nLMF2Root was running for %lg seconds.\n",stop_time-start_time);
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\nLMF2Root was running for %lg seconds.\n",stop_time-start_time); fflush(Ueber->logfile_handle);}

	// do calibration:
	for (__int32 detector = 0;detector < 3;++detector) {
		Det_struct * detstruct;
		std::string det_name;
		if (detector == 0) {detstruct = Ueber->proj; det_name = "proj";}
		if (detector == 1) {detstruct = Ueber->rec;  det_name = "rec";}
		if (detector == 2) {detstruct = Ueber->elec; det_name = "elec";}
		if (detstruct->auto_calibration && (!detstruct->use_reconstruction) && detstruct->sorter) {
			detstruct->sorter->do_calibration();

	/*	if (detstruct->pos_walk_calibrator) { // XXX
				profile_class * prof = detstruct->pos_walk_calibrator->pos_profile[0][0];
				double f = Ueber->proj_Sorter->fu;
				double bx = f*0.5*(prof->center_of_right_bin - prof->center_of_left_bin)/(prof->number_of_columns-1);
				double by = f*0.5*(prof->center_of_upper_bin - prof->center_of_lower_bin)/(prof->number_of_rows-1);
				Ueber->Hist->fill(555,"Non-linearity_layer_u",-1000.,-1000.,"x [mm]","dx [mm]",prof->number_of_columns,prof->center_of_left_bin*f-bx,prof->center_of_right_bin*f+bx,prof->number_of_rows,prof->center_of_lower_bin*f-by,prof->center_of_upper_bin*f+by);
				for (__int32 ix = 0;ix < prof->number_of_columns;++ix) {
					for (__int32 iy = 0;iy < prof->number_of_rows;++iy) {
						double c = prof->double_source_matrix[ix][iy];
						if (c != 0.) {
							__int32 sadsad = 0;
						}
						Ueber->Hist->Fill(555,prof->get_bin_center_x(ix)*f,prof->get_bin_center_y(iy)*f,c*f);
						//Ueber->Hist->SetBinContentAt(555,prof->get_bin_center_x(ix),prof->get_bin_center_y(iy),c);
					}
				}
				f = Ueber->proj_Sorter->fv;
				prof = detstruct->pos_walk_calibrator->pos_profile[1][0];
				bx = f*0.5*(prof->center_of_right_bin - prof->center_of_left_bin)/(prof->number_of_columns-1);
				by = f*0.5*(prof->center_of_upper_bin - prof->center_of_lower_bin)/(prof->number_of_rows-1);
				Ueber->Hist->fill(556,"Non-linearity_layer_v",-1000.,-1000.,"x [mm]","dx [mm]",prof->number_of_columns,prof->center_of_left_bin*f-bx,prof->center_of_right_bin*f+bx,prof->number_of_rows,prof->center_of_lower_bin*f-by,prof->center_of_upper_bin*f+by);
				for (__int32 ix = 0;ix < prof->number_of_columns;++ix) {
					for (__int32 iy = 0;iy < prof->number_of_rows;++iy) {
						double c = prof->double_source_matrix[ix][iy];
						if (c != 0.) {
							__int32 sadsad = 0;
						}
						Ueber->Hist->Fill(556,prof->get_bin_center_x(ix)*f,prof->get_bin_center_y(iy)*f,c*f);
						//Ueber->Hist->SetBinContentAt(555,prof->get_bin_center_x(ix),prof->get_bin_center_y(iy),c);
					}
				}
				f = Ueber->proj_Sorter->fw;
				prof = detstruct->pos_walk_calibrator->pos_profile[2][0];
				bx = f*0.5*(prof->center_of_right_bin - prof->center_of_left_bin)/(prof->number_of_columns-1);
				by = f*0.5*(prof->center_of_upper_bin - prof->center_of_lower_bin)/(prof->number_of_rows-1);
				Ueber->Hist->fill(557,"Non-linearity_layer_w",-1000.,-1000.,"x [mm]","dx [mm]",prof->number_of_columns,prof->center_of_left_bin*f-bx,prof->center_of_right_bin*f+bx,prof->number_of_rows,prof->center_of_lower_bin*f-by,prof->center_of_upper_bin*f+by);
				for (__int32 ix = 0;ix < prof->number_of_columns;++ix) {
					for (__int32 iy = 0;iy < prof->number_of_rows;++iy) {
						double c = prof->double_source_matrix[ix][iy];
						if (c != 0.) {
							__int32 sadsad = 0;
						}
						Ueber->Hist->Fill(557,prof->get_bin_center_x(ix)*f,prof->get_bin_center_y(iy)*f,c*f);
						//Ueber->Hist->SetBinContentAt(555,prof->get_bin_center_x(ix),prof->get_bin_center_y(iy),c);
					}
				}
			}
*/

			write_correction_tables(Ueber, det_name, detstruct);
	
			if (detstruct->sorter->scalefactors_calibrator) {
				printf("Hex_%s: Good scale-factors (relative to f_U) are:\nf_V=%lf\nf_W=%lf\nOffset on layer W=%lf\nThese values must be changed in the config-file.\nAfter this please repeat the calibration until the values do not change any more.",det_name.c_str(),detstruct->sorter->scalefactors_calibrator->best_fv,detstruct->sorter->scalefactors_calibrator->best_fw,detstruct->sorter->scalefactors_calibrator->best_w_offset);
				if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"Hex_%s: Good scale-factors (relative to f_U) are:\nf_V=%lf\nf_W=%lf\nOffset on layer W=%lf\nThese values must be changed in the config-file.\nAfter this please repeat the calibration until the values do not change any more.",det_name.c_str(),detstruct->sorter->scalefactors_calibrator->best_fv,detstruct->sorter->scalefactors_calibrator->best_fw,detstruct->sorter->scalefactors_calibrator->best_w_offset); fflush(Ueber->logfile_handle);}
			}
		}
	}

L666:
	// I think this is where we end up after all files were read or IPA is invoked...
	if(Ueber->ipa_enable)
		Ueber->IPA->main_loop();

	if (Ueber->outputRootFile_handle) {

		if(Ueber->CH) {
			printf("\nWriting ColAHelL histograms to the root file ...");
			//write histograms to the root file

			// Presorter histograms first, if we have some...
			int maxHist = Ueber->CH->GetMaxPrsHists();
			if(maxHist>0) {
				histo_handler* hPrs = Ueber->CH->GetPrsHists();
				for ( int j=0;j < maxHist; j++ ) {	
					if(hPrs->H1d_vector[j] != 0){
						Ueber->rt->add_hist(hPrs->H1d_vector[j],Ueber->outputRootFile_handle);
					}
					if(hPrs->H2d_vector[j] != 0){
						Ueber->rt->add_hist(hPrs->H2d_vector[j],Ueber->outputRootFile_handle);
					}
					if(hPrs->H3d_vector[j] != 0){
						Ueber->rt->add_hist(hPrs->H3d_vector[j],Ueber->outputRootFile_handle);
					}
				}
			}
			// ColAHelL histograms...
			maxHist = Ueber->CH->GetMaxHists();
			if(maxHist>0) {
				histo_handler* hCHe = Ueber->CH->GetHists(EL);
				histo_handler* hCHion = Ueber->CH->GetHists(IO);
				histo_handler* hCHuser = Ueber->CH->GetHists(US);

				for ( int j=0;j < maxHist; j++ ) {	
					if(hCHe->H1d_vector[j] != 0)
						Ueber->rt->add_hist(hCHe->H1d_vector[j],Ueber->outputRootFile_handle);
					if(hCHe->H2d_vector[j] != 0)
						Ueber->rt->add_hist(hCHe->H2d_vector[j],Ueber->outputRootFile_handle);
					if(hCHe->H3d_vector[j] != 0)
						Ueber->rt->add_hist(hCHe->H3d_vector[j],Ueber->outputRootFile_handle);
					if(hCHion->H1d_vector[j] != 0)
						Ueber->rt->add_hist(hCHion->H1d_vector[j],Ueber->outputRootFile_handle);
					if(hCHion->H2d_vector[j] != 0)
						Ueber->rt->add_hist(hCHion->H2d_vector[j],Ueber->outputRootFile_handle);
					if(hCHion->H3d_vector[j] != 0)
						Ueber->rt->add_hist(hCHion->H3d_vector[j],Ueber->outputRootFile_handle);
				}
			}
			// User histograms...
			maxHist = Ueber->CH->GetMaxUserHists();
			if(maxHist>0) {
				printf("\nWriting user histograms to the root file ...");
				histo_handler* hCHuser = Ueber->CH->GetHists(US);

				for ( int j=0;j < maxHist; j++ ) {	
					if(hCHuser->H1d_vector[j] != 0)
						Ueber->rt->add_hist(hCHuser->H1d_vector[j],Ueber->outputRootFile_handle);
					if(hCHuser->H2d_vector[j] != 0)
						Ueber->rt->add_hist(hCHuser->H2d_vector[j],Ueber->outputRootFile_handle);
					if(hCHuser->H3d_vector[j] != 0)
						Ueber->rt->add_hist(hCHuser->H3d_vector[j],Ueber->outputRootFile_handle);
				}
			}

		}
		printf("\nwriting root file... ");
		Ueber->outputRootFile_handle->Write();
		printf(" written\n");
	}
	
	if (Ueber->outputRootFile_handle) {
		printf("closing root file... ");
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"closing root file... "); fflush(Ueber->logfile_handle);}
		Ueber->Hist->Reset();
		Ueber->outputRootFile_handle->Close();
		Ueber->outputRootFile_handle = 0;
		printf(" closed\n");
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle," closed\n"); fflush(Ueber->logfile_handle);}
	}
	
	if (Ueber->CH) {delete Ueber->CH; Ueber->CH = 0;} 

	if (Ueber->config_info) {delete Ueber->config_info; Ueber->config_info = 0;}

	if (Ueber->parser) {delete Ueber->parser; Ueber->parser = 0;}

	if (parameter) {delete[] parameter;  parameter = 0;}

	if (Ueber->LMF_output) {
		printf("closing output LMF file... ");
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"closing output LMF file... "); fflush(Ueber->logfile_handle);}
		Ueber->LMF_output->CloseOutputLMF();
		Ueber->LMF_output->CloseInputLMF();
		delete Ueber->LMF_output;   Ueber->LMF_output = 0;
		printf("closed\n");
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"closed\n"); fflush(Ueber->logfile_handle);}
	}
	if (Ueber->rt) {delete Ueber->rt; Ueber->rt = 0;}

	if (Ueber->peak_tracker1) {delete Ueber->peak_tracker1; Ueber->peak_tracker1 = 0;}
	if (Ueber->peak_tracker2) {delete Ueber->peak_tracker2; Ueber->peak_tracker2 = 0;}

	if (Ueber->error) {
		Red(true);
		printf("\nError %i occurred.\n",Ueber->error);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\nError %i occurred.\n",Ueber->error); fflush(Ueber->logfile_handle);}
		Red(false);
	}

	#ifdef dont_use_MFC
		#ifndef LINUX
			printf("\n\nplease check your compiler settings.\nYou have not used the MFC-libraries\n");
			printf("This is only necessary if you are using the (free) Express Version of the Visual Studio\n");
		#endif
	#endif

	if (Ueber->logfile_handle) {
		fclose(Ueber->logfile_handle);
		Ueber->logfile_handle = 0;
	}

	#ifndef dont_use_MFC
		if (Ueber->parallel_sorter) {
			parallel_sort_class * p = (parallel_sort_class *)Ueber->parallel_sorter;
			delete p; Ueber->parallel_sorter = 0;}
	#endif

	for (__int32 detector = 0;detector < 3;++detector) {
		Det_struct * detstruct = 0;
	
		if (detector == 0) {detstruct = Ueber->proj;}
		if (detector == 1) {detstruct = Ueber->rec;}
		if (detector == 2) {detstruct = Ueber->elec;}

		if (!detstruct) continue;

		if (detstruct->sumu_watchdog) {delete detstruct->sumu_watchdog; detstruct->sumu_watchdog = 0;}
		if (detstruct->sumv_watchdog) {delete detstruct->sumv_watchdog; detstruct->sumv_watchdog = 0;}
		if (detstruct->sumw_watchdog) {delete detstruct->sumw_watchdog; detstruct->sumw_watchdog = 0;}
	}

	if (Ueber->proj_Sorter) {delete Ueber->proj_Sorter; Ueber->proj->sorter = Ueber->proj_Sorter = 0;}
	if (Ueber->rec_Sorter ) {delete Ueber->rec_Sorter;  Ueber->rec->sorter  = Ueber->rec_Sorter  = 0;}
	if (Ueber->elec_Sorter) {delete Ueber->elec_Sorter; Ueber->elec->sorter = Ueber->elec_Sorter = 0;}

	if (Ueber->proj) {delete Ueber->proj; Ueber->proj = 0;}
	if (Ueber->rec)  {delete Ueber->rec;  Ueber->rec = 0; }
	if (Ueber->elec) {delete Ueber->elec; Ueber->elec = 0;}

	if (Ueber->Hist) {delete Ueber->Hist; Ueber->Hist = 0;}

	delete_file("bot_config.txt");
	delete_file("I_am_busy.txt");

	if (Ueber->processed_data_file_names) {
		for (__int32 i=0;i<5000;i++) {
			if (Ueber->processed_data_file_names[i]) {
				delete[] Ueber->processed_data_file_names[i]; 
				Ueber->processed_data_file_names[i]=0; 
				continue;
			}
			break;
		}
	}

	if (Ueber->LMF_input) {
		delete Ueber->LMF_input;	Ueber->LMF_input = 0;
	}

	if (Ueber->HDF5_input) {
		delete Ueber->HDF5_input; Ueber->HDF5_input = 0;
	}

	if (Ueber->bot_mode) {
			char command[300];
			sprintf(command,"del %s",Ueber->config_file_name);
			system(command);
			goto L100;
	}

	if (Ueber->processed_data_file_names) {
		for (__int32 i=0;i<5000;i++) {
			if (Ueber->processed_data_file_names[i]) {delete[] Ueber->processed_data_file_names[i]; Ueber->processed_data_file_names[i]=0; continue;}
			break;
		}
		delete[] Ueber->processed_data_file_names;
		Ueber->processed_data_file_names = 0;
	}

	if (Ueber) {delete Ueber; Ueber = 0;}
	delete_file("I_am_alive.txt");
	delete_file("are_you_alive.txt");

#ifdef _DEBUG
	//_CrtDumpMemoryLeaks();
#endif

	return 1;
}










/////////////////////////////////////////////////////////////////////////////
int lmf2root_main(char * args_line)
/////////////////////////////////////////////////////////////////////////////
{
	char * inp = new char[strlen(args_line)+1]; inp[0] = 0;
	char * temp = new char[strlen(args_line)+1]; temp[0] = 0;

	sprintf(inp, args_line);

	char * argv[50];
	for (__int32 i=0;i<50;i++) argv[i] = 0;

	__int32 nof_ars = 0;

	__int32 index = 0;

	for (unsigned __int32 i=0;i<strlen(args_line);i++) {
		if (inp[i] != ' ' && inp[i] != '\t') {
			temp[index] = inp[i];
			temp[index+1] = 0;
			index++;
			continue;
		} else {
			while (true) {
				if (inp[i] == ' ' || inp[i] == '\t') {
					i++;
					continue;
				}
				i--;
				break;
			}
			nof_ars++;
			argv[nof_ars-1] = new char[index+1];
			sprintf(argv[nof_ars-1],temp);
			index = 0;
		}
	}
	if (index > 0) {
		nof_ars++;
		argv[nof_ars-1] = new char[index+1];
		sprintf(argv[nof_ars-1],temp);
	}

	delete[] inp; inp = 0;
	delete[] temp; temp = 0;
	int return_value = lmf2root_main(nof_ars, argv);
	for (__int32 i=0;i<50;i++) {
		if (argv[i]) {delete[] argv[i]; argv[i] = 0;}
	}
	return return_value;
}


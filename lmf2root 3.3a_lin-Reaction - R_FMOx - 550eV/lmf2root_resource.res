        ��  ��                  0� P   I D S _ S O U R C E _ C P P   L M F 2 R O O T _ C P P       0         
#define LMF2ROOTVERSION (3.00)

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
#include "Ueberstruct.h"
#include "../ADC/ADC_meta.h"


#include "../ipa/ipa.h"

//tbb
//#include "Event.h"
//#include "ADC_analysis.h"


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
		printf("\n \n %s does not exist!\n\n",name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s does not exist!\n\n",name); fflush(Ueber->logfile_handle);}
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
			cout << "Output file name is " << root_output_file_name_with_extension << endl;
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


	if (Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT) Ueber->stop_reading_files = false;  	// ADC - lmf2root loophole since no AGAT files are supported..







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
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8HPRAW ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8HP ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_TDC8 ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_2TDC8 ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_RAW32BIT ||
			Ueber->LMF_input->DAQ_ID == DAQ_ID_SIMPLE) {
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

//// ADC - ADC - ADC - ADC - ADC - ADC - ADC - ADC -ADC - ADC -ADC - ADC
		//OutputDebugString("\nLMF_read_event_loop1");

		ADC_meta * adc = (ADC_meta*)Ueber->adc_meta;
					
		if(adc && (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4 || Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8 || Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT)){
			adc->GetTDCandCNTArray(tdc_ns,cnt);
			}

		//if (Ueber->parameter[850] > 0.5) {
		//	std::cout << "\nLMA file mode selected..";
		//}
//// ADC - ADC - ADC - ADC - ADC - ADC - ADC - ADC -ADC - ADC -ADC - ADC

			
	


		// start of fADC8-stuff
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8) {
			printf("please check the source code. This part is commented out.");
			Ueber->stop_reading_files = true;
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
		if (parameter[15] > 0.5) {
			for (__int32 i=0;i<NUM_CHANNELS;++i) {
				for (__int32 j=0;j<cnt[i];++j) {	
					if(i!=Ueber->scan_channel) {
						tdc_ns[i*NUM_CHANNELS + j] += (double(rand())/RAND_MAX-0.5)*Ueber->LMF_input->tdcresolution;
					}
				}
			}
		}

		if (Ueber->stop_reading_files) return true;
		
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
					printf("event %I64i of %I64i\n",Ueber->LMF_input->uint64_number_of_read_events, Ueber->LMF_input->uint64_Numberofevents);
					if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"event %I64i of %I64i\n",Ueber->LMF_input->uint64_number_of_read_events, Ueber->LMF_input->uint64_Numberofevents); fflush(Ueber->logfile_handle);}
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
				if (Ueber->LMF_input->uint64_Numberofevents > 0) {
					stars_in_progress_bar = __int32(52. *Ueber->LMF_input->uint64_number_of_read_events / Ueber->LMF_input->uint64_Numberofevents);
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
								sprintf(name,"%s_%03d.lmf",Ueber->input_lmf_file_name, Ueber->output_LMF_file_index);
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
			if (Ueber->LMF_output && Ueber->LMF_output->DAQ_ID_output == DAQ_ID_FADC8) {
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



				if (Ueber->LMF_output && (Ueber->LMF_output->DAQ_ID_output != DAQ_ID_FADC8) && (Ueber->LMF_output->DAQ_ID_output != DAQ_ID_FADC4)) {
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
							double a = tdc_ns[i*NUM_IONS+j]/Ueber->LMF_input->tdcresolution;
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
	} // end while ...

	return true;
}

void ProcessLMAFile(char *input_LMA_file_name, Ueberstruct * Ueber, double tds_ns[][NUM_IONS])
{
	//lma2root *lma = new lma2root(input_LMA_file_name);
	//lma->ReadLMAHeader();
}















////////////////////////////////////////////////////////////////////
void ProcessLMFFile(char * input_LMF_file_name, Ueberstruct * Ueber)
////////////////////////////////////////////////////////////////////
{
	sprintf(Ueber->input_lmf_file_name, input_LMF_file_name);
	Ueber->reverse_time_direction_in_event = 0;
	if (Ueber->LMF_output_file_name_without_extension[0] == 0) Ueber->output_LMF_file_index = 0;
	double * parameter = Ueber->parameter;
	__int32 * cnt = Ueber->cnt;
	if (Ueber->stop_reading_files) return;

	bool write_sorted_file = ((parameter[49] > 0.5) ? true:false);

	Ueber->LMF_input = new LMF_IO(NUM_CHANNELS,NUM_IONS);
	if (Ueber->LMF_input) {
		if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4) {
			Ueber->LMF_input->CloseInputLMF();
			delete Ueber->LMF_input;	Ueber->LMF_input = 0;
			//extract_mean_pulse_fADC4(input_LMF_file_name, Ueber, tdc_ns); // xxx Robert
			Ueber->LMF_input = new LMF_IO(NUM_CHANNELS,NUM_IONS);
		}
	}
	Ueber->fast_mode = parameter[11] > 0.5 ? 1:0;

	Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions = __int32(parameter[9] + 0.1);
	if (parameter[9] < -0.1) Ueber->LMF_input->TDC8HP.channel_offset_for_rising_transitions--;

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
				return;
			}
		}
	}

	if ( parameter[58] > 0 ) printf("\n\n%G events left (out of %G)\n",parameter[58] - Ueber->eventcounter,parameter[58]);

	printf("\n\nInput file is: %s\n",Ueber->input_lmf_file_name);
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n\nreading input file %s\n",Ueber->input_lmf_file_name); fflush(Ueber->logfile_handle);}
	if(!File_is_readable(Ueber->input_lmf_file_name)) {
		printf("\n \n %s could not be opened!\n\n",Ueber->input_lmf_file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_lmf_file_name); fflush(Ueber->logfile_handle);}
		if (Ueber->LMF_input) {
			Ueber->LMF_input->CloseInputLMF();
			delete Ueber->LMF_input;	Ueber->LMF_input = 0;
		}
		return;
	}

	if (Ueber->LMF_input->DAQ_ID != DAQ_ID_AGAT && !Ueber->LMF_input->OpenInputLMF(Ueber->input_lmf_file_name)) {	// AGAT files are handled separate..
		printf("\n \n %s could not be opened!\n\n",Ueber->input_lmf_file_name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\n \n %s could not be opened!\n\n",Ueber->input_lmf_file_name); fflush(Ueber->logfile_handle);}
		Ueber->stop_reading_files = true;
		if (Ueber->LMF_input) {
			Ueber->LMF_input->CloseInputLMF();
			delete Ueber->LMF_input;	Ueber->LMF_input = 0;
		}
		return;
	}


	if (write_sorted_file) {
		char name[552];
		if (strlen(Ueber->LMF_output_file_name_without_extension) > 0) {
			sprintf(name,"%s_%03d.lmf",Ueber->LMF_output_file_name_without_extension,Ueber->output_LMF_file_index);
		} else {
			sprintf(name,"%s_%03d.lmf",Ueber->input_lmf_file_name,Ueber->output_LMF_file_index);
		}
		printf("\nwriting new LMF-file\n\"%s\"\n",name);
		if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"\nwriting new LMF-file\n\"%s\"\n",name); fflush(Ueber->logfile_handle);}
	}

	Ueber->number_of_channels	= ((__int32(parameter[54]+0.1) == 0) ? (Ueber->LMF_input->number_of_channels + Ueber->LMF_input->number_of_channels2) : __int32(parameter[54]+0.1));
	__int32 number_of_hits		= ((__int32(parameter[55]+0.1) == 0) ? (Ueber->LMF_input->max_number_of_hits) : __int32(parameter[55]+0.1));
	
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


	for (__int32 detector = 0; detector < 3; ++detector) {
		__int32 para_offset;
		sort_class * sorter		= 0;
		Det_struct * det_struct = 0;
		if (detector == 0) {para_offset = 100; sorter = Ueber->proj_Sorter; det_struct = Ueber->proj;}
		if (detector == 1) {para_offset = 200; sorter = Ueber->rec_Sorter;  det_struct = Ueber->rec;}
		if (detector == 2) {para_offset = 300; sorter = Ueber->elec_Sorter; det_struct = Ueber->elec;}
		if (sorter) {
			sorter->TDC_resolution_ns = Ueber->LMF_input->tdcresolution;
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
						return;
					}
					parallel->start_sort_threads();
			} else parameter[10] = 0.;
			}
	}
#endif

	if (parameter[2] > 0.5) {
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
				printf("\n            Events in file  :%11I64i ",Ueber->LMF_input->uint64_Numberofevents);
		#endif
		printf("\n            LMVersionString : %s    ",Ueber->LMF_input->Versionstring.c_str());
	}

	printf("\n0%%           25%%          50%%          75%%        100%%\n");
	printf("|------------|------------|------------|-----------|\n");


	ADC_meta * adc = 0;
	if (Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC4 || Ueber->LMF_input->DAQ_ID == DAQ_ID_FADC8 || Ueber->LMF_input->DAQ_ID == DAQ_ID_AGAT) { // is this a LMF with ADC-data?
		// ADC stuff	 
	  if(!Ueber->adc_meta) Ueber->adc_meta = ADC_meta::instance(Ueber, NUM_CHANNELS, MAX_NBR_PULSES, MAX_NBR_PEAKS);//(void*)new ADC_analysis(Ueber);
			adc = (ADC_meta*)Ueber->adc_meta;


		if (parameter[950] > 0.5 && adc) {
		//	adc->inspectDPA(parameter[951]);
		}	
	}


	OutputDebugString("\ngo into LMF_read_event_loop1");
	LMF_read_event_loop(Ueber);
	
	if (adc) {adc->deleteInstance(); adc = 0;}

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

	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"finished reading LMF-file.\n"); fflush(Ueber->logfile_handle);}

	
/*	if (Ueber->bot_mode) {
		char command[300]; command[0] = 0;
		sprintf(command,"del %s", input_lmf_file_name);
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


	std::string tags_ascii[24] = {"X","Y","TOF","PX","PY","PZ","P","PHI_DEG","PHIPOS","THETA_DEG","ENERGY",
								  "PSUMX","PSUMY","PSUMZ","PSUM","PSUM_PHI_DEG","PSUM_THETA_DEG",
								  "PRELX","PRELY","PRELZ","PREL","PREL_PHI_DEG","PREL_THETA_DEG","KER"};

	Ueberstruct * Ueber = new Ueberstruct();
	Ueber->config_file_name[0] = 0;
	Ueber->logfilepathname[0] = 0;
	Ueber->proj = Ueber->rec = Ueber->elec = 0;
	Ueber->rt = 0;
	Ueber->config_version = 0;
	Ueber->Hist = 0;
	Ueber->current_proj_Sorter = Ueber->current_elec_Sorter = Ueber->current_rec_Sorter = 0;
	Ueber->proj_Sorter = Ueber->rec_Sorter = Ueber->elec_Sorter = 0;
	Ueber->parser = 0;
	Ueber->processed_data_file_names = 0;
	Ueber->input_lmf_file_name[0] = 0;
	Ueber->bot_mode = 0;
	Ueber->scan_channel = -1;
	Ueber->scan_in_data = false;
	Ueber->scan_num_vals = 0;
	Ueber->scan_factor = 1;

	Ueber->adc_meta = 0;
	Ueber->CH = 0;

	//split_pathname(argv[0],Ueber->applicationfilepath,0);
	get_current_directory(Ueber->applicationfilepath, 300);
	
	#ifndef dont_use_MFC
		printf("running in folder\n%s\n",Ueber->applicationfilepath);
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
			printf("\nenter \"lmf2root -help\" to display help the screen\n");
			printf("Now I will try to read \"config.txt\"...\n");
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

		if (Ueber->parser->external_command == "readLMF") {
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
						if (get_next_new_data_file_name(new_name, main_name, Ueber)) new_file_found = true;
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
						ProcessLMFFile((char *)(Ueber->parser->convert_backslash_to_slash(new_name), Ueber);
					#else
						ProcessLMFFile((char *)(new_name), Ueber);
					#endif
					delete_file("I_am_busy.txt");

					if (single_file_mode) break;
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
				Ueber->IPA->enable_wiggletool(200, -10.0, 200.0);
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
// end of interactive paramter adjustment stuff

//--------------------------------------------------------------------------------------
// ColAHelL definitions: ask T.Jahnke for support
//--------------------------------------------------------------------------------------
		if (Ueber->parser->external_command == "use_ColAHelL" || Ueber->parser->external_command == "use_ColaHell" || Ueber->parser->external_command == "use_colahell" ) {
			
			White(true);
			printf("You're using ColAHelL! Good job! ColAHelL is your friend!\n");
			White(false);

			if (Ueber->parser->number_of_external_arguments == 0) {

				Ueber->CH = new CH_Int();
				Ueber->use_CH = true;
				int err = Ueber->CH->ReadConfig();
				if(err<0) {
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
				printf("Reading '%s'...\n\n",Ueber->parser->external_arguments[0]);
				White(true);
				int err = Ueber->CH->ReadConfig(Ueber->parser->external_arguments[0].c_str());
				if(err<0) {
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
	printf("lmf2root was running for %lg seconds.\n",stop_time-start_time);
	if (Ueber->logfile_handle) {fprintf(Ueber->logfile_handle,"lmf2root was running for %lg seconds.\n",stop_time-start_time); fflush(Ueber->logfile_handle);}

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
		printf("\nWriting histograms to the root file ...");
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
			for ( int j=0;j < maxHist; j++ ) {	
				if(hCHe->H1d_vector[j] != 0){
					Ueber->rt->add_hist(hCHe->H1d_vector[j],Ueber->outputRootFile_handle);
				}
				if(hCHe->H2d_vector[j] != 0){
					Ueber->rt->add_hist(hCHe->H2d_vector[j],Ueber->outputRootFile_handle);
				}
				if(hCHe->H3d_vector[j] != 0){
					Ueber->rt->add_hist(hCHe->H3d_vector[j],Ueber->outputRootFile_handle);
				}
				if(hCHion->H1d_vector[j] != 0){
					Ueber->rt->add_hist(hCHion->H1d_vector[j],Ueber->outputRootFile_handle);
				}
				if(hCHion->H2d_vector[j] != 0){
					Ueber->rt->add_hist(hCHion->H2d_vector[j],Ueber->outputRootFile_handle);
				}
				if(hCHion->H3d_vector[j] != 0){
					Ueber->rt->add_hist(hCHion->H3d_vector[j],Ueber->outputRootFile_handle);
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

=" L   I D S _ S O U R C E _ C P P   L M F _ I O _ C P P       0         
#include "LMF_IO.h"

#ifndef LINUX
	#define WINVER 0x0501
	#pragma warning(disable : 4996)
#endif

#define LMF_IO_CLASS_VERSION (2018)


#define DAQ_SOURCE_CODE		0x80000000
#define DAN_SOURCE_CODE		0x40000000
#define CCF_HISTORY_CODE	0x20000000



void MyFILE::seek(unsigned __int64 pos)
{
	__int32 rval = _fseeki64(file,  pos,SEEK_SET);
	if (rval == 0) this->position = pos; else error = 1;


/*
	if (pos < 1000000000) {
		f_stream.seekg((unsigned int)pos, std::ios::beg);
		this->position = pos;
		return;
	}

	unsigned __int64 real_pos = 0;
	f_stream.seekg(0, std::ios::beg);
	while (true) {
		
//		f_stream.seekg(1000000000,std::ios::cur);
		real_pos += 1000000000;
		__int64 diff = pos - real_pos;
		if (diff < 1000000000) break;
	}
	__int64 diff = pos - real_pos;
	__int32 diff1 = __int32(diff);
	f_stream.seekg(diff1,std::ios::cur);*/

}




/*
void MyFILE::seek_to_end()
{
	f_stream.seekg(0, std::ios::end);
}
*/





#define INT_MIN_ (-2147483647-1)




void LMF_IO::write_times(MyFILE * out_file,time_t _Starttime,time_t _Stoptime) {
	unsigned __int32 dummy;
	if (CTime_version_output > 2003) {
		dummy = INT_MIN_+10; *out_file << dummy;
		dummy = (unsigned __int32)_Starttime; *out_file << dummy;
		dummy =0;			*out_file << dummy;
		dummy =INT_MIN_+10;	*out_file << dummy;
		dummy = (unsigned __int32)_Stoptime;	*out_file << dummy;
		dummy = 0;			*out_file << dummy;
	} else {
		dummy = (unsigned __int32)_Starttime;	*out_file << dummy;
		dummy = (unsigned __int32)_Stoptime;	*out_file << dummy;
	}
}









unsigned __int32 ReadCStringLength(MyFILE &in_file)
{
	unsigned __int64 qwLength;
	unsigned __int32 dwLength;
	unsigned __int16 wLength;
	unsigned __int8  bLength;

	__int32 nCharSize = sizeof(__int8);

	// First, try to read a one-byte length
	in_file >> bLength;

	if (bLength < 0xff) 
		return bLength;

	// Try a two-byte length
	in_file >> wLength;
	if (wLength == 0xfffe)
	{
		// Unicode string.  Start over at 1-byte length
		nCharSize = sizeof(wchar_t);

		in_file >> bLength;
		if (bLength < 0xff)
			return bLength;

		// Two-byte length
		in_file >> wLength;
		// Fall through to continue on same branch as ANSI string
	}
	if (wLength < 0xffff)
		return wLength;

	// 4-byte length
	in_file >> dwLength;
	if (dwLength < 0xffffffff)
		return dwLength;

	// 8-byte length
	in_file >> qwLength;

	return (unsigned __int32)qwLength;
}




void Read_CString_as_StdString(MyFILE &in_file,std::string &stdstring)
{
	__int32 length = ReadCStringLength(in_file);
	__int8 * temp_string = new __int8[length+1];
	in_file.read(temp_string,length);
	temp_string[length] = 0;
	stdstring = temp_string;
	delete[] temp_string;
	temp_string = 0;
}



void WriteCStringLength(MyFILE &out_file, unsigned __int32 nLength)
{
	unsigned __int8 dummy_uint8;
	unsigned __int16 dummy_uint16;
	unsigned __int32 dummy_uint32;
	unsigned __int64 dummy_uint64;

	if (nLength < 255)
	{
		dummy_uint8 = nLength; out_file << dummy_uint8;
	}
	else if (nLength < 0xfffe)
	{
		dummy_uint8 = 0xff; out_file << dummy_uint8;
		dummy_uint16 = nLength; out_file << dummy_uint16;
	}
	else if (nLength < 0xffffffff)
	{
		dummy_uint8 = 0xff; out_file << dummy_uint8;
		dummy_uint16 = 0xffff; out_file << dummy_uint16;
		dummy_uint32 = nLength; out_file << dummy_uint32;
	}
	else
	{
		dummy_uint8 = 0xff; out_file << dummy_uint8;
		dummy_uint16 = 0xffff; out_file << dummy_uint16;
		dummy_uint32 = 0xffffffff; out_file << dummy_uint32;
		dummy_uint64 = nLength; out_file << dummy_uint64;
	}
}



void Write_StdString_as_CString(MyFILE &out_file,std::string &stdstring)
{
	unsigned __int32 length = (unsigned __int32)(stdstring.length());
	WriteCStringLength(out_file,length);
	out_file.write(stdstring.c_str(),length);
}







/////////////////////////////////////////////////////////////////
__int32 LMF_IO::GetVersionNumber()
/////////////////////////////////////////////////////////////////
{
	return LMF_IO_CLASS_VERSION;
}











/////////////////////////////////////////////////////////////////
LMF_IO::LMF_IO(__int32 _num_channels, __int32 _num_ions)
/////////////////////////////////////////////////////////////////
{
	num_channels = _num_channels;
	num_ions = _num_ions;
	Initialize();
}






/////////////////////////////////////////////////////////////////
LMF_IO::~LMF_IO()
/////////////////////////////////////////////////////////////////
{
	if (InputFileIsOpen)  CloseInputLMF();
	if (OutputFileIsOpen) CloseOutputLMF();
	if (CAMAC_Data)		{delete[] CAMAC_Data; CAMAC_Data = 0;}
	if (i32TDC)			{delete[] i32TDC; i32TDC = 0;}
	if (us16TDC)		{delete[] us16TDC; us16TDC = 0;}
	if (dTDC)			{delete[] dTDC; dTDC = 0;}
	if (number_of_hits) {delete[] number_of_hits; number_of_hits = 0;}

	for (__int32 i=0;i<3;++i) {
		if (TDC8HP.TDC_info[i]) {
			if (TDC8HP.TDC_info[i]->INLCorrection)	{delete[] TDC8HP.TDC_info[i]->INLCorrection;	TDC8HP.TDC_info[i]->INLCorrection = 0;}
			if (TDC8HP.TDC_info[i]->DNLData)		{delete[] TDC8HP.TDC_info[i]->DNLData;			TDC8HP.TDC_info[i]->DNLData = 0;}
			if (TDC8HP.TDC_info[i]) {delete TDC8HP.TDC_info[i]; TDC8HP.TDC_info[i] = 0;}
		}
	}
	if (ui32buffer) {delete[] ui32buffer; ui32buffer = 0; this->ui32buffer_size = 0;}

	

	if (CCFHistory_strings) {
		for (__int32 i=0;i<number_of_CCFHistory_strings;++i) {
			if (CCFHistory_strings[i]) {delete CCFHistory_strings[i]; CCFHistory_strings[i] = 0;}
		}
		delete[] CCFHistory_strings; CCFHistory_strings = 0;
	}
	if (DAQ_source_strings) {
		for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
			if (DAQ_source_strings[i]) {delete DAQ_source_strings[i]; DAQ_source_strings[i] = 0;}
		}
		delete[] DAQ_source_strings; DAQ_source_strings = 0;
	}
	if (DAN_source_strings) {
		for (__int32 i=0;i<number_of_DAN_source_strings;++i) {
			if (DAN_source_strings[i]) {delete DAN_source_strings[i]; DAN_source_strings[i] = 0;}
		}
		delete[] DAN_source_strings; DAN_source_strings = 0;
	}

	if (CCFHistory_strings_output) {
		for (__int32 i=0;i<number_of_CCFHistory_strings_output;++i) {
			if (CCFHistory_strings_output[i]) {
				delete CCFHistory_strings_output[i]; 
				CCFHistory_strings_output[i] = 0;
			}
		}
		delete[] CCFHistory_strings_output; CCFHistory_strings_output = 0;
	}
	if (DAQ_source_strings_output) {
		for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {
			if (DAQ_source_strings_output[i]) {delete DAQ_source_strings_output[i]; DAQ_source_strings_output[i] = 0;}
		}
		delete[] DAQ_source_strings_output; DAQ_source_strings_output = 0;
	}
	if (DAN_source_strings_output) {
		for (__int32 i=0;i<number_of_DAN_source_strings_output;++i) {
			if (DAN_source_strings_output[i]) {delete DAN_source_strings_output[i]; DAN_source_strings_output[i] = 0;}
		}
		delete[] DAN_source_strings_output; DAN_source_strings_output = 0;
	}

	if (Parameter) {delete[] Parameter; Parameter = 0;}
	if (Parameter_old) {delete[] Parameter_old; Parameter_old = 0;}
}





/////////////////////////////////////////////////////////////////
void LMF_IO::Initialize()
/////////////////////////////////////////////////////////////////
{
	Parameter = new double[10000];
	Parameter_old = new double[10000];
	uint64_LMF_EventCounter = -1;
	not_Cobold_LMF = false;
	errorflag = 0;
	iLMFcompression = 0;
	number_of_channels2 = 0;
	InputFileIsOpen = false;
	OutputFileIsOpen = false;
	SIMPLE_DAQ_ID_Orignial = 0;
	input_lmf = 0;
	output_lmf = 0;
	DAQ_ID = 0;
	frequency = 1.;
	common_mode = 0;
	timerange = 0;
	must_read_first = true;
	uint64_number_of_read_events = 0;
	uint64_number_of_written_events = 0;
	number_of_channels_output = 0;
	number_of_channels2_output = -1;
	number_of_channels = 0;
	max_number_of_hits = 0;
	max_number_of_hits_output = 0;
	max_number_of_hits2_output = -1;
	data_format_in_userheader_output = -2;
	output_byte_counter = 0;
	timestamp_format = 0;
	timestamp_format_output = -1;
	DOUBLE_timestamp = 0.;
	ui64_timestamp = 0;
	skip_header = false;
	time_reference_output = 0;
	system_timeout_output = -1;
	Numberofcoordinates = -2;
	Numberofcoordinates_output = -2;
	DAQ_ID_output = 0;
	User_header_size_output = 0;
	CAMAC_Data = 0;
	changed_mask_read = 0;
	tdcresolution_output = -1.;
	tdcresolution = 0.;
	TDC8HP.SyncValidationChannel = 0;
	TDC8HP.exotic_file_type = 0;
	User_header_size = 0;
	TDC8HP.UserHeaderVersion = 7; // 4 = Cobold 2008 first release in 2008
	                              // 5 = Cobold 2008 R2 in August 2009
								  // 6 = Cobold 2009?
								  // 7 = Cobold 2011 R1 + R2

	HM1.use_normal_method = true;
	TDC8PCI2.use_normal_method_2nd_card = true;
	TDC8PCI2.use_normal_method = true;
	TDC8HP.VHR_25ps = true;
	common_mode_output = -1;
	DAQ_info = "";
	DAQ_info_Length = 0;
	Versionstring = "";
	FilePathName = "";
	OutputFilePathName = "";
	Comment = "";
	Comment_output = "";
	DAQ_info = "";
	Camac_CIF = "";
	TDC8HP.variable_event_length = 0;
	DAQVersion_output = -1;
	LMF_Version = -1;
	LMF_Version_output = -1;
	TDC8HP.csConfigFile = "";
	TDC8HP.csINLFile = "";
	TDC8HP.csDNLFile = "";
	TDC8HP.GroupingEnable_p66		 = false;
	TDC8HP.GroupingEnable_p66_output = false;

	DAQ_SOURCE_CODE_bitmasked  = 0;
	DAN_SOURCE_CODE_bitmasked  = 0;
	CCF_HISTORY_CODE_bitmasked = 0;
	number_of_CCFHistory_strings = 0;
	number_of_DAN_source_strings = 0;
	number_of_DAQ_source_strings = 0;
	number_of_CCFHistory_strings_output = 0;
	number_of_DAN_source_strings_output = 0;
	number_of_DAQ_source_strings_output = 0;
	CCFHistory_strings = 0;
	DAN_source_strings = 0;
	DAQ_source_strings = 0;
	CCFHistory_strings_output = 0;
	DAN_source_strings_output = 0;
	DAQ_source_strings_output = 0;

	TDC8HP.number_of_bools   = 0;
	TDC8HP.number_of_doubles = 0;
	TDC8HP.number_of_doubles = 0;

	TDC8PCI2.variable_event_length = 0;

	i32TDC =	new __int32[num_channels*num_ions];
	us16TDC =	new unsigned __int16[num_channels*num_ions];
	dTDC =		new double [num_channels*num_ions];
	number_of_hits = new unsigned __int32[num_channels];

	memset(number_of_hits,0,num_channels*sizeof(__int32));

	TDC8HP.DMAEnable = true;
	TDC8HP.SSEEnable = false;
	TDC8HP.MMXEnable = true;
	TDC8HP.GroupTimeOut = 0.1;

	TDC8HP.channel_offset_for_rising_transitions = 0;

//	TDC8HP.bdummy = false;
//	TDC8HP.idummy = 0;
//	TDC8HP.ddummy = 0.;

	TDC8HP.i32NumberOfDAQLoops  = 1;
	TDC8HP.TDC8HP_DriverVersion = 0x00000000;
	TDC8HP.iTriggerChannelMask  = 0;
	TDC8HP.iTime_zero_channel   = 0;

	TDC8HP.Number_of_TDCs = 0;

	time_t osBinaryTime;
	time( &osBinaryTime );
	time_t time_dummy(osBinaryTime);
	Starttime = time_dummy;
	Stoptime = time_dummy;
	Starttime_output = 0;
	Stoptime_output  = 0;
	CTime_version    = 0;
	CTime_version_output = 0;
	number_of_DAQ_source_strings_output = -1;
	number_of_CCFHistory_strings_output = -1;
	number_of_DAN_source_strings_output = -1;

	number_of_bytes_in_PostEventData = 0;

	LMF_Header_version = 476759;
	ui64LevelInfo = 0;

	Cobold_Header_version = 2002;
	Cobold_Header_version_output = 0;

	TDC8HP.OffsetTimeZeroChannel_s = 0.;
	TDC8HP.BinsizeType = 0;

	for (__int32 i=0;i<3;++i) {
		TDC8HP.TDC_info[i] = 0;
		TDC8HP.TDC_info[i] = new TDC8HP_info_struct;
		TDC8HP.TDC_info[i]->INLCorrection = 0;
		TDC8HP.TDC_info[i]->INLCorrection	= new __int32[8*1024];
		TDC8HP.TDC_info[i]->DNLData = 0;
		TDC8HP.TDC_info[i]->DNLData			= new unsigned __int16[8*1024];
	}



	fADC8.driver_version = 0;
	fADC8.bReadCustomData = false;
	fADC8.i32NumberOfDAQLoops = 0;
	fADC8.number_of_bools = 0;
	fADC8.number_of_int32s = 0;
	fADC8.number_of_uint32s = 0;
	fADC8.number_of_doubles = 0;
	fADC8.GroupEndMarker = 0;
	fADC8.i32NumberOfADCmodules = 0;
	fADC8.iEnableGroupMode = 1;
	fADC8.iTriggerChannel = 0;
	fADC8.iPreSamplings_in_4800ps_units = 0;
	fADC8.iPostSamplings_in_9600ps_units = 0;
	fADC8.iEnableTDCinputs = 0;
	fADC8.veto_gate_length = 0;
	fADC8.veto_delay_length = 0;
	fADC8.veto_mask = 0;
	fADC8.dGroupRangeStart = 0.;
	fADC8.at_least_1_signal_was_written = false;
	fADC8.dGroupRangeEnd = 0.;
	for (__int32 mod = 0;mod<8;mod++) {
		fADC8.firmware_version[mod]	= 0;
		fADC8.serial_number[mod]		= 0;
		for (__int32 ch = 0;ch<8;ch++) {
			fADC8.GND_level[mod][ch] = 0;
			fADC8.iThreshold_GT[mod][ch] = 0;
			fADC8.iThreshold_LT[mod][ch] = 0;
		}
		fADC8.iChannelMode [mod][0] = fADC8.iChannelMode [mod][1] = 0;
		fADC8.iSynchronMode[mod][0] = fADC8.iSynchronMode[mod][1] = 0;
		fADC8.dSyncTimeOffset[mod][0] = fADC8.dSyncTimeOffset[mod][1] = 0.;
	}

///////////////////////////////////


	fADC4.packet_count = -1;
	fADC4.number_of_bools = 0;
	fADC4.number_of_int32s = 0;
	fADC4.number_of_uint32s = 0;
	fADC4.number_of_doubles = 0;
	fADC4.GroupEndMarker = 12345678;
	fADC4.driver_version = 0;
	fADC4.i32NumberOfDAQLoops = 1;
	fADC4.bReadCustomData = true;
	fADC4.i32NumberOfADCmodules = 0;
	fADC4.iTriggerChannel = 0;
	fADC4.dGroupRangeStart = 0.;
	fADC4.dGroupRangeEnd = 0.;
	fADC4.csConfigFile = fADC4.csINLFile = fADC4.csDNLFile = "";
	for (__int32 i=0;i<20;i++) fADC4.bits_per_mVolt[i] = 0.;
	for (__int32 i=0;i<80;i++) fADC4.GND_level[i] = 0.;

	ui32buffer_size = 0;
	ui32buffer = 0;

	error_text[0] = (char*)"no error";
	error_text[1] = (char*)"error reading timestamp";
	error_text[2] = (char*)"error reading data";
	error_text[3] = (char*)"input file is already open";
	error_text[4] = (char*)"could not open input file";
	error_text[5] = (char*)"could not connect CAchrive to input file";
	error_text[6] = (char*)"error reading header";
	error_text[7] = (char*)"LMF not data of TDC8PCI2 or 2TDC8PCI2 or TDC8HP or CAMAC";
	error_text[8] = (char*)"file format not supported (only unsigned __int16 16bit and signed integer 32)";
	error_text[9] = (char*)"input file not open";
	error_text[10] = (char*)"output file not open";
	error_text[11] = (char*)"could not open output file";
	error_text[12] = (char*)"output file is already open";
	error_text[13] = (char*)"could not connect CAchrive to output file";
	error_text[14] = (char*)"some parameters are not initialized";
	error_text[15] = (char*)"CAMAC data tried to read with wrong function";
	error_text[16] = (char*)"seek does not work with non-fixed event lengths";
	error_text[17] = (char*)"writing file with non-fixed event length dan DAQVersion < 2008 no possible";
	error_text[18] = (char*)"end of input file";
	error_text[19] = (char*)"more channels in file than specified at new LMF_IO()";
	error_text[20] = (char*)"more hits per channel in file than specified at new LMF_IO()";
	error_text[21] = (char*)"more bytes after event than are reserved in LMF_IO source code";
}






/////////////////////////////////////////////////////////////////
void LMF_IO::CloseInputLMF()
/////////////////////////////////////////////////////////////////
{
	if (input_lmf)	{input_lmf->close(); delete input_lmf; input_lmf = 0;}
	InputFileIsOpen = false;
}







/////////////////////////////////////////////////////////////////
unsigned __int64 LMF_IO::GetLastLevelInfo()
/////////////////////////////////////////////////////////////////
{
	return ui64LevelInfo;
}








/////////////////////////////////////////////////////////////////
bool LMF_IO::OpenNonCoboldFile(void)
/////////////////////////////////////////////////////////////////
{
	input_lmf->seek(0);

	if (skip_header) {
		DAQ_ID = DAQ_ID_RAW32BIT;
		data_format_in_userheader = 10;
		return true;
	}

	if (DAQ_ID == DAQ_ID_RAW32BIT) return true;

	DAQ_ID = DAQ_ID_SIMPLE;

	*input_lmf >> SIMPLE_DAQ_ID_Orignial;
	unsigned tempi;
	*input_lmf >> tempi; uint64_Numberofevents = tempi;
	*input_lmf >> data_format_in_userheader;

	User_header_size = 3*sizeof(__int32);
	Headersize = 0;

	return true;
}








/////////////////////////////////////////////////////////////////
bool LMF_IO::OpenInputLMF(std::string Filename)
/////////////////////////////////////////////////////////////////
{
	return OpenInputLMF((__int8*)Filename.c_str());
}







/////////////////////////////////////////////////////////////////
bool LMF_IO::OpenInputLMF(__int8 * LMF_Filename)
/////////////////////////////////////////////////////////////////
{
	unsigned __int32	unsigned_int_Dummy;
	__int32				data_format_in_header;
	__int8				byte_Dummy;
	__int32				byte_counter;

	if (Parameter) memset(Parameter,0,10000*sizeof(double));
	if (Parameter_old) memset(Parameter_old,0,10000*sizeof(double));
	for (__int32 i=901;i<=932;i++) Parameter_old[i] = -1.e201;

	if (InputFileIsOpen) {
		errorflag = 3; // file is already open
		return false;
	}
	input_lmf = new MyFILE(true);

	TDC8HP.UserHeaderVersion = 0; // yes, 0 is ok here and 2 in LMF_IO::initialization is also ok

	input_lmf->open(LMF_Filename);

	if (input_lmf->error) {
		errorflag = 4; // could not open file
		input_lmf = 0;
		return false;
	}

//L10:

//  READ LMF-HEADER
	ArchiveFlag = 0;

	*input_lmf >> ArchiveFlag;

	DAQ_SOURCE_CODE_bitmasked  =  ArchiveFlag & DAQ_SOURCE_CODE;
	DAN_SOURCE_CODE_bitmasked  =  ArchiveFlag & DAN_SOURCE_CODE;
	CCF_HISTORY_CODE_bitmasked =  ArchiveFlag & CCF_HISTORY_CODE;

	ArchiveFlag = ArchiveFlag & 0x1fffffff;
	
	if (ArchiveFlag != 476758 && ArchiveFlag != 476759) { // is this not a Cobold list mode file?
		not_Cobold_LMF = true;
		if (!OpenNonCoboldFile()) return false;
		errorflag = 0; // no error
		InputFileIsOpen = true;
		return true;
	}

	if (ArchiveFlag == 476758) Cobold_Header_version = 2002;
	if (ArchiveFlag == 476759) Cobold_Header_version = 2008;

	*input_lmf >> data_format_in_header;

	data_format_in_userheader = data_format_in_header;
	if (data_format_in_header != LM_USERDEF && data_format_in_header != LM_SHORT && data_format_in_header != LM_SLONG && data_format_in_header != LM_DOUBLE && data_format_in_header != LM_CAMAC) {
		errorflag = 8;
		CloseInputLMF();
		return false;
	}

	if (Cobold_Header_version <= 2002)	*input_lmf >> Numberofcoordinates;
	if (Cobold_Header_version >= 2008)	{
		unsigned __int64 temp;
		*input_lmf >> temp;
		Numberofcoordinates =__int32(temp);
	}
	
	if (Numberofcoordinates >= 0) CAMAC_Data = new unsigned __int32[Numberofcoordinates];

	if (Cobold_Header_version <= 2002)	*input_lmf >> Headersize;
	if (Cobold_Header_version >= 2008)	{
		unsigned __int64 temp;
		*input_lmf >> temp;
		Headersize =__int32(temp);
	}


	if (Cobold_Header_version <= 2002)	*input_lmf >> User_header_size;
	if (Cobold_Header_version >= 2008)	{
		unsigned __int64 temp;
		*input_lmf >> temp;
		User_header_size =__int32(temp);
	}

	if (skip_header) {
		__int32 backstep;
		if (Cobold_Header_version <= 2002) backstep = sizeof(__int32)*5;
		if (Cobold_Header_version >= 2008) backstep = sizeof(__int32)*2 + sizeof(__int64)*3;
		for (unsigned __int32 i=0;i<Headersize-backstep;++i) *input_lmf >> byte_Dummy;
		for (unsigned __int32 i=0;i<User_header_size;++i) *input_lmf >> byte_Dummy;
		errorflag = 0; // no error
		InputFileIsOpen = true;
		goto L666;
	}

	if (Cobold_Header_version >= 2008) *input_lmf >> uint64_Numberofevents;
	if (Cobold_Header_version <= 2002) {
		__int32 temp;
		*input_lmf >> temp;
		uint64_Numberofevents = temp;
	}

	// get CTime version:
	if (!CTime_version) {
		unsigned __int32 dummy_uint32;
		unsigned __int64 pos = (unsigned __int64)(input_lmf->tell()); 
		*input_lmf >> dummy_uint32;
		if (dummy_uint32 == INT_MIN_+10) CTime_version = 2005; else CTime_version = 2003;
		input_lmf->seek(pos);
	}
	
	CTime_version_output = CTime_version;

	unsigned __int32 dummy_uint32;
	if(CTime_version >= 2005) {
		*input_lmf >> dummy_uint32;
		*input_lmf >> dummy_uint32; Starttime = dummy_uint32;
		*input_lmf >> dummy_uint32;
	} else {*input_lmf >> dummy_uint32; Starttime = dummy_uint32;}

	if(CTime_version >= 2005) {
		*input_lmf >> dummy_uint32;
		*input_lmf >> dummy_uint32; Stoptime = dummy_uint32;
		*input_lmf >> dummy_uint32;
	} else {
		*input_lmf >> dummy_uint32; Stoptime = dummy_uint32;
	}

	Read_CString_as_StdString(*input_lmf,Versionstring);
	Read_CString_as_StdString(*input_lmf,FilePathName);
	Read_CString_as_StdString(*input_lmf,Comment);
	Comment_output = Comment;

	byte_counter = 0;

	if (CCF_HISTORY_CODE_bitmasked) {
		*input_lmf >> number_of_CCFHistory_strings;
		CCFHistory_strings = new std::string*[number_of_CCFHistory_strings];
		memset(CCFHistory_strings,0,sizeof(std::string*)*number_of_CCFHistory_strings);
		for (__int32 i=0;i<number_of_CCFHistory_strings;++i) {
			__int32 string_len;
			*input_lmf >> string_len;
			CCFHistory_strings[i] = new std::string();
			CCFHistory_strings[i]->reserve(string_len);
			while (string_len > 0) {
				__int8 c;
				*input_lmf >> c;	
				*CCFHistory_strings[i] += c;
				--string_len;
			}
		}
	}

	if (DAN_SOURCE_CODE_bitmasked) {
		*input_lmf >> number_of_DAN_source_strings;
		DAN_source_strings = new std::string*[number_of_DAN_source_strings];
		memset(DAN_source_strings,0,sizeof(std::string*)*number_of_DAN_source_strings);
		for (__int32 i=0;i<number_of_DAN_source_strings;++i) {
			__int32 string_len;
			*input_lmf >> string_len;
			DAN_source_strings[i] = new std::string();
			DAN_source_strings[i]->reserve(string_len);
			while (string_len > 0) {
				__int8 c;
				*input_lmf >> c; 
				*DAN_source_strings[i] += c;
				--string_len;
			}
		}
	}	



	if (User_header_size == 0) {
		errorflag = 6;
		InputFileIsOpen = false;
		goto L666;
	}

	if (__int32(input_lmf->tell()) != __int32(Headersize)) {
		errorflag = 6;
		InputFileIsOpen = false;
		goto L666;
	}



//  READ USER-HEADER
	if (Cobold_Header_version >= 2008) {
		*input_lmf >> LMF_Header_version;   byte_counter += sizeof(__int32);
	}

	if (Cobold_Header_version <= 2002) {*input_lmf >> unsigned_int_Dummy;	byte_counter += sizeof(__int32);}

	if (Cobold_Header_version >= 2008) {
		unsigned __int64 temp;
		*input_lmf >> temp;	unsigned_int_Dummy = (unsigned __int32)(temp);
		byte_counter += sizeof(unsigned __int64);
	}


	if (unsigned_int_Dummy != User_header_size) {
		errorflag = 6; // error reading header
		CloseInputLMF();
		return false;
	}

	*input_lmf >> DAQVersion;	byte_counter += sizeof(__int32);	// Version is always 2nd value

	*input_lmf >> DAQ_ID;		byte_counter += sizeof(__int32);	// DAQ_ID is always 3ed value

	if (DAQ_ID == DAQ_ID_TDC8)		goto L100;
	if (DAQ_ID == DAQ_ID_2TDC8)		goto L100;
	if (DAQ_ID == DAQ_ID_TDC8HP)	goto L100;
	if (DAQ_ID == DAQ_ID_TDC8HPRAW)	goto L100;
	if (DAQ_ID == DAQ_ID_HM1)		goto L100;
	if (DAQ_ID == DAQ_ID_CAMAC)		goto L100;
	if (DAQ_ID == DAQ_ID_HM1_ABM)	goto L100;
	if (DAQ_ID == DAQ_ID_TCPIP)		goto L100;
	if (DAQ_ID == DAQ_ID_FADC8)		goto L100;
	if (DAQ_ID == DAQ_ID_FADC4)		goto L100;

	errorflag = 7; // LMF not data of TDC8PCI2 or 2TDC8PCI2 or TDC8HP or CAMAC or HM1 or TCPIP
	CloseInputLMF();
	return false;

L100:
	if (DAQ_ID == DAQ_ID_TDC8)	 byte_counter += ReadTDC8PCI2Header();
	if (DAQ_ID == DAQ_ID_2TDC8)	 byte_counter += Read2TDC8PCI2Header();
	if (DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) byte_counter += ReadTDC8HPHeader_LMFV_1_to_7(byte_counter);
	if (DAQ_ID == DAQ_ID_HM1 || DAQ_ID == DAQ_ID_HM1_ABM)     byte_counter += ReadHM1Header();
	if (DAQ_ID == DAQ_ID_CAMAC)  byte_counter += ReadCAMACHeader();
	if (DAQ_ID == DAQ_ID_TCPIP)  byte_counter += ReadTCPIPHeader();
	if (DAQ_ID == DAQ_ID_FADC8)  byte_counter += ReadfADC8Header();
	if (DAQ_ID == DAQ_ID_FADC4)  byte_counter += ReadfADC4Header_up_to_v11();

	if ((__int32(User_header_size) != byte_counter) || (data_format_in_userheader != data_format_in_header)) {
		if (!(DAQ_ID == DAQ_ID_TDC8 && this->LMF_Version == 0x8)
			&& !(this->LMF_Version == 7 && (DAQ_ID == DAQ_ID_HM1 || DAQ_ID == DAQ_ID_HM1_ABM) && DAQVersion == 20080507 && User_header_size > 100000)) {
			errorflag = 6; // error reading header
			CloseInputLMF();
			return false;
		}
	}

	if (data_format_in_header != LM_USERDEF && data_format_in_userheader != LM_SHORT && data_format_in_userheader != LM_SLONG && data_format_in_userheader != LM_DOUBLE && data_format_in_userheader != LM_CAMAC) {
		errorflag = 6; // error reading header
		CloseInputLMF();
		return false;
	}

	errorflag = 0; // no error
	InputFileIsOpen = true;

L666:

	return true;
}












/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteTDC8PCI2Header()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;
	__int32 int_Dummy = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	int_Dummy = __int32(DAQ_info.length());
	*output_lmf << int_Dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	if ((DAQVersion_output >= 20020408 && TDC8PCI2.use_normal_method)) {
		*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);
	}

	if (LMF_Version_output >= 0x9) {
		if ((number_of_DAQ_source_strings_output < 0) || (!DAQ_source_strings_output)) number_of_DAQ_source_strings_output = 0;
		*output_lmf << number_of_DAQ_source_strings_output;  byte_counter += sizeof(__int32);

		for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {
			unsigned __int32 unsigned_int_Dummy = (unsigned __int32)(DAQ_source_strings_output[i]->length());
			*output_lmf << unsigned_int_Dummy;   byte_counter += sizeof(__int32);
			output_lmf->write(DAQ_source_strings_output[i]->c_str(),__int32(DAQ_source_strings_output[i]->length()));
			byte_counter += (unsigned __int32)(DAQ_source_strings_output[i]->length());
		}
	}

	*output_lmf << system_timeout_output; byte_counter += sizeof(__int32);		//   system time-out
	*output_lmf << time_reference_output; byte_counter += sizeof(__int32);
	*output_lmf << common_mode_output; byte_counter += sizeof(__int32);		//   0 common start    1 common stop
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	if (DAQVersion_output >= 20020408 && TDC8PCI2.use_normal_method) {*output_lmf << TDCDataType; byte_counter += sizeof(__int32);}
	*output_lmf << timerange; byte_counter += sizeof(__int32);	// time range of the tdc in microseconds

	if (DAQVersion_output < 20080507) {
		*output_lmf << number_of_channels_output; byte_counter += sizeof(__int32);			// number of channels
		*output_lmf << max_number_of_hits_output; byte_counter += sizeof(__int32);			// number of hits
	} else {
		__int64 i64_temp = number_of_channels_output;
		*output_lmf << i64_temp; byte_counter += sizeof(__int64);			// number of channels
		i64_temp = max_number_of_hits_output;
		*output_lmf << i64_temp; byte_counter += sizeof(__int64);			// number of hits
	}

	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	*output_lmf << module_2nd;	byte_counter += sizeof(__int32);	// indicator for 2nd module data
	
	if (DAQVersion_output >= 20020408 && TDC8PCI2.use_normal_method && (DAQ_ID_output == DAQ_ID_TDC8)) {
		*output_lmf << TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
		*output_lmf << TDC8PCI2.OpenTime_1st_card;			byte_counter += sizeof(__int32); // open time 1st card
		*output_lmf << TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
		*output_lmf << TDC8PCI2.TriggerFalling_1st_card;	byte_counter += sizeof(__int32); // trigger falling edge 1st card
		*output_lmf << TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card
		*output_lmf << TDC8PCI2.EmptyCounter_1st_card;		byte_counter += sizeof(__int32); // EmptyCounter 1st card
		*output_lmf << TDC8PCI2.EmptyCounter_since_last_Event_1st_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 1st card
	} 

	return byte_counter;
}




















/////////////////////////////////////////////////////////////////
void LMF_IO::CloseOutputLMF()
/////////////////////////////////////////////////////////////////
{
	if (!output_lmf) return;

	if (DAQ_ID_output == DAQ_ID_RAW32BIT) {
		output_lmf->close(); delete output_lmf; output_lmf = 0;
		OutputFileIsOpen = false;
		return;
	}

	if (Cobold_Header_version_output == 0) Cobold_Header_version_output = Cobold_Header_version;
	
	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		output_lmf->seek(0);

		*output_lmf << SIMPLE_DAQ_ID_Orignial;
		unsigned __int32 dummy = __int32(uint64_number_of_written_events);
		*output_lmf << dummy;
		*output_lmf << data_format_in_userheader_output;

		output_lmf->close(); delete output_lmf; output_lmf = 0;
		OutputFileIsOpen = false;
		return;
	}

//	if (out_ar) {
		output_lmf->seek(0);

		WriteFirstHeader();

		output_lmf->flush();
		Headersize_output = (unsigned __int32)(output_lmf->tell());
		unsigned __int64 seek_value;
		if (Cobold_Header_version_output <= 2002) seek_value = 3*sizeof(unsigned __int32);
		if (Cobold_Header_version_output >= 2008) seek_value = 2*sizeof(unsigned __int32) + sizeof(unsigned __int64);

		output_lmf->seek(seek_value);
		if (Cobold_Header_version_output <= 2002) *output_lmf << Headersize_output;
		if (Cobold_Header_version_output >= 2008) {
			unsigned __int64 temp = Headersize_output;
			*output_lmf << temp;
		}
		output_lmf->flush();
		output_lmf->seek(Headersize_output);

		if (Cobold_Header_version_output >= 2008 || DAQVersion_output >= 20080000) {
			*output_lmf << LMF_Header_version;
		}

		if (Cobold_Header_version_output <= 2002) *output_lmf << User_header_size_output;
		if (Cobold_Header_version_output >= 2008) {
			unsigned __int64 temp = User_header_size_output;
			*output_lmf << temp;
		}
		
		//out_ar->Close(); out_ar=0;
//	}
	
	if (output_lmf) {output_lmf->close(); delete output_lmf; output_lmf = 0;}
	OutputFileIsOpen = false;
}











/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadTDC8PCI2Header()
/////////////////////////////////////////////////////////////////
{
	__int32 byte_counter;
	byte_counter = 0;



	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit
	
	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;

//	input_lmf->flush();
	unsigned __int64 StartPosition = input_lmf->tell();
	__int32 old_byte_counter = byte_counter;
L50:
	byte_counter = old_byte_counter;

//	TRY
	if (DAQVersion >= 20020408 && TDC8PCI2.use_normal_method) {
		*input_lmf >> LMF_Version;
		byte_counter += sizeof(__int32);
	}

	if (DAQVersion >= 20080507) {
		if (LMF_Version >= 0x8 && data_format_in_userheader == -1) TDC8PCI2.variable_event_length = 1;
		if (LMF_Version >= 0x8) {
			*input_lmf >> number_of_DAQ_source_strings;    byte_counter += sizeof(__int32);
			DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
			memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
			for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
				__int32 string_len;
				*input_lmf >> string_len;    byte_counter += sizeof(__int32);
				DAQ_source_strings[i] = new std::string();
				DAQ_source_strings[i]->reserve(string_len);
				while (string_len > 0) {
					__int8 c;
					*input_lmf >> c;     byte_counter += sizeof(__int8);
					*DAQ_source_strings[i] += c;
					--string_len;
				}
			}
		}
	}

	*input_lmf >> system_timeout;	byte_counter += sizeof(__int32);		//   system time-out
	*input_lmf >> time_reference;	byte_counter += sizeof(__int32);
	*input_lmf >> common_mode;		byte_counter += sizeof(__int32);			//   0 common start    1 common stop
	*input_lmf >> tdcresolution;	byte_counter += sizeof(double);	// tdc resolution in ns


	TDCDataType = 1;
	if (DAQVersion >= 20020408 && TDC8PCI2.use_normal_method) {*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);}
	*input_lmf >> timerange; byte_counter += sizeof(__int32);			// time range of the tdc in microseconds

	if (this->DAQVersion < 20080507) {
		*input_lmf >> number_of_channels; byte_counter += sizeof(__int32);	// number of channels
		*input_lmf >> max_number_of_hits; byte_counter += sizeof(__int32);	// number of hits
	} else {
		__int64 i64_temp;
		*input_lmf >> i64_temp; number_of_channels = (unsigned __int32)(i64_temp); byte_counter += sizeof(__int64);	// number of channels
		*input_lmf >> i64_temp; max_number_of_hits = (unsigned __int32)(i64_temp); byte_counter += sizeof(__int64);	// number of hits
	}

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}
	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	*input_lmf >> module_2nd;	byte_counter += sizeof(__int32);		// indicator for 2nd module data
	
	if (byte_counter == __int32(User_header_size - 12)) return byte_counter;

	if (DAQVersion >= 20020408 && TDC8PCI2.use_normal_method) {
		*input_lmf >> TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
		*input_lmf >> TDC8PCI2.OpenTime_1st_card;			byte_counter += sizeof(__int32); // open time 1st card
		*input_lmf >> TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
		*input_lmf >> TDC8PCI2.TriggerFalling_1st_card;		byte_counter += sizeof(__int32); // trigger falling edge 1st card
		*input_lmf >> TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card
		*input_lmf >> TDC8PCI2.EmptyCounter_1st_card;		byte_counter += sizeof(__int32); // EmptyCounter 1st card
		*input_lmf >> TDC8PCI2.EmptyCounter_since_last_Event_1st_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 1st card
	}

/*
	CATCH( CArchiveException, e )
		if (!TDC8PCI2.use_normal_method) return 0;
		TDC8PCI2.use_normal_method = false;
		in_ar->Close(); delete in_ar; in_ar = 0;
		input_lmf->seek(StartPosition);
		in_ar = new CArchive(input_lmf,CArchive::load);
		goto L50;
	END_CATCH
*/

	if (byte_counter != __int32(User_header_size - 12) && DAQVersion < 20080507) {
		if (!TDC8PCI2.use_normal_method) return 0;
		TDC8PCI2.use_normal_method = false;
		input_lmf->seek(StartPosition);
		goto L50;
	}

	if (LMF_Version == 0x8) {
		input_lmf->flush();
		input_lmf->seek((unsigned __int64)(this->Headersize+this->User_header_size));
		input_lmf->flush();
	}

	return byte_counter;
}








/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadCAMACHeader()
/////////////////////////////////////////////////////////////////
{
	__int32 byte_counter;
	byte_counter = 0;

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;

	*input_lmf >> Camac_CIF_Length;	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[Camac_CIF_Length+1];
	input_lmf->read(__int8_temp,Camac_CIF_Length);
	__int8_temp[Camac_CIF_Length] = 0;
	Camac_CIF = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += Camac_CIF_Length;

	*input_lmf >> system_timeout; byte_counter += sizeof(__int32);		// system time-out
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);		// data format (2=short integer)

	return byte_counter;
}





/////////////////////////////////////////////////////////////////
__int32	LMF_IO::Read2TDC8PCI2Header()
/////////////////////////////////////////////////////////////////
{
	unsigned __int64 StartPosition;
	__int32 old_byte_counter;
	bool desperate_mode;

	TDC8PCI2.variable_event_length = 0;
	__int32 byte_counter;
	byte_counter = 0;

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;

	*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	*input_lmf >> system_timeout; byte_counter += sizeof(__int32);		// system time-out
	
	if (LMF_Version >= 9 && DAQVersion >= 20110208) {				// handle new CStringArray information
		unsigned __int32 ui32CStringCount = system_timeout;			// # of defined CStrings

		__int32 iDummy;
		__int8 cDummy;
		for(unsigned __int32 ui32Count=0;ui32Count<ui32CStringCount;ui32Count++) {	// now read every CString
			*input_lmf >> iDummy; byte_counter += sizeof(__int32);					// how many characters to read
			for(unsigned __int32 ui32Count2=0;ui32Count2<(unsigned int)iDummy;ui32Count2++)		// read character by character
				input_lmf->read(&cDummy,1);							// skip over DAq-source code
			byte_counter += iDummy;
		}
		*input_lmf >> system_timeout; byte_counter += sizeof(__int32);		// again system time-out
	}
	
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> common_mode; byte_counter += sizeof(__int32);			// 0 common start    1 common stop
	*input_lmf >> tdcresolution; byte_counter += sizeof(double);	// tdc resolution in ns

	*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);
	*input_lmf >> timerange; byte_counter += sizeof(__int32);			// time range of the tdc in microseconds
	if (DAQVersion >= 20080507 && DAQVersion < 20110208) {
		TDC8PCI2.variable_event_length = 1;
		__int32 iDummy;
		unsigned __int64 i64temp;
		*input_lmf >> i64temp;		number_of_channels =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of channels
		*input_lmf >> i64temp;		max_number_of_hits =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of hits
		if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
		if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

		*input_lmf >> i64temp;		number_of_channels2 =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of channels2
		*input_lmf >> i64temp;		max_number_of_hits2 =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of hits2
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);

		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);

		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);		// Sync Mode    (parameter 60)
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);		// IO address 2 (parameter 61)
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);
		*input_lmf >> iDummy;			byte_counter += sizeof(__int32);

		goto L200;
	}

	if (DAQVersion >= 20110208) {
		TDC8PCI2.variable_event_length = 1;
		__int32 iDummy;
		unsigned __int64 i64temp;
		*input_lmf >> i64temp;		number_of_channels =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of channels
		*input_lmf >> i64temp;		max_number_of_hits =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of hits
		if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
		if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

		*input_lmf >> i64temp;		number_of_channels2 =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of channels2
		*input_lmf >> i64temp;		max_number_of_hits2 =__int32(i64temp);	byte_counter += sizeof(unsigned __int64);			// number of hits2
		
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // indicator 2nd module
		
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
		
		
		
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // gate delay
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // gate open
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // write empty events
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // trigger at falling edge
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // trigger at rising edge
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
		
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32);

		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // gate delay
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // gate open
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // write empty events
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // trigger at falling edge
		*input_lmf >> iDummy;		byte_counter += sizeof(__int32); // trigger at rising edge
		if (byte_counter < __int32(User_header_size) - 20) {
			*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
			*input_lmf >> iDummy;		byte_counter += sizeof(__int32);
		}

		goto L200;
	}
	

	*input_lmf >> number_of_channels; byte_counter += sizeof(__int32);	// number of channels
	*input_lmf >> max_number_of_hits; byte_counter += sizeof(__int32);	// number of hits
	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}
	*input_lmf >> number_of_channels2; byte_counter += sizeof(__int32);	// number of channels2
	*input_lmf >> max_number_of_hits2; byte_counter += sizeof(__int32);	// number of hits2
	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);				// data format (2=short integer)

//	input_lmf->flush();
	StartPosition = input_lmf->tell();
	old_byte_counter = byte_counter;
	desperate_mode = false;
L50:
	byte_counter = old_byte_counter;

//	TRY

	if (TDC8PCI2.use_normal_method_2nd_card) {
		if (DAQVersion >= 20020408) {
			*input_lmf >> TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
			*input_lmf >> TDC8PCI2.OpenTime_1st_card ;			byte_counter += sizeof(__int32); // open time 1st card
			*input_lmf >> TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
			*input_lmf >> TDC8PCI2.TriggerFalling_1st_card;		byte_counter += sizeof(__int32); // trigger falling edge 1st card
			*input_lmf >> TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card
			*input_lmf >> TDC8PCI2.EmptyCounter_1st_card;		byte_counter += sizeof(__int32); // EmptyCounter 1st card
			*input_lmf >> TDC8PCI2.EmptyCounter_since_last_Event_1st_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 1st card
		}
		*input_lmf >> TDC8PCI2.sync_test_on_off;			byte_counter += sizeof(__int32); // sync test on/off
		*input_lmf >> TDC8PCI2.io_address_2nd_card;			byte_counter += sizeof(__int32); // io address 2nd card
		*input_lmf >> TDC8PCI2.GateDelay_2nd_card;			byte_counter += sizeof(__int32); // gate delay 2nd card
		*input_lmf >> TDC8PCI2.OpenTime_2nd_card;			byte_counter += sizeof(__int32); // open time 2nd card
		*input_lmf >> TDC8PCI2.WriteEmptyEvents_2nd_card;	byte_counter += sizeof(__int32); // write empty events 2nd card
		*input_lmf >> TDC8PCI2.TriggerFallingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger falling edge 2nd card
		*input_lmf >> TDC8PCI2.TriggerRisingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger rising edge 2nd card
		*input_lmf >> TDC8PCI2.EmptyCounter_2nd_card;		byte_counter += sizeof(__int32); // EmptyCounter 2nd card
		*input_lmf >> TDC8PCI2.EmptyCounter_since_last_Event_2nd_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 2nd card
	} else {
		*input_lmf >> module_2nd;							byte_counter += sizeof(__int32);	// indicator for 2nd module data
		*input_lmf >> TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
		*input_lmf >> TDC8PCI2.OpenTime_1st_card ;			byte_counter += sizeof(__int32); // open time 1st card
		*input_lmf >> TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
		*input_lmf >> TDC8PCI2.TriggerFalling_1st_card;		byte_counter += sizeof(__int32); // trigger falling edge 1st card
		*input_lmf >> TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card

		*input_lmf >> TDC8PCI2.GateDelay_2nd_card;			byte_counter += sizeof(__int32); // gate delay 2nd card
		if (!desperate_mode) { // this is only a quick fix.
			*input_lmf >> TDC8PCI2.OpenTime_2nd_card;			byte_counter += sizeof(__int32); // open time 2nd card
			*input_lmf >> TDC8PCI2.WriteEmptyEvents_2nd_card;	byte_counter += sizeof(__int32); // write empty events 2nd card
			*input_lmf >> TDC8PCI2.TriggerFallingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger falling edge 2nd card
			*input_lmf >> TDC8PCI2.TriggerRisingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger rising edge 2nd card
		}
	}

/*
	CATCH( CArchiveException, e )
		if (!TDC8PCI2.use_normal_method_2nd_card) return 0;
		TDC8PCI2.use_normal_method_2nd_card = false;
		in_ar->Close(); delete in_ar; in_ar = 0;
		input_lmf->seek(StartPosition);
		in_ar = new CArchive(input_lmf,CArchive::load);
		goto L50;
	END_CATCH
*/

L200:

	if (DAQVersion < 20080507) {
		if (byte_counter != __int32(User_header_size - 12)) {
			if (desperate_mode) return 0;
			if (!TDC8PCI2.use_normal_method_2nd_card) {
				desperate_mode = true;
			}
			TDC8PCI2.use_normal_method_2nd_card = false;
			input_lmf->seek(StartPosition);
			goto L50;
		}
	}

	if (DAQVersion >= 20080507 && DAQVersion < 20110208) {
		if (byte_counter != __int32(User_header_size - 20)) {
			byte_counter = -1000; // XXX (this line is okay. I have put the XXX just to bring it to attention)
		}
	}

	return byte_counter;
}






/////////////////////////////////////////////////////////////////
bool LMF_IO::ReadNextfADC4packet(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words)
/////////////////////////////////////////////////////////////////
{
	bEnd_of_group_detected = false;
	if (fADC4.packet_count < 1) {
		if (this->LMF_Version > 11) {
		__int32 event_end_marker;
		*input_lmf >> event_end_marker;
			if (event_end_marker != fADC4.GroupEndMarker) {
			errorflag = 2;
			return false;
		}
		*input_lmf >> ui64_timestamp;
		} else {
			ui64_timestamp = 0;
		}
		DOUBLE_timestamp = ui64_timestamp * 1.e-12; // convert to seconds
		*input_lmf >> fADC4.packet_count;
	}
	
	if (errorflag) return false;

	while (fADC4.packet_count > 0) {
		fADC4.packet_count--;

		*input_lmf >> packet->channel;		if (errorflag) return false;
		*input_lmf >> packet->card;			if (errorflag) return false;
		*input_lmf >> packet->type;			if (errorflag) return false;
		*input_lmf >> packet->flags;		if (errorflag) return false;
		*input_lmf >> packet->length;		if (errorflag) return false;
		*input_lmf >> packet->timestamp;	if (errorflag) return false;
		
		//if (packet->flags) continue;
		
		if (packet->card >10) continue;

		if (packet->type == CRONO_PACKET_TYPE_16_BIT_SIGNED) {
			__int32 length = __int32(packet->length*4) < buffersize_in_16bit_words ? packet->length*4 : buffersize_in_16bit_words;
			input_lmf->read((char*)i16buffer,length * sizeof(__int16));
			if (errorflag) return false;
			for (unsigned __int32 k=0;k<(packet->length*4 - length);k++) {
				__int16 dummy16;
				*input_lmf >> dummy16;
				if (errorflag) return false;
			}
			packet->length = length >> 2;
		}

		if (   packet->type != CRONO_PACKET_TYPE_16_BIT_SIGNED
		    && packet->type != CRONO_PACKET_TYPE_TDC_RISING
		    && packet->type != CRONO_PACKET_TYPE_TDC_FALLING
		    && packet->type != CRONO_PACKET_TYPE_TDC_DATA
		    && packet->type != CRONO_PACKET_TYPE_TIMESTAMP_ONLY) {
					continue;
		}

		if (fADC4.packet_count > 0) return true;
	}


	if (this->LMF_Version > 12) {
		unsigned __int32 changed_mask;
		*input_lmf >> changed_mask;
		if (changed_mask)
		{
			for (__int32 i = 0; i<32; i++)
			{
				if (changed_mask & 0x1)
				{
					double double_temp;
					*input_lmf >> double_temp;
					Parameter[901 + i] = double_temp;
				}
				changed_mask >>= 1;
			}
		}
	}



	bEnd_of_group_detected = true;
	number_of_bytes_in_PostEventData = 0;

	if (fADC4.bReadCustomData) {
		*input_lmf >> number_of_bytes_in_PostEventData;
		if (number_of_bytes_in_PostEventData == 12345678) {
			unsigned __int64 pos = input_lmf->tell();
			input_lmf->seek(pos-4);
			number_of_bytes_in_PostEventData = 0;
		}
		if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
			input_lmf->error = 21;
			return false;
		}
		if (number_of_bytes_in_PostEventData < 0) {
			input_lmf->error = 21;
			return false;
		}
		for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
			unsigned __int8 byte_dummy;
			*input_lmf >> byte_dummy;
			ui8_PostEventData[i] = byte_dummy;
		}
	}

	++uint64_number_of_read_events;

	return true;
}






/////////////////////////////////////////////////////////////////
bool LMF_IO::ReadNextfADC8Signal(fADC8_signal_info_struct& signal_info, bool& bEnd_of_group_detected, unsigned __int32 * ui32buffer, __int32 buffersize_in_32bit_words)
/////////////////////////////////////////////////////////////////
{
	if (!input_lmf) {
		errorflag = 9;
		return false;
	}
	if (errorflag) return false;

//L100:

	input_lmf->read((char*)(&signal_info),sizeof(signal_info));

	if (input_lmf->error) {
		if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1;
		return false;
	}

	if (signal_info.adc_channel == -1 && signal_info.ModuleIndex == -1) {
		goto L200;
	}

	if (signal_info.signal_type != 8) {
		__int32 leng = buffersize_in_32bit_words > signal_info.signallength_including_header_in_32bit_words ? signal_info.signallength_including_header_in_32bit_words : buffersize_in_32bit_words;
		input_lmf->read(ui32buffer, leng*4);
		if (input_lmf->error) {
			if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1;
			return false;
		}
		if (leng < signal_info.signallength_including_header_in_32bit_words) {
			__int8 byte_dummy;
			for (__int32 i=0;i<(signal_info.signallength_including_header_in_32bit_words-leng)*4;i++) input_lmf->read(&byte_dummy,1);
		}
	}


L200:

	__int32 header_word;
	*input_lmf >> header_word;

	if (header_word != 12345678 && header_word != 0) {
		errorflag = 2;
		return false;
	}

	if (header_word == fADC8.GroupEndMarker) { // end of group
		bEnd_of_group_detected = true;
		++uint64_number_of_read_events;
	} else bEnd_of_group_detected = false;

	if (bEnd_of_group_detected && fADC8.bReadCustomData) {
		*input_lmf >> number_of_bytes_in_PostEventData;
		if (input_lmf->error) {
			if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1;
			return false;
		}
		if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
			input_lmf->error = 21;
			return false;
		}
		for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
			unsigned __int8 byte_dummy;
			*input_lmf >> byte_dummy;
			ui8_PostEventData[i] = byte_dummy;
		}
	}

	if (input_lmf->error) {
		if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1;
		return false;
	}

	return true;
}





/////////////////////////////////////////////////////////////////
bool LMF_IO::WriteNextfADC8Signal(fADC8_signal_info_struct* signal_info, unsigned __int32 * ui32buffer)
/////////////////////////////////////////////////////////////////
{
	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return false;
	}
	if (fADC8.at_least_1_signal_was_written) *output_lmf << __int32(0);

	output_lmf->write((char*)(signal_info),sizeof(fADC8_signal_info_struct));

	if (signal_info->signal_type != 8) {
		output_lmf->write((char*)ui32buffer, signal_info->signallength_including_header_in_32bit_words*4);
	}

	fADC8.at_least_1_signal_was_written = true;
	
	return true;
}




/////////////////////////////////////////////////////////////////
bool LMF_IO::WritefADC8EndGroupMarker()
/////////////////////////////////////////////////////////////////
{
	if (fADC8.at_least_1_signal_was_written) *output_lmf << __int32(12345678);

	if (fADC8.bReadCustomData) {
		*output_lmf << number_of_bytes_in_PostEventData;
		for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
			*output_lmf << ui8_PostEventData[i];
		}
	}

	fADC8.at_least_1_signal_was_written = false;
	++uint64_number_of_written_events;

	return true;
}







__int32	LMF_IO::ReadfADC8_header_LMFversion9()
{
	__int32 byte_counter = 0;
	*input_lmf >> this->number_of_DAQ_source_strings;  byte_counter += sizeof(__int32);

	DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
	memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
		__int32 string_len;
		*input_lmf >> string_len;	byte_counter += sizeof(__int32);
		DAQ_source_strings[i] = new std::string();
		DAQ_source_strings[i]->reserve(string_len);
		while (string_len > 0) {
			__int8 c;
			*input_lmf >> c;		byte_counter += sizeof(__int8);
			*DAQ_source_strings[i] += c;
			--string_len;
		}
	}

	*input_lmf >> time_reference;	byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution;	byte_counter += sizeof(double);					// tdc resolution in ns
	*input_lmf >> TDCDataType;		byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64;
	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of channels
	number_of_channels =__int32(temp_uint64);

	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of hits
	max_number_of_hits =__int32(temp_uint64);

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)



	// bools
	__int32 counter = 0;
	*input_lmf >> fADC8.number_of_bools;		byte_counter += sizeof(__int32);
	for (__int32 i=counter;i<fADC8.number_of_bools;++i) {
		bool bool_dummy;
		*input_lmf >> bool_dummy;				byte_counter += sizeof(bool);
	}


	// unsigned 32bit
	counter = 0;
	*input_lmf >> fADC8.number_of_uint32s;		byte_counter += sizeof(__int32);
	
	*input_lmf >> fADC8.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.i32NumberOfADCmodules;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iEnableGroupMode;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iTriggerChannel;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iPreSamplings_in_4800ps_units;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iPostSamplings_in_9600ps_units;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iEnableTDCinputs;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[0][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[0][1];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[1][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[1][1];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[2][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[2][1];		byte_counter += sizeof(__int32); counter++;
	for (__int32 i=counter;i<fADC8.number_of_uint32s;++i) {
		unsigned __int32 uint_dummy;
		*input_lmf >> uint_dummy;				byte_counter += sizeof(unsigned __int32);
	}


	// signed 32bit
	counter = 0;
	*input_lmf >> fADC8.number_of_int32s;		byte_counter += sizeof(__int32);

	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*input_lmf >> fADC8.iThreshold_GT[imod][ich];  byte_counter += sizeof(__int32);  counter++;
		}
	}
	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*input_lmf >> fADC8.iThreshold_LT[imod][ich];  byte_counter += sizeof(__int32);  counter++;
		}
	}
	*input_lmf >> fADC8.iSynchronMode[0][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[0][1];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[1][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[1][1];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[2][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[2][1];	byte_counter += sizeof(__int32); counter++;
	
	
	for (__int32 i=counter;i<fADC8.number_of_int32s;++i) {
		__int32 int_dummy;
		*input_lmf >> int_dummy;				byte_counter += sizeof(__int32);
	}

	// doubles
	counter = 0;
	*input_lmf >> fADC8.number_of_doubles;		byte_counter += sizeof(__int32);
	*input_lmf >> fADC8.dGroupRangeStart;		byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dGroupRangeEnd;			byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[0][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[0][1];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[1][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[1][1];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[2][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[2][1];	byte_counter += sizeof(double); counter++;

	for (__int32 i=counter;i<fADC8.number_of_doubles;++i) {
		double double_dummy;
		*input_lmf >> double_dummy;				byte_counter += sizeof(double);
	}

	*input_lmf >> fADC8.GroupEndMarker;			byte_counter += sizeof(__int32); counter++;

	return byte_counter;
}







__int32	LMF_IO::ReadfADC8_header_LMFversion10()
{
	__int32 byte_counter = 0;
	*input_lmf >> this->number_of_DAQ_source_strings;  byte_counter += sizeof(__int32);

	DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
	memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
		__int32 string_len;
		*input_lmf >> string_len;	byte_counter += sizeof(__int32);
		DAQ_source_strings[i] = new std::string();
		DAQ_source_strings[i]->reserve(string_len);
		while (string_len > 0) {
			__int8 c;
			*input_lmf >> c;	byte_counter += sizeof(__int8);
			*DAQ_source_strings[i] += c;
			--string_len;
		}
	}

	*input_lmf >> time_reference;	byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution;	byte_counter += sizeof(double);		// tdc resolution in ns
	*input_lmf >> TDCDataType;		byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64;
	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of channels
	number_of_channels =__int32(temp_uint64);

	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of hits
	max_number_of_hits =__int32(temp_uint64);

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)


	// bools
	__int32 counter = 0;
	*input_lmf >> fADC8.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=counter;i<fADC8.number_of_bools;++i) {
		bool bool_dummy;
		*input_lmf >> bool_dummy;	byte_counter += sizeof(bool);
	}

	// unsigned 32bit
	counter = 0;
	*input_lmf >> fADC8.number_of_uint32s;		byte_counter += sizeof(__int32);
	
	*input_lmf >> fADC8.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.i32NumberOfADCmodules;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iEnableGroupMode;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iTriggerChannel;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iPreSamplings_in_4800ps_units;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iPostSamplings_in_9600ps_units;	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iEnableTDCinputs;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[0][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[0][1];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[1][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[1][1];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[2][0];		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iChannelMode[2][1];		byte_counter += sizeof(__int32); counter++;
	
	for (__int32 i=0;i<3;i++) {
		*input_lmf >> fADC8.firmware_version[i];	; byte_counter += sizeof(__int32); counter++; // adc module firmware
		*input_lmf >> fADC8.serial_number[i];		; byte_counter += sizeof(__int32); counter++; // adc module serial number
	}
	
	*input_lmf >> fADC8.driver_version; byte_counter += sizeof(__int32); counter++; // adc driver version number
	
	for (__int32 i=counter;i<fADC8.number_of_uint32s;++i) {
		unsigned __int32 uint_dummy;
		*input_lmf >> uint_dummy;				byte_counter += sizeof(unsigned __int32);
	}

	// signed 32bit
	counter = 0;
	*input_lmf >> fADC8.number_of_int32s;		byte_counter += sizeof(__int32);

	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*input_lmf >> fADC8.iThreshold_GT[imod][ich];	byte_counter += sizeof(__int32);  counter++;
		}
	}
	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*input_lmf >> fADC8.iThreshold_LT[imod][ich];	byte_counter += sizeof(__int32);  counter++;
		}
	}
	*input_lmf >> fADC8.iSynchronMode[0][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[0][1];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[1][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[1][1];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[2][0];	byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.iSynchronMode[2][1];	byte_counter += sizeof(__int32); counter++;
	
	*input_lmf >> fADC8.GroupEndMarker;			byte_counter += sizeof(__int32); counter++;
	
	if (fADC8.GroupEndMarker != 12345678) return 0;
	
	*input_lmf >> fADC8.veto_gate_length;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.veto_delay_length;		byte_counter += sizeof(__int32); counter++;
	*input_lmf >> fADC8.veto_mask;				byte_counter += sizeof(__int32); counter++;

	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*input_lmf >> fADC8.GND_level[imod][ich]; byte_counter += sizeof(__int32); counter++;
		}
	}

	__int32 int_dummy;
	fADC8.bReadCustomData = false;
	if (counter < fADC8.number_of_int32s) {
		*input_lmf >> int_dummy;			byte_counter += sizeof(__int32);
		fADC8.bReadCustomData = int_dummy ? true:false;
		counter++;
	}
	iLMFcompression = 0;
	if (counter < fADC8.number_of_int32s) {
		*input_lmf >> iLMFcompression;			byte_counter += sizeof(__int32);
		counter++;
	}
	if (counter < fADC8.number_of_int32s) {
		for (__int32 i=counter;i<fADC8.number_of_int32s;++i) {
			*input_lmf >> int_dummy;			byte_counter += sizeof(__int32);
			counter++;
		}
	}

	// doubles
	counter = 0;
	*input_lmf >> fADC8.number_of_doubles;		byte_counter += sizeof(__int32);
	*input_lmf >> fADC8.dGroupRangeStart;		byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dGroupRangeEnd;			byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[0][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[0][1];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[1][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[1][1];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[2][0];	byte_counter += sizeof(double); counter++;
	*input_lmf >> fADC8.dSyncTimeOffset[2][1];	byte_counter += sizeof(double); counter++;
	for (__int32 i=counter;i<fADC8.number_of_doubles;++i) {
		double double_dummy;
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
	}

	*input_lmf >> int_dummy;	byte_counter += sizeof(__int32);
	if (int_dummy != 123456789) return 0;
	
	return byte_counter;
}















__int32	LMF_IO::WritefADC8_header_LMFversion10()
{
	unsigned __int32 unsigned_int_Dummy;

	__int32 byte_counter = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);			// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value

	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);

	if ((number_of_DAQ_source_strings_output < 0) || (!DAQ_source_strings_output)) number_of_DAQ_source_strings_output = 0;
	*output_lmf << number_of_DAQ_source_strings_output;  byte_counter += sizeof(__int32);

	for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {
		unsigned_int_Dummy = (unsigned __int32)(DAQ_source_strings_output[i]->length());
		*output_lmf << unsigned_int_Dummy;   byte_counter += sizeof(__int32);
		output_lmf->write(DAQ_source_strings_output[i]->c_str(),__int32(DAQ_source_strings_output[i]->length()));
		byte_counter += (unsigned __int32)(DAQ_source_strings_output[i]->length());
	}

	*output_lmf <<  time_reference_output;	byte_counter += sizeof(__int32);
	*output_lmf <<  tdcresolution_output;	byte_counter += sizeof(double);		// tdc resolution in ns
	*output_lmf <<  TDCDataType;		byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64 = number_of_channels_output;
	*output_lmf <<  temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of channels

	temp_uint64 = max_number_of_hits;
	*output_lmf <<  temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of hits

	*output_lmf <<  data_format_in_userheader_output;	byte_counter += sizeof(__int32);	// data format (2=short integer)

	// bools
	__int32 counter = 0;
	*output_lmf <<  fADC8.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=counter;i<fADC8.number_of_bools;++i) {
		bool bool_dummy = false;
		*output_lmf <<  bool_dummy;	byte_counter += sizeof(bool);
	}

	// unsigned 32bit
	counter = 0;
	*output_lmf <<  fADC8.number_of_uint32s;		byte_counter += sizeof(__int32);
	
	*output_lmf <<  fADC8.i32NumberOfDAQLoops;		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.i32NumberOfADCmodules;	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iEnableGroupMode;			byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iTriggerChannel;			byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iPreSamplings_in_4800ps_units;	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iPostSamplings_in_9600ps_units;	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iEnableTDCinputs;			byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[0][0];		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[0][1];		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[1][0];		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[1][1];		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[2][0];		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iChannelMode[2][1];		byte_counter += sizeof(__int32); counter++;
	
	for (__int32 i=0;i<3;i++) {
		*output_lmf <<  fADC8.firmware_version[i];	byte_counter += sizeof(__int32); counter++; // adc module firmware
		*output_lmf <<  fADC8.serial_number[i];		byte_counter += sizeof(__int32); counter++; // adc module serial number
	}
	
	*output_lmf <<  fADC8.driver_version; byte_counter += sizeof(__int32); counter++; // adc driver version number
	
	for (__int32 i=counter;i<fADC8.number_of_uint32s;++i) {
		unsigned __int32 uint_dummy = 0;
		*output_lmf <<  uint_dummy;				byte_counter += sizeof(unsigned __int32);
	}

	// signed 32bit
	counter = 0;
	*output_lmf <<  fADC8.number_of_int32s;		byte_counter += sizeof(__int32);

	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*output_lmf <<  fADC8.iThreshold_GT[imod][ich];	byte_counter += sizeof(__int32);  counter++;
		}
	}
	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*output_lmf <<  fADC8.iThreshold_LT[imod][ich];	byte_counter += sizeof(__int32);  counter++;
		}
	}
	*output_lmf <<  fADC8.iSynchronMode[0][0];	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iSynchronMode[0][1];	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iSynchronMode[1][0];	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iSynchronMode[1][1];	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iSynchronMode[2][0];	byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.iSynchronMode[2][1];	byte_counter += sizeof(__int32); counter++;
	
	*output_lmf <<  fADC8.GroupEndMarker;			byte_counter += sizeof(__int32); counter++;
	
	*output_lmf <<  fADC8.veto_gate_length;		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.veto_delay_length;		byte_counter += sizeof(__int32); counter++;
	*output_lmf <<  fADC8.veto_mask;				byte_counter += sizeof(__int32); counter++;

	for (__int32 imod = 0; imod < 3; imod++) {
		for (__int32 ich = 0; ich < 8; ich++) {
			*output_lmf <<  fADC8.GND_level[imod][ich]; byte_counter += sizeof(__int32); counter++;
		}
	}

	if (counter < fADC8.number_of_int32s) {
		for (__int32 i=counter;i<fADC8.number_of_int32s;++i) {
			*output_lmf <<  __int32(fADC8.bReadCustomData ? 1 : 0);	byte_counter += sizeof(__int32);
		}
	}

	// doubles
	counter = 0;
	*output_lmf <<  fADC8.number_of_doubles;		byte_counter += sizeof(__int32);
	*output_lmf <<  fADC8.dGroupRangeStart;			byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dGroupRangeEnd;			byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[0][0];	byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[0][1];	byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[1][0];	byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[1][1];	byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[2][0];	byte_counter += sizeof(double); counter++;
	*output_lmf <<  fADC8.dSyncTimeOffset[2][1];	byte_counter += sizeof(double); counter++;
	for (__int32 i=counter;i<fADC8.number_of_doubles;++i) {
		double double_dummy = 0.;
		*output_lmf <<  double_dummy; byte_counter += sizeof(double);
	}

	__int32 int_dummy = 123456789;
	*output_lmf <<  int_dummy;	byte_counter += sizeof(__int32);

	return byte_counter;
}














/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadfADC8Header()
/////////////////////////////////////////////////////////////////
{
	__int32	byte_counter = 0;

	*input_lmf >> frequency;	byte_counter += sizeof(double);			// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);	// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	byte_counter += DAQ_info_Length;  // DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;

	*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	
	if (LMF_Version == 10) byte_counter += ReadfADC8_header_LMFversion10();
	if (LMF_Version == 9)  byte_counter += ReadfADC8_header_LMFversion9();

	return byte_counter;
}










__int32	LMF_IO::ReadfADC4Header_up_to_v11()
{
	__int32 i32Dummy = 0;
	unsigned __int64 ui64Dummy = 0;
	double dDummy = 0.;
	bool bDummy = false;
	char cDummy = 0;

	__int32 byte_counter = 0;
	
	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;
	
	*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	
	*input_lmf >> this->number_of_DAQ_source_strings;  byte_counter += sizeof(__int32);

	DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
	memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
		__int32 string_len;
		*input_lmf >> string_len;	byte_counter += sizeof(__int32);
		DAQ_source_strings[i] = new std::string();
		DAQ_source_strings[i]->reserve(string_len);
		while (string_len > 0) {
			__int8 c;
			*input_lmf >> c;	byte_counter += sizeof(__int8);
			*DAQ_source_strings[i] += c;
			--string_len;
		}
	}

	*input_lmf >> time_reference;	byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution;	byte_counter += sizeof(double);					// tdc resolution in ns
	*input_lmf >> TDCDataType;		byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64;
	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of channels
	number_of_channels =__int32(temp_uint64);

	*input_lmf >> temp_uint64;		byte_counter += sizeof(unsigned __int64);		// number of hits
	max_number_of_hits =__int32(temp_uint64);

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)

	*input_lmf >> fADC4.i32NumberOfADCmodules;	byte_counter += sizeof(__int32);


	__int32 check_value;
	*input_lmf >> check_value; byte_counter += sizeof(__int32);

	if (check_value == 12345678) {
		for (__int32 i=0;i<fADC4.i32NumberOfADCmodules;i++) {
			*input_lmf >> fADC4.ndigo_parameters[i].bandwidth;		byte_counter += sizeof(double);
			*input_lmf >> fADC4.ndigo_parameters[i].board_id;		byte_counter += sizeof(__int32);
			*input_lmf >> fADC4.ndigo_parameters[i].channel_mask;	byte_counter += sizeof(__int32);
			*input_lmf >> fADC4.ndigo_parameters[i].channels;		byte_counter += sizeof(__int32);
			*input_lmf >> fADC4.ndigo_parameters[i].sample_period;	byte_counter += sizeof(double);
			*input_lmf >> fADC4.ndigo_parameters[i].sample_rate;	byte_counter += sizeof(double);
			*input_lmf >> fADC4.ndigo_parameters[i].version;		byte_counter += sizeof(__int32);

			if (DAQVersion >= 20140617) { // DAQ_VERSION13FADC4
				ndigo_static_info * ndsi = &fADC4.ndigo_info[i];	
				*input_lmf >> ndsi->size;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->version;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->board_id;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->driver_revision;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->firmware_revision;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->board_revision;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->board_configuration;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->adc_resolution;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->nominal_sample_rate;				byte_counter += sizeof(double);
				*input_lmf >> ndsi->analog_bandwidth;			byte_counter += sizeof(double);
				*input_lmf >> ndsi->chip_id;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->board_serial;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->flash_serial_low;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->flash_serial_high;			byte_counter += sizeof(__int32);
				*input_lmf >> ndsi->flash_valid;			byte_counter += sizeof(__int32);
				*input_lmf >> i32Dummy;			byte_counter += sizeof(__int32);
				ndsi->dc_coupled = (i32Dummy != 0) ? true : false;
				*input_lmf >> ndsi->subversion_revision;			byte_counter += sizeof(__int32);
				for (__int32 k=0;k<20;k++) {
					*input_lmf >> ndsi->calibration_date[k];			byte_counter += sizeof(__int8);
				}
			} else {
				memset(&fADC4.ndigo_info[i], 0 ,sizeof(ndigo_static_info));
			}

		}
	} else {
			unsigned __int64 pos = input_lmf->tell(); input_lmf->seek(pos-4); byte_counter -= sizeof(__int32);
	for (__int32 i=0;i<fADC4.i32NumberOfADCmodules;i++) {
		input_lmf->read((char*)&fADC4.ndigo_parameters[i],sizeof(ndigo_param_info)); byte_counter += sizeof(ndigo_param_info);
				__int32 size_diff = sizeof(ndigo_param_info) - fADC4.ndigo_parameters[i].size;
				if (size_diff) {
					pos = input_lmf->tell(); input_lmf->seek(pos-size_diff); byte_counter -= size_diff;
				}
	}

	}

	



	for (unsigned __int32 i = 0; i < number_of_channels; i++) {
		*input_lmf >> fADC4.threshold_GT[i];	byte_counter += sizeof(double);
		*input_lmf >> fADC4.threshold_LT[i];	byte_counter += sizeof(double);
	}

	for (unsigned __int32 i = 0; i < number_of_channels; i++) {
		*input_lmf >> fADC4.set_DC_offset[i];	byte_counter += sizeof(double);
	}

	if (LMF_Version > 10) {
		for (unsigned __int32 i = 0; i < number_of_channels; i++) {
			*input_lmf >> fADC4.GND_level[i];	byte_counter += sizeof(double);
		}
	} else memset(fADC4.GND_level,0,80*sizeof(double));

	__int32 num_temp = fADC4.i32NumberOfADCmodules < 4 ? 4 : fADC4.i32NumberOfADCmodules;
	for (__int32 i=0;i<num_temp;i++) { // write ADC mode 4x1.25Gs or 2x2.5Gs or 1x5Gs)
		*input_lmf >> fADC4.sampling_mode[i];	byte_counter += sizeof(__int32);
	}

	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32); // precursor in steps of 3.2 ns
	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32); // length after end of signal in steps of 3.2 ns
	for (unsigned __int32 i=2;i<number_of_channels;i++) {
		*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32); // read precursor and length for all other channels. Currently only for future use
		*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32);
	}

	*input_lmf >> ui64Dummy;	byte_counter += sizeof(__int64);				// ui64RisingEnable
	*input_lmf >> ui64Dummy;	byte_counter += sizeof(__int64);				// ui64FallingEnable
	*input_lmf >> i32Dummy;		byte_counter += sizeof(__int32);				// i32TriggerEdge
	*input_lmf >> fADC4.iTriggerChannel;	byte_counter += sizeof(__int32);	// i32TriggerChannel
	fADC4.iTriggerChannel--;

	*input_lmf >> bDummy;	byte_counter += sizeof(bool);					// bOutputLevel
	*input_lmf >> fADC4.dGroupRangeStart;	byte_counter += sizeof(double);	// dGroupRangeStart
	*input_lmf >> fADC4.dGroupRangeEnd;		byte_counter += sizeof(double);	// dGroupRangeEnd

	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32);				// read Config file
	fADC4.csConfigFile = "";
	for(__int32 i32Count=0;i32Count<i32Dummy;i32Count++)
	{
		input_lmf->read(&cDummy,1);	byte_counter += sizeof(char);
		fADC4.csConfigFile += cDummy;
	}
	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32);				// read INL file
	fADC4.csINLFile = "";
	for(__int32 i32Count=0;i32Count<i32Dummy;i32Count++)
	{
		input_lmf->read(&cDummy,1);	byte_counter += sizeof(char);
		fADC4.csINLFile += cDummy;
	}
	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32);				// read DNL file
	fADC4.csDNLFile = "";
	for(__int32 i32Count=0;i32Count<i32Dummy;i32Count++)
	{
		input_lmf->read(&cDummy,1);	byte_counter += sizeof(char);
		fADC4.csDNLFile += cDummy;
	}

	*input_lmf >> i32Dummy;	byte_counter += sizeof(__int32);	// dummy
	*input_lmf >> dDummy;	byte_counter += sizeof(double);		// dGroupTimeout

	// read TDCInfo
	*input_lmf >> TDC8HP.Number_of_TDCs;  byte_counter += sizeof(__int32);
	for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
		*input_lmf >> TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelCount;		byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelStart;		byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
		*input_lmf >> TDC8HP.TDC_info[iCount]->serialNumber;		byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->version;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->fifoSize;			byte_counter += sizeof(__int32);
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
		*input_lmf >> TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
	}

	// now handle additional bools
	*input_lmf >> fADC4.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=0;i<fADC4.number_of_bools;++i) {
		*input_lmf >> bDummy;	byte_counter += sizeof(bool);
	}


	// now handle additional 32Bit values (signed or unsigned)
	*input_lmf >> fADC4.number_of_int32s;		byte_counter += sizeof(__int32);
	if (fADC4.number_of_int32s == 3 && this->LMF_Version == 11) this->LMF_Version = 12; // dirty work around for old LMFs
	
	__int32 counter = fADC4.number_of_int32s;
	if (counter > 0) {*input_lmf >> fADC4.i32NumberOfDAQLoops;		counter--;	byte_counter += sizeof(__int32);}
	if (counter > 0) {*input_lmf >> TDC8HP.TDC8HP_DriverVersion;	counter--;	byte_counter += sizeof(__int32);}
	while (counter > 0) {
		*input_lmf >> i32Dummy;
		byte_counter += sizeof(__int32);
		counter--;
	}

	*input_lmf >> fADC4.number_of_doubles;		byte_counter += sizeof(__int32);
	counter = fADC4.number_of_doubles;
	if (LMF_Version > 10) {
		for (__int32 i=0;i<fADC4.i32NumberOfADCmodules;i++) {
			*input_lmf >> fADC4.bits_per_mVolt[i];	byte_counter += sizeof(double);
			counter--;
		}
	} else {
		*input_lmf >> dDummy; byte_counter += sizeof(double);	counter--;
		for (__int32 i=0;i<fADC4.i32NumberOfADCmodules;i++) fADC4.bits_per_mVolt[i] = dDummy;
	}
	while (counter >0)	{*input_lmf >> dDummy;	byte_counter += sizeof(double);	counter--;}

	*input_lmf >> i32Dummy;			byte_counter += sizeof(__int32);
	if (i32Dummy != 1234567) return 0;
	
	return byte_counter;
}










/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadTDC8HPHeader_LMFV_1_to_7(__int32 byte_counter_external)
/////////////////////////////////////////////////////////////////
{
	this->TDC8HP.ui32oldRollOver = 0;
	this->TDC8HP.ui64RollOvers = 0;
	this->TDC8HP.ui32AbsoluteTimeStamp = 0;
	this->TDC8HP.ui64TDC8HP_AbsoluteTimeStamp = 0;

	__int32		byte_counter = 0;

	unsigned __int64 StartPosition1 = input_lmf->tell();

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;

//	input_lmf->flush();
	unsigned __int64 StartPosition = input_lmf->tell();
	__int32 old_byte_counter = byte_counter;
L50:
	byte_counter = old_byte_counter;

	if (DAQVersion > 20080000) TDC8HP.UserHeaderVersion = 4;

	if (TDC8HP.UserHeaderVersion >= 1) {
		*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
		if (LMF_Version == 8) TDC8HP.UserHeaderVersion = 5;
		if (LMF_Version >= 9) TDC8HP.UserHeaderVersion = 6;
		if (LMF_Version >= 10) TDC8HP.UserHeaderVersion = 7;
		if (LMF_Version >= 11) TDC8HP.UserHeaderVersion = 8;
	}
	
	if (TDC8HP.UserHeaderVersion >= 5) {
		input_lmf->seek(StartPosition1);
		if (TDC8HP.UserHeaderVersion < 7) return ReadTDC8HPHeader_LMFV_8_to_9();
		return ReadTDC8HPHeader_LMFV_10();
	}

	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution; byte_counter += sizeof(double);		// tdc resolution in ns
	if (tdcresolution < 0.0001 || tdcresolution > 100.) {
		if (TDC8HP.UserHeaderVersion != 0) return 0;
		TDC8HP.UserHeaderVersion = 1;

		input_lmf->seek(StartPosition);
		goto L50;
	}
	//tdcresolution = double(__int32(tdcresolution*1000000+0.01))/1000000.;
	if (TDC8HP.UserHeaderVersion >= 1) {
		*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);
	}


	unsigned __int64 TempPosition = input_lmf->tell();
	bool bRead32bit = false;

L100:

	if (DAQVersion < 20080000 || bRead32bit) {
		*input_lmf >> number_of_channels; byte_counter += sizeof(__int32);			// number of channels
	} else {
		unsigned __int64 temp;
		*input_lmf >> temp;
		if (temp > 4000000000) {
			input_lmf->seek(TempPosition);
			bRead32bit = true;
			goto L100;
		}
		number_of_channels =__int32(temp);
		byte_counter += sizeof(unsigned __int64);			// number of channels
	}

	if (number_of_channels < 1 || number_of_channels > 100) {
		if (TDC8HP.UserHeaderVersion != 0) return 0;
//		in_ar->Close(); delete in_ar; in_ar = 0;
		input_lmf->seek(StartPosition);
//		in_ar = new CArchive(input_lmf,CArchive::load);
		goto L50;
	}


	if (DAQVersion < 20080000 || bRead32bit) {
		*input_lmf >> max_number_of_hits; byte_counter += sizeof(__int32);			// number of hits
	} else {
		unsigned __int64 temp;
		*input_lmf >> temp;
		max_number_of_hits =__int32(temp);
		byte_counter += sizeof(unsigned __int64);			// number of hits
	}

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	if (max_number_of_hits < 1 || max_number_of_hits > 100) {
		if (TDC8HP.UserHeaderVersion != 0) return 0;
//		in_ar->Close(); delete in_ar; in_ar = 0;
		input_lmf->seek(StartPosition);
//		in_ar = new CArchive(input_lmf,CArchive::load);
		goto L50;
	}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)
	*input_lmf >> TDC8HP.no_config_file_read;	byte_counter += sizeof(__int32);	// parameter 60-1
	if (DAQVersion < 20080000 || bRead32bit) {
		unsigned __int32 temp;
		*input_lmf >> temp; TDC8HP.RisingEnable_p61 = temp;		byte_counter += sizeof(__int32);	// parameter 61-1
		*input_lmf >> temp; TDC8HP.FallingEnable_p62 = temp;		byte_counter += sizeof(__int32);	// parameter 62-1
	} else {
		*input_lmf >> TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
		*input_lmf >> TDC8HP.FallingEnable_p62;		byte_counter += sizeof(__int64);	// parameter 62-1
	}
	*input_lmf >> TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*input_lmf >> TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1
	*input_lmf >> TDC8HP.OutputLevel_p65;		byte_counter += sizeof(__int32);	// parameter 65-1
	*input_lmf >> TDC8HP.GroupingEnable_p66;	byte_counter += sizeof(__int32);	// parameter 66-1
	*input_lmf >> TDC8HP.AllowOverlap_p67;		byte_counter += sizeof(__int32);	// parameter 67-1
	*input_lmf >> TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);		// parameter 68-1
	*input_lmf >> TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);		// parameter 69-1
	*input_lmf >> TDC8HP.GroupRangeEnd_p70;		byte_counter += sizeof(double);		// parameter 70-1
	*input_lmf >> TDC8HP.ExternalClock_p71;		byte_counter += sizeof(__int32);	// parameter 71-1
	*input_lmf >> TDC8HP.OutputRollOvers_p72;	byte_counter += sizeof(__int32);	// parameter 72-1
	*input_lmf >> TDC8HP.DelayTap0_p73;			byte_counter += sizeof(__int32);	// parameter 73-1
	*input_lmf >> TDC8HP.DelayTap1_p74;			byte_counter += sizeof(__int32);	// parameter 74-1
	*input_lmf >> TDC8HP.DelayTap2_p75;			byte_counter += sizeof(__int32);	// parameter 75-1
	*input_lmf >> TDC8HP.DelayTap3_p76;			byte_counter += sizeof(__int32);	// parameter 76-1
	*input_lmf >> TDC8HP.INL_p80;				byte_counter += sizeof(__int32);	// parameter 80-1
	*input_lmf >> TDC8HP.DNL_p81;				byte_counter += sizeof(__int32);	// parameter 81-1

	*input_lmf >> TDC8HP.csConfigFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csConfigFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csConfigFile_Length);
	__int8_temp[TDC8HP.csConfigFile_Length] = 0;
	TDC8HP.csConfigFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csConfigFile_Length;

	*input_lmf >> TDC8HP.csINLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csINLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csINLFile_Length);
	__int8_temp[TDC8HP.csINLFile_Length] = 0;
	TDC8HP.csINLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csINLFile_Length;

	*input_lmf >> TDC8HP.csDNLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csDNLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csDNLFile_Length);
	__int8_temp[TDC8HP.csDNLFile_Length] = 0;
	TDC8HP.csDNLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csDNLFile_Length;

	if (DAQVersion < 20080000) {
		if (byte_counter == __int32(User_header_size - 12 - 4)) {  // Cobold 2002 v11
			if (TDC8HP.UserHeaderVersion < 2) TDC8HP.UserHeaderVersion = 2;
			*input_lmf >> TDC8HP.SyncValidationChannel;  byte_counter += sizeof(__int32); 	// parameter 77-1
		}
		if (byte_counter == __int32(User_header_size - 12 - 4 - 4)) { // never used in official Cobold releases
			TDC8HP.UserHeaderVersion = 3;
			*input_lmf >> TDC8HP.SyncValidationChannel;  byte_counter += sizeof(__int32);
			*input_lmf >> TDC8HP.VHR_25ps;				 byte_counter += sizeof(bool);
		}
	} else {
		if (this->data_format_in_userheader == -1) TDC8HP.variable_event_length = 1;
		*input_lmf >> TDC8HP.SyncValidationChannel;  byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.VHR_25ps;  byte_counter += sizeof(bool);

		*input_lmf >> TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);

		*input_lmf >> TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
		*input_lmf >> TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
		*input_lmf >> TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

		unsigned __int64 TempPosition = input_lmf->tell();
		bool bNoTDCInfoRead = false;
		__int32 old_byte_counter2 = byte_counter;

L200:
		byte_counter = old_byte_counter2;
		input_lmf->seek(TempPosition);

		if (!bNoTDCInfoRead) {
			// read TDCInfo
			*input_lmf >> TDC8HP.Number_of_TDCs;  byte_counter += sizeof(__int32);
				if (TDC8HP.Number_of_TDCs<0 || TDC8HP.Number_of_TDCs>3) {bNoTDCInfoRead = true;	goto L200;}
			for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
				*input_lmf >> TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->channelCount;		byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->channelStart;		byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
				*input_lmf >> TDC8HP.TDC_info[iCount]->serialNumber;		byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->version;				byte_counter += sizeof(__int32);
				*input_lmf >> TDC8HP.TDC_info[iCount]->fifoSize;			byte_counter += sizeof(__int32);
				input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);	byte_counter += sizeof(__int32)*8*1024;
				input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
				*input_lmf >> TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
			}
		}

		bool bool_dummy;
		*input_lmf >> bool_dummy; byte_counter += sizeof(bool);
		*input_lmf >> bool_dummy; byte_counter += sizeof(bool);
		*input_lmf >> bool_dummy; byte_counter += sizeof(bool);
		*input_lmf >> bool_dummy; byte_counter += sizeof(bool);
		*input_lmf >> bool_dummy; byte_counter += sizeof(bool);

		*input_lmf >> TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
		if (TDC8HP.i32NumberOfDAQLoops == 0x12345678) TDC8HP.i32NumberOfDAQLoops = 1;
		*input_lmf >> TDC8HP.TDC8HP_DriverVersion;	byte_counter += sizeof(__int32);
		if (TDC8HP.TDC8HP_DriverVersion == 0x12345678) TDC8HP.TDC8HP_DriverVersion = 0;
		*input_lmf >> TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
		if (TDC8HP.iTriggerChannelMask == 0x12345678) TDC8HP.iTriggerChannelMask = 0;
		*input_lmf >> TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32);
		if (TDC8HP.iTime_zero_channel == 0x12345678) TDC8HP.iTime_zero_channel = -1;

		if (TDC8HP.i32NumberOfDAQLoops != 0x12345678) {
			if (TDC8HP.i32NumberOfDAQLoops<-2 || TDC8HP.i32NumberOfDAQLoops > 2e6) {bNoTDCInfoRead = true;	goto L200;}
		}
		
		if (TDC8HP.iTime_zero_channel != 0x12345678) {
			if (TDC8HP.iTime_zero_channel <-1 || TDC8HP.iTime_zero_channel  > 50)  {bNoTDCInfoRead = true;	goto L200;}
		}

		__int32 int_dummy;
		*input_lmf >> int_dummy; byte_counter += sizeof(__int32);

		double double_dummy;
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
		
		if (!bNoTDCInfoRead && (byte_counter_external + byte_counter == __int32(User_header_size+4))) {bNoTDCInfoRead = true;	goto L200;}
	}
	
	if (bRead32bit && LMF_Version == 6) {
		TDC8HP.exotic_file_type = 1;
	}

	return byte_counter;
}











/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadTDC8HPHeader_LMFV_8_to_9()
// reads LMF headers from Cobold 2008 R2 (release August 2009)
/////////////////////////////////////////////////////////////////
{
	this->TDC8HP.ui32oldRollOver = 0;
	this->TDC8HP.ui64RollOvers = 0;
	this->TDC8HP.ui32AbsoluteTimeStamp = 0;
	this->TDC8HP.ui64TDC8HP_AbsoluteTimeStamp = 0;

	__int32		byte_counter = 0;

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;

//	input_lmf->flush();
	//unsigned __int64 StartPosition = input_lmf->tell();
	__int32 old_byte_counter = byte_counter;

	byte_counter = old_byte_counter;

	if (DAQVersion > 20080000) TDC8HP.UserHeaderVersion = 4;

	if (TDC8HP.UserHeaderVersion >= 1) {
		*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
		if (LMF_Version == 8) TDC8HP.UserHeaderVersion = 5;
		if (LMF_Version >= 9) TDC8HP.UserHeaderVersion = 6;
	}

	*input_lmf >> this->number_of_DAQ_source_strings;  byte_counter += sizeof(__int32);

	DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
	memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
		__int32 string_len;
		*input_lmf >> string_len;    byte_counter += sizeof(__int32);
		DAQ_source_strings[i] = new std::string();
		DAQ_source_strings[i]->reserve(string_len);
		while (string_len > 0) {
			__int8 c;
			*input_lmf >> c;     byte_counter += sizeof(__int8);
			*DAQ_source_strings[i] += c;
			--string_len;
		}
	}
	
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution; byte_counter += sizeof(double);		// tdc resolution in ns
	*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64;
	*input_lmf >> temp_uint64;
	number_of_channels =__int32(temp_uint64);
	byte_counter += sizeof(unsigned __int64);			// number of channels

	*input_lmf >> temp_uint64;
	max_number_of_hits =__int32(temp_uint64);
	byte_counter += sizeof(unsigned __int64);			// number of hits

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)

	bool temp_bool;

	*input_lmf >> temp_bool; TDC8HP.no_config_file_read = temp_bool; byte_counter += sizeof(bool);	// parameter 60-1
		
	*input_lmf >> TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
	*input_lmf >> TDC8HP.FallingEnable_p62;		byte_counter += sizeof(__int64);	// parameter 62-1
	*input_lmf >> TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*input_lmf >> TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1

	*input_lmf >> temp_bool; TDC8HP.OutputLevel_p65    = temp_bool;	byte_counter += sizeof(bool);	// parameter 65-1
	*input_lmf >> temp_bool; TDC8HP.GroupingEnable_p66 = temp_bool;	byte_counter += sizeof(bool);	// parameter 66-1
	*input_lmf >> temp_bool; TDC8HP.AllowOverlap_p67   = temp_bool;	byte_counter += sizeof(bool);	// parameter 67-1

	*input_lmf >> TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);		// parameter 68-1
	*input_lmf >> TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);		// parameter 69-1
	*input_lmf >> TDC8HP.GroupRangeEnd_p70;		byte_counter += sizeof(double);		// parameter 70-1

	*input_lmf >> temp_bool; TDC8HP.ExternalClock_p71   = temp_bool;	byte_counter += sizeof(bool);	// parameter 71-1
	*input_lmf >> temp_bool; TDC8HP.OutputRollOvers_p72 = temp_bool;	byte_counter += sizeof(bool);	// parameter 72-1

	*input_lmf >> TDC8HP.DelayTap0_p73;			byte_counter += sizeof(__int32);	// parameter 73-1
	*input_lmf >> TDC8HP.DelayTap1_p74;			byte_counter += sizeof(__int32);	// parameter 74-1
	*input_lmf >> TDC8HP.DelayTap2_p75;			byte_counter += sizeof(__int32);	// parameter 75-1
	*input_lmf >> TDC8HP.DelayTap3_p76;			byte_counter += sizeof(__int32);	// parameter 76-1

	*input_lmf >> temp_bool; TDC8HP.INL_p80 = temp_bool;	byte_counter += sizeof(bool);	// parameter 80-1
	*input_lmf >> temp_bool; TDC8HP.DNL_p81 = temp_bool;	byte_counter += sizeof(bool);	// parameter 81-1

	*input_lmf >> TDC8HP.csConfigFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csConfigFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csConfigFile_Length);
	__int8_temp[TDC8HP.csConfigFile_Length] = 0;
	TDC8HP.csConfigFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csConfigFile_Length;

	*input_lmf >> TDC8HP.csINLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csINLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csINLFile_Length);
	__int8_temp[TDC8HP.csINLFile_Length] = 0;
	TDC8HP.csINLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csINLFile_Length;

	*input_lmf >> TDC8HP.csDNLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csDNLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csDNLFile_Length);
	__int8_temp[TDC8HP.csDNLFile_Length] = 0;
	TDC8HP.csDNLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csDNLFile_Length;

	if (this->data_format_in_userheader == -1) TDC8HP.variable_event_length = 1;
	*input_lmf >> TDC8HP.SyncValidationChannel;  byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.VHR_25ps;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);
	*input_lmf >> TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

	//// read TDCInfo
	*input_lmf >> TDC8HP.Number_of_TDCs;  byte_counter += sizeof(__int32);
	for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
		*input_lmf >> TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelCount;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelStart;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
		*input_lmf >> TDC8HP.TDC_info[iCount]->serialNumber;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->version;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->fifoSize;				byte_counter += sizeof(__int32);
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
		*input_lmf >> TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
	}

	*input_lmf >> TDC8HP.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=0;i<TDC8HP.number_of_bools;++i) {
		bool bool_dummy;
		*input_lmf >> bool_dummy;	byte_counter += sizeof(bool);
	}

	*input_lmf >> TDC8HP.number_of_int32s; byte_counter += sizeof(__int32);

	*input_lmf >> TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.TDC8HP_DriverVersion;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32); // 1 is first channel
	
	if (TDC8HP.UserHeaderVersion >= 6) {
		*input_lmf >> TDC8HP.BinsizeType;	byte_counter += sizeof(__int32); // 1 is first channel
		for (__int32 i=5;i<TDC8HP.number_of_int32s;++i) {
			__int32 int_dummy;
			*input_lmf >> int_dummy;	byte_counter += sizeof(__int32);
		}
	} else {
		for (__int32 i=4;i<TDC8HP.number_of_int32s;++i) {
			__int32 int_dummy;
			*input_lmf >> int_dummy;	byte_counter += sizeof(__int32);
		}
	}
	

	

	*input_lmf >> TDC8HP.number_of_doubles; byte_counter += sizeof(__int32);
	if (TDC8HP.UserHeaderVersion == 5) {
		for (__int32 i=0;i<TDC8HP.number_of_doubles;++i) {
			double double_dummy;
			*input_lmf >> double_dummy; byte_counter += sizeof(double);
		}
	}
	if (TDC8HP.UserHeaderVersion >= 6) {
		*input_lmf >> TDC8HP.OffsetTimeZeroChannel_s;  byte_counter += sizeof(double);
		for (__int32 i=1;i<TDC8HP.number_of_doubles;++i) {
			double double_dummy;
			*input_lmf >> double_dummy; byte_counter += sizeof(double);
		}
	}

	return byte_counter;
}





/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadTDC8HPHeader_LMFV_10()
// reads LMF headers from Cobold 2011 R1+R2
/////////////////////////////////////////////////////////////////
{
	this->TDC8HP.ui32oldRollOver = 0;
	this->TDC8HP.ui64RollOvers = 0;
	this->TDC8HP.ui32AbsoluteTimeStamp = 0;
	this->TDC8HP.ui64TDC8HP_AbsoluteTimeStamp = 0;

	__int32		byte_counter = 0;

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;



//	input_lmf->flush();
	//unsigned __int64 StartPosition = input_lmf->tell();
	__int32 old_byte_counter = byte_counter;

	byte_counter = old_byte_counter;

	*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	TDC8HP.UserHeaderVersion = 7;

	*input_lmf >> this->number_of_DAQ_source_strings;  byte_counter += sizeof(__int32);

	DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
	memset(DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {
		__int32 string_len;
		*input_lmf >> string_len;    byte_counter += sizeof(__int32);
		DAQ_source_strings[i] = new std::string();
		DAQ_source_strings[i]->reserve(string_len);
		while (string_len > 0) {
			__int8 c;
			*input_lmf >> c;     byte_counter += sizeof(__int8);
			*DAQ_source_strings[i] += c;
			--string_len;
		}
	}
	
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution; byte_counter += sizeof(double);		// tdc resolution in ns
	*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);

	unsigned __int64 temp_uint64;
	*input_lmf >> temp_uint64;
	number_of_channels =__int32(temp_uint64);
	byte_counter += sizeof(unsigned __int64);			// number of channels

	*input_lmf >> temp_uint64;
	max_number_of_hits =__int32(temp_uint64);
	byte_counter += sizeof(unsigned __int64);			// number of hits

	if (__int32(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (__int32(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);	// data format (2=short integer)

	bool temp_bool;

	*input_lmf >> temp_bool; TDC8HP.no_config_file_read = temp_bool; byte_counter += sizeof(bool);	// parameter 60-1
		
	*input_lmf >> TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
	*input_lmf >> TDC8HP.FallingEnable_p62;		byte_counter += sizeof(__int64);	// parameter 62-1
	*input_lmf >> TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*input_lmf >> TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1

	*input_lmf >> temp_bool; TDC8HP.OutputLevel_p65    = temp_bool;	byte_counter += sizeof(bool);	// parameter 65-1
	*input_lmf >> temp_bool; TDC8HP.GroupingEnable_p66 = temp_bool;	byte_counter += sizeof(bool);	// parameter 66-1
	*input_lmf >> temp_bool; TDC8HP.AllowOverlap_p67   = temp_bool;	byte_counter += sizeof(bool);	// parameter 67-1

	*input_lmf >> TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);		// parameter 68-1
	*input_lmf >> TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);		// parameter 69-1
	*input_lmf >> TDC8HP.GroupRangeEnd_p70;		byte_counter += sizeof(double);		// parameter 70-1

	*input_lmf >> temp_bool; TDC8HP.ExternalClock_p71   = temp_bool;	byte_counter += sizeof(bool);	// parameter 71-1
	*input_lmf >> temp_bool; TDC8HP.OutputRollOvers_p72 = temp_bool;	byte_counter += sizeof(bool);	// parameter 72-1

	*input_lmf >> TDC8HP.DelayTap0_p73;			byte_counter += sizeof(__int32);	// parameter 73-1
	*input_lmf >> TDC8HP.DelayTap1_p74;			byte_counter += sizeof(__int32);	// parameter 74-1
	*input_lmf >> TDC8HP.DelayTap2_p75;			byte_counter += sizeof(__int32);	// parameter 75-1
	*input_lmf >> TDC8HP.DelayTap3_p76;			byte_counter += sizeof(__int32);	// parameter 76-1

	*input_lmf >> temp_bool; TDC8HP.INL_p80 = temp_bool;	byte_counter += sizeof(bool);	// parameter 80-1
	*input_lmf >> temp_bool; TDC8HP.DNL_p81 = temp_bool;	byte_counter += sizeof(bool);	// parameter 81-1

	*input_lmf >> TDC8HP.csConfigFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csConfigFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csConfigFile_Length);
	__int8_temp[TDC8HP.csConfigFile_Length] = 0;
	TDC8HP.csConfigFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csConfigFile_Length;

	*input_lmf >> TDC8HP.csINLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csINLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csINLFile_Length);
	__int8_temp[TDC8HP.csINLFile_Length] = 0;
	TDC8HP.csINLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csINLFile_Length;

	*input_lmf >> TDC8HP.csDNLFile_Length;
	byte_counter += sizeof(__int32);
	__int8_temp = new __int8[TDC8HP.csDNLFile_Length+1];
	input_lmf->read(__int8_temp,TDC8HP.csDNLFile_Length);
	__int8_temp[TDC8HP.csDNLFile_Length] = 0;
	TDC8HP.csDNLFile = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += TDC8HP.csDNLFile_Length;

	if (this->data_format_in_userheader == -1) TDC8HP.variable_event_length = 1;
	*input_lmf >> TDC8HP.SyncValidationChannel;  byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.VHR_25ps;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);
	*input_lmf >> TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
	*input_lmf >> TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

	//// read TDCInfo
	*input_lmf >> TDC8HP.Number_of_TDCs;  byte_counter += sizeof(__int32);
	for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
		*input_lmf >> TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelCount;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->channelStart;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
		*input_lmf >> TDC8HP.TDC_info[iCount]->serialNumber;			byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->version;				byte_counter += sizeof(__int32);
		*input_lmf >> TDC8HP.TDC_info[iCount]->fifoSize;				byte_counter += sizeof(__int32);
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
		input_lmf->read((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
		*input_lmf >> TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
	}

	*input_lmf >> TDC8HP.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=0;i<TDC8HP.number_of_bools;++i) {
		bool bool_dummy;
		*input_lmf >> bool_dummy;	byte_counter += sizeof(bool);
	}

	*input_lmf >> TDC8HP.number_of_int32s; byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.TDC8HP_DriverVersion;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32); // 1 means: first channel
	*input_lmf >> TDC8HP.BinsizeType;	byte_counter += sizeof(__int32); // 1 means: first channel
	for (__int32 i=5;i<TDC8HP.number_of_int32s;++i) {
		__int32 int_dummy;
		*input_lmf >> int_dummy;	byte_counter += sizeof(__int32);
	}

	*input_lmf >> TDC8HP.number_of_doubles; byte_counter += sizeof(__int32);
	*input_lmf >> TDC8HP.OffsetTimeZeroChannel_s;  byte_counter += sizeof(double);
	for (__int32 i=1;i<TDC8HP.number_of_doubles;++i) {
		double double_dummy;
		*input_lmf >> double_dummy; byte_counter += sizeof(double);
	}

	return byte_counter;
}













/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadTCPIPHeader()
/////////////////////////////////////////////////////////////////
{
	__int32	byte_counter = 0;
 

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;
	*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);

	*input_lmf >> data_format_in_userheader;	byte_counter   += sizeof(__int32);
	*input_lmf >> number_of_channels;			byte_counter   += sizeof(__int32);
	*input_lmf >> max_number_of_hits;			byte_counter   += sizeof(__int32);

	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}

	tdcresolution = 1.;

	return byte_counter;
}














/////////////////////////////////////////////////////////////////
__int32	LMF_IO::ReadHM1Header()
/////////////////////////////////////////////////////////////////
{
	__int32		byte_counter;
	byte_counter = 0;

	HM1.use_normal_method = false;

	*input_lmf >> frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*input_lmf >> IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*input_lmf >> timestamp_format;	byte_counter += sizeof(__int32);	// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	*input_lmf >> DAQ_info_Length;	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	__int8 * __int8_temp = new __int8[DAQ_info_Length+1];
	input_lmf->read(__int8_temp,DAQ_info_Length);	// DAQInfo always 8th value
	__int8_temp[DAQ_info_Length] = 0;
	DAQ_info = __int8_temp;
	delete[] __int8_temp; __int8_temp = 0;
	byte_counter += DAQ_info_Length;


	__int32 nominalHeaderLength;
	nominalHeaderLength = sizeof(__int32)*21 + sizeof(double)*2 + DAQ_info_Length;	// size of user defined header
	if (DAQ_ID == DAQ_ID_HM1_ABM) nominalHeaderLength += 24*sizeof(__int32);

	if (DAQVersion >= 20080507) HM1.use_normal_method = true;

	if (nominalHeaderLength == __int32(User_header_size)) HM1.use_normal_method = true;

	if (DAQVersion >= 20020408 && HM1.use_normal_method) {
		*input_lmf >> LMF_Version; byte_counter += sizeof(__int32);
	}

	*input_lmf >> system_timeout; byte_counter += sizeof(__int32);		//   system time-out
	*input_lmf >> time_reference; byte_counter += sizeof(__int32);
	*input_lmf >> HM1.FAK_DLL_Value; byte_counter += sizeof(__int32);
	*input_lmf >> HM1.Resolution_Flag; byte_counter += sizeof(__int32);
	*input_lmf >> HM1.trigger_mode_for_start; byte_counter += sizeof(__int32);
	*input_lmf >> HM1.trigger_mode_for_stop; byte_counter += sizeof(__int32);
	*input_lmf >> tdcresolution; byte_counter += sizeof(double);		// tdc resolution in ns

	if (DAQVersion >= 20020408 && HM1.use_normal_method) {
		*input_lmf >> TDCDataType; byte_counter += sizeof(__int32);
	}

	*input_lmf >> HM1.Even_open_time; byte_counter += sizeof(__int32);
	*input_lmf >> HM1.Auto_Trigger; byte_counter += sizeof(__int32);

	if (DAQVersion >= 20020408) {
		__int64 i64_dummy;
		*input_lmf >> i64_dummy; number_of_channels = (unsigned __int32)i64_dummy; byte_counter += sizeof(__int64);			// number of channels
		*input_lmf >> i64_dummy; max_number_of_hits = (unsigned __int32)i64_dummy; byte_counter += sizeof(__int64);			// number of hits
	} else {
	*input_lmf >> number_of_channels; byte_counter += sizeof(__int32);			// number of channels
	*input_lmf >> max_number_of_hits; byte_counter += sizeof(__int32);			// number of hits
	}
	if (int(number_of_channels) > num_channels)	{errorflag = 19;	return -100000;}
	if (int(max_number_of_hits) > num_ions)		{errorflag = 20;	return -100000;}
	*input_lmf >> HM1.set_bits_for_GP1; byte_counter += sizeof(__int32);
	*input_lmf >> data_format_in_userheader;	byte_counter += sizeof(__int32);				// data format (2=short integer)
	*input_lmf >> module_2nd;	byte_counter += sizeof(__int32);		// indicator for 2nd module data

	if (DAQ_ID == DAQ_ID_2HM1) {
		*input_lmf >> DAQSubVersion;						byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_FAK_DLL_Value;				byte_counter += sizeof(__int32);  // parameter 10-1
		*input_lmf >> HM1.TWOHM1_Resolution_Flag;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_trigger_mode_for_start;	byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_trigger_mode_for_stop;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_res_adjust;				byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_tdcresolution;				byte_counter += sizeof(double);
		*input_lmf >> HM1.TWOHM1_test_overflow;				byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_number_of_channels;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_number_of_hits;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_set_bits_for_GP1;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_HM1_ID_1;					byte_counter += sizeof(__int32);
		*input_lmf >> HM1.TWOHM1_HM1_ID_2;					byte_counter += sizeof(__int32);
	}

	if (DAQ_ID == DAQ_ID_HM1_ABM) {
		max_number_of_hits = 1;
		*input_lmf >> HM1.ABM_m_xFrom;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_xTo;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_yFrom;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_yTo;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_xMin;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_xMax;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_yMin;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_yMax;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_xOffset;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_yOffset;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_m_zOffset;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_Mode;				byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_OsziDarkInvert;	byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_ErrorHisto;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_XShift;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_YShift;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_ZShift;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_ozShift;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_wdShift;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_ucLevelXY;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_ucLevelZ;			byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_uiABMXShift;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_uiABMYShift;		byte_counter += sizeof(__int32);
		*input_lmf >> HM1.ABM_uiABMZShift;		byte_counter += sizeof(__int32);
		if (DAQVersion > 20080507) {
			__int32 i32Dummy;
			*input_lmf >> i32Dummy;		byte_counter += sizeof(__int32);
		}
	}

	if (DAQ_ID == 1 && DAQVersion == 20080507 && LMF_Version == 7) {
		input_lmf->flush();
		input_lmf->seek(Headersize+User_header_size);
		input_lmf->flush();
	}

	return byte_counter;
}





























/////////////////////////////////////////////////////////////////
__int32 LMF_IO::WriteTCPIPHeader()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value

	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);
	*output_lmf << time_reference; byte_counter += sizeof(__int32);


	*output_lmf << data_format_in_userheader_output;					byte_counter   += sizeof(__int32);
	*output_lmf << number_of_channels_output;			byte_counter   += sizeof(__int32);
	*output_lmf << max_number_of_hits_output;			byte_counter   += sizeof(__int32);

	return byte_counter;
}

















/////////////////////////////////////////////////////////////////
__int32	LMF_IO::Write2TDC8PCI2Header()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;
	//__int32 int_Dummy = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += __int32(DAQ_info.length());

	*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);

	*output_lmf << system_timeout_output; byte_counter += sizeof(__int32);		//   system time-out
	*output_lmf << time_reference_output; byte_counter += sizeof(__int32);
	*output_lmf << common_mode_output; byte_counter += sizeof(__int32);		//   0 common start    1 common stop
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	*output_lmf << TDCDataType; byte_counter += sizeof(__int32);
	*output_lmf << timerange; byte_counter += sizeof(__int32);	// time range of the tdc in microseconds

	*output_lmf << number_of_channels_output; byte_counter += sizeof(__int32);			// number of channels
	*output_lmf << max_number_of_hits_output; byte_counter += sizeof(__int32);			// number of hits
	*output_lmf << number_of_channels2_output; byte_counter += sizeof(__int32);	// number of channels2
	*output_lmf << max_number_of_hits2_output; byte_counter += sizeof(__int32);	// number of hits2
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	if (TDC8PCI2.use_normal_method_2nd_card) {
		if (DAQVersion_output >= 20020408) {
			*output_lmf << TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
			*output_lmf << TDC8PCI2.OpenTime_1st_card ;			byte_counter += sizeof(__int32); // open time 1st card
			*output_lmf << TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
			*output_lmf << TDC8PCI2.TriggerFalling_1st_card;	byte_counter += sizeof(__int32); // trigger falling edge 1st card
			*output_lmf << TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card
			*output_lmf << TDC8PCI2.EmptyCounter_1st_card;		byte_counter += sizeof(__int32); // EmptyCounter 1st card
			*output_lmf << TDC8PCI2.EmptyCounter_since_last_Event_1st_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 1st card
		}
		*output_lmf << TDC8PCI2.sync_test_on_off;			byte_counter += sizeof(__int32); // sync test on/off
		*output_lmf << TDC8PCI2.io_address_2nd_card;		byte_counter += sizeof(__int32); // io address 2nd card
		*output_lmf << TDC8PCI2.GateDelay_2nd_card;			byte_counter += sizeof(__int32); // gate delay 2nd card
		*output_lmf << TDC8PCI2.OpenTime_2nd_card;			byte_counter += sizeof(__int32); // open time 2nd card
		*output_lmf << TDC8PCI2.WriteEmptyEvents_2nd_card;	byte_counter += sizeof(__int32); // write empty events 2nd card
		*output_lmf << TDC8PCI2.TriggerFallingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger falling edge 2nd card
		*output_lmf << TDC8PCI2.TriggerRisingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger rising edge 2nd card
		*output_lmf << TDC8PCI2.EmptyCounter_2nd_card;		byte_counter += sizeof(__int32); // EmptyCounter 2nd card
		*output_lmf << TDC8PCI2.EmptyCounter_since_last_Event_2nd_card;	byte_counter += sizeof(__int32); // Empty Counter since last event 2nd card
	} else {
		*output_lmf << module_2nd;							byte_counter += sizeof(__int32); // indicator for 2nd module data
		*output_lmf << TDC8PCI2.GateDelay_1st_card;			byte_counter += sizeof(__int32); // gate delay 1st card
		*output_lmf << TDC8PCI2.OpenTime_1st_card ;			byte_counter += sizeof(__int32); // open time 1st card
		*output_lmf << TDC8PCI2.WriteEmptyEvents_1st_card;	byte_counter += sizeof(__int32); // write empty events 1st card
		*output_lmf << TDC8PCI2.TriggerFalling_1st_card;	byte_counter += sizeof(__int32); // trigger falling edge 1st card
		*output_lmf << TDC8PCI2.TriggerRising_1st_card;		byte_counter += sizeof(__int32); // trigger rising edge 1st card

		*output_lmf << TDC8PCI2.GateDelay_2nd_card;			byte_counter += sizeof(__int32); // gate delay 2nd card
		*output_lmf << TDC8PCI2.OpenTime_2nd_card;			byte_counter += sizeof(__int32); // open time 2nd card
		*output_lmf << TDC8PCI2.WriteEmptyEvents_2nd_card;	byte_counter += sizeof(__int32); // write empty events 2nd card
		*output_lmf << TDC8PCI2.TriggerFallingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger falling edge 2nd card
		*output_lmf << TDC8PCI2.TriggerRisingEdge_2nd_card;	byte_counter += sizeof(__int32); // trigger rising edge 2nd card
	}
	return byte_counter;
}






/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteTDC8HPHeader_LMFV_1_to_7()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;
	double	double_Dummy = 0.;
	__int32		int_Dummy = 0;
	//unsigned __int32 unsigned_int_Dummy = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value

	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	if ((DAQVersion_output >= 20020408 && TDC8HP.UserHeaderVersion >=1) || TDC8HP.UserHeaderVersion >= 4) {
		*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);
	}

	*output_lmf << time_reference; byte_counter += sizeof(__int32);
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	if ((DAQVersion_output >= 20020408 && TDC8HP.UserHeaderVersion >= 1) || TDC8HP.UserHeaderVersion >= 4) {*output_lmf << TDCDataType; byte_counter += sizeof(__int32);}

	if (DAQVersion_output < 20080000) {
		*output_lmf << number_of_channels_output; byte_counter += sizeof(__int32);			// number of channels
		*output_lmf << max_number_of_hits_output; byte_counter += sizeof(__int32);			// number of hits
	}
	if (DAQVersion_output >= 20080000) {
		unsigned __int64 temp;
		temp = number_of_channels_output;	*output_lmf << temp; byte_counter += sizeof(unsigned __int64);			// number of channels
		temp = max_number_of_hits_output;	*output_lmf << temp; byte_counter += sizeof(unsigned __int64);			// number of hits
	}
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	*output_lmf << TDC8HP.no_config_file_read;	byte_counter += sizeof(__int32);	// parameter 60-1
	if (DAQVersion_output < 20080000) {
		unsigned __int32 temp;
		temp =__int32(TDC8HP.RisingEnable_p61);  *output_lmf << temp;	byte_counter += sizeof(__int32);	// parameter 61-1
		temp =__int32(TDC8HP.FallingEnable_p62); *output_lmf << temp;	byte_counter += sizeof(__int32);	// parameter 62-1
	}
	if (DAQVersion_output >= 20080000) {
		*output_lmf << TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
		*output_lmf << TDC8HP.FallingEnable_p62;	byte_counter += sizeof(__int64);	// parameter 62-1
	}
	*output_lmf << TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*output_lmf << TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1
	*output_lmf << TDC8HP.OutputLevel_p65;		byte_counter += sizeof(__int32);	// parameter 65-1
	*output_lmf << TDC8HP.GroupingEnable_p66_output;	byte_counter += sizeof(__int32);	// parameter 66-1
	*output_lmf << TDC8HP.AllowOverlap_p67;		byte_counter += sizeof(__int32);	// parameter 67-1
	*output_lmf << TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);	// parameter 68-1
	*output_lmf << TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);	// parameter 69-1
	*output_lmf << TDC8HP.GroupRangeEnd_p70;	byte_counter += sizeof(double);	// parameter 70-1
	*output_lmf << TDC8HP.ExternalClock_p71;	byte_counter += sizeof(__int32);	// parameter 71-1
	*output_lmf << TDC8HP.OutputRollOvers_p72;	byte_counter += sizeof(__int32);	// parameter 72-1
	*output_lmf << TDC8HP.DelayTap0_p73;		byte_counter += sizeof(__int32);	// parameter 73-1
	*output_lmf << TDC8HP.DelayTap1_p74;		byte_counter += sizeof(__int32);	// parameter 74-1
	*output_lmf << TDC8HP.DelayTap2_p75;		byte_counter += sizeof(__int32);	// parameter 75-1
	*output_lmf << TDC8HP.DelayTap3_p76;		byte_counter += sizeof(__int32);	// parameter 76-1
	*output_lmf << TDC8HP.INL_p80;				byte_counter += sizeof(__int32);	// parameter 80-1
	*output_lmf << TDC8HP.DNL_p81;				byte_counter += sizeof(__int32);	// parameter 81-1

	dummy = __int32(TDC8HP.csConfigFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csConfigFile.c_str(), __int32(TDC8HP.csConfigFile.length()));
	byte_counter += __int32(TDC8HP.csConfigFile.length());

	dummy = __int32(TDC8HP.csINLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csINLFile.c_str(), __int32(TDC8HP.csINLFile.length()));
	byte_counter += __int32(TDC8HP.csINLFile.length());

	dummy = __int32(TDC8HP.csDNLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csDNLFile.c_str(), __int32(TDC8HP.csDNLFile.length()));
	byte_counter += __int32(TDC8HP.csDNLFile.length());

	if (TDC8HP.UserHeaderVersion >= 2 ) {
		*output_lmf << TDC8HP.SyncValidationChannel; byte_counter += sizeof(__int32);
	}
	if (TDC8HP.UserHeaderVersion >= 3 ) {
		*output_lmf << TDC8HP.VHR_25ps; byte_counter += sizeof(bool);
	}

	if (TDC8HP.UserHeaderVersion >= 4) {

		*output_lmf << TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);

		*output_lmf << TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
		*output_lmf << TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
		*output_lmf << TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

		//// write TDCInfo
		*output_lmf << TDC8HP.Number_of_TDCs;	byte_counter += sizeof(__int32);
		for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
			*output_lmf << TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->channelCount;		byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->channelStart;		byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
			*output_lmf << TDC8HP.TDC_info[iCount]->serialNumber;		byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->version;			byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.TDC_info[iCount]->fifoSize;			byte_counter += sizeof(__int32);
			output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
			output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
			*output_lmf << TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
		}

		bool bool_dummy = false;
		*output_lmf << bool_dummy; byte_counter += sizeof(bool);
		*output_lmf << bool_dummy; byte_counter += sizeof(bool);
		*output_lmf << bool_dummy; byte_counter += sizeof(bool);
		*output_lmf << bool_dummy; byte_counter += sizeof(bool);
		*output_lmf << bool_dummy; byte_counter += sizeof(bool);

		*output_lmf << TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC8HP_DriverVersion; byte_counter += sizeof(__int32);	
		*output_lmf << TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32);
		*output_lmf << int_Dummy;				byte_counter += sizeof(__int32);

		*output_lmf << double_Dummy; byte_counter += sizeof(double);
		*output_lmf << double_Dummy; byte_counter += sizeof(double);
		*output_lmf << double_Dummy; byte_counter += sizeof(double);
		*output_lmf << double_Dummy; byte_counter += sizeof(double);
		*output_lmf << double_Dummy; byte_counter += sizeof(double);
	}
	
	return byte_counter;
}












/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteTDC8HPHeader_LMFV_8_to_9()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter	= 0;

	bool			 bool_dummy = false;
	//double			 double_Dummy	= 0.;
	//__int32			 int_Dummy		= 0;
	unsigned __int32 unsigned_int_Dummy = 0;
	__int64			 int64_dummy	= 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value

	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);

	if ((number_of_DAQ_source_strings_output < 0) || (!DAQ_source_strings_output)) number_of_DAQ_source_strings_output = 0;
	*output_lmf << number_of_DAQ_source_strings_output;  byte_counter += sizeof(__int32);

	for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {
		unsigned_int_Dummy = (unsigned __int32)(DAQ_source_strings_output[i]->length());
		*output_lmf << unsigned_int_Dummy;   byte_counter += sizeof(__int32);
		output_lmf->write(DAQ_source_strings_output[i]->c_str(),__int32(DAQ_source_strings_output[i]->length()));
		byte_counter += (unsigned __int32)(DAQ_source_strings_output[i]->length());
	}

	*output_lmf << time_reference; byte_counter += sizeof(__int32);
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	*output_lmf << TDCDataType; byte_counter += sizeof(__int32);

	int64_dummy = number_of_channels_output;
	*output_lmf << int64_dummy; byte_counter += sizeof(__int64);			// number of channels
	int64_dummy = max_number_of_hits_output;
	*output_lmf << int64_dummy; byte_counter += sizeof(__int64);			// number of hits
	
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	bool_dummy = TDC8HP.no_config_file_read ? true: false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 60-1
	*output_lmf << TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
	*output_lmf << TDC8HP.FallingEnable_p62;	byte_counter += sizeof(__int64);	// parameter 62-1
	
	*output_lmf << TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*output_lmf << TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1
	bool_dummy = TDC8HP.OutputLevel_p65 ? true:false;
	*output_lmf << bool_dummy;		byte_counter += sizeof(bool);	// parameter 65-1
	bool_dummy = TDC8HP.GroupingEnable_p66_output ? true : false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 66-1
	bool_dummy = TDC8HP.AllowOverlap_p67 ? true: false;
	*output_lmf << bool_dummy;		byte_counter += sizeof(bool);	// parameter 67-1
	*output_lmf << TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);	// parameter 68-1
	*output_lmf << TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);	// parameter 69-1
	*output_lmf << TDC8HP.GroupRangeEnd_p70;	byte_counter += sizeof(double);	// parameter 70-1
	bool_dummy = TDC8HP.ExternalClock_p71 ? true:false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 71-1
	bool_dummy = TDC8HP.OutputRollOvers_p72 ? true:false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 72-1
	*output_lmf << TDC8HP.DelayTap0_p73;		byte_counter += sizeof(__int32);	// parameter 73-1
	*output_lmf << TDC8HP.DelayTap1_p74;		byte_counter += sizeof(__int32);	// parameter 74-1
	*output_lmf << TDC8HP.DelayTap2_p75;		byte_counter += sizeof(__int32);	// parameter 75-1
	*output_lmf << TDC8HP.DelayTap3_p76;		byte_counter += sizeof(__int32);	// parameter 76-1
	bool_dummy = TDC8HP.INL_p80 ? true : false;
	*output_lmf << bool_dummy;				byte_counter += sizeof(bool);	// parameter 80-1
	bool_dummy = TDC8HP.DNL_p81 ? true:false;
	*output_lmf << bool_dummy;				byte_counter += sizeof(bool);	// parameter 81-1

	dummy = __int32(TDC8HP.csConfigFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csConfigFile.c_str(), __int32(TDC8HP.csConfigFile.length()));
	byte_counter += __int32(TDC8HP.csConfigFile.length());

	dummy = __int32(TDC8HP.csINLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csINLFile.c_str(), __int32(TDC8HP.csINLFile.length()));
	byte_counter += __int32(TDC8HP.csINLFile.length());

	dummy = __int32(TDC8HP.csDNLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csDNLFile.c_str(), __int32(TDC8HP.csDNLFile.length()));
	byte_counter += __int32(TDC8HP.csDNLFile.length());

	*output_lmf << TDC8HP.SyncValidationChannel; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.VHR_25ps; byte_counter += sizeof(bool);
	
	*output_lmf << TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);

	*output_lmf << TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
	*output_lmf << TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
	*output_lmf << TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

	//// write TDCInfo
	*output_lmf << TDC8HP.Number_of_TDCs;	byte_counter += sizeof(__int32);
	for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
		*output_lmf << TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->channelCount;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->channelStart;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
		*output_lmf << TDC8HP.TDC_info[iCount]->serialNumber;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->version;			byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->fifoSize;			byte_counter += sizeof(__int32);
		output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
		output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
		*output_lmf << TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
	}


	*output_lmf << TDC8HP.number_of_bools; byte_counter += sizeof(__int32);
	for (__int32 i=4; i<TDC8HP.number_of_bools; ++i) {
		bool bdummy = false;
		*output_lmf << bdummy; byte_counter += sizeof(bool);
	}

	if (TDC8HP.number_of_int32s < 4) TDC8HP.number_of_int32s = 4;
	if (TDC8HP.UserHeaderVersion >= 6) TDC8HP.number_of_int32s = 5;
	*output_lmf << TDC8HP.number_of_int32s; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.TDC8HP_DriverVersion; byte_counter += sizeof(__int32);	
	*output_lmf << TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32);

	if (TDC8HP.UserHeaderVersion == 5) {
		for (__int32 i=4; i<TDC8HP.number_of_int32s; ++i) {
			__int32 idummy = 0;
			*output_lmf << idummy; byte_counter += sizeof(__int32);
		}
	} else {
		*output_lmf << TDC8HP.BinsizeType;	byte_counter += sizeof(__int32);
		for (__int32 i=5; i<TDC8HP.number_of_int32s; ++i) {
			__int32 idummy = 0;
			*output_lmf << idummy; byte_counter += sizeof(__int32);
		}
	}

	if (TDC8HP.UserHeaderVersion == 5) {
		*output_lmf << TDC8HP.number_of_doubles; byte_counter += sizeof(__int32);
		for (__int32 i=0; i<TDC8HP.number_of_doubles; ++i) {
			 double ddummy = 0.;
			*output_lmf << ddummy; byte_counter += sizeof(double);
		}
	}
	if (TDC8HP.UserHeaderVersion >= 6) {
		if (TDC8HP.number_of_doubles == 0) {
			*output_lmf << __int32(1); byte_counter += sizeof(__int32);
			*output_lmf << double(0.); byte_counter += sizeof(double);
		} else {
			*output_lmf << TDC8HP.number_of_doubles; byte_counter += sizeof(__int32);
			*output_lmf << TDC8HP.OffsetTimeZeroChannel_s; byte_counter += sizeof(double);
			for (__int32 i=1; i<TDC8HP.number_of_doubles; ++i) {
				 double ddummy = 0.;
				*output_lmf << ddummy; byte_counter += sizeof(double);
			}
		}
	}
	
	return byte_counter;
}







/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteTDC8HPHeader_LMFV_10_to_12()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter	= 0;

	bool			 bool_dummy = false;
	//double			 double_Dummy	= 0.;
	//__int32			 int_Dummy		= 0;
	unsigned __int32 unsigned_int_Dummy = 0;
	__int64			 int64_dummy	= 0;

	number_of_bytes_in_PostEventData = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value

	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += (unsigned __int32)(DAQ_info.length());

	*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);

	if ((number_of_DAQ_source_strings_output < 0) || (!DAQ_source_strings_output)) number_of_DAQ_source_strings_output = 0;
	*output_lmf << number_of_DAQ_source_strings_output;  byte_counter += sizeof(__int32);

	for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {
		unsigned_int_Dummy = (unsigned __int32)(DAQ_source_strings_output[i]->length());
		*output_lmf << unsigned_int_Dummy;   byte_counter += sizeof(__int32);
		output_lmf->write(DAQ_source_strings_output[i]->c_str(),__int32(DAQ_source_strings_output[i]->length()));
		byte_counter += (unsigned __int32)(DAQ_source_strings_output[i]->length());
	}

	*output_lmf << time_reference; byte_counter += sizeof(__int32);
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	*output_lmf << TDCDataType; byte_counter += sizeof(__int32);

	int64_dummy = number_of_channels_output;
	*output_lmf << int64_dummy; byte_counter += sizeof(__int64);			// number of channels
	int64_dummy = max_number_of_hits_output;
	*output_lmf << int64_dummy; byte_counter += sizeof(__int64);			// number of hits
	
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	bool_dummy = TDC8HP.no_config_file_read ? true: false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 60-1
	*output_lmf << TDC8HP.RisingEnable_p61;		byte_counter += sizeof(__int64);	// parameter 61-1
	*output_lmf << TDC8HP.FallingEnable_p62;	byte_counter += sizeof(__int64);	// parameter 62-1
	
	*output_lmf << TDC8HP.TriggerEdge_p63;		byte_counter += sizeof(__int32);	// parameter 63-1
	*output_lmf << TDC8HP.TriggerChannel_p64;	byte_counter += sizeof(__int32);	// parameter 64-1
	bool_dummy = TDC8HP.OutputLevel_p65 ? true:false;
	*output_lmf << bool_dummy;		byte_counter += sizeof(bool);	// parameter 65-1
	bool_dummy = TDC8HP.GroupingEnable_p66_output ? true : false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 66-1
	bool_dummy = TDC8HP.AllowOverlap_p67 ? true: false;
	*output_lmf << bool_dummy;		byte_counter += sizeof(bool);	// parameter 67-1
	*output_lmf << TDC8HP.TriggerDeadTime_p68;	byte_counter += sizeof(double);	// parameter 68-1
	*output_lmf << TDC8HP.GroupRangeStart_p69;	byte_counter += sizeof(double);	// parameter 69-1
	*output_lmf << TDC8HP.GroupRangeEnd_p70;	byte_counter += sizeof(double);	// parameter 70-1
	bool_dummy = TDC8HP.ExternalClock_p71 ? true:false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 71-1
	bool_dummy = TDC8HP.OutputRollOvers_p72 ? true:false;
	*output_lmf << bool_dummy;	byte_counter += sizeof(bool);	// parameter 72-1
	*output_lmf << TDC8HP.DelayTap0_p73;		byte_counter += sizeof(__int32);	// parameter 73-1
	*output_lmf << TDC8HP.DelayTap1_p74;		byte_counter += sizeof(__int32);	// parameter 74-1
	*output_lmf << TDC8HP.DelayTap2_p75;		byte_counter += sizeof(__int32);	// parameter 75-1
	*output_lmf << TDC8HP.DelayTap3_p76;		byte_counter += sizeof(__int32);	// parameter 76-1
	bool_dummy = TDC8HP.INL_p80 ? true : false;
	*output_lmf << bool_dummy;				byte_counter += sizeof(bool);	// parameter 80-1
	bool_dummy = TDC8HP.DNL_p81 ? true:false;
	*output_lmf << bool_dummy;				byte_counter += sizeof(bool);	// parameter 81-1

	dummy = __int32(TDC8HP.csConfigFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csConfigFile.c_str(), __int32(TDC8HP.csConfigFile.length()));
	byte_counter += __int32(TDC8HP.csConfigFile.length());

	dummy = __int32(TDC8HP.csINLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csINLFile.c_str(), __int32(TDC8HP.csINLFile.length()));
	byte_counter += __int32(TDC8HP.csINLFile.length());

	dummy = __int32(TDC8HP.csDNLFile.length());
	*output_lmf << dummy;	byte_counter += sizeof(__int32);
	output_lmf->write(TDC8HP.csDNLFile.c_str(), __int32(TDC8HP.csDNLFile.length()));
	byte_counter += __int32(TDC8HP.csDNLFile.length());

	*output_lmf << TDC8HP.SyncValidationChannel; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.VHR_25ps; byte_counter += sizeof(bool);
	
	*output_lmf << TDC8HP.GroupTimeOut;  byte_counter += sizeof(double);

	*output_lmf << TDC8HP.SSEEnable;  byte_counter += sizeof(bool);
	*output_lmf << TDC8HP.MMXEnable;  byte_counter += sizeof(bool);
	*output_lmf << TDC8HP.DMAEnable;  byte_counter += sizeof(bool);

	//// write TDCInfo
	*output_lmf << TDC8HP.Number_of_TDCs;	byte_counter += sizeof(__int32);
	for(__int32 iCount=0;iCount<TDC8HP.Number_of_TDCs;++iCount) {
		*output_lmf << TDC8HP.TDC_info[iCount]->index;				byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->channelCount;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->channelStart;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelCount;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->highResChannelStart;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelCount;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->lowResChannelStart;	byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->resolution;			byte_counter += sizeof(double);
		*output_lmf << TDC8HP.TDC_info[iCount]->serialNumber;		byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->version;			byte_counter += sizeof(__int32);
		*output_lmf << TDC8HP.TDC_info[iCount]->fifoSize;			byte_counter += sizeof(__int32);
		output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->INLCorrection,sizeof(__int32)*8*1024);		byte_counter += sizeof(__int32)*8*1024;
		output_lmf->write((__int8*)TDC8HP.TDC_info[iCount]->DNLData,sizeof(unsigned __int16)*8*1024);	byte_counter += sizeof(unsigned __int16)*8*1024;
		*output_lmf << TDC8HP.TDC_info[iCount]->flashValid;			byte_counter += sizeof(bool);
	}


	TDC8HP.number_of_bools = 0;
	*output_lmf << TDC8HP.number_of_bools; byte_counter += sizeof(__int32);

	TDC8HP.number_of_int32s = 5;
	*output_lmf << TDC8HP.number_of_int32s; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.i32NumberOfDAQLoops;	byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.TDC8HP_DriverVersion; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.iTriggerChannelMask;	byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.iTime_zero_channel;	byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.BinsizeType;	byte_counter += sizeof(__int32);

	TDC8HP.number_of_doubles = 1;
	*output_lmf << TDC8HP.number_of_doubles; byte_counter += sizeof(__int32);
	*output_lmf << TDC8HP.OffsetTimeZeroChannel_s; byte_counter += sizeof(double);
	for (__int32 i=1; i<TDC8HP.number_of_doubles; ++i) {
		 double ddummy = 0.;
		*output_lmf << ddummy; byte_counter += sizeof(double);
	}

	return byte_counter;
}











/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteHM1Header()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;
	//__int32 int_Dummy = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = (unsigned __int32)(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += __int32(DAQ_info.length());

	if (DAQVersion_output >= 20020408 && HM1.use_normal_method) {*output_lmf << LMF_Version_output; byte_counter += sizeof(__int32);}

	*output_lmf << system_timeout_output; byte_counter += sizeof(__int32);		//   system time-out
	*output_lmf << time_reference_output; byte_counter += sizeof(__int32);
	if (DAQ_ID_output == DAQ_ID_HM1 || DAQ_ID_output == DAQ_ID_HM1_ABM) {
		*output_lmf << HM1.FAK_DLL_Value; byte_counter += sizeof(__int32);
		*output_lmf << HM1.Resolution_Flag; byte_counter += sizeof(__int32);
		*output_lmf << HM1.trigger_mode_for_start; byte_counter += sizeof(__int32);
		*output_lmf << HM1.trigger_mode_for_stop; byte_counter += sizeof(__int32);
	}
	*output_lmf << tdcresolution_output; byte_counter += sizeof(double);		// tdc resolution in ns

	TDCDataType = 1;
	if (DAQVersion_output >= 20020408 && HM1.use_normal_method) {*output_lmf << TDCDataType; byte_counter += sizeof(__int32);}
	
	if (DAQ_ID_output == DAQ_ID_HM1 || DAQ_ID_output == DAQ_ID_HM1_ABM) {
		*output_lmf << HM1.Even_open_time; byte_counter += sizeof(__int32);
		*output_lmf << HM1.Auto_Trigger; byte_counter += sizeof(__int32);
	}

	if (DAQVersion_output >= 20080507) {
		unsigned __int64 dummy = number_of_channels_output;
		*output_lmf << dummy; byte_counter += sizeof(__int64);			// number of channels
		if (DAQ_ID_output == DAQ_ID_HM1_ABM) {
			dummy = 0;
			*output_lmf << dummy; byte_counter += sizeof(__int64);
		} else { 
			dummy = max_number_of_hits_output;
			*output_lmf << dummy; byte_counter += sizeof(__int64);		// number of hits
		}
	} else {
	*output_lmf << number_of_channels_output; byte_counter += sizeof(__int32);			// number of channels
	if (DAQ_ID_output == DAQ_ID_HM1_ABM) {
			*output_lmf << __int32(0); byte_counter += sizeof(__int32);
	} else { 
		*output_lmf << max_number_of_hits_output; byte_counter += sizeof(__int32);			// number of hits
	}
	}
	
	if (DAQ_ID_output == DAQ_ID_HM1 || DAQ_ID_output == DAQ_ID_HM1_ABM) {
		*output_lmf << HM1.set_bits_for_GP1; byte_counter += sizeof(__int32);
	}
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);				// data format (2=short integer)

	if (DAQ_ID_output == DAQ_ID_HM1 || DAQ_ID_output == DAQ_ID_HM1_ABM) {*output_lmf << module_2nd;	byte_counter += sizeof(__int32);}	// indicator for 2nd module data

	if (DAQ_ID_output == DAQ_ID_2HM1) {
		*output_lmf << DAQSubVersion;						byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_FAK_DLL_Value;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_Resolution_Flag;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_trigger_mode_for_start;	byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_trigger_mode_for_stop;	byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_res_adjust;				byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_tdcresolution;			byte_counter += sizeof(double);
		*output_lmf << HM1.TWOHM1_test_overflow;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_number_of_channels;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_number_of_hits;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_set_bits_for_GP1;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_HM1_ID_1;					byte_counter += sizeof(__int32);
		*output_lmf << HM1.TWOHM1_HM1_ID_2;					byte_counter += sizeof(__int32);
	}

	if (DAQ_ID_output == DAQ_ID_HM1_ABM) {
		*output_lmf << HM1.ABM_m_xFrom;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_xTo;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_yFrom;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_yTo;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_xMin;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_xMax;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_yMin;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_yMax;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_xOffset;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_yOffset;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_m_zOffset;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_Mode;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_OsziDarkInvert;	byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_ErrorHisto;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_XShift;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_YShift;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_ZShift;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_ozShift;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_wdShift;			byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_ucLevelXY;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_ucLevelZ;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_uiABMXShift;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_uiABMYShift;		byte_counter += sizeof(__int32);
		*output_lmf << HM1.ABM_uiABMZShift;		byte_counter += sizeof(__int32);
	}

	return byte_counter;
}







/////////////////////////////////////////////////////////////////
__int32	LMF_IO::WriteCAMACHeader()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32 byte_counter;
	byte_counter = 0;

	*output_lmf << frequency;	byte_counter += sizeof(double);		// frequency is always 4th value
	*output_lmf << IOaddress;	byte_counter += sizeof(__int32);		// IO address (parameter 1) always 5th value
	*output_lmf << timestamp_format_output;	byte_counter += sizeof(__int32);		// TimeInfo (parameter 2) always 6th value  (0,1,2)*32Bit

	unsigned __int32 dummy = __int32(DAQ_info.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);		// Length of DAQInfo always 7th value
	output_lmf->write(DAQ_info.c_str(), __int32(DAQ_info.length()));	// DAQInfo always 8th value
	byte_counter += __int32(DAQ_info.length());

	dummy = __int32(Camac_CIF.length());
	*output_lmf << dummy;
	byte_counter += sizeof(__int32);
	output_lmf->write(Camac_CIF.c_str(), __int32(Camac_CIF.length()));
	byte_counter += __int32(Camac_CIF.length());

	*output_lmf << system_timeout_output; byte_counter += sizeof(__int32);		// system time-out
	*output_lmf << time_reference_output; byte_counter += sizeof(__int32);
	*output_lmf << data_format_in_userheader_output;	byte_counter += sizeof(__int32);		// data format (2=short integer)

	return byte_counter;
}





bool LMF_IO::OpenOutputLMF(std::string LMF_Filename)
{
	return OpenOutputLMF((__int8*)LMF_Filename.c_str());
}



/////////////////////////////////////////////////////////////////
bool LMF_IO::OpenOutputLMF(__int8 * LMF_Filename)
/////////////////////////////////////////////////////////////////
{
	//double				double_Dummy = 0.;
	//unsigned __int32	unsigned_int_Dummy = 0;
	//__int32				int_Dummy = 0;

	if (Parameter_old) memset(Parameter_old,0,10000*sizeof(double));
	for (__int32 i=901;i<=932;i++) Parameter_old[i] = -1.e201;

	if (OutputFileIsOpen) {
		errorflag = 12; // file is already open
		return false;
	}
	output_lmf = new MyFILE(false);

	fADC8.at_least_1_signal_was_written = false;

	output_lmf->open(LMF_Filename);

	if (output_lmf->error) {
		errorflag = 11; // could not open output file
		return false;
	}

//	out_ar = new CArchive(output_lmf,CArchive::store);
//	if (!out_ar) {
//		errorflag = 13; // could not connect CAchrive to output file
//		output_lmf->Close(); output_lmf = 0;
//		return false;
//	}

	if (number_of_DAQ_source_strings_output == -1 && number_of_DAQ_source_strings > 0) {
		number_of_DAQ_source_strings_output = number_of_DAQ_source_strings;
		DAQ_source_strings_output = new std::string*[number_of_DAQ_source_strings];
		memset(DAQ_source_strings_output,0,sizeof(std::string*)*number_of_DAQ_source_strings_output);
		for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {DAQ_source_strings_output[i] = new std::string(); *DAQ_source_strings_output[i] = *DAQ_source_strings[i];}
	} else number_of_DAQ_source_strings_output = 0;

	if (number_of_CCFHistory_strings_output == -1 && number_of_CCFHistory_strings > 0) {
		number_of_CCFHistory_strings_output = number_of_CCFHistory_strings;
		CCFHistory_strings_output = new std::string*[number_of_CCFHistory_strings];
		memset(CCFHistory_strings_output,0,sizeof(std::string*)*number_of_CCFHistory_strings_output);
		for (__int32 i=0;i<number_of_CCFHistory_strings;++i)  {CCFHistory_strings_output[i] = new std::string(); *CCFHistory_strings_output[i] = *CCFHistory_strings[i];}
	} else number_of_CCFHistory_strings_output = 0;

	if (number_of_DAN_source_strings_output == -1 && number_of_DAN_source_strings > 0) {
		number_of_DAN_source_strings_output = number_of_DAN_source_strings;
		DAN_source_strings_output = new std::string*[number_of_DAN_source_strings];
		memset(DAN_source_strings_output,0,sizeof(std::string*)*number_of_DAN_source_strings_output);
		for (__int32 i=0;i<number_of_DAN_source_strings;++i) {DAN_source_strings_output[i] = new std::string(); *DAN_source_strings_output[i] = *DAN_source_strings[i];}
	} else number_of_DAN_source_strings_output = 0;

	if (DAQ_ID_output == DAQ_ID_RAW32BIT) {
		if (number_of_channels_output  == 0) number_of_channels_output  = number_of_channels;
		if (max_number_of_hits_output  == 0) max_number_of_hits_output  = max_number_of_hits;
		if (number_of_channels_output == 0 || max_number_of_hits_output == 0) {errorflag =14; return false;}
		Numberofcoordinates_output = number_of_channels_output*(1+max_number_of_hits_output);
		data_format_in_userheader_output = 10;

		errorflag = 0; // no error
		OutputFileIsOpen = true;
		return true;
	}

	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		if (number_of_channels_output  == 0) number_of_channels_output  = number_of_channels;
		if (max_number_of_hits_output  == 0) max_number_of_hits_output  = max_number_of_hits;
		if (number_of_channels_output == 0 || max_number_of_hits_output == 0) {errorflag =14; return false;}
		Numberofcoordinates_output = number_of_channels_output*(1+max_number_of_hits_output);
		if (SIMPLE_DAQ_ID_Orignial == 0) SIMPLE_DAQ_ID_Orignial = DAQ_ID;
		*output_lmf << SIMPLE_DAQ_ID_Orignial;
		unsigned __int32 dummy = (unsigned __int32)(uint64_number_of_written_events);
		*output_lmf << dummy;
		*output_lmf << data_format_in_userheader_output;

		errorflag = 0; // no error
		OutputFileIsOpen = true;
		return true;
	}


	
// Preparing to write LMF-header:
// -----------------------------------
	OutputFilePathName = LMF_Filename;

	Headersize_output = 0;

	if (DAQVersion_output == -1) DAQVersion_output = DAQVersion;

	if (DAQVersion_output >= 20080000 && Cobold_Header_version_output == 0) Cobold_Header_version_output = 2008;
	if (Cobold_Header_version_output == 0) Cobold_Header_version_output = Cobold_Header_version;

	if (tdcresolution_output < 0.) tdcresolution_output = tdcresolution;
	if (Starttime_output == 0)	{
		Starttime_output = Starttime;
	}
	if (Stoptime_output == 0)	{
		Stoptime_output = Stoptime;
	}
	if (system_timeout_output == -1) system_timeout_output = system_timeout;
	if (time_reference_output ==  0) time_reference_output = time_reference;
	if (common_mode_output    == -1) common_mode_output    = common_mode;

	if (number_of_channels_output  == 0) number_of_channels_output  = number_of_channels;
	if (max_number_of_hits_output  == 0) max_number_of_hits_output  = max_number_of_hits;
	
	if (data_format_in_userheader_output == -2) data_format_in_userheader_output = data_format_in_userheader;
	if (DAQ_ID_output == 0) DAQ_ID_output = DAQ_ID;
	if (DAQ_ID_output == DAQ_ID_TDC8HPRAW) DAQ_ID_output = DAQ_ID_TDC8HP;

	if (DAQ_ID_output == DAQ_ID_2TDC8) {
		if (number_of_channels2_output == -1) number_of_channels2_output = number_of_channels2;
		if (max_number_of_hits2_output == -1) max_number_of_hits2_output = max_number_of_hits2;
	}

	if (max_number_of_hits2_output < __int32(0) || max_number_of_hits2_output > __int32(100000)) max_number_of_hits2_output = max_number_of_hits_output;

	if (timestamp_format_output == -1) timestamp_format_output = timestamp_format;

	if (Numberofcoordinates_output == -2) {
		if (data_format_in_userheader_output == LM_SHORT)  Numberofcoordinates_output = timestamp_format_output * 2;
		if (data_format_in_userheader_output == LM_DOUBLE) Numberofcoordinates_output = timestamp_format_output == 0 ? 0 : 1;
		if (data_format_in_userheader_output == LM_SLONG)  Numberofcoordinates_output = timestamp_format_output;
		if (data_format_in_userheader_output == LM_SHORT || data_format_in_userheader_output == LM_DOUBLE || data_format_in_userheader_output == LM_SLONG) {
			Numberofcoordinates_output += number_of_channels_output  * (max_number_of_hits_output +1);
			Numberofcoordinates_output += number_of_channels2_output * (max_number_of_hits2_output+1);
		}
		if (data_format_in_userheader_output == LM_CAMAC) Numberofcoordinates_output = Numberofcoordinates;

		if (data_format_in_userheader_output == LM_USERDEF) {
			if (DAQVersion_output < 20080000) {this->errorflag = 17; return false;}
			if (DAQ_ID_output == DAQ_ID_TDC8HP) {
				//Numberofcoordinates_output = 2 + timestamp_format_output + number_of_channels_output*(1+max_number_of_hits_output); // old
				Numberofcoordinates_output = number_of_channels_output*(1+max_number_of_hits_output);
				if (LMF_Version >=9) Numberofcoordinates_output++; // for level info
			}
			if (DAQ_ID_output == DAQ_ID_TDC8 || DAQ_ID_output == DAQ_ID_2TDC8) {
				//Numberofcoordinates_output = 2 + timestamp_format_output*2 + (number_of_channels_output + number_of_channels2_output)*(1+max_number_of_hits_output); // old
				Numberofcoordinates_output = (number_of_channels_output + number_of_channels2_output)*(1+max_number_of_hits_output);
			}
			if (DAQ_ID_output == DAQ_ID_HM1) {
				//Numberofcoordinates_output = 2 + timestamp_format_output*2 + (number_of_channels_output + number_of_channels2_output)*(1+max_number_of_hits_output); // old
				Numberofcoordinates_output = (number_of_channels_output + number_of_channels2_output)*(1+max_number_of_hits_output);
			}
			if (DAQ_ID_output == DAQ_ID_FADC8) {
				Numberofcoordinates_output = 3*number_of_channels_output*(1+max_number_of_hits_output);
			}
			if (DAQ_ID_output == DAQ_ID_FADC4) {
				return false;
			}
		}
	}



//  WRITE LMF-HEADER:


	WriteFirstHeader();

	output_lmf->flush();

	Headersize_output = (unsigned __int32)(output_lmf->tell());

	unsigned __int64 seek_value;
	if (Cobold_Header_version_output <= 2002) seek_value = 3*sizeof(unsigned __int32);
	if (Cobold_Header_version_output >= 2008) seek_value = 2*sizeof(unsigned __int32) + sizeof(unsigned __int64);
	output_lmf->seek(seek_value);


	if (Cobold_Header_version_output <= 2002) *output_lmf << Headersize_output;
	if (Cobold_Header_version_output >= 2008) {
		unsigned __int64 temp = Headersize_output;
		*output_lmf << temp;
	}

	output_lmf->flush();
	output_lmf->seek(Headersize_output);

	if (LMF_Version_output == -1) {
		LMF_Version_output = LMF_Version;
		if (LMF_Version_output == -1) LMF_Version_output = 12; // XXX if necessary: modify to latest LMF version number
	}

	output_byte_counter = 0;

//  WRITE USER-HEADER
	if (Cobold_Header_version_output >= 2008 || DAQVersion_output >= 20080000) {
		*output_lmf << LMF_Header_version;
		output_byte_counter += sizeof(__int32);
	}
	if (Cobold_Header_version_output <= 2002) {
		*output_lmf << User_header_size_output;
		output_byte_counter += sizeof(__int32);
	}
	if (Cobold_Header_version_output >= 2008) {
		unsigned __int64 temp = User_header_size_output;
		*output_lmf << temp;
		output_byte_counter += sizeof(unsigned __int64);
	}

	*output_lmf << DAQVersion_output;			output_byte_counter += sizeof(__int32);	// Version is always 2nd value
	*output_lmf << DAQ_ID_output;		output_byte_counter += sizeof(__int32);	// DAQ_ID is always 3ed value

	if (DAQ_ID_output == DAQ_ID_TDC8)	 output_byte_counter += WriteTDC8PCI2Header();
	if (DAQ_ID_output == DAQ_ID_2TDC8)	 output_byte_counter += Write2TDC8PCI2Header();
	if (DAQ_ID_output == DAQ_ID_TDC8HP || DAQ_ID_output == DAQ_ID_TDC8HPRAW)	 {
		if (this->LMF_Version_output <  8) output_byte_counter += WriteTDC8HPHeader_LMFV_1_to_7();
		if (this->LMF_Version_output >= 8 && this->LMF_Version_output <= 9) output_byte_counter += WriteTDC8HPHeader_LMFV_8_to_9();
		if (this->LMF_Version_output >= 10) output_byte_counter += WriteTDC8HPHeader_LMFV_10_to_12();
	}

	if (DAQ_ID_output == DAQ_ID_HM1)	 output_byte_counter += WriteHM1Header();
	if (DAQ_ID_output == DAQ_ID_HM1_ABM) output_byte_counter += WriteHM1Header();
	if (DAQ_ID_output == DAQ_ID_CAMAC)   output_byte_counter += WriteCAMACHeader();
	if (DAQ_ID_output == DAQ_ID_TCPIP)   output_byte_counter += WriteTCPIPHeader();
	if (DAQ_ID_output == DAQ_ID_FADC8)   output_byte_counter += WritefADC8_header_LMFversion10();

	User_header_size_output = output_byte_counter;

	errorflag = 0; // no error
	OutputFileIsOpen = true;
	return true;
}











/////////////////////////////////////////////////////////////////
void LMF_IO::WriteFirstHeader()
/////////////////////////////////////////////////////////////////
{
	if (Cobold_Header_version_output >= 2008) {
		unsigned __int32 ArchiveFlagtemp = 476759;

		ArchiveFlagtemp = ArchiveFlagtemp | DAN_SOURCE_CODE;
		ArchiveFlagtemp = ArchiveFlagtemp | DAQ_SOURCE_CODE; 
		ArchiveFlagtemp = ArchiveFlagtemp | CCF_HISTORY_CODE; 

		*output_lmf << ArchiveFlagtemp;
		*output_lmf << data_format_in_userheader_output;

		unsigned __int64 temp;
		temp = Numberofcoordinates_output;  *output_lmf << temp;
		temp = Headersize_output;			*output_lmf << temp;
		temp = User_header_size_output;		*output_lmf << temp;
		*output_lmf << uint64_number_of_written_events;
	}

	if (Cobold_Header_version_output <= 2002) {
		unsigned __int32 ArchiveFlagtemp = 476758;
		*output_lmf << ArchiveFlagtemp;
		*output_lmf << data_format_in_userheader_output;

		*output_lmf << Numberofcoordinates_output;
		*output_lmf << Headersize_output;
		*output_lmf << User_header_size_output;
		unsigned __int32 dummy = (unsigned __int32)(uint64_number_of_written_events);
		*output_lmf << dummy;
	}

	write_times(output_lmf,Starttime_output,Stoptime_output);

	Write_StdString_as_CString(*output_lmf,Versionstring);
	Write_StdString_as_CString(*output_lmf,OutputFilePathName);
	Write_StdString_as_CString(*output_lmf,Comment_output);

	if (CCFHistory_strings_output) {
		if (number_of_CCFHistory_strings_output > 0) {
			*output_lmf << number_of_CCFHistory_strings_output;
			for (__int32 i=0;i<number_of_CCFHistory_strings_output;++i) {
			unsigned __int32 unsigned_int_Dummy = (unsigned __int32)(CCFHistory_strings_output[i]->length());
				*output_lmf << unsigned_int_Dummy;					
			output_lmf->write(CCFHistory_strings_output[i]->c_str(),__int32(CCFHistory_strings_output[i]->length()));
			}
		} else *output_lmf << __int32(0);
	} else {
		number_of_CCFHistory_strings_output = 0;
		if (Cobold_Header_version_output >= 2008) *output_lmf << number_of_CCFHistory_strings_output;
	}

	if (DAN_source_strings_output) {
		if (number_of_DAN_source_strings_output > 0) {
			*output_lmf << number_of_DAN_source_strings_output;
			for (__int32 i=0;i<number_of_DAN_source_strings_output;++i) {
			unsigned __int32 unsigned_int_Dummy = (unsigned __int32)(DAN_source_strings_output[i]->length());
				*output_lmf << unsigned_int_Dummy;
			output_lmf->write(DAN_source_strings_output[i]->c_str(),__int32(DAN_source_strings_output[i]->length()));
			}
		} else *output_lmf << __int32(0);
	} else {
		number_of_DAN_source_strings_output = 0;
		if (Cobold_Header_version_output >= 2008) *output_lmf << number_of_DAN_source_strings_output;
	}
}














/////////////////////////////////////////////////////////////////
void LMF_IO::WriteEventHeader(unsigned __int64 timestamp, unsigned __int32 cnt[])
/////////////////////////////////////////////////////////////////
{
	unsigned __int64 HeaderLength = 0;
	// EventLength information in 64Bit
	if (DAQVersion_output >= 20080000 || data_format_in_userheader_output == LM_USERDEF) {
		HeaderLength = timestamp_format_output + 1 + number_of_channels_output;	// +1 for HeaderLength itself

		HeaderLength = timestamp_format_output*sizeof(__int32) + 2*sizeof(__int64);	// +2 for HeaderLength itself, +2 for EventCounter (size in __int32)
		for(__int32 iCount=0;iCount<number_of_channels_output;++iCount) HeaderLength += cnt[iCount]*sizeof(__int32);
		HeaderLength += number_of_channels_output*sizeof(__int16);
#ifdef LINUX
		HeaderLength = (HeaderLength & 0x00ffffffffffffffLL) | 0xff00000000000000LL;	// set 0xff in bits 63..56 as EventMarker
#else
		HeaderLength = (HeaderLength & 0x00ffffffffffffff) | 0xff00000000000000;	// set 0xff in bits 63..56 as EventMarker
#endif
		*output_lmf << HeaderLength;
		*output_lmf << uint64_number_of_written_events; 
	}

	myLARGE_INTEGER LARGE_Timestamp;
	LARGE_Timestamp.QuadPart = timestamp;

	if (timestamp_format_output > 0) {
		*output_lmf << LARGE_Timestamp.LowPart;
		if (timestamp_format_output == 2) *output_lmf << LARGE_Timestamp.HighPart;
	}
}













/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(unsigned __int64 timestamp, unsigned __int32 cnt[], __int32 * i32TDC)
/////////////////////////////////////////////////////////////////
{
	unsigned __int16 dummy_uint16;
	__int32 dummy_int32;
	double dummy_double;

	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	
	++uint64_number_of_written_events;

	WriteEventHeader(timestamp,cnt);

	__int32 i,j;
	if (DAQ_ID_output != DAQ_ID_SIMPLE) {
		for (i=0;i<number_of_channels_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits_output) hits = max_number_of_hits_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(i32TDC[i*num_ions+j]);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = (unsigned __int16)(0);
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) {
					dummy_double = double(i32TDC[i*num_ions+j]);
					*output_lmf << dummy_double;
				}
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				*output_lmf << __int32(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));
				for (j=0;j<hits;++j) *output_lmf << i32TDC[i*num_ions+j];
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_int32;
			}
			if (data_format_in_userheader_output == LM_USERDEF) {
				dummy_uint16 = (unsigned __int16)(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));
				*output_lmf << dummy_uint16;
				if (DAQ_ID_output == DAQ_ID_2TDC8 || DAQ_ID_output == DAQ_ID_TDC8 || DAQ_ID_output == DAQ_ID_HM1) {
					for (j=0;j<hits;++j) {
						dummy_uint16 = (unsigned __int16)(i32TDC[i*num_ions+j]);
						*output_lmf << dummy_uint16;
					}
				} else {
				for (j=0;j<hits;++j) *output_lmf << i32TDC[i*num_ions+j];
				}
			}
		}
				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 9) {
					*output_lmf << ui64LevelInfo;
				}

		if (this->LMF_Version_output >= 11) {
			unsigned __int32 changed_mask = 0;
			//__int32 max_par_index = 932;
			for (__int32 i=0; i<32; i++) {
				if (Parameter_old[i+901] != Parameter[i+901]) {
					changed_mask += (1 << i);
				}
			}

			*output_lmf << changed_mask;

			for (__int32 i=0; i<32; i++) {
				if (Parameter_old[i+901] != Parameter[i+901]) {
					*output_lmf << Parameter[i+901];
					Parameter_old[i+901] = Parameter[i+901];
				}
			}
		}

				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 10) {
					*output_lmf << number_of_bytes_in_PostEventData;
					for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) *output_lmf << ui8_PostEventData[i];
		}
	}
	if (DAQ_ID_output == DAQ_ID_2TDC8 && number_of_channels2_output > 0) { // here we write only the stuff for the second TDC (because the first TDC was written above)
		for (i=number_of_channels_output;i<number_of_channels2_output+number_of_channels_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits2_output) hits = max_number_of_hits2_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));  *output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(i32TDC[i*num_ions+j]);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));		*output_lmf << dummy_double;
				for (j=0;j<hits;++j) {
					dummy_double = double(i32TDC[i*num_ions+j]);
					*output_lmf << dummy_double;
			}
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				dummy_int32 = __int32(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));		*output_lmf << dummy_int32;
				for (j=0;j<hits;++j) *output_lmf << i32TDC[i*num_ions+j];
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_int32;
			}
			if (data_format_in_userheader_output == LM_USERDEF) {
				dummy_uint16 = (unsigned __int16)(hits + (DAQ_ID == DAQ_ID_HM1 ? 1 : 0));
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(i32TDC[i*num_ions+j]);
					*output_lmf << dummy_uint16;
				}
			}
		}
	}

	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		unsigned __int32 channel;
		unsigned __int32 i;
		i = 0;
		for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) i = i + cnt[channel] + (cnt[channel]>0 ? 1 : 0);
		if (data_format_in_userheader_output == 2)  {dummy_uint16 = (unsigned __int16)i; *output_lmf << dummy_uint16;}
		if (data_format_in_userheader_output == 10) *output_lmf << i;
		if (data_format_in_userheader_output == 2) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_uint16 = (unsigned __int16)((channel << 8) + cnt[channel]);
					*output_lmf << dummy_uint16;
					for (i=0;i<cnt[channel];++i) {
						dummy_uint16 = (unsigned __int16)(i32TDC[channel*num_ions+i]);
						*output_lmf << dummy_uint16;
					}
				}
			}
		}
		if (data_format_in_userheader_output == 10) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_int32 = __int32((channel << 24) + cnt[channel]);
					*output_lmf << dummy_int32;
					for (i=0;i<cnt[channel];++i) *output_lmf << i32TDC[channel*num_ions+i];
				}
			}
		}
	} // end if (DAQ_ID_output == DAQ_ID_SIMPLE)

	return;
}



/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(double timestamp, unsigned __int32 cnt[], __int32 *i32TDC)
/////////////////////////////////////////////////////////////////
{
	if (!output_lmf|| !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}

	unsigned __int64 new_timestamp = (unsigned __int64)(timestamp * frequency + 0.001);

	WriteTDCData(new_timestamp, cnt, i32TDC);

	return;
}












/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(unsigned __int64 timestamp, unsigned __int32 cnt[], double * d64TDC)
/////////////////////////////////////////////////////////////////
{
	unsigned __int16 dummy_uint16;
	double dummy_double;
	__int32 dummy_int32;

	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	
	++uint64_number_of_written_events;

	WriteEventHeader(timestamp,cnt);

	__int32 i,j;
	__int32 ii;
	if (DAQ_ID_output != DAQ_ID_SIMPLE) {
		for (i=0;i<number_of_channels_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits_output) hits = max_number_of_hits_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(d64TDC[i*num_ions+j]+1.e-6);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) 	*output_lmf << d64TDC[i*num_ions+j];
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				dummy_int32 = __int32(hits);
				*output_lmf << dummy_int32;
				for (j=0;j<hits;++j) {
					if (d64TDC[i*num_ions+j] >= 0.) ii =__int32(d64TDC[i*num_ions+j]+1.e-6);
					if (d64TDC[i*num_ions+j] <  0.) ii =__int32(d64TDC[i*num_ions+j]-1.e-6);
					*output_lmf << ii;
				}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_int32;
			}
			if (data_format_in_userheader_output == LM_USERDEF) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					if (d64TDC[i*num_ions+j] >= 0.) ii =__int32(d64TDC[i*num_ions+j]+1.e-6);
					if (d64TDC[i*num_ions+j] <  0.) ii =__int32(d64TDC[i*num_ions+j]-1.e-6);
					if (DAQ_ID_output != DAQ_ID_HM1 && DAQ_ID_output != DAQ_ID_TDC8 && DAQ_ID_output != DAQ_ID_2TDC8) *output_lmf << ii;
						else *output_lmf << (unsigned __int16)(ii);
				}
			}
		}
				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 9) {
					*output_lmf << ui64LevelInfo;
				}
				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 10) {
					*output_lmf << number_of_bytes_in_PostEventData;
					for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) *output_lmf << ui8_PostEventData[i];
		}
	}
	if (DAQ_ID_output == DAQ_ID_2TDC8 && number_of_channels2_output > 0) {
		for (i=number_of_channels_output;i<number_of_channels_output+number_of_channels2_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits2_output) hits = max_number_of_hits2_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(d64TDC[i*num_ions+j]+1.e-6);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) *output_lmf << d64TDC[i*num_ions+j];
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				*output_lmf <<__int32(hits);
				for (j=0;j<hits;++j) {
					if (d64TDC[i*num_ions+j] >= 0.) ii =__int32(d64TDC[i*num_ions+j]+1.e-6);
					if (d64TDC[i*num_ions+j] <  0.) ii =__int32(d64TDC[i*num_ions+j]-1.e-6);
					*output_lmf << ii;
				}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_int32;
			}
		}
	}

	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		unsigned __int32 channel;
		unsigned __int32 i;
		i = 0;
		for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) i = i + cnt[channel] + (cnt[channel]>0 ? 1 : 0);
		if (data_format_in_userheader_output == 2)  {
			dummy_uint16 = (unsigned __int16)i;
			*output_lmf << dummy_uint16;
		}
		if (data_format_in_userheader_output == 10) *output_lmf << i;
		if (data_format_in_userheader_output == 2) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_uint16 = (unsigned __int16)((channel << 8) + cnt[channel]);
					*output_lmf << dummy_uint16;
					for (i=0;i < cnt[channel]; ++i) {
						dummy_uint16 = (unsigned __int16)(d64TDC[channel*num_ions+i]);
						*output_lmf << dummy_uint16;
					}
				}
			}
		}
		if (data_format_in_userheader_output == 10) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_int32 = __int32((channel << 24) + cnt[channel]);
					*output_lmf << dummy_int32;
					for (i=0;i<cnt[channel];++i) {
						dummy_int32 = __int32(d64TDC[channel*num_ions+i]);
						*output_lmf << dummy_int32;
					}
				}
			}
		}
	} // end if (DAQ_ID_output == DAQ_ID_SIMPLE)

	return;
}

/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(double timestamp, unsigned __int32 cnt[], double *dtdc)
/////////////////////////////////////////////////////////////////
{
	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	unsigned __int64 new_timestamp = (unsigned __int64)(timestamp * frequency + 0.001);

	WriteTDCData(new_timestamp, cnt, dtdc);

	return;
}






/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(unsigned __int64 timestamp, unsigned __int32 cnt[], __int64 * i64TDC)
	/////////////////////////////////////////////////////////////////
{
	unsigned __int16 dummy_uint16;
	double dummy_double;
	__int32 dummy_int32;

	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}

	++uint64_number_of_written_events;

	WriteEventHeader(timestamp,cnt);

	__int32 i,j;
	__int32 ii;
	if (DAQ_ID_output != DAQ_ID_SIMPLE) {
		for (i=0;i<number_of_channels_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits_output) hits = max_number_of_hits_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(i64TDC[i*num_ions+j]+1.e-6);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) 	*output_lmf << double(i64TDC[i*num_ions+j]);
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				dummy_int32 = __int32(hits);
				*output_lmf << dummy_int32;
				for (j=0;j<hits;++j) {
					if (i64TDC[i*num_ions+j] >= 0.) ii =__int32(i64TDC[i*num_ions+j]+1.e-6);
					if (i64TDC[i*num_ions+j] <  0.) ii =__int32(i64TDC[i*num_ions+j]-1.e-6);
					*output_lmf << ii;
				}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_int32;
			}
			if (data_format_in_userheader_output == LM_USERDEF) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					if (i64TDC[i*num_ions+j] >= 0.) ii =__int32(i64TDC[i*num_ions+j]+1.e-6);
					if (i64TDC[i*num_ions+j] <  0.) ii =__int32(i64TDC[i*num_ions+j]-1.e-6);
					if (DAQ_ID_output != DAQ_ID_HM1 && DAQ_ID_output != DAQ_ID_TDC8 && DAQ_ID_output != DAQ_ID_2TDC8) *output_lmf << ii;
					else *output_lmf << (unsigned __int16)(ii);
				}
			}
		}
		if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 9) {
			*output_lmf << ui64LevelInfo;
		}
		if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 10) {
			*output_lmf << number_of_bytes_in_PostEventData;
			for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) *output_lmf << ui8_PostEventData[i];
		}
	}
	if (DAQ_ID_output == DAQ_ID_2TDC8 && number_of_channels2_output > 0) {
		for (i=number_of_channels_output;i<number_of_channels_output+number_of_channels2_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits2_output) hits = max_number_of_hits2_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_uint16 = (unsigned __int16)(i64TDC[i*num_ions+j]+1.e-6);
					*output_lmf << dummy_uint16;
				}
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) *output_lmf << i64TDC[i*num_ions+j];
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				*output_lmf <<__int32(hits);
				for (j=0;j<hits;++j) {
					if (i64TDC[i*num_ions+j] >= 0.) ii =__int32(i64TDC[i*num_ions+j]+1.e-6);
					if (i64TDC[i*num_ions+j] <  0.) ii =__int32(i64TDC[i*num_ions+j]-1.e-6);
					*output_lmf << ii;
				}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_int32;
			}
		}
	}

	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		unsigned __int32 channel;
		unsigned __int32 i;
		i = 0;
		for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) i = i + cnt[channel] + (cnt[channel]>0 ? 1 : 0);
		if (data_format_in_userheader_output == 2)  {
			dummy_uint16 = (unsigned __int16)i;
			*output_lmf << dummy_uint16;
		}
		if (data_format_in_userheader_output == 10) *output_lmf << i;
		if (data_format_in_userheader_output == 2) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_uint16 = (unsigned __int16)((channel << 8) + cnt[channel]);
					*output_lmf << dummy_uint16;
					for (i=0;i < cnt[channel]; ++i) {
						dummy_uint16 = (unsigned __int16)(i64TDC[channel*num_ions+i]);
						*output_lmf << dummy_uint16;
					}
				}
			}
		}
		if (data_format_in_userheader_output == 10) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_int32 = __int32((channel << 24) + cnt[channel]);
					*output_lmf << dummy_int32;
					for (i=0;i<cnt[channel];++i) {
						dummy_int32 = __int32(i64TDC[channel*num_ions+i]);
						*output_lmf << dummy_int32;
					}
				}
			}
		}
	} // end if (DAQ_ID_output == DAQ_ID_SIMPLE)

	return;
}

/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(double timestamp, unsigned __int32 cnt[], __int64 *i64tdc)
	/////////////////////////////////////////////////////////////////
{
	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	unsigned __int64 new_timestamp = (unsigned __int64)(timestamp * frequency + 0.001);

	WriteTDCData(new_timestamp, cnt, i64tdc);

	return;
}










/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(unsigned __int64 timestamp, unsigned __int32 cnt[], unsigned __int16 * us16TDC)
/////////////////////////////////////////////////////////////////
{
	unsigned __int16 dummy_uint16;
	double dummy_double;
	__int32 dummy_int32;


	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	
	++uint64_number_of_written_events;

	WriteEventHeader(timestamp,cnt);

	__int32 i,j;

	if (DAQ_ID_output == DAQ_ID_HM1_ABM) {
		for (i=0;i<number_of_channels_output;++i) {
			if (data_format_in_userheader_output == 2)  *output_lmf << us16TDC[i*num_ions];
			if (data_format_in_userheader_output == 5)  {dummy_double = double(us16TDC[i*num_ions]); *output_lmf << dummy_double;}
			if (data_format_in_userheader_output == 10) {dummy_int32 = __int32(us16TDC[i*num_ions]); *output_lmf << dummy_int32;}
		}
	}
	if (DAQ_ID_output != DAQ_ID_SIMPLE && DAQ_ID_output != DAQ_ID_HM1_ABM) {
		for (i=0;i<number_of_channels_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits_output) hits = max_number_of_hits_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) *output_lmf << us16TDC[i*num_ions+j];
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) {dummy_double = double(us16TDC[i*num_ions+j]); *output_lmf << dummy_double;}
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				*output_lmf <<__int32(hits);
				for (j=0;j<hits;++j) {dummy_int32 = __int32(us16TDC[i*num_ions+j]); *output_lmf << dummy_int32;}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits_output;++j) *output_lmf << dummy_int32;
			}
			if (data_format_in_userheader_output == LM_USERDEF) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) {
					dummy_int32 = __int32(us16TDC[i*num_ions+j]);
					if (DAQ_ID_output != DAQ_ID_HM1 && DAQ_ID_output != DAQ_ID_TDC8 && DAQ_ID_output != DAQ_ID_2TDC8) *output_lmf << dummy_int32;
						else *output_lmf << (unsigned __int16)(dummy_int32);
				}
			}
		}
				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 9) {
					*output_lmf << ui64LevelInfo;
				}
				if (DAQ_ID_output == DAQ_ID_TDC8HP && this->LMF_Version_output >= 10) {
					*output_lmf << number_of_bytes_in_PostEventData;
					for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) *output_lmf << ui8_PostEventData[i];
		}
	}

	if (DAQ_ID_output == DAQ_ID_2TDC8 && number_of_channels2_output > 0) {
		for (i=number_of_channels_output;i<number_of_channels_output+number_of_channels2_output;++i) {
			__int32 hits = cnt[i];
			if (hits > max_number_of_hits2_output) hits = max_number_of_hits2_output;
			if (data_format_in_userheader_output == 2) {
				dummy_uint16 = (unsigned __int16)(hits);
				*output_lmf << dummy_uint16;
				for (j=0;j<hits;++j) *output_lmf << us16TDC[i*num_ions+j];
				dummy_uint16 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_uint16;
			}
			if (data_format_in_userheader_output == 5) {
				dummy_double = double(hits);
				*output_lmf << dummy_double;
				for (j=0;j<hits;++j) {dummy_double = double(us16TDC[i*num_ions+j]); *output_lmf << dummy_double;}
				dummy_double = 0.;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_double;
			}
			if (data_format_in_userheader_output == 10) {
				dummy_int32 = __int32(hits);
				*output_lmf << dummy_int32;
				for (j=0;j<hits;++j) {dummy_int32 = __int32(us16TDC[i*num_ions+j]); *output_lmf << dummy_int32;}
				dummy_int32 = 0;
				for (j=hits;j<max_number_of_hits2_output;++j) *output_lmf << dummy_int32;
			}
		}
	}

	if (DAQ_ID_output == DAQ_ID_SIMPLE) {
		unsigned __int32 channel;
		unsigned __int32 i = 0;

		for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) i = i + cnt[channel] + (cnt[channel]>0 ? 1 : 0);
		if (data_format_in_userheader_output == 2)  {dummy_uint16 = (unsigned __int16)i; *output_lmf << dummy_uint16;}
		if (data_format_in_userheader_output == 10) *output_lmf << i;
		if (data_format_in_userheader_output == 2) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_uint16 = (unsigned __int16)((channel << 8) + cnt[channel]);
					*output_lmf << dummy_uint16;
					for (i=0;i<cnt[channel];++i) *output_lmf << us16TDC[channel*num_ions+i];
				}
			}
		}
		if (data_format_in_userheader_output == 10) {
			for (channel=0;channel < (unsigned __int32)number_of_channels_output;++channel) {
				if (cnt[channel]>0) {
					dummy_int32 = __int32((channel << 24) + cnt[channel]);
					*output_lmf << dummy_int32;
					for (i=0;i<cnt[channel];++i) {dummy_int32 = __int32(us16TDC[channel*num_ions+i]); *output_lmf << dummy_int32;}
				}
			}
		}
	} // end if (DAQ_ID_output == DAQ_ID_SIMPLE)

	return;
}







/////////////////////////////////////////////////////////////////
void LMF_IO::WriteTDCData(double timestamp, unsigned __int32 cnt[], unsigned __int16 * us16TDC)
/////////////////////////////////////////////////////////////////
{
	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	unsigned __int64 new_timestamp = (unsigned __int64)(timestamp * frequency);

	WriteTDCData(new_timestamp, cnt, us16TDC);

	return;
}








/////////////////////////////////////////////////////////////////
__int32 LMF_IO::GetErrorStatus()
/////////////////////////////////////////////////////////////////
{
	return errorflag;
}







/////////////////////////////////////////////////////////////////
bool LMF_IO::SeekToEventNumber(unsigned __int64 target_number)
/////////////////////////////////////////////////////////////////
{
	if (DAQ_ID == DAQ_ID_SIMPLE) return false;

	if (target_number == 0) {
		input_lmf->seek((unsigned __int64)(Headersize + User_header_size));
		must_read_first = true;
		errorflag = 0;
		input_lmf->error = 0;
		return true;
	}

	if (data_format_in_userheader == LM_USERDEF) {errorflag = 16; return false;}

	if (!input_lmf) {
		errorflag = 9;
		return false;
	}


/*	unsigned __int64 filesize;
	unsigned __int64 pos = input_lmf->tell();
	input_lmf->seek_to_end();
	filesize = input_lmf->tell();
	input_lmf->seek(pos);
*/

	

	if (target_number < 0) return false;
	if (target_number > uint64_Numberofevents) return false;
	__int32 eventsize;
	if (data_format_in_userheader == 2 ) eventsize = 2 * Numberofcoordinates;
	if (data_format_in_userheader == 5 ) eventsize = 8 * Numberofcoordinates;
	if (data_format_in_userheader == 10) eventsize = 4 * Numberofcoordinates;

	if (DAQ_ID == DAQ_ID_RAW32BIT) eventsize = 4 * number_of_channels * (max_number_of_hits+1);

	if (input_lmf->filesize < (unsigned __int64)(eventsize)*target_number + (unsigned  __int64)(Headersize + User_header_size)) return false;

	unsigned __int64 new_position = (unsigned __int64)(eventsize)*target_number + (unsigned __int64)(Headersize + User_header_size);

	input_lmf->seek(new_position);



	uint64_number_of_read_events = target_number;
	must_read_first = true;
	errorflag = 0;
	input_lmf->error = 0;
	return true;
}







///////////////////////////////////////////////////////////////////////////////////
__int32 LMF_IO::PCIGetTDC_TDC8HP_25psGroupMode(unsigned __int64 &ref_ui64TDC8HPAbsoluteTimeStamp, __int32 count, unsigned __int32 * Buffer)
///////////////////////////////////////////////////////////////////////////////////
{
	memset(number_of_hits,0,num_channels*sizeof(__int32));		// clear the hit-counts values in _TDC array

	unsigned __int32 ui32DataWord;
	bool bOKFlag = false;
	unsigned __int8 ucTDCChannel;

	for(__int32 i = 0; i < count ; ++i)
	{
		ui32DataWord = Buffer[i];
		if ((ui32DataWord & 0xf8000000) == 0x18000000) // handle output level info
		{
			__int32 n = ui32DataWord & 0x7e00000;
			n >>= 21;
			unsigned __int64 ui64_temp_LevelInfo = ui32DataWord & 0x1fffff;
			if (n > 20) continue;
			if (n < 9) ui64_temp_LevelInfo >>= (9-n); else ui64_temp_LevelInfo <<= (n-9);
			
			n-=9;
			if (n<0) n = 0;

			unsigned __int64 ui64_tempL_LevelInfo = ui64LevelInfo >> (n+21);
			ui64_tempL_LevelInfo <<= (n+21);
			unsigned __int64 ui64_tempR_LevelInfo = n!=0 ? ui64LevelInfo << (64-n) : 0;
			ui64_tempR_LevelInfo = n!=0 ? ui64_tempR_LevelInfo >> (64-n) : 0;

			//unsigned __int64 old = ui64LevelInfo;
			ui64LevelInfo = ui64_tempL_LevelInfo | ui64_tempR_LevelInfo | ui64_temp_LevelInfo;

			continue;
		}
		if( (ui32DataWord&0xC0000000)>0x40000000)		// valid data only if rising or falling trigger indicated
		{
			__int32 lTDCData = (ui32DataWord&0x00FFFFFF);
			if(lTDCData & 0x00800000)				// detect 24 bit signed flag
				lTDCData |= 0xff000000;				// if detected extend negative value to 32 bit
			if(!this->TDC8HP.VHR_25ps) 				// correct for 100ps if necessary
				lTDCData >>= 2;
			
			ucTDCChannel = (unsigned __int8)((ui32DataWord&0x3F000000)>>24);		// extract channel information
			// calculate TDC channel to _TDC channel
			bool valid = false;
			if((ucTDCChannel > 41) && (ucTDCChannel < 51)) {
				ucTDCChannel -= 25;
				valid = true;
			}
			else if((ucTDCChannel > 20) && (ucTDCChannel < 30)) {
				ucTDCChannel -= 12;
				valid = true;
			} else if ((ucTDCChannel >= 0) && (ucTDCChannel < 9)) {
				valid = true;
			}
			
			if (valid) {
			bool bIsFalling = true;
			if ((ui32DataWord&0xC0000000) == 0xC0000000) bIsFalling = false;

			if (!bIsFalling) {
				ucTDCChannel += TDC8HP.channel_offset_for_rising_transitions;
			}

			if(ucTDCChannel < num_channels)	// if detected channel fits into TDC array then sort
			{
				++number_of_hits[ucTDCChannel];
				__int32 cnt = number_of_hits[ucTDCChannel];
				// increase Hit Counter;
				
				// test for oversized Hits
				if(cnt > num_ions) {
					--number_of_hits[ucTDCChannel];
					--cnt;
				}
				else			
					// if Hit # ok then store it
					i32TDC[ucTDCChannel*num_ions+cnt-1] = lTDCData;

				bOKFlag = true;
			}
		} 
		} 
		else
		{
			if ((ui32DataWord & 0xf0000000) == 0x00000000) {			// GroupWord detected
				this->TDC8HP.ui32AbsoluteTimeStamp = ui32DataWord & 0x00ffffff;
		}
			else if ((ui32DataWord & 0x10000000) == 0x10000000) {			// RollOverWord detected ?
				unsigned __int32 ui32newRollOver = (ui32DataWord & 0x00ffffff);
				if (ui32newRollOver > this->TDC8HP.ui32oldRollOver) {
					this->TDC8HP.ui64RollOvers += ui32newRollOver - this->TDC8HP.ui32oldRollOver;
				} else if (ui32newRollOver < this->TDC8HP.ui32oldRollOver) {
					this->TDC8HP.ui64RollOvers += ui32newRollOver;
					this->TDC8HP.ui64RollOvers += 1;
					this->TDC8HP.ui64RollOvers += (unsigned __int32)(0x00ffffff) - this->TDC8HP.ui32oldRollOver;
				}
				this->TDC8HP.ui32oldRollOver = ui32newRollOver;
			}
			//	only for debugging:
#ifdef _DEBUG
			else if (((ui32DataWord & 0xc0000000)>>30) == 0x00000001)			// ErrorWord detected ?
			{
				__int32 channel = (ui32DataWord & 0x3f000000)>>24;
				__int32 error = (ui32DataWord   & 0x00ff0000)>>16;
				__int32 count = ui32DataWord    & 0x0000ffff;
			}
#endif

		}
	}

	if (bOKFlag)
	{
		ref_ui64TDC8HPAbsoluteTimeStamp  = this->TDC8HP.ui64RollOvers * (unsigned __int64)(0x0000000001000000);
		ref_ui64TDC8HPAbsoluteTimeStamp += (unsigned __int64)(this->TDC8HP.ui32AbsoluteTimeStamp);
		this->TDC8HP.ui64TDC8HP_AbsoluteTimeStamp = ref_ui64TDC8HPAbsoluteTimeStamp;
	}
	
	return bOKFlag;
}







/////////////////////////////////////////////////////////////////
bool LMF_IO::Read_TDC8HP_raw_format(unsigned __int64 &ui64TDC8HP_AbsoluteTimeStamp_)
/////////////////////////////////////////////////////////////////
{
	__int32 count;
	*input_lmf >> count;
	if (input_lmf->error) return false;
	if (!count) return false;
	if (ui32buffer_size < count) {
		if (ui32buffer) {delete[] ui32buffer; ui32buffer = 0;}
		ui32buffer_size = count + 5000;
		ui32buffer = new unsigned __int32[ui32buffer_size];
	}
	input_lmf->read(this->ui32buffer,count*sizeof(__int32));
	if (input_lmf->error) return false;
	if (!PCIGetTDC_TDC8HP_25psGroupMode(ui64TDC8HP_AbsoluteTimeStamp_, count, this->ui32buffer)) return false;
	return true;
}












/////////////////////////////////////////////////////////////////
bool LMF_IO::ReadNextEvent()
/////////////////////////////////////////////////////////////////
{
	unsigned __int32	i,j;
	__int32				int_Dummy;
	unsigned __int16	unsigned_short_Dummy;
	double			double_Dummy;

	if (!input_lmf) {
		errorflag = 9;
		return false;
	}
	if (DAQ_ID != DAQ_ID_SIMPLE) {
		if ((max_number_of_hits == 0 || number_of_channels == 0) && (DAQ_ID != DAQ_ID_CAMAC)) {
			errorflag = 14;
			return false;
		}
		if ((data_format_in_userheader == LM_CAMAC) && (DAQ_ID == DAQ_ID_CAMAC)) {
			return ReadNextCAMACEvent();
		}
	}

	if (input_lmf->error) {
		if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1;
		return false;
	}

	if (data_format_in_userheader ==  2) memset(us16TDC,0,num_channels*num_ions*2);
	if (data_format_in_userheader ==  5) memset(dTDC,0,num_channels*num_ions*8);
	if (data_format_in_userheader == 10) memset(i32TDC,0,num_channels*num_ions*4);
	if (TDC8HP.variable_event_length == 1) memset(i32TDC,0,num_channels*num_ions*4);
	if (TDC8PCI2.variable_event_length == 1) memset(us16TDC,0,num_channels*num_ions*2);

	unsigned __int64 HPTDC_event_length = 0;
	unsigned __int64 TDC8PCI2_event_length = 0;

	DOUBLE_timestamp = 0.;
	ui64_timestamp = 0;

	if (TDC8HP.variable_event_length == 1) {
		if (this->TDC8HP.UserHeaderVersion >= 5 && this->TDC8HP.GroupingEnable_p66) {
			
			while (!Read_TDC8HP_raw_format(ui64_timestamp)) {
				if (input_lmf->error) break;
			}
			if (this->LMF_Version >= 11) {
				changed_mask_read = 0;
				*input_lmf >> changed_mask_read;
				unsigned __int32 changed_mask = changed_mask_read;
				if (changed_mask) {
					for (__int32 i=0;i<32;i++) {
						if (changed_mask & 0x1) {
							double double_temp;
							*input_lmf >> double_temp;
							Parameter[901+i] = double_temp;
						}
						changed_mask >>= 1;
					}
				}
			}
			number_of_bytes_in_PostEventData = 0;
			if (this->TDC8HP.UserHeaderVersion >= 7) {
				*input_lmf >> number_of_bytes_in_PostEventData;
				if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
					input_lmf->error = 21;
					return false;
				}
				for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
					unsigned __int8 byte_dummy;
					*input_lmf>>byte_dummy;
					ui8_PostEventData[i] = byte_dummy;
				}
			}

			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
			++uint64_number_of_read_events;
			DOUBLE_timestamp = double(ui64_timestamp)/frequency;  // time stamp in seconds.
			must_read_first = false;
			return true;
		} 
		if (!TDC8HP.exotic_file_type) {
			*input_lmf >> HPTDC_event_length;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}

#ifdef LINUX
			if((HPTDC_event_length & 0xff00000000000000LL) != 0xff00000000000000LL) {
#else
			if((HPTDC_event_length & 0xff00000000000000) != 0xff00000000000000) {
#endif
				this->errorflag = 2;
				return false;
			}

#ifdef LINUX
			HPTDC_event_length = HPTDC_event_length & 0x00ffffffffffffffLL;
#else
			HPTDC_event_length = HPTDC_event_length & 0x00ffffffffffffff;
#endif

			*input_lmf >> uint64_LMF_EventCounter;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
		} else {
			__int32 i32HPTDC_event_length;
			*input_lmf >> i32HPTDC_event_length;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
		}
	}

	if (TDC8PCI2.variable_event_length == 1) {
			*input_lmf >> TDC8PCI2_event_length;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}

#ifdef LINUX
			if((TDC8PCI2_event_length & 0xff00000000000000LL) != 0xff00000000000000LL) {
#else
			if((TDC8PCI2_event_length & 0xff00000000000000) != 0xff00000000000000) {
#endif
				this->errorflag = 2;
				return false;
			}

#ifdef LINUX
			TDC8PCI2_event_length = TDC8PCI2_event_length & 0x00ffffffffffffffLL;
#else
			TDC8PCI2_event_length = TDC8PCI2_event_length & 0x00ffffffffffffff;
#endif

			*input_lmf >> uint64_LMF_EventCounter;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
	}


	if (DAQ_ID == DAQ_ID_HM1 && DAQVersion >= 20080507) {
			unsigned __int64 event_length;
			*input_lmf >> event_length;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}

#ifdef LINUX
			if((event_length & 0xff00000000000000LL) != 0xff00000000000000LL) {
				this->errorflag = 2;
				return false;
			}
			event_length = event_length & 0x00ffffffffffffffLL;
#else
			if((event_length & 0xff00000000000000) != 0xff00000000000000) {
				this->errorflag = 2;
				return false;
			}
			event_length = event_length & 0x00ffffffffffffff;
#endif

			*input_lmf >> uint64_LMF_EventCounter;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
	}

	++uint64_number_of_read_events;

	//-------------------------------
	//  Read Time Stamp
	
	if (timestamp_format > 0) {

		if (timestamp_format > 0 ) {
			myLARGE_INTEGER aaa;
			aaa.QuadPart = 0;
			*input_lmf >> aaa.LowPart;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
			if (timestamp_format == 2 ) *input_lmf >> aaa.HighPart;
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 1; return false;}
			ui64_timestamp = aaa.QuadPart;
		}
		DOUBLE_timestamp = double(ui64_timestamp)/frequency;  // time stamp in seconds.
	}


	if (DAQ_ID != DAQ_ID_SIMPLE && !TDC8HP.variable_event_length && !TDC8PCI2.variable_event_length) {
		for (i=0;i<number_of_channels+number_of_channels2;++i) {
				if (DAQ_ID == DAQ_ID_HM1_ABM) {
					number_of_hits[i] = 1;
					if (data_format_in_userheader ==  2) *input_lmf >> us16TDC[i*num_ions];
					if (data_format_in_userheader ==  5) *input_lmf >> dTDC[i*num_ions];
					if (data_format_in_userheader == 10) *input_lmf >> i32TDC[i*num_ions];
					if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
				}
				if (DAQ_ID != DAQ_ID_HM1_ABM) {
					if (data_format_in_userheader ==  2) {
						*input_lmf >> unsigned_short_Dummy;
						if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
						if (DAQ_ID == DAQ_ID_HM1) unsigned_short_Dummy = (unsigned_short_Dummy & 0x0007) - 1;
						number_of_hits[i] = (__int32 )unsigned_short_Dummy;
						for (j=0;j<max_number_of_hits;++j)  *input_lmf >> us16TDC[i*num_ions+j];
					}
					if (data_format_in_userheader ==  5) {
						*input_lmf >> double_Dummy;
						if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
						number_of_hits[i] =__int32(double_Dummy+0.1);
						for (j=0;j<max_number_of_hits;++j)  *input_lmf >> dTDC[i*num_ions+j];
					}
					if (data_format_in_userheader == 10) {
						*input_lmf >> int_Dummy;
						if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
						number_of_hits[i] = (__int32 )int_Dummy;
						for (j=0;j<max_number_of_hits;++j)  *input_lmf >> i32TDC[i*num_ions+j];
					}
					if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
				}
		} // for i
	}

	if (DAQ_ID == DAQ_ID_HM1 && data_format_in_userheader == -1) {
		for (unsigned __int32 channel=0;channel<number_of_channels+number_of_channels2;++channel) {
			unsigned  __int16 nn;
			*input_lmf >> nn;		// store hits for this channel

			__int32 n = nn;
			n = (n & 0x07)-1;		// now hits only
			n = n < 0 ? 0 : n;		// avoid "negative" hits
			if (n > __int32(max_number_of_hits)) n = max_number_of_hits;

			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
			number_of_hits[channel] = n;
			for(__int32 i=0;i<(__int32 )n;++i) {	// transfer selected hits
				unsigned __int16 us16data;
				*input_lmf >> us16data;
				us16TDC[channel*num_ions+i] = us16data;
			}

			if (DAQVersion > 20080507) {
				__int16 level_info;
				*input_lmf >> level_info;				// Data for LevelInfo compatibility
				*input_lmf >> number_of_bytes_in_PostEventData;
				if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
					input_lmf->error = 21;
					return false;
				}
				for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
					unsigned __int8 byte_dummy;
					*input_lmf>>byte_dummy;
					ui8_PostEventData[i] = byte_dummy;
				}
			}

			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
		}
	}

	if ((DAQ_ID == DAQ_ID_TDC8 || DAQ_ID == DAQ_ID_2TDC8) && TDC8PCI2.variable_event_length == 1) {
		for (unsigned __int32 channel=0;channel<number_of_channels+number_of_channels2;++channel) {
			unsigned  __int16 n;
			*input_lmf >> n;		// store hits for this channel

			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
			number_of_hits[channel] = n;
			for(__int32 i=0;i<(__int32 )n;++i) {	// transfer selected hits
				unsigned __int16 us16data;
				*input_lmf >> us16data;
				us16TDC[channel*num_ions+i] = us16data;
			}
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
		}
	}

	if ((DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) && TDC8HP.variable_event_length == 1) {
		for (unsigned __int32 channel=0;channel<number_of_channels;++channel) {
			unsigned  __int16 n16;
			unsigned  __int32 n32;
			__int32 n;
			if (TDC8HP.exotic_file_type) {
				*input_lmf >> n32;
				n = n32;
			} else {
				*input_lmf >> n16;		// store hits for this channel
				n = n16;
			}
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
			number_of_hits[channel] = n;
			for(__int32 i=0;i<(__int32 )n;++i) {	// transfer selected hits
				__int32 i32data;
				*input_lmf >> i32data;
				i32TDC[channel*num_ions+i] = i32data;
			}
			if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
		}
		if (DAQ_ID == DAQ_ID_TDC8HP && this->LMF_Version >= 9) {
			*input_lmf >> ui64LevelInfo;
		}
		if (this->LMF_Version >= 11) {
			changed_mask_read = 0;
			*input_lmf >> changed_mask_read;
			unsigned __int32 changed_mask = changed_mask_read;
			if (changed_mask) {
				for (__int32 i=0;i<32;i++) {
					if (changed_mask & 0x1) {
						double double_temp;
						*input_lmf >> double_temp;
						Parameter[901+i] = double_temp;
					}
					changed_mask >>= 1;
				}
			}
		}
		number_of_bytes_in_PostEventData = 0;
		if (DAQ_ID == DAQ_ID_TDC8HP && this->TDC8HP.UserHeaderVersion >= 7) {
			*input_lmf >> number_of_bytes_in_PostEventData;
			if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
				input_lmf->error = 21;
				return false;
			}
			for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
				unsigned __int8 byte_dummy;
				*input_lmf>>byte_dummy;
				ui8_PostEventData[i] = byte_dummy;
			}
		}
	}

	if (DAQ_ID == DAQ_ID_SIMPLE) {
		unsigned __int16 us16_Dummy;
		__int32 i32_Dummy;
		__int32 number_of_words;
		__int32 channel;
		number_of_words = 0;
		if (data_format_in_userheader == 2 ) {*input_lmf >> us16_Dummy; number_of_words = us16_Dummy;}
		if (data_format_in_userheader == 10) *input_lmf >> number_of_words;
		if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
		for (i=0;i<number_of_channels;++i) number_of_hits[i] = 0;

		bool read_channel_marker;
		read_channel_marker = true;
		if (data_format_in_userheader == 2) {
			while (number_of_words > 0) {
				number_of_words--;
				*input_lmf >> us16_Dummy;
				if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
				if (read_channel_marker) {
					read_channel_marker = false;
					channel = 0;
					channel =__int32((us16_Dummy & 0xff00)  >> 8);
					number_of_hits[channel] =__int32(us16_Dummy & 0x00ff);
					i=0;
				} else {
					us16TDC[channel*num_ions+i] = us16_Dummy;
					++i;
					if (i == number_of_hits[channel]) read_channel_marker = true;
				}
			}
		}
		if (data_format_in_userheader == 10) {
			while (number_of_words > 0) {
				number_of_words--;
				*input_lmf >> i32_Dummy;
				if (input_lmf->error) {if (input_lmf->eof) this->errorflag = 18; else this->errorflag = 2; return false;}
				if (read_channel_marker) {
					read_channel_marker = false;
					channel = 0;
					channel =__int32((i32_Dummy & 0xff000000)  >> 24);
					number_of_hits[channel] =__int32(i32_Dummy & 0x000000ff);
					i=0;
				} else {
					i32TDC[channel*num_ions+i] = i32_Dummy;
					++i;
					if (i == number_of_hits[channel]) read_channel_marker = true;
				}
			}
		}
	}
	
	if (LMF_Version >= 9 && DAQ_ID == DAQ_ID_2TDC8 && DAQVersion >= 20110208) {
		number_of_bytes_in_PostEventData = 0;
		*input_lmf >> number_of_bytes_in_PostEventData;
		if (number_of_bytes_in_PostEventData > MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA) {
			input_lmf->error = 21;
			return false;
		}
		for (__int32 i=0;i<number_of_bytes_in_PostEventData;i++) {
			unsigned __int8 byte_dummy;
			*input_lmf>>byte_dummy;
			ui8_PostEventData[i] = byte_dummy;
		}
	}

	must_read_first = false;
	return true;
}








/////////////////////////////////////////////////////////////////
void LMF_IO::WriteCAMACArray(double timestamp, unsigned __int32 data[])
/////////////////////////////////////////////////////////////////
{
	unsigned __int16 dummy_uint16;
	unsigned __int8  dummy_uint8;


	if (!output_lmf || !OutputFileIsOpen) {
		errorflag = 10;
		return;
	}
	
	++uint64_number_of_written_events;

	myLARGE_INTEGER LARGE_Timestamp;
	LARGE_Timestamp.QuadPart = (__int64)(timestamp * frequency);

	if (timestamp_format_output >= 1) {
		if (data_format_in_userheader_output == 2) {
			dummy_uint16 = (unsigned __int16)(LARGE_Timestamp.LowPart & 0x0000ffff); *output_lmf << dummy_uint16;			// 32 Bit Low part, lower  16 Bit
			dummy_uint16 = (unsigned __int16)(LARGE_Timestamp.LowPart & 0xffff0000); *output_lmf << dummy_uint16;			// 32 Bit Low part, higher 16 Bit
		}
		if (data_format_in_userheader_output == 6) {
			dummy_uint8 = (unsigned __int8)(LARGE_Timestamp.LowPart & 0x000000ff); *output_lmf << dummy_uint8;			// 32 Bit Low part, lower 16 Bit, lower 8 bit
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.LowPart >> 8) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit Low part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)0; *output_lmf << dummy_uint8;	
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.LowPart >> 16) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit Low part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.LowPart >> 24) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit Low part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)0; *output_lmf << dummy_uint8;
		}
	}
	if (timestamp_format_output == 2) {
		if (data_format_in_userheader_output == 2) {
			dummy_uint16 = (unsigned __int16)(LARGE_Timestamp.HighPart & 0x0000ffff); *output_lmf << dummy_uint16;			// 32 Bit High part, lower  16 Bit
			dummy_uint16 = (unsigned __int16)(LARGE_Timestamp.HighPart & 0xffff0000); *output_lmf << dummy_uint16;			// 32 Bit High part, higher 16 Bit
		}
		if (data_format_in_userheader_output == 6) {
			dummy_uint8 = (unsigned __int8)(LARGE_Timestamp.HighPart & 0x000000ff); *output_lmf << dummy_uint8;			// 32 Bit High part, lower 16 Bit, lower 8 bit
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.HighPart >> 8) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit High part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)0; *output_lmf << dummy_uint8;	
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.HighPart >> 16) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit High part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)((LARGE_Timestamp.HighPart >> 24) & 0x000000ff); *output_lmf << dummy_uint8;	// 32 Bit High part, lower 16 Bit, high 8 bit
			dummy_uint8 = (unsigned __int8)0; *output_lmf << dummy_uint8;
		}
	}

	__int32 i;

	for (i=0;i<Numberofcoordinates - timestamp_format_output * 2;++i) {
		if (data_format_in_userheader_output == 2) {
			dummy_uint16 = (unsigned __int16)( data[i]        & 0x0000ffff); *output_lmf << dummy_uint16;
		}
		if (data_format_in_userheader_output == 6) {
			dummy_uint8 = (unsigned __int8)( data[i]        & 0x000000ff); *output_lmf << dummy_uint8;
			dummy_uint8 = (unsigned __int8)((data[i] >>  8) & 0x000000ff); *output_lmf << dummy_uint8;
			dummy_uint8 = (unsigned __int8)((data[i] >> 16) & 0x000000ff); *output_lmf << dummy_uint8;
		}
	}

	return;
}






/////////////////////////////////////////////////////////////////
bool LMF_IO::ReadNextCAMACEvent()
/////////////////////////////////////////////////////////////////
{
	__int32 i;

	if (!input_lmf) {
		errorflag = 9;
		return false;
	}

	++uint64_number_of_read_events;

	//-------------------------------
	//  Read Time Stamp
	DOUBLE_timestamp = 0.;
	
	unsigned __int32 time_temp;
	unsigned __int8 byte_1,byte_2,byte_3;
	unsigned __int16 unsigned_short;

	if (timestamp_format > 0) {
		//TRY
			myLARGE_INTEGER LARGE_timestamp;
			LARGE_timestamp.QuadPart = 0;
			if (timestamp_format >= 1 )
			{
				if (data_format_in_userheader == 2) {
					*input_lmf >> unsigned_short;
					time_temp = unsigned_short;
					*input_lmf >> unsigned_short;
					LARGE_timestamp.LowPart = time_temp + unsigned_short*256*256;
				}
				if (data_format_in_userheader == 6) {
					*input_lmf >> byte_1; *input_lmf >> byte_2; *input_lmf >> byte_3;
					time_temp = byte_1 + byte_2 * 256 + byte_3 * 256 * 256;
					LARGE_timestamp.LowPart = time_temp;

					*input_lmf >> byte_1; *input_lmf >> byte_2; *input_lmf >> byte_3;
					time_temp = byte_1 + byte_2 * 256 + byte_3 * 256 * 256;
					LARGE_timestamp.LowPart += time_temp * 256*256;
				}
			}
			if (timestamp_format == 2 ) {
				if (data_format_in_userheader == 2) {
					*input_lmf >> unsigned_short;
					time_temp = unsigned_short;
					*input_lmf >> unsigned_short;
					LARGE_timestamp.HighPart = time_temp + unsigned_short*256*256;
				}
				if (data_format_in_userheader == 6) {
					*input_lmf >> byte_1; *input_lmf >> byte_2; *input_lmf >> byte_3;
					time_temp = byte_1 + byte_2 * 256 + byte_3 * 256 * 256;
					LARGE_timestamp.HighPart = time_temp;

					*input_lmf >> byte_1; *input_lmf >> byte_2; *input_lmf >> byte_3;
					time_temp = byte_1 + byte_2 * 256 + byte_3 * 256 * 256;
					LARGE_timestamp.HighPart += time_temp * 256*256;
				}
			}
			ui64_timestamp = LARGE_timestamp.QuadPart;
/*		CATCH(CArchiveException,e)
			errorflag = 1;	// error reading timestamp
			return false;
		END_CATCH
		*/
		DOUBLE_timestamp = double(ui64_timestamp)/frequency;  // time stamp in seconds.
	}

//	TRY
		for (i=0;i<Numberofcoordinates - timestamp_format * 2;++i) {
			if (data_format_in_userheader == 6) {
				*input_lmf >> byte_1; *input_lmf >> byte_2; *input_lmf >> byte_3;
				CAMAC_Data[i] = byte_1 + byte_2 * 256 + byte_3 * 256 * 256;
			}
			if (data_format_in_userheader == 2) {
				*input_lmf >> unsigned_short;
				CAMAC_Data[i] = unsigned_short;
			}
		} // for i
/*	CATCH(CArchiveException,e)
		errorflag = 2; // error reading data
		return false;
	END_CATCH
	*/

	must_read_first = false;
	return true;
}





/////////////////////////////////////////////
void LMF_IO::GetCAMACArray(unsigned __int32 data[])
/////////////////////////////////////////////
{
	if (must_read_first) {
		if (!ReadNextCAMACEvent()) return;
	}
	for (__int32 i=0;i<Numberofcoordinates - timestamp_format * 2;++i) data[i] = CAMAC_Data[i];
}





/////////////////////////////////////////////////////////////////
void LMF_IO::GetTDCDataArray(__int32 *tdc)
/////////////////////////////////////////////////////////////////
{
	__int32 i,j;
	__int32 ii;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	//__int32 max_channel = (number_of_channels+number_of_channels2 < num_channels) ? (number_of_channels+number_of_channels2) : num_channels;
	//__int32 max_hits = (max_number_of_hits < num_ions) ? max_number_of_hits : num_ions;
	__int32 max_channel = num_channels;
	__int32 max_hits = num_ions;

	if (data_format_in_userheader == LM_USERDEF) {
		if (DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = i32TDC[i*num_ions+j];
			}
		}
		if (DAQ_ID == DAQ_ID_TDC8 || DAQ_ID == DAQ_ID_2TDC8 || DAQ_ID == DAQ_ID_HM1) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = us16TDC[i*num_ions+j];
			}
		}
	}
	if (data_format_in_userheader == 10) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = i32TDC[i*num_ions+j];
		}
	}
	if (data_format_in_userheader == 2) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] =__int32(us16TDC[i*num_ions+j]);
		}
	}
	if (data_format_in_userheader == 5) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) {
				if (dTDC[i*num_ions+j] >= 0.) ii =__int32(dTDC[i*num_ions+j]+1.e-19);
				if (dTDC[i*num_ions+j] <  0.) ii =__int32(dTDC[i*num_ions+j]-1.e-19);
				tdc[i*num_ions+j] = ii;
			}
		}
	}
}




/////////////////////////////////////////////////////////////////
void LMF_IO::GetTDCDataArray(__int64 *tdc)
/////////////////////////////////////////////////////////////////
{
	__int32 i,j;
	__int32 ii;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	//__int32 max_channel = (number_of_channels+number_of_channels2 < num_channels) ? (number_of_channels+number_of_channels2) : num_channels;
	//__int32 max_hits = (max_number_of_hits < num_ions) ? max_number_of_hits : num_ions;
	__int32 max_channel = num_channels;
	__int32 max_hits = num_ions;

	if (data_format_in_userheader == LM_USERDEF) {
		if (DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = i32TDC[i*num_ions+j];
			}
		}
		if (DAQ_ID == DAQ_ID_TDC8 || DAQ_ID == DAQ_ID_2TDC8 || DAQ_ID == DAQ_ID_HM1) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = us16TDC[i*num_ions+j];
			}
		}
	}
	if (data_format_in_userheader == 10) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = i32TDC[i*num_ions+j];
		}
	}
	if (data_format_in_userheader == 2) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] =__int32(us16TDC[i*num_ions+j]);
		}
	}
	if (data_format_in_userheader == 5) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) {
				if (dTDC[i*num_ions+j] >= 0.) ii =__int32(dTDC[i*num_ions+j]+1.e-19);
				if (dTDC[i*num_ions+j] <  0.) ii =__int32(dTDC[i*num_ions+j]-1.e-19);
				tdc[i*num_ions+j] = ii;
			}
		}
	}
}




/////////////////////////////////////////////////////////////////
void LMF_IO::GetTDCDataArray(unsigned __int16 *tdc)
/////////////////////////////////////////////////////////////////
{
	__int32 i,j;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	//__int32 max_channel = (number_of_channels+number_of_channels2 < num_channels) ? (number_of_channels+number_of_channels2) : num_channels;
	//__int32 max_hits = (max_number_of_hits < num_ions) ? max_number_of_hits : num_ions;
	__int32 max_channel = num_channels;
	__int32 max_hits = num_ions;

	if (data_format_in_userheader == LM_USERDEF) {
		if (DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = (unsigned __int16)(i32TDC[i*num_ions+j]);
			}
		}
		if (DAQ_ID == DAQ_ID_TDC8 || DAQ_ID == DAQ_ID_2TDC8 || DAQ_ID == DAQ_ID_HM1) {
			for (i=0;i<max_channel;++i) {
				for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = (unsigned __int16)(us16TDC[i*num_ions+j]);
			}
		}
	}
	if (data_format_in_userheader == 10) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = (unsigned __int16)(i32TDC[i*num_ions+j]);
		}
	}
	if (data_format_in_userheader == 2) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = us16TDC[i*num_ions+j];
		}
	}
	if (data_format_in_userheader == 5) {
		for (i=0;i<max_channel;++i) {
			for (j=0;j<max_hits;++j) tdc[i*num_ions+j] = (unsigned __int16)(dTDC[i*num_ions+j]+1e-7);
		}
	}
}



/////////////////////////////////////////////////////////////////
void LMF_IO::GetTDCDataArray(double *tdc)
/////////////////////////////////////////////////////////////////
{
	__int32 i,j;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	__int32 max_channel = (__int32(number_of_channels+number_of_channels2) < num_channels) ? __int32(number_of_channels+number_of_channels2) : num_channels;
	__int32 max_hits = (__int32(max_number_of_hits) < num_ions) ? __int32(max_number_of_hits) : num_ions;
	//__int32 max_channel = num_channels;
	//__int32 max_hits = num_ions;

	if (data_format_in_userheader == LM_USERDEF) {
		if (DAQ_ID == DAQ_ID_TDC8HP || DAQ_ID == DAQ_ID_TDC8HPRAW) {
			for (i=0;i<max_channel;++i) {
				__int32 n = number_of_hits[i];
				for (j=0;j<n;++j) tdc[i*num_ions+j] = double(i32TDC[i*num_ions+j]);
			}
		}
		if (DAQ_ID == DAQ_ID_TDC8 || DAQ_ID == DAQ_ID_2TDC8 || DAQ_ID == DAQ_ID_HM1) {
			for (i=0;i<max_channel;++i) {
				__int32 n = number_of_hits[i];
				for (j=0;j<n;++j) tdc[i*num_ions+j] = double(us16TDC[i*num_ions+j]);
			}
		}
	}
	if (data_format_in_userheader == 10) {
		for (i=0;i<max_channel;++i) {
			__int32 n = number_of_hits[i];
			for (j=0;j<n;++j) tdc[i*num_ions+j] = double(i32TDC[i*num_ions+j]);
		}
	}
	if (data_format_in_userheader == 2) {
		for (i=0;i<max_channel;++i) {
			__int32 n = number_of_hits[i];
			for (j=0;j<n;++j) tdc[i*num_ions+j] = double(us16TDC[i*num_ions+j]);
		}
	}
	if (data_format_in_userheader == 5) {
		for (i=0;i<max_channel;++i) {
			__int32 n = number_of_hits[i];
			for (j=0;j<n;++j) tdc[i*num_ions+j] = dTDC[i*num_ions+j];
		}
	}
}





/////////////////////////////////////////////////////////////////
void LMF_IO::GetNumberOfHitsArray(unsigned __int32 cnt[]) {
/////////////////////////////////////////////////////////////////
	__int32 i;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	for (i=0;i<num_channels;++i) cnt[i] = (number_of_hits[i] < (unsigned __int32)(num_ions)) ? number_of_hits[i] : num_ions;
}




/////////////////////////////////////////////////////////////////
void LMF_IO::GetNumberOfHitsArray(__int32 cnt[]) {
/////////////////////////////////////////////////////////////////
	__int32 i;

	if (must_read_first) {
		if (!ReadNextEvent()) return;
	}

	for (i=0;i<num_channels;++i)cnt[i] = (number_of_hits[i] < (unsigned __int32)(num_ions)) ? number_of_hits[i] : num_ions;
}




/////////////////////////////////////////////////////////////////
const char* LMF_IO::GetErrorText(__int32 error_id)
/////////////////////////////////////////////////////////////////
{
	return error_text[error_id];
}


/////////////////////////////////////////////////////////////////
void LMF_IO::GetErrorText(__int32 error_code, __int8 text[])
/////////////////////////////////////////////////////////////////
{
	sprintf(text,"%s",error_text[error_code]);
	return;
}


/////////////////////////////////////////////////////////////////
void LMF_IO::GetErrorText(__int8 text[])
/////////////////////////////////////////////////////////////////
{
	GetErrorText(errorflag,text);
	return;
}


/////////////////////////////////////////////////////////////////
unsigned __int64 LMF_IO::GetEventNumber()
/////////////////////////////////////////////////////////////////
{
	return uint64_number_of_read_events;
}

/////////////////////////////////////////////////////////////////
unsigned __int32 LMF_IO::GetNumberOfChannels()
/////////////////////////////////////////////////////////////////
{
	return number_of_channels+number_of_channels2;
}

/////////////////////////////////////////////////////////////////
unsigned __int32 LMF_IO::GetMaxNumberOfHits()
/////////////////////////////////////////////////////////////////
{
	return max_number_of_hits;
}


/////////////////////////////////////////////////////////////////
double LMF_IO::GetDoubleTimeStamp()
/////////////////////////////////////////////////////////////////
{
	if (must_read_first) {
		if (!ReadNextEvent()) return 0.;
	}
	return DOUBLE_timestamp;
}

/////////////////////////////////////////////////////////////////
unsigned __int64 LMF_IO::Getuint64TimeStamp()
/////////////////////////////////////////////////////////////////
{
	if (must_read_first) {
		if (!ReadNextEvent()) return ui64_timestamp;
	}
	return ui64_timestamp;
}




/////////////////////////////////////////////////////////////////
bool LMF_IO::Clone(LMF_IO * clone)
/////////////////////////////////////////////////////////////////
{
	if (!clone) return 0;

	clone->Versionstring			= this->Versionstring;
	clone->FilePathName				= this->FilePathName;
	clone->OutputFilePathName		= this->OutputFilePathName;
	clone->Comment					= this->Comment;
	clone->Comment_output			= this->Comment_output;
	clone->DAQ_info					= this->DAQ_info;
	clone->Camac_CIF				= this->Camac_CIF;
	
	clone->iLMFcompression			= this->iLMFcompression;

	clone->Starttime				= this->Starttime;
	clone->Stoptime					= this->Stoptime;
	clone->Starttime_output			= this->Starttime_output;
	clone->Stoptime_output			= this->Stoptime_output;
	
	clone->time_reference			= this->time_reference;
	clone->time_reference_output	= this->time_reference_output;
	
	clone->ArchiveFlag				= this->ArchiveFlag;
	clone->Cobold_Header_version	= this->Cobold_Header_version;
	clone->Cobold_Header_version_output	= this->Cobold_Header_version_output;
	
	clone->uint64_LMF_EventCounter	= this->uint64_LMF_EventCounter;
	clone->uint64_number_of_read_events	= this->uint64_number_of_read_events;
	clone->uint64_Numberofevents	= this->uint64_Numberofevents;
	
	clone->Numberofcoordinates		= this->Numberofcoordinates;
	clone->CTime_version			= this->CTime_version;
	clone->CTime_version_output		= this->CTime_version_output;
	clone->CTime_version_output		= this->CTime_version_output;
	clone->SIMPLE_DAQ_ID_Orignial	= this->SIMPLE_DAQ_ID_Orignial;
	clone->DAQVersion				= this->DAQVersion;
	clone->DAQVersion_output		= this->DAQVersion_output;
	clone->DAQ_ID					= this->DAQ_ID;
	clone->DAQ_ID_output			= this->DAQ_ID_output;
	clone->data_format_in_userheader	= this->data_format_in_userheader;
	clone->data_format_in_userheader_output	= this->data_format_in_userheader_output;
	
	clone->Headersize				= this->Headersize;
	clone->User_header_size			= this->User_header_size;
	clone->User_header_size_output	= this->User_header_size_output;
	
	clone->IOaddress				= this->IOaddress;
	clone->timestamp_format			= this->timestamp_format;
	clone->timestamp_format_output	= this->timestamp_format_output;
	clone->timerange				= this->timerange;
	
	clone->number_of_channels		= this->number_of_channels;
	clone->number_of_channels2		= this->number_of_channels2;
	clone->max_number_of_hits		= this->max_number_of_hits;
	clone->max_number_of_hits2		= this->max_number_of_hits2;
	
	clone->number_of_channels_output	= this->number_of_channels_output;
	clone->number_of_channels2_output	= this->number_of_channels2_output;
	clone->max_number_of_hits_output	= this->max_number_of_hits_output;
	clone->max_number_of_hits2_output	= this->max_number_of_hits2_output;
	
	clone->DAQSubVersion			= this->DAQSubVersion;
	clone->module_2nd				= this->module_2nd;
	clone->system_timeout			= this->system_timeout;
	clone->system_timeout_output	= this->system_timeout_output;
	clone->common_mode				= this->common_mode;
	clone->common_mode_output		= this->common_mode_output;
	clone->DAQ_info_Length			= this->DAQ_info_Length;
	clone->Camac_CIF_Length			= this->Camac_CIF_Length;
	clone->LMF_Version				= this->LMF_Version;
	clone->LMF_Version_output		= this->LMF_Version_output;
	clone->TDCDataType				= this->TDCDataType;
	
	clone->LMF_Header_version		= this->LMF_Header_version;
	
	clone->number_of_bytes_in_PostEventData = this->number_of_bytes_in_PostEventData;
	
	clone->tdcresolution			= this->tdcresolution;
	clone->tdcresolution_output		= this->tdcresolution_output;
	clone->frequency				= this->frequency;
	clone->DOUBLE_timestamp			= this->DOUBLE_timestamp;
	clone->ui64_timestamp			= this->ui64_timestamp;

	clone->number_of_CCFHistory_strings = this->number_of_CCFHistory_strings;
	clone->number_of_DAN_source_strings = this->number_of_DAN_source_strings;
	clone->number_of_DAQ_source_strings = this->number_of_DAQ_source_strings;
	clone->number_of_CCFHistory_strings_output = this->number_of_CCFHistory_strings_output;
	clone->number_of_DAN_source_strings_output = this->number_of_DAN_source_strings_output;
	clone->number_of_DAQ_source_strings_output = this->number_of_DAQ_source_strings_output;

	if (number_of_CCFHistory_strings >= 0) {
		clone->CCFHistory_strings = new std::string*[number_of_CCFHistory_strings];
		memset(clone->CCFHistory_strings,0,sizeof(std::string*)*number_of_CCFHistory_strings);
	}
	if (number_of_DAN_source_strings >= 0) {
		clone->DAN_source_strings = new std::string*[number_of_DAN_source_strings];
		memset(clone->DAN_source_strings,0,sizeof(std::string*)*number_of_DAN_source_strings);
	}
	if (number_of_DAQ_source_strings >= 0) {
		clone->DAQ_source_strings = new std::string*[number_of_DAQ_source_strings];
		memset(clone->DAQ_source_strings,0,sizeof(std::string*)*number_of_DAQ_source_strings);
	}
	if (this->CCFHistory_strings) {
		for (__int32 i=0;i<number_of_CCFHistory_strings;++i) {clone->CCFHistory_strings[i] = new std::string(); *clone->CCFHistory_strings[i] = *this->CCFHistory_strings[i];}
	}
	if (this->DAN_source_strings) {
		for (__int32 i=0;i<number_of_DAN_source_strings;++i) {clone->DAN_source_strings[i] = new std::string(); *clone->DAN_source_strings[i] = *this->DAN_source_strings[i];}
	}
	if (this->DAQ_source_strings) {
		for (__int32 i=0;i<number_of_DAQ_source_strings;++i) {clone->DAQ_source_strings[i] = new std::string(); *clone->DAQ_source_strings[i] = *this->DAQ_source_strings[i];}
	}

	if (number_of_CCFHistory_strings_output>=0) {
		clone->CCFHistory_strings_output = new std::string*[number_of_CCFHistory_strings_output];
		memset(clone->CCFHistory_strings_output,0,sizeof(std::string*)*number_of_CCFHistory_strings_output);
	}
	if (number_of_DAN_source_strings_output>=0) {
		clone->DAN_source_strings_output = new std::string*[number_of_DAN_source_strings_output];
		memset(clone->DAN_source_strings_output,0,sizeof(std::string*)*number_of_DAN_source_strings_output);
	}
	if (number_of_DAQ_source_strings_output>=0) {
		clone->DAQ_source_strings_output = new std::string*[number_of_DAQ_source_strings_output];
		memset(clone->DAQ_source_strings_output,0,sizeof(std::string*)*number_of_DAQ_source_strings_output);
	}
	if (this->CCFHistory_strings_output) {
		for (__int32 i=0;i<number_of_CCFHistory_strings_output;++i) {clone->CCFHistory_strings_output[i] = new std::string(); *clone->CCFHistory_strings_output[i] = *this->CCFHistory_strings_output[i];}
	}
	if (this->DAN_source_strings_output) {
		for (__int32 i=0;i<number_of_DAN_source_strings_output;++i) {clone->DAN_source_strings_output[i] = new std::string(); *clone->DAN_source_strings_output[i] = *this->DAN_source_strings_output[i];}
	}
	if (this->DAQ_source_strings_output) {
		for (__int32 i=0;i<number_of_DAQ_source_strings_output;++i) {clone->DAQ_source_strings_output[i] = new std::string(); *clone->DAQ_source_strings_output[i] = *this->DAQ_source_strings_output[i];}
	}

	clone->CCF_HISTORY_CODE_bitmasked = this->CCF_HISTORY_CODE_bitmasked;
	clone->DAN_SOURCE_CODE_bitmasked = this->DAN_SOURCE_CODE_bitmasked;
	clone->DAQ_SOURCE_CODE_bitmasked = this->DAQ_SOURCE_CODE_bitmasked;
	
	clone->errorflag				= this->errorflag;
	clone->skip_header				= this->skip_header;
	
	clone->uint64_number_of_written_events	= this->uint64_number_of_written_events;
	
	clone->not_Cobold_LMF			= this->not_Cobold_LMF;
	clone->Headersize_output		= this->Headersize_output;
	clone->output_byte_counter		= this->output_byte_counter;
	clone->Numberofcoordinates_output	= this->Numberofcoordinates_output;
	clone->must_read_first			= this->must_read_first;

	clone->num_channels				= this->num_channels;
	clone->num_ions					= this->num_ions;



// pointers and other special stuff:


//clone->InputFileIsOpen		= this->InputFileIsOpen;
//clone->in_ar					= this->in_ar;
//clone->input_lmf				= this->input_lmf;

//clone->OutputFileIsOpen			= this->OutputFileIsOpen;
//clone->out_ar					= this->out_ar;
//clone->output_lmf				= this->output_lmf;


	// TDC8PCI2
	clone->TDC8PCI2.GateDelay_1st_card			= this->TDC8PCI2.GateDelay_1st_card;
	clone->TDC8PCI2.OpenTime_1st_card			= this->TDC8PCI2.OpenTime_1st_card;
	clone->TDC8PCI2.WriteEmptyEvents_1st_card	= this->TDC8PCI2.WriteEmptyEvents_1st_card;
	clone->TDC8PCI2.TriggerFalling_1st_card		= this->TDC8PCI2.TriggerFalling_1st_card;
	clone->TDC8PCI2.TriggerRising_1st_card		= this->TDC8PCI2.TriggerRising_1st_card;
	clone->TDC8PCI2.EmptyCounter_1st_card		= this->TDC8PCI2.EmptyCounter_1st_card;
	clone->TDC8PCI2.EmptyCounter_since_last_Event_1st_card = this->TDC8PCI2.EmptyCounter_since_last_Event_1st_card;
	clone->TDC8PCI2.use_normal_method			= this->TDC8PCI2.use_normal_method;
	clone->TDC8PCI2.use_normal_method_2nd_card	= this->TDC8PCI2.use_normal_method_2nd_card;
	clone->TDC8PCI2.sync_test_on_off			= this->TDC8PCI2.sync_test_on_off;
	clone->TDC8PCI2.io_address_2nd_card			= this->TDC8PCI2.io_address_2nd_card;
	clone->TDC8PCI2.GateDelay_2nd_card			= this->TDC8PCI2.GateDelay_2nd_card;
	clone->TDC8PCI2.OpenTime_2nd_card			= this->TDC8PCI2.OpenTime_2nd_card;
	clone->TDC8PCI2.WriteEmptyEvents_2nd_card	= this->TDC8PCI2.WriteEmptyEvents_2nd_card;
	clone->TDC8PCI2.TriggerFallingEdge_2nd_card = this->TDC8PCI2.TriggerFallingEdge_2nd_card;
	clone->TDC8PCI2.TriggerRisingEdge_2nd_card	= this->TDC8PCI2.TriggerRisingEdge_2nd_card;
	clone->TDC8PCI2.EmptyCounter_2nd_card		= this->TDC8PCI2.EmptyCounter_2nd_card;
	clone->TDC8PCI2.EmptyCounter_since_last_Event_2nd_card = this->TDC8PCI2.EmptyCounter_since_last_Event_2nd_card;
	clone->TDC8PCI2.variable_event_length		= this->TDC8PCI2.variable_event_length;


	// HM1
	clone->HM1.FAK_DLL_Value			= this->HM1.FAK_DLL_Value;
	clone->HM1.Resolution_Flag			= this->HM1.Resolution_Flag;
	clone->HM1.trigger_mode_for_start	= this->HM1.trigger_mode_for_start;
	clone->HM1.trigger_mode_for_stop	= this->HM1.trigger_mode_for_stop;
	clone->HM1.Even_open_time			= this->HM1.Even_open_time;
	clone->HM1.Auto_Trigger				= this->HM1.Auto_Trigger;
	clone->HM1.set_bits_for_GP1			= this->HM1.set_bits_for_GP1;
	clone->HM1.ABM_m_xFrom				= this->HM1.ABM_m_xFrom;
	clone->HM1.ABM_m_xTo				= this->HM1.ABM_m_xTo;
	clone->HM1.ABM_m_yFrom				= this->HM1.ABM_m_yFrom;
	clone->HM1.ABM_m_yTo				= this->HM1.ABM_m_yTo;
	clone->HM1.ABM_m_xMin				= this->HM1.ABM_m_xMin;
	clone->HM1.ABM_m_xMax				= this->HM1.ABM_m_xMax;
	clone->HM1.ABM_m_yMin				= this->HM1.ABM_m_yMin;
	clone->HM1.ABM_m_yMax				= this->HM1.ABM_m_yMax;
	clone->HM1.ABM_m_xOffset			= this->HM1.ABM_m_xOffset;
	clone->HM1.ABM_m_yOffset			= this->HM1.ABM_m_yOffset;
	clone->HM1.ABM_m_zOffset			= this->HM1.ABM_m_zOffset;
	clone->HM1.ABM_Mode					= this->HM1.ABM_Mode;
	clone->HM1.ABM_OsziDarkInvert		= this->HM1.ABM_OsziDarkInvert;
	clone->HM1.ABM_ErrorHisto			= this->HM1.ABM_ErrorHisto;
	clone->HM1.ABM_XShift				= this->HM1.ABM_XShift;
	clone->HM1.ABM_YShift				= this->HM1.ABM_YShift;
	clone->HM1.ABM_ZShift				= this->HM1.ABM_ZShift;
	clone->HM1.ABM_ozShift				= this->HM1.ABM_ozShift;
	clone->HM1.ABM_wdShift				= this->HM1.ABM_wdShift;
	clone->HM1.ABM_ucLevelXY			= this->HM1.ABM_ucLevelXY;
	clone->HM1.ABM_ucLevelZ				= this->HM1.ABM_ucLevelZ;
	clone->HM1.ABM_uiABMXShift			= this->HM1.ABM_uiABMXShift;
	clone->HM1.ABM_uiABMYShift			= this->HM1.ABM_uiABMYShift;
	clone->HM1.ABM_uiABMZShift			= this->HM1.ABM_uiABMZShift;
	clone->HM1.use_normal_method		= this->HM1.use_normal_method;

	clone->HM1.TWOHM1_FAK_DLL_Value		= this->HM1.TWOHM1_FAK_DLL_Value;
	clone->HM1.TWOHM1_Resolution_Flag	= this->HM1.TWOHM1_Resolution_Flag;
	clone->HM1.TWOHM1_trigger_mode_for_start	= this->HM1.TWOHM1_trigger_mode_for_start;
	clone->HM1.TWOHM1_trigger_mode_for_stop		= this->HM1.TWOHM1_trigger_mode_for_stop;
	clone->HM1.TWOHM1_res_adjust		= this->HM1.TWOHM1_res_adjust;
	clone->HM1.TWOHM1_tdcresolution		= this->HM1.TWOHM1_tdcresolution;
	clone->HM1.TWOHM1_test_overflow		= this->HM1.TWOHM1_test_overflow;
	clone->HM1.TWOHM1_number_of_channels	= this->HM1.TWOHM1_number_of_channels;
	clone->HM1.TWOHM1_number_of_hits	= this->HM1.TWOHM1_number_of_hits;
	clone->HM1.TWOHM1_set_bits_for_GP1	= this->HM1.TWOHM1_set_bits_for_GP1;
	clone->HM1.TWOHM1_HM1_ID_1			= this->HM1.TWOHM1_HM1_ID_1;
	clone->HM1.TWOHM1_HM1_ID_2			= this->HM1.TWOHM1_HM1_ID_2;

	clone->TDC8HP.no_config_file_read	= this->TDC8HP.no_config_file_read;
	clone->TDC8HP.RisingEnable_p61		= this->TDC8HP.RisingEnable_p61;
	clone->TDC8HP.FallingEnable_p62		= this->TDC8HP.FallingEnable_p62;
	clone->TDC8HP.TriggerEdge_p63		= this->TDC8HP.TriggerEdge_p63;
	clone->TDC8HP.TriggerChannel_p64		= this->TDC8HP.TriggerChannel_p64;
	clone->TDC8HP.OutputLevel_p65		= this->TDC8HP.OutputLevel_p65;
	clone->TDC8HP.GroupingEnable_p66		= this->TDC8HP.GroupingEnable_p66;
	clone->TDC8HP.GroupingEnable_p66_output	= this->TDC8HP.GroupingEnable_p66_output;
	clone->TDC8HP.AllowOverlap_p67		= this->TDC8HP.AllowOverlap_p67;
	clone->TDC8HP.TriggerDeadTime_p68	= this->TDC8HP.TriggerDeadTime_p68;
	clone->TDC8HP.GroupRangeStart_p69	= this->TDC8HP.GroupRangeStart_p69;
	clone->TDC8HP.GroupRangeEnd_p70		= this->TDC8HP.GroupRangeEnd_p70;
	clone->TDC8HP.ExternalClock_p71		= this->TDC8HP.ExternalClock_p71;
	clone->TDC8HP.OutputRollOvers_p72	= this->TDC8HP.OutputRollOvers_p72;
	clone->TDC8HP.DelayTap0_p73			= this->TDC8HP.DelayTap0_p73;
	clone->TDC8HP.DelayTap1_p74			= this->TDC8HP.DelayTap1_p74;
	clone->TDC8HP.DelayTap2_p75			= this->TDC8HP.DelayTap2_p75;
	clone->TDC8HP.DelayTap3_p76			= this->TDC8HP.DelayTap3_p76;
	clone->TDC8HP.INL_p80				= this->TDC8HP.INL_p80;
	clone->TDC8HP.DNL_p81				= this->TDC8HP.DNL_p81;
	clone->TDC8HP.csConfigFile			= this->TDC8HP.csConfigFile;
	clone->TDC8HP.csINLFile				= this->TDC8HP.csINLFile;
	clone->TDC8HP.csDNLFile				= this->TDC8HP.csDNLFile;
	clone->TDC8HP.csConfigFile_Length	= this->TDC8HP.csConfigFile_Length;
	clone->TDC8HP.csINLFile_Length		= this->TDC8HP.csINLFile_Length;
	clone->TDC8HP.csDNLFile_Length		= this->TDC8HP.csDNLFile_Length;
	clone->TDC8HP.UserHeaderVersion		= this->TDC8HP.UserHeaderVersion;
	clone->TDC8HP.VHR_25ps				= this->TDC8HP.VHR_25ps;
	clone->TDC8HP.SyncValidationChannel	= this->TDC8HP.SyncValidationChannel;
	clone->TDC8HP.variable_event_length	= this->TDC8HP.variable_event_length;
	clone->TDC8HP.SSEEnable				= this->TDC8HP.SSEEnable;
	clone->TDC8HP.MMXEnable				= this->TDC8HP.MMXEnable;
	clone->TDC8HP.DMAEnable				= this->TDC8HP.DMAEnable;
	clone->TDC8HP.GroupTimeOut			= this->TDC8HP.GroupTimeOut;

	clone->TDC8HP.i32NumberOfDAQLoops   = this->TDC8HP.i32NumberOfDAQLoops;   
	clone->TDC8HP.TDC8HP_DriverVersion 	= this->TDC8HP.TDC8HP_DriverVersion; 	
	clone->TDC8HP.iTriggerChannelMask	= this->TDC8HP.iTriggerChannelMask;	
	clone->TDC8HP.iTime_zero_channel	= this->TDC8HP.iTime_zero_channel;	

	clone->TDC8HP.Number_of_TDCs = this->TDC8HP.Number_of_TDCs;
	for (__int32 i = 0;i<3;++i) {
		if (this->TDC8HP.TDC_info[i]) {
			clone->TDC8HP.TDC_info[i]->index				= this->TDC8HP.TDC_info[i]->index;
			clone->TDC8HP.TDC_info[i]->channelCount			= this->TDC8HP.TDC_info[i]->channelCount;
			clone->TDC8HP.TDC_info[i]->channelStart			= this->TDC8HP.TDC_info[i]->channelStart;
			clone->TDC8HP.TDC_info[i]->highResChannelCount	= this->TDC8HP.TDC_info[i]->highResChannelCount;
			clone->TDC8HP.TDC_info[i]->highResChannelStart	= this->TDC8HP.TDC_info[i]->highResChannelStart;
			clone->TDC8HP.TDC_info[i]->lowResChannelCount	= this->TDC8HP.TDC_info[i]->lowResChannelCount;
			clone->TDC8HP.TDC_info[i]->lowResChannelStart	= this->TDC8HP.TDC_info[i]->lowResChannelStart;
			clone->TDC8HP.TDC_info[i]->resolution			= this->TDC8HP.TDC_info[i]->resolution;
			clone->TDC8HP.TDC_info[i]->serialNumber			= this->TDC8HP.TDC_info[i]->serialNumber;
			clone->TDC8HP.TDC_info[i]->version				= this->TDC8HP.TDC_info[i]->version;
			clone->TDC8HP.TDC_info[i]->fifoSize				= this->TDC8HP.TDC_info[i]->fifoSize;
			clone->TDC8HP.TDC_info[i]->flashValid			= this->TDC8HP.TDC_info[i]->flashValid;
			memcpy(clone->TDC8HP.TDC_info[i]->INLCorrection,this->TDC8HP.TDC_info[i]->INLCorrection,sizeof(__int32)*8*1024);
			memcpy(clone->TDC8HP.TDC_info[i]->DNLData,this->TDC8HP.TDC_info[i]->DNLData,sizeof(__int16)*8*1024);
		}
	}




	clone->fADC8.driver_version					=	this->fADC8.driver_version;
	clone->fADC8.i32NumberOfDAQLoops			=	this->fADC8.i32NumberOfDAQLoops;
	clone->fADC8.number_of_bools				=	this->fADC8.number_of_bools;
	clone->fADC8.number_of_int32s				=	this->fADC8.number_of_int32s;
	clone->fADC8.number_of_uint32s				=	this->fADC8.number_of_uint32s;
	clone->fADC8.number_of_doubles				=	this->fADC8.number_of_doubles;
	clone->fADC8.GroupEndMarker					=	this->fADC8.GroupEndMarker;
	clone->fADC8.i32NumberOfADCmodules			=	this->fADC8.i32NumberOfADCmodules;
	clone->fADC8.iEnableGroupMode				=	this->fADC8.iEnableGroupMode;
	clone->fADC8.iTriggerChannel				=	this->fADC8.iTriggerChannel;
	clone->fADC8.iPreSamplings_in_4800ps_units	=	this->fADC8.iPreSamplings_in_4800ps_units;
	clone->fADC8.iPostSamplings_in_9600ps_units	=	this->fADC8.iPostSamplings_in_9600ps_units;
	clone->fADC8.iEnableTDCinputs				=	this->fADC8.iEnableTDCinputs;
	clone->fADC8.bReadCustomData				=	this->fADC8.bReadCustomData;
	clone->fADC8.veto_gate_length				=	this->fADC8.veto_gate_length;
	clone->fADC8.veto_delay_length				=	this->fADC8.veto_delay_length;
	clone->fADC8.veto_mask						=	this->fADC8.veto_mask;
	clone->fADC8.dGroupRangeStart				=	this->fADC8.dGroupRangeStart;
	clone->fADC8.dGroupRangeEnd					=	this->fADC8.dGroupRangeEnd;
	clone->fADC8.at_least_1_signal_was_written	=	this->fADC8.at_least_1_signal_was_written;
	for (__int32 m = 0;m<8;m++) {
		clone->fADC8.firmware_version[m]	= this->fADC8.firmware_version[m];
		clone->fADC8.serial_number[m]		= this->fADC8.serial_number[m];
		for (__int32 adc = 0;adc<2;adc++) {
			clone->fADC8.dSyncTimeOffset[m][adc] = this->fADC8.dSyncTimeOffset[m][adc];
			clone->fADC8.iChannelMode[m][adc] =	   this->fADC8.iChannelMode[m][adc];								 // 0 = 1.25Gs, 1 = 2.5Gs, 2 = 5Gs
			clone->fADC8.iThreshold_GT[m][m] =		this->fADC8.iThreshold_GT[m][m];
			clone->fADC8.iThreshold_LT[m][m] =		this->fADC8.iThreshold_LT[m][m];
			clone->fADC8.iSynchronMode [m][adc] =	 this ->fADC8.iSynchronMode [m][adc];
		}
		for (__int32 ch = 0;ch<10;ch++) clone->fADC8.GND_level[m][ch] = this->fADC8.GND_level[m][ch];
	}


	clone->fADC4.packet_count =				this->fADC4.packet_count;
	clone->fADC4.number_of_bools =			this->fADC4.number_of_bools;
	clone->fADC4.number_of_int32s =			this->fADC4.number_of_int32s;
	clone->fADC4.number_of_uint32s =		this->fADC4.number_of_uint32s;
	clone->fADC4.number_of_doubles =		this->fADC4.number_of_doubles;
	clone->fADC4.GroupEndMarker =			this->fADC4.GroupEndMarker;
	clone->fADC4.driver_version =			this->fADC4.driver_version;
	clone->fADC4.i32NumberOfDAQLoops =		this->fADC4.i32NumberOfDAQLoops;
	clone->fADC4.bReadCustomData =			this->fADC4.bReadCustomData;
	clone->fADC4.i32NumberOfADCmodules =	this->fADC4.i32NumberOfADCmodules;
	clone->fADC4.iTriggerChannel =			this->fADC4.iTriggerChannel;
	clone->fADC4.dGroupRangeStart =			this->fADC4.dGroupRangeStart;
	clone->fADC4.dGroupRangeEnd =			this->fADC4.dGroupRangeEnd;
	clone->fADC4.csConfigFile =				this->fADC4.csConfigFile;
	clone->fADC4.csINLFile =				this->fADC4.csINLFile;
	clone->fADC4.csDNLFile =				this->fADC4.csDNLFile;
	for (__int32 i=0;i<20;i++) clone->fADC4.bits_per_mVolt[i] = this->fADC4.bits_per_mVolt[i];

	return true;
}






//////////////////////////////////////////////////////////////////
// decompress_asynchronous_fADC8_signal()
//////////////////////////////////////////////////////////////////

bool LMF_IO::decompress_asynchronous_fADC8_signal(fADC8_signal_info_struct &signal_info, unsigned __int32 source_buffer[], unsigned __int32 source_buffer_size_in_32bit_words, __int16 i16bit_target_buffer[], __int32 target_buffer_size_in_16bit_words, __int32 &number_of_filled_16bit_words)
{
	number_of_filled_16bit_words = -signal_info.timestamp_subindex - 1;

	if ((signal_info.signal_type == FADC8_HEADER_EVENT_ID_TDC) || (signal_info.signal_type == FADC8_HEADER_EVENT_ID_END_MARKER)) {
		number_of_filled_16bit_words = 0;
		return false;
	}

	unsigned __int32 limit_;
	if (signal_info.signallength_including_header_in_32bit_words - 4 < __int32(source_buffer_size_in_32bit_words)) {
		limit_ = signal_info.signallength_including_header_in_32bit_words - 4;
	} else {
		limit_ = source_buffer_size_in_32bit_words;
	}

	if (   (signal_info.signal_type == FADC8_HEADER_EVENT_ID_1_25G_ADC1)
		|| (signal_info.signal_type == FADC8_HEADER_EVENT_ID_1_25G_ADC2)
		|| (signal_info.signal_type == FADC8_HEADER_EVENT_ID_1_25G_ADC3)
		|| (signal_info.signal_type == FADC8_HEADER_EVENT_ID_1_25G_ADC4)) {

		for (unsigned __int32 i=4;i<limit_ + 4;i+=4)  {
			if (target_buffer_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*3) {
				break;
			}
			for (__int32 j=0; j<4; j++) {
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+j] >> 20) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+j] >> 10) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+j]      ) & 0x3ff;
			}
		}

		number_of_filled_16bit_words++;
		//trace_filter(i16bit_target_buffer,number_of_filled_16bit_words, signal_info);
		return true;
	} // endif 1.25 GS mode


	if (   (signal_info.signal_type == FADC8_HEADER_EVENT_ID_2_5G_ADC12)
		|| (signal_info.signal_type == FADC8_HEADER_EVENT_ID_2_5G_ADC34)) {
	
		for (unsigned __int32 i=4;i<limit_ + 4;i+=8)  {
			if (target_buffer_size_in_16bit_words < number_of_filled_16bit_words +1 + 2*4*3) {
				break;
			}

			for (__int32 j=0; j<4; j++) {
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0+j] >> 20) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4+j] >> 20) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0+j] >> 10) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4+j] >> 10) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0+j]      ) & 0x3ff;
				i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4+j]      ) & 0x3ff;
			}
			
		}

		number_of_filled_16bit_words++;
		//trace_filter(i16bit_target_buffer,number_of_filled_16bit_words, signal_info);
		return true;
	} // endif 2.5 Gs mode

	if (signal_info.signal_type == FADC8_HEADER_EVENT_ID_5G_ADC1234) {

		for (unsigned __int32 i=4;i<limit_ + 4;i+=16)  {
			if (target_buffer_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*4*3) {
				break;
			}

			for (__int32 j=0; j<4; j++) {
				for (__int32 shift_ = 20;shift_ >= 0;shift_-=10) {
					i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+8 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+12+j] >> shift_) & 0x3ff;
				}
			}

		}

		number_of_filled_16bit_words++;
		//trace_filter(i16bit_target_buffer,number_of_filled_16bit_words, signal_info);
		return true;
	} // endif 5 Gs mode

	return false;
}



















//////////////////////////////////////////////////////////////////
// decompress_synchronous_fADC8_signal()
//////////////////////////////////////////////////////////////////

bool LMF_IO::decompress_synchronous_fADC8_signal(fADC8_signal_info_struct &signal_info, unsigned __int32 source_buffer[], unsigned __int32 source_buffer_size_in_32bit_words, __int32 &number_of_filled_16bit_words,
												__int16 i16bit_target_buffer1[], __int32 target_buffer1_size_in_16bit_words,
												__int16 i16bit_target_buffer2[], __int32 target_buffer2_size_in_16bit_words,
												__int16 i16bit_target_buffer3[], __int32 target_buffer3_size_in_16bit_words,
												__int16 i16bit_target_buffer4[], __int32 target_buffer4_size_in_16bit_words)
{

	number_of_filled_16bit_words = -1;

	if (signal_info.signal_type != FADC8_HEADER_EVENT_ID_5G_ADC1234) {
		number_of_filled_16bit_words = 0;
		return false;
	}

	unsigned __int32 limit_;
	if (signal_info.signallength_including_header_in_32bit_words - 4 < __int32(source_buffer_size_in_32bit_words)) {
		limit_ = signal_info.signallength_including_header_in_32bit_words - 4;
	} else {
		limit_ = source_buffer_size_in_32bit_words;
	}


	if (signal_info.megasamples_per_second == 5000) {

		for (unsigned __int32 i=4;i<limit_ + 4-16;i+=16)  { // xxx stimmt -16?
			if (target_buffer1_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*4*3) break;

			for (__int32 j=0; j<4; j++) {
				for (__int32 shift_ = 20;shift_ >= 0;shift_-=10) {
					i16bit_target_buffer1[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer1[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+8 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer1[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer1[(++number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+12+j] >> shift_) & 0x3ff;
				}
			}

		}
		number_of_filled_16bit_words++;
		//trace_filter(i16bit_target_buffer1,number_of_filled_16bit_words, signal_info);
		return true;
	}





	if (signal_info.megasamples_per_second == 2500) {
		for (unsigned __int32 i=4;i<limit_ + 4-16;i+=16)  {
			if (target_buffer1_size_in_16bit_words < number_of_filled_16bit_words +1 + 2*4*3) break;
			if (target_buffer3_size_in_16bit_words < number_of_filled_16bit_words +1 + 2*4*3) break;

			for (__int32 j=0; j<4; j++) {
				for (__int32 shift_ = 20;shift_ >= 0;shift_-=10) {
					++number_of_filled_16bit_words;
					i16bit_target_buffer1[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer3[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+8 +j] >> shift_) & 0x3ff;

					++number_of_filled_16bit_words;
					i16bit_target_buffer1[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4 +j] >> shift_) & 0x3ff;					
					i16bit_target_buffer3[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+12+j] >> shift_) & 0x3ff;
				}	
			}

		}
		//trace_filter(i16bit_target_buffer1,number_of_filled_16bit_words, signal_info);
		//trace_filter(i16bit_target_buffer3,number_of_filled_16bit_words, signal_info);
		number_of_filled_16bit_words++;
		return true;
	}


		if (signal_info.megasamples_per_second == 1250) {
		for (unsigned __int32 i=4;i<limit_ + 4-16;i+=16)  {
			if (target_buffer1_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*3) break;
			if (target_buffer2_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*3) break;
			if (target_buffer3_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*3) break;
			if (target_buffer4_size_in_16bit_words < number_of_filled_16bit_words +1 + 4*3) break;

			for (__int32 j=0; j<4; j++) {
				for (__int32 shift_ = 20;shift_ >= 0;shift_-=10) {
					++number_of_filled_16bit_words;
					i16bit_target_buffer1[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+0 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer2[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+4 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer3[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+8 +j] >> shift_) & 0x3ff;
					i16bit_target_buffer4[(number_of_filled_16bit_words) > 0 ? number_of_filled_16bit_words : 0] = (source_buffer[i+12+j] >> shift_) & 0x3ff;
				}
			}

		}
		//trace_filter(i16bit_target_buffer1,number_of_filled_16bit_words, signal_info);
		//trace_filter(i16bit_target_buffer2,number_of_filled_16bit_words, signal_info);
		//trace_filter(i16bit_target_buffer3,number_of_filled_16bit_words, signal_info);
		//trace_filter(i16bit_target_buffer4,number_of_filled_16bit_words, signal_info);
		number_of_filled_16bit_words++;
		return true;
	}

	number_of_filled_16bit_words = 0;
	return false;
}



   �n  H   I D S _ S O U R C E _ C P P   L M F _ I O _ H       0         #pragma once

#ifndef _LMF_IO_
	#define _LMF_IO_
#include "fstream"
//#include "stdio.h"
#include "time.h"

#ifndef LINUX
	#ifndef WIN32
		#ifndef WIN64
			#define LINUX
		#endif
	#endif
#endif


#ifdef LINUX
	#include "string.h"
	#define _fseeki64 fseeko
	#define _ftelli64 ftello

	#ifndef __int32_IS_DEFINED
		#define __int32_IS_DEFINED
		#define __int32 int
		#define __int16 short
		#define __int64 long long
		#define __int8 char
	#endif
#endif

#ifndef LINUX
#pragma warning(disable : 4996)
#endif





#define CRONO_OK 0
#define CRONO_WINDRIVER_NOT_FOUND 1
#define CRONO_DEVICE_NOT_FOUND 2
#define CRONO_NOT_INITIALIZED 3
	/** No card with the matching board id was found in the device list*/
#define CRONO_ERROR_TRIGGER_CARD_NOT_FOUND 4
/* hen a capture on a closed card is called*/
#define CRONO_WRONG_STATE 4
/* The pointer given to a d125 function was not a valid pointer*/
#define CRONO_INVALID_DEVICE 5
#define CRONO_SYNC_NDIGO_MASTER_REQUIRED 6
#define CRONO_SYNC_INVALID_CONFIG 100


#define CRONO_PACKET_TYPE_16_BIT_SIGNED 1
#define CRONO_PACKET_TYPE_TIMESTAMP_ONLY 128
#define CRONO_PACKET_TYPE_TDC_RISING 144
#define CRONO_PACKET_TYPE_TDC_FALLING 145
#define CRONO_PACKET_TYPE_END_OF_BUFFER 129
#define CRONO_PACKET_TYPE_TDC_DATA	8
#define CRONO_READ_OK 0
#define CRONO_READ_NO_DATA 1
#define CRONO_READ_INTERNAL_ERROR 2


struct ndigo_packet {
	unsigned char channel;  // channel number ADC: 0-3 , 4=TDC, 5="Gate"-input
	unsigned char card;		// card number
	unsigned char type;		// signal type (adc, TDC or G)
	unsigned char flags;	// error handling .. overflows usw.
	unsigned int length;	// number of samples divided by 4.
	unsigned __int64 timestamp;		// pico sec. absolut time stamp of first bin
	unsigned __int64 data[1];		// raw data of ADC trace
};

typedef struct {

	/* The number of bytes occupied by the structure*/
	int size;
	/*	A version number that is increased when the definition of the structure is changed. 
		The increment can be larger than one to match driver version numbers or similar. Set to 0 for all versions up to first release.
	*/
	int version;

	/* Bandwidth. Currently fixed at 3G, might later be configurable to 1G (fullBW=false).*/
	double bandwidth;

	/*	Actual sample rate of currently sampled data. This is affected by the number of channels in use
		in the current ADC mode and possibly also by changes to the clock.
	*/
	double sample_rate;
	/* The period one sample in the data represent in picoseconds */
	double sample_period;

	/* the ID the board should use to identify itself in the output data stream. 0 to 255. */
	int board_id;
	/* Number of channels in the current mode.*/
	int channels;

	/*

	*/
	int channel_mask;

	/** The total amount of buffer in bytes*/
	__int64 total_buffer;

	__int64 free_buffer;
} ndigo_param_info;



struct ndigo_static_info {
	/*! \brief The number of bytes occupied by the structure
	 */
	int size;

	/*! \brief A version number
	 * 
	 * that is increased when the definition of the structure is changed. 
	 * The increment can be larger than one to match driver version numbers or similar. 
	 * Set to 0 for all versions up to first release.
	 */
	int version;
	
	/*! \brief Index of the board as passed to the constructor
	 *
	 * or set via int @link conffuncts ndigo_set_board_id(ndigo_device *device, int board_id) @endlink.
	 */
	int board_id;
	
	/*! \brief driver revision number
	 *
	 * The lower three bytes contain a triple level hierarchy of version numbers. * E.g. 0x010103 codes
	 * version 1.1.3. A change in the first digit generally requires a recompilation of user applications.
	 * Change in the second digit denote significant improvements or changes that don't break compatibility
	 * and the third digit changes with minor bugfixes and the like.
	 */
	int driver_revision;
	
	/*! \brief Revision number of the FPGA configuration.
	 *
	 * This can be read from a register.
	 */
	int firmware_revision;
	
	/*! \brief board revision number
	 *
	 * This can be read from a register. It is a four bit number that changes when the schematic of the
	 * board is changed. 
	 *  - 0: Experimental first board Version. Labeled Rev. 1
	 *  - 2: First commercial Version. Labeled "Rev. 2"
	 *  - 3: for the version produced starting in 2011 labeled "Rev. 3"
	 */
	int board_revision;
	
	/*! \brief The same board schematic can be populated in multiple Variants.
	 *
	 * This is a four bit code that	can be read from a register.
	 * - For board revision 0 this always reads 0. 
	 * - For board revision 2 
	 *		- bit 0 determines whether this board contains an 8-bit or 10-bit ADC
	 *		- bit 1 determines whether differential inputs are used
	 *		- bit 2 determines whether the tdc-oscillator is present
	 *		- bit 3 = 1 signifies a special version of the board.
	 * - For board revision 3
	 *		- bit 2 determines input connectors (0: single ende, 1: differential)
	 *		- for the other bits see user guide
	*/
	int board_configuration;
	
	/*! \brief Number of bits of the ADC
	 *
	 * set to 0 if unknown. Should be 14.
	 */
	int adc_resolution;
	
	/*! \brief Maximum sample rate.
	 *	- 2:5e8 = 250Msps for the Ndigo250M,
	 *	- 1:25e8 = 125Msps for the Ndigo125M 
	 *	- 5e9 = 5Gsps for the Ndigo5G
	 */
	double nominal_sample_rate;	
	
	/*! \brief analog bandwidth
	 *
	 *	- 1.25e8 for 125MHz for Ndigo250M 
	 *	- 3e9 for 3Ghz for Ndigo5G
	 */
	double analog_bandwidth;
	
	/*! \brief chipID as read from the 16 bit adc chip id register at SPI address 0.
	 *
	 * his value should be cached.
	 */
	int chip_id;
	
	/*! \brief Serial number with year and running number in 8.24 format.
	 *
	 * The number is identical to the one printed on the silvery sticker on the board.
	 */
	int board_serial;
	
	/*! \brief 64bit serial number from the configuration flash chip
	 */
	int flash_serial_low;

	/*! \brief 64bit serial number from the configuration flash chip
	 */
	int flash_serial_high;
	
	/*! \brief If not 0 the driver found valid calibration data in the flash on the board and is using it.
	 */
	int flash_valid;
	
	/*! \brief Returns false for the standard AC coupled Ndigo5G.
	 * 
	 * Returns true for the Ndigo250M.
	 */
	unsigned char dc_coupled;
	
	/*! \brief Subversion revision id of the FPGA configuration.
	 *
	 * This can be read from a register.
	 */
	int subversion_revision;
	
	/*! \brief calibration date
	 *
	 * DIN EN ISO 8601 string YYYY-MM-DD HH:DD describing the time when the card was calibrated.
	 */
	char calibration_date[20];
};






#define FADC8_HEADER_EVENT_ID_1_25G_ADC1		0x0
#define FADC8_HEADER_EVENT_ID_1_25G_ADC2		0x1
#define FADC8_HEADER_EVENT_ID_1_25G_ADC3		0x2
#define FADC8_HEADER_EVENT_ID_1_25G_ADC4		0x3
#define FADC8_HEADER_EVENT_ID_2_5G_ADC12		0x4
#define FADC8_HEADER_EVENT_ID_2_5G_ADC34		0x5
#define FADC8_HEADER_EVENT_ID_5G_ADC1234		0x7
#define FADC8_HEADER_EVENT_ID_TDC				0x8
#define FADC8_HEADER_EVENT_ID_DIRECT_MEMORY_START		0xc	  // add 1.2.2011
#define FADC8_HEADER_EVENT_ID_DIRECT_MEMORY_STOP		0xd	  // add 1.2.2011
#define FADC8_HEADER_EVENT_ID_ANALYZED_SIGNAL	0xe
#define FADC8_HEADER_EVENT_ID_END_MARKER		0xf



#ifndef FADC8_MANAGER_CLASS_IS_DEFINED
	struct fADC8_signal_info_struct {
		__int8	ModuleIndex;
		__int8	ADCChipIndex;
		__int64 timestamp_in_ps; // unit = 1ps
		__int64	coarse_timestamp; // units: 4.8 ns
		__int8	signal_type;
		unsigned __int32 counter;
		double  TDC_value_ps;
		__int32 signallength_including_header_in_32bit_words;
		__int8	timestamp_subindex;
		__int8	transition_type;
		__int16	adc_channel;
		__int8	TDC_subcounter;
		__int16 megasamples_per_second;
		__int8	is_a_synchronous_signal;
	//	__int32 i32dummy_for_future_use;
	//	unsigned __int32 header[4];
	};
#endif


#define MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA 4000

class MyFILE
{
public:
	MyFILE(bool mode_reading_) {error = 0; eof = false; mode_reading = mode_reading_; file = 0; position = 0; filesize = 0;}
	~MyFILE() {close(); error = 0; eof = false;}

	FILE * file;

	__int64 get_position() {
		if (!file) return 0;
		return position;
	}

	bool open(__int8* name) {
		if (file) {error = 1; return false;}
		eof = false;
		if (mode_reading) {
			file = fopen(name,"rb");
			if (!file) {error = 1; return false;}
			if (ferror(file)) {fclose(file); error = 1; return false;}
			__int8 dummy;
			if (!fread(&dummy,1,1,file)) {fclose(file); error = 1; return false;}
			if (ferror(file)) {fclose(file); error = 1; return false;}
			_fseeki64(file,  (unsigned __int64)0, SEEK_END);
			filesize = _ftelli64(file);
			_fseeki64(file,  (unsigned __int64)0, SEEK_SET);
		} else file = fopen(name,"wb");
		if (file) return true; else {error = 1; return false;}
	}

	void close() {
		if (file) {fclose(file);  file = 0;} else error = 1;
		position = 0; filesize = 0; eof = false;
	}

	unsigned __int64 tell() {return position;}
	void seek(unsigned __int64 pos);

	void read(__int8* string,__int32 length_bytes) {
		unsigned __int32 read_bytes = (unsigned __int32)(fread(string,1,length_bytes,file));
		if (__int32(read_bytes) != length_bytes) {
			error = 1;
			if (feof(file)) eof = true;
		}
		position += length_bytes;
	}

	void read(unsigned __int32 * dest,__int32 length_bytes) {
		unsigned __int32 read_bytes = (unsigned __int32)(fread(dest,1,length_bytes,file));
		if (__int32(read_bytes) != length_bytes) {
			error = 1;
			if (feof(file)) eof = true;
		}
		position += length_bytes;
	}

	void write(const __int8* string,__int32 length) {
		fwrite(string,1,length,file);
		position += length;
		if (filesize < position) filesize = position;
	}

	void flush() {fflush(file);}

	MyFILE & operator>>(unsigned __int8 &c)		{read((__int8*)&c,sizeof(unsigned __int8));		return *this;}
	MyFILE & operator>>(__int8 &c)				{read((__int8*)&c,sizeof(__int8));				return *this;}
	MyFILE & operator>>(unsigned __int16 &l)	{read((__int8*)&l,sizeof(unsigned __int16));	return *this;}
	MyFILE & operator>>(unsigned __int32 &l)	{read((__int8*)&l,sizeof(unsigned __int32));	return *this;}
	MyFILE & operator>>(unsigned __int64 &l)	{read((__int8*)&l,sizeof(unsigned __int64));	return *this;}
	MyFILE & operator>>(__int16 &s)				{read((__int8*)&s,sizeof(__int16));				return *this;}
	MyFILE & operator>>(__int32 &l)				{read((__int8*)&l,sizeof(__int32));				return *this;}
	MyFILE & operator>>(__int64 &l)				{read((__int8*)&l,sizeof(__int64));				return *this;}
	MyFILE & operator>>(double &d)				{read((__int8*)&d,sizeof(double));				return *this;}
	MyFILE & operator>>(bool &d)				{read((__int8*)&d,sizeof(bool));				return *this;}

	MyFILE & operator<<(unsigned __int8 c)		{write((__int8*)&c,sizeof(unsigned __int8));	return *this;}
	MyFILE & operator<<(__int8 c)				{write((__int8*)&c,sizeof(__int8));				return *this;}
	MyFILE & operator<<(unsigned __int16 l)		{write((__int8*)&l,sizeof(unsigned __int16));	return *this;}
	MyFILE & operator<<(unsigned __int32 l)		{write((__int8*)&l,sizeof(unsigned __int32));	return *this;}
	MyFILE & operator<<(unsigned __int64 l)		{write((__int8*)&l,sizeof(unsigned __int64));	return *this;}
	MyFILE & operator<<(__int16 s)				{write((__int8*)&s,sizeof(__int16));			return *this;}
	MyFILE & operator<<(__int32 l)				{write((__int8*)&l,sizeof(__int32));			return *this;}
	MyFILE & operator<<(__int64 l)				{write((__int8*)&l,sizeof(__int64));			return *this;}
	MyFILE & operator<<(double d)				{write((__int8*)&d,sizeof(double));				return *this;}
	MyFILE & operator<<(bool d)					{write((__int8*)&d,sizeof(bool));				return *this;}

	__int32 error;
	bool eof;
	unsigned __int64 filesize;

private:
	bool mode_reading;
	unsigned __int64 position;
};






#ifndef _WINNT_
typedef union _myLARGE_INTEGER {
	struct {
		unsigned __int32 LowPart;
		__int32 HighPart;
	}; 
	struct {
		unsigned __int32 LowPart;
		__int32 HighPart;
	} u;
	__int64 QuadPart;
} myLARGE_INTEGER,  *PmyLARGE_INTEGER;
#endif


#define DAQ_ID_HM1     0x000001
#define DAQ_ID_TDC8	   0x000002
#define DAQ_ID_CAMAC   0x000003
#define DAQ_ID_2HM1	   0x000004
#define DAQ_ID_2TDC8   0x000005
#define DAQ_ID_HM1_ABM 0x000006
#define DAQ_ID_TDC8HP  0x000008
#define DAQ_ID_TCPIP   0x000009
#define DAQ_ID_TDC8HPRAW 0x000010
#define DAQ_ID_FADC8	0x000011		// for fADC8
#define DAQ_ID_FADC4	0x000012		// for fADC4
#define DAQ_ID_AGAT  0x000013		//recorded with acquiris (agat) 

#define DAQ_ID_RAW32BIT 100
#define DAQ_ID_SIMPLE 101


#define LM_BYTE					1	//  8bit integer
#define LM_SHORT				2	// 16bit integer
#define LM_LONG					3	// 32bit integer
#define	LM_FLOAT				4   // 32bit IEEE float
#define LM_DOUBLE				5	// 64bit IEEE float
#define LM_CAMAC				6	// 24bit integer
#define LM_DOUBLELONG			7	// 64bit integer
#define LM_SBYTE				8	// signed 8bit integer
#define LM_SSHORT				9	// signed 16bit integer
#define LM_SLONG				10	// signed 32bit integer
#define LM_SDOUBLELONG			11	// signed 64bit integer
#define LM_LASTKNOWNDATAFORMAT	LM_SDOUBLELONG
#define LM_USERDEF				-1	// user will handle the reading 




struct TDC8PCI2_struct
{
	__int32 GateDelay_1st_card;
	__int32 OpenTime_1st_card;
	__int32 WriteEmptyEvents_1st_card;
	__int32 TriggerFalling_1st_card;
	__int32 TriggerRising_1st_card;
	__int32 EmptyCounter_1st_card;
	__int32 EmptyCounter_since_last_Event_1st_card;
	bool use_normal_method;
	bool use_normal_method_2nd_card;
	__int32 sync_test_on_off;
	__int32 io_address_2nd_card;
	__int32 GateDelay_2nd_card;
	__int32 OpenTime_2nd_card;
	__int32 WriteEmptyEvents_2nd_card;
	__int32 TriggerFallingEdge_2nd_card;
	__int32 TriggerRisingEdge_2nd_card;
	__int32 EmptyCounter_2nd_card;
	__int32 EmptyCounter_since_last_Event_2nd_card;
	__int32 variable_event_length;
};






struct TDC8HP_info_struct
{
	__int32 index;
	__int32 channelCount;
	__int32 channelStart;
	__int32 highResChannelCount;
	__int32 highResChannelStart;
	__int32 lowResChannelCount;
	__int32 lowResChannelStart;
	double resolution;
	__int32 serialNumber;
	__int32 version;
	__int32 fifoSize;
	__int32 *INLCorrection;
	unsigned __int16 *DNLData;
	bool flashValid;
};





struct TDC8HP_struct
{
	__int32 exotic_file_type;
	__int32 no_config_file_read;
	unsigned __int64 RisingEnable_p61;
	unsigned __int64 FallingEnable_p62;
	__int32 TriggerEdge_p63;
	__int32 TriggerChannel_p64;
	__int32 OutputLevel_p65;
	__int32 GroupingEnable_p66;
	__int32 GroupingEnable_p66_output;
	__int32 AllowOverlap_p67;
	double TriggerDeadTime_p68;
	double GroupRangeStart_p69;
	double GroupRangeEnd_p70;
	__int32 ExternalClock_p71;
	__int32 OutputRollOvers_p72;
	__int32 DelayTap0_p73;
	__int32 DelayTap1_p74;
	__int32 DelayTap2_p75;
	__int32 DelayTap3_p76;
	__int32 INL_p80;
	__int32 DNL_p81;
	double  OffsetTimeZeroChannel_s;
	__int32 BinsizeType;

	std::string	csConfigFile, csINLFile, csDNLFile;

	__int32	csConfigFile_Length, csINLFile_Length, csDNLFile_Length;
	__int32	UserHeaderVersion;
	bool VHR_25ps;
	__int32	SyncValidationChannel;
	__int32 variable_event_length;
	bool SSEEnable, MMXEnable, DMAEnable;
	double GroupTimeOut;

	__int32		Number_of_TDCs;
	TDC8HP_info_struct * TDC_info[3];

//	bool	bdummy;
//	__int32	idummy;
//	double	ddummy;


	__int32 i32NumberOfDAQLoops;
	unsigned __int32 TDC8HP_DriverVersion;
	__int32 iTriggerChannelMask;
	__int32 iTime_zero_channel;

	
	__int32 number_of_bools;
	__int32 number_of_int32s;
	__int32 number_of_doubles;

	unsigned __int32 ui32oldRollOver;
	unsigned __int64 ui64RollOvers;
	unsigned __int32 ui32AbsoluteTimeStamp;
//	unsigned __int64 ui64OldTimeStamp;
	unsigned __int64 ui64TDC8HP_AbsoluteTimeStamp;

	__int32 channel_offset_for_rising_transitions;
};





struct fADC4_struct
{
	unsigned __int32 driver_version;
	__int32 i32NumberOfADCmodules;
	__int32 i32NumberOfDAQLoops;
	bool	bReadCustomData;
	__int32 iTriggerChannel; // first has index 0
	__int32 number_of_bools;
	__int32 number_of_int32s;
	__int32 number_of_uint32s;
	__int32 number_of_doubles;
	
	double dGroupRangeStart;
	double dGroupRangeEnd;

	__int32 GroupEndMarker;

	__int32 packet_count;
	ndigo_param_info ndigo_parameters[20];
	double threshold_GT[80];
	double threshold_LT[80];
	double GND_level[80];
	double set_DC_offset[80];
	__int32 sampling_mode[20];
	double bits_per_mVolt[20];
	std::string	csConfigFile, csINLFile, csDNLFile;

	ndigo_static_info ndigo_info[10];
};




struct fADC8_struct
{
	unsigned __int32 driver_version;

	__int32 i32NumberOfDAQLoops;
	
	__int32 number_of_bools;
	__int32 number_of_int32s;
	__int32 number_of_uint32s;
	__int32 number_of_doubles;
	
	__int32 GroupEndMarker;
	
	__int32 i32NumberOfADCmodules;
	__int32 iEnableGroupMode;
	__int32 iTriggerChannel; // first has index 0
	__int32 iPreSamplings_in_4800ps_units;
	__int32 iPostSamplings_in_9600ps_units;
	__int32 iEnableTDCinputs;
	__int32 iChannelMode[8][2]; // 0 = 1.25Gs, 1 = 2.5Gs, 2 = 5Gs
	__int32 iThreshold_GT[8][8];
	__int32 iThreshold_LT[8][8];
	__int32 iSynchronMode [8][2];
	bool	bReadCustomData;

	__int32 firmware_version[8];
	__int32 serial_number[8];

	__int32 veto_gate_length;
	__int32 veto_delay_length;
	__int32 veto_mask;

	__int32 GND_level[8][10];
	
	double dGroupRangeStart;
	double dGroupRangeEnd;
	double dSyncTimeOffset[8][2];
	bool at_least_1_signal_was_written;
};








struct HM1_struct
{
	__int32 FAK_DLL_Value;
	__int32 Resolution_Flag;
	__int32 trigger_mode_for_start;
	__int32 trigger_mode_for_stop;
	__int32 Even_open_time;
	__int32 Auto_Trigger;
	__int32 set_bits_for_GP1;
	__int32 ABM_m_xFrom;	
	__int32 ABM_m_xTo;
	__int32 ABM_m_yFrom;
	__int32 ABM_m_yTo;
	__int32 ABM_m_xMin;
	__int32 ABM_m_xMax;
	__int32 ABM_m_yMin;
	__int32 ABM_m_yMax;
	__int32 ABM_m_xOffset;
	__int32 ABM_m_yOffset;
	__int32 ABM_m_zOffset;
	__int32 ABM_Mode;
	__int32 ABM_OsziDarkInvert;
	__int32 ABM_ErrorHisto;
	__int32 ABM_XShift;
	__int32 ABM_YShift;
	__int32 ABM_ZShift;
	__int32 ABM_ozShift;
	__int32 ABM_wdShift;
	__int32 ABM_ucLevelXY;
	__int32 ABM_ucLevelZ;
	__int32 ABM_uiABMXShift;
	__int32 ABM_uiABMYShift;
	__int32 ABM_uiABMZShift;
	bool use_normal_method;

	__int32 TWOHM1_FAK_DLL_Value;
	__int32 TWOHM1_Resolution_Flag;
	__int32 TWOHM1_trigger_mode_for_start;
	__int32 TWOHM1_trigger_mode_for_stop;
	__int32 TWOHM1_res_adjust;
	double TWOHM1_tdcresolution;
	__int32 TWOHM1_test_overflow;
	__int32 TWOHM1_number_of_channels;
	__int32 TWOHM1_number_of_hits;
	__int32 TWOHM1_set_bits_for_GP1;
	__int32 TWOHM1_HM1_ID_1;
	__int32 TWOHM1_HM1_ID_2;
};


















class LMF_IO
{
public:

	LMF_IO(__int32 Number_of_Channels, __int32 Number_of_Hits);
    ~LMF_IO();

	static __int32			GetVersionNumber();
	bool			OpenInputLMF(__int8* Filename);
	bool			OpenInputLMF(std::string Filename);
	void			CloseInputLMF();

	bool			OpenOutputLMF(__int8* Filename);
	bool			OpenOutputLMF(std::string Filename);
	void			CloseOutputLMF();

	unsigned __int64 GetLastLevelInfo();

	bool			Clone(LMF_IO * clone);

	void			WriteTDCData(double timestamp,unsigned __int32 cnt[],__int32 *tdc);
	void			WriteTDCData(unsigned __int64 timestamp,unsigned __int32 cnt[],__int32 *tdc);
	void			WriteTDCData(double timestamp,unsigned __int32 cnt[],double *tdc);
	void			WriteTDCData(unsigned __int64 timestamp,unsigned __int32 cnt[],double *tdc);
	void			WriteTDCData(double timestamp,unsigned __int32 cnt[],__int64 *tdc);
	void			WriteTDCData(unsigned __int64 timestamp,unsigned __int32 cnt[],__int64 *tdc);
	void			WriteTDCData(double timestamp,unsigned __int32 cnt[],unsigned __int16 *tdc);
	void			WriteTDCData(unsigned __int64 timestamp,unsigned __int32 cnt[],unsigned __int16 *tdc);

	bool			ReadNextfADC4packet(ndigo_packet * packet, bool &bEnd_of_group_detected, __int16 * i16buffer, __int32 buffersize_in_16bit_words);
	
	bool			ReadNextfADC8Signal(fADC8_signal_info_struct& signal_info, bool& bEnd_of_group_detected, unsigned __int32 * ui32buffer, __int32 buffersize_in_32bit_words);
	bool			WriteNextfADC8Signal(fADC8_signal_info_struct* signal_info, unsigned __int32 * ui32buffer);
	bool			WritefADC8EndGroupMarker();

	void			WriteFirstHeader();
	bool			ReadNextEvent();
	bool			ReadNextCAMACEvent();
	void			GetNumberOfHitsArray(unsigned __int32 cnt[]);
	void			GetNumberOfHitsArray(__int32 cnt[]);
	void			GetTDCDataArray(__int32 *tdc);
	void			GetTDCDataArray(__int64 *tdc);
	void			GetTDCDataArray(double *tdc);
	void			GetTDCDataArray(unsigned __int16 * tdc);

	void			GetCAMACArray(unsigned __int32 []);
	void			WriteCAMACArray(double, unsigned int[]);

	unsigned __int64		GetEventNumber();
	double			GetDoubleTimeStamp();
	unsigned __int64 Getuint64TimeStamp();
	unsigned __int32	GetNumberOfChannels();
	unsigned __int32	GetMaxNumberOfHits(); 
	bool			SeekToEventNumber(unsigned __int64 Eventnumber);

	const char *	GetErrorText(__int32 error_id);
	void			GetErrorText(__int32 error_id, __int8 char_buffer[]);
	void			GetErrorText(__int8 char_buffer[]);

	__int32				GetErrorStatus();


	bool decompress_synchronous_fADC8_signal(fADC8_signal_info_struct &signal_info, unsigned __int32 source_buffer[], unsigned __int32 source_buffer_size_in_32bit_words, __int32 &number_of_filled_16bit_words,
												__int16 i16bit_target_buffer1[], __int32 target_buffer1_size_in_16bit_words,
												__int16 i16bit_target_buffer2[], __int32 target_buffer2_size_in_16bit_words,
												__int16 i16bit_target_buffer3[], __int32 target_buffer3_size_in_16bit_words,
												__int16 i16bit_target_buffer4[], __int32 target_buffer4_size_in_16bit_words);

	bool decompress_asynchronous_fADC8_signal(fADC8_signal_info_struct &signal_info, unsigned __int32 source_buffer[], unsigned __int32 source_buffer_size_in_32bit_words, __int16 i16bit_target_buffer[], __int32 target_buffer_size_in_16bit_words, __int32 &number_of_filled_16bit_words);

	void prepare_Cobold2008b_TDC8HP_header_output() { // Cobold2008 release August 2009
		//TDC8HP.TriggerDeadTime_p68;	
		//TDC8HP.GroupRangeStart_p69;	
		//TDC8HP.GroupRangeEnd_p70;
		//Starttime_output;
		//Stoptime_output;
		TDC8HP.FallingEnable_p62 = 1044991; // all 9 channels on 2 cards
		TDC8HP.TriggerChannel_p64 = 1;
		frequency = 1.e12; // 1ps
		tdcresolution_output = 0.025; // or 0.016 ?
		TDC8HP.VHR_25ps = true;
		TDC8HP.UserHeaderVersion = 7;
		DAQVersion_output = 20080507;
		Cobold_Header_version_output = 2008;
		system_timeout_output = 5; // obsolete
		time_reference_output = (__int32)int(Starttime_output);
		common_mode_output = 0; // 0 = common start
		data_format_in_userheader_output = -1;
		DAQ_ID_output = 8; // 8 = TDC8HP
		timestamp_format_output = 2;
		LMF_Version_output = 10;
		IOaddress = 0;
		number_of_DAQ_source_strings_output = 0;
		TDC8HP.OutputRollOvers_p72 = 1;
		TDC8HP.Number_of_TDCs = 0;
		TDC8HP.number_of_int32s = 5;
		TDC8HP.number_of_bools = 0;
		TDC8HP.number_of_doubles = 1;
		TDC8HP.GroupingEnable_p66_output = false;
	}

private:
	void			Initialize();

	void			write_times(MyFILE *,time_t,time_t);

	bool			OpenNonCoboldFile(void);

	__int32			ReadCAMACHeader();
	__int32			ReadTDC8PCI2Header();
	__int32			Read2TDC8PCI2Header();
	__int32			ReadTDC8HPHeader_LMFV_1_to_7(__int32 byte_counter_external);
	__int32			ReadTDC8HPHeader_LMFV_8_to_9();
	__int32			ReadTDC8HPHeader_LMFV_10();
	__int32			ReadHM1Header();
	__int32			ReadTCPIPHeader();
	__int32			ReadfADC8Header();
	__int32			WritefADC8_header_LMFversion10();
	__int32			ReadfADC8_header_LMFversion10();
	__int32			ReadfADC8_header_LMFversion9();

	__int32			ReadfADC4Header_up_to_v11();

	void			WriteEventHeader(unsigned __int64 timestamp, unsigned __int32 cnt[]);
	__int32			WriteCAMACHeader();
	__int32			WriteTDC8PCI2Header();
	__int32			Write2TDC8PCI2Header();
	__int32			WriteTDC8HPHeader_LMFV_1_to_7();
	__int32			WriteTDC8HPHeader_LMFV_8_to_9();
	__int32			WriteTDC8HPHeader_LMFV_10_to_12();
	__int32			WriteHM1Header();
	__int32			WriteTCPIPHeader();
	bool			Read_TDC8HP_raw_format(unsigned __int64 &ui64TDC8HP_AbsoluteTimeStamp);
	__int32			PCIGetTDC_TDC8HP_25psGroupMode(unsigned __int64 &ui64TDC8HPAbsoluteTimeStamp, __int32 count, unsigned __int32 * Buffer);


public:
	std::string			Versionstring;
	std::string			FilePathName;
	std::string			OutputFilePathName;
	std::string			Comment;
	std::string			Comment_output;
	std::string			DAQ_info;
	std::string			Camac_CIF;

	char *error_text[40];

	double *Parameter;
	double *Parameter_old;

	time_t			Starttime;
	time_t			Stoptime;
	time_t			Starttime_output;
	time_t			Stoptime_output;

	__int32				time_reference;
	__int32				time_reference_output;
	
	unsigned __int32	ArchiveFlag;
	__int32				Cobold_Header_version;
	__int32				Cobold_Header_version_output;

	unsigned __int64	uint64_LMF_EventCounter;
	unsigned __int64	uint64_number_of_read_events;
	unsigned __int64	uint64_Numberofevents;

	__int32				Numberofcoordinates;
	__int32				CTime_version,CTime_version_output;
	unsigned __int32	SIMPLE_DAQ_ID_Orignial;
	__int32	DAQVersion;
	__int32	DAQVersion_output;
	unsigned __int32	DAQ_ID;
	unsigned __int32	DAQ_ID_output;
	__int32				data_format_in_userheader;
	__int32				data_format_in_userheader_output;

	unsigned __int32	Headersize;
	unsigned __int32	User_header_size;
	unsigned __int32	User_header_size_output;

	__int32				IOaddress;
	__int32				timestamp_format;
	__int32				timestamp_format_output;
	__int32				timerange;

	unsigned __int32	number_of_channels;
	unsigned __int32	number_of_channels2;
	unsigned __int32	max_number_of_hits;
	unsigned __int32	max_number_of_hits2;

	__int32				number_of_channels_output;
	__int32				number_of_channels2_output;
	__int32				max_number_of_hits_output;
	__int32				max_number_of_hits2_output;

	unsigned __int64	ui64LevelInfo;

	unsigned __int32	changed_mask_read;

	unsigned __int8		ui8_PostEventData[MAX_NUMBER_OF_BYTES_IN_POSTEVENTDATA];
	__int32				number_of_bytes_in_PostEventData;

	__int32				DAQSubVersion;
	__int32				module_2nd;
	__int32				system_timeout;
	__int32				system_timeout_output;
	__int32				common_mode;
	__int32				common_mode_output;
	__int32				DAQ_info_Length;
	__int32				Camac_CIF_Length;
	__int32				LMF_Version;
	__int32				LMF_Version_output;
	__int32				TDCDataType;
	
	__int32				iLMFcompression;
	unsigned __int32	LMF_Header_version;

	double			tdcresolution;
	double			tdcresolution_output;
	double			frequency;
	double			DOUBLE_timestamp;
	unsigned __int64	ui64_timestamp;

	__int32				errorflag;
	bool			skip_header;

	unsigned __int32 DAQ_SOURCE_CODE_bitmasked;
	unsigned __int32 DAN_SOURCE_CODE_bitmasked;
	unsigned __int32 CCF_HISTORY_CODE_bitmasked;

	__int32	number_of_CCFHistory_strings;
	__int32	number_of_CCFHistory_strings_output;
	std::string ** CCFHistory_strings;
	std::string ** CCFHistory_strings_output;
	__int32	number_of_DAN_source_strings;
	__int32	number_of_DAN_source_strings_output;
	std::string ** DAN_source_strings;
	std::string ** DAN_source_strings_output;
	__int32	number_of_DAQ_source_strings;
	__int32 number_of_DAQ_source_strings_output;
	std::string ** DAQ_source_strings;
	std::string ** DAQ_source_strings_output;

	HM1_struct		HM1;
	TDC8HP_struct	TDC8HP;
	TDC8PCI2_struct TDC8PCI2;
	fADC8_struct	fADC8;
	fADC4_struct	fADC4;

	unsigned __int64 uint64_number_of_written_events;

	MyFILE *	input_lmf;
	MyFILE *	output_lmf;

	bool			InputFileIsOpen;
	bool			OutputFileIsOpen;

private:
	unsigned __int32 *	ui32buffer;
	__int32				ui32buffer_size;
	bool				not_Cobold_LMF;
	unsigned __int32	Headersize_output;
	__int32				output_byte_counter;
	__int32				Numberofcoordinates_output;
	bool				must_read_first;
	unsigned __int32	* number_of_hits;
	__int32				*i32TDC;
	unsigned __int16	*us16TDC;
	double				*dTDC;

	__int32				num_channels;
	__int32				num_ions;


	unsigned __int32	* CAMAC_Data;
};
#endif �1  P   I D S _ S O U R C E _ C P P   A N A L Y S I S _ C P P       0         #pragma warning(disable : 4800)
#include "OS_Version.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtupleD.h"

#include <math.h>
#include "rootstuff.h"
#include "Histo.h"
#include "TF1.h"
//#include "TMinuit.h"

#include "functions.h"
#include "Ueberstruct.h"

//#include "polynomial_TOF_to_P.hpp"
//#include "resort64c.h"

#define AUTO __LINE__+1000 //returns the line number in the current file. Is used to make life easier for histogram numbering.

char ntuple_identifier_A[2048];
bool create_ntuple_identifier_A;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int analysis(__int64 eventcounter, double  parameter[], TTree * Data, Ueberstruct * Ueber)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	
	// do some renaming for less typing...
	Histo * Hist = Ueber->Hist;

	Ueber->start_new_root_file = false;

	// variable used when avoiding ColAHelL
	double reaction=-1., ehit=-1., rhit=-1., phit=-1., bunchmarker=-1., timestamp =-1., LMFStart=-1., LMFStop=-1.;
	double ex[64],ey[64],etof[64],etime[64],emethod[64];
	double rx[64],ry[64],rtof[64],rtime[64],rmethod[64];
	double px[64],py[64],ptof[64],ptime[64],pmethod[64];
	double scan_val[16];

	double pdet_size	= 5. + parameter[122];				// size of projectile detector (for displaying spectra properly)
	double rdet_size	= 5. + parameter[222];				// size of recoil detector (for displaying spectra properly)
	double edet_size	= 5. + parameter[322];				// size of electron detector (for displaying spectra properly)

	double NTupleData[250];
	bool WriteNTuple = false;
	char tmp[200];;

	int NRmax = (int)(parameter[61]+0.01);
	int NEmax = (int)(parameter[62]+0.01);
	int NPmax = (int)(parameter[60]+0.01);

	if (Ueber->error) {
		Ueber->stop_reading_files = true;
		return 0;
	}

	if(eventcounter == 0) {
		Ueber->EntriesInFile = 0;
		Ueber->eventswritten = 0;
		create_ntuple_identifier_A = true;
	}


	if(Ueber->EntriesInFile == 0) {

		int MaxRec = 0;
		int MaxElec = 0;
		int MaxProj = 0;

		char ntuple_identifier[500];
		bool Check = 0;
		// set branch addresses...

		Data->SetBranchAddress("reaction",&reaction);
		Data->SetBranchAddress("ehit",&ehit);
		Data->SetBranchAddress("rhit",&rhit);
		Data->SetBranchAddress("phit",&phit);
		Data->SetBranchAddress("bunchmarker",&bunchmarker);
		Data->SetBranchAddress("timestamp",&timestamp);
		Data->SetBranchAddress("LMFStart",&LMFStart);
		Data->SetBranchAddress("LMFStop",&LMFStop);


		// find out if/how many scan parameters are in Ntuple
		for(int i=15;i>=0;i--) {
			sprintf(ntuple_identifier,"scanval%i",i);
			if(Data->GetBranch(ntuple_identifier)) {
				Ueber->scan_in_data = true;
				Ueber->scan_num_vals = i+1;
				break;
			}
		}

		//check whether we have a scan value in the ntuple
		if(Ueber->scan_in_data) {
			for(int i=0;i<Ueber->scan_num_vals;i++) {
				sprintf(ntuple_identifier,"scanval%i",i);
				Data->SetBranchAddress(ntuple_identifier,&scan_val[i]);
			}
		}

		// find out how many recoils are in Ntuple
		do {
			sprintf(ntuple_identifier,"r%ix",++MaxRec);
		} while(Data->GetBranch(ntuple_identifier));
		MaxRec = (MaxRec<17 ? MaxRec : 16);
		Ueber->RecInNTuple = MaxRec;

		for(int i=0;i<MaxRec-1;i++) {
			sprintf(ntuple_identifier,"r%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rx[i]);
			sprintf(ntuple_identifier,"r%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ry[i]);
			sprintf(ntuple_identifier,"r%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rtime[i]);
			sprintf(ntuple_identifier,"r%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rtof[i]);
			sprintf(ntuple_identifier,"r%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&rmethod[i]);
		}

		// find out how many electrons are in Ntuple
		do {
			sprintf(ntuple_identifier,"e%ix",++MaxElec);
		} while(Data->GetBranch(ntuple_identifier));
		MaxElec = (MaxElec<17 ? MaxElec : 16);
		Ueber->ElecInNTuple = MaxElec;

		for(int i=0;i<MaxElec-1;i++) {
			sprintf(ntuple_identifier,"e%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ex[i]);
			sprintf(ntuple_identifier,"e%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ey[i]);
			sprintf(ntuple_identifier,"e%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&etime[i]);
			sprintf(ntuple_identifier,"e%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&etof[i]);
			sprintf(ntuple_identifier,"e%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&emethod[i]);
		}

		// find out how many projectyles are in Ntuple
		do {
			sprintf(ntuple_identifier,"p%ix",++MaxProj);
		} while(Data->GetBranch(ntuple_identifier));
		MaxProj = (MaxProj<17 ? MaxProj : 16);
		Ueber->ProjInNTuple = MaxProj;
		
		for(int i=0;i<MaxProj-1;i++) {
			sprintf(ntuple_identifier,"p%ix",i+1);
			Data->SetBranchAddress(ntuple_identifier,&px[i]);
			sprintf(ntuple_identifier,"p%iy",i+1);
			Data->SetBranchAddress(ntuple_identifier,&py[i]);
			sprintf(ntuple_identifier,"p%imcp",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ptime[i]);
			sprintf(ntuple_identifier,"p%itof",i+1);
			Data->SetBranchAddress(ntuple_identifier,&ptof[i]);
			sprintf(ntuple_identifier,"p%iflag",i+1);
			Data->SetBranchAddress(ntuple_identifier,&pmethod[i]);
		}
	}

	Data->GetEntry(Ueber->EntriesInFile);

	if(Ueber->EntriesInFile < Data->GetEntries()-1) {
		++Ueber->EntriesInFile;
	} else {
		Ueber->EntriesInFile = 0;
	}

	CH_event_struct * CH_evt;
	if(Ueber->use_CH) {
		double ds = -1.0;
		ds =  Ueber->CH->GetPDetSize();
		if(ds>0.0)
			pdet_size = ds + 5.0;
		ds =  Ueber->CH->GetIDetSize();
		if(ds>0.0)
			rdet_size = ds + 5.0;
		ds =  Ueber->CH->GetEDetSize();
		if(ds>0.0)
			edet_size = ds + 5.0;

		//copy data to ColAHelL
		Ueber->CH->ResetEvent();
		Ueber->CH->SetEvent_reaction(reaction);
		Ueber->CH->SetEvent_bm(bunchmarker);
		for(int i=0;i<ehit;i++)
			Ueber->CH->SetEvent_xyttof(EL,ex[i],ey[i],etime[i],etof[i],emethod[i]);
		for(int i=0;i<rhit;i++)
			Ueber->CH->SetEvent_xyttof(IO,rx[i],ry[i],rtime[i],rtof[i],rmethod[i]);
		for(int i=0;i<phit;i++)
			Ueber->CH->SetEvent_xyttof(PR,px[i],py[i],ptime[i],ptof[i],pmethod[i]);
		for(int i=0;i<Ueber->scan_num_vals;i++)
			Ueber->CH->SetEvent_scanvals(scan_val[i]);

		CH_evt = Ueber->CH->GetEvent();
// IPA needs food... Feeding is done here...
// -----------------------------------------------------------
		if(Ueber->ipa_enable) {
			Ueber->IPA->append_data(CH_evt);
			return 0;
		}

// Process the data and plot histograms. We might have several different
// reaction definitions belonging to the same "channel" (e.g. in case
// the "ANY" presorter was used), so we might need to loop over event 
// several times...

		// Let's see if we find a reaction definition for the current channel..
		// -----------------------------------------------------------	

		int EvtBelongsToReactions = Ueber->CH->EvtBelongsToReactions((int)CH_evt->reaction);
		// We need to loop, as there might be several reactions using the same channel number...
		// (e.g. when people use the ANY presorter...)
		for(int i=0;i<EvtBelongsToReactions;i++) {
//			char reaction_dir[80];
//			strcpy(reaction_dir,Ueber->Proc->reaction_list[cur_reac]->name);

			// ColAHelL kicks in here! ...
			// -----------------------------------------------------------
			// process particles...
			Ueber->CH->ProcessEvent(i);
			// fill std. ColAHelL histograms
			Ueber->CH->FillHistograms();
//------------------------------------------------------------------------------------------------
//
// Include your analysis here.
//
//	you may use:
//			electrons: e[n]->...
//			ions: r[n]->...
//			diatomic molecule: mol->...
//
//	information on the reaction occured for this event (if identified, i.e. found_current_reaction == true) is in: 
//			cur_reaction->...
//
//	a char string with the name of the current reaction is in:
//			reaction_dir
//
//  "raw" event data is located in CH_evt.
//
//	"AUTO" can now be used instead of histogram numbers
//	Hist->fill(AUTO, "ehit_ntuple", ehit, 1., "Electron hits (ntuple)", 43, -1.25, 20.25, "ehit", "ntuple");
//
//------------------------------------------------------------------------------------------------

//	__int64 TotalNumEvents = (__int64)(Data->GetEntries()); //total number of events in the root file that is open


		
	



		

//------------------------------------------------------------------------------------------------


			if(WriteNTuple) {
				int NRec = (int)rhit;
				int NElec = (int)ehit;
				int NPro = (int)phit;
		
				if(NRec > NRmax)
					NRec = NRmax;
				if(NElec > NEmax)
					NElec = NEmax;
				if(NPro > NPmax)
					NPro = NPmax;

				NTupleData[0] = reaction;
				NTupleData[1] = double(NElec);
				NTupleData[2] = double(NRec);
				NTupleData[3] = double(NPro);
				NTupleData[4] = bunchmarker;
				NTupleData[5] = timestamp;
				NTupleData[6] = (double)Ueber->LMF_input->Starttime;
				NTupleData[7] = (double)Ueber->LMF_input->Stoptime; 

				int absind = 8;		

				if(Ueber->scan_in_data) {
					for(int i=0;i<Ueber->scan_num_vals;i++)
					{
						NTupleData[absind+i]=scan_val[i];
					}
					absind += Ueber->scan_num_vals;
				}

				for(int i=0;i<NRec;i++)
				{	// begin for: write all Recoils to NTuple
						NTupleData[absind+i*5] = CH_evt->r.x[i];
						NTupleData[absind+i*5+1] = CH_evt->r.y[i];
						NTupleData[absind+i*5+2] = CH_evt->r.time[i];
						NTupleData[absind+i*5+3] = CH_evt->r.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->r.method[i];
				}	// end for: write all Recoils to NTuple

				absind += NRmax*5;

				for(int i=0;i<NElec;i++)
				{	// begin for: write all Electrons to NTuple
						NTupleData[absind+i*5] = CH_evt->e.x[i];
						NTupleData[absind+i*5+1] = CH_evt->e.y[i];
						NTupleData[absind+i*5+2] = CH_evt->e.time[i];
						NTupleData[absind+i*5+3] = CH_evt->e.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->e.method[i];
				}	// end for: write all Electrons to NTuple

				absind += NEmax*5;

				for(int i=0;i<NPro;i++)
				{	// begin for: write all Projectiles to NTuple
						NTupleData[absind+i*5] = CH_evt->p.x[i];
						NTupleData[absind+i*5+1] = CH_evt->p.y[i];
						NTupleData[absind+i*5+2] = CH_evt->p.time[i];
						NTupleData[absind+i*5+3] = CH_evt->p.tof[i];
						NTupleData[absind+i*5+4] = CH_evt->p.method[i];
				}	// end for: write all Projectiles to NTuple

			//////////////////////////////////////////////////////////////////////////
			//	write data to ntuple
			//////////////////////////////////////////////////////////////////////////

				if(create_ntuple_identifier_A) {
					create_ntuple_identifier_A = false;
			
					sprintf(ntuple_identifier_A,"reaction:ehit:rhit:phit:bunchmarker:timestamp:LMFStart:LMFStop");

					if(Ueber->scan_in_data) {

					for(int i=0;i<Ueber->scan_num_vals;i++)
					{
						sprintf(tmp,":scanval%d",i);
						strcat(ntuple_identifier_A,tmp);
					}
				}

				for(int i=1;i<=NRmax;i++)
				{
					sprintf(tmp,":r%dx:r%dy:r%dmcp:r%dtof:r%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
				}
				for(int i=1;i<=NEmax;i++)
				{
					sprintf(tmp,":e%dx:e%dy:e%dmcp:e%dtof:e%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
				}
				for(int i=1;i<=NPmax;i++)
				{
					sprintf(tmp,":p%dx:p%dy:p%dmcp:p%dtof:p%dflag",i,i,i,i,i);
					strcat(ntuple_identifier_A,tmp);
					}
				}

				Hist->NTupleD(0,"Data","BESSY2010",ntuple_identifier_A, 32000, NTupleData);
				++Ueber->eventswritten;

				if(parameter[57]>0.5) {
					unsigned __int64 max_events = (unsigned __int64)(parameter[56]+0.1);
					if(Ueber->eventswritten > (__int64)max_events && max_events > 0) {
						Ueber->start_new_root_file = true;
						Ueber->eventswritten = 0;
						Hist->Reset();
					}
				}
			}
		}
	} //end loop over event

// A few very basic histograms
	Hist->fill(AUTO,"reaction_ntuple",reaction,1.,"Reaction channel flag (ntuple)",203,-1.25,100.25,"reaction_flag","ntuple");
	Hist->fill(AUTO,"rhit_ntuple",rhit,1.,"Recoil hits (ntuple)",43,-1.25,20.25,"rhit","ntuple");
	Hist->fill(AUTO,"phit_ntuple",phit,1.,"Projectile hits (ntuple)",43,-1.25,20.25,"phit","ntuple");
	Hist->fill(AUTO,"ehit_ntuple",ehit,1.,"Electron hits (ntuple)",43,-1.25,20.25,"ehit","ntuple");

	return 0;
} �Q  T   I D S _ S O U R C E _ C P P   F U N C T I O N S _ C P P         0         #include "OS_Version.h"
#include <math.h>

#include "Ueberstruct.h"

#include "TCanvas.h"    //Needed for Tof_to_P
#include "TGraph.h"		//Needed for Tof_to_P
#include "TF1.h"		//Needed for Tof_to_P
//#include <vector>
//#include <string>
//#include <iostream>
//#include <sstream>
//
//std::vector<std::vector<double>> read_array_file(std::string filename){
//	std::vector<std::vector<double>> results;
//	
//	ifstream array_file(filename);
//	std::string line;
//	std::string comment = "//";
//	std::string cur_open = "{";
//	std::string cur_close = "}";
//	std::string equ = "=";
//	std::string comma = ",";
//	std::string semi = ";";
//	std::string open_comment = "/*";
//	std::string close_comment = "*/";
//	std::string space = " ";
//	std::string tab = "\t";
//	std::string tofs = "tofs";
//	std::string mons = "mons";
//
//	size_t found_tofs;
//	size_t found_mons;
//	size_t found_comment;
//	size_t found_semi;
//	size_t found_cur_open;
//	size_t found_cur_close;
//	size_t found_space;
//	size_t found_tab;
//	size_t found_Open_comment;
//	size_t found_Close_comment;
//	size_t found_comma;
//
//	bool multi_line_comment = false;
//	bool inside_tofs = false;
//	bool inside_mons = false;
//
//	if (array_file.is_open())
//	{
//		while (!array_file.eof())
//		{
//			getline(array_file, line);
//
//			if (multi_line_comment) {
//				found_Close_comment = line.find(close_comment);
//				if (found_Close_comment != std::string::npos) {
//					multi_line_comment = false;
//					line = line.substr(found_Close_comment + 2);
//					//cout << "found close comment. multi_line_comment= " << multi_line_comment<< endl;
//				}
//			}
//			if (multi_line_comment == false) {
//				found_Open_comment = line.find(open_comment);
//				if (found_Open_comment != std::string::npos) {
//					line = line.substr(0, found_Open_comment);
//					multi_line_comment = true;
//					//cout << "found open comment. multi_line_comment= " << multi_line_comment<< endl;
//				}
//
//				// remove the comments
//				found_comment = line.find(comment);
//				if (found_comment != std::string::npos)
//					line = line.substr(0, found_comment);
//
//				// replace the tabs with spaces
//				found_tab = line.find(tab);
//				while (found_tab != std::string::npos) {
//					line = line.replace(found_tab, 1, " ");
//					found_tab = line.find(tab);
//				}
//
//				// replace the { with spaces
//				found_cur_open = line.find(cur_open);
//				while (found_cur_open != std::string::npos) {
//					line = line.replace(found_cur_open, 1, " ");
//					found_cur_open = line.find(cur_open);
//				}
//
//				// replace the } with spaces
//				found_cur_close = line.find(cur_close);
//				while (found_cur_close != std::string::npos) {
//					line = line.replace(found_cur_close, 1, " ");
//					found_cur_close = line.find(cur_close);
//				}
//
//				// replace the comma with spaces
//				found_comma = line.find(cur_close);
//				while (found_comma != std::string::npos) {
//					line = line.replace(found_comma, 1, " ");
//					found_comma = line.find(found_comma);
//				}
//
//
//				found_tofs = line.find(tofs);
//				if (found_tofs != std::string::npos) {
//					inside_tofs = true;
//					line = line.erase(found_tofs, 4);
//				}
//				
//				if (inside_tofs) {
//					
//					found_semi	= line.find(semi);
//					if (found_semi != std::string::npos) {
//						line = line.replace(found_semi, 1, " ");
//						inside_tofs = false;
//					}
//
//					std::stringstream line_stream(line);
//					while (!line_stream.eof()) {
//						double temp;
//						line_stream >> temp;
//						results[0].push_back(temp);
//					}
//				}
//
//
//				found_mons = line.find(mons);
//				if (found_mons != std::string::npos) {
//					inside_mons = true;
//					line = line.erase(found_mons, 4);
//				}
//
//				if (inside_mons) {
//
//					found_semi = line.find(semi);
//					if (found_semi != std::string::npos) {
//						line = line.replace(found_semi, 1, " ");
//						inside_mons = false;
//					}
//
//					std::stringstream line_stream(line);
//					while (!line_stream.eof()) {
//						double temp;
//						line_stream >> temp;
//						results[1].push_back(temp);
//					}
//				}
//			
//			
//			}
//
//		}
//	}
//
//	std::cout << "\nTofs = ";
//	for (int i = 0; i > results[0].size(); i++) {
//		std::cout << results[0][i] << "  ";
//	}
//
//	std::cout << "\nMons = ";
//	for (int i = 0; i > results[1].size(); i++) {
//		std::cout << results[1][i] << "  ";
//	}
//
//	return results;
//}


//Takes at set of tof and momentum points as imput and returns a TF1 fuction that is fitted to these points.
//This TF1 can then be used to convert tof to momentum.  The funciton displays the results of the fit in a window while LMF2Root is running.
//It uses a polynomial function A*x^4+B*x^3+C*x^2+D*x+E, where A,B,C,D,E are fitted parameters.  
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Example (this code is placed in Analysis.cpp):
//		static TF1 * H_plus;
//		if (eventcounter == 0) { //this code runs only once!!
//		
//			double tofs[] = { 1602.28,	1600.11,	1597.91,	1595.69,	1181.80,	1179.73 }; 
//			double moms[] = { -38.51,	-38.12,		-37.73,		-37.33,		38.12,		38.51	};
//		
//			H_plus = Tof_to_P("H_plus_", 5, tofs, moms);
//		}
//////////////////////////////////////////////////////////////////////////////////////////////////////
// to get the value from the fitted function: H_plus->Eval(CH_evt->r.tof[0])
//
// By Joshua Williams
TF1 *  Tof_to_P(char * name, int num_data_points, double *tofs, double *moms)
{
	TCanvas *c1 = new TCanvas(name, name, 200, 10, 700, 500);
	c1->SetGrid();
	gPad->SetRightMargin(Float_t(0.12));
	gPad->SetLeftMargin(Float_t(0.11));
	gPad->SetTopMargin(Float_t(0.0825));
	gPad->SetBottomMargin(Float_t(0.087));

	TGraph * gr1 = new TGraph(num_data_points, tofs, moms);
	gr1->SetTitle(name);
	gr1->SetMarkerColor(kBlue);
	gr1->SetMarkerStyle(21);
	gr1->SetLineWidth(0);
	gr1->SetLineColor(0);
	gr1->Draw();

	double uxmin = gPad->GetUxmin();
	double uxmax = gPad->GetUxmax();

	char fitname[256];
	sprintf(fitname, "%s-myfit", name);

	TF1 * myfit = new TF1(fitname, "[4]*x**4+[3]*x**3+[2]*x**2+[1]*x+[0]", uxmin, uxmax);
	myfit->SetParName(0, "E");
	myfit->SetParName(1, "D");
	myfit->SetParName(2, "C");
	myfit->SetParName(3, "B");
	myfit->SetParName(4, "A");
	myfit->SetParameter(0, 0);
	myfit->SetParameter(1, 0);
	myfit->SetParameter(2, 0);
	myfit->SetParameter(3, 0);
	myfit->SetParameter(4, 0);
	std::cout << "Fit parameters for " << name << std::endl;
	gr1->Fit(fitname, "M");


	gPad->Modified();
	gPad->Update();
	gPad->Modified();
	gPad->Update();

	//printf("0%           25%          50%          75%        100%\r\n");
	//printf("|------------|------------|------------|-----------|\r\n");

	return myfit;
}

// **************************************************************
// **
// ** momentum calculation from t_0 (2:1 time-focusing geometry)
// **
// **************************************************************

double tof2mom21(double tof_ns, double acc_mm, double drift_mm, double Efeld_Vpcm, double mass_au, double charge_au, double  parameter[])
{
	double s1 = acc_mm * 1e-3;
	double s2 = drift_mm * 1e-3;
	double t = tof_ns * 1e-9;
	double q = charge_au * 1.6022e-19;
	double mass = mass_au * 9.1095e-31;
	
	double a = Efeld_Vpcm*100. * q/mass;

	double v0 = -t*a + sqrt(2.*a*s1) + s2*sqrt(a/2.* s1);
	for (__int32 i=0;i<10;++i)
	{
		double test1Tof = -v0/a + sqrt(v0*v0 + 2.*a*s1)/a + s2/sqrt(v0*v0 + 2.*a*s1);
		double v1 = 1.01*v0;
		double test2Tof = -v1/a + sqrt(v1*v1 + 2.*a*s1)/a + s2/sqrt(v1*v1 + 2.*a*s1);
		double dtdvo = (test2Tof - test1Tof) / (v1-v0);
		v0 = v0 + 0.7*(t - test1Tof)/dtdvo;
	}

	double mom = v0 * mass/ (9.1095e-31*2.1877e6);
	return mom;
}


// **************************************************************
// **
// **  momentum calculation from t_0 (homogeneous acceleration)
// **
// **************************************************************

double tof2mom(double tof_ns, double acc_mm, double Efield_Vpcm, double mass_au, double charge_au, double parameter[])
{
	double s1 = acc_mm * 1e-3;
	double t = tof_ns * 1e-9;
	double vau = 2.1877e+6;                              // unity velocity in atomic units [m/s]
	double mom;

	mom = (s1 / t - 0.5 * 1.7588e13*Efield_Vpcm * charge_au / mass_au * t ) / vau * mass_au;

	return mom;
}
	



// **************************************************************
// **
// ** calculate momenta for tof-direction /w McLaren geometry
// **
// **************************************************************
double toftomom(double tof_ns, double acc_mm, double drift_mm, double Efeld_Vpcm, double mass_au, double charge_au, double  parameter[]) {

// *******************************************
// this is lothars subroutine
// ********************************************

	double s1 = acc_mm * 1e-3;
	double s2 = drift_mm* 1e-3;

	double a = 1.7588E13 * parameter[1070];	// electrons

	double etof = tof_ns * 1.e-9;


	double vo = - etof*a + sqrt(2.*a*s1) + s2*sqrt(a/2./s1);

	double dt = 10.;               

	for(__int32 i=1;i<7;i++) {
		
		if((vo*vo + 2.*a*s1) > 1e-9) {
    
			double test1tof = -vo/a +  sqrt(vo*vo + 2.*a*s1 )/a + s2/sqrt(vo*vo + 2.*a*s1);
			double v1 = 1.01 * vo;
			double test2tof = -v1/a + sqrt(v1*v1 + 2.*a*s1 )/a + s2/sqrt(v1*v1 + 2.*a*s1);
			double dtdvo = (test2tof - test1tof) / (v1 - vo);
			vo = vo + 0.7 * (etof - test1tof) / dtdvo;
		}
	}
	
	return vo/2.18E6;
}




	
// **************************************************************
// **
// ** position correction for electron (E-B-drift)
// **
// **************************************************************

double elec_pos_corr(double e1x_mm, double e1y_mm, double e1t_ns, char direction, double  parameter[])
{
	double outputvalue=754653.;

	
	// position correction for different times (E-B-drift)
	double e1x_corr		= -1. * ((e1x_mm - ( parameter[355] + (parameter[360] - parameter[355]) * (e1t_ns - parameter[357]) / (parameter[362] - parameter[357]) ) ));
	double e1y_corr		= (e1y_mm - ( parameter[356] + (parameter[361] - parameter[356]) * (e1t_ns - parameter[357]) / (parameter[362] - parameter[357]) ) );
	
	if (direction=='x') {outputvalue=e1x_corr;}
	else {if (direction=='y') {outputvalue=e1y_corr;} else {outputvalue=85734.;} }

	return outputvalue;
}

	
// **************************************************************
// **
// ** calculate momenta for electron (magnetic field)
// **
// **************************************************************


double elec2mom(double e1x_mm, double e1y_mm, double e1t_ns, char direction, double  parameter[])
{
	double e_beta, epx1, epy1;
	double pi = 3.14159265359;
	double outputvalue=754653.;

	// parameter 1206 = width how much around a node is thrown away [cos]
	// parameter 1209 = gyration period [ns]
	// parameter 1213 = additional detector rotation [deg]

	// magnetic field, nodes etc.
	double b_field		= 357.238755341 / fabs (parameter[1081]);										// magnetic field in Gauss from wiggle [ns]
	double e_alpha		= fmod((( (e1t_ns - parameter[1081]) * 360. / parameter[1081]) + 1080.), 360.);	// angle in gyration phase
	double rad_callib	= sqrt(2. - 2. * cos(e_alpha * pi / 180) );

	if (rad_callib >= parameter[1206])
		{						// Begin IF: not in the node
		e_beta		= -1 * fmod( ((e_alpha/2.) +360.), 180. ); // + parameter[1213];
			rad_callib	= 1/rad_callib * b_field / 124.38406963662;

		epx1	= ( cos(e_beta * pi / 180.) * (e1x_mm) + sin(e_beta * pi / 180.) * e1y_mm) * rad_callib;
		epy1	= ( -sin(e_beta * pi / 180.) * (e1x_mm) + cos(e_beta * pi / 180.) * e1y_mm) * rad_callib;
		}						// End IF: not in the node
	else
	{
		epx1 = -34693498;
		epy1 = -9276134903;
	}


	
	if (direction=='x') {outputvalue=epx1;}
	else {if (direction=='y') {outputvalue=epy1;} else {outputvalue=85734.;} }


	return outputvalue;
}



// **************************************************************
// **
// ** find channels in PIPICO: t1 [ns], m1,m2 [au], q1,q2 [e], s[mm], fieldE [V/cm] 
// **
// **************************************************************

double t2(double t1,double m1,double m2, double q1, double q2, double s, double fieldE, double  parameter[], double shift_ns = 0.0) {
	
	return 0.0; //CH_t2(t1, m1, m2, q1, q2, s, fieldE, shift_ns);
}


double t3(double t1, double t2, double m1,double m2, double m3, double q1, double q2, double q3, double s, double fieldE) 
{
	//convert to mks

	m1 = m1 * 1.661e-27;
	m2 = m2 * 1.661e-27;
	m3 = m3 * 1.661e-27;
	
	q1 = q1 * 1.60322e-19;
	q2 = q2 * 1.60322e-19;
	q3 = q3 * 1.60322e-19;
  
	t1 = t1 * 1e-9;
	t2 = t2 * 1e-9;

	s = s / 1000.0; 
	fieldE = fieldE * 100.0;

	double t3 = 1./(2.*fieldE*t1*t2*q3)*(-fieldE*t1*t1*t2*q1-fieldE*t1*t2*t2*q2+2.*m1*t2*s+2.*m2*t1*s+sqrt((-fieldE*t1*t1*t2*q1-fieldE*t1*t2*t2*q2+2.*m1*t2*s+2.*m2*t1*s)*(-fieldE*t1*t1*t2*q1-fieldE*t1*t2*t2*q2+2.*m1*t2*s+2.*m2*t1*s)+8*fieldE*m3*t1*t1*t2*t2*q3*s));
	//printf("%f\n",t3*1e+9);
	return t3 * 1e+9;
}


// **************************************************************
// **
// ** calculate momenta for electron x,y direction (magnetic field)
// ** using Mirko's functions.
// **
// **************************************************************

double electron_px(double etof, double xe, double ye, double  parameter[]) {

	double pex,pey,qe,w,me,a,b,pau,wigglepos;

	me = 9.1083E-31;
	qe = 1.602E-19;
	pau = me*300.e6/137.;

	double fieldB = parameter[1080] / 10000.;
	double fieldBns = parameter[1081];


	wigglepos =  fieldBns * 1e-9;
	if(1==0)fieldB = 2.*me*3.14152 / (qe * wigglepos);

	w = qe / me * fieldB;
	a = (1. - cos(w * etof*1.e-9)) / w;
	b = (sin(w * etof*1.e-9)) / w;

	pey = me * (-xe/1000. * a - b*ye/1000.) / (a*a + b*b);
	pex = me * (xe/1000. * b - a*ye/1000.) / (a*a + b*b);
	pex = pex / pau;
	pey = pey / pau;
	return pex;
}
double electron_py(double etof, double xe, double ye, double  parameter[]) {

	double pex,pey,qe,w,me,a,b,pau,wigglepos;

	me = 9.1083E-31;
	qe = 1.602E-19;
	pau = me*300.e6/137.;

	double fieldB = parameter[1080] / 10000.;
	double fieldBns = parameter[1081];

	wigglepos =  fieldBns * 1e-9;
	if(1==0)fieldB = 2.*me*3.14152 / (qe * wigglepos);

	w = qe / me * fieldB;
	a = (1. - cos(w * etof*1.e-9)) / w;
	b = (sin(w * etof*1.e-9)) / w;

	pey = me * (-xe/1000. * a - b*ye/1000.) / (a*a + b*b);
	pex = me * (xe/1000. * b - a*ye/1000.) / (a*a + b*b);
	pex = pex / pau;
	pey = pey / pau;
	return pey;
}




// ********************************************************************
// **
// ** calculate angle of Delta_Phi between Recoil and projectile
// **
// ********************************************************************

double deltaphi(double phi_projectile, double phi_molecule)
{
	double phi_p = phi_projectile + 180.;
	double phi_m = phi_molecule + 180.;
	double relphi;

	relphi = phi_p - phi_m;
		if (phi_p > phi_m) {
			if (relphi > 180.) {
				relphi = relphi - 360.;
			}
		}else {
			if (phi_p <= phi_m) {
				if (relphi < -180.) {
					relphi = relphi + 360.;
				}
			}
		}
			
	return relphi;
}


// ********************************************************************
// **
// ** Labframe transformation - 3 vectors are transformed into a new
// ** frame of reference
// **
// ********************************************************************

void labframe_transformation(__int32 direction, double a[], double b[], double c[]) {

//   Author: Achim Czasch, email: czasch@atom.uni-frankfurt.de
//
//   Function: the 3 vectors are transformed into a new frame of reference.
//             The new frame of reference is defined in the following way:
//             Vector  a  becomes the x-axis (|a|,0,0) of the new frame of reference.
//             The y-axis of the new system lies in the half plane which is spanned
//             by a and b. (It is not necessarily parallel to b)
//             The orientation of the new system is in accordance to the right-hand-rule:
//             Thumb is vector a. The new y-axis is in the half-plane which is spanned
//             by thumb and index-finger. The middle finger points in the direction of the
//             the z-axis.
//             The first argument (direction )can be +1 or -1. In case of -1 the inverse transformation
//             is performed.


	double norm_vector_x[3];
	double norm_vector_y[3];
	double norm_vector_z[3];

	double n;
	__int32 i;

	double transformation_matrix[3][3];
	
	n=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	norm_vector_x[0] = a[0]/n;
	norm_vector_x[1] = a[1]/n;
	norm_vector_x[2] = a[2]/n;

	norm_vector_z[0] = -b[1]*a[2]+b[2]*a[1];
	norm_vector_z[1] = -b[2]*a[0]+b[0]*a[2];
	norm_vector_z[2] = -b[0]*a[1]+b[1]*a[0];

	n=sqrt(norm_vector_z[0]*norm_vector_z[0] + norm_vector_z[1]*norm_vector_z[1] + norm_vector_z[2]*norm_vector_z[2]);
	norm_vector_z[0] = norm_vector_z[0]/n;
	norm_vector_z[1] = norm_vector_z[1]/n;
	norm_vector_z[2] = norm_vector_z[2]/n;

	norm_vector_y[0] = -a[1]*norm_vector_z[2]+a[2]*norm_vector_z[1];
	norm_vector_y[1] = -a[2]*norm_vector_z[0]+a[0]*norm_vector_z[2];
	norm_vector_y[2] = -a[0]*norm_vector_z[1]+a[1]*norm_vector_z[0];

	n=sqrt(norm_vector_y[0]*norm_vector_y[0] + norm_vector_y[1]*norm_vector_y[1] + norm_vector_y[2]*norm_vector_y[2]);
	norm_vector_y[0] = norm_vector_y[0]/n;
	norm_vector_y[1] = norm_vector_y[1]/n;
	norm_vector_y[2] = norm_vector_y[2]/n;

	for (i=0;i<3;i++) transformation_matrix[0][i] = norm_vector_x[i];
	for (i=0;i<3;i++) transformation_matrix[1][i] = norm_vector_y[i];
	for (i=0;i<3;i++) transformation_matrix[2][i] = norm_vector_z[i];


//         If direction =  +1: This transforms into the new system
//         If direction =  -1: Reverse transformation

	
	double vec_new[3];

	if (direction==-1) {
			vec_new[0]=a[0]*transformation_matrix[0][0]+a[1]*transformation_matrix[1][0]+a[2]*transformation_matrix[2][0];
			vec_new[1]=a[0]*transformation_matrix[0][1]+a[1]*transformation_matrix[1][1]+a[2]*transformation_matrix[2][1];
			vec_new[2]=a[0]*transformation_matrix[0][2]+a[1]*transformation_matrix[1][2]+a[2]*transformation_matrix[2][2];
			a[0] = vec_new[0];	a[1] = vec_new[1];	a[2] = vec_new[2];

			vec_new[0]=b[0]*transformation_matrix[0][0]+b[1]*transformation_matrix[1][0]+b[2]*transformation_matrix[2][0];
			vec_new[1]=b[0]*transformation_matrix[0][1]+b[1]*transformation_matrix[1][1]+b[2]*transformation_matrix[2][1];
			vec_new[2]=b[0]*transformation_matrix[0][2]+b[1]*transformation_matrix[1][2]+b[2]*transformation_matrix[2][2];
			b[0] = vec_new[0];	b[1] = vec_new[1];	b[2] = vec_new[2];

			vec_new[0]=c[0]*transformation_matrix[0][0]+c[1]*transformation_matrix[1][0]+c[2]*transformation_matrix[2][0];
			vec_new[1]=c[0]*transformation_matrix[0][1]+c[1]*transformation_matrix[1][1]+c[2]*transformation_matrix[2][1];
			vec_new[2]=c[0]*transformation_matrix[0][2]+c[1]*transformation_matrix[1][2]+c[2]*transformation_matrix[2][2];
			c[0] = vec_new[0];	c[1] = vec_new[1];	c[2] = vec_new[2];
	} else {
		if (direction==+1) {
				vec_new[0]=a[0]*transformation_matrix[0][0]+a[1]*transformation_matrix[0][1]+a[2]*transformation_matrix[0][2];
				vec_new[1]=a[0]*transformation_matrix[1][0]+a[1]*transformation_matrix[1][1]+a[2]*transformation_matrix[1][2];
				vec_new[2]=a[0]*transformation_matrix[2][0]+a[1]*transformation_matrix[2][1]+a[2]*transformation_matrix[2][2];
				a[0] = vec_new[0];	a[1] = vec_new[1];	a[2] = vec_new[2];

				vec_new[0]=b[0]*transformation_matrix[0][0]+b[1]*transformation_matrix[0][1]+b[2]*transformation_matrix[0][2];
				vec_new[1]=b[0]*transformation_matrix[1][0]+b[1]*transformation_matrix[1][1]+b[2]*transformation_matrix[1][2];
				vec_new[2]=b[0]*transformation_matrix[2][0]+b[1]*transformation_matrix[2][1]+b[2]*transformation_matrix[2][2];
				b[0] = vec_new[0];	b[1] = vec_new[1];	b[2] = vec_new[2];

				vec_new[0]=c[0]*transformation_matrix[0][0]+c[1]*transformation_matrix[0][1]+c[2]*transformation_matrix[0][2];
				vec_new[1]=c[0]*transformation_matrix[1][0]+c[1]*transformation_matrix[1][1]+c[2]*transformation_matrix[1][2];
				vec_new[2]=c[0]*transformation_matrix[2][0]+c[1]*transformation_matrix[2][1]+c[2]*transformation_matrix[2][2];
				c[0] = vec_new[0];	c[1] = vec_new[1];	c[2] = vec_new[2];
		} else {
				a[0] = 0.; a[1] = 0.; a[2] = 0.;
				b[0] = 0.; b[1] = 0.; b[2] = 0.;
				c[0] = 0.; c[1] = 0.; c[2] = 0.;
		}
	}
}
   �v  l   I D S _ S O U R C E _ C P P   S O R T _ A N D _ W R I T E _ N T U P L E _ C P P         0         #include "OS_Version.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TFile.h"

#include <math.h>
#include "rootstuff.h"
#include "Histo.h"
#include "TF1.h"
//#include "TMinuit.h"
#include "TNtupleD.h"

#include "functions.h"
#include "Ueberstruct.h"


//#include "ADC_analysis.h"
#include "../ADC/ADC_meta.h"

#include "console.h"

__int32 detector_stuff(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[]);

char ntuple_identifier[2048];
bool create_ntuple_identifier;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
__int32 sort_and_write_NTuple(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[], double timestamp)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
#ifdef _DEBUG
	assert(Ueber);
#endif
	
	double NTupleData[250]; 
	double scan_val[16];

	bool writeNTuple = false;

	if (Ueber->error) {
		Ueber->stop_reading_files = true;
		return 0;
	}

	Ueber->write_this_event_into_sorted_LMF = false;
	Ueber->start_new_root_file = false;

	int NRec  = 0;
	int NElec = 0;
	int NPro  = 0;

	int NRmax = (int)(parameter[61]+0.01);
	int NEmax = (int)(parameter[62]+0.01);
	int NPmax = (int)(parameter[60]+0.01);

	if(eventcounter == 0) {
		Ueber->eventswritten = 0;
		create_ntuple_identifier = true;
	}

	Histo * Hist = Ueber->Hist;

	if (!detector_stuff(eventcounter,Ueber,parameter)) return 0;

	double * tdc_ns = Ueber->tdc_ns; // this line must be below   if (!detectorstuff()...
	__int32 * cnt = Ueber->cnt; // this line must be below   if (!detectorstuff()...

	double bunchspacing = 0.;
	int bunchmarker_channel = -1;
	bool useBM = false;

	double pdet_size	= 5. + parameter[122];				// size of projectile detector (for displaying spectra properly)
	double rdet_size	= 5. + parameter[222];				// size of recoil detector (for displaying spectra properly)
	double edet_size	= 5. + parameter[322];				// size of electron detector (for displaying spectra properly)

	if(Ueber->scan_channel>-1) {
		Ueber->scan_in_data = true;
		for(int i=0;i<Ueber->scan_num_vals;i++)
			scan_val[i] = (double)tdc_ns[Ueber->scan_channel*NUM_IONS+i]/(double)Ueber->scan_factor;
	}

	// Include your analysis of the TDC data here.
//------------------------------------------------------------------------------------------------
	// Access the results from detector reconstruction using the detector struct.
	// The name of the struct is "proj", "rec" or "elec" 
	// depending on the detector. The struct is created as a part of "Ueber". 
	// Thus it is addressed using "Ueber->".
	//
	// For example the electron detector is accessed by "Ueber->elec".
	//
	// See Ueberstruct.h for the definition of the detector struct. However, the following 
	// are the important variables:
	//
	//		__int32		number_of_reconstructed_hits	- number of hits or this event
	//		__int32		method[NUM_HIT]					- reconstruction method that was used
	//		double	x[NUM_HIT]						- x position of impact
	//		double	y[NUM_HIT]						- y position of impact
	//		double	time[NUM_HIT]					- MCP time 
	//
	// The first hit is NUM_HIT=0.
	//
	// Another example: in order to write the x-position of the second hit on a 
	// recoil detector to xpos use:
	//		double xpos = Ueber->rec->x[1];
	//
//------------------------------------------------------------------------------------------------



// Feed things into ColAHelL if enabled...	
// ("evt->"-struct is only available when using ColAHelL!)
// -----------------------------------------------------------
	CH::CH_event_struct* evt = 0;

	if(Ueber->use_CH) {
		useBM = Ueber->CH->GetUseBM();		
		if(useBM) {
			bunchspacing = Ueber->CH->GetBunchSpacing();
			bunchmarker_channel = Ueber->CH->GetBMChannel();
		}
		double ds = -1.0;
		ds =  Ueber->CH->GetPDetSize();
		if(ds>0.0)
			pdet_size = ds + 5.0;
		ds =  Ueber->CH->GetIDetSize();
		if(ds>0.0)
			rdet_size = ds + 5.0;
		ds =  Ueber->CH->GetEDetSize();
		if(ds>0.0)
			edet_size = ds + 5.0;

		// Limit data to max. 64 hits 
		Ueber->elec->number_of_reconstructed_hits = (Ueber->elec->number_of_reconstructed_hits<65) ? Ueber->elec->number_of_reconstructed_hits : 64;
		Ueber->rec->number_of_reconstructed_hits = (Ueber->rec->number_of_reconstructed_hits<65) ? Ueber->rec->number_of_reconstructed_hits : 64;
		Ueber->proj->number_of_reconstructed_hits = (Ueber->proj->number_of_reconstructed_hits<65) ? Ueber->proj->number_of_reconstructed_hits : 64;
		
		// copy event data to ColAHelL
		Ueber->CH->ResetEvent();
		for(int i=0;i<Ueber->elec->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(EL, Ueber->elec->x[i], Ueber->elec->y[i], Ueber->elec->time[i], (double)Ueber->elec->method[i]); 	
		for(int i=0;i<Ueber->rec->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(IO, Ueber->rec->x[i], Ueber->rec->y[i], Ueber->rec->time[i], (double)Ueber->rec->method[i]); 
		for(int i=0;i<Ueber->proj->number_of_reconstructed_hits;i++) 
			Ueber->CH->SetEvent_xyt(PR, Ueber->proj->x[i], Ueber->proj->y[i], Ueber->proj->time[i], (double)Ueber->proj->method[i]); 
		// get the bunchmarker time from the tdc array...
		Ueber->CH->SetBMTime(tdc_ns, cnt);
		// copy scan data to ColAHelL
		for(int i=0;i<Ueber->scan_num_vals;i++)
			Ueber->CH->SetEvent_scanvals(scan_val[i]);

		// calculate TOFs using ColAHelL definitions
		Ueber->CH->CalculateTOFs();
		evt = Ueber->CH->GetEvent();

// IPA needs food, as well... Feeding is done here...
// (IPA works only together with ColAHelL!)
// -----------------------------------------------------------
		if(Ueber->ipa_enable && Ueber->IPA) { 
			if(Ueber->IPA->GetChannel() == -1) {
				Ueber->IPA->append_data(evt);
				return 0;
			} else {
				std::vector<CH::CH_event_struct> IPAPrsEvts = Ueber->CH->PresortEvent();
				for(int i=0;i<IPAPrsEvts.size();i++) {
					Ueber->IPA->append_data(&IPAPrsEvts[i]);
				}
			}
			return 0;
		}
	}

//---------------------------------------------------------------------------------------------------------------------------
// some typical spectra ....
//---------------------------------------------------------------------------------------------------------------------------

	if(useBM)
	{
		Hist->fill1(100,"bunchmarker1",tdc_ns[bunchmarker_channel*NUM_IONS+0],1.,"Time Bunchmarker #1",10000,-10000.,10000.,"bunchmarker [ns]","raw");
		Hist->fill1(101,"bunchmarker_hit",cnt[bunchmarker_channel],1.,"Bunchmarker hits",63,-1.25,30.25,"bunchmarker hits","raw");
	}

	if(Ueber->scan_channel>-1) {
		Hist->fill1(108,"scan_hits",cnt[Ueber->scan_channel],1.,"number of hits scan channel",1000,-1.,25.,"cnt","raw");
		Hist->fill1(109,"scan value",scan_val[0],1.,"scan value",1000,0.,100.,"cnt","raw");
	}

	if(Ueber->use_CH) {
		Hist->fill1(105,"phit",evt->p.num_hits,1.,"projectile number of reconstructed hits",1000,-1.,25.,"evt->p.num_hits","raw");
		Hist->fill1(106,"ehit",evt->e.num_hits,1.,"electron number of reconstructed hits",1000,-1.,25.,"evt->e.num_hits","raw");
		Hist->fill1(107,"rhit",evt->r.num_hits,1.,"recoil number of reconstructed hits",1000,-1.,25.,"evt->r.num_hits","raw");


		if (evt->p.num_hits>0 )
		{	// begin if: one projectile
			Hist->fill2(110,"proj1xy",evt->p.x[0],evt->p.y[0],1.,"Projectile position #1",400,-1.*pdet_size,pdet_size,"proj1 x [mm]",400,-1.*pdet_size,pdet_size,"proj1 y [mm]","raw");
		}	// end if: one projectile


		if (evt->r.num_hits>0 )
		{	//  begin if: one recoil
			Hist->fill2(140,"rec1xy",evt->r.x[0],evt->r.y[0],1.,"Recoil position #1",400,-1.*rdet_size,rdet_size,"rec1 x [mm]",400,-1.*rdet_size,rdet_size,"rec1 y [mm]","raw");

			Hist->fill1(145,"rec1tof",evt->r.tof[0],1.,"Recoil TOF #1",10000,-10000.,60000.,"rec1 TOF [ns]","raw");
			Hist->fill1(1451,"rec1tof_short",evt->r.tof[0],1.,"Recoil TOF #1",10000,-2500.,7500.,"rec1 TOF [ns]","raw");
			Hist->fill2(146,"rec1wiggle",evt->r.tof[0],(sqrt(evt->r.x[0]*evt->r.x[0]+evt->r.y[0]*evt->r.y[0])),1.,"Recoil wiggle #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.,rdet_size,"rec1 radius [mm]","raw");
			Hist->fill2(147,"rec1xfish",evt->r.tof[0],evt->r.x[0],1.,"Recoil x-fish #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 x-pos [mm]","raw");
			Hist->fill2(1471,"rec1xfish_short",evt->r.tof[0],evt->r.x[0],1.,"Recoil x-fish #1",500,-2500.,7500.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 x-pos [mm]","raw");
			Hist->fill2(148,"rec1yfish",evt->r.tof[0],evt->r.y[0],1.,"Recoil y-fish #1",500,0.,60000.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 y-pos [mm]","raw");
			Hist->fill2(1481,"rec1yfish_short",evt->r.tof[0],evt->r.y[0],1.,"Recoil y-fish #1",500,-2500.,7500.,"rec1 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec1 y-pos [mm]","raw");
		}	//  end if: one recoil


		if (evt->r.num_hits>1 )
		{	//  begin if: two recoils
			Hist->fill2(150,"rec2xy",evt->r.x[1],evt->r.y[1],1.,"Recoil position #2",400,-1.*rdet_size,rdet_size,"rec2 x [mm]",400,-1.*rdet_size,rdet_size,"rec2 y [mm]","raw");

			Hist->fill1(155,"rec2tof",evt->r.tof[1],1.,"Recoil TOF #2",1000,-10000.,60000.,"rec2 TOF [ns]","raw");
			Hist->fill1(1551,"rec2tof_short",evt->r.tof[1],1.,"Recoil TOF #2",1000,-2500.,7500.,"rec2 TOF [ns]","raw");
			Hist->fill2(156,"rec2wiggle",evt->r.tof[1],(sqrt(evt->r.x[1]*evt->r.x[1]+evt->r.y[1]*evt->r.y[1])),1.,"Recoil wiggle #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.,rdet_size,"rec2 radius [mm]","raw");
			Hist->fill2(157,"rec2xfish",evt->r.tof[1],evt->r.x[1],1.,"Recoil x-fish #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 x-pos [mm]","raw");
			Hist->fill2(1571,"rec2xfish_short",evt->r.tof[1],evt->r.x[1],1.,"Recoil x-fish #2",500,-2500.,7500.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 x-pos [mm]","raw");
			Hist->fill2(158,"rec2yfish",evt->r.tof[1],evt->r.y[1],1.,"Recoil y-fish #2",500,0.,60000.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 y-pos [mm]","raw");
			Hist->fill2(1581,"rec2yfish_short",evt->r.tof[1],evt->r.y[1],1.,"Recoil y-fish #2",500,-2500.,7500.,"rec2 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec2 y-pos [mm]","raw");
		
			Hist->fill2(160,"recsumdiff",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff",1000,10000.,40000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1601,"recsumdiff_short",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff",1000,-5000.,15000.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(161,"recsumdiff2",2.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (2*r1 + r2)",800,0.,50000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1611,"recsumdiff2_short",2.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (2*r1 + r2)",800,-7500.,22500.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(162,"recsumdiff3",3.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (3*r1 + r2)",800,0.,60000.,"rec TOF sum [ns]",610,-6000.,100.,"rec TOF diff [ns]","raw");
			Hist->fill2(1621,"recsumdiff3_short",3.*evt->r.tof[0]+evt->r.tof[1],evt->r.tof[0]-evt->r.tof[1],1.,"Recoil TOF sum vs diff (3*r1 + r2)",800,-10000.,30000.,"rec TOF sum [ns]",410,-4000.,100.,"rec TOF diff [ns]","raw");
		
			Hist->fill2(163,"pipico",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum",500,0.,50000.,"rec1 TOF [ns]",500,0.,50000.,"rec2 TOF [ns]","raw");
			Hist->fill2(1631,"pipico_short",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum",800,-2500.,7500.,"rec1 TOF [ns]",800,-2500.,7500.,"rec2 TOF [ns]","raw");
		
			double phi0 = atan2(evt->r.x[0],evt->r.y[0])*180./3.14152;
			double phi1 = atan2(evt->r.x[1],evt->r.y[1])*180./3.14152;

			if(fabs(deltaphi(phi0,phi1)-180.0)<10.0)
			{
				Hist->fill2(165,"pipico_momcon",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum (coarse) /w simple mom. check",500,0.,50000.,"rec1 TOF [ns]",500,0.,50000.,"rec2 TOF [ns]","raw");
				Hist->fill2(1651,"pipico_momcon_short",evt->r.tof[0],evt->r.tof[1],1.,"PIPICO spectrum (fine) /w simple mom. check",800,-2500.,7500.,"rec1 TOF [ns]",800,-2500.,7500.,"rec2 TOF [ns]","raw");
			}
		}	//  end if: two recoils
	

		if ( evt->r.num_hits>2)
		{	// begin if: three recoils
			Hist->fill2(170,"rec3xy",evt->r.x[2],evt->r.y[2],1.,"Recoil position #3",400,-1.*rdet_size,rdet_size,"rec3 x [mm]",400,-1.*rdet_size,rdet_size,"rec3 y [mm]","raw");

			Hist->fill1(175,"rec3tof",evt->r.tof[2],1.,"Recoil TOF #3",1000,-10000.,60000.,"rec3 TOF [ns]","raw");
			Hist->fill1(1751,"rec3tof_short",evt->r.tof[2],1.,"Recoil TOF #3",1000,-2500.,7500.,"rec3 TOF [ns]","raw");
			Hist->fill2(176,"rec3wiggle",evt->r.tof[2],(sqrt(evt->r.x[2]*evt->r.x[2]+evt->r.y[2]*evt->r.y[2])),1.,"Recoil wiggle #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.,rdet_size,"rec3 radius [mm]","raw");
			Hist->fill2(177,"rec3xfish",evt->r.tof[2],evt->r.x[2],1.,"Recoil x-fish #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 x-pos [mm]","raw");
			Hist->fill2(1771,"rec3xfish_short",evt->r.tof[2],evt->r.x[2],1.,"Recoil x-fish #3",500,-2500.,7500.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 x-pos [mm]","raw");
			Hist->fill2(178,"rec3yfish",evt->r.tof[2],evt->r.y[2],1.,"Recoil y-fish #3",500,0.,60000.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 y-pos [mm]","raw");
			Hist->fill2(1781,"rec3yfish_short",evt->r.tof[2],evt->r.y[2],1.,"Recoil y-fish #3",500,-2500.,7500.,"rec3 TOF [ns]",200,-1.*rdet_size,rdet_size,"rec3 y-pos [mm]","raw");
			Hist->fill2(179,"pipipico",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[2],1.,"PIPIPICO spectrum (coarse)",800,0.,80000.,"rec1+rec2 TOF [ns]",500,0.,50000.,"rec3 TOF [ns]","raw");
			Hist->fill2(1791,"pipipico_short",evt->r.tof[0]+evt->r.tof[1],evt->r.tof[2],1.,"PIPIPICO spectrum (fine)",500,-5000.,15000.,"rec1+rec2 TOF [ns]",500,-2500.,7500.,"rec3 TOF [ns]","raw");
		}	// end if: three recoils


		if ( evt->e.num_hits>0 )
		{	// begin if: one electron
			Hist->fill2(180,"elec1xy",evt->e.x[0],evt->e.y[0],1.,"Electron position #1",400,-1.*edet_size,edet_size,"elec x [mm]",400,-1.*edet_size,edet_size,"elec y [mm]","raw");

			if(useBM && (bunchspacing > 0.0))
			{
				Hist->fill1(185,"elec1tof",evt->e.tof[0],1.,"Electron TOF #1",int(6*bunchspacing),-1.5*bunchspacing,1.5*bunchspacing,"elec1 TOF [ns]","raw");
				Hist->fill2(186,"elec1wiggle",evt->e.tof[0],(sqrt(evt->e.x[0]*evt->e.x[0]+evt->e.y[0]*evt->e.y[0])),1.,"Electron wiggle #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.,edet_size,"elec1 radius [mm]","raw");
				Hist->fill2(187,"elec1xfish",evt->e.tof[0],evt->e.x[0],1.,"Electron x-fish #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 x-pos [mm]","raw");
				Hist->fill2(188,"elec1yfish",evt->e.tof[0],evt->e.y[0],1.,"Electron y-fish #1",500,-5.,1.1*bunchspacing,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 y-pos [mm]","raw");
			}
			else 
			{
				Hist->fill1(185,"elec1tof",evt->e.tof[0],1.,"Electron TOF #1",5000,-2500.,2500.,"elec1 TOF [ns]","raw");
				Hist->fill2(186,"elec1wiggle",evt->e.tof[0],(sqrt(evt->e.x[0]*evt->e.x[0]+evt->e.y[0]*evt->e.y[0])),1.,"Electron wiggle #1",500,0.,250.,"elec1 TOF [ns]",200,-1.,edet_size,"elec1 radius [mm]","raw");
				Hist->fill2(187,"elec1xfish",evt->e.tof[0],evt->e.x[0],1.,"Electron x-fish #1",500,0.,250.,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 x-pos [mm]","raw");
				Hist->fill2(188,"elec1yfish",evt->e.tof[0],evt->e.y[0],1.,"Electron y-fish #1",500,0.,250.,"elec1 TOF [ns]",200,-1.*edet_size,edet_size,"elec1 y-pos [mm]","raw");
			}

		}	// end if: one electron

		if ( evt->e.num_hits>1 )
		{	// begin if: two electrons
			Hist->fill2(190,"elec2xy",evt->e.x[1],evt->e.y[1],1.,"Electron position #2",400,-1.*edet_size,edet_size,"elec2 x [mm]",400,-1.*edet_size,edet_size,"elec2 y [mm]","raw");
			double dr = sqrt((evt->e.x[0] - evt->e.x[1])*(evt->e.x[0] - evt->e.x[1]) + (evt->e.y[0] - evt->e.y[1])*(evt->e.y[0] - evt->e.y[1]));
			//Hist->fill2(199, "dtof_vs_dr", fabs(evt->e.tof[0] - evt->e.tof[1]),dr,1.,"dToF vs dR",100,-50,50, )


			if(useBM && (bunchspacing > 0.0))
			{
				Hist->fill1(195,"elec2tof",evt->e.tof[1],1.,"Electron TOF #2",int(6*bunchspacing),-1.5*bunchspacing,1.5*bunchspacing,"elec2 TOF [ns]","raw");
				Hist->fill2(196,"elec2wiggle",evt->e.tof[1],(sqrt(evt->e.x[1]*evt->e.x[1]+evt->e.y[1]*evt->e.y[1])),1.,"Electron wiggle #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.,edet_size,"elec2 radius [mm]","raw");
				Hist->fill2(197,"elec2xfish",evt->e.tof[1],evt->e.x[1],1.,"Electron x-fish #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 x-pos [mm]","raw");
				Hist->fill2(198,"elec2yfish",evt->e.tof[1],evt->e.y[1],1.,"Electron y-fish #2",500,-5.,1.1*bunchspacing,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 y-pos [mm]","raw");
			}
			else
			{
				Hist->fill1(195,"elec2tof",evt->e.tof[1],1.,"Electron TOF #2",2000,-1000.,1000.,"elec2 TOF [ns]","raw");
				Hist->fill2(196,"elec2wiggle",evt->e.tof[1],(sqrt(evt->e.x[1]*evt->e.x[1]+evt->e.y[1]*evt->e.y[1])),1.,"Electron wiggle #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.,edet_size,"elec2 radius [mm]","raw");
				Hist->fill2(197,"elec2xfish",evt->e.tof[1],evt->e.x[1],1.,"Electron x-fish #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 x-pos [mm]","raw");
				Hist->fill2(198,"elec2yfish",evt->e.tof[1],evt->e.y[1],1.,"Electron y-fish #2",500,-1000.,1000.,"elec2 TOF [ns]",200,-1.*edet_size,edet_size,"elec2 y-pos [mm]","raw");
			}
		}	// end if: two electrons

		if (evt->e.num_hits == 2) {
			double dpos_mm = sqrt((evt->e.x[0] - evt->e.x[1])*(evt->e.x[0] - evt->e.x[1]) + (evt->e.y[0] - evt->e.y[1])*(evt->e.y[0] - evt->e.y[1]));
			double dtof_ns = fabs(evt->e.time[0] - evt->e.time[1]);

			Hist->fill1(80, "dpos", dpos_mm, 1., "dpos_mm", 2000, -50., 50., "dpos [mm]", "DPA");
			Hist->fill1(81, "dtof", dtof_ns, 1., "dtof_mm", 2000, 0., 1000., "dtof [ns]", "DPA");
			Hist->fill2(82, "dtof_vs_dpos", dtof_ns, fabs(dpos_mm), 1., "dtof_vs_dpos", 250, 0, 1000., "dtof [ns]", 200, 0, 50, "dpos [mm]", "DPA");
		}
	} // end if: "use CH"

//---------------------------------------------------------------------------------------------------------------------------
// add missing spectra etc. here....
//---------------------------------------------------------------------------------------------------------------------------

	ADC_meta * adc = (ADC_meta*)(Ueber->adc_meta);

	if (adc->m_event->channels[0]) {
		if (adc->m_event->channels[0]->NbrPulses) {
			if (adc->m_event->channels[0]->Pulses[0].NbrPeaks) {
				//...
			}
		}
	}



//---------------------------------------------------------------------------------------------------------------------------
//	                            P R E S O R T E R                                      
//---------------------------------------------------------------------------------------------------------------------------

	// Dublicate event in case it fits to several presorters... (default: true)
	bool dublicate = true;

	// presort event
	std::vector<CH::CH_event_struct> PrsEvts = Ueber->CH->PresortEvent();
		
	// write resulting data to file 
	for(unsigned int prsnum=0;prsnum<PrsEvts.size();prsnum++) {
		
		writeNTuple = true;

		NRec = PrsEvts[prsnum].r.num_hits;
		if(NRec > NRmax)
			NRec = NRmax;
		
		NElec = PrsEvts[prsnum].e.num_hits;
		if(NElec > NEmax)
			NElec = NEmax;

		NPro = PrsEvts[prsnum].p.num_hits;
		if(NPro > NPmax)
			NPro = NPmax;

		NTupleData[0] = double(PrsEvts[prsnum].reaction);
		NTupleData[1] = double(NElec);
		NTupleData[2] = double(NRec);
		NTupleData[3] = double(NPro);
		NTupleData[4] = PrsEvts[prsnum].bunchmarker;
        NTupleData[5] = timestamp;
        NTupleData[6] = (double)Ueber->LMF_input->Starttime;
        NTupleData[7] = (double)Ueber->LMF_input->Stoptime; 

		int absind = 8;

		if(Ueber->scan_channel>-1) {
			for(int i=0;i<Ueber->scan_num_vals;i++)
			{
				NTupleData[absind+i]=scan_val[i];
			}
			absind += Ueber->scan_num_vals;
		}

		for(int i=0;i<NRec;i++)
		{	
			NTupleData[absind+i*5]   = double(PrsEvts[prsnum].r.x[i]);
			NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].r.y[i]);
			NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].r.time[i]);
			NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].r.tof[i]);
			NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].r.method[i]);
		}	
		absind += NRmax*5;

		for(int i=0;i<NElec;i++)
		{	
			NTupleData[absind+i*5]   = double(PrsEvts[prsnum].e.x[i]);
			NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].e.y[i]);
			NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].e.time[i]);
			NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].e.tof[i]);
			NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].e.method[i]);
		}	
		absind += NEmax*5;

		for(int i=0;i<NPro;i++)
		{	
			NTupleData[absind+i*5]   = double(PrsEvts[prsnum].p.x[i]);
			NTupleData[absind+i*5+1] = double(PrsEvts[prsnum].p.y[i]);
			NTupleData[absind+i*5+2] = double(PrsEvts[prsnum].p.time[i]);
			NTupleData[absind+i*5+3] = double(PrsEvts[prsnum].p.tof[i]);
			NTupleData[absind+i*5+4] = double(PrsEvts[prsnum].p.method[i]);
		}	
		absind += NPmax*5;

	//////////////////////////////////////////////////////////////////////////
	//	write data to ntuple
	//////////////////////////////////////////////////////////////////////////

		char tmp[200];

		if(create_ntuple_identifier) {
			create_ntuple_identifier = false;
	
		sprintf(ntuple_identifier,"reaction:ehit:rhit:phit:bunchmarker:timestamp:LMFStart:LMFStop");

		if(Ueber->scan_channel>-1) {
			for(int i=0;i<Ueber->scan_num_vals;i++)
			{
				sprintf(tmp,":scanval%d",i);
				strcat(ntuple_identifier,tmp);
			}
		}

		for(int i=1;i<=NRmax;i++)
		{
			sprintf(tmp,":r%dx:r%dy:r%dmcp:r%dtof:r%dflag",i,i,i,i,i);
			strcat(ntuple_identifier,tmp);
		}
		for(int i=1;i<=NEmax;i++)
		{
			sprintf(tmp,":e%dx:e%dy:e%dmcp:e%dtof:e%dflag",i,i,i,i,i);
			strcat(ntuple_identifier,tmp);
		}
		for(int i=1;i<=NPmax;i++)
		{
			sprintf(tmp,":p%dx:p%dy:p%dmcp:p%dtof:p%dflag",i,i,i,i,i);
			strcat(ntuple_identifier,tmp);
			}
		}

		Hist->NTupleD(0,"Data","BESSY2012",ntuple_identifier, 32000, NTupleData);
		++Ueber->eventswritten;

		// plot some histogramms of all presorted events
		Hist->fill1(200,"reaction_ntuple",PrsEvts[prsnum].reaction,1.,"Reaction channel flag",203,-1.25,100.25,"reaction_flag","ntuple");
	
		Hist->fill1(205,"phit",PrsEvts[prsnum].p.num_hits,1.,"projectile number of reconstructed hits",1000,-1.,25.,"evt->p.num_hits","ntuple");
		Hist->fill1(206,"ehit",PrsEvts[prsnum].e.num_hits,1.,"electron number of reconstructed hits",1000,-1.,25.,"evt->e.num_hits","ntuple");
		Hist->fill1(207,"rhit",PrsEvts[prsnum].r.num_hits,1.,"recoil number of reconstructed hits",1000,-1.,25.,"evt->r.num_hits","ntuple");
			

		if (NElec>0)
		{	
			Hist->fill2(210,"elec1xy_ntuple",PrsEvts[prsnum].e.x[0],PrsEvts[prsnum].e.y[0],1.,"Electron position #1 (ntuple)",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[0] [mm]",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
			Hist->fill1(215,"elec1tof_ntuple",PrsEvts[prsnum].e.tof[0],1.,"Electron TOF #1 (ntuple)",1000,-500.0,500.0,"PrsEvts[prsnum].e.tof[0] [ns]","ntuple");
			Hist->fill2(216,"elec1wiggle_ntuple",PrsEvts[prsnum].e.tof[0],(sqrt(PrsEvts[prsnum].e.x[0]*PrsEvts[prsnum].e.x[0]+PrsEvts[prsnum].e.y[0]*PrsEvts[prsnum].e.y[0])),1.,"Electron wiggle #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.,edet_size,"elec_radius [mm]","ntuple");
			Hist->fill2(217,"elec1xfish_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.x[0],1.,"Electron x-fish #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[0] [mm]","ntuple");
			Hist->fill2(218,"elec1yfish_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.y[0],1.,"Electron y-fish #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
			if ( fabs(PrsEvts[prsnum].e.x[0])<10. ) Hist->fill2(219,"elec1yfishcut_ntuple",PrsEvts[prsnum].e.tof[0],PrsEvts[prsnum].e.y[0],1.,"Electron y-fish (x<10) #1 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[0] [ns]",20,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[0] [mm]","ntuple");
		}	

		if (NElec>1)
		{	
			Hist->fill2(220,"elec2xy_ntuple",PrsEvts[prsnum].e.x[1],PrsEvts[prsnum].e.y[1],1.,"Electron position #2 (ntuple)",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[1] [mm]",400,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[1] [mm]","ntuple");
			Hist->fill1(225,"elec2tof_ntuple",PrsEvts[prsnum].e.tof[1],1.,"Electron TOF #2 (ntuple)",1000,-500.,500.,"PrsEvts[prsnum].e.tof[1] [ns]","ntuple");
			Hist->fill2(226,"elec2wiggle_ntuple",PrsEvts[prsnum].e.tof[1],(sqrt(PrsEvts[prsnum].e.x[1]*PrsEvts[prsnum].e.x[1]+PrsEvts[prsnum].e.y[1]*PrsEvts[prsnum].e.y[1])),1.,"Electron wiggle #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.,edet_size,"elec2_radius [mm]","ntuple");
			Hist->fill2(227,"elec2xfish_ntuple",PrsEvts[prsnum].e.tof[1],PrsEvts[prsnum].e.x[1],1.,"Electron x-fish #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.x[1] [mm]","ntuple");
			Hist->fill2(228,"elec2yfish_ntuple",PrsEvts[prsnum].e.tof[1],PrsEvts[prsnum].e.y[1],1.,"Electron y-fish #2 (ntuple)",250,0.,250.,"PrsEvts[prsnum].e.tof[1] [ns]",200,-1.*edet_size,edet_size,"PrsEvts[prsnum].e.y[1] [mm]","ntuple");
		}	

		if (NRec>0)
		{	
			Hist->fill2(250,"rec1xy_ntuple",PrsEvts[prsnum].r.x[0],PrsEvts[prsnum].r.y[0],1.,"Recoil position #1 (ntuple)",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[0] [mm]",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[0] [mm]","ntuple");
			Hist->fill1(255,"rec1tof_ntuple",PrsEvts[prsnum].r.tof[0],1.,"Recoil TOF #1 (ntuple)",30100,-100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]","ntuple");
			Hist->fill2(256,"rec1wiggle_ntuple",PrsEvts[prsnum].r.tof[0],(sqrt(PrsEvts[prsnum].r.x[0]*PrsEvts[prsnum].r.x[0]+PrsEvts[prsnum].r.y[0]*PrsEvts[prsnum].r.y[0])),1.,"Recoil wiggle #1 (ntuple)",500,10000.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.,rdet_size,"rec1_radius [mm]","ntuple");
			Hist->fill2(257,"rec1xfish_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.x[0],1.,"Recoil x-fish #1 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[0] [mm]","ntuple");
			Hist->fill2(258,"rec1yfish_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.y[0],1.,"Recoil y-fish #1 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[0] [mm]","ntuple");
			Hist->fill2(259,"rec1tof_flag_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].reaction,1.,"Recoil TOF Reaction channel #1 (ntuple)",100,5000.,30000.,"PrsEvts[prsnum].r.tof[0] [ns]",405,-1.125,100.125,"reaction_flag","ntuple");
		}	

		if (NRec>1)
		{	
			Hist->fill2(260,"rec2xy_ntuple",PrsEvts[prsnum].r.x[1],PrsEvts[prsnum].r.y[1],1.,"Recoil position #2 (ntuple)",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[1] [mm]",400,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[1] [mm]","ntuple");
			Hist->fill1(265,"rec2tof_ntuple",PrsEvts[prsnum].r.tof[1],1.,"Recoil TOF #2 (ntuple)",30100,-100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
			Hist->fill2(266,"rec2wiggle_ntuple",PrsEvts[prsnum].r.tof[1],(sqrt(PrsEvts[prsnum].r.x[1]*PrsEvts[prsnum].r.x[1]+PrsEvts[prsnum].r.y[1]*PrsEvts[prsnum].r.y[1])),1.,"Recoil wiggle #2 (ntuple)",500,10000.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.,rdet_size,"rec2_radius [mm]","ntuple");
			Hist->fill2(267,"rec2xfish_ntuple",PrsEvts[prsnum].r.tof[1],PrsEvts[prsnum].r.x[1],1.,"Recoil x-fish #2 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.x[1] [mm]","ntuple");
			Hist->fill2(268,"rec2yfish_ntuple",PrsEvts[prsnum].r.tof[1],PrsEvts[prsnum].r.y[1],1.,"Recoil y-fish #2 (ntuple)",5000,100.,30000.,"PrsEvts[prsnum].r.tof[1] [ns]",200,-1.*rdet_size,rdet_size,"PrsEvts[prsnum].r.y[1] [mm]","ntuple");
			Hist->fill2(270,"pipico_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.tof[1],1.,"PIPICO spectrum (ntuple)",500,0.,50000.,"PrsEvts[prsnum].r.tof[0] [ns]",500,0.,50000.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
			Hist->fill2(271,"pipico_short_ntuple",PrsEvts[prsnum].r.tof[0],PrsEvts[prsnum].r.tof[1],1.,"PIPICO spectrum fine (ntuple)",800,-2500.,7500.,"PrsEvts[prsnum].r.tof[0] [ns]",800,-2500.,7500.,"PrsEvts[prsnum].r.tof[1] [ns]","ntuple");
		}		

		if (NPro>0)
		{	
			Hist->fill2(290,"proj1xy_ntuple",PrsEvts[prsnum].p.x[0],PrsEvts[prsnum].p.y[0],1.,"Projectile position #1 (ntuple)",400,-1.*pdet_size,pdet_size,"PrsEvts[prsnum].p.x[0] [mm]",400,-1.*pdet_size,pdet_size,"PrsEvts[prsnum].p.y[0] [mm]","ntuple");
		}	

		// ColAHelL momentum stuff kicks in here! 
		// (comment in if you don't mind processing time and need momenta etc. already at this point...)
		// -----------------------------------------------------------
		/*
		if(Ueber->use_CH)  {
			int EvtBelongsToReactions = Ueber->CH->EvtBelongsToReactions((int)PrsEvts[prsnum].reaction);
			// We need to loop, as there might be several reactions using the same channel number...
			// (e.g. when people use the ANY presorter...)
			for(int i=0;i<EvtBelongsToReactions;i++) {
				// process particles...
				Ueber->CH->ProcessEvent(&PrsEvts[prsnum],i);
				// fill std. ColAHelL histograms
				Ueber->CH->FillHistograms();
			}
		}
		*/
		if(dublicate == false)
			break;	
	}
				
	if(parameter[57]>0.5) {
		unsigned __int32 max_events = (unsigned __int32)(parameter[56]+0.1);
		if(Ueber->eventswritten > max_events && max_events > 0) {
			Ueber->start_new_root_file = true;
			Ueber->eventswritten = 0;
		}
	}

	return 1;
}
  K�  \   I D S _ S O U R C E _ C P P   D E T E C T O R _ S T U F F _ C P P       0         #include "OS_Version.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TApplication.h"
#include "TFile.h"


#include <math.h>
#include "rootstuff.h"
#include "Histo.h"
#include "TF1.h"
//#include "TMinuit.h"

#include "Ueberstruct.h"
//#include "ADC_analysis.h"
#include "../ADC/ADC_meta.h"

#ifndef dont_use_MFC
	#include "parallel_sort_class.h"
#endif




std::string add(std::string a, std::string b)
{
	a.append(b);
	return a;
}

std::string add(std::string a, std::string b, std::string c)
{
	a.append(b);
	a.append(c);
	return a;
}







void get_results_after_sorting(const char * name, __int32 param_offset,__int32 hist_offset, sort_class * sorter,Det_struct * det, double  parameter[], Histo * Hist, Ueberstruct * Ueber)
{
	std::string det_name = name;
	if ((parameter[param_offset] > 0.5) && (parameter[param_offset] < 2.5)) {
		bool isDLD = true;
		if((parameter[param_offset] > 1.5) && (parameter[param_offset] < 2.5)) isDLD = false;
		
		// write results to detector struct.
		//
		// The sort_and_write_NTuple function will be able to access for each hit:
		//		x position [mm]
		//		y position [mm]
		//		MCP time [ns]
		//		reconstruction method
		//		number of hits on detector
		
		for (__int32 i=0;i<det->number_of_reconstructed_hits;++i) {
			det->x[i] = sorter->output_hit_array[i]->x; // xxx
			det->y[i] = sorter->output_hit_array[i]->y;
			det->method[i] = sorter->output_hit_array[i]->method;
			det->time[i] = sorter->output_hit_array[i]->time; 
		}

		if (!Ueber->fast_mode) {

			if (det->number_of_reconstructed_hits > 0) {
				if (Hist->is_defined(hist_offset+3))
					Hist->fill2(hist_offset+3,det->x[0],det->y[0],1.);
				else
					Hist->fill2(hist_offset+3,add(det_name,std::string("_xy_mm_sorted")),det->x[0],det->y[0],1.,"",200,-100.,100.,"",200,-100.,+100,"",det_name.c_str()); 
			}

			if (det->number_of_reconstructed_hits > 0) {
				if (Hist->is_defined(hist_offset+0)) Hist->fill1(hist_offset+0,det->method[0]); else Hist->fill1(hist_offset+0,add(det_name,std::string("_method_hit_1")),det->method[0],1.,add(std::string("method_hit_1 "),det_name,std::string(" method_hit_1")),25,-0.5,24.5,add(det_name,std::string("_method_hit_1")),det_name.c_str());
			}
			if (det->number_of_reconstructed_hits > 1) {
				if (Hist->is_defined(hist_offset+1)) Hist->fill1(hist_offset+1,det->method[1]); else Hist->fill1(hist_offset+1,add(det_name,std::string("_method_hit_2")),det->method[1],1.,add(std::string("method_hit_2 "),det_name,std::string(" method_hit_2")),25,-0.5,24.5,add(det_name,std::string("_method_hit_2")),det_name.c_str());
			}
		}

	}
}








void detector_stuff_before_sorting(const char * name,__int32 param_offset,__int32 hist_offset,sort_class * sorter,Det_struct * det, double  parameter[], Histo * Hist, Ueberstruct * Ueber)
{
	double * tdc_ns = Ueber->tdc_ns;
	__int32 * cnt = Ueber->cnt;
	std::string det_name = name;
	__int32 hist_index_counter = -1;

	if ((parameter[param_offset] < 0.5) || (parameter[param_offset] > 2.5)) return;
	
	double u_ns = 10000.;
	double v_ns = 10000.;
	double w_ns = 10000.;
	double sumu_ns = 5000.;
	double sumv_ns = 5000.;
	double sumw_ns = 5000.;
	double Xuv_mm = 10000.; 
	double Xuw_mm = 10000.;
	double Xvw_mm = 10000.; 
	double Yuv_mm = 10000.; 
	double Yuw_mm = 10000.; 
	double Yvw_mm = 10000.;

	bool isDLD = true;
	if((parameter[param_offset] > 1.5) && (parameter[param_offset] < 2.5)) isDLD = false;

	// shift time sums to 0ns
	if(isDLD) {
		sorter->shift_sums(+1,det->offset_timesum_U,det->offset_timesum_V);
	} else {
		sorter->shift_sums(+1,det->offset_timesum_U,det->offset_timesum_V,det->offset_timesum_W);
	}

	// shift layer w so that all 2 center lines of the layers meet in one point
	if(!isDLD) {
		sorter->shift_layer_w(+1,det->hex_offset_w);
	}

	// shift all layers so that the position picture is centered around X=0,Y=0
	sorter->shift_position_origin(+1,det->center_X,det->center_Y);
	// plot number of hits before reconstruction
	if (!Ueber->fast_mode) {
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cmcp]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_mcp_hit")),cnt[sorter->Cmcp],1.,add(std::string("Num hits "),det_name,std::string(" MCP")),34,-0.75,16.25,add(det_name,std::string("_mcp")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cu1]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_u1_hit")),cnt[sorter->Cu1],1.,add(std::string("Num hits "),det_name,std::string(" u1")),34,-0.75,16.25,add(det_name,std::string("_u1")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cu2]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_u2_hit")),cnt[sorter->Cu2],1.,add(std::string("Num hits "),det_name,std::string(" u2")),34,-0.75,16.25,add(det_name,std::string("_u2")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cv1]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_v1_hit")),cnt[sorter->Cv1],1.,add(std::string("Num hits "),det_name,std::string(" v1")),34,-0.75,16.25,add(det_name,std::string("_v1")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cv2]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_v2_hit")),cnt[sorter->Cv2],1.,add(std::string("Num hits "),det_name,std::string(" v2")),34,-0.75,16.25,add(det_name,std::string("_v2")),det_name.c_str());
		if (!isDLD) {		
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cw1]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_w1_hit")),cnt[sorter->Cw1],1.,add(std::string("Num hits "),det_name,std::string(" w1")),34,-0.75,16.25,add(det_name,std::string("_w1")),det_name.c_str());
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,cnt[sorter->Cw2]); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_w2_hit")),cnt[sorter->Cw2],1.,add(std::string("Num hits "),det_name,std::string(" w2")),34,-0.75,16.25,add(det_name,std::string("_w2")),det_name.c_str());
		} else hist_index_counter += 2;
	} else hist_index_counter += 7;

	double mcp_signal = 0.;
	if (sorter->use_MCP) {
		if (cnt[sorter->Cmcp]>0) mcp_signal = tdc_ns[sorter->Cmcp*NUM_IONS+0]; else mcp_signal = -1.e100;
	}

	// calculate time sum u-layer
	if (sorter->use_u1 && sorter->use_u2) {
		if (cnt[sorter->Cu1]>0 && cnt[sorter->Cu2]>0) {
			u_ns	= sorter->signal_corrector->correct_pos(tdc_ns[sorter->Cu1*NUM_IONS+0],tdc_ns[sorter->Cu2*NUM_IONS+0],0);
			sumu_ns = sorter->signal_corrector->correct_sum(tdc_ns[sorter->Cu1*NUM_IONS+0],tdc_ns[sorter->Cu2*NUM_IONS+0],0) - 2.*mcp_signal;
		}
	}

	// calculate time sum v-layer
	if (sorter->use_v1 && sorter->use_v2) {
		if (cnt[sorter->Cv1]>0 && cnt[sorter->Cv2]>0) {
			v_ns	= sorter->signal_corrector->correct_pos(tdc_ns[sorter->Cv1*NUM_IONS+0],tdc_ns[sorter->Cv2*NUM_IONS+0],1);
			sumv_ns = sorter->signal_corrector->correct_sum(tdc_ns[sorter->Cv1*NUM_IONS+0],tdc_ns[sorter->Cv2*NUM_IONS+0],1) - 2.*mcp_signal;
		}
	}

	// calculate time sum w-layer
	if (!isDLD) {
		if (sorter->use_w1 && sorter->use_w2) {
			if (cnt[sorter->Cw1]>0 && cnt[sorter->Cw2]>0) {
				w_ns	= sorter->signal_corrector->correct_pos(tdc_ns[sorter->Cw1*NUM_IONS+0],tdc_ns[sorter->Cw2*NUM_IONS+0],2);
				sumw_ns = sorter->signal_corrector->correct_sum(tdc_ns[sorter->Cw1*NUM_IONS+0],tdc_ns[sorter->Cw2*NUM_IONS+0],2) - 2.*mcp_signal;
			}
		}
	}

	// calculate position on u and v from raw data
	++hist_index_counter;
	if (fabs(sumu_ns) < 10.) {
		det->sumu_watchdog->track(sumu_ns);
		if (!det->auto_calibration) {
			if (!Ueber->fast_mode) {
				if (det->sumu_watchdog->this_is_the_first_run && det->sumu_watchdog->number_of_resets == 1) {
					if (Hist->is_defined(hist_offset+hist_index_counter)) Hist->fill1(hist_offset+hist_index_counter,sumu_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_ns_offset_for_peak_tracker")),sumu_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum offset_for_peak_tracker"),det_name,std::string(" sum_u [ns]")),det_name.c_str());
				}
			}
			if (det->use_sum_tracker_layer_u && det->sumu_watchdog->is_new_result_available()) det->offset_timesum_U -= det->sumu_watchdog->last_return_value;
		}
	}
	++hist_index_counter;
	if (fabs(sumv_ns) < 10.) {
		det->sumv_watchdog->track(sumv_ns);
		if (!det->auto_calibration) {
			if (!Ueber->fast_mode) {
				if (det->sumv_watchdog->this_is_the_first_run && det->sumv_watchdog->number_of_resets == 1) {
					if (Hist->is_defined(hist_offset+hist_index_counter)) Hist->fill1(hist_offset+hist_index_counter,sumv_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumv_ns_offset_for_peak_tracker")),sumv_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum offset_for_peak_tracker"),det_name,std::string(" sum_v [ns]")),det_name.c_str());
				}
			}
				
			if (det->use_sum_tracker_layer_v && det->sumv_watchdog->is_new_result_available()) det->offset_timesum_V -= det->sumv_watchdog->last_return_value;
		}
	}
	++hist_index_counter;
	if (!isDLD && fabs(sumw_ns) < 10.) {
		det->sumw_watchdog->track(sumw_ns);
		if (!det->auto_calibration) {
			if (!Ueber->fast_mode) {
				if (det->sumw_watchdog->this_is_the_first_run && det->sumw_watchdog->number_of_resets == 1) {
					if (Hist->is_defined(hist_offset+hist_index_counter)) Hist->fill1(hist_offset+hist_index_counter,sumw_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumw_ns_offset_for_peak_tracker")),sumw_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum offset_for_peak_tracker"),det_name,std::string(" sum_w [ns]")),det_name.c_str());
				}
			}
			if (det->use_sum_tracker_layer_w && det->sumw_watchdog->is_new_result_available()) det->offset_timesum_W -= det->sumw_watchdog->last_return_value;
		}
	}

	// calculate positions using different layers [mm]
	if(isDLD) {
		Xuv_mm	=	u_ns * sorter->fu;
		Yuv_mm	=	v_ns * sorter->fv;
	} else {
		Xuv_mm	=	u_ns * sorter->fu;
		Yuv_mm	=	(u_ns * sorter->fu - 2.*v_ns * sorter->fv) / sqrt(3.);
		Xuw_mm	=	Xuv_mm;
		Yuw_mm	=	(2.*w_ns * sorter->fw  -  u_ns * sorter->fu) / sqrt(3.);
		Xvw_mm	=	(v_ns * sorter->fv  +  w_ns * sorter->fw);
		Yvw_mm	=	(w_ns * sorter->fw - v_ns * sorter->fv) / sqrt(3.);
	}

	double detector_NL_map_fill = 0.;
	double detector_resolution_map_fill = 0.;
	if (!isDLD && det->sorter->scalefactors_calibrator) {
		det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_up_to_now = -1.;
		det->sorter->scalefactors_calibrator->detector_map_devi_value_up_to_now = -1.e10;
	}

	if((!Ueber->fast_mode) || det->auto_calibration) {
		// create map of disaster
		double sumu_ns = -1.e60;
		double sumv_ns = -2.e60;
		double sumw_ns = -3.e60;
		__int32 * cnt = Ueber->cnt;
		if (sorter->use_u1 && sorter->use_u2) {
			if (cnt[sorter->Cu1]>0 && cnt[sorter->Cu2]>0) sumu_ns = tdc_ns[sorter->Cu1*NUM_IONS+0] + tdc_ns[sorter->Cu2*NUM_IONS+0];
		}
		if (sorter->use_v1 && sorter->use_v2) {
			if (cnt[sorter->Cv1]>0 && cnt[sorter->Cv2]>0) sumv_ns = tdc_ns[sorter->Cv1*NUM_IONS+0] + tdc_ns[sorter->Cv2*NUM_IONS+0];
		}
		if (!isDLD) {
			if (sorter->use_w1 && sorter->use_w2) {
				if (cnt[sorter->Cw1]>0 && cnt[sorter->Cw2]>0) sumw_ns = tdc_ns[sorter->Cw1*NUM_IONS+0] + tdc_ns[sorter->Cw2*NUM_IONS+0];
			}
		}
		if (sorter->use_MCP) {
			double mcp_signal = -1.e100;
			if (cnt[sorter->Cmcp]>0) mcp_signal = tdc_ns[sorter->Cmcp*NUM_IONS+0];
			sumu_ns -= 2.*mcp_signal;
			sumv_ns -= 2.*mcp_signal;
			if (!isDLD) sumw_ns -= 2.*mcp_signal;
		}

		
		double tu1_ns = 1.e50;
		double tu2_ns = 1.e60;
		double tv1_ns = 1.e70;
		double tv2_ns = 1.e80;
		double tw1_ns = 1.e90;
		double tw2_ns = 1.e110;
		if (sorter->use_u1) if (cnt[sorter->Cu1] > 0) tu1_ns = tdc_ns[sorter->Cu1*NUM_IONS+0];
		if (sorter->use_u2) if (cnt[sorter->Cu2] > 0) tu2_ns = tdc_ns[sorter->Cu2*NUM_IONS+0];
		if (sorter->use_v1) if (cnt[sorter->Cv1] > 0) tv1_ns = tdc_ns[sorter->Cv1*NUM_IONS+0];
		if (sorter->use_v2) if (cnt[sorter->Cv2] > 0) tv2_ns = tdc_ns[sorter->Cv2*NUM_IONS+0];
		if (!isDLD) {
			if (sorter->use_w1) if (cnt[sorter->Cw1] > 0) tw1_ns = tdc_ns[sorter->Cw1*NUM_IONS+0];
			if (sorter->use_w1) if (cnt[sorter->Cw2] > 0) tw2_ns = tdc_ns[sorter->Cw2*NUM_IONS+0];
		}

		double u_ns = tu1_ns - tu2_ns;
		double v_ns = tv1_ns - tv2_ns;
		double w_ns = tw1_ns - tw2_ns;

		if ((!det->use_reconstruction) && (!det->auto_calibration) && sorter->use_pos_correction) {
			if (sorter->use_u1 && sorter->use_v2) u_ns = sorter->signal_corrector->correct_pos(tu1_ns,tu2_ns,0);
			if (sorter->use_v1 && sorter->use_v2) v_ns = sorter->signal_corrector->correct_pos(tv1_ns,tv2_ns,1);
			if (!isDLD) {
				if (sorter->use_w1 && sorter->use_w2) w_ns = sorter->signal_corrector->correct_pos(tw1_ns,tw2_ns,2);
			}
		}

		if (fabs(sumu_ns) > 10.) u_ns = -1.e90;
		if (fabs(sumv_ns) > 10.) v_ns = -1.e92;
		if ((!isDLD) && fabs(sumw_ns) > 10.) w_ns = -1.e93;

		bool feed = false;
		if (!det->use_reconstruction) feed = true;
		if (det->use_reconstruction && det->number_of_reconstructed_hits>0) {
			if (det->method[0] == 0) feed = true;
		}
		if (feed) {
			det->sorter->feed_calibration_data(det->auto_calibration, det->hex_offset_w);
			if (!isDLD && det->sorter->scalefactors_calibrator) {
				if (det->auto_calibration == 1) { // == 1 is necessary here
					if (det->sorter->scalefactors_calibrator->map_is_full_enough()) {
						++det->auto_calibration;
						printf("\nCalibration finished on detector \"%s\"\n",name);
						Ueber->stop_reading_files = true;
					}
				}
				detector_NL_map_fill = det->sorter->scalefactors_calibrator->detector_map_devi_fill;
				detector_resolution_map_fill = det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_fill;
			}
		}
		
	}


	// plot histograms
	if (!Ueber->fast_mode) {
		double Yuv_ns,Yuw_ns,Yvw_ns,Xvw_ns;
		if(isDLD) {
			Yuv_ns	=	v_ns;
		} else {
			Yuv_ns	=	(u_ns - 2.*v_ns) / sqrt(3.);
			Yuw_ns	=	(2.*w_ns  -  u_ns) / sqrt(3.);
			Xvw_ns	=	v_ns + w_ns;
			Yvw_ns	=	(w_ns - v_ns) / sqrt(3.);
		}
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,u_ns,Yuv_ns); else if (!isDLD) Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uv_ns")),u_ns,Yuv_ns,1.,"U-V-Layer without shifts etc. [ns]",350,-175.,175.,"x_uv [ns]",350,-175.,175.,"y_uv [ns]",det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,u_ns,Yuw_ns); else if (!isDLD) Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uw_ns")),u_ns,Yuw_ns,1.,"U-W-Layer without shifts etc. [ns]",350,-175.,175.,"x_uw [ns]",350,-175.,175.,"y_uw [ns]",det_name.c_str());
		if (!isDLD) {
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xvw_ns,Yvw_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_vw_ns")),Xvw_ns,Yvw_ns,1.,"V-W-Layer without shifts etc. [ns]",350,-175.,175.,"x_vw [ns]",350,-175.,175.,"y_vw [ns]",det_name.c_str());
		} else ++hist_index_counter;

		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumu_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_ns")),sumu_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_u [ns]")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumv_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumv_ns")),sumv_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_v [ns]")),det_name.c_str());
		if (!isDLD) {
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumw_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumw_ns")),sumw_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_w [ns]")),det_name.c_str());
		} else ++hist_index_counter;
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,sumu_ns,sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sum_uv_ns")),sumu_ns,sumv_ns,1.,"",150,-150.,150.,"sum u [ns]",150,-150.,150.,"sum v [ns]",det_name.c_str());

		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,u_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_u_ns")),u_ns,1.,"",500,-250.,250.,add(std::string("u-layer "),det_name,std::string(" [ns]")),det_name.c_str());
		
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,v_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_v_ns")),v_ns,1.,"",500,-250.,250.,add(std::string("v-layer "),det_name,std::string(" [ns]")),det_name.c_str());
		if (!isDLD) {
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,w_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_w_ns")),w_ns,1.,"",500,-250.,250.,add(std::string("w-layer "),det_name,std::string(" [ns]")),det_name.c_str());
		} else ++hist_index_counter;

		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,u_ns,sumu_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_u_ns")),u_ns,sumu_ns,1.,"",400,-100.,100.,add(std::string("u-layer "),det_name,std::string(" [ns]")),160,-5.,5.,add(std::string("Time sum "),det_name,std::string(" u [ns]")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,v_ns,sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumv_v_ns")),v_ns,sumv_ns,1.,"",400,-100.,100.,add(std::string("v-layer "),det_name,std::string(" [ns]")),160,-5.,5.,add(std::string("Time sum "),det_name,std::string(" v [ns]")),det_name.c_str());
		if (!isDLD) {
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,w_ns,sumw_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumw_w_ns")),w_ns,sumw_ns,1.,"",400,-100.,100.,add(std::string("w-layer "),det_name,std::string(" [ns]")),160,-5.,5.,add(std::string("Time sum "),det_name,std::string(" w [ns]")),det_name.c_str());
		} else ++hist_index_counter;

		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xuv_mm,Yuv_mm); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uv_mm")),Xuv_mm,Yuv_mm,1.,"",280,-70.,70.,"Pos x UV [mm]",280,-70.,70.,"Pos y UV [mm]",det_name.c_str());

		if (!isDLD) {
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xuw_mm,Yuw_mm); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uw_mm")),Xuw_mm,Yuw_mm,1.,"",280,-70.,70.,"Pos x UW [mm]",280,-70.,70.,"Pos y UW [mm]",det_name.c_str());
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xvw_mm,Yvw_mm); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_vw_mm")),Xvw_mm,Yvw_mm,1.,"",280,-70.,70.,"Pos x VW [mm]",280,-70.,70.,"Pos y VW [mm]",det_name.c_str());
			double x=-1.e100;
			double y=-1.e100;
			if (fabs(sumu_ns) < 10. && fabs(sumv_ns) < 10.) {
				x = Xuv_mm;	y = Yuv_mm;
			} else if (fabs(sumu_ns) < 10. && fabs(sumw_ns) < 10.) {
				x = Xuw_mm;	y = Yuw_mm;
			} else if (fabs(sumv_ns) < 10. && fabs(sumw_ns) < 10.) {
				x = Xvw_mm;	y = Yvw_mm;
			}
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,x,y); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_xy_mm")),x,y,1.,"xy raw [mm]",160,-80.,80,"x [mm]",160,-80.,80.,"y [mm]",det_name.c_str());
		} else hist_index_counter += 3;

		if (det->sorter) {
			if (det->sorter->scalefactors_calibrator) {
				__int32 n = __int32(det->sorter->scalefactors_calibrator->detector_map_size * 0.5 + 0.1);
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) {
					//Hist->fill2(hist_offset+hist_index_counter,det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,detector_NL_map_fill);
					if (det->sorter->scalefactors_calibrator->detector_map_devi_value_up_to_now > -1.e9) Hist->SetBinContentAt(hist_offset+hist_index_counter,det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,det->sorter->scalefactors_calibrator->detector_map_devi_value_up_to_now);
				} else {
					//Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_non_linearity")),det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,detector_NL_map_fill,"",2*n,-n-0.5,+n-0.5,"",2*n,-n-0.5,+n-0.5,"",det_name.c_str()); 
					if (det->sorter->scalefactors_calibrator->detector_map_devi_value_up_to_now > -1.e9) Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_non_linearity")),det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,det->sorter->scalefactors_calibrator->detector_map_devi_value_up_to_now,"",2*n,-n-0.5,+n-0.5,"",2*n,-n-0.5,+n-0.5,"",det_name.c_str()); 
				}
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) {
					//Hist->fill2(hist_offset+hist_index_counter,det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,detector_resolution_map_fill);
					if (det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_up_to_now >= 0.) Hist->SetBinContentAt(hist_offset+hist_index_counter,det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_up_to_now);
				} else {
					//Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_resolution_FWHM")),det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,detector_resolution_map_fill,"",2*n,-n-0.5,+n-0.5,"",2*n,-n-0.5,+n-0.5,"",det_name.c_str()); 
					if (det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_up_to_now >= 0.) Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_resolution_FWHM")),det->sorter->scalefactors_calibrator->binx-n,det->sorter->scalefactors_calibrator->biny-n,det->sorter->scalefactors_calibrator->detector_map_resol_FWHM_up_to_now,"",2*n,-n-0.5,+n-0.5,"",2*n,-n-0.5,+n-0.5,"",det_name.c_str()); 
				}
			}
		}
	}


	ADC_meta * adc = (ADC_meta*)(Ueber->adc_meta);
if (adc) {
	// tbb tests



	// tobuaer
	if (det->number_of_reconstructed_hits > 1) {
		if (Hist->is_defined(hist_offset + (++hist_index_counter))) Hist->fill2(hist_offset + hist_index_counter, sumu_ns, sumv_ns); else Hist->fill2(hist_offset + hist_index_counter, add(det_name, std::string("_sumu_sumv_ns_2x2")), sumu_ns, sumv_ns, 1., "", 400, -4., 4., add(det_name, std::string("sumu [ns]")), 400, -4., 4., add(det_name, std::string(" sumv [ns]")), det_name.c_str());
	}else hist_index_counter += 1;

	if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,sumu_ns,sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_2x2")),sumu_ns,sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());	
	if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_ROT_2x2mm")),0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());

	if(adc->m_event->channels[2] && adc->m_event->channels[2]->NbrPulses && adc->m_event->channels[2]->Pulses[0].NbrPeaks){
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch2")),adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch0")),adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
		
		if(adc->m_event->channels[2]->Pulses[0].Peaks[0].Height >= 120. && adc->m_event->channels[2]->Pulses[0].Peaks[0].Height <= 140.){	
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch2_hight")),adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch0_hight")),adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
		}
	}


	if(Xuv_mm < -5. && Xuv_mm > - 8 && Yuv_mm < 8. && Yuv_mm > 5.){
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumu_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_ns_tbb")),sumu_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_u [ns]")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumv_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumv_ns_tbb")),sumv_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_v [ns]")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumw_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumw_ns_tbb")),sumw_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_w [ns]")),det_name.c_str());
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xuv_mm,Yuv_mm); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uv_mm_tbb")),Xuv_mm,Yuv_mm,1.,"",280,-70.,70.,"Pos x UV [mm]",280,-70.,70.,"Pos y UV [mm]",det_name.c_str());

		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,sumu_ns,sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_tbb")),sumu_ns,sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());	
		if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_ROT_tbb")),0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());

	
		if(adc->m_event->channels[2] && adc->m_event->channels[2]->NbrPulses && adc->m_event->channels[2]->Pulses[0].NbrPeaks){

			
			// poshalf cfd .. 
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch2_3x3")),adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
			if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch0_3x3")),adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
			
			
			if(adc->m_event->channels[2]->Pulses[0].Peaks[0].Height >= 120. && adc->m_event->channels[2]->Pulses[0].Peaks[0].Height <= 140.){		
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch2_3x3_hight")),adc->m_event->channels[2]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[2]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_ch0_3x3_hight")),adc->m_event->channels[0]->Pulses[0].Peaks[0].PosHalfLeft,adc->m_event->channels[0]->Pulses[0].Peaks[0].Cfd,1.,"",400,0.,40.,"PosHalfLeft [stepsize]",400,0.,40.,"cfd [stepsize]","PosHalfLeft_vs_Cfd");

				
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumu_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_ns_tbb_hight")),sumu_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_u [ns]")),det_name.c_str());
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumv_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumv_ns_tbb_hight")),sumv_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_v [ns]")),det_name.c_str());
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill1(hist_offset+hist_index_counter,sumw_ns); else Hist->fill1(hist_offset+hist_index_counter,add(det_name,std::string("_sumw_ns_tbb_hight")),sumw_ns,1.,"",20*400,-200.,200.,add(std::string("Time sum "),det_name,std::string(" sum_w [ns]")),det_name.c_str());
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,Xuv_mm,Yuv_mm); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_uv_mm_tbb_hight")),Xuv_mm,Yuv_mm,1.,"",280,-70.,70.,"Pos x UV [mm]",280,-70.,70.,"Pos y UV [mm]",det_name.c_str());

				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,sumu_ns,sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_tbb_hight")),sumu_ns,sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());	
				if (Hist->is_defined(hist_offset+(++hist_index_counter))) Hist->fill2(hist_offset+hist_index_counter,0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns); else Hist->fill2(hist_offset+hist_index_counter,add(det_name,std::string("_sumu_sumv_ns_ROT_tbb_hight")),0.7071*sumu_ns-0.7071*sumv_ns,0.7071*sumu_ns+0.7071*sumv_ns,1.,"",400,-4.,4.,add(det_name,std::string("sumu [ns]")),400,-4.,4.,add(det_name,std::string(" sumv [ns]")),det_name.c_str());
			}
			//std::cout << "\n " << adc->m_event->channels[2]->Pulses[0].Peaks[0].Height;
		}
	
		
	
	
	}
	
	






	//std::cout << "\n" << adc->m_adc_info[1].delay;
	//std::cout << "\n" << Ueber->LMF_input->DAQ_ID;
	// tbb tests
}







	// everything above happened without changing the TDC array! 
}










///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int detector_stuff(__int64 eventcounter, Ueberstruct * Ueber, double  parameter[])
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	Histo * Hist = Ueber->Hist;
	const double pi = 3.14159265359;

//--------------------------------------------------------------------------------------------------------------------------
//  Detector calibration and stuff using Achim's sort routine...
//--------------------------------------------------------------------------------------------------------------------------

// use raw data only to create histograms

	detector_stuff_before_sorting("proj",100,10100,Ueber->current_proj_Sorter,Ueber->proj, parameter, Hist, Ueber);
	detector_stuff_before_sorting("rec" ,200,10200,Ueber->current_rec_Sorter ,Ueber->rec , parameter, Hist, Ueber);
	detector_stuff_before_sorting("elec",300,10300,Ueber->current_elec_Sorter,Ueber->elec, parameter, Hist, Ueber);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// do sorting !
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Ueber->proj->number_of_reconstructed_hits = 0;
	Ueber->rec->number_of_reconstructed_hits = 0;
	Ueber->elec->number_of_reconstructed_hits = 0;

	if (parameter[10] < 0.5) { // 0 means: do not use multi-core processing

		if (Ueber->current_proj_Sorter) {
			if (Ueber->proj->use_reconstruction) Ueber->current_proj_Sorter->sort();
			else Ueber->current_proj_Sorter->run_without_sorting();

			Ueber->proj->number_of_reconstructed_hits = Ueber->current_proj_Sorter->output_hit_array_counter < NUM_IONS ? Ueber->current_proj_Sorter->output_hit_array_counter : NUM_IONS;
			}
		if (Ueber->current_rec_Sorter) {
			if (Ueber->rec->use_reconstruction) Ueber->current_rec_Sorter->sort();
			else Ueber->current_rec_Sorter->run_without_sorting();

			Ueber->rec->number_of_reconstructed_hits = Ueber->current_rec_Sorter->output_hit_array_counter < NUM_IONS ? Ueber->current_rec_Sorter->output_hit_array_counter : NUM_IONS;

			}
		if (Ueber->current_elec_Sorter) {
			if (Ueber->elec->use_reconstruction) Ueber->current_elec_Sorter->sort();
			else Ueber->current_elec_Sorter->run_without_sorting();

			Ueber->elec->number_of_reconstructed_hits = Ueber->current_elec_Sorter->output_hit_array_counter < NUM_IONS ? Ueber->current_elec_Sorter->output_hit_array_counter : NUM_IONS;
		}

	} else {   // use multi-core processing?
#ifndef dont_use_MFC
		parallel_sort_class * parallel = (parallel_sort_class*)Ueber->parallel_sorter;

		__int32 ID = parallel->next_thread_ID_to_be_sorted;
		parallel->do_sort_with_thread_ID(ID);
		parallel->next_thread_ID_to_be_sorted++;
		if (parallel->next_thread_ID_to_be_sorted == parallel->number_of_threads) parallel->next_thread_ID_to_be_sorted = 0;

		if (parallel->filled_threads < parallel->number_of_threads) return 0;

		ID = parallel->next_thread_ID_to_be_sorted;
		parallel->wait_for_thread_with_ID(ID);

		Ueber->tdc_ns = parallel->events[ID]->tdc_ns;
		Ueber->cnt = parallel->events[ID]->cnt;
		Ueber->current_proj_Sorter = parallel->events[ID]->sorter_proj;
		Ueber->current_rec_Sorter = parallel->events[ID]->sorter_rec;
		Ueber->current_elec_Sorter = parallel->events[ID]->sorter_elec;

	
		if (Ueber->current_proj_Sorter)	{
			Ueber->proj->number_of_reconstructed_hits = Ueber->current_proj_Sorter->output_hit_array_counter < NUM_IONS ? Ueber->current_proj_Sorter->output_hit_array_counter : NUM_IONS;
		}
		if (Ueber->current_rec_Sorter)	{
			Ueber->rec->number_of_reconstructed_hits  = Ueber->current_rec_Sorter->output_hit_array_counter  < NUM_IONS ? Ueber->current_rec_Sorter->output_hit_array_counter : NUM_IONS; 
		}
		if (Ueber->current_elec_Sorter)	{
			Ueber->elec->number_of_reconstructed_hits = Ueber->current_elec_Sorter->output_hit_array_counter < NUM_IONS ? Ueber->current_elec_Sorter->output_hit_array_counter : NUM_IONS;
		}

		if (Ueber->elec->number_of_reconstructed_hits != Ueber->current_elec_Sorter->output_hit_array_counter && Ueber->current_elec_Sorter->output_hit_array_counter < NUM_IONS) {
			int xxx = 0;
		}
#else
		printf("\nplease disable multi-core processing in the config file.\n");
		printf("(not supported if lmf2root was not compiled with MFC-libs.)\n");
		Ueber->stop_reading_files = true;
#endif		
	}

	if (Ueber->proj->number_of_reconstructed_hits < 0) {
		printf("proj: Some parameters for the reconstruction routine\nare not correctly initialized\n");
		Ueber->proj->number_of_reconstructed_hits = 0;
	}
	if (Ueber->rec->number_of_reconstructed_hits < 0) {
		printf("rec: Some parameters for the reconstruction routine\nare not correctly initialized\n");
		Ueber->rec->number_of_reconstructed_hits = 0;
	}
	if (Ueber->elec->number_of_reconstructed_hits < 0) {
		printf("elec: Some parameters for the reconstruction routine\nare not correctly initialized\n");
		Ueber->elec->number_of_reconstructed_hits = 0;
	}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// get results
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	get_results_after_sorting("proj",100,10150,Ueber->current_proj_Sorter,Ueber->proj, parameter, Hist, Ueber);
	get_results_after_sorting("rec" ,200,10250,Ueber->current_rec_Sorter ,Ueber->rec , parameter, Hist, Ueber);
	get_results_after_sorting("elec",300,10350,Ueber->current_elec_Sorter,Ueber->elec, parameter, Hist, Ueber);

	return 1;
} 
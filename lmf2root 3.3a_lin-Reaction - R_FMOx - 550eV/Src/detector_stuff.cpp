#include "OS_Version.h"

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
if (!adc) {
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
}
#ifndef CH_PRESORTER_ALREADY_INCLUDED
	#define CH_PRESORTER_ALREADY_INCLUDED
	
	#define PRESORTERVERSION (0.2)

	#include "CH_Basics.h"
	#include "../CH_Event.h"
	#include "CH_Spectrometer.h"
	#include "CH_Tof.h"

	#include <math.h>

	#include "CH_Hist.h"

	// Maximum number of presorter hists per channel.
	#define NHISTS 50

	typedef bool (CH::presorter_class::*prsp)();

	namespace CH {

		class presorter_class
		{	
			public:
				CH_event_struct evti;
				CH_event_struct evto;
				// These are presorter histos only!
				histo_handler * H;

			private:
				int Flag;
				char dir[80];
				int HistoBlock;

				spectrometer_class * spect;
				tof_calc_class * ctof;

				// array of presorter functions to run
				prsp run_prs_array[64];
				int num_of_presorters;

				// presorter results
				bool result;
				bool e_hit_valid_array[64];
				bool r_hit_valid_array[64];
				bool p_hit_valid_array[64];
				std::vector<CH_event_struct> PrsEvtOut;

				// conditions on number of hits
				int ehit_min;
				int ehit_max;
				int rhit_min;
				int rhit_max;
				int phit_min;
				int phit_max;

				// single tof variables
				double etof_min;
				double etof_max;
				double epos_center_x;
				double epos_center_y;
				double epos_radius;

				double rtof_min;
				double rtof_max;
				double rpos_center_x;
				double rpos_center_y;
				double rpos_radius;

				double ptof_min;
				double ptof_max;
				double ppos_center_x;
				double ppos_center_y;
				double ppos_radius;

				// additional pipico variables
				double mass1_amu;
				double charge1_au;
				double mass2_amu;
				double charge2_au;
				double efield1_Vcm;
				double acceleration1_mm;
				double efield2_Vcm;
				double acceleration2_mm;
				double efield3_Vcm;
				double acceleration3_mm;
				double pipico_width_ns;
				double shift_ns;

				// user presorter parameters
				double usr_par[16][32]; 
				int num_of_usr_par[32]; 

				//sumdiff presorter parameters
				int number_of_recoils_to_take;
				double rtof_sum_center;
				double rtof_sum_width;
				double rtof_diff_center;
				double rtof_diff_width;

				// polyatomic presorter variables
				int min_number_of_recoils; 
				double tofmin[16];
				double tofmax[16];
				int ions_defined;


				bool tofgood(CH_det_struct* particle, int i, double min, double max); 
				bool posgood(CH_det_struct* particle, int i, double center_x, double center_y, double radius); 

			public:

				presorter_class(int Flag, histo_handler * Hist, spectrometer_class * spect, tof_calc_class * ctof);
				~presorter_class();

				// function pointer handling
				void append_presorter_to_list(prsp prs);

				// Get presorter channel
				int Get_flag();
				
				// Set histogram block 
				void SetHistoBlock(int num);

				// run presorter
				bool Run();
				
				// retrieve results
				bool is_evt_valid();
				bool *Get_hit_valid_array_e();
				bool *Get_hit_valid_array_r();
				bool *Get_hit_valid_array_p();

				///// presorter init
				// parameters for hit conditions
				void set_hits_e(int min, int max);
				void set_hits_r(int min, int max);
				void set_hits_p(int min, int max);

				// parameters for pos and tof conditions
				void set_tof_e(double tof_min, double tof_max);
				void set_pos_e(double pos_center_x, double pos_center_y, double pos_radius);
				void set_tof_r(double tof_min, double tof_max);
				void set_pos_r(double pos_center_x, double pos_center_y, double pos_radius);
				void set_tof_p(double tof_min, double tof_max);
				void set_pos_p(double pos_center_x, double pos_center_y, double pos_radius);

				// parameters for pipico
				void set_pipico(double mass1_amu, double charge1_au, double mass2_amu, double charge2_au, double pipico_width);

				// parameters for sumdiff
				void set_sumdiff(int number_of_recoils_to_take, double rtof_sum_center, double rtof_sum_width, double rtof_diff_center, double rtof_diff_width);

				// set polyatomic parameters
				void set_polyatomic(int min_number_of_recoils);
				void polyatomic_add_ion(double tmin, double tmax);

				// set user parameters
				void add_user_parameter(int prs_num, double val);
				
				// presorter functions
				bool hits_e();
				bool hits_r();
				bool hits_p();

				bool tof_e();
				bool tof_r();
				bool tof_p();

				bool pos_e();
				bool pos_r();
				bool pos_p();

				bool pipico();

				bool sumdiff();

				bool polyatomic();

				// User defined presorters
				bool user_prs_0();
				bool user_prs_1();
				bool user_prs_2();
				bool user_prs_3();
				bool user_prs_4();
				bool user_prs_5();
				bool user_prs_6();
				bool user_prs_7();

		};	// end presorter_class
	}

#endif;
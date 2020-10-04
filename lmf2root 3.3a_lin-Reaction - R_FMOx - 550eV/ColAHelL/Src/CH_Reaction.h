#ifndef CH_REACTION_ALREADY_INCLUDED
	#define CH_REACTION_ALREADY_INCLUDED
	
	#include "CH_Basics.h"
//	#include "CH_Ranges.h"
	#include <vector>

	using namespace std;

	namespace CH
	{
		// reaction types (more to be added.. Even values are reactions including electrons!)
		#define RTYPE_ELEC					0
		#define RTYPE_ION					1
		#define RTYPE_ION_ELEC				2
		#define RTYPE_DIATOMIC				3
		#define RTYPE_DIATOMIC_ELEC			4
		#define RTYPE_POLYATOMIC			5
		#define RTYPE_POLYATOMIC_ELEC		6
		#define RTYPE_DIATOMIC_DISS			7
		#define RTYPE_DIATOMIC_DISS_ELEC	8
/*
		#define TAGTYPE_ELEC			1
		#define TAGTYPE_ION				2
		#define TAGTYPE_PROJ			3
		#define TAGTYPE_DIATOMIC		10
*/

		#define FIRST_DIATOMIC_TAG	14
		#define NUMBER_OF_TAGS		27

		// types for tagging and (in)validating of particles 
		#define RTAG_X				0
		#define RTAG_Y				1
		#define RTAG_TOF			2
		#define RTAG_PX				3
		#define RTAG_PY				4
		#define RTAG_PZ				5
		#define RTAG_PXY			6
		#define RTAG_PXZ			7
		#define RTAG_PYZ			8
		#define RTAG_P				9
		#define RTAG_PHI_DEG		10
		#define RTAG_PHIPOS			11
		#define RTAG_THETA_DEG		12
		#define RTAG_ENERGY			13

		//additional types for diatomics
		#define RTAG_PSUMX			14
		#define RTAG_PSUMY			15
		#define RTAG_PSUMZ			16
		#define RTAG_PSUM			17
		#define RTAG_PSUM_PHI_DEG	18
		#define RTAG_PSUM_THETA_DEG	19
		#define RTAG_PRELX			20
		#define RTAG_PRELY			21
		#define RTAG_PRELZ			22
		#define RTAG_PREL			23
		#define RTAG_PREL_PHI_DEG	24
		#define RTAG_PREL_THETA_DEG	25
		#define RTAG_KER			26

		// MFPAD types
		#define NONE -1;
		#define LIN 0;
		#define CIRC 1;
		#define AUGER 2;

		struct rtag_struct {
			short particle_num; // Apply condition to particle PARTICLE_NUM. (-1: all)
			int type;			// Reaction type (DIATOMIC, SINGLEION, ELEC etc.)
			int tag_type;		// Variable which is restricted (X, Y, KER etc.)
			int channel;
			double min;
			double max;
			bool invalidate;
			bool validate;
			bool tag;
			bool negate;
		};

		struct expected_range_struct {
			int rtype;
			int rtag;
			double min;
			double max;
		};

		struct MFPAD_cond_struct {
			int type;
			double rmin;
			double rmax;
			double emin;
			double emax;
		};

		struct LFPAD_cond_struct {
			int type;
			double emin;
			double emax;
		};

		struct reaction_struct {
			char name[180];
			int channel;
			int type;
			int ID;

			int diatomic_prel;

			int reac_num; // this is reaction definition number "reac_num" found in the config file... 

			vector <double> mass;
			vector <double> charge;
			vector <double> t_mean;

			double t_width[16]; //for Polyatomic
			int number_of_ions; //for Polyatomic
			int use_ion_matrix; //for Polyatomic 
			double max_target_value; //for Polyatomic
			double ambiguity_parameter; //for Polyatomic
			std::string interpolation_filename[16];		
			std::vector<std::vector<std::vector<double>>>    interpolation_points;
			std::vector<std::vector<double>>    interpolation_adjustment_variables;
			bool incomplete; // for Polyatomic

			//vector <rtag_struct> tag;
			rtag_struct *tag[16];
			int num_tags;
			expected_range_struct expctd[16];

			cor_param *e_fac;
			mom_cor_param *e_mom_fac;

			cor_param *r_fac;
			mom_cor_param *r_mom_fac;

			cor_param *p_fac;
			mom_cor_param *p_mom_fac;

			bool randomize_ions;
			int rnd_array[16];
			short rnd_count;

			MFPAD_cond_struct *MF_cond;
			LFPAD_cond_struct *LF_cond;

			int Dalitz_array[3];
			int Newton_array[3];

			bool two_elec;
			double e_master_min;
			double e_master_max;

			bool photon_scan;
			int ph_scan_channel;
			double ph_min;
			double ph_max;
			double ph_step;
		};
	}
#endif;
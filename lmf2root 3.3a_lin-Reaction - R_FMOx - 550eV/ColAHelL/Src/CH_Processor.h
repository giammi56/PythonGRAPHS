#pragma once

#include "CH_Basics.h"
#include "CH_Reaction.h"
#include "..\CH_Event.h"
#include "CH_Electron.h"
#include "CH_Ion.h"
#include "CH_Diatomic.h"
#include "CH_Polyatomic.h"

#define MAX_PARTICLES_TO_PROCESS 64
#define MAX_NUM_REACTIONS 64


namespace CH
{
	class processor_class {  
		
	private:

		particle_raw r_raw[MAX_PARTICLES_TO_PROCESS];
		particle_raw e_raw[MAX_PARTICLES_TO_PROCESS];
		particle_raw p_raw[MAX_PARTICLES_TO_PROCESS];
		int cur_reac_num;

	public:

		bool at_least_one_r_valid;
		bool at_least_one_e_valid;

		void new_reaction();
//		void reset_reaction(int reac_num);
		void new_rtag(__int32 edit_reaction_num);

		reaction_struct *reaction_list[MAX_NUM_REACTIONS];
		int number_of_reactions_defined;

		spectrometer_class *Spect;
		
		bool cut_wiggle;
		double cut_wiggle_width;

		// These are the default parameters applied independent of the reaction
		//--------------------------------
		cor_param *e_fac;
		mom_cor_param *e_mom_fac;

		cor_param *r_fac;
		mom_cor_param *r_mom_fac;

		cor_param *p_fac;
		mom_cor_param *p_mom_fac;
		//--------------------------------

		processor_class();
		~processor_class();

		int find_reaction(int channel);
		int find_reaction(int channel, int already_processed);
		int EvtBelongsToReactions(int chan);
		reaction_struct * get_current_reaction();

		bool process_all(CH_event_struct *evti, electron_class ** el, ion_class ** rec, diatomic_class * molec, polyatomic_class * big_molec, int which_reaction_def);
		bool process_electrons(reaction_struct *cur_reaction, CH_event_struct *evti, electron_class ** el);
		bool process_ions(reaction_struct *cur_reaction, CH_event_struct *evti, ion_class ** rec, diatomic_class * molec, polyatomic_class * big_molec);

	}; // end processor class
}






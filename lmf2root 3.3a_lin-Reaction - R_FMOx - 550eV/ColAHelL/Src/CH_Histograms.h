#pragma once

#pragma warning(disable : 4996)

#include "CH_Basics.h"
#include "CH_Reaction.h"
#include "..\CH_Event.h"
#include "CH_Electron.h"
#include "CH_Ion.h"
#include "CH_Diatomic.h"
#include "CH_Polyatomic.h"
#include "CH_Hist.h"
#include "Math.h"

// Defines how many histograms per channel and ion (or electron) can be used...
// Smaller values speed up histogram writing when closing the file.. 
#define HISTS_PER_CHANNEL 200

namespace CH
{
	class histograms_class { 

	private:
		reaction_struct *cur_reaction;
		ion_class * r[64];
		
		char CH_master_folder[180];
		bool use_master_folder;

	public:
		histograms_class();
		~histograms_class();

		void plot(reaction_struct *reac, electron_class ** e, ion_class ** r, diatomic_class * mol, polyatomic_class * big_mol, CH_event_struct *evt);

		void plot_electrons(reaction_struct *cur_reaction, electron_class ** e, int ehit, double *scan_val);
		void plot_ions(reaction_struct *cur_reaction, ion_class ** r, int rhit, electron_class ** e, int ehit, double *scan_val);
		void plot_diatomic(reaction_struct *cur_reaction, diatomic_class * mol, electron_class ** e, int ehit, double *scan_val);
		void plot_polyatomic(reaction_struct *cur_reaction, polyatomic_class * big_mol, electron_class ** e, int ehit, double *scan_val);
		
		void SetMasterFolder(const char* folder);

		void UserAnalysis(reaction_struct *reac, electron_class ** e, ion_class ** r, diatomic_class * mol, polyatomic_class * big_mol, CH_event_struct *evt);
		int MaxUserHists;

		histo_handler * Hist_e;
		histo_handler * Hist_ions;
		histo_handler * Hist_User;


		double rdet_size;
		double edet_size;
		double pdet_size;

		//these two are used to flag whether we had at least one valid ion or electron after processing... 
		bool electron_valid;
		bool ion_valid;
	};
}
#include "CH_Histograms.h"
#include "assert.h"
#include "Math.h"
#include <vector>

namespace CH
{
	histograms_class::histograms_class() 
	{
		Hist_e = new histo_handler();
		Hist_ions = new histo_handler();
		Hist_User =  new histo_handler();
		MaxUserHists = 0;

		cur_reaction = 0;

		rdet_size = 55.0;
		edet_size = 55.0;
		pdet_size = 55.0;

		this->electron_valid = false;
		this->ion_valid = false;

		this->use_master_folder = false;
	}
	
	histograms_class::~histograms_class() 
	{
		if (Hist_e) delete Hist_e;
		Hist_e = 0;
		if (Hist_ions) delete Hist_ions;
		Hist_ions = 0;
		if (Hist_User) delete Hist_User;
		Hist_User = 0;
		//if (cur_reaction) delete cur_reaction; cur_reaction = 0;
	}

	void histograms_class::SetMasterFolder(const char* folder) {
		this->use_master_folder = true;
		strncpy(this->CH_master_folder,folder,180);
		strcat(this->CH_master_folder,"/");
	}

	void histograms_class::plot(reaction_struct *reac, electron_class ** e, ion_class ** r, diatomic_class * mol, polyatomic_class * big_mol, CH_event_struct *evt) {
		
		assert(reac);
		cur_reaction = reac;
		int ehit = evt->e.num_hits;
		int rhit = evt->r.num_hits;
		// get the pointer of the vector!
		double* scan_val = 0;
		if(evt->scanval.size()>0)
			scan_val = &evt->scanval[0];

		// plot heavy stuff
		switch (cur_reaction->type) {
			
			case(RTYPE_ELEC):
				if(this->electron_valid)
					plot_electrons(cur_reaction, e, ehit, scan_val);
				break;
			case(RTYPE_ION):
				if(this->ion_valid)
					plot_ions(cur_reaction,r,rhit, e, ehit, scan_val);
				break;
			case(RTYPE_ION_ELEC):
				if(this->ion_valid && this->electron_valid) {
					plot_ions(cur_reaction,r,rhit, e, ehit, scan_val);
					plot_electrons(cur_reaction, e, ehit, scan_val);
				}
				break;
			case(RTYPE_DIATOMIC):
			case(RTYPE_DIATOMIC_DISS):
				if(mol->valid)
					plot_diatomic(cur_reaction, mol, e, ehit, scan_val);
				break;
			case(RTYPE_DIATOMIC_ELEC):
			case(RTYPE_DIATOMIC_DISS_ELEC):
				if(mol->valid && this->electron_valid) {
					plot_diatomic(cur_reaction, mol, e, ehit, scan_val);
					plot_electrons(cur_reaction, e, ehit, scan_val);
				}
				break;			
			case(RTYPE_POLYATOMIC):
				if(big_mol->valid)
					plot_polyatomic(cur_reaction, big_mol, e, ehit, scan_val);
				break;
			case(RTYPE_POLYATOMIC_ELEC):
				if(big_mol->valid && this->electron_valid) {
					plot_polyatomic(cur_reaction, big_mol, e, ehit, scan_val);
					plot_electrons(cur_reaction, e, ehit, scan_val);
				}
				break;
		}
	}
}
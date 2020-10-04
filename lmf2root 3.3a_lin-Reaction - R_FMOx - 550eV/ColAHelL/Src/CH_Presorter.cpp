#pragma warning(disable : 4996)
#include "CH_Presorter.h"
#include <math.h>
#include "CH_Tof.h"
#include "CH_Spectrometer.h"
#include "CH_FUN_Lowlevel.h"



namespace CH 
{
	presorter_class::presorter_class(int Flag, histo_handler * Hist, spectrometer_class * spect, tof_calc_class * ctof)
	{
		this->num_of_presorters = 0;

		this->Flag = Flag;
		this->H = Hist;
		this->ctof = ctof;
		this->spect = spect;
		
		for(int i=0;i<16;i++) {
			this->num_of_usr_par[i] = 0;
			for(int j=0;j<32;j++)
				this->usr_par[i][j];
		}

		sprintf_s(dir,80,"_presorter check/channel %d",Flag);
	}

	presorter_class::~presorter_class() {
	}

	void presorter_class::SetHistoBlock(int num) {
		this->HistoBlock = num;
	}

	int presorter_class::Get_flag() {
		return Flag;
	}

	void presorter_class::append_presorter_to_list(prsp prs) {
		run_prs_array[num_of_presorters] = prs;
		this->num_of_presorters++;
	}
	
	bool presorter_class::is_evt_valid() {
		return this->result;
	}

	bool *presorter_class::Get_hit_valid_array_e() {
		return this->e_hit_valid_array;
	}

	bool *presorter_class::Get_hit_valid_array_r() {
		return this->r_hit_valid_array;
	}

	bool *presorter_class::Get_hit_valid_array_p() {
		return this->p_hit_valid_array;
	}

	bool presorter_class::Run()
	{
		// Oh no, depending on the tof-calc method we need to loop over all reference particles!
		int num_ref_particles = 1;
		bool use_particle_as_ref = true;
		bool *ref_hit_valid_array; 
		
		switch(this->ctof->GetTOFType()) {
			case 0: // BM
				use_particle_as_ref = false;
				if(this->ctof->GetBMspacing()>0.1) {
					use_particle_as_ref = true;
					num_ref_particles = this->evti.e.num_hits;
					ref_hit_valid_array = e_hit_valid_array;
				} else {
					num_ref_particles = 1;
				}
				break;
			case 1: // PROJ
				use_particle_as_ref = true;
				num_ref_particles = this->evti.p.num_hits;
				ref_hit_valid_array = p_hit_valid_array;
				break;
			case 2: // ELEC
				use_particle_as_ref = true;
				num_ref_particles = this->evti.e.num_hits;
				ref_hit_valid_array = e_hit_valid_array;
				break;
			default:
				break;
		}

		bool event_valid = false;

		int ref;
		for(ref=0;ref<num_ref_particles;ref++) {
			event_valid = true;

			//prepare arrays for tagging hits
			for(int i=0;i<this->evti.e.num_hits;i++) 
				e_hit_valid_array[i] = true;
			for(int i=0;i<this->evti.r.num_hits;i++) 
				r_hit_valid_array[i] = true;
			for(int i=0;i<this->evti.p.num_hits;i++) 
				p_hit_valid_array[i] = true;

			ctof->calc(&this->evti,ref);
			
			// loop over all presorters, that are requested for this event
			for(int i=0;i<this->num_of_presorters;i++) {
				// if one of the presorters tags the event as non valid the event needs to be discarded 
				if(!(this->*run_prs_array[i])())
					event_valid = false;
			}
			
			// if the reference particle is flagged as non valid the event needs to be discarded, as well.
			if((use_particle_as_ref) && (!ref_hit_valid_array[ref]))
				event_valid = false;
			
			// check if we found a valid event with the current reference particle
			if(event_valid) {
				break;
			}
		}

		this->result = event_valid;

		// create presorted event
		if(event_valid) {

			// We need to cleanup in case we used a reference particle and multihit:
			// Reference particle will be particle 0
			for(int i=0;i<ref;i++)
				ref_hit_valid_array[i]=false;


			this->evto.reaction = this->Flag;
			this->evto.scanval = this->evti.scanval;
			this->evto.e.num_hits = 0;
			this->evto.r.num_hits = 0;
			this->evto.p.num_hits = 0;

			this->evto.e.x.clear();
			this->evto.e.y.clear();
			this->evto.e.time.clear();
			this->evto.e.tof.clear();
			this->evto.e.method.clear();
			this->evto.r.x.clear();
			this->evto.r.y.clear();
			this->evto.r.time.clear();
			this->evto.r.tof.clear();
			this->evto.r.method.clear();
			this->evto.p.x.clear();
			this->evto.p.y.clear();
			this->evto.p.time.clear();
			this->evto.p.tof.clear();
			this->evto.p.method.clear();

			// Deal with electron multihit, which might be corrupt when using 
			// modulo BM-spacing...
			if(this->ctof->usesBM() && (this->ctof->GetBMspacing()<1000.0) && (this->ctof->GetBMspacing()>0.5)) {
				double max_spread = this->ctof->GetBMspacing();
				double min_time;
				bool min_valid = false;

				for(int i=0;i<this->evti.e.num_hits;i++) {
					if(this->e_hit_valid_array[i]) {
						if(min_valid) {
							if(fabs(evti.e.time[i] - min_time) > max_spread)
								this->e_hit_valid_array[i] = false;
						} else {
							min_time = evti.e.time[i]; 
							min_valid = true;
						}
					}
				}
			}

			for(int i=0; i<this->evti.e.num_hits; i++) {
				if(this->e_hit_valid_array[i]) {
					this->evto.e.x.push_back(this->evti.e.x[i]);
					this->evto.e.y.push_back(this->evti.e.y[i]);
					this->evto.e.time.push_back(this->evti.e.time[i]);
					this->evto.e.tof.push_back(this->evti.e.tof[i]);
					this->evto.e.method.push_back(this->evti.e.method[i]);
					this->evto.e.num_hits++;
				}
			}

			for(int i=0; i<this->evti.r.num_hits; i++) {
				if(this->r_hit_valid_array[i]) {
					this->evto.r.x.push_back(this->evti.r.x[i]);
					this->evto.r.y.push_back(this->evti.r.y[i]);
					this->evto.r.time.push_back(this->evti.r.time[i]);
					this->evto.r.tof.push_back(this->evti.r.tof[i]);
					this->evto.r.method.push_back(this->evti.r.method[i]);
					this->evto.r.num_hits++;
				}
			}

			for(int i=0; i<this->evti.p.num_hits; i++) {
				if(this->p_hit_valid_array[i]) {
					this->evto.p.x.push_back(this->evti.p.x[i]);
					this->evto.p.y.push_back(this->evti.p.y[i]);
					this->evto.p.time.push_back(this->evti.p.time[i]);
					this->evto.p.tof.push_back(this->evti.p.tof[i]);
					this->evto.p.method.push_back(this->evti.p.method[i]);
					this->evto.p.num_hits++;
				}
			}
		}
		return event_valid;
	}

	bool presorter_class::tofgood(CH_det_struct* particle, int i, double min, double max) { 
		return((particle->tof[i] >= min) && (particle->tof[i] < max));
	}	

	bool presorter_class::posgood(CH_det_struct* particle, int i, double center_x, double center_y, double radius) { 
		return(sqrt((particle->x[i] -center_x)*(particle->x[i] - center_x) + (particle->y[i] - center_y)*(particle->y[i] - center_y)) < radius);
	}
	void presorter_class::add_user_parameter(int prs_num, double val) {
		this->usr_par[prs_num][this->num_of_usr_par[prs_num]] = val;
		this->num_of_usr_par[prs_num]++;
	}

	void presorter_class::set_hits_e(int min, int max) {
		this->ehit_min = min;
		this->ehit_max = max;
	}

	void presorter_class::set_hits_r(int min, int max) {
		this->rhit_min = min;
		this->rhit_max = max;
	}

	void presorter_class::set_hits_p(int min, int max) {
		this->phit_min = min;
		this->phit_max = max;
	}

	void presorter_class::set_pos_e(double pos_center_x, double pos_center_y, double pos_radius) {
		this->epos_center_x  = pos_center_x;
		this->epos_center_y  = pos_center_y;
		this->epos_radius  = pos_radius;
	}

	void presorter_class::set_tof_e(double tof_min, double tof_max) {
		this->etof_min = tof_min;
		this->etof_max  = tof_max;
	}

	void presorter_class::set_pos_r(double pos_center_x, double pos_center_y, double pos_radius) {
		this->rpos_center_x  = pos_center_x;
		this->rpos_center_y  = pos_center_y;
		this->rpos_radius  = pos_radius;
	}

	void presorter_class::set_tof_r(double tof_min, double tof_max) {
		this->rtof_min = tof_min;
		this->rtof_max  = tof_max;
	}

	void presorter_class::set_pos_p(double pos_center_x, double pos_center_y, double pos_radius) {
		this->ppos_center_x  = pos_center_x;
		this->ppos_center_y  = pos_center_y;
		this->ppos_radius  = pos_radius;
	}

	void presorter_class::set_tof_p(double tof_min, double tof_max) {
		this->ptof_min = tof_min;
		this->ptof_max  = tof_max;
	}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
	//////////////////////////   N U M    H I T S   ///////////////////////////////////////

	bool presorter_class::hits_e() {
		int num_hits = 0;
		for(int i=0;i< this->evti.e.num_hits;i++) {
			if(num_hits >= this->ehit_max)
				this->e_hit_valid_array[i] = false;
			if(this->e_hit_valid_array[i])
				num_hits++;
		}
		if(num_hits < this->ehit_min) {
			return false;
		}
		return true;
	}

	bool presorter_class::hits_r() {
		int num_hits = 0;
		for(int i=0;i< this->evti.r.num_hits;i++) {
			if(num_hits >= this->rhit_max)
				this->r_hit_valid_array[i] = false;
			if(this->r_hit_valid_array[i])
				num_hits++;
		}
		if(num_hits < this->rhit_min) {
			return false;
		}
		return true;
	}

	bool presorter_class::hits_p() {
		int num_hits = 0;
		for(int i=0;i< this->evti.p.num_hits;i++) {
			if(num_hits >= this->phit_max)
				this->p_hit_valid_array[i] = false;
			if(this->p_hit_valid_array[i])
				num_hits++;
		}
		if(num_hits < this->phit_min) {
			return false;
		}
		return true;
	}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
	//////////////////////////   T O F   ///////////////////////////////////////

	bool presorter_class::tof_e() {
		int num_hits = 0;
		for(int i=0;i<this->evti.e.num_hits;i++) {
			H->fill1(0 + this->HistoBlock*NHISTS,"electron tof(no prs)",evti.e.tof[i],1.,"electron TOF",200,0.,this->etof_max*1.1,"TOF [ns]",dir);
			if(this->e_hit_valid_array[i]) {
				if(tofgood(&this->evti.e, i, this->etof_min, this->etof_max)) {
					H->fill1(1 + this->HistoBlock*NHISTS,"electron tof",evti.e.tof[i],1.,"electron TOF",200,0.,this->etof_max*1.1,"TOF [ns]",dir);
					num_hits++;			
				} else {
					e_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0) 
			return true;
		else
			return false;
	}

	bool presorter_class::tof_r() {
		int num_hits = 0;
		for(int i=0;i<this->evti.r.num_hits;i++) {
			H->fill1(2 + this->HistoBlock*NHISTS,"ion tof(no prs)",evti.r.tof[i],1.,"ion TOF",2000,0.,rtof_max*1.1,"TOF [ns]",dir);
			if(this->r_hit_valid_array[i]) {
				if(tofgood(&this->evti.r, i, this->rtof_min, this->rtof_max)) {
					H->fill1(3 + this->HistoBlock*NHISTS,"ion tof",evti.r.tof[i],1.,"ion TOF",2000,0.,rtof_max*1.1,"TOF [ns]",dir);
					num_hits++;			
				} else {
					r_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0)
			return true;
		else
			return false;
	}

	bool presorter_class::tof_p() {
		int num_hits = 0;
		for(int i=0;i<this->evti.p.num_hits;i++) {
			H->fill1(4 + this->HistoBlock*NHISTS,"projectile tof(no prs)",evti.p.tof[i],1.,"projectile TOF",200,0.,ptof_max*1.1,"TOF [ns]",dir);
			if(this->p_hit_valid_array[i]) {
				if(tofgood(&this->evti.p, i, this->ptof_min, this->ptof_max)) {
					H->fill1(5 + this->HistoBlock*NHISTS,"projectile tof",evti.p.tof[i],1.,"projectile TOF",200,0.,ptof_max*1.1,"TOF [ns]",dir);
					num_hits++;			
				} else {
					p_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0)
			return true;
		else
			return false;
	}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------
	//////////////////////////   P O S   ///////////////////////////////////////

	bool presorter_class::pos_e() {
		int num_hits = 0;
		for(int i=0;i<this->evti.e.num_hits;i++) {
			H->fill2(6 + this->HistoBlock*NHISTS,"electron pos(no prs)",evti.e.x[i],evti.e.y[i],1.,"electron position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
			if(this->e_hit_valid_array[i]) {
				if(posgood(&this->evti.e, i, this->epos_center_x, this->epos_center_y, this->epos_radius)) {
					H->fill2(7 + this->HistoBlock*NHISTS,"electron pos",evti.e.x[i],evti.e.y[i],1.,"electron position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
					num_hits++;			
				} else {
					e_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0)
			return true;
		else
			return false;
	}

	bool presorter_class::pos_r() {
		int num_hits = 0;
		for(int i=0;i<this->evti.r.num_hits;i++) {
			H->fill2(8 + this->HistoBlock*NHISTS,"ion pos(no prs)",evti.r.x[i],evti.r.y[i],1.,"ion position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
			if(this->r_hit_valid_array[i]) {
				if(posgood(&this->evti.r, i, this->rpos_center_x, this->rpos_center_y, this->rpos_radius)) {
					H->fill2(9 + this->HistoBlock*NHISTS,"ion pos",evti.r.x[i],evti.r.y[i],1.,"ion position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
					num_hits++;			
				} else {
					r_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0)
			return true;
		else
			return false;
	}

	bool presorter_class::pos_p() {
		int num_hits = 0;
		for(int i=0;i<this->evti.p.num_hits;i++) {
			H->fill2(10 + this->HistoBlock*NHISTS,"projectile pos(no prs)",evti.p.x[i],evti.p.y[i],1.,"projectile position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
			if(this->p_hit_valid_array[i]) {
				if(posgood(&this->evti.p, i, this->ppos_center_x, this->ppos_center_y, this->ppos_radius)) {
					H->fill2(11 + this->HistoBlock*NHISTS,"projectile pos",evti.p.x[i],evti.p.y[i],1.,"projectile position",200,-60.0,60.0,"x [mm]",200,-60.0,60.0,"y [mm]",dir);
					num_hits++;			
				} else {
					p_hit_valid_array[i] = false;
				}
			}
		}
		if(num_hits > 0)
			return true;
		else
			return false;
	}


//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	//////////////////////////   P I P I C O   ///////////////////////////////////////
	// using spectrometer struct
	void presorter_class::set_pipico(double mass1_amu, double charge1_au, double mass2_amu, double charge2_au, double pipico_width) 
	{
		this->pipico_width_ns = pipico_width;
		this->mass1_amu = mass1_amu;
		this->charge1_au = charge1_au;
		this->mass2_amu = mass2_amu;
		this->charge2_au = charge2_au;
		this->pipico_width_ns = pipico_width_ns;
		this->acceleration1_mm = (spect->ion_side->number_of_regions>0 ? spect->ion_side->lengths[0] : 0.0);
		this->acceleration2_mm = (spect->ion_side->number_of_regions>1 ? spect->ion_side->lengths[1] : 0.0);
		this->acceleration3_mm = (spect->ion_side->number_of_regions>2 ? spect->ion_side->lengths[2] : 0.0);
		this->efield1_Vcm = (spect->ion_side->number_of_regions>0 ? spect->ion_side->Efields[0] : 0.0);
		this->efield2_Vcm = (spect->ion_side->number_of_regions>1 ? spect->ion_side->Efields[1] : 0.0);
		this->efield3_Vcm = (spect->ion_side->number_of_regions>2 ? spect->ion_side->Efields[2] : 0.0);
	 }

	bool presorter_class::pipico()
	{	
		if(evti.r.num_hits<2)
			return false;

		// assume there is no valid ion pair
		for (int i=0;i<(evti.r.num_hits);i++) {
			r_hit_valid_array[i] = false;
		}

		H->fill2(15 + this->HistoBlock*NHISTS,"pipico big range(no prs, first two hits)",evti.r.tof[0],evti.r.tof[1],1.,"PIPICO spectrum (coarse)",500,0.,50000.,"ion1 TOF [ns]",500,0.,50000.,"ion2 TOF [ns]",dir);
		H->fill2(16 + this->HistoBlock*NHISTS,"pipico small range(no prs, first two hits)",evti.r.tof[0],evti.r.tof[1],1.,"PIPICO spectrum (fine)",800,80.,20000.,"ion1 TOF [ns]",800,80.,20000.,"ion2 TOF [ns]",dir);
		
		for (int j=0;j<(evti.r.num_hits);j++) {
			for (int k=j+1;k<(evti.r.num_hits);k++) {

				H->fill2(17 + this->HistoBlock*NHISTS,"pipico big range(no prs, all hits)",evti.r.tof[j],evti.r.tof[k],1.,"PIPICO spectrum (coarse)",500,0.,50000.,"ion1 TOF [ns]",500,0.,50000.,"ion2 TOF [ns]",dir);
				H->fill2(18 + this->HistoBlock*NHISTS,"pipico small range(no prs, all hits)",evti.r.tof[j],evti.r.tof[k],1.,"PIPICO spectrum (fine)",800,80.,20000.,"ion1 TOF [ns]",800,80.,20000.,"ion2 TOF [ns]",dir);

				double t2 = t2_3accel(evti.r.tof[j], this->acceleration1_mm, this->acceleration2_mm, this->acceleration3_mm, this->efield1_Vcm, this->efield2_Vcm, this->efield3_Vcm, this->charge1_au, this->charge2_au, this->mass1_amu, this->mass2_amu);
				double t2swap = t2_3accel(evti.r.tof[j], this->acceleration1_mm, this->acceleration2_mm, this->acceleration3_mm, this->efield1_Vcm, this->efield2_Vcm, this->efield3_Vcm, this->charge2_au, this->charge1_au, this->mass2_amu, this->mass1_amu);

				if(fabs(evti.r.tof[k] - t2) < this->pipico_width_ns || fabs(evti.r.tof[k] - t2swap) < this->pipico_width_ns) {
					if(evti.r.tof[j]>0.0 && evti.r.tof[k]>0.0) {
						// found a pair, set the corresponding ions to 'true'
						r_hit_valid_array[j] = true;
						r_hit_valid_array[k] = true;
						H->fill2(19 + this->HistoBlock*NHISTS,"pipico big range",evti.r.tof[j],evti.r.tof[k],1.,"PIPICO spectrum (coarse)",500,0.,50000.,"ion1 TOF [ns]",500,0.,50000.,"ion2 TOF [ns]",dir);
						H->fill2(20 + this->HistoBlock*NHISTS,"pipico small range",evti.r.tof[j],evti.r.tof[k],1.,"PIPICO spectrum (fine)",800,80.,20000.,"ion1 TOF [ns]",800,80.,20000.,"ion2 TOF [ns]",dir);
						return true;
					}
				}
			}
		}

		// we looped over all ions, but did not find any pair that fits...
		return false;
	}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	//////////////////////////   S U M - D I F F   T O F s   ///////////////////////////////////////
	// electron & recoil (sum/diff) TOF for n particles
	void presorter_class::set_sumdiff(int number_of_recoils_to_take, double rtof_sum_center, double rtof_sum_width, double rtof_diff_center, double rtof_diff_width)
	{
		this->number_of_recoils_to_take = number_of_recoils_to_take; 
		this->rtof_sum_center = rtof_sum_center; 
		this->rtof_sum_width  = rtof_sum_width;  
		this->rtof_diff_center = rtof_diff_center; 
		this->rtof_diff_width  = rtof_diff_width;
	}

	bool presorter_class::sumdiff()
	{ 
		if((evti.r.num_hits < number_of_recoils_to_take) )
			return false;

		// assume there is no valid ion combination
		for (int i=0;i<(evti.r.num_hits);i++) {
			r_hit_valid_array[i] = false;
		}

		double rtof_sum, rtof_diff, rtof_sum_a, rtof_sum_b;
		double rtof_sum_found, rtof_diff_found, rtof_sum_a_found, rtof_sum_b_found;

		rtof_sum=0.;
		rtof_sum_a=0.;  // first ions
		rtof_sum_b=0.;  // second ions
		rtof_diff=0.;

		rtof_sum_found=0.;
		rtof_sum_a_found=0.;  // first ions
		rtof_sum_b_found=0.;  // second ions
		rtof_diff_found=0.;

		for (int j=0;j<(evti.r.num_hits-number_of_recoils_to_take+0.5);j++) //xxxm
		{ // begin for: j - first "fitting" ion
			int i1 = 0;
			rtof_sum=0.;
			rtof_sum_a=0.;  // first ions
			rtof_sum_b=0.;  // second ions
			rtof_diff=0.;

			rtof_sum_found=0.;
			rtof_sum_a_found=0.;  // first ions
			rtof_sum_b_found=0.;  // second ions
			rtof_diff_found=0.;

			for (int l=0;l<number_of_recoils_to_take;l++)
			{ // begin for: recoils j...j+l
				//rtof_sum = rtof_sum + evti->r.tof[j+l];
				if (l>((number_of_recoils_to_take-0.5)/2.)) {rtof_sum_b = rtof_sum_b + evti.r.tof[j+l];}  // (3+4)
				else { rtof_sum_a = rtof_sum_a + evti.r.tof[j+l];i1++;}    // (0+1+2)
				rtof_sum = rtof_sum_a + rtof_sum_b;
				rtof_diff = rtof_sum_b - rtof_sum_a;
			} // end for: recoils j...j+l
			char title1[50];
			char xtitle1[50];
			char ytitle1[50];
			char title2[50];
			char xtitle2[50];
			char ytitle2[50];
			sprintf(title1,"%s%li%s","PI",number_of_recoils_to_take,"CO (coarse)");
			sprintf(xtitle1,"%s%li%s","Sum 1...",i1," ions");
			sprintf(ytitle1,"%s%li%s%li%s","Sum ",i1+1,"...",number_of_recoils_to_take," ions");
			sprintf(title2,"%s%li%s","Sum TOF vs. diff TOF (",number_of_recoils_to_take," ions)");
			sprintf(xtitle2,"%s%li%s","Sum TOF ",number_of_recoils_to_take," ions");
			sprintf(ytitle2,"%s%li%s","Diff TOF ",number_of_recoils_to_take," ions");
				
			H->fill2(15 + this->HistoBlock*NHISTS,"PInCO coarse (no prs, just first hits)",rtof_sum_a,rtof_sum_b,1.,title1,1000,0.,20000.,xtitle1,1000,0.,20000.,ytitle1,dir);
			H->fill2(16 + this->HistoBlock*NHISTS,"sum/diff coarse (no prs, just first hits)",rtof_sum,rtof_diff,1.,title2,1500,0.,30000.,xtitle2,750,-5000.,15000.,ytitle2,dir);
			H->fill2(17 + this->HistoBlock*NHISTS,"sum/diff fine (no prs, just first hits)",rtof_sum,rtof_diff,1.,title2,500,(rtof_sum_center-3.*rtof_sum_width),(rtof_sum_center+3.*rtof_sum_width),xtitle2,500,(rtof_diff_center-3.*rtof_diff_width),(rtof_diff_center+3.*rtof_diff_width),ytitle2,dir);
			if( (fabs(rtof_sum-rtof_sum_center)<rtof_sum_width) && (fabs(rtof_diff-rtof_diff_center)<rtof_diff_width) )
			{ // begin if: good recoil sum/diff tof (j)
				for (int k=0;k<number_of_recoils_to_take;k++)
				{ // begin for: store recoil hits (k)
					r_hit_valid_array[j+k] = true;
				} // end for: store recoil hits (k)

				H->fill2(10504+this->Flag*1000,"sum/diff fine (prs, initial found hits)",rtof_sum,rtof_diff,1.,title2,500,0.,20000.,"rtofsum",500,-7000.,7000.,"rtof_diff",dir);
				// calculate sum, sum_a, sum_b and diff for found ions!
				int ion = 0;
				for (int k=0;k<number_of_recoils_to_take;k++)
				{ // begin for: recoils 0...k
					if(r_hit_valid_array[k]) {
						if (ion>((number_of_recoils_to_take-0.5)/2.)) {
							rtof_sum_b_found = rtof_sum_b_found + evti.r.tof[k];
						}  // (3+4)
						else {
							rtof_sum_a_found = rtof_sum_a_found + evti.r.tof[k];
						}    // (0+1+2)
						rtof_sum_found = rtof_sum_b_found + rtof_sum_a_found;
						rtof_diff_found = rtof_sum_b_found - rtof_sum_a_found;
						ion++;
					}
				} // end for: recoils 0...k
				H->fill2(18 + this->HistoBlock*NHISTS,"sum/diff big range coarse (good recoils)",rtof_sum_found,rtof_diff_found,1.,title1,250,0.,30000.,xtitle1,250,-5000.,15000.,ytitle1,dir);
				H->fill2(19 + this->HistoBlock*NHISTS,"sum/diff big range fine (good recoil)",rtof_sum_found,rtof_diff_found,1.,title2,500,(rtof_sum_center-1.5*rtof_sum_width),(rtof_sum_center+1.5*rtof_sum_width),xtitle2,500,(rtof_diff_center-1.5*rtof_diff_width),(rtof_diff_center+1.5*rtof_diff_width),ytitle2,dir);
				return true;
			} // end if: good recoil sum/diff tof (with loop over j)
		} // end for: recoil hits (j)
		return false;
	} 

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------------------------------------

	//////////////////////////   P O L Y A T O M I C   ///////////////////////////////////////
	void presorter_class::set_polyatomic(int min_number_of_recoils)
	{
		this->min_number_of_recoils = min_number_of_recoils;	
		this->ions_defined=0;
	}

	void presorter_class::polyatomic_add_ion(double tmin, double tmax)
	{
		this->tofmin[this->ions_defined]=tmin;
		this->tofmax[this->ions_defined++]=tmax;
	}

	bool presorter_class::polyatomic()
	{	
		if(evti.r.num_hits < min_number_of_recoils )
			return false;
		
		// assume there is no valid ion combination
		for (int i=0;i<(evti.r.num_hits);i++) {
			r_hit_valid_array[i] = false;
		}

		char title[50];
		char xtitle[50];
		char ytitle[50];
		
		// plot Martin's histograms!
		double rtof_sum1 = 0.;
		double rtof_sum2 = 0.;	
		int i1, i2;
		for (int i = 0; i < (min_number_of_recoils-0.01)/2; i++)
           			{rtof_sum1+= evti.r.tof[i]; i1=i;}
		i2 = int((min_number_of_recoils-0.01)/2) + 1;
		for (int i = i2; i < min_number_of_recoils; i++)
         			rtof_sum2+= evti.r.tof[i];

		sprintf_s(title,50,"%s%li%s","PI",min_number_of_recoils,"CO(no prs)");
		sprintf_s(xtitle,50,"%s%li%s","Sum 1...",i1+1," ions");
		sprintf_s(ytitle,50,"%s%li%s%li%s","Sum ",i2+1,"...",min_number_of_recoils," ions");

		H->fill2(15 + this->HistoBlock*NHISTS,title,rtof_sum1,rtof_sum2,1.,title,500,0.,10000.,xtitle,500,0.,15000.,ytitle,dir);
						
		if(ions_defined>0) {
		// We have some tof-min/max definitions, so try to find these ions!
			int current_ion = 0;

			for(int i=0; i<evti.r.num_hits; i++) {
				if(evti.r.tof[i] > this->tofmin[current_ion] && evti.r.tof[i] < this->tofmax[current_ion]) {
					r_hit_valid_array[i] = true;
					current_ion++;
				}
				if(current_ion > ions_defined)
					break;
			}

			if(current_ion>=min_number_of_recoils) {
				// plot histograms
				int ion=0;
				rtof_sum1 = 0.;
				rtof_sum2 = 0.;
				for(int i=0; i<evti.r.num_hits; i++) {
					if(r_hit_valid_array[i] == true) {
						sprintf_s(title,50,"Ion %i fish X",ion + 1);
						H->fill2(17 + this->HistoBlock*NHISTS+(ion*2),title,evti.r.tof[i],evti.r.x[i],1.,title,int((this->tofmax[ion]-this->tofmin[ion])/2.0),this->tofmin[ion]*0.9,this->tofmax[ion]*1.1,"TOF [ns]",200,-1.*60.0,60.0,"ion x-pos [mm]",dir);
						sprintf_s(title,50,"Ion %i fish Y",ion + 1);
						H->fill2(17 + this->HistoBlock*NHISTS+(ion*2)+1,title,evti.r.tof[i],evti.r.y[i],1.,title,int((this->tofmax[ion]-this->tofmin[ion])/2.0),this->tofmin[ion]*0.9,this->tofmax[ion]*1.1,"TOF [ns]",200,-1.*60.0,60.0,"ion y-pos [mm]",dir);
						
						if(ion<(min_number_of_recoils-0.01)/2)
							rtof_sum1 += evti.r.tof[i];
						else 
							rtof_sum2 += evti.r.tof[i];

						ion++;
					}
				}
				sprintf_s(title,50,"%s%li%s","PI",min_number_of_recoils,"CO");
				H->fill2(16 + this->HistoBlock*NHISTS,title,rtof_sum1,rtof_sum2,1.,title,500,0.,10000.,xtitle,500,0.,15000.,ytitle,dir);

				return true;
			} else 
				return false;
		}
		// There were no ion definitions but event has enough ions, so flag  
		// first ions as ok and return "true"...
		for (int i=0;i<min_number_of_recoils;i++) {
			r_hit_valid_array[i] = true;
		}
		return true;
	}
}


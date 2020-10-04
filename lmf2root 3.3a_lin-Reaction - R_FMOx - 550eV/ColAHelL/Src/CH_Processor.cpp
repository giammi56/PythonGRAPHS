#include "CH_Processor.h"

#include "Math.h"


#include "assert.h"

namespace CH
{
	processor_class::processor_class() 
	{
		this->Spect = 0;

		cut_wiggle = false;
		cut_wiggle_width = 0.0;

		e_fac = new cor_param();
		e_fac->dt = 0.0;
		e_fac->dx = 0.0;
		e_fac->dy = 0.0;
		e_fac->overall_stretch = 1.0;
		e_fac->x_stretch = 1.0;
		e_fac->y_stretch = 1.0;
		e_fac->mir_x = false;
		e_fac->mir_y = false;
		e_fac->EBx = 0.0;
		e_fac->EBy = 0.0;
		e_fac->rot_ang = 0.0;
		e_fac->raw_order = true;

		e_mom_fac = new mom_cor_param();
		e_mom_fac->dz = 0.0;
		e_mom_fac->dx = 0.0;
		e_mom_fac->dy = 0.0;
		e_mom_fac->overall_stretch = 1.0;
		e_mom_fac->x_stretch = 1.0;
		e_mom_fac->y_stretch = 1.0;
		e_mom_fac->z_stretch = 1.0;
		e_mom_fac->mir_x = false;
		e_mom_fac->mir_y = false;
		e_mom_fac->mir_z = false;
		e_mom_fac->rot_ang = 0.0;

		r_fac = new cor_param();
		r_fac->dt = 0.0;
		r_fac->dx = 0.0;
		r_fac->dy = 0.0;
		r_fac->overall_stretch = 1.0;
		r_fac->x_stretch = 1.0;
		r_fac->y_stretch = 1.0;
		r_fac->mir_x = false;
		r_fac->mir_y = false;
		r_fac->EBx = 0.0;
		r_fac->EBy = 0.0;
		r_fac->rot_ang = 0.0;
		r_fac->raw_order = true;

		r_mom_fac = new mom_cor_param();
		r_mom_fac->dz = 0.0;
		r_mom_fac->dx = 0.0;
		r_mom_fac->dy = 0.0;
		r_mom_fac->overall_stretch = 1.0;
		r_mom_fac->x_stretch = 1.0;
		r_mom_fac->y_stretch = 1.0;
		r_mom_fac->z_stretch = 1.0;
		r_mom_fac->mir_x = false;
		r_mom_fac->mir_y = false;
		r_mom_fac->mir_z = false;
		r_mom_fac->rot_ang = 0.0;

		p_fac = new cor_param();
		p_fac->dt = 0.0;
		p_fac->dx = 0.0;
		p_fac->dy = 0.0;
		p_fac->overall_stretch = 1.0;
		p_fac->x_stretch = 1.0;
		p_fac->y_stretch = 1.0;
		p_fac->mir_x = false;
		p_fac->mir_y = false;
		p_fac->EBx = 0.0;
		p_fac->EBy = 0.0;
		p_fac->rot_ang = 0.0;
		r_fac->raw_order = true;

		p_mom_fac = new mom_cor_param();
		p_mom_fac->dz = 0.0;
		p_mom_fac->dx = 0.0;
		p_mom_fac->dy = 0.0;
		p_mom_fac->overall_stretch = 1.0;
		p_mom_fac->x_stretch = 1.0;
		p_mom_fac->y_stretch = 1.0;
		p_mom_fac->z_stretch = 1.0;
		p_mom_fac->mir_x = false;
		p_mom_fac->mir_y = false;
		p_mom_fac->mir_z = false;
		p_mom_fac->rot_ang = 0.0;

		for (__int32 i=0;i<MAX_NUM_REACTIONS;i++) reaction_list[i] = 0;

		number_of_reactions_defined = 0;
		cur_reac_num = 0;
	}
	
	processor_class::~processor_class() 
	{

		if (e_fac) delete e_fac;
		if (e_mom_fac) delete e_mom_fac;
		
		if (r_fac) delete r_fac;
		if (r_mom_fac) delete r_mom_fac;
		
		if (p_fac) delete p_fac;
		if (p_mom_fac) delete p_mom_fac;

		for (__int32 i=0;i<number_of_reactions_defined;i++) {
			assert(reaction_list[i]);
			if (reaction_list[i]) {
				if (reaction_list[i]->e_fac) {delete reaction_list[i]->e_fac; reaction_list[i]->e_fac = 0;}
				if (reaction_list[i]->r_fac) {delete reaction_list[i]->r_fac; reaction_list[i]->r_fac = 0;}
				if (reaction_list[i]->p_fac) {delete reaction_list[i]->p_fac; reaction_list[i]->p_fac = 0;}
				if (reaction_list[i]->e_mom_fac) {delete reaction_list[i]->e_mom_fac; reaction_list[i]->e_mom_fac = 0;}
				if (reaction_list[i]->r_mom_fac) {delete reaction_list[i]->r_mom_fac; reaction_list[i]->r_mom_fac = 0;}
				if (reaction_list[i]->p_mom_fac) {delete reaction_list[i]->p_mom_fac; reaction_list[i]->p_mom_fac = 0;}
				if (reaction_list[i]->MF_cond) {delete reaction_list[i]->MF_cond; reaction_list[i]->MF_cond = 0;}
				if (reaction_list[i]->LF_cond) {delete reaction_list[i]->LF_cond; reaction_list[i]->LF_cond = 0;}

				for (__int32 k=0;k<16;k++) reaction_list[i]->interpolation_filename->clear();
				reaction_list[i]->interpolation_points.clear();
				reaction_list[i]->interpolation_adjustment_variables.clear();

				for (__int32 k=0;k<reaction_list[i]->num_tags;k++) {
					if (reaction_list[i]->tag[k]) {
						delete reaction_list[i]->tag[k];
						reaction_list[i]->tag[k] = 0;
					}
				}

				delete reaction_list[i];
				reaction_list[i] = 0;
			}
		}
		number_of_reactions_defined = 0;
	}

	// Will give you only the first one, even if several reactions share same channel
	int processor_class::find_reaction(int channel) {
		return find_reaction(channel, 0);
	}

	// Will give you the (n+1)th reaction if several reactions share same channel
	int processor_class::find_reaction(int channel, int skip_n_reacs) {
		int skipped = 0;
		for(int i=0;i<this->number_of_reactions_defined;i++) {
			if(this->reaction_list[i]->channel == channel) {
				if(skipped == skip_n_reacs)
					return i;
				else
					skipped++;
			}
		}
		return -1;
	}

	// reports how many reactions are defined for channel CHAN
	int processor_class::EvtBelongsToReactions(int chan) {
		int num_reacs = 0;
		for(int i=0;i<this->number_of_reactions_defined;i++) {
			if(this->reaction_list[i]->channel == chan) {
				num_reacs++;
			}
		}
		return num_reacs;
	}
	// return the reaction used for previous event processing
	reaction_struct * processor_class::get_current_reaction() {
		return this->reaction_list[cur_reac_num];
	}

	bool processor_class::process_all(CH_event_struct *evti, electron_class ** el, ion_class ** rec, diatomic_class * molec, polyatomic_class * big_molec, int which_reaction_def) {
		cur_reac_num = find_reaction((int)evti->reaction, which_reaction_def);
		if(cur_reac_num>-1) { 
			// process ions (diatomic, polyatomic etc.) 
			bool valid_ions = process_ions(this->reaction_list[cur_reac_num], evti, rec, molec, big_molec);
			// Process electrons
			// only even reaction types include electrons!
			bool valid_electrons = false;
			if(this->reaction_list[cur_reac_num]->type % 2 == 0)
				valid_electrons = process_electrons(this->reaction_list[cur_reac_num], evti, el);	
			return valid_ions && valid_electrons;
		}	
		return false;
	}
				
	bool processor_class::process_electrons(reaction_struct *cur_reaction, CH_event_struct *evti, electron_class ** e) 
	{
		
		int ehit = evti->e.num_hits;
		this->at_least_one_e_valid = false;

		for(int i=0;i<(int)ehit;i++) {
				
			e_raw[i].method = evti->e.method[i];
			e_raw[i].data.x = evti->e.x[i];
			e_raw[i].data.y = evti->e.y[i];
			e_raw[i].data.time = evti->e.time[i];
			e_raw[i].data.tof = evti->e.tof[i];
			e_raw[i].m = 1.0 / MASSAU;
			e_raw[i].q = -1.0;

			e[i]->set_raw(e_raw[i]);
			e[i]->valid = true;

			if(this->e_fac->raw_order) {
// GN				e[i]->shift_stretch_raw(this->e_fac);
				e[i]->shift_stretch_raw(cur_reaction->e_fac);
				e[i]->rotate_raw(cur_reaction->e_fac->rot_ang);
			} else {
				e[i]->rotate_raw(cur_reaction->e_fac->rot_ang);
//	GN_end			e[i]->shift_stretch_raw(this->e_fac);
				e[i]->shift_stretch_raw(cur_reaction->e_fac);
			}
			e[i]->calc_phi_pos();

			//cannot process electrons any further without spectrometer def.
			if(this->Spect==0)
				return true;
			e[i]->process(this->Spect);
			
			//e[i]->shift_stretch_mom(this->e_mom_fac);
			e[i]->shift_stretch_mom(cur_reaction->e_mom_fac);

			if(this->cut_wiggle)
				e[i]->cut_wiggle(this->Spect, this->cut_wiggle_width);	
	
//			e[i]->autotune();

			// set electrons as not valid as requested by user
			for(int k=0;k<cur_reaction->num_tags;k++) {
				if(cur_reaction->tag[k]->type == RTYPE_ELEC) {
					e[i]->invalidate(cur_reaction->tag[k]);
				}
			}
			if(e[i]->valid)
				this->at_least_one_e_valid = true;
		}

		// has no meaning so far.. 
		return this->at_least_one_e_valid;
	}

	void processor_class::new_rtag(__int32 edit_reaction_num)
	{
		reaction_list[edit_reaction_num]->tag[reaction_list[edit_reaction_num]->num_tags++] = new rtag_struct();
	}

	void processor_class::new_reaction()
	{
		reaction_list[number_of_reactions_defined] = new reaction_struct();
		reaction_list[number_of_reactions_defined]->reac_num = number_of_reactions_defined;
	
		memset(reaction_list[number_of_reactions_defined]->tag,0,16*sizeof(rtag_struct *));

		reaction_list[number_of_reactions_defined]->two_elec = false;
		reaction_list[number_of_reactions_defined]->e_master_min = -1.0;
		reaction_list[number_of_reactions_defined]->e_master_max = 10000000.0;

		reaction_list[number_of_reactions_defined]->photon_scan = false;
		reaction_list[number_of_reactions_defined]->ph_scan_channel = -1;
		reaction_list[number_of_reactions_defined]->ph_min = -1.0;
		reaction_list[number_of_reactions_defined]->ph_max = -1.0;
		reaction_list[number_of_reactions_defined]->ph_step = -1.0;

		reaction_list[number_of_reactions_defined]->randomize_ions = false;
		reaction_list[number_of_reactions_defined]->rnd_count = 0;

		reaction_list[number_of_reactions_defined]->Dalitz_array[0] = 0;
		reaction_list[number_of_reactions_defined]->Dalitz_array[1] = 1;
		reaction_list[number_of_reactions_defined]->Dalitz_array[2] = 2;

		reaction_list[number_of_reactions_defined]->Newton_array[0] = 0;
		reaction_list[number_of_reactions_defined]->Newton_array[1] = 1;
		reaction_list[number_of_reactions_defined]->Newton_array[2] = 2;

		reaction_list[number_of_reactions_defined]->num_tags = 0;
		
		reaction_list[number_of_reactions_defined]->MF_cond = new MFPAD_cond_struct();
		reaction_list[number_of_reactions_defined]->MF_cond->type = -1;
		reaction_list[number_of_reactions_defined]->MF_cond->emin = 0.0;
		reaction_list[number_of_reactions_defined]->MF_cond->emax = 10000.0;
		reaction_list[number_of_reactions_defined]->MF_cond->rmin = 0.0;
		reaction_list[number_of_reactions_defined]->MF_cond->rmax = 10000.0;

		reaction_list[number_of_reactions_defined]->LF_cond = new LFPAD_cond_struct();
		reaction_list[number_of_reactions_defined]->LF_cond->type = -1;
		reaction_list[number_of_reactions_defined]->LF_cond->emin = 0.0;
		reaction_list[number_of_reactions_defined]->LF_cond->emax = 10000.0;

		reaction_list[number_of_reactions_defined]->e_fac = new cor_param();
		reaction_list[number_of_reactions_defined]->e_fac->dt = this->e_fac->dt; //0.0;
		reaction_list[number_of_reactions_defined]->e_fac->dx = this->e_fac->dx; //0.0;
		reaction_list[number_of_reactions_defined]->e_fac->dy = this->e_fac->dy; //0.0;
		reaction_list[number_of_reactions_defined]->e_fac->overall_stretch = this->e_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_fac->x_stretch = this->e_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_fac->y_stretch = this->e_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_fac->mir_x = this->e_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->e_fac->mir_y = this->e_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->e_fac->raw_order = this->e_fac->raw_order; //false;

		reaction_list[number_of_reactions_defined]->e_mom_fac = new mom_cor_param();
		reaction_list[number_of_reactions_defined]->e_mom_fac->dx = this->e_mom_fac->dx ; //0.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->dy = this->e_mom_fac->dy ; //0.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->dz = this->e_mom_fac->dz ; //0.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->overall_stretch = this->e_mom_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->x_stretch = this->e_mom_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->y_stretch = this->e_mom_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->z_stretch = this->e_mom_fac->z_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->e_mom_fac->mir_x = this->e_mom_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->e_mom_fac->mir_y = this->e_mom_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->e_mom_fac->mir_z = this->e_mom_fac->mir_z; //false;

		reaction_list[number_of_reactions_defined]->e_fac->EBx = this->e_fac->EBx; //0.0;
		reaction_list[number_of_reactions_defined]->e_fac->EBy = this->e_fac->EBy; //0.0;
		reaction_list[number_of_reactions_defined]->e_fac->rot_ang = this->e_fac->rot_ang; //0.0;

		reaction_list[number_of_reactions_defined]->r_fac = new cor_param();
		reaction_list[number_of_reactions_defined]->r_fac->dt = this->r_fac->dt; //0.0;
		reaction_list[number_of_reactions_defined]->r_fac->dx = this->r_fac->dx; //0.0;
		reaction_list[number_of_reactions_defined]->r_fac->dy = this->r_fac->dy; //0.0;
		reaction_list[number_of_reactions_defined]->r_fac->overall_stretch = this->r_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_fac->x_stretch = this->r_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_fac->y_stretch = this->r_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_fac->mir_x = this->r_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->r_fac->mir_y = this->r_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->r_fac->raw_order = this->r_fac->raw_order; //false;

		reaction_list[number_of_reactions_defined]->r_mom_fac = new mom_cor_param();
		reaction_list[number_of_reactions_defined]->r_mom_fac->dx = this->r_mom_fac->dx ; //0.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->dy = this->r_mom_fac->dy ; //0.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->dz = this->r_mom_fac->dz ; //0.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->overall_stretch = this->r_mom_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->x_stretch = this->r_mom_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->y_stretch = this->r_mom_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->z_stretch = this->r_mom_fac->z_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->r_mom_fac->mir_x = this->r_mom_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->r_mom_fac->mir_y = this->r_mom_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->r_mom_fac->mir_z = this->r_mom_fac->mir_z; //false;

		reaction_list[number_of_reactions_defined]->r_fac->EBx = this->r_fac->EBx; //0.0;
		reaction_list[number_of_reactions_defined]->r_fac->EBy = this->r_fac->EBy; //0.0;
		reaction_list[number_of_reactions_defined]->r_fac->rot_ang = this->r_fac->rot_ang; //0.0;

		reaction_list[number_of_reactions_defined]->p_fac = new cor_param();
		reaction_list[number_of_reactions_defined]->p_fac->dt = this->p_fac->dt; //0.0;
		reaction_list[number_of_reactions_defined]->p_fac->dx = this->p_fac->dx; //0.0;
		reaction_list[number_of_reactions_defined]->p_fac->dy = this->p_fac->dy; //0.0;
		reaction_list[number_of_reactions_defined]->p_fac->overall_stretch = this->p_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_fac->x_stretch = this->p_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_fac->y_stretch = this->p_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_fac->mir_x = this->p_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->p_fac->mir_y = this->p_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->p_fac->raw_order = this->p_fac->raw_order; //false;

		reaction_list[number_of_reactions_defined]->p_mom_fac = new mom_cor_param();
		reaction_list[number_of_reactions_defined]->p_mom_fac->dx = this->p_mom_fac->dx ; //0.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->dy = this->p_mom_fac->dy ; //0.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->dz = this->p_mom_fac->dz ; //0.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->overall_stretch = this->p_mom_fac->overall_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->x_stretch = this->p_mom_fac->x_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->y_stretch = this->p_mom_fac->y_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->z_stretch = this->p_mom_fac->z_stretch; //1.0;
		reaction_list[number_of_reactions_defined]->p_mom_fac->mir_x = this->p_mom_fac->mir_x; //false;
		reaction_list[number_of_reactions_defined]->p_mom_fac->mir_y = this->p_mom_fac->mir_y; //false;
		reaction_list[number_of_reactions_defined]->p_mom_fac->mir_z = this->p_mom_fac->mir_z; //false;

		reaction_list[number_of_reactions_defined]->p_fac->EBx = this->p_fac->EBx; //0.0;
		reaction_list[number_of_reactions_defined]->p_fac->EBy = this->p_fac->EBy; //0.0;
		reaction_list[number_of_reactions_defined]->p_fac->rot_ang = this->p_fac->rot_ang; //0.0;

		number_of_reactions_defined++;
	}





	bool processor_class::process_ions(reaction_struct *cur_reaction, CH_event_struct *evti, ion_class ** r, diatomic_class * mol, polyatomic_class * big_mol)
	{
		// less typing...
		int rhit = evti->r.num_hits;
		
		this->at_least_one_r_valid = false;

				
		if(rhit>0) {
			for(int i=0;i<evti->r.num_hits;i++) {
				r_raw[i].method = evti->r.method[i];
				r_raw[i].data.x = evti->r.x[i];
				r_raw[i].data.y = evti->r.y[i];
				r_raw[i].data.time = evti->r.time[i];
				r_raw[i].data.tof = evti->r.tof[i];
			}/*	
				// we do not know mass and charge, set it to proton to avoid possible crashes...
				r_raw[i].m = 1;
				r_raw[i].q = 1;

				r[i]->set_raw(&r_raw[i]);
				r[i]->set_channel(0,false);
					
				r[i]->set_t_mean(cur_reaction->t_mean[0]);

				if(this->r_fac->raw_order) {
					r[i]->shift_stretch_raw(this->r_fac);		
					r[i]->rotate_raw(this->r_fac->rot_ang);
				} else {
					r[i]->rotate_raw(this->r_fac->rot_ang);
					r[i]->shift_stretch_raw(this->r_fac);		
				}
				r[i]->calc_phi_pos();				
			}	*/			
		}

		//cannot process ions any further without spectrometer def.
		if(this->Spect==0)
			return true;

		switch (cur_reaction->type) {
			case(RTYPE_ELEC):
				break;
			case(RTYPE_ION):
			case(RTYPE_ION_ELEC):
			case(RTYPE_DIATOMIC_DISS):
			case(RTYPE_DIATOMIC_DISS_ELEC):

				for(int i=0;i<(int)rhit;i++) {

					r_raw[i].m = cur_reaction->mass[0];
					r_raw[i].q = cur_reaction->charge[0];

					r[i]->set_raw(&r_raw[i]);
					r[i]->set_channel(cur_reaction->channel,true);
					r[i]->set_t_mean(cur_reaction->t_mean[0]);

					if(this->r_fac->raw_order) {
						//r[i]->shift_stretch_raw(this->r_fac);
						r[i]->shift_stretch_raw(cur_reaction->r_fac);
						r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);
					} else {						
						r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);
						//r[i]->shift_stretch_raw(this->r_fac);
						r[i]->shift_stretch_raw(cur_reaction->r_fac);
					}
					r[i]->calc_phi_pos();

				//	if (this->Spect->ion_side->interpolation_approximation)
				//		r[i]->process(this->Spect, cur_reaction->interpolation_points[0], cur_reaction->interpolation_adjustment_variables[0]);
				//	else
						r[i]->process(this->Spect);

					//r[i]->shift_stretch_mom(this->r_mom_fac);
					r[i]->shift_stretch_mom(cur_reaction->r_mom_fac);

					// set ions as not valid as requested by user
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_ION) {
							if((cur_reaction->tag[k]->particle_num<0) || ( i == cur_reaction->tag[k]->particle_num)) {
								r[i]->invalidate(cur_reaction->tag[k]);
							}
						}
					}
					if(r[i]->valid)
						at_least_one_r_valid = true;
				}

				// if reaction is just an atomic ion : stop here
				if(cur_reaction->type!=RTYPE_DIATOMIC_DISS && cur_reaction->type!=RTYPE_DIATOMIC_DISS_ELEC)
					break;
				
				// Dissociation of diatomic into charged + neutral
				// Create a fake molecule from single ion already processed before...
					
				r[1]->mom.x = -r[0]->mom.x;
				r[1]->mom.y = -r[0]->mom.y;
				r[1]->mom.z = -r[0]->mom.z;

				r[1]->raw.m = cur_reaction->mass[1];
				r[1]->raw.q = cur_reaction->charge[1];
				r[1]->valid = r[0]->valid;

				mol->set_ions(r[0],r[1]);
				mol->set_channel(cur_reaction->channel,true);

				// Get rid of possible overlapping doubly charged channel...
				if(rhit>1) {
					mol->valid = false;
					break;
				}
				//calc rel_mom and cm_mom
				mol->process_diatomic_only(this->Spect);
						
				// set molecules as not valid as requested by user
				for(int k=0;k<cur_reaction->num_tags;k++) {
					if(cur_reaction->tag[k]->type == RTYPE_DIATOMIC) {
						mol->invalidate(cur_reaction->tag[k]);
					}
				}
				// set ions as not valid as requested by user
				for(int k=0;k<cur_reaction->num_tags;k++) {
					if(cur_reaction->tag[k]->type == RTYPE_ION) {
						for(int m=0;m<2;m++) {
							if((cur_reaction->tag[k]->particle_num<0) || ( m == cur_reaction->tag[k]->particle_num)) {
								mol->ion[m]->invalidate(cur_reaction->tag[k]);
							}
						}

						if(!mol->ion[0]->valid || !mol->ion[1]->valid) {
							mol->set_valid(false);									
						}
					}
				}

				break;

			case(RTYPE_DIATOMIC):
			case(RTYPE_DIATOMIC_ELEC):

				if(rhit > 1) {

					mol->set_mom_rel_calc(cur_reaction->diatomic_prel);

					for(int i=0;i<2;i++) {
						r_raw[i].m = cur_reaction->mass[i];
						r_raw[i].q = cur_reaction->charge[i];
						r[i]->set_raw(&r_raw[i]);
						r[i]->set_channel(cur_reaction->channel,true);
						r[i]->set_t_mean(cur_reaction->t_mean[i]);
						if(this->r_fac->raw_order) {							
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);
							r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);				
						} else {
							r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);				
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);
						}

						r[i]->calc_phi_pos();

					}
					mol->set_ions(r[0],r[1]);
					mol->set_channel(cur_reaction->channel,true);

					//mol->set_mom_fac(this->r_mom_fac);
					mol->set_mom_fac(cur_reaction->r_mom_fac);

/*					if (this->Spect->ion_side->interpolation_approximation) {
						if ((cur_reaction->mass[0] != cur_reaction->mass[1]) || (cur_reaction->charge[0] != cur_reaction->charge[1]))
							mol->sort_ions_and_process(Spect, cur_reaction->interpolation_points, cur_reaction->interpolation_adjustment_variables[0], cur_reaction->interpolation_adjustment_variables[1]);
						else
							mol->process(Spect, cur_reaction->interpolation_points[0], cur_reaction->interpolation_points[1], cur_reaction->interpolation_adjustment_variables[0], cur_reaction->interpolation_adjustment_variables[1]);
					}
					else {*/
						if ((cur_reaction->mass[0] != cur_reaction->mass[1]) || (cur_reaction->charge[0] != cur_reaction->charge[1]))
							mol->sort_ions_and_process(Spect);
						else
							mol->process(Spect);
				//	}

					// set molecules as not valid as requested by user
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_DIATOMIC) {
							mol->invalidate(cur_reaction->tag[k]);
						}
					}
					// set ions as not valid as requested by user
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_ION) {
							for(int m=0;m<2;m++) {
								if((cur_reaction->tag[k]->particle_num<0) || ( m == cur_reaction->tag[k]->particle_num)) {
									mol->ion[m]->invalidate(cur_reaction->tag[k]);
								}
							}

							if(!mol->ion[0]->valid || !mol->ion[1]->valid) {
								mol->set_valid(false);								
							}
						}
					}

					// swap ions if requested
					if(cur_reaction->randomize_ions) {
						if((rand() / (float)RAND_MAX)>0.5) {
							ion_class *tmp = mol->ion[1];
							mol->ion[1] = mol->ion[0];
							mol->ion[0] = tmp;
							mol->process_diatomic_only(Spect);
						}
					}

				} else { 
					mol->set_valid(false);
				}
				
				break;

			case(RTYPE_POLYATOMIC):
			case(RTYPE_POLYATOMIC_ELEC):
				
				if (rhit > (cur_reaction->number_of_ions - 0.5) ) //reaction
				{

					big_mol->reset();

					for(int i=0;i<cur_reaction->number_of_ions;i++) {
						r_raw[i].m = cur_reaction->mass[i];
						r_raw[i].q = cur_reaction->charge[i];
						r[i]->set_channel(cur_reaction->channel,true);
						r[i]->set_raw(&r_raw[i]);
						r[i]->set_t_mean(cur_reaction->t_mean[i]);
						r[i]->set_t_width(cur_reaction->t_width[i]);

						if(this->r_fac->raw_order) {				
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);
							r[i]->rotate_raw(this->r_fac->rot_ang);				
						} else {
							r[i]->rotate_raw(this->r_fac->rot_ang);				
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);
						}
						r[i]->calc_phi_pos();
					}
					big_mol->set_ions(cur_reaction->number_of_ions,r);
					big_mol->set_channel(cur_reaction->channel,true);
					if(cur_reaction->incomplete) {
						big_mol->incomplete = true;
						big_mol->ion[cur_reaction->number_of_ions]->channel = cur_reaction->channel;
						big_mol->ion[cur_reaction->number_of_ions]->valid = true;
						big_mol->ion[cur_reaction->number_of_ions]->raw.m = cur_reaction->mass[cur_reaction->number_of_ions];
						big_mol->ion[cur_reaction->number_of_ions]->raw.q = cur_reaction->charge[cur_reaction->number_of_ions];
					}

					if(cur_reaction->use_ion_matrix == 1)
					{
						big_mol->valid = big_mol->run_ion_matrix(cur_reaction,evti,this->Spect);
					}
					else
					{
						big_mol->process(this->Spect);  
						big_mol->valid = true;
					}

					// set molecules as not valid as requested by user
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_POLYATOMIC) {
							big_mol->invalidate(cur_reaction->tag[k]);
								
						}
					}
					int num_i = cur_reaction->number_of_ions;
					if(cur_reaction->incomplete)
						num_i++;
					// set ions as not valid as requested by user	
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_ION) {
							for (int m = 0; m < num_i; m++)
							{
								if((cur_reaction->tag[k]->particle_num<0) || ( m == cur_reaction->tag[k]->particle_num)) {
									big_mol->ion[m]->invalidate(cur_reaction->tag[k]);					
									if (!big_mol->ion[m]->valid)
									{
										big_mol->set_valid(false);
										break;
									}
								}
							}
						}
					}

					// swap ions if requested
					if(cur_reaction->randomize_ions) {
						for (int i=0; i<cur_reaction->rnd_count; i++) {
							int r = rand() % cur_reaction->rnd_count;  // generate a random position in rnd_array
							ion_class *tmp = big_mol->ion[cur_reaction->rnd_array[i]]; 
							big_mol->ion[cur_reaction->rnd_array[i]] =  big_mol->ion[cur_reaction->rnd_array[r]];
							big_mol->ion[cur_reaction->rnd_array[r]] = tmp;
						}
					}

					Coordinate_System molframe = Coordinate_System(big_mol->ion[0]->mom, big_mol->ion[1]->mom);
					big_mol->trafo_to_molframe(molframe);

				} else {
					big_mol->valid = false;
				}
			
				break;//end case(POLYATOMIC)
/*			
				case(RTYPE_POLYATOMIC_ELEC):
				
				if (rhit > (cur_reaction->number_of_ions - 0.5) ) //reaction 
				{

					big_mol->reset();
					for(int i=0;i<cur_reaction->number_of_ions;i++) {
						r_raw[i].m = cur_reaction->mass[i];
						r_raw[i].q = cur_reaction->charge[i];
						r[i]->set_raw(&r_raw[i]);
						r[i]->set_channel(cur_reaction->channel,true);
						r[i]->set_t_mean(cur_reaction->t_mean[i]);
						r[i]->set_t_width(cur_reaction->t_width[i]);
			
						if(this->r_fac->raw_order) {
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);					
							r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);				
						} else {
							r[i]->rotate_raw(cur_reaction->r_fac->rot_ang);				
							//r[i]->shift_stretch_raw(this->r_fac);
							r[i]->shift_stretch_raw(cur_reaction->r_fac);					
						}
						r[i]->calc_phi_pos();
					}
				
					big_mol->set_ions(cur_reaction->number_of_ions,r);
					big_mol->set_channel(cur_reaction->channel,true);
					if(cur_reaction->use_ion_matrix == 1)
					{
						big_mol->valid = big_mol->run_ion_matrix(cur_reaction,evti,this->Spect);
					}
					else
					{
						big_mol->process(this->Spect);  
						big_mol->valid = true;
					}
			
					Coordinate_System molframe = Coordinate_System(big_mol->ion[0]->mom, big_mol->ion[1]->mom);
					big_mol->trafo_to_molframe(molframe);
				// set molecules as not valid as requested by user
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_POLYATOMIC) {
							big_mol->invalidate(cur_reaction->tag[k]);
						}
					}
				// set ions as not valid as requested by user
				//not yet specific for each ion
						
					for(int k=0;k<cur_reaction->num_tags;k++) {
						if(cur_reaction->tag[k]->type == RTYPE_ION) {
							for (int m = 0; m < cur_reaction->number_of_ions; m++)
							{
								big_mol->ion[m]->invalidate(cur_reaction->tag[k]);
														
								if (!big_mol->ion[m]->valid)
								{
									big_mol->set_valid(false);	
									break;
								}
							}
						}
					}
				
				}
				break;//end case(POLYATOMIC_ELEC) */
			default:
				break;
		}
		// has no meaning, so far...
		return true;
	}

}
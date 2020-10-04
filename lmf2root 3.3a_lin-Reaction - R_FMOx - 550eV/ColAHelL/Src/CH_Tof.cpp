#include "OS_Version.h"
#include "CH_Tof.h"
#include "..\CH_Event.h"

#include "Math.h"

namespace CH 
{
	tof_calc_class::tof_calc_class()
	{
		this->tofType = 0;
			
		this->reference_hit = 0;

		this->UseBM = false;
		this->BMchannel = 0;
		this->BMspacing = 0.0;

		this->mod_shift_e = 0.0;
		this->mod_shift_r = 0.0;
		this->mod_shift_p = 0.0;
	}
	tof_calc_class::tof_calc_class(__int32 Type)
	{
		this->tofType = Type;
			
		this->reference_hit = 0;

		this->UseBM = false;
		this->BMchannel = 0;
		this->BMspacing = 0.0;

		this->mod_shift_e = 0.0;
		this->mod_shift_r = 0.0;
		this->mod_shift_p = 0.0;
	}
	tof_calc_class::tof_calc_class(__int32 Type, __int32 BMchannel)
	{
		this->tofType = Type;
			
		this->reference_hit = 0;

		this->UseBM = true;
		this->BMchannel = BMchannel;
		this->BMspacing = 0.0;

		this->mod_shift_e = 0.0;
		this->mod_shift_r = 0.0;
		this->mod_shift_p = 0.0;
	}
	tof_calc_class::tof_calc_class(__int32 Type, __int32 BMchannel, double BMspacing)
	{
		this->tofType = Type;
			
		this->reference_hit = 0;

		this->UseBM = true;
		this->BMchannel = BMchannel;
		this->BMspacing = BMspacing;

		this->mod_shift_e = 0.0;
		this->mod_shift_r = 0.0;
		this->mod_shift_p = 0.0;
	}		

	tof_calc_class::~tof_calc_class()
	{
	}

	double tof_calc_class::tof_of(double a, double b) {
		return a - b;
	}

	double tof_calc_class::mod_tof_of(double a, double b) {
		double res = fmod((a - b),this->BMspacing); 
		if(res<0)res+=this->BMspacing;
		return res;
	}

	void tof_calc_class::set_reference_hit(int hit) {
		this->reference_hit = hit;
	}

	void tof_calc_class::set_tof_type(__int32 type) {
		this->tofType = type;
	}
	void tof_calc_class::set_mod_shift_e(double shift) {
		this->mod_shift_e = shift;
	}
	void tof_calc_class::set_mod_shift_p(double shift) {
		this->mod_shift_p = shift;
	}
	bool tof_calc_class::usesBM() {
		return this->UseBM;
	}
	int tof_calc_class::GetTOFType() {
		return this->tofType;
	}
	double tof_calc_class::GetBMspacing() {
		return this->BMspacing;
	}
	int tof_calc_class::GetBMchannel() {
		return this->BMchannel;
	}
	double tof_calc_class::GetBMtime(double * tdc_ns, __int32 * cnt) {
		if(this->UseBM) {
			if(cnt[this->BMchannel]>1 && this->BMspacing > 0.0)
				return tdc_ns[this->BMchannel*NUM_IONS+1];
			else
				if(cnt[this->BMchannel]>0)
					return tdc_ns[this->BMchannel*NUM_IONS + 0];
		}
		return 0.0;
	}

	void tof_calc_class::calc(CH_event_struct *evt) {
		calc(evt, this->reference_hit);
	}

	double tof_calc_class::get_t_0_shift() {
		if(this->UseBM) {
			if(this->BMspacing>0.0) {
				return this->mod_shift_e;
			}else{
				return 0.0;
			}
		} else {
			return this->mod_shift_p;
		}
	}

	void tof_calc_class::calc(CH_event_struct *evt, int reference_electron) {
		if(this->UseBM) {
			if(this->BMspacing>0.0) {

				// calc electron tofs
				for(int i=0;i<evt->e.num_hits;i++) {
					evt->e.tof[i] = mod_tof_of(evt->e.time[i] + this->mod_shift_e, evt->bunchmarker + 1000.0 * this->BMspacing);
				}
				//calc recoil tofs
				for(int i=0;i<evt->r.num_hits;i++) {
					if(evt->e.num_hits>reference_electron)
						evt->r.tof[i] = tof_of(evt->r.time[i],evt->e.time[reference_electron]) + evt->e.tof[reference_electron];
					else
						evt->r.tof[i] = -10000000000.0;
				}
				//calc proj tofs
				for(int i=0;i<evt->p.num_hits;i++) {
					if(evt->e.num_hits>reference_electron)
						evt->p.tof[i] = tof_of(evt->p.time[i],evt->e.time[reference_electron]) + evt->e.tof[reference_electron];
					else
						evt->p.tof[i] = -10000000000.0;
				}

			} else {

				// calc electron tofs
				for(int i=0;i<evt->e.num_hits;i++) {
					evt->e.tof[i] = tof_of(evt->e.time[i], evt->bunchmarker);
				}
				//calc recoil tofs
				for(int i=0;i<evt->r.num_hits;i++) {
					evt->r.tof[i] = tof_of(evt->r.time[i], evt->bunchmarker);
				}
				//calc proj tofs
				for(int i=0;i<evt->p.num_hits;i++) {
					evt->p.tof[i] = tof_of(evt->p.time[i], evt->bunchmarker);
				}
			
			}

		} else {
			//PROJ is reference
			if(this->tofType==1) { 
				// calc electron tofs
				for(int i=0;i<evt->e.num_hits;i++) {
					evt->e.tof[i] = tof_of(evt->e.time[i], evt->p.time[reference_electron] + this->mod_shift_p);
				}
				//calc recoil tofs
				for(int i=0;i<evt->r.num_hits;i++) {
					evt->r.tof[i] = tof_of(evt->r.time[i], evt->p.time[reference_electron] + this->mod_shift_p);
				}
				//calc proj tofs
				for(int i=0;i<evt->p.num_hits;i++) {
					evt->p.tof[i] = tof_of(evt->p.time[i], evt->p.time[reference_electron] + this->mod_shift_p);
				}
			}
			//ELEC is reference
			if(this->tofType==2) { 
				// calc electron tofs
				for(int i=0;i<evt->e.num_hits;i++) {
					evt->e.tof[i] = tof_of(evt->e.time[i], evt->e.time[reference_electron] + this->mod_shift_e);
				}
				//calc recoil tofs
				for(int i=0;i<evt->r.num_hits;i++) {
					evt->r.tof[i] = tof_of(evt->r.time[i], evt->e.time[reference_electron] + this->mod_shift_e);
				}
				//calc proj tofs
				for(int i=0;i<evt->p.num_hits;i++) {
					evt->p.tof[i] = tof_of(evt->p.time[i], evt->e.time[reference_electron] + this->mod_shift_e);
				}
			}
		}
	}
}
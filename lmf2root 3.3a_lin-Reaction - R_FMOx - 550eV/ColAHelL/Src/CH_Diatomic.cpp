#include "CH_Diatomic.h"
#include "CH_Basics.h"
#include "CH_Ion.h"

#include "Math.h"
namespace CH
{

	diatomic_class::diatomic_class()
	{
		this->tmp_ion = 0; //new ion_class();

		this->mfac = new mom_cor_param();

		this->mfac->dx = 0.0;
		this->mfac->dy = 0.0;
		this->mfac->dz = 0.0;

		this->mfac->x_stretch = 1.0;
		this->mfac->y_stretch = 1.0;
		this->mfac->z_stretch = 1.0;
		this->mfac->overall_stretch = 1.0;

		this->mfac->mir_x = false;
		this->mfac->mir_y = false;
		this->mfac->mir_z = false;
	}

	diatomic_class::~diatomic_class()
	{
		//delete tmp_ion;
		if (mfac) delete mfac;
		mfac = 0;
	}

	void diatomic_class::set_ions(ion_class *ion1, ion_class *ion2) {
		this->ion[0] = ion1;
		this->ion[1] = ion2;
		calc_mu();  
	}

	void diatomic_class::calc_mu() {
		this->reduced_mass = this->ion[0]->raw.m * MASSAU * this->ion[1]->raw.m  * MASSAU / (this->ion[0]->raw.m * MASSAU + this->ion[1]->raw.m * MASSAU);  	
	}

	void diatomic_class::set_mom_rel_calc(int calctype) {
		this->RelMomType = calctype;
	}

	void diatomic_class::set_channel(int channel) {
		this->channel = channel;
	}

	void diatomic_class::set_channel(int channel, bool valid) {
		this->channel = channel;
		this->valid = valid;
	}

	void diatomic_class::set_valid(bool valid) {
		this->valid = valid;
	}

	void diatomic_class::set_mom_fac(mom_cor_param * fac) {
		this->mfac->dx = fac->dx;
		this->mfac->dy = fac->dy;
		this->mfac->dz = fac->dz;

		this->mfac->x_stretch = fac->x_stretch;
		this->mfac->y_stretch = fac->y_stretch;
		this->mfac->z_stretch = fac->z_stretch;
		this->mfac->overall_stretch = fac->overall_stretch;

		this->mfac->mir_x = fac->mir_x;
		this->mfac->mir_y = fac->mir_y;
		this->mfac->mir_z = fac->mir_z;
	}

	void diatomic_class::sort_ions_and_process(spectrometer_class *spect) { //, const std::vector<std::vector<std::vector<double>>> &interpolation_points, const std::vector<double> &adjustment_coeff1, const std::vector<double> &adjustment_coeff2) {
		// find right m/q by checking back to back emission (i.e. cm_mom = 0)

		double temp_m1 = this->ion[0]->raw.m;
		double temp_q1 = this->ion[0]->raw.q;
		double temp_m2 = this->ion[1]->raw.m;
		double temp_q2 = this->ion[1]->raw.q;

		// process first combination
		process(spect); //, interpolation_points[0], interpolation_points[1], adjustment_coeff1, adjustment_coeff2);
		double cm_mom_Ax = this->mom_cm.x;
		double cm_mom_Ay = this->mom_cm.y;
		double cm_mom_Az = this->mom_cm.z;


		// swap particle properties

		this->ion[0]->raw.m = temp_m2;
		this->ion[0]->raw.q = temp_q2;
		this->ion[1]->raw.m = temp_m1;
		this->ion[1]->raw.q = temp_q1;

		// process second combination
		process(spect); //, interpolation_points[1], interpolation_points[0], adjustment_coeff2, adjustment_coeff1);
		double cm_mom_Bx = this->mom_cm.x;
		double cm_mom_By = this->mom_cm.y;
		double cm_mom_Bz = this->mom_cm.z;

		///if(fabs(cm_mom_Ax) < fabs(cm_mom_Bx) 
		if((fabs(cm_mom_Ax) + fabs(cm_mom_Ay) + fabs(cm_mom_Az)) < (fabs(cm_mom_Bx) + fabs(cm_mom_By) + fabs(cm_mom_Bz))) {
			
			// restore initial parameter set, inital asignment was ok
			this->ion[0]->raw.m = temp_m1;
			this->ion[0]->raw.q = temp_q1;
			this->ion[1]->raw.m = temp_m2;
			this->ion[1]->raw.q = temp_q2;

			process(spect); //, interpolation_points[0], interpolation_points[1], adjustment_coeff1, adjustment_coeff2);
		
		} else {
			// swap ions, as we want the ion with inital mass m1 to be ion[0]
			tmp_ion = ion[1];
			ion[1] = ion[0];
			ion[0] = tmp_ion;
			process_diatomic_only(spect);
		}
	}

	void diatomic_class::process(spectrometer_class *spect) { //, const std::vector<std::vector<double>> &tofs_and_moms1, const std::vector<std::vector<double>> &tofs_and_moms2, const std::vector<double> &adjustment_coeff1, const std::vector<double> &adjustment_coeff2) {
		this->ion[0]->process(spect); //, tofs_and_moms1, adjustment_coeff1);
		this->ion[0]->shift_stretch_mom(this->mfac);

		this->ion[1]->process(spect); //, tofs_and_moms2, adjustment_coeff2);
		this->ion[1]->shift_stretch_mom(this->mfac);

		process_diatomic_only(spect);
	}

	void diatomic_class::process_diatomic_only(spectrometer_class *spect) {

		// relative momenta		
		if(spect->ion_side->linear_approximation || spect->ion_side->number_of_regions>2) {
			this->mom_rel = (this->ion[0]->mom - this->ion[1]->mom)/2.0;
		} else { 

			double sr = spect->ion_side->lengths[0] / 1000.0; // in m
			double E = spect->ion_side->Efields[0] * 100; // in V/m
			double m0 = this->ion[0]->raw.m * MUKG; // in kg
			double m1 = this->ion[1]->raw.m * MUKG;
			double q0 = this->ion[0]->raw.q * COULOMB; // in As
			double q1 = this->ion[1]->raw.q * COULOMB;
			double t0 = this->ion[0]->raw.data.tof * 1e-9; // in s
			double t1 = this->ion[1]->raw.data.tof * 1e-9;
			double x0 = this->ion[0]->raw.data.x/1000.0; // in m
			double x1 = this->ion[1]->raw.data.x/1000.0; // in m
			double y0 = this->ion[0]->raw.data.y/1000.0; // in m
			double y1 = this->ion[1]->raw.data.y/1000.0; // in m

			switch(this->RelMomType) {
			case 0:
				this->mom_rel = (this->ion[0]->mom - this->ion[1]->mom)/2.0;
				break;
			case 1:
				this->mom_rel.x = this->reduced_mass * (x0/t0 - x1/t1) * MPSAU;
				this->mom_rel.y = this->reduced_mass * (y0/t0 - y1/t1) * MPSAU;
				this->mom_rel.z = this->reduced_mass * (sr/t0 - sr/t1 - E/2.0 * (t0*q0/m0 - t1*q1/m1)) * MPSAU;
				break;
			case 2:
				this->mom_rel.x = m0*m1*(x0 - x1)/(t1*m0 + t0*m1) * MPSAU;
				this->mom_rel.y = m0*m1*(y0 - y1)/(t1*m0 + t0*m1) * MPSAU;
				this->mom_rel.z = E/2.0 * (pow(t1,2)*q1*m0 - pow(t0,2)*q0*m1) / (t0*m1 + t1*m0) * MPSAU;
				break;
			}
		}
		
		// cm momenta
		this->mom_cm = (this->ion[0]->mom + this->ion[1]->mom);
	}

	double diatomic_class::KER() {
		return this->mom_rel.Mag() * this->mom_rel.Mag() / (2.0 * this->reduced_mass) * EVAU;
	}

	bool diatomic_class::check_mom_cm_x(double psumx_width, int channel, bool invalidate)
	{
		if(fabs(this->mom_cm.x)>psumx_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		return true;
	}

	bool diatomic_class::check_mom_cm_y(double psumy_width, int channel, bool invalidate)
	{
		if(fabs(this->mom_cm.y)>psumy_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		return true;
	}

	bool diatomic_class::check_mom_cm_z(double psumz_width, int channel, bool invalidate)
	{
		if(fabs(this->mom_cm.z)>psumz_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		return true;
	}

	bool diatomic_class::check_mom_cm(double psumx_width, double psumy_width, double psumz_width, int channel, bool invalidate)
	{	
		if(!(check_mom_cm_x(psumx_width,channel,invalidate)))
			return false;

		if(!(check_mom_cm_y(psumy_width,channel,invalidate)))
			return false;

		if(!(check_mom_cm_z(psumz_width,channel,invalidate)))
			return false;

		if(channel>=0)
			this->channel = channel;

		return true;
	}


	bool diatomic_class::check_KER(double energy_from, double energy_to, int channel, bool invalidate)
	{
		if((this->KER() < energy_from) || (this->KER() > energy_to))
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(channel>=0)
			this->channel = channel;

		return true;
	}

	bool diatomic_class::check_validity(double var, double min, double max, bool negate)
	{
		if(negate) {
			if(var < min || var >= max)
				return false;
			else 
				return this->valid;
		} else {
			if(var > min && var <= max)
				return false;
			else 
				return this->valid;
		}
	}

	void diatomic_class::invalidate(rtag_struct* rtag) 
	{
		switch (rtag->tag_type) {
		
		case RTAG_PSUMX:
			this->valid = check_validity(this->mom_cm.x,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUMY:
			this->valid = check_validity(this->mom_cm.y,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUMZ:
			this->valid = check_validity(this->mom_cm.z,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUM:
			this->valid = check_validity(this->mom_cm.Mag(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUM_PHI_DEG:
			this->valid = check_validity(this->mom_cm.Phi_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PSUM_THETA_DEG:
			this->valid = check_validity(this->mom_cm.Theta_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PRELX:
			this->valid = check_validity(this->mom_rel.x,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PRELY:
			this->valid = check_validity(this->mom_rel.y,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PRELZ:
			this->valid = check_validity(this->mom_rel.z,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PREL:
			this->valid = check_validity(this->mom_rel.Mag(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PREL_PHI_DEG:
			this->valid = check_validity(this->mom_rel.Phi_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PREL_THETA_DEG:
			this->valid = check_validity(this->mom_rel.Theta_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_KER:
			this->valid = check_validity(this->KER(),rtag->min,rtag->max,rtag->negate);
			break;

		default:
			break;
		}

	}
}
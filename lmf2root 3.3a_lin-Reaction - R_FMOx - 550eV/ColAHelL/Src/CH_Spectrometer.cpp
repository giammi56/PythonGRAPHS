#include "CH_Spectrometer.h"

namespace CH
{

	spectrometer_class::spectrometer_class(){
		this->ion_side = new spec_arm();

		this->ion_side->linear_approximation = false;
//		this->ion_side->interpolation_approximation = false;
		this->ion_side->number_of_regions = 0;
		
		for(int i=0;i<MAX_NUM_REGIONS;i++) {
			this->ion_side->Efields[i] = 0.0;
			this->ion_side->lengths[i] = 0.0;
		}

		this->electron_side = new spec_arm();

		this->electron_side->linear_approximation = false;
//		this->ion_side->interpolation_approximation = false;
		this->electron_side->number_of_regions = 0;

		for(int i=0;i<MAX_NUM_REGIONS;i++) {
			this->electron_side->Efields[i] = 0.0;
			this->electron_side->lengths[i] = 0.0;
		}

		this->Bfield_G = 0;
		this->Bfield_clockwise = false;

		this->VJet = 0.0;
		this->AngJet = 0.0;
		
		this->MeanTOFe = 0.0;
	}

	spectrometer_class::~spectrometer_class(){
		if (ion_side) delete ion_side; ion_side = 0;
		if (electron_side) delete electron_side; electron_side = 0;
	}

	void spectrometer_class::set_Bfield(double B_val, bool B_clockwise, bool B_in_ns) {
		this->Bfield_clockwise = B_clockwise;
		if(B_in_ns) {
			this->Bfield_ns = B_val;
			this->Bfield_G = 357.238755341 / B_val; 
		} else {
			if(B_val <1e-6) // catch very tiny (i.e. zero) B-fields
				this->Bfield_ns = 1e+8;
			else
				this->Bfield_ns = 357.238755341 / B_val ;
			this->Bfield_G = B_val;
		}
	}

	void spectrometer_class::set_VJet(double vjet, double angjet) {
		this->VJet = vjet;
		this->AngJet = angjet;
	}

	void spectrometer_class::set_ion_arm_lin(double E_val) {	
		this->ion_side->linear_approximation = true;
		this->ion_side->Efields[0] = E_val;
	}

/*	void spectrometer_class::set_ion_arm_interpolation() {
		this->ion_side->interpolation_approximation = true;
	}
*/
	void spectrometer_class::set_ion_arm(unsigned short num_regions, double *lengths, double *E_val) {
		this->ion_side->linear_approximation = false;
		this->ion_side->number_of_regions = num_regions;

		for(int i=0;i<num_regions;i++) {
			this->ion_side->Efields[i] = E_val[i];
			this->ion_side->lengths[i] = lengths[i];
		}
	}

	void spectrometer_class::set_electron_arm(unsigned short num_regions, double *lengths, double *E_val) {
		this->electron_side->linear_approximation = false;
		this->electron_side->number_of_regions = num_regions;

		for(int i=0;i<num_regions;i++) {
			this->electron_side->Efields[i] = E_val[i];
			this->electron_side->lengths[i] = lengths[i];
		}
	}

	void spectrometer_class::set_electron_arm_lin(double E_val) {	
		this->electron_side->linear_approximation = true;
		this->electron_side->Efields[0] = E_val;
	}

	void spectrometer_class::set_mean_tof_e(double mean_tof) {	
		this->MeanTOFe = mean_tof;
	}

}
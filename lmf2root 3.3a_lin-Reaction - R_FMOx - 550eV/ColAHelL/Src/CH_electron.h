#ifndef CH_ELEC_ALREADY_INCLUDED
	#define CH_ELEC_ALREADY_INCLUDED

	#include "CH_Particle.h"
	#include "CH_Spectrometer.h"
	#include "CH_Basics.h"

	namespace CH
	{
		class electron_class : public particle_class {  
			
		public:

			electron_class();
			electron_class(particle_raw rawdata);
			~electron_class();

			void set_raw(particle_raw rawdata);

			void cut_wiggle(spectrometer_class * spect, double width);
			void within_range_of_other(electron_class *ref_elec, double tof_from, double tof_to, int channel, bool invalidate);
			void EB_pos_corr(double ex_k1, double ey_k1, double et_k1,double ex_k2, double ey_k2, double et_k2);

		}; // end electron class
	}

#endif;
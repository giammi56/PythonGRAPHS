#ifndef CH_PARTICLE_ALREADY_INCLUDED
	#define CH_PARTICLE_ALREADY_INCLUDED

	#include "CH_Basics.h"
	#include "CH_Spectrometer.h"
	#include "CH_Reaction.h"

	namespace CH
	{
		
		class particle_class { 

		public:
			bool valid;
			int channel;
			
			particle_raw raw;
			CH_vector mom;       // in a.u.
			double energy();

			cor_param cor;

		public:

			particle_class();
			~particle_class();
			
			void set_raw(particle_raw *rawdata);
			void set_raw(double x, double y, double t, int method = 0);
			void calc_phi_pos();

			void set_channel(int channel);
			void set_channel(int channel, bool valid);
			void set_valid(bool valid);
			void set_t_mean(double t);
			void set_t_width(double t);

			void shift_stretch_raw(cor_param * fac);
			void rotate_raw(double ang);

			void shift_stretch_mom(mom_cor_param * fac);

			void process(spectrometer_class *spect); //, const std::vector<std::vector<double>> &tofs_and_moms = std::vector<std::vector<double>>(), const std::vector<double> & adjustment_coeff = std::vector<double>(5,0) ); //adjustment_coeff are initlized to zero (so they will not change pz)
			void particle_class::process(spectrometer_class *spect, double x, double y, double t); //, const std::vector<std::vector<double>> &tofs_and_moms = std::vector<std::vector<double>>(), const std::vector<double> & adjustment_coeff = std::vector<double>(5,0) ); // for ion matrix (polyatomic)
			void autotune();

			bool check_tof(double tof_from, double tof_to, int channel, bool invalidate);
			bool check_mom(double px_width, double py_width, double pz_width, int channel, bool invalidate);
			bool check_energy(double energy_from, double energy_to, int channel, bool invalidate);
			void invalidate(rtag_struct* rtag);

		private:
			bool check_validity(double var, double min, double max, bool negate);
			double t_mean;
			double t_width;

			friend class electron_class;
			friend class ion_class;
			friend class polyatomic_class;
		}; // end particle class
	}
#endif;
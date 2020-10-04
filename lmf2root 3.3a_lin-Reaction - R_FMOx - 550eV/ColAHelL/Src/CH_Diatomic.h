#ifndef CH_DIATOMIC_ALREADY_INCLUDED
	#define CH_DIATOMIC_ALREADY_INCLUDED

	#include "CH_Basics.h"
	#include "CH_Ion.h"
	
	#define RELMOMSTD 0;
	#define RELMOMNOP 1;
	#define RELMOMNOR 2;

	namespace CH
	{
		class diatomic_class {  

		private:

			double reduced_mass;
			mom_cor_param * mfac;
			ion_class * tmp_ion;
			void calc_mu();
			bool check_validity(double var, double min, double max, bool negate);
			int RelMomType;

		public:

			diatomic_class();
			~diatomic_class();
			
			bool valid;
			int channel;

			ion_class *ion[2];

			CH_vector mom_rel;
			CH_vector mom_cm;
			double KER();

			void set_ions(ion_class *ion1, ion_class *ion2);

			void set_channel(int channel);
			void set_channel(int channel, bool valid);
			void set_valid(bool valid);

			void shift_stretch_rel(cor_param * fac);

			void set_mom_fac(mom_cor_param * fac);

			void set_mom_rel_calc(int calctype);

			void sort_ions_and_process(spectrometer_class *spect); //, const std::vector<std::vector<std::vector<double>>> &interpolation_points = std::vector<std::vector<std::vector<double>>>(0, std::vector<std::vector<double>>() ), const std::vector<double> & adjustment_coeff1 = std::vector<double>(5, 0), const std::vector<double> & adjustment_coeff2 = std::vector<double>(5, 0));
			void process(spectrometer_class *spect); //, const std::vector<std::vector<double>> &tofs_and_moms1 = std::vector<std::vector<double>>(), const std::vector<std::vector<double>> &tofs_and_moms2 = std::vector<std::vector<double>>(), const std::vector<double> & adjustment_coeff1 = std::vector<double>(5, 0), const std::vector<double> & adjustment_coeff2 = std::vector<double>(5, 0));
			void process_diatomic_only(spectrometer_class *spect);
		
			bool check_mom_cm_x(double psumx_width, int channel, bool invalidate);
			bool check_mom_cm_y(double psumy_width, int channel, bool invalidate);
			bool check_mom_cm_z(double psumz_width, int channel, bool invalidate);
			bool check_mom_cm(double psumx_width, double psumy_width, double psumz_width, int channel, bool invalidate);

			bool check_KER(double energy_from, double energy_to, int channel, bool invalidate);

			void invalidate(rtag_struct* rtag);
			
		}; // end diatomic class
	}
#endif;



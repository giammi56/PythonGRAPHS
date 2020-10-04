#ifndef CH_SPEC_ALREADY_INCLUDED
	#define CH_SPEC_ALREADY_INCLUDED

	namespace CH
	{

		#define MAX_NUM_REGIONS 3
		
		struct spec_arm{
			bool linear_approximation;
//			bool interpolation_approximation;

		#ifdef FLOAT_ONLY
			float number_of_regions;
			float lengths[MAX_NUM_REGIONS];
			float Efields[MAX_NUM_REGIONS];
		#else
			double number_of_regions;
			double lengths[MAX_NUM_REGIONS];
			double Efields[MAX_NUM_REGIONS];
		#endif

		};

		class spectrometer_class{

		public:

			spec_arm * ion_side;
			spec_arm * electron_side;

			#ifdef FLOAT_ONLY
				float Bfield_G;
				float Bfield_ns;
			#else
				double Bfield_G;
				double Bfield_ns;
			#endif

			bool Bfield_clockwise;
			double VJet;
			double AngJet;
			double MeanTOFe;

		public:
				
			spectrometer_class();
			~spectrometer_class();

			void set_Bfield(double B_val, bool B_clockwise, bool B_in_ns);
			void set_VJet(double vjet, double angjet);
			void set_ion_arm_lin(double E_val);
//			void set_ion_arm_interpolation();
			void set_ion_arm(unsigned short num_regions, double *lengths, double *E_val);
			void set_electron_arm(unsigned short num_regions, double *lengths, double *E_val);
			void set_electron_arm_lin(double E_val);
			void set_mean_tof_e(double mean_tof);

			friend class particle_class;
			friend class electron_class;
			friend class presorter_class;
			
		}; // end spectrometer class
	}

#endif;
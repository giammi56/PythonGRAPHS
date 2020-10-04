#ifndef CH_ION_ALREADY_INCLUDED
	#define CH_ION_ALREADY_INCLUDED

	#include "CH_Particle.h"
	#include "CH_Spectrometer.h"
	#include "CH_Basics.h"
	

	namespace CH
	{
		
		class ion_class : public particle_class { 
			
		public:

			ion_class();
			~ion_class();
		
			void clone_from(ion_class * ion);

		friend class diatomic_class;
		friend class polyatomic_class;
		}; // end ion class
	}

#endif;
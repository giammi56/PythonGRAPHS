#ifndef CH_RANGE_ALREADY_INCLUDED
	#define CH_RANGE_ALREADY_INCLUDED
	
	#include "CH_Basics.h"
	#include <string>	
	#include <vector>

	using namespace std;

	namespace CH
	{
		
		// particle types for defining default ranges
		#define NUM_PARTICLE_TYPES 4 
		#define RA_ELECTRON		0
		#define RA_ION			1
		#define RA_PROJECTILE	2
		#define RA_DIATOMIC		3

		// variables for single particles (ions, elecs...)
		#define NUM_SINGLE_VARS 6 
		#define TOF		0
		#define P		1
		#define PX		2
		#define PY		3
		#define PZ		4
		#define ENERGY	5

		// variables for diatomic molecules
		#define NUM_DIATOMIC_VARS 9 
		#define PREL	0
		#define PRELX	1
		#define PRELY	2
		#define PRELZ	3
		#define PCM		4
		#define PCMX	5
		#define PCMY	6
		#define PCMZ	7
		#define KINER   8
		
		class range_class { 
		
		private:
			double *range_map_from;
			double *range_map_to;
			int array_size[NUM_PARTICLE_TYPES];
			void init_default_values();
			void prepare_tags();

		public:
			range_class();
			~range_class();

			vector <string> particle_types;	
			vector <string> single_vars;
			vector <string> diatomic_vars;

			double get_from(int particle_type, int var);
			double get_to(int particle_type, int var);
			bool set_from(int particle_type, int var, double val);
			bool set_to(int particle_type, int var, double val);

		};

	}
#endif;
#ifndef CH_BASICS_ALREADY_INCLUDED
	#define CH_BASICS_ALREADY_INCLUDED

	#include <vector>
	#include "CH_vector.h"

	#define PAU 1.99285E-24				// Momentum in AU
	#define MASSAU 1822.888				// Mass in AU
	#define VAU 2.1877e+6				// Velocity in AU
	#define EVAU 27.2114				// eV in AU
	#define EFIELDAU 0.00805			// Electric field (V/cm) in AU
	#define PI 3.14159265				// PI
	#define MEKG 9.1093897E-31			// Electron Mass in Kg
	#define MUKG 1.66053886e-27			// u in Kg
	#define COULOMB	1.60217733E-19		// Elementary charge
	#define MPSAU 4.5706e-7				// 1 m/s in AU

	#define MIR_NONE 0;
	#define MIR_X 1;
	#define MIR_Y 2;
	#define MIR_XY 3;

	namespace CH
	{

		struct cor_param{ 

		#ifdef FLOAT_ONLY
			float dx;
			float dy;
			float dt;

			float x_stretch;
			float y_stretch;
			float overall_stretch;

			bool mir_x;
			bool mir_y;

			float rot_ang;
			float EBx;
			float EBy;

			bool raw_order;
		#else
			double dx;
			double dy;
			double dt;

			double x_stretch;
			double y_stretch;
			double overall_stretch;

			bool mir_x;
			bool mir_y;

			double rot_ang;
			double EBx;
			double EBy;

			bool raw_order;
		#endif
		};

		struct mom_cor_param{ 

		#ifdef FLOAT_ONLY
			float dx;
			float dy;
			float dz;

			float x_stretch;
			float y_stretch;
			float z_stretch;
			float overall_stretch;

			bool mir_x;
			bool mir_y;
			bool mir_z;

			float rot_ang;
		#else
			double dx;
			double dy;
			double dz;

			double x_stretch;
			double y_stretch;
			double z_stretch;
			double overall_stretch;

			bool mir_x;
			bool mir_y;
			bool mir_z;

			double rot_ang;
		#endif
		};


		struct xyt{

		#ifdef FLOAT_ONLY
			float x; // in mm
			float y; // in mm
			float t; // in ns
		#else
			double x; // in mm
			double y; // in mm
			double t; // in ns
		#endif

		};

		struct CH_data_struct{
			double x;
			double y;
			double time;
			double tof;
			double mcp;
		};

		struct particle_raw{
			double method;    // this is a double to make it work with root!

		#ifdef FLOAT_ONLY
			float m; // in a.u.
			float q; // in a.u.
			CH_data_struct data; // tof and position in ns and mm
			
			float phi; // angle on detector
		#else
			double m; // in a.u.
			double q; // in a.u.
			CH_data_struct data; // tof and position in ns and mm
			double phi; // angle on detector
		#endif

		};

	}
#endif;


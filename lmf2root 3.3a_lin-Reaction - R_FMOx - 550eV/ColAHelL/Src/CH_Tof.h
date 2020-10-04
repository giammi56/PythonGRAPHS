#ifndef CH_TOF_ALREADY_INCLUDED
	#define CH_TOF_ALREADY_INCLUDED


	#include "CH_Basics.h"
	#include "..\CH_Event.h"

	namespace CH 
	{
		class tof_calc_class{  // find out how to inherit from electron struct!

		private:
			__int32 tofType;    // 0: BM, 1:PROJ. 2:ELEC
			
			__int32 reference_hit;

			bool UseBM;
			__int32	BMchannel;
			double BMspacing;

			double mod_shift_e;
			double mod_shift_r;
			double mod_shift_p;

			double tof_of(double a, double b);
			double mod_tof_of(double a, double b);

		public:
			tof_calc_class();
			tof_calc_class(__int32 Type);
			tof_calc_class(__int32 Type, __int32 BMchannel);
			tof_calc_class(__int32 Type, __int32 BMchannel, double BMspacing);
			~tof_calc_class();

			void set_reference_hit(__int32 hit);
			void set_tof_type(__int32 type);
			void set_mod_shift_e(double shift);
			void set_mod_shift_p(double shift);
			bool usesBM();
			int GetTOFType();
			double GetBMspacing();
			int GetBMchannel();
			double GetBMtime(double * tdc_ns, __int32 *cnt);
			double get_t_0_shift();

			void calc(CH_event_struct *evt);
			void calc(CH_event_struct *evt, int reference_electron);

			friend class presorter_class;
		};
	}
#endif;
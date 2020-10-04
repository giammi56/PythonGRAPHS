#ifndef CH_LOWLEVEL_ALREADY_INCLUDED
	#define CH_LOWLEVEL_ALREADY_INCLUDED

	namespace CH
	{
		double ps_get_y_at_x(double * y, __int32 number_of_points, double x);

//		double interpolation2mom(double tof, const std::vector<std::vector<double>> &tofs_and_moms);

//		std::vector<std::vector<double>> read_array_file(std::string filename);

		double NEW_v(double acc1, double acc2, double acc3, double acc1_length, double acc2_length, double acc3_length, double v, double t);

		double tof2mom_3accel(double tof_ns, double acc1mm, double acc2mm, double acc3mm, double Efield1_Vpcm, double Efield2_Vpcm, double Efield3_Vpcm, double charge_au, double mass_amu);

		double calc_px(double tof_ns, double x_mm, double y_mm, double mass_amu, double charge_au, double BField_ns, bool BField_clockwise = true);

		double calc_py(double tof_ns, double x_mm, double y_mm, double mass_amu, double charge_au, double BField_ns, bool BField_clockwise = true);

		double t2_3accel(double tof_ns, double acc1mm, double acc2mm, double acc3mm, double Efield1_Vpcm, double Efield2_Vpcm, double Efield3_Vpcm, double charge1_au, double charge2_au, double mass1_amu, double mass2_amu);
	}

#endif
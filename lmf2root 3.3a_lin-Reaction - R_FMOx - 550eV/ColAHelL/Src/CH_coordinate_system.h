#pragma once
#include "CH_vector.h"


class Coordinate_System{
public:
	CH_vector x_axis;
	CH_vector y_axis;
	CH_vector z_axis;

	Coordinate_System();
	Coordinate_System(CH_vector v1, CH_vector v2);
	Coordinate_System(CH_vector x, CH_vector y, CH_vector z);
	//GN
	Coordinate_System(CH_vector v1, CH_vector v2, int);
	//GN_!
	//~Coordinate_System();


	CH_vector Coordinate_System::project_vector( CH_vector v);
};
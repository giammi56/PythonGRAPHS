#pragma once
#include "CH_coordinate_system.h"
#include "CH_vector.h"

Coordinate_System::Coordinate_System(void){
}

//Coordinate_System::~Coordinate_System(void)
//{
//	delete x_axis;
//}

//Molecular frame transform
//Give it the momentum vectors of two particles and it will return three unit vectors that are defined in the molecules frame.
//Z will be exactly align with vector 1.
//Y will be determined by crossing vector 1 and vector 2. 
//X will be determined by crossing Y and vector 1.
Coordinate_System::Coordinate_System(CH_vector v1, CH_vector v2){
	this->z_axis = v1.Norm();
	this->y_axis = v1.Cross(v2).Norm();
	this->x_axis =  this->y_axis.Cross(v1).Norm();
}

Coordinate_System::Coordinate_System(CH_vector x, CH_vector y, CH_vector z){
	x_axis=x; 
	y_axis=y; 
	z_axis=z;
}

//GN
//Molecular frame transform according to Kilian: this system has the modified cos(theta) and phi fucntions
//Give it the momentum vectors of two particles and it will return three unit vectors that are defined in the molecules frame.
//X will be exactly align with vector 1.
//Z will be determined by crossing vector 1 and vector 2. 
//Y will be determined by crossing vector 1 and Z.
//This procedure leads to the same results as labframe_transformation when used in combiantion with cos(theta)_mod and phi_mod
Coordinate_System::Coordinate_System(CH_vector v1, CH_vector v2, int) {
	this->x_axis = v1.Norm();
	this->z_axis = v1.Cross(v2).Norm();
	//this->y_axis = v1.Cross(this->z_axis).Norm();
	//The cross product is anticommutative (i.e., a x b = - b x a)/
	this->y_axis = v1.Cross(this->z_axis).Norm();
	//GN_!
}

CH_vector Coordinate_System::project_vector( CH_vector v){
	CH_vector pro;
	pro.x=v.Dot(this->x_axis);
	pro.y=v.Dot(this->y_axis);
	pro.z=v.Dot(this->z_axis);
	return pro;
}

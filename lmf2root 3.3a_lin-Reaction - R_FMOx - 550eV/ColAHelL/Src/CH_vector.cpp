#include "CH_vector.h"
#include <math.h>
#include <stdio.h>
#include <iostream> 
#include <string>
#include <sstream> //double to string

#define PI (3.14159265358979323846264338327950288419716939937510)

CH_vector::CH_vector(void)
{
}

CH_vector::CH_vector(double X, double Y, double Z){
	x=X;
	y=Y;
	z=Z;
};
CH_vector::~CH_vector(void)
{
}

//returns the cross product of the two vectors
CH_vector CH_vector::Cross(CH_vector b){
	//returns the cross product of the two vectors
	CH_vector c = CH_vector();
	c.x=(this->y*b.z-this->z*b.y);
	c.y=(this->z*b.x-this->x*b.z);
	c.z=(this->x*b.y-this->y*b.x);
	return c;
}

//returns the lenght of the vector
double CH_vector::Mag(){
	return sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
}


//find how much of a[2] is perpendicular to b[2]
double CH_vector::Perp(CH_vector b){
//	|axb|	|a||b|sin(alpha)
//	----- = ---------------- = a_perp
//	 |b|		  |b|
//	
//	   |   a
//	   |  /		
//	   | /	
//	   |/
//	   |---------b	
//

	CH_vector cross_product = this->Cross(b);
	return cross_product.Mag()/b.Mag();
}

//returns the dot product of the two vectors
double CH_vector::Dot(CH_vector b){
	return (this->x*b.x + this->y*b.y + this->z*b.z);
}

//returns a unit vector
CH_vector CH_vector::Norm(){
	CH_vector c;
	c.x=this->x/this->Mag();
	c.y=this->y/this->Mag();
	c.z=this->z/this->Mag();
	return c;
}

double CH_vector::Parallel(CH_vector b){
//	|a.b|	|a||b|cos(alpha)
//	----- = ---------------- = a_parallel
//	 |b|		  |b|
//	
//	   |   a
//	   |  /		
//	   | /	
//	   |/
//	   |---------b		
//
	return this->Dot(b)/b.Mag();
}


//returns the phi angle of the vector (phi is defined from x and is in the xy plane, range +-PI, x=0)
double CH_vector::Phi(){
	double a_phi = atan2(this->y, this->x); //radians
	//if(a_phi<0.0)	a_phi =2.0*PI + a_phi;
	return a_phi;
}

//returns the phi angle of the vector in deg (phi is defined from x and is in the xy plane, range +-180°, x=0°)
double CH_vector::Phi_deg(){
	double a_phi = atan2(this->y, this->x); //radians
	//if(a_phi<0.0)	a_phi =2.0*PI + a_phi;
	return a_phi / PI * 180.0;
}

double CH_vector::Phi_deg_mod(){
	double a_phi = atan2((this->y)*-1, this->z); //radians
	//if(a_phi<0.0)	a_phi =2.0*PI + a_phi;
	return a_phi / PI * 180.0;
}

double CH_vector::Phi_deg_mod_mirr(){
	double a_phi = atan2((this->y)*-1, this->z); //radians
    if(a_phi>180.)	a_phi =a_phi-360.;
	return a_phi / PI * 180.0;
}

//returns the theta angle of the vector (theta is defined from z, range -1..1) 
double CH_vector::Cos_Theta(){
	//return fabs(atan(this->z / sqrt(this->x*this->x + this->y*this->y)) - PI/2);	//radians;
	return this->z / sqrt(this->x*this->x + this->y*this->y + this->z*this->z);	//radians;
}

//GN
double CH_vector::Cos_Theta_mod(){
	//return fabs(atan(this->z / sqrt(this->x*this->x + this->y*this->y)) - PI/2);	//radians;
	return this->x / sqrt(this->x*this->x + this->y*this->y + this->z*this->z)*-1;	//radians;
}
//GN_!

//returns the theta angle of the vector (theta is defined from z, range 0..PI) 
double CH_vector::Theta(){
	//return fabs(atan(this->z / sqrt(this->x*this->x + this->y*this->y)) - PI/2);	//radians;
	return acos(this->z / sqrt(this->x*this->x + this->y*this->y + this->z*this->z));	//radians;
}

//returns the theta angle of the vector (theta is defined from z, range 0°..180°) 
double CH_vector::Theta_deg(){
	//return fabs(atan(this->z / sqrt(this->x*this->x + this->y*this->y)) - PI/2);	//radians;
	return acos(this->z / sqrt(this->x*this->x + this->y*this->y + this->z*this->z)) / PI * 180.0;	//radians;
}

//the angle between two vectors in radians
double CH_vector::Angle(CH_vector b){
	double cos = this->Dot(b) / (this->Mag()*b.Mag() ); 
	if ( cos >  1.0 ) {
		return  acos( 1. );
	}
	if ( cos <  -1.0 ) {
		return  acos( -1. );
	}
	return  acos( cos );
}

//the angle between two vectors in degrees
double CH_vector::Angle_deg(CH_vector b){
	double cos = this->Dot(b) / ( this->Mag()*b.Mag() ); 
	if ( cos >  1.0 ) {
		return  acos( 1. )*180./PI;
	}
	if ( cos <  -1.0 ) {
		return  acos( -1. )*180./PI;
	}
	return  acos( cos )*180./PI;
}


//Rodrigues' rotation
//rotates v about k by angle alpha ( alpha in radians)
CH_vector CH_vector::rotate_about_k(CH_vector k_axis, double alpha){ 
	//v_rot = v cos(alpha) + (k_unit x v)sin(alpha) + k_unit(k_unit.v)(1-cos(alpha))
		
	CH_vector k_unit = k_axis.Norm();
	
	CH_vector result = *this * cos(alpha) + k_unit.Cross(*this)*sin(alpha) + k_unit * (k_unit.Dot(*this)) * (1-cos(alpha));

	return result;	
}


//adds to vectors together and returns a vector
CH_vector CH_vector::operator+( const CH_vector& other ) const{
	CH_vector result=CH_vector();
	result.x=this->x + other.x;
	result.y=this->y + other.y;
	result.z=this->z + other.z;
	return result;
}

//subtracts to vectors and returns a vector
CH_vector CH_vector::operator-( const CH_vector& other ) const{
	CH_vector result=CH_vector();
	result.x=this->x - other.x;
	result.y=this->y - other.y;
	result.z=this->z - other.z;
	return result;
}

//multiplies a vectors by a scaler together and returns a vector
CH_vector CH_vector::operator*( double scaler ) const{
	CH_vector result=CH_vector();
	result.x=this->x * scaler;
	result.y=this->y * scaler;
	result.z=this->z * scaler;
	return result;
}

//divides a vectors by a scaler together and returns a vector
CH_vector CH_vector::operator/( double scaler ) const{
	CH_vector result=CH_vector();
	result.x=this->x / scaler;
	result.y=this->y / scaler;
	result.z=this->z / scaler;
	return result;
}


/*void CH_vector::operator+=( const CH_vector& other ) {
	
	this->x += other.x;
	this->y += other.y;
	this->z += other.z;
	
}
*/
void CH_vector::print(){
	std::cout << "X=" << this->x << ", Y="<<this->y << ", Z="<<this->z << std::endl;
}


std::string CH_vector::to_string(){
	std::ostringstream s;
	s << "X=" << this->x << ", Y="<<this->y << ", Z="<<this->z << std::endl;
	return s.str();
}



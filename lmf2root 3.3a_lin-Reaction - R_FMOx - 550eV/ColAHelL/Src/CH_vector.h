#pragma once
#include <string>
//#include "CH_coordinate_system.h"




//lowlevel vector class
class CH_vector
{
public:
	double x;
	double y;
	double z;
	CH_vector(void);
	CH_vector(double X, double Y, double Z);
	~CH_vector(void);

	CH_vector Cross(CH_vector b); 
	double CH_vector::Mag();
	double CH_vector::Perp( CH_vector b);
	
	double CH_vector::Dot(CH_vector b);
	CH_vector CH_vector::Norm();
	double CH_vector::Parallel(CH_vector b);
	double CH_vector::Phi();
	double CH_vector::Phi_deg();
	//GN
	double CH_vector::Phi_deg_mod();
	double CH_vector::Phi_deg_mod_mirr();
	double CH_vector::Cos_Theta();
	double CH_vector::Cos_Theta_mod();
	//GN_!
	double CH_vector::Theta();
	double CH_vector::Theta_deg();
	double CH_vector::Angle(CH_vector b);
	double CH_vector::Angle_deg(CH_vector b);

	CH_vector CH_vector::rotate_about_k(CH_vector k_axis, double alpha);


	CH_vector operator+( const CH_vector& other ) const;
	CH_vector operator-( const CH_vector& other ) const;

	CH_vector operator*( double scaler ) const;
	CH_vector operator/( double scaler ) const;
	//CH_vector operator+( const CH_vector& other ) const;
	//CH_vector operator*(  const CH_vector& other) const;

	void CH_vector::print();
	std::string CH_vector::to_string();


};


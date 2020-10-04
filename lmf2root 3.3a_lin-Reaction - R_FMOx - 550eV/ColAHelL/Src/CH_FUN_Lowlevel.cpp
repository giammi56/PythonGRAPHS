#include "Math.h"
#include "float.h"
#include "CH_Basics.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


namespace CH
{

// Achim's magic "parabola spline"...
// --------------------------------------------------------------------------------
	double ps_parabola_fit(double y1, double y2, double y3, double x_pos)
	{
		double a = 0.5*(y1+y3) - y2;
		double b = 2*y2 - 1.5*y1 - 0.5*y3;
		return (a*x_pos + b)*x_pos + y1;
	}

	double ps_get_y_at_x(double * y, __int32 number_of_points, double x)
	{
		if (x < 0. || x > number_of_points) return 1.e200;
		if (number_of_points < 3) return 1.e200;

		if (x <= 1.) {
			return ps_parabola_fit(y[0], y[1], y[2], x);
		}
		if (x >= number_of_points-2) {
			return ps_parabola_fit(y[number_of_points-3], y[number_of_points-2], y[number_of_points-1], x-number_of_points+3);
		}
		__int32 index = __int32(x);
		double dx = x-index;
		double yi = y[index];
		double yip1 = y[index+1];
		double y1 = ps_parabola_fit(y[index-1], yi, yip1, dx+1.);
		double y2 = ps_parabola_fit(yi, yip1, y[index+2], dx);

		return y1*(1.-dx) + y2*dx;
	} 
// --------------------------------------------------------------------------------



//	double interpolation2mom(double tof, const std::vector<std::vector<double>> &tofs_and_moms) {
//		int a = 0;
//		int b = tofs_and_moms[0].size()-1;
//		int mid = (b-a)/2;
//		double result = -123456789.;
//		if (tof > tofs_and_moms[0][0] && tof < tofs_and_moms[0][b]) {
//			while ((b - a) != 1) { //search by halfs
//				mid = (b - a) / 2 + a;
//				if (tof > tofs_and_moms[0][mid] && tof < tofs_and_moms[0][b])
//					a = mid;
//				else
//					b = mid;
//			}
//			result = (tof- tofs_and_moms[0][a])*(tofs_and_moms[1][b] - tofs_and_moms[1][a]) /(tofs_and_moms[0][b] - tofs_and_moms[0][a]) + tofs_and_moms[1][a];
//		}//if it is not inside the range of points the we will us a linear approximation
//		else if (tofs_and_moms[0].size() > 2){
//			result = (tof - tofs_and_moms[0][a])*(tofs_and_moms[1][b] - tofs_and_moms[1][a]) / (tofs_and_moms[0][b] - tofs_and_moms[0][a]) + tofs_and_moms[1][a];
//		}
//		return result;
//	}
//
//
//	std::vector<std::vector<double>> read_array_file(std::string filename) {
//		std::vector<std::vector<double>> results;
//		std::vector<double> push_temp;
//		results.push_back(push_temp);
//		results.push_back(push_temp);
//
//		std::ifstream array_file(filename);
//		std::string line;
//		std::string comment = "//";
//		std::string cur_open = "{";
//		std::string cur_close = "}";
//		std::string sq_open = "[";
//		std::string sq_close = "]";
//		std::string equ = "=";
//		std::string comma = ",";
//		std::string semi = ";";
//		std::string open_comment = "/*";
//		std::string close_comment = "*/";
//		std::string space = " ";
//		std::string tab = "\t";
//		std::string tofs = "tofs";
//		std::string moms = "moms";
//
//		size_t found_tofs;
//		size_t found_moms;
//		size_t found_comment;
//		size_t found_semi;
//		size_t found_cur_open;
//		size_t found_cur_close;
//		size_t found_sq_open;
//		size_t found_sq_close;
////		size_t found_space;
//		size_t found_tab;
//		size_t found_Open_comment;
//		size_t found_Close_comment;
//		size_t found_comma;
//		size_t found_equ;
//
//		bool multi_line_comment = false;
//		bool inside_tofs = false;
//		bool inside_moms = false;
//
//		if (array_file.is_open())
//		{
//			while (!array_file.eof())
//			{
//				getline(array_file, line);
//
//				if (multi_line_comment) {
//					found_Close_comment = line.find(close_comment);
//					if (found_Close_comment != std::string::npos) {
//						multi_line_comment = false;
//						line = line.substr(found_Close_comment + 2);
//						//cout << "found close comment. multi_line_comment= " << multi_line_comment<< endl;
//					}
//				}
//				if (multi_line_comment == false) {
//					found_Open_comment = line.find(open_comment);
//					if (found_Open_comment != std::string::npos) {
//						line = line.substr(0, found_Open_comment);
//						multi_line_comment = true;
//						//cout << "found open comment. multi_line_comment= " << multi_line_comment<< endl;
//					}
//
//					// remove the comments
//					found_comment = line.find(comment);
//					if (found_comment != std::string::npos)
//						line = line.substr(0, found_comment);
//
//					// replace the tabs with spaces
//					found_tab = line.find(tab);
//					while (found_tab != std::string::npos) {
//						line = line.replace(found_tab, 1, " ");
//						found_tab = line.find(tab);
//					}
//
//					// replace the { with spaces
//					found_cur_open = line.find(cur_open);
//					while (found_cur_open != std::string::npos) {
//						line = line.replace(found_cur_open, 1, " ");
//						found_cur_open = line.find(cur_open);
//					}
//
//					// replace the } with spaces
//					found_cur_close = line.find(cur_close);
//					while (found_cur_close != std::string::npos) {
//						line = line.replace(found_cur_close, 1, " ");
//						found_cur_close = line.find(cur_close);
//					}
//
//					// replace the [ with spaces
//					found_sq_open = line.find(sq_open);
//					while (found_sq_open != std::string::npos) {
//						line = line.replace(found_sq_open, 1, " ");
//						found_sq_open = line.find(sq_open);
//					}
//
//					// replace the ] with spaces
//					found_sq_close = line.find(sq_close);
//					while (found_sq_close != std::string::npos) {
//						line = line.replace(found_sq_close, 1, " ");
//						found_sq_close = line.find(sq_close);
//					}
//
//					// replace the comma with spaces
//					found_comma = line.find(comma);
//					while (found_comma != std::string::npos) {
//						line = line.replace(found_comma, 1, " ");
//						found_comma = line.find(comma);
//					}
//
//					// replace the = with spaces
//					found_equ = line.find(equ);
//					while (found_equ != std::string::npos) {
//						line = line.replace(found_equ, 1, " ");
//						found_equ = line.find(equ);
//					}
//
//					found_tofs = line.find(tofs);
//					if (found_tofs != std::string::npos) {
//						inside_tofs = true;
//						line = line.erase(found_tofs, 4);
//					}
//
//					if (inside_tofs) {
//
//						found_semi = line.find(semi);
//						if (found_semi != std::string::npos) {
//							line = line.replace(found_semi, 1, " ");
//							inside_tofs = false;
//						}
//
//						std::stringstream line_stream(line);
//						while (!line_stream.eof()) {
//							double temp;
//							line_stream >> temp;
//							results[0].push_back(temp);
//						}
//					}
//
//
//					found_moms = line.find(moms);
//					if (found_moms != std::string::npos) {
//						inside_moms = true;
//						line = line.erase(found_moms, 4);
//					}
//
//					if (inside_moms) {
//
//						found_semi = line.find(semi);
//						if (found_semi != std::string::npos) {
//							line = line.replace(found_semi, 1, " ");
//							inside_moms = false;
//						}
//
//						std::stringstream line_stream(line);
//						while (!line_stream.eof()) {
//							double temp;
//							line_stream >> temp;
//							results[1].push_back(temp);
//						}
//					}
//
//
//				}
//
//			}
//		}
//
////		std::cout << "\nTofs = ";
////		for (int i = 0; i < results[0].size(); i++) {
////			std::cout << results[0][i] << "  ";
////		}
////
////		std::cout << "\nmoms = ";
////		for (int i = 0; i < results[1].size(); i++) {
////			std::cout << results[1][i] << "  ";
////		}
//
//		if (results[0].size() != results[1].size())
//			printf("\nERROR: interpolation arraies must be the same size!\n");
//
//		if (results[0][0] - results[0][1] > 0) { //check to see if it assending or decending
//			std::vector<std::vector<double>> flip(2, std::vector<double>());
//			int end = results[0].size()-1;
//			for (int i = 0; i <= end; i++) {
//				flip[0].push_back(results[0][end - i]);
//				flip[1].push_back(results[1][end - i]);
//			}
//			return flip;
//		}
//
//		return results;
//	}



	double NEW_v(double acc1, double acc2, double acc3, double acc1_length, double acc2_length, double acc3_length, double v, double t){

				double t_function = 0.;
				double Dt_function_dv = 0.;
				
				
				double a2_term = 0.0;
				double a3_term = 0.0;
				double v1 = 0.0;
				
					// Newton's method x1 = x0 - f(x0)/f'(x0)
					/*
					t_function =	
									-t + (-v + sqrt(2.0*acc1*acc1_length + v*v))/acc1 + 
									( (-sqrt(2.0*acc1*acc1_length + v*v) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v))/acc2 )+ 
									( (-sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + 2.0*acc3*acc3_length + v*v))/acc3 );

					Dt_function_dv = 
									(-1.0 + v/sqrt(2.0*acc1*acc1_length + v*v))/acc1 + 
									( (v*(-(1.0/sqrt(2.0*acc1*acc1_length + v*v)) + 1.0/sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v)))/acc2 )+ 
									( (v*(-(1.0/sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v)) + 1.0/sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + 2.0*acc3*acc3_length + v*v)))/acc3 );		
					*/
					//I must do it in parts because otherwise it will blow up.
					
					if (acc2 != 0. &&  acc2_length != 0.){
						a2_term = (-sqrt(2.0*acc1*acc1_length + v*v) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v))/acc2;
					}
					else if (acc2_length != 0. && acc2 == 0.){
						a2_term = acc2_length/sqrt(2.*acc1*acc1_length + v*v); 
					}
					else{
						a2_term = 0.0;
					}
					if (acc3 != 0. &&  acc3_length != 0.){
						a3_term = (-sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v*v) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + 2.0*acc3*acc3_length + v*v))/acc3;
					}
					else if (acc3_length != 0. && acc3 == 0.) {
						a3_term = acc3_length/sqrt(2.*acc1*acc1_length + 2.*acc2*acc2_length + v*v);
					}
					else{
						a3_term = 0.0;
					}
					
					//printf("a2_term = %5.5e | a3_term = %5.5e \n", a2_term, a3_term);

					t_function =	-t + (-v + sqrt(2.0*acc1*acc1_length + v*v))/acc1 + a2_term + a3_term; 
					
					
					v1 = v*(1.01);
					
					
					if (acc2 != 0. &&  acc2_length != 0.){
						a2_term = (-sqrt(2.0*acc1*acc1_length + v1*v1) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v1*v1))/acc2;
					}
					else if (acc2_length != 0. && acc2 == 0.){
						a2_term =acc2_length/sqrt(2.*acc1*acc1_length + v1*v1); 
					}
					else{
						a2_term = 0.0;
					}
					
					if (acc3 != 0. &&  acc3_length != 0.){
						a3_term = (-sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + v1*v1) + sqrt(2.0*acc1*acc1_length + 2.0*acc2*acc2_length + 2.0*acc3*acc3_length + v1*v1))/acc3;
					}
					else if (acc3_length != 0. && acc3 == 0.) {
						a3_term = acc3_length/sqrt(2.*acc1*acc1_length + 2.*acc2*acc2_length + v1*v1);
					}
					else{
						a3_term = 0.0;
					}
					
					
					Dt_function_dv =( (-t + (-v1 + sqrt(2.0*acc1*acc1_length + v1*v1))/acc1 + a2_term + a3_term) - t_function )/ (v1 - v); 

					
					
		
					return  v -  t_function / Dt_function_dv;;


	}

	//// **************************************************************
	//// **
	//// ** momentum calculation for time of flight direction for 3 
	//// ** accel. regions of any size.  This function uses newton's
	//// ** method to solve the problem.  Returns P in Au.  
	//// ** 
	//// ** By Joshua Williams
	//// **************************************************************
	double tof2mom_3accel(double tof_ns, double acc1mm, double acc2mm, double acc3mm, 
						  double Efield1_Vpcm, double Efield2_Vpcm, double Efield3_Vpcm, double charge_au, double mass_amu){

			
			//don't try to find P if the time of flight is negative
			if (tof_ns > 0.1){ 

				
				double E_Field_au_over_si =  1. / 5.1421e11; // 									au * m / V
				double length_au_over_m =  1. / 5.2917720859e-11;  // 						au / m
				double time_au_over_s = 1. / ( 2.418884326505e-17 ); // 						au / s
				double AUmass_over_AMU = 1.660538782e-27 / 9.10938215e-31 ;//		(au / kg) * (kg / amu)		
				double MperS_over_auVelocity = 2.1876912633e6;//								m / ( s * au )
				double Kg_over_AMU = 1.660538782e-27;			
				//double Psi_to_Pau = 1. / 1.992851565e-24;
				
				
				//convert to AU units
				double t = (tof_ns * 1e-9) * time_au_over_s  ;
				double mass = mass_amu *  AUmass_over_AMU;
				double q = charge_au;
				double acc1_length = acc1mm * 0.001 * length_au_over_m;
				double acc2_length = acc2mm * 0.001 * length_au_over_m;
				double acc3_length = acc3mm * 0.001 * length_au_over_m;
				double Efield1_Vau = Efield1_Vpcm * 100.0 * E_Field_au_over_si;
				double Efield2_Vau = Efield2_Vpcm * 100.0 * E_Field_au_over_si;
				double Efield3_Vau = Efield3_Vpcm * 100.0 * E_Field_au_over_si;
				double v=1;


				//find accel. for each region
				double acc1 = (Efield1_Vau * q) / mass;
				double acc2 = (Efield2_Vau * q) / mass;
				double acc3 = (Efield3_Vau * q) / mass;
				
				//find the an initial value for v (this is done through the a first order seires expansion of the t(v) about zero)
				double acc1_temp = 0.001, acc2_temp = 0.001, acc3_temp = 0.001;
				if(acc1 != 0.) acc1_temp = acc1;
				if(acc2 != 0.) acc2_temp = acc2;
				if(acc3 != 0.) acc3_temp = acc3;
				//YOU SHOULD CHANGE V IF YOU KNOW A BETTER GUESS 
				v = sqrt(2.)*sqrt(acc1_temp*acc1_length) - (sqrt(2.)*acc1_temp*sqrt(acc1_temp*acc1_length))/acc2_temp + (sqrt(2.)*acc1_temp*sqrt(acc1_temp*acc1_length + acc2_temp*acc2_length))/acc2_temp - (sqrt(2.)*acc1_temp*sqrt(acc1_temp*acc1_length + acc2_temp*acc2_length))/acc3_temp + (sqrt(2.)*acc1_temp*sqrt(acc1_temp*acc1_length + acc2_temp*acc2_length + acc3_temp*acc3_length))/acc3_temp - acc1_temp*t;
				
				if ( _finite(v) == 0.0 ) {
					//**********
					v=5;		//if you get here this proabably means you are using a decceration field, so I will guess a large value
					//*********
				}

				double new_v = 0.;
				
				
				//int limit = 10;
				int shift_counter = 1;
				int loop_counter = 0;
				for( int i=0;i<50;i++){
					loop_counter++;
					new_v = NEW_v( acc1,  acc2,  acc3,  acc1_length,  acc2_length,  acc3_length,  v, t);
					
					if ( _finite(new_v) == 0.0 ) {
						new_v = v * (1. +  pow(2., (double)shift_counter) / 1000.);   //move it, at first just a little bit but increase the size of the shift until you get a finite value for the new_v
						shift_counter++;
					}
					else{
						shift_counter = 1;  // reset the shift counter when you get a finite value for new_v
					}
					
					if ( fabs((new_v * mass  - v * mass )/ (v * mass )) < 1e-6) {
						i = 500; //thus exiting the loop
					}  
					v = new_v;
				}
				
				if ( loop_counter >= 50) { return (1e50 ); } // the function did not converge (but the value might still be pretty to the correct value) 
				return v * mass;
							  
				}

			else{
				return -100000.; 
				//if negative tof then the data point is crap so return crap
			}

	}

	// **************************************************************
	// **
	// ** calculate momenta for electron x,y direction (magnetic field)
	// ** using Mirko's functions.
	// **
	// **************************************************************

	double calc_px(double tof_ns, double x_mm, double y_mm, double mass_amu, double charge_au, double BField_ns, bool BField_clockwise = true) {

		double px;

		if(mass_amu<1.0) {	
			double w,a,b;

			double m = mass_amu * MASSAU * MEKG;
			double q = charge_au * COULOMB;
			double pau = m*300.e6/137.;

			double fieldB = 2.*m*3.14152 / (q * BField_ns * 1e-9);

			if(BField_clockwise)
				fieldB = -fieldB;

			w = q / m * fieldB;
			a = (1. - cos(w * tof_ns*1.e-9)) / w;
			b = (sin(w * tof_ns*1.e-9)) / w;

			px = m * (x_mm/1000. * b - a*y_mm/1000.) / (a*a + b*b);
			px = px / pau;
		} else {
			double vau = 2.1877e+6;
			px = x_mm/1000. / ((tof_ns)*1e-9) / vau * mass_amu * MASSAU;
		}
		return px;
	}

	double calc_py(double tof_ns, double x_mm, double y_mm, double mass_amu, double charge_au, double BField_ns, bool BField_clockwise = true) {

		double py;

		if(mass_amu<1.0) {	
			
			double w,a,b;

			double m = mass_amu * MASSAU * MEKG;
			double q = charge_au * COULOMB;
			double pau = m*300.e6/137.;

			double fieldB = 2.*m*3.14152 / (q * BField_ns * 1e-9);

			if(BField_clockwise)
				fieldB = -fieldB;

			w = q / m * fieldB;
			a = (1. - cos(w * tof_ns*1.e-9)) / w;
			b = (sin(w * tof_ns*1.e-9)) / w;

			py = m * (-x_mm/1000. * a - b*y_mm/1000.) / (a*a + b*b);
			py = py / pau;

		} else {
			double vau = 2.1877e+6;
			py = y_mm/1000. / ((tof_ns)*1e-9) / vau * mass_amu * MASSAU;
		}
		
		return py;
	}

	// **************************************************************
	// **
	// ** This function can calculate the time of flight for the
	// ** second recoil.  This function uses tof2mom_3accel
	// ** and assumes that it is a two body break up. 
	// **
	// ** 
	// ** By Joshua Williams
	// **************************************************************
	double t2_3accel(double tof_ns, double acc1mm, double acc2mm, double acc3mm, 
						  double Efield1_Vpcm, double Efield2_Vpcm, double Efield3_Vpcm, double charge1_au, double charge2_au, double mass1_amu, double mass2_amu){
			double t = tof_ns * 1e-9;
			double mass = mass2_amu * 1.660538782e-27;
			double q = charge2_au * 1.602176487e-19;;
			double acc1m = acc1mm * .001;
			double acc2m = acc2mm * .001;
			double acc3m = acc3mm * .001;
			double Efield1_Vpm = Efield1_Vpcm * 100.0;
			double Efield2_Vpm = Efield2_Vpcm * 100.0;
			double Efield3_Vpm = Efield3_Vpcm * 100.0;

			//find accel. for each region
			double acc1 = (Efield1_Vpm * q)/ mass;
			double acc2 = (Efield2_Vpm * q)/ mass;
			double acc3 = (Efield3_Vpm * q)/ mass;
			
			double  acc2_zero_bool = 0.;
			double  acc3_zero_bool = 0.;
				

				
			double a2_term = 0.0;
			double a3_term = 0.0;


			double r1_px = - tof2mom_3accel(tof_ns, acc1mm, acc2mm, acc3mm, 
							  Efield1_Vpcm, Efield2_Vpcm, Efield3_Vpcm, charge1_au, mass1_amu);
			
			double v = r1_px / (  mass / (9.1093826e-31*2.1876912633e6) );
			

			if (acc2 != 0. &&  acc2m != 0.){
				a2_term = (-sqrt(2.0*acc1*acc1m + v*v) + sqrt(2.0*acc1*acc1m + 2.0*acc2*acc2m + v*v))/acc2;
			}
			else if (acc2m != 0. && acc2 == 0.){
				a2_term = acc2m/sqrt(2.*acc1*acc1m + v*v); 
			}
			else{
				a2_term = 0.0;
			}
			
			if (acc3 != 0. &&  acc3m != 0.){
				a3_term = (-sqrt(2.0*acc1*acc1m + 2.0*acc2*acc2m + v*v) + sqrt(2.0*acc1*acc1m + 2.0*acc2*acc2m + 2.0*acc3*acc3m + v*v))/acc3;
			}
			else if (acc3m != 0. && acc3 == 0.) {
				a3_term = acc3m/sqrt(2.*acc1*acc1m + 2.*acc2*acc2m + v*v);
			}
			else{
				a3_term = 0.0;
			}
			
			/*if (counter < 10){
				printf("a2_term = %e, a3_term = %e, " , a2_term, a3_term);
			}*/	
			double t2 = (-v + sqrt(2.0*acc1*acc1m + v*v))/acc1 + a2_term + a3_term; 
			
			
			
			return (t2 * 1.e+9);
		
	}
}
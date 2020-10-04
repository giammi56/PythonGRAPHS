#pragma once

#ifndef FUNCTIONS_ALREADY_INCLUDED
	#define FUNCTIONS_ALREADY_INCLUDED

	#include "OS_Version.h"

	#include "TF1.h"		//Needed for Tof_to_P
	#include <vector>
	#include <string>



//std::vector<std::vector<double>> * read_array_file(std::string filename);


	//Takes at set of tof and momentum points as imput and returns a TF1 fuction that is fitted to these points.
	//This TF1 can then be used to convert tof to momentum.  The funciton displays the results of the fit in a window while LMF2Root is running.
	//It uses a polynomial function A*x^4+B*x^3+C*x^2+D*x+C, where A,B,C,D,E are fitted parameters.  
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//Example (this code is placed in Analysis.cpp):
	//		static TF1 * H_plus;
	//		if (eventcounter == 0) { //this code runs only once!!
	//		
	//			double tofs[] = { 1602.28,	1600.11,	1597.91,	1595.69,	1181.80,	1179.73 }; 
	//			double moms[] = { -38.51,	-38.12,		-37.73,		-37.33,		38.12,		38.51	};
	//		
	//			H_plus = Tof_to_P("H_plus_", 5, tofs, moms);
	//		}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// to get the value from the fitted function: H_plus->Eval(CH_evt->r.tof[0])
	//
	// By Joshua Williams
	TF1 *  Tof_to_P(char * name, int num_data_points, double *tofs, double *moms);


	// **************************************************************
	// **
	// **  momentum calculation from t_0 (two field regions)
	// **
	// **************************************************************
	double tof2mom21(double tof_ns, double acc_mm, double drift_mm, double Efeld_Vpcm, double mass_au, double charge_au, double  parameter[]);


	// **************************************************************
	// **
	// **  momentum calculation from t_0 (homogeneous acceleration)
	// **
	// **************************************************************
	double tof2mom(double tof_ns, double acc_mm, double Efield_Vpcm, double mass_au, double charge_au, double parameter[]);


	// **************************************************************
	// **
	// ** calculate momenta for tof-direction /w McLaren geometry
	// **
	// **************************************************************
	double toftomom(double tof_ns, double acc_mm, double drift_mm, double Efeld_Vpcm, double mass_au, double charge_au, double  parameter[]);
	

	// **************************************************************
	// **
	// ** position correction for electron (E-B-drift)
	// **
	// **************************************************************
	double elec_pos_corr(double e1x_mm, double e1y_mm, double e1t_ns, char direction, double  parameter[]);	


	// **************************************************************
	// **
	// ** calculate momenta for electron (magnetic field)
	// **
	// **************************************************************
	double elec2mom(double e1x_mm, double e1y_mm, double e1t_ns, char direction, double  parameter[]);


	// **************************************************************
	// **
	// ** find channels in PIPICO: t1 [ns], m1,m2 [au], q1,q2 [e], s[mm], fieldE [V/cm] 
	// **
	// **************************************************************
	double t2(double t1,double m1,double m2, double q1, double q2, double s, double fieldE, double  parameter[]); 


	// **************************************************************
	// **
	// ** find channels in PIPIPICO: t1,t2 [ns], m1,m2,m2 [u], q1,q2,q3 [e], s[mm], fieldE [V/cm] 
	// **
	// **************************************************************
	double t3(double t1, double t2, double m1,double m2, double m3, double q1, double q2, double q3, double s, double fieldE);

	// **************************************************************
	// **
	// ** calculate momenta for electron x,y direction (magnetic field)
	// ** using Mirko Hattass' functions.
	// **
	// **************************************************************

	double electron_px(double etof, double xe, double ye, double  parameter[]);
	double electron_py(double etof, double xe, double ye, double  parameter[]);


	// ********************************************************************
	// **
	// ** calculate angle of Delta_Phi between Recoil and projectile
	// **
	// ********************************************************************
	double deltaphi(double phi_projectile, double phi_molecule);


	// **************************************************************
	// **
	// ** Labframe transformation - 3 vectors are transformed into a new
	// ** frame of reference
	// **
	// **************************************************************
	void labframe_transformation(__int32 direction, double a[], double b[], double c[]);

#endif
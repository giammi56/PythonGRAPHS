#pragma warning(disable : 4996)
#include "../src/CH_Presorter.h"
#include <math.h>
#include "../src/CH_Tof.h"
#include "../src/CH_Spectrometer.h"
#include "../src/CH_FUN_Lowlevel.h"

namespace CH 
{

	// User defined presorter functions:
	// ------------------------------------------------------------------------------
	// Return 'false' if the whole event is not valid.
	// Use the 'hit_valid_arrays' to check which hits are already discarded by other presorters,
	// and to set which hits are still valid after your presorter was used:
	//
	// bool e_hit_valid_array[64];
	// bool r_hit_valid_array[64];
	// bool p_hit_valid_array[64];
	//
	// -You can access the event as 'evti.', e.g. to access the number of hits of the electrons: 'evti->e.num_hits'.
	// -The parameters set in the XML-File can be accessed as 'par[x]'.
	// -If you want to plot a histogram, use: H->fill1(HISTNUM + this->HistoBlock*NHISTS,...,dir); 
	//  (use 50< HISTNM < 100. If you need more histograms ask Till..
	//
	// Enjoy!
	// ------------------------------------------------------------------------------

	bool presorter_class::user_prs_0() {
		double *par;
		par = this->usr_par[0];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_1() {
		double *par;
		par = this->usr_par[1];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_2() {
		double *par;
		par = this->usr_par[2];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_3() {
		double *par;
		par = this->usr_par[3];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_4() {
		double *par;
		par = this->usr_par[4];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_5() {
		double *par;
		par = this->usr_par[5];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_6() {
		double *par;
		par = this->usr_par[6];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

	bool presorter_class::user_prs_7() {
		double *par;
		par = this->usr_par[7];
		// add your code here
		//--------------------------------------------
	

		return true;
	}

}
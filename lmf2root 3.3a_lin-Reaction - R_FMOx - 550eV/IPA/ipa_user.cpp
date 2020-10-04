// ------------------------------------------------------------------------------------------
// Template for tuning of the sum momentum of a diatomic.
//
// Simply copy functions to "user.cpp" replacing the old ones...
// ------------------------------------------------------------------------------------------

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TProfile.h"
#include "TCutG.h"
#include "TPolyLine.h"
#include "TObjArray.h"
#include "TF1.h"
#include "ipa.h"
void ipa_class::init_user_histograms() {
	// Due to the multithreading the histograms need to be defined as a double array 
	// with first index being "[th]".
	// Example: "this->H1D[th][0]" is the first 1D histogram.

	double prel = 150.0;
	int prel_bins = 200;

	double psum = 50.0;
	int psum_bins = 50;

	this->myCanvas[0] = root->newCanvas("IPA","Interactive Parameter Adjustment - The Canvas",10,300,330*3,330*2);
	this->myCanvas[0]->Divide(3,2);
	char s[20];
	for(int th=0;th<MAXTHREADS;th++) {
		sprintf_s(s,"px1x2_%i",th);
		this->H2D[th][0] = root->newTH2D(s,"",prel_bins,-prel,prel,"p_rel_x1",prel_bins,-prel,prel,"p_rel_x2","COLZ");
		sprintf_s(s,"py1y2_%i",th);
		this->H2D[th][1] = root->newTH2D(s,"",prel_bins,-prel,prel,"p_rel_y1",prel_bins,-prel,prel,"p_rel_y2","COLZ");
		sprintf_s(s,"pz1z2_%i",th);
		this->H2D[th][2] = root->newTH2D(s,"",prel_bins,-prel,prel,"p_rel_z1",prel_bins,-prel,prel,"p_rel_z2","COLZ");
		sprintf_s(s,"psumxy_%i",th);
		this->H2D[th][3] = root->newTH2D(s,"",psum_bins,-psum,psum,"p_sum_x",psum_bins,-psum,psum,"p_sum_y","COLZ");
		sprintf_s(s,"psumxz_%i",th);
		this->H2D[th][4] = root->newTH2D(s,"",psum_bins,-psum,psum,"p_sum_x",psum_bins,-psum,psum,"p_sum_z","COLZ");
		sprintf_s(s,"psumyz_%i",th);
		this->H2D[th][5] = root->newTH2D(s,"",psum_bins,-psum,psum,"p_sum_y",psum_bins,-psum,psum,"p_sum_z","COLZ");
	}
}

void ipa_class::user_function(CH_event_struct *evti, int th) {
	// Due to the multithreading you need to add "[th]" to e, r, mol and big_mol!
	// Example: second electron is "e[th][1]->", a diatomic is "mol[th]->" etc.

	H2D[th][0]->Fill(mol[th]->ion[0]->mom.x,mol[th]->ion[1]->mom.x,1.0);
	H2D[th][1]->Fill(mol[th]->ion[0]->mom.y,mol[th]->ion[1]->mom.y,1.0);
	H2D[th][2]->Fill(mol[th]->ion[0]->mom.z,mol[th]->ion[1]->mom.z,1.0);

	H2D[th][3]->Fill(mol[th]->mom_cm.x,mol[th]->mom_cm.y,1.0);
	H2D[th][4]->Fill(mol[th]->mom_cm.x,mol[th]->mom_cm.z,1.0);
	H2D[th][5]->Fill(mol[th]->mom_cm.y,mol[th]->mom_cm.z,1.0);
}

void ipa_class::user_draw() {
	// Plotting the histos: first array index is always "[0]", second is the one of 
	// the histogram you want to plot. E.g. H2D[0][2]->Draw() will plot your histogram
	// number two.

	myCanvas[0]->cd(1);
	H2D[0][0]->Draw();
	myCanvas[0]->cd(2);
	H2D[0][1]->Draw();
	myCanvas[0]->cd(3);
	H2D[0][2]->Draw();

	myCanvas[0]->cd(4);
	H2D[0][3]->Draw();	
	myCanvas[0]->cd(5);
	H2D[0][4]->Draw();
	myCanvas[0]->cd(6);
	H2D[0][5]->Draw();

	myCanvas[0]->Modified();
	myCanvas[0]->Update();
}
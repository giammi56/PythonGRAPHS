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
#include "TLine.h"
#include "TObjArray.h"
#include "TF1.h"
#include "ipa.h"

void ipa_class::enable_wiggletool(int bins, double min, double max) {
	this->NumPars[PPCOPAGE]= 1;
	this->wiggletool_enabled = true;
	this->wiggletool_min=min;
	this->wiggletool_max=max;

	this->H_wiggletool_x = root->newTH2D("fish_x","The Fish, x-direction",bins,min,max,"time-of-flight [ns]",100,-60.0,60.0,"x [mm]","COLZ");
	this->H_wiggletool_y = root->newTH2D("fish_y","The Fish, y-direction",bins,min,max,"time-of-flight [ns]",100,-60.0,60.0,"y [mm]","COLZ");
	this->H_wiggletool_r = root->newTH2D("wiggle","The Wiggle",bins,min,max,"time-of-flight [ns]",100,0.0,60.0,"r [mm]","COLZ");

}

void ipa_class::init_wiggletool() {

	this->myCanvas2 = root->newCanvas("IPA_PPCO","IPA - the wiggle tool (press space to update)",40,20,330*2,330*3);	
	this->myCanvas2->SetFixedAspectRatio(false);
	this->myCanvas2->Divide(1,3);

	this->wiggletool_t_0_elec = IPA_Proc[0]->e_fac->dt;
	this->wiggletool_t_0_tof_calc = IPA_Tof->get_t_0_shift();

	add_par(this->wiggletool_t_0_tof_calc,-1000.0,1000.0,1.0,"Overall t0 shift",70);
}

void ipa_class::Fill_wiggletool(CH_event_struct *evnt) 
{
	if(evnt->e.num_hits>0) {

		e_raw[0]->method = evnt->e.method[0];
		e_raw[0]->data.x = evnt->e.x[0];
		e_raw[0]->data.y = evnt->e.y[0];
		e_raw[0]->data.time = evnt->e.time[0];
		e_raw[0]->data.tof = evnt->e.tof[0];
	
		e[0][0]->set_raw(*e_raw[0]);
		e[0][0]->shift_stretch_raw(IPA_Proc[0]->e_fac);
		e[0][0]->rotate_raw(IPA_Proc[0]->e_fac->rot_ang);
		e[0][0]->calc_phi_pos();		

		this->H_wiggletool_x->Fill(e[0][0]->raw.data.tof,e[0][0]->raw.data.x);
		this->H_wiggletool_y->Fill(e[0][0]->raw.data.tof,e[0][0]->raw.data.y);
		this->H_wiggletool_r->Fill(e[0][0]->raw.data.tof,sqrt(e[0][0]->raw.data.x*e[0][0]->raw.data.x+e[0][0]->raw.data.y*e[0][0]->raw.data.y));
		num_events_read++;
	}
}

void ipa_class::draw_wiggletool() {
	this->myCanvas2->cd(1);
	this->H_wiggletool_x->Draw("COLZ");
	draw_wiggletool_markers();

	this->myCanvas2->cd(2);
	this->H_wiggletool_y->Draw("COLZ");
	draw_wiggletool_markers();

	this->myCanvas2->cd(3);
	this->H_wiggletool_r->Draw("COLZ");
	draw_wiggletool_markers(true);

	this->myCanvas2->Modified();
	this->myCanvas2->Update();
}

void ipa_class::draw_wiggletool_markers(bool r) {
	double B_ns= this->IPA_Proc[0]->Spect->Bfield_ns; 
	B_ns-= (IPA_Proc[0]->e_fac->dt - this->wiggletool_t_0_elec);
	B_ns+= (this->wiggletool_t_0_tof_calc - par[70].value);

	while(wiggletool_min<B_ns && B_ns<wiggletool_max) {
		double ymin=-60.0;
		if(r)
			ymin=-0.0;
		
		TLine *line = new TLine(B_ns,ymin,B_ns,60.0);
		line->SetLineColor(kRed);
	    line->SetLineWidth(1);
		line->Draw();

		B_ns += this->IPA_Proc[0]->Spect->Bfield_ns;
	}
}
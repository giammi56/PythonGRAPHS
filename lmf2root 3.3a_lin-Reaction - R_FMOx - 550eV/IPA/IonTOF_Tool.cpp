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
#include "TText.h"
#include "TObjArray.h"
#include "TF1.h"
#include "ipa.h"

void ipa_class::enable_iontoftool(int bins, double min, double max) {
	this->NumPars[PPCOPAGE]= 9;
	this->iontoftool_enabled = true;
	this->iontoftool_min=min;
	this->iontoftool_max=max;

	this->H_iontoftool_x = root->newTH2D("fish_x","",bins,min,max,"time-of-flight [ns]",100,-60.0,60.0,"x [mm]","COLZ");
	this->H_iontoftool_y = root->newTH2D("fish_y","",bins,min,max,"time-of-flight [ns]",100,-60.0,60.0,"y [mm]","COLZ");
}

void ipa_class::init_iontoftool() {
	this->myCanvas2 = root->newCanvas("IPA_PPCO","IPA - the ion tof tool (press space to update)",40,20,330*3,330*2);	
	this->myCanvas2->SetFixedAspectRatio(false);
	this->myCanvas2->Divide(1,2);

	// add ions
	add_par(41,1.0,400.0,1.0,"m/q ion 1",70);
	add_par(42,1.0,400.0,1.0,"m/q ion 2",71);
	add_par(43,1.0,400.0,1.0,"m/q ion 3",72);
	add_par(44,1.0,400.0,1.0,"m/q ion 4",73);
	add_par(68,1.0,400.0,1.0,"m/q ion 5",74);
	add_par(69,1.0,400.0,1.0,"m/q ion 6",75);
	add_par(70,1.0,400.0,1.0,"m/q ion 7",76);
	add_par(14,1.0,400.0,1.0,"m/q ref ion",77);
	add_par(2999.0,1.0,100000.0,100.0,"tof ref ion",78);

	// store intial position offsets as they are already included in the histo
	this->dx = IPA_Proc[0]->r_fac->dx;
	this->dy = IPA_Proc[0]->r_fac->dy;
	this->dt = IPA_Proc[0]->r_fac->dt;
	this->drot = IPA_Proc[0]->r_fac->rot_ang; 
}

void ipa_class::Fill_iontoftool(CH_event_struct *evnt) {
	for(int i=0;i<evnt->r.num_hits;i++) {
		r_raw[0]->method = evnt->r.method[i];
		r_raw[0]->data.x = evnt->r.x[i];
		r_raw[0]->data.y = evnt->r.y[i];
		r_raw[0]->data.time = evnt->r.time[i];
		r_raw[0]->data.tof = evnt->r.tof[i];
	
		r[0][0]->set_raw(r_raw[0]);
		r[0][0]->shift_stretch_raw(IPA_Proc[0]->r_fac);
		r[0][0]->rotate_raw(IPA_Proc[0]->r_fac->rot_ang);
		r[0][0]->calc_phi_pos();		

		this->H_iontoftool_x->Fill(r[0][0]->raw.data.tof,r[0][0]->raw.data.x);
		this->H_iontoftool_y->Fill(r[0][0]->raw.data.tof,r[0][0]->raw.data.y);
	}
	num_events_read++;

}

void ipa_class::draw_iontoftool() {
	this->myCanvas2->cd(1);
	this->H_iontoftool_x->Draw("COLZ");
	draw_iontoftool_markers();
	draw_iontoftool_jet_markers(false);

	this->myCanvas2->cd(2);
	this->H_iontoftool_y->Draw("COLZ");
	draw_iontoftool_markers();
	draw_iontoftool_jet_markers(true);

	this->myCanvas2->Modified();
	this->myCanvas2->Update();
}

void ipa_class::draw_iontoftool_jet_markers(bool y) {
	double y1, y0;
	double vel_mmns=IPA_Proc[0]->Spect->VJet/1e+9*1000.0;
	double rot = (IPA_Proc[0]->Spect->AngJet + (IPA_Proc[0]->r_fac->rot_ang-this->drot)) /180.0*PI;
	if(y) {
		y0 = (IPA_Proc[0]->r_fac->dy - this->dy) + vel_mmns * cos(rot)*(this->iontoftool_min);
		y1 = y0 + vel_mmns * cos(rot)*(this->iontoftool_max-this->iontoftool_min);
	} else {	
		y0 = (IPA_Proc[0]->r_fac->dx - this->dx) + vel_mmns * sin(rot)*(this->iontoftool_min);
		y1 = y0 + vel_mmns * sin(rot)*(this->iontoftool_max-this->iontoftool_min);
	}
	TLine *line = new TLine(this->iontoftool_min,y0,this->iontoftool_max,y1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->SetLineStyle(2);
	line->Draw();
}

void ipa_class::draw_iontoftool_markers() {
	double tofs[16];

	if(this->IPA_Proc[0]->Spect->ion_side->linear_approximation) {
		// tof are realtive to reference ion
		double reftof = par[78].value;
		double refmq = par[77].value;

		TLine *line = new TLine(reftof,-60.0,reftof,60.0);
		line->SetLineColor(kBlack);
		line->SetLineWidth(2);
		line->Draw();

		for(int i=0;i<7;i++) {
			tofs[i] = reftof/sqrt(refmq/par[70+i].value); 
		}
	} else {
		// get spectrometer voltages
		double U[4];
		U[0]=0.0;
		for(int j=0;j<IPA_Proc[0]->Spect->ion_side->number_of_regions;j++) {
			U[j+1]=U[j]+IPA_Proc[0]->Spect->ion_side->Efields[j]*IPA_Proc[0]->Spect->ion_side->lengths[j]/10.0;
		}
		// get tofs from spectrometer geometry
		for(int i=0;i<7;i++) {
			tofs[i] = 0.0; 
			double V0 = 0.;
			double V1 = 0.;
			double mq = par[70+i].value*MASSAU; 

			for(int j=0;j<IPA_Proc[0]->Spect->ion_side->number_of_regions;j++) {		
				V0=2.188e-1*sqrt((U[j])/mq/27.21*2.); 
				V1=2.188e-1*sqrt((U[j+1])/mq/27.21*2.); 
			
				tofs[i] = tofs[i]+(IPA_Proc[0]->Spect->ion_side->lengths[j]/10.0)*2./(V0+V1);
			}
		}

	}

	// add ion t_0
	for(int i=0;i<7;i++) {
		tofs[i]+= (IPA_Proc[0]->r_fac->dt - this->dt);
	}

	for(int i=0;i<7;i++) {
		char tag[16];
		if(tofs[i]<this->iontoftool_max) {
			TLine *line = new TLine(tofs[i],-60.0,tofs[i],60.0);
			line->SetLineColor(kRed);
			line->SetLineWidth(1);
			line->Draw();
			sprintf(tag,"m/q=%i",(int)par[70+i].value);
			TText *th1 = new TText(tofs[i],60,tag);
			th1->SetTextAlign(11); 
			th1->SetTextSize(0.02f);
			th1->SetTextAngle(45);
			th1->Draw();

		}
	}
}
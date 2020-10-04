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

void ipa_class::enable_PPCO(int bins, double min, double max, bool momcon) {
	this->NumPars[PPCOPAGE]= 12;
	this->PPCO_enabled = true;
	this->PPCO_momcon = momcon;
	this->pco_bins = bins;
	this->pco_min = min;
	this->pco_max = max;
	this->H_PIPICO = root->newTH2D("PIPICO","Photoion/photoion-coincidence plot",pco_bins,pco_min,pco_max,"time-of-flight ion 1 [ns]",pco_bins,pco_min,pco_max,"time-of-flight ion 1 [ns]","COLZ");
}

void ipa_class::init_PIPICO() {

	this->myCanvas2 = root->newCanvas("IPA_PPCO","IPA - the PIPICO tool press space to update)",40,320,330*2,330*2);	
	this->myCanvas2->SetFixedAspectRatio(true);

	add_par(14.0,1.0,1000.0,1.0,"ch 1:  mass 1 [amu]",70);
	add_par(14.0,1.0,1000.0,1.0,"ch 1:  mass 2 [amu]",71);
	add_par(1.0,1.0,1000.0,1.0,"ch 1: charge 1 [au]",72);
	add_par(1.0,1.0,1000.0,1.0,"ch 1: charge 2 [au]",73);

	add_par(16.0,1.0,1000.0,1.0,"ch 2:  mass 1 [amu]",74);
	add_par(16.0,1.0,1000.0,1.0,"ch 2:  mass 2 [amu]",75);
	add_par(1.0,1.0,1000.0,1.0,"ch 2: charge 1 [au]",76);
	add_par(1.0,1.0,1000.0,1.0,"ch 2: charge 2 [au]",77);

	add_par(17.0,1.0,1000.0,1.0,"ch 3:  mass 1 [amu]",78);
	add_par(1.0,1.0,1000.0,1.0,"ch 3:  mass 2 [amu]",79);
	add_par(1.0,1.0,1000.0,1.0,"ch 3: charge 1 [au]",80);
	add_par(1.0,1.0,1000.0,1.0,"ch 3: charge 2 [au]",81);
/*
	add_par(1.0,1.0,1000.0,1.0,"ch 4:  mass 1 [amu]",82);
	add_par(1.0,1.0,1000.0,1.0,"ch 4:  mass 2 [amu]",83);
	add_par(1.0,1.0,1000.0,1.0,"ch 4: charge 1 [au]",84);
	add_par(1.0,1.0,1000.0,1.0,"ch 4: charge 2 [au]",85);
*/
}

void ipa_class::FillPPCO(CH_event_struct *evnt) 
{
	if(evnt->r.num_hits>1) {
		for(int i=0;i<2;i++) {
			r_raw[i]->method = evnt->r.method[i];
			r_raw[i]->data.x = evnt->r.x[i];
			r_raw[i]->data.y = evnt->r.y[i];
			r_raw[i]->data.time = evnt->r.time[i];
			r_raw[i]->data.tof = evnt->r.tof[i];
		
			r[0][i]->set_raw(r_raw[i]);
			r[0][i]->shift_stretch_raw(IPA_Proc[0]->r_fac);
			r[0][i]->rotate_raw(IPA_Proc[0]->r_fac->rot_ang);
			r[0][i]->calc_phi_pos();		
		}
		if(this->PPCO_momcon) {
			double phi_p = r[0][0]->raw.phi + 180.;
			double phi_m = r[0][1]->raw.phi + 180.;
			double relphi;

			relphi = phi_p - phi_m;
				if (phi_p > phi_m) {
					if (relphi > 180.) {
						relphi = relphi - 360.;
					}
				} else {
					if (phi_p <= phi_m) {
						if (relphi < -180.) {
							relphi = relphi + 360.;
						}
					}
				}
			
			if((fabs(relphi - 180.0) < 10.0)) {
				this->H_PIPICO->Fill(r[0][0]->raw.data.tof,r[0][1]->raw.data.tof);
				num_events_read++;
			}
		} else {
			this->H_PIPICO->Fill(r[0][0]->raw.data.tof,r[0][1]->raw.data.tof);
			num_events_read++;
		}
	}
}

void ipa_class::draw_PIPICO() {
	this->myCanvas2->cd(0);
	this->H_PIPICO->Draw();

	for(int i=0;i<3;i++) {
		draw_PIPICOline(i+1, this->ppco_m1[i], this->ppco_m2[i], this->ppco_q1[i], this->ppco_q2[i]);
	}
	this->myCanvas2->Modified();
	this->myCanvas2->Update();
	
}

void ipa_class::draw_PIPICOline(int color, double MassP1_amu, double MassP2_amu, double ChargeP1_au, double ChargeP2_au) {

    double tof1[1000];
    double tof2[1000];
    double t1;

    double tof1_2[1000];
    double tof2_2[1000];

    int MinTof,MaxTof;
    int YMinTof,YMaxTof;

	double LenReg1_mm = this->IPA_Proc[0]->Spect->ion_side->lengths[0];
	double LenReg2_mm = 0.0;
	double LenReg3_mm = 0.0;
	if(this->IPA_Proc[0]->Spect->ion_side->number_of_regions>1)
		LenReg2_mm = this->IPA_Proc[0]->Spect->ion_side->lengths[1];
	if(this->IPA_Proc[0]->Spect->ion_side->number_of_regions>2)
		LenReg3_mm = this->IPA_Proc[0]->Spect->ion_side->lengths[2];

	double EField1_Vpcm = this->IPA_Proc[0]->Spect->ion_side->Efields[0];
	double EField2_Vpcm = 0.0;
	double EField3_Vpcm = 0.0;
	if(this->IPA_Proc[0]->Spect->ion_side->number_of_regions>1)
		EField2_Vpcm = this->IPA_Proc[0]->Spect->ion_side->Efields[1];
	if(this->IPA_Proc[0]->Spect->ion_side->number_of_regions>2)
		EField3_Vpcm = this->IPA_Proc[0]->Spect->ion_side->Efields[2];


// get histogram
	TH2D *hist1 = this->H_PIPICO;

	MinTof = (int)pco_min;
    MaxTof = (int)pco_max;

    YMinTof = (int)pco_min;
    YMaxTof = (int)pco_max;

    t1 = double(MinTof);
    for (int i=0;i<1000;++i)
    {

        tof1[i] = t1;
        tof2[i] = t2_3accel(t1+IPA_Proc[0]->r_fac->dt, LenReg1_mm, LenReg2_mm, LenReg3_mm, EField1_Vpcm, EField2_Vpcm, EField3_Vpcm,  ChargeP1_au,  ChargeP2_au, MassP1_amu, MassP2_amu);
       
        if(tof2[i]>YMaxTof)
            tof2[i]=YMaxTof;
        if(tof2[i]<YMinTof)
            tof2[i]=YMinTof;
               
        t1 = t1 + double((MaxTof-MinTof)/1000.);
    }

    TPolyLine *pline = new TPolyLine(1000,tof1,tof2);
    pline->SetLineColor(color);
    pline->SetLineWidth(1);
    pline->Draw();

	// Asym. breakup? We need a second line!
    if((ChargeP1_au!=ChargeP2_au) || (MassP1_amu!=MassP2_amu)) {

        t1 = double(MinTof);

		for (int i=0;i<1000;++i)
        {
            tof1_2[i] = t1;
            tof2_2[i] = t2_3accel(t1+IPA_Proc[0]->r_fac->dt, LenReg1_mm, LenReg2_mm, LenReg3_mm, EField1_Vpcm, EField2_Vpcm, EField3_Vpcm, ChargeP2_au, ChargeP1_au, MassP2_amu, MassP1_amu);
           
            if(tof2_2[i]>YMaxTof)
                tof2_2[i]=YMaxTof;
            if(tof2_2[i]<YMinTof)
                tof2_2[i]=YMinTof;
            t1 = t1 + double((MaxTof-MinTof)/1000.);
        }

        TPolyLine *pline = new TPolyLine(1000,tof1_2,tof2_2);
        pline->SetLineColor(color);
        pline->SetLineWidth(1);
        pline->Draw();
    }

}

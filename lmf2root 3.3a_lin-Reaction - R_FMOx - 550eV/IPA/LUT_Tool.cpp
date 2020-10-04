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

void ipa_class::enable_LUT(double mom, int bins_mom, int bins_phi, int bins_ctheta, int electrons) {
	this->NumPars[PPCOPAGE]= 7;
	this->LUT_enabled = true;
	this->LUT_electrons = electrons;
	this->LUT_mom = mom;
	this->LUT_bins_mom = bins_mom;
	this->LUT_bins_phi = bins_phi;
	this->LUT_bins_ctheta = bins_ctheta;
	this->LUT_Profile = 0;

	for(int i=0;i<bins_phi;i++) {
		char s[10];
		for(int j=0;j<MAXTHREADS;j++) {
			sprintf(s,"phi_%i_%i",i,j);//(double)i*(360.0/bins_phi));
			this->ctpeArray[j][i] = new TH2D(s,"",bins_ctheta,-1.0,1.0,bins_mom,mom*0.6,mom*1.4);
		}
		this->aSlices[i] = new TObjArray();
		this->Fit[i] = new TH1D();
	}
}

void ipa_class::write_correctiontable() {

	FILE *f = fopen("corr.txt","wt");

	fprintf(f,"int bins_phi = %d;\n",this->LUT_bins_phi);
	fprintf(f,"int bins_ctheta = %d;\n",this->LUT_bins_ctheta);

	fprintf(f,"double corr_tab[%d][%d] = { ",this->LUT_bins_phi,this->LUT_bins_ctheta);
	
	for(int phiang=0;phiang<this->LUT_bins_phi;phiang++) {

		fprintf(f,"{ ");
		for (int i = 1;i<=this->LUT_bins_ctheta;++i) {

			fprintf(f,"%g ",this->LUT_mom/this->Fit[phiang]->GetBinContent(i));

			if(i<this->LUT_bins_ctheta)
				fprintf(f,",");
		}
		if(phiang<this->LUT_bins_phi - 1)
			fprintf(f,"},\n    ");
		else
			fprintf(f,"}\n    ");
	}

	fprintf(f," };\n");
	fclose(f);
}

void ipa_class::init_extract_LUT() {
	this->myCanvas2 = root->newCanvas("IPA_PPCO","IPA - the Lookuptable extractor tool (press 'V' (shift + 'v'!) to update)",0,0,1680,1080);	
	this->myCanvas2->SetFixedAspectRatio(false);
	int hdiv=this->LUT_bins_phi / 6; 
	this->myCanvas2->Divide(6,hdiv);

	add_par(0.0,0.0,1000.0,1.0,"zero socket",70);
	add_par(1.0,1.0,4.0,1.0,"merge N bins",71);
	add_par(1.0,0.0,1.0,1.0,"sliding merge (yes=1)",72);
	add_par(this->LUT_mom*0.8,0.0,1000.0,0.1,"mom range min",73);
	add_par(this->LUT_mom*1.2,0.0,1000.0,0.1,"mom range max",74);
	add_par(this->LUT_mom,1.0,1000.0,1.0,"expected mom",75);
	add_par(0.0,0.0,1.0,1.0,"Gauss/Profile (G=0)",76);
}

void ipa_class::fill_extract_LUT(CH_event_struct *evti, int th) {
	int ehit = evti->e.num_hits;
	for(int i=0;i<(int)ehit;i++) {
		if(e[th][i]->valid) {
			int phiarr = (int)((e[0][i]->mom.Phi_deg()+180.0)/(360.0/(double)(this->LUT_bins_phi)));
			if(phiarr<this->LUT_bins_phi)
				ctpeArray[th][phiarr]->Fill(cos(e[th][i]->mom.Theta()),e[th][i]->mom.Mag());
		}
	}
}

void ipa_class::draw_extract_LUT() {

	char s[10];
	sprintf(s,"QNR");
	if(this->LUT_mergeNbins>1) {
			sprintf(s,"QG%i",this->LUT_mergeNbins);
		if(this->LUT_slidingmerge>0)
			sprintf(s,"QG%iS",this->LUT_mergeNbins);
	}

	//TF1 *F = new TF1("G","gaus(3)",this->LUT_momrange_min,this->LUT_momrange_max);
//	TF1 *F = 0;
//	F->SetRange(this->LUT_momrange_min,this->LUT_momrange_max);
	// do fitting + draw result
	for(int i=0;i<this->LUT_bins_phi;i++) {

		char nm[10];

		if(this->LUT_Profile>0) {
			sprintf(nm,"phi_%i_0_px",i);
			TCutG *cutg = new TCutG("mycut",3);
			cutg->SetVarX("y");
			cutg->SetVarY("x");
			cutg->SetPoint(0,-1.0,this->LUT_momrange_min);
			cutg->SetPoint(1,1.0,this->LUT_momrange_min);
			cutg->SetPoint(2,1.0,this->LUT_momrange_max);
			cutg->SetPoint(3,-1.0,this->LUT_momrange_max);
			ctpeArray[0][i]->ProfileX(nm,1,-1,"[mycut]");
			this->Fit[i] = (TH1D*)gDirectory->Get(nm);
		} else {
			sprintf(nm,"phi_%i_0_1",i);
			ctpeArray[0][i]->FitSlicesY(0,0,-1,this->LUT_zero,s);
			this->Fit[i] = (TH1D*)gDirectory->Get(nm);
		}
			
		this->Fit[i]->SetLineWidth(Width_t(3.0));
		myCanvas2->cd(i+1);
		ctpeArray[0][i]->Draw("COLZ");
		Fit[i]->Draw("same");
	}
	myCanvas2->Modified();
	myCanvas2->Update();

}
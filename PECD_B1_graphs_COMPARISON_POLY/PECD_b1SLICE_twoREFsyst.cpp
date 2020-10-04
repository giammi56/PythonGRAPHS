// PECD_graphs.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <stdio.h>
#include <atlstr.h>

#include <TApplication.h>

#include <TDirectory.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraph2D.h>
#include <THistPainter.h>
#include <THStack.h>
#include <TH2.h>
#include <TH1.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TPad.h>
#include <TPaveLabel.h>
#include <TRandom.h>
#include <TRefArray.H>
#include <TStyle.h>
#include <TString.h>
#include <TSystem.h>
#include <TText.h>
#include <TTree.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TError.h>

#include <iostream>
#include <vector>
#define NOPRINTSAVE 1
#define NOPRINTOUT 1
#define NOPRINTINLOAD 1
#define NOPRINTOUTLOAD 0
#define NOPRINTOUTLOAD2D 1

using namespace std;


void DrawTextInPad(double xCoord, double yCoord, char* text, double size) {
	TText* t1 = new TText(xCoord, yCoord, text);
	t1->SetNDC();
	t1->SetTextColor(1);
	t1->SetTextSize(size);
	t1->SetTextFont(62);
	t1->Draw();
}

void DividegPad(Int_t nx, Int_t ny,
	Float_t l, Float_t r, Float_t t, Float_t b)
{
	Int_t ix, iy, n=0;
	Double_t x1, x2, y1, y2;
	Double_t dx = ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l);
	Double_t dl = dx/(1-l);
	Double_t dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
	Double_t db = dy/(1-b);
	char *name  = new char [strlen(gPad->GetName())+6];

	y1 = 0;
	y2 = db;
	for (iy=0; iy<ny; iy++) {
		x1 = 0;
		x2 = dl;
		for (ix=0;ix<nx;ix++) {
			if (x1 > x2) continue;
			n++;
			sprintf(name,"%s_%d",gPad->GetName(),n);
			TPad* pad = new TPad(name,name,x1,y1,x2,y2,0);
			if (ix==0)    pad->SetLeftMargin(l);
			if (ix==nx-1) pad->SetRightMargin(r);
			if (iy==ny-1) pad->SetTopMargin(t);
			if (iy==0)    pad->SetBottomMargin(b);
			x1 = x2;
			if (ix==nx-2) x2 = 1;
			else          x2 = x1+dx;
			pad->SetNumber(n);
			pad->Draw();
		}
		y1 = y2;
		if (iy==ny-2) y2 = 1;
		else          y2 = y1+dy;
	}
}

void centrelabelaxis(TH2D* graph) {
	graph->GetXaxis()->SetTitleOffset(1);
	graph->GetYaxis()->SetTitleOffset(1);
	graph->GetZaxis()->SetTitleOffset(1);
}

void centrelabelaxis(TGraphAsymmErrors* graph) {
	graph->GetXaxis()->SetTitleOffset(1);
	graph->GetYaxis()->SetTitleOffset(1);
	//graph->GetZaxis()->SetTitleOffset(1);
}

void normalize_TH2D(TH2D* TH2D_Array) {

	Double_t norm = 1. / TH2D_Array->Integral();
	TH2D_Array->Scale(norm);
	//std::cout << TH2D_Array->GetName() << ": TH2D normalized \n" << endl << endl;

	return;
}

void savegraphas(char* tempjk, TH2D* asymmTH2D) {

	//#if NOPRINTSAVE == 1
	//	return;
	//#endif

	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	std::sprintf(tempjk, "%s.png", tempjk);
	cdummy->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
							  //#if NOPRINTOUT == 1
							  //	cdummy->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
							  //#endif
	cdummy->cd();
	centrelabelaxis(asymmTH2D);
	//asymmTH2D->GetZaxis()->SetLimits(-0.2,0.2);
	//asymmTH2D->GetZaxis()->SetRangeUser(-0.2,0.2);
	asymmTH2D->Draw();
	cdummy->SaveAs(tempjk);
	if (cdummy) { 
		cdummy->Close();
		gSystem->ProcessEvents();
	}

	//std::cout << asymmTH2D->GetName() << ": Exported as .png \n" << endl << endl;

	return;
}

void savegraphas(char* tempjk, TGraphAsymmErrors* ratio_genarray, TF1* fit1) {

#if NOPRINTSAVE == 1
	return;
#endif

	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	//TLegend *leg = cdummy->cd()->BuildLegend();
	//leg->SetTextFont(30);
	//leg->SetY1NDC(2);
	std::sprintf(tempjk, "%s.png", tempjk);
	cdummy->SetBatch(kFALSE);
#if NOPRINTOUT == 1
	cdummy->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
	cdummy->cd();
	ratio_genarray->SetTitle("");
	ratio_genarray->Draw("AP");
	fit1->Draw("SAME");
	cdummy->SaveAs(tempjk);
	if (cdummy) { 
		cdummy->Close();
		gSystem->ProcessEvents();
	}

	//std::cout << ratio_genarray->GetName() << ": Exported as .png \n" << endl << endl;

	return;
}

void savegraphas(char* tempjk, TH1D* ratioTH1D) {

#if NOPRINTSAVE == 1
	return;
#endif

	//ratioTH1D->SetTitle(tempjk);
	//ratioTH1D->SetTitle("");
	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	std::sprintf(tempjk, "%s.png", tempjk);
	//cdummy->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOPRINTOUT == 1
	cdummy->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
	cdummy->cd();
	ratioTH1D->SetMarkerStyle(20);
	ratioTH1D->SetMarkerSize(1);
	ratioTH1D->SetStats(0);
	ratioTH1D->GetYaxis()->SetTitle("b1");
	ratioTH1D->GetYaxis()->SetRangeUser(-0.3, 0.3);
	//ratioTH1D->Draw();
	ratioTH1D->Draw("E1");
	// (*aSlice)[0]->Draw();
	cdummy->SaveAs(tempjk);
	if (cdummy) { 
		cdummy->Close();
		gSystem->ProcessEvents();
	}

	//std::cout << asymmTH2D->GetName() << ": Exported as .png \n" << endl << endl;

	return;
}

//alternative method to define the fitting function
Double_t buildfunction(Double_t* x, Double_t* par) {
	Double_t arg = 0;
	arg = (x[0] - par[1])/par[2];
	Double_t f1new = par[0] * x[0] / (1 - 0.125*par[1] * x[0] * x[0]) + 1;
	return f1new;
}

TGraphAsymmErrors* ratio_TH1D(TH1D* Working_Array_1, TH1D* Working_Array_2) {
	int a = Working_Array_1->GetNbinsX();
	//std::cout << "a = " << a << endl << endl;

	// Asymmetry = (h1 - h2)/(h1 + h2)  where h1 = this
	auto ratioTH1D = new TGraphAsymmErrors(a);
	// https://en.wikipedia.org/wiki/Poisson_distribution#Assumptions:_When_is_the_Poisson_distribution_an_appropriate_model? //
	// second option in Divide() midp = lancaster mid-p_value, n = normal propagation
	ratioTH1D->Divide(Working_Array_1, Working_Array_2, "pois n"); //avarage value at 1
																   //ratioTH1D = (TGraphAsymmErrors*)Working_Array_1->GetAsymmetry(Working_Array_2);
																   //std::printf("Division of %s and %s in %s: DONE! Size ratio vector = %d \n",Working_Array_1->GetTitle(), Working_Array_2->GetTitle(), cname, ratioTH1D->Sizeof()/sizeof(int));

	return ratioTH1D;
}

TH1D* asymm_TH1D(TH1D* &Working_Array_1, TH1D* &Working_Array_2, const char* &cname) {

	TH1D* asymmTH1D = new TH1D();
	asymmTH1D = (TH1D*)Working_Array_1->GetAsymmetry(Working_Array_2); //GetAsymmetry is TH1*, cast to TH1D*
																	   //std::printf("Asymmetry of %s and %s in %s: DONE! Size ratio vector = %d \n",Working_Array_1->GetTitle(), Working_Array_2->GetTitle(), cname, asymmTH1D->Sizeof()/sizeof(int));

	return asymmTH1D;
}

TH2D* asymm_TH2D(TH2D* Working_Array2D_1, TH2D* Working_Array2D_2, const char* &cname) {

	//Working_Array2D_1->Sumw2();	
	normalize_TH2D(Working_Array2D_1);
	//Working_Array2D_2->Sumw2();
	normalize_TH2D(Working_Array2D_2);
	double asym_hist1_counts = Working_Array2D_1->GetEntries();
	double asym_hist2_counts = Working_Array2D_2->GetEntries();
	TH2D* asymmTH2D = (TH2D*)Working_Array2D_1->GetAsymmetry(Working_Array2D_2, asym_hist1_counts/asym_hist2_counts);

	////TEST
	//TH2D* top = (TH2D*)top->Add(Working_Array2D_1, Working_Array2D_2);
	//TH2D* bottom = (TH2D*)bottom->Add(Working_Array2D_1, Working_Array2D_2);
	//TH2D* asymmTH2D = (TH2D*)asymmTH2D->Divide(top,bottom);
	//cout << "Total number of entries = " << asymmTH2D->GetEntries() << endl;
	//cout << "Size of asymmTH2D: " << (asymmTH2D->Sizeof())/sizeof(int) << endl << endl;

	cout << "Division TH2D DONE!" << endl << endl;

	return asymmTH2D;
}

TH2D* adaptentries(TH2D* asymmTH2D) {

	//higer entries to allow FitSliceY()
	int nx = asymmTH2D->GetXaxis()->GetNbins();
	int ny = asymmTH2D->GetYaxis()->GetNbins();
	for(int x=0; x<nx; x++){
		for(int y=0; y<ny; y++){
			Double_t temp = asymmTH2D->GetBinContent(x, y);
			asymmTH2D->SetBinContent(x,y,temp * 10000);
			//cout << "New bin content: old = " << temp << " new = " << asymmTH2D->GetBinContent(x, y) << endl;
		}
	}

	//cout << "Total number of entries = " << asymmTH2D->GetEntries() << endl;
	//cout << "Size of asymmTH2D: " << (asymmTH2D->Sizeof())/sizeof(int) << endl << endl;
	//cout << "Division TH2D DONE!" << endl << endl;

	return asymmTH2D;
}

TH1D* adaptentries_red(TH1D* ratioTH1D) {

	int nx = ratioTH1D->GetXaxis()->GetNbins();
	for(int x=0; x<nx; x++){
		Double_t temp = ratioTH1D->GetBinContent(x);
		Double_t err = ratioTH1D->GetBinError(x);
		ratioTH1D->SetBinError(x,err / 100); // sqrt(10000)
		ratioTH1D->SetBinContent(x,temp / 10000);
		//cout << "New bin content: old = " << temp << " new = " << asymmTH2D->GetBinContent(x, y) << endl;
	}
	//cout << "Total number of entries = " << asymmTH2D->GetEntries() << endl;
	//cout << "Size of asymmTH2D: " << (asymmTH2D->Sizeof())/sizeof(int) << endl << endl;
	//cout << "Division TH2D DONE!" << endl << endl;

	return ratioTH1D;
}

// TAKE CARE THAT THE FEEDING PARAMETERS ARE IN TEH SAME ORDER!!!
void load_TH1D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH1D** S_LCP_TH1D_Array, TH1D** S_RCP_TH1D_Array, TH1D** R_LCP_TH1D_Array, TH1D** R_RCP_TH1D_Array, TH1D** general_Array, TString* dirlist, TString* histlist1D, const char dir[]) {
	int counter = 0;
	char name[200];
	for (int i = 0; i < 2; i++) {
		//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
		TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory(dir+dirlist[i]));
		TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory(dir+dirlist[i]));
		TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory(dir+dirlist[i]));
		TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory(dir+dirlist[i]));
		gDirectory->pwd();
		//for (int j = 0; j < LCPhistlist->Sizeof(); j++) {

		//12x6
		for (int k = 0; k < 72; k++) {
			std::sprintf(name, "S_RCP_%d_%d",k,i);
			S_RCP_TH1D_Array[k+i*72] = (TH1D*)S_Rdir->FindObjectAny(histlist1D[k]);
			S_RCP_TH1D_Array[k+i*72]->SetName(name);
			S_RCP_TH1D_Array[k+i*72]->SetTitle(name);
			S_RCP_TH1D_Array[k+i*72]->Sumw2();
			S_RCP_TH1D_Array[k+i*72]->Scale(1. / S_RCP_TH1D_Array[k]->Integral());
			//std::cout << "Loading and normalisation of " << S_RCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << S_RCP_TH1D_Array[k]->GetEntries() << endl << endl;
			general_Array[counter++] = S_RCP_TH1D_Array[k+i*72];

			std::sprintf(name, "S_LCP_%d_%d",k,i);
			S_LCP_TH1D_Array[k+i*72] = (TH1D*)S_Ldir->FindObjectAny(histlist1D[k]);
			S_LCP_TH1D_Array[k+i*72]->SetName(name);
			S_LCP_TH1D_Array[k+i*72]->SetTitle(name);
			S_LCP_TH1D_Array[k+i*72]->Sumw2();
			S_LCP_TH1D_Array[k+i*72]->Scale(1. / S_LCP_TH1D_Array[k+i*72]->Integral());
			//std::cout << "Loading and normalisation of " << S_LCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << S_LCP_TH1D_Array[k]->GetEntries() << endl << endl;
			general_Array[counter++] = S_LCP_TH1D_Array[k+i*72];

			std::sprintf(name, "R_RCP_%d_%d",k,i);
			R_RCP_TH1D_Array[k+i*72] = (TH1D*)R_Rdir->FindObjectAny(histlist1D[k]);
			R_RCP_TH1D_Array[k+i*72]->SetName(name);
			R_RCP_TH1D_Array[k+i*72]->SetTitle(name);
			R_RCP_TH1D_Array[k+i*72]->Sumw2();
			R_RCP_TH1D_Array[k+i*72]->Scale(1. / R_RCP_TH1D_Array[k+i*72]->Integral());
			//std::cout << "Loading and normalisation of " << R_RCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << R_RCP_TH1D_Array[k]->GetEntries() << endl << endl;
			general_Array[counter++] = R_RCP_TH1D_Array[k+i*72];

			std::sprintf(name, "R_LCP_%d_%d",k,i);
			R_LCP_TH1D_Array[k+i*72] = (TH1D*)R_Ldir->FindObjectAny(histlist1D[k]);
			R_LCP_TH1D_Array[k+i*72]->SetName(name);
			R_LCP_TH1D_Array[k+i*72]->SetTitle(name);
			R_LCP_TH1D_Array[k+i*72]->Sumw2();
			R_LCP_TH1D_Array[k+i*72]->Scale(1. / R_LCP_TH1D_Array[k+i*72]->Integral());
			//std::cout << "Loading and normalisation of " << R_LCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << R_LCP_TH1D_Array[k]->GetEntries() << endl << endl;
			general_Array[counter++] = R_LCP_TH1D_Array[k+i*72];
		}

		std::cout << "Total number of TH1D histos = " << counter << endl;

		//Autoput of all loaded graphs
		const char *cnameinternal[2][4] = {{"cos(theta)_S_RCP","cos(theta)_S_LCP","cos(theta)_R_RCP","cos(theta)_R_LCP"},
		{"cos(theta)_S_RCP_flip","cos(theta)_S_LCP_flip","cos(theta)_R_RCP_flip","cos(theta)_R_LCP_flip"}};

		//plotting all cos(theta) plots
		for (int j = 0; j < 4; j++) {
			//char canvastitle[100];
			//std::sprintf(canvastitle, "all entries = %d", counter / 4);
			char canvasname[200];
			std::sprintf(canvasname, "All %s plots", cnameinternal[i][j]);
			TCanvas* cplot = new TCanvas(canvasname, canvasname);
			cplot->SetCanvasSize(1920, 1080);
			cplot->SetWindowSize(1920 + 4, 1080 + 28);
			cplot->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOPRINTINLOAD == 1
			cplot->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif 
									//12x6
			cplot->Divide(12, 6);
			for (int k = 0; k < 72; k++) {
				cplot->cd(k+1);
				//general_Array[k]->SetOption("E");
				general_Array[k*4+j+i*288]->Draw("E1");
				//gPad->Modified();
				gPad->Update();
			}
			cplot->Write();
			cplot->Close();
			gSystem->ProcessEvents();
		}
	}
	std::cout << endl;
	return;
}

void load_TH2D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH2D** S_LCP_TH2D_Array, TH2D** S_RCP_TH2D_Array, TH2D** R_LCP_TH2D_Array, TH2D** R_RCP_TH2D_Array, TH2D** general_Array2D, TString* dirlist, TString* histlist2D, const char dir[]) {
	int counter = 0;
	char name[200];
	for (int i = 0; i < 1; i++) {
		//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
		TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory(dir+dirlist[i]));
		TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory(dir+dirlist[i]));
		TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory(dir+dirlist[i]));
		TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory(dir+dirlist[i]));

		gDirectory->pwd();

		for (int k = 0; k < 72; k++) {
			std::sprintf(name, "MFPAD_S_RCP_%d_%d",k,i);
			S_RCP_TH2D_Array[k] = (TH2D*)S_Rdir->FindObjectAny(histlist2D[k]);
			S_RCP_TH2D_Array[k]->SetName(name);
			S_RCP_TH2D_Array[k]->SetTitle(name);
			//S_RCP_TH2D_Array[k+i*3]->Sumw2();													 
			//normalize_TH2D(S_RCP_TH2D_Array[k]);
			general_Array2D[counter++] = S_RCP_TH2D_Array[k];

			std::sprintf(name, "MFPAD_S_LCP_%d_%d",k,i);
			S_LCP_TH2D_Array[k] = (TH2D*)S_Ldir->FindObjectAny(histlist2D[k]);
			S_LCP_TH2D_Array[k]->SetName(name);
			S_LCP_TH2D_Array[k]->SetTitle(name);
			//S_LCP_TH2D_Array[k+i*3]->Sumw2();
			//normalize_TH2D(S_LCP_TH2D_Array[k]);
			general_Array2D[counter++] = S_LCP_TH2D_Array[k];

			std::sprintf(name, "MFPAD_R_RCP_%d_%d",k,i);
			R_RCP_TH2D_Array[k] = (TH2D*)R_Rdir->FindObjectAny(histlist2D[k]);
			R_RCP_TH2D_Array[k]->SetName(name);
			R_RCP_TH2D_Array[k]->SetTitle(name);
			//R_RCP_TH2D_Array[k+i*3]->Sumw2();													 
			//normalize_TH2D(R_RCP_TH2D_Array[k]);
			general_Array2D[counter++] = R_RCP_TH2D_Array[k];

			std::sprintf(name, "MFPAD_R_LCP_%d_%d",k,i);
			R_LCP_TH2D_Array[k] = (TH2D*)R_Ldir->FindObjectAny(histlist2D[k]);
			R_LCP_TH2D_Array[k]->SetName(name);
			R_LCP_TH2D_Array[k]->SetTitle(name);
			//R_LCP_TH2D_Array[k+i*3]->Sumw2();													 
			//normalize_TH2D(R_LCP_TH2D_Array[k]);
			general_Array2D[counter++] = R_LCP_TH2D_Array[k];
		}

		std::cout << "Total number of TH2D histos = " << counter << endl;

		const char *cnameinternal[2][4] = {{"MFPAD_S_RCP","MFPAD_S_LCP","MFPAD_R_RCP","MFPAD_R_LCP"},
		{"MFPAD_S_RCP_flip","MFPAD_S_LCP_flip","MFPAD_R_RCP_flip","MFPAD_R_LCP_flip"}};

		for (int j = 0; j < 4; j++) {
			char canvasname[200];
			std::sprintf(canvasname, "All %s plots", cnameinternal[i][j]);
			TCanvas* cplot = new TCanvas(canvasname, canvasname);
			cplot->SetCanvasSize(1920, 1080);
			cplot->SetWindowSize(1920 + 4, 1080 + 28);
			cplot->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOPRINTINLOAD == 1
			cplot->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
			cplot->Divide(12,6);
			for (int k = 0; k < 72; k++) {
				cplot->cd(k+1);
				//general_Array[k]->SetOption("E");
				general_Array2D[k*4+j+i*288]->Draw("COLZ");
				//gPad->Modified();
				gPad->Update();
			}
			cplot->Write();
			cplot->Close();
			gSystem->ProcessEvents();
		}
	}

	std::cout << endl;

	return;

}

//int main(int argc, char** argv)
int main(int argc, char* argv[])

{
	TApplication theApp("App", &argc, argv);

	// single TH1D arrays: adapt to the requested dimension i*j
	int dim = 576; //12x6x4x2
	TH1D** S_LCP_TH1D_Array;
	S_LCP_TH1D_Array = new TH1D*[dim/4];

	TH1D** R_LCP_TH1D_Array;
	R_LCP_TH1D_Array = new TH1D*[dim/4];

	TH1D** S_RCP_TH1D_Array;
	S_RCP_TH1D_Array = new TH1D*[dim/4];

	TH1D** R_RCP_TH1D_Array;
	R_RCP_TH1D_Array = new TH1D*[dim/4];

	TH1D** Working_Array_1;
	Working_Array_1 = new TH1D*[dim];

	TH1D** Working_Array_2;
	Working_Array_2 = new TH1D*[dim];

	TH1D** general_Array;
	general_Array = new TH1D*[dim];

	TGraphAsymmErrors** ratio_genarray;
	ratio_genarray = new TGraphAsymmErrors*[dim];

	TH1D** asymm_genarray;
	asymm_genarray = new TH1D*[dim];

	// TH2D arrays
	//12x6x4x2
	TH2D** S_LCP_TH2D_Array;
	S_LCP_TH2D_Array = new TH2D*[dim/4];

	TH2D** R_LCP_TH2D_Array;
	R_LCP_TH2D_Array = new TH2D*[dim/4];

	TH2D** S_RCP_TH2D_Array;
	S_RCP_TH2D_Array = new TH2D*[dim/4];

	TH2D** R_RCP_TH2D_Array;
	R_RCP_TH2D_Array = new TH2D*[dim/4];

	TH2D** Working_Array2D_1;
	Working_Array2D_1 = new TH2D*[dim];

	TH2D** Working_Array2D_2;
	Working_Array2D_2 = new TH2D*[dim];

	TH2D** general_Array2D;
	general_Array2D = new TH2D*[dim];

	TH2D** asymmTH2D;
	asymmTH2D = new TH2D*[dim];

	TH1D** copyaSliceTH1D;
	copyaSliceTH1D = new TH1D*[dim];

	// Fititign parameters
	TF1** fit1;
	fit1 = new TF1*[];
	//PECD fitting
	TF1 *f1 = new TF1("PECD","([0]*x) / (1 + 0.5*[1]*(3*x*x-1)) +1", -0.9, 0.9);; //1D PECD fitting mod (+1 for the divide) NB. the coeffiecient on b2 could be wrong (0.5)
	f1->SetParameters(0,0);
	f1->SetParLimits(0,-1,1);
	//f1->SetParLimits(1,-6.5,6.5);
	f1->SetParNames("b1","b2");

	TF1 *f2 = new TF1("PECD4","([0]*x + 0.5*[2]*(5*x*x*x-3*x)) / (1 + 0.5*[1]*(3*x*x-1) + 0.125*[3]*(35*x*x*x*x-30*x*x+3)) +1", -0.9, 0.9);; //1D PECD fitting mod (+1 for the divide)
	f2->SetParameters(0,0,0,0);
	f2->SetParLimits(0,-1,1);
	f2->SetParLimits(1,-1.5,1.5);
	f2->SetParLimits(2,-1.5,1.5);
	f2->SetParLimits(3,-1.5,1.5);
	f2->SetParNames("b1","b2","b3","b4");

	TF1 *f3 = new TF1("PECD6","([0]*x + 0.5*[2]*(5*x*x*x-3*x) +0.125*[4]*(63*x*x*x*x*x-70*x*x*x+15*x)) / (1 + 0.5*[1]*(3*x*x-1) + 0.125*[3]*(35*x*x*x*x-30*x*x+3) + 0.0625*[5]*(231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5)) +1", -0.9, 0.9);; //1D PECD fitting mod (+1 for the divide)
	f3->SetParameters(0,0,0,0,0,0);
	f3->SetParLimits(0,-1,1);
	f3->SetParLimits(1,-1,1);
	f3->SetParLimits(2,-1,1);
	f3->SetParLimits(3,-1,1);
	f3->SetParLimits(4,-1,1);
	f3->SetParLimits(5,-1,1);
	f3->SetParNames("b1","b2","b3","b4","b5","b6");
	//Double_t* x;
	//x = new Double_t[];
	//Double_t* par;
	//par = new Double_t[];

	//Adapt dimensio to request
	Double_t b1[288];
	Double_t b1_err[288];

	double cos_theta[288/4];
	double phi[288/4];

	TGraph2D** b1_map;
	b1_map = new TGraph2D*[3];

	TH2D** b1_map_copy;
	b1_map_copy = new TH2D*[3];

	Int_t MyPalette[200];
	Double_t r[]    = {0.0, 1.0, 1.0, 1.0};
	Double_t g[]    = {0.0, 1.0, 1.0, 0.0};
	Double_t b[]    = {1.0, 1.0, 0.0, 0.0};
	Double_t stop[] = {0.0, 0.5, 0.75, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(4, stop, r, g, b, 200);
	for (int i=0;i<200;i++) MyPalette[i] = FI+i;
	gStyle->SetPalette(200,MyPalette);

	//// DO NOT CHANGE THE ORDER OF THESE INPUTS R-C3H3F3O_546eV_CL_10800-2350ns_multiCH11_MFPAD_30.root
	////550 eV CH11
	//TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	//if (S_LFile->IsOpen()) printf("%s: File opened successfully!\n", S_LFile->GetName());
	//TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	//if (S_RFile->IsOpen()) printf("%s: File opened successfully!\n", S_RFile->GetName());
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	//if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	//if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//550eV CH9
	TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CL_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root", "READ");
	if (S_LFile->IsOpen()) printf("%s: File opened successfully!\n", S_LFile->GetName());
	TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CR_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root", "READ");
	if (S_RFile->IsOpen()) printf("%s: File opened successfully!\n", S_RFile->GetName());
	TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root", "READ");
	if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_9600-3700ns_newEfield_multiCH9_MFPAD_30_t0.root", "READ");
	if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//// DO NOT CHANGE THE ORDER OF THESE INPUTS
	//TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/S-C3H3F3O_poly_30-69m_550eV_CL_reaction_POLY_3x3_test.root", "READ");
	//if (S_LFile->IsOpen()) printf("%s: File opened successfully!\n", S_LFile->GetName());
	//TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/S-C3H3F3O_poly_30-69m_550eV_CR_reaction_POLY_3x3_test.root", "READ");
	//if (S_RFile->IsOpen()) printf("%s: File opened successfully!\n", S_RFile->GetName());
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CL_reaction_POLY_3x3_test.root", "READ");
	//if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CR_reaction_POLY_3x3_test.root", "READ");
	//if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//File to save outout b1
	FILE *fp = fopen("b1_550eV_ID_ALL_CR-CL_CH9m69_12x8.dat","w+");
	//FILE *fp = fopen("b1_ID_ALL_m69_3x3_test.dat","w");
	fprintf(fp,"cos(theta) \t phi \t p1 \t p1_err \n");

	//FILE NAME
	// POLYATOMIC - THREE CHARGED FRAGMENTS
	//TFile *f = new TFile("FMtox_ID_ALL_m69_3x3_test.root","RECREATE");
	//TFile *f = new TFile("FMtox_ID3_m41_69.root", "RECREATE");
	//POLIATOMIC - NEUTRAL
	TFile *f = new TFile("FMtox_550eV_ID_ALL_CR-CL_CH9m69_12x6.root","RECREATE");

	if (f->IsOpen()) {
		printf("File CREATED successfully\n");
	}

	//const char dir[] = {"angular_distr_el/CH11/ID_ALL_mol_e0_valid"};
	const char dir[] = {"angular_distr_el/CH9/ID_ALL_mol_e0_valid"};
	TString dirlist[] = {"/EN_gate/MFPADs_multinew","/EN_gate/MFPADs_multinew_flipped"};
	//12x6 = 72
	TString histlist1D[] = { "cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-0.67_phi_-180","cos(theta)_e[0]_costheta_-0.33_phi_-180","cos(theta)_e[0]_costheta_-0.00_phi_-180","cos(theta)_e[0]_costheta_0.33_phi_-180","cos(theta)_e[0]_costheta_0.67_phi_-180",
		"cos(theta)_e[0]_costheta_-1.00_phi_-150","cos(theta)_e[0]_costheta_-0.67_phi_-150","cos(theta)_e[0]_costheta_-0.33_phi_-150","cos(theta)_e[0]_costheta_-0.00_phi_-150","cos(theta)_e[0]_costheta_0.33_phi_-150","cos(theta)_e[0]_costheta_0.67_phi_-150",
		"cos(theta)_e[0]_costheta_-1.00_phi_-120","cos(theta)_e[0]_costheta_-0.67_phi_-120","cos(theta)_e[0]_costheta_-0.33_phi_-120","cos(theta)_e[0]_costheta_-0.00_phi_-120","cos(theta)_e[0]_costheta_0.33_phi_-120","cos(theta)_e[0]_costheta_0.67_phi_-120",
		"cos(theta)_e[0]_costheta_-1.00_phi_-90", "cos(theta)_e[0]_costheta_-0.67_phi_-90", "cos(theta)_e[0]_costheta_-0.33_phi_-90", "cos(theta)_e[0]_costheta_-0.00_phi_-90", "cos(theta)_e[0]_costheta_0.33_phi_-90", "cos(theta)_e[0]_costheta_0.67_phi_-90",
		"cos(theta)_e[0]_costheta_-1.00_phi_-60", "cos(theta)_e[0]_costheta_-0.67_phi_-60", "cos(theta)_e[0]_costheta_-0.33_phi_-60", "cos(theta)_e[0]_costheta_-0.00_phi_-60", "cos(theta)_e[0]_costheta_0.33_phi_-60", "cos(theta)_e[0]_costheta_0.67_phi_-60",
		"cos(theta)_e[0]_costheta_-1.00_phi_-30", "cos(theta)_e[0]_costheta_-0.67_phi_-30", "cos(theta)_e[0]_costheta_-0.33_phi_-30", "cos(theta)_e[0]_costheta_-0.00_phi_-30", "cos(theta)_e[0]_costheta_0.33_phi_-30", "cos(theta)_e[0]_costheta_0.67_phi_-30",
		"cos(theta)_e[0]_costheta_-1.00_phi_0",   "cos(theta)_e[0]_costheta_-0.67_phi_0",   "cos(theta)_e[0]_costheta_-0.33_phi_0",   "cos(theta)_e[0]_costheta_-0.00_phi_0",   "cos(theta)_e[0]_costheta_0.33_phi_0",   "cos(theta)_e[0]_costheta_0.67_phi_0",
		"cos(theta)_e[0]_costheta_-1.00_phi_30",  "cos(theta)_e[0]_costheta_-0.67_phi_30",  "cos(theta)_e[0]_costheta_-0.33_phi_30",  "cos(theta)_e[0]_costheta_-0.00_phi_30",  "cos(theta)_e[0]_costheta_0.33_phi_30",  "cos(theta)_e[0]_costheta_0.67_phi_30",
		"cos(theta)_e[0]_costheta_-1.00_phi_60",  "cos(theta)_e[0]_costheta_-0.67_phi_60",  "cos(theta)_e[0]_costheta_-0.33_phi_60",  "cos(theta)_e[0]_costheta_-0.00_phi_60",  "cos(theta)_e[0]_costheta_0.33_phi_60",  "cos(theta)_e[0]_costheta_0.67_phi_60",
		"cos(theta)_e[0]_costheta_-1.00_phi_90",  "cos(theta)_e[0]_costheta_-0.67_phi_90",  "cos(theta)_e[0]_costheta_-0.33_phi_90",  "cos(theta)_e[0]_costheta_-0.00_phi_90",  "cos(theta)_e[0]_costheta_0.33_phi_90",  "cos(theta)_e[0]_costheta_0.67_phi_90",
		"cos(theta)_e[0]_costheta_-1.00_phi_120", "cos(theta)_e[0]_costheta_-0.67_phi_120", "cos(theta)_e[0]_costheta_-0.33_phi_120", "cos(theta)_e[0]_costheta_-0.00_phi_120", "cos(theta)_e[0]_costheta_0.33_phi_120", "cos(theta)_e[0]_costheta_0.67_phi_120",
		"cos(theta)_e[0]_costheta_-1.00_phi_150", "cos(theta)_e[0]_costheta_-0.67_phi_150", "cos(theta)_e[0]_costheta_-0.33_phi_150", "cos(theta)_e[0]_costheta_-0.00_phi_150", "cos(theta)_e[0]_costheta_0.33_phi_150", "cos(theta)_e[0]_costheta_0.67_phi_150" };

	////9x8
	//TString histlist1D[72] = {  "cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-1.00_phi_-135","cos(theta)_e[0]_costheta_-1.00_phi_-90","cos(theta)_e[0]_costheta_-1.00_phi_-45","cos(theta)_e[0]_costheta_-1.00_phi_0","cos(theta)_e[0]_costheta_-1.00_phi_45","cos(theta)_e[0]_costheta_-1.00_phi_90","cos(theta)_e[0]_costheta_-1.00_phi_135",
	//							  "cos(theta)_e[0]_costheta_-0.78_phi_-180","cos(theta)_e[0]_costheta_-0.78_phi_-135","cos(theta)_e[0]_costheta_-0.78_phi_-90","cos(theta)_e[0]_costheta_-0.78_phi_-45","cos(theta)_e[0]_costheta_-0.78_phi_0","cos(theta)_e[0]_costheta_-0.78_phi_45","cos(theta)_e[0]_costheta_-0.78_phi_90","cos(theta)_e[0]_costheta_-0.78_phi_135",
	//							  "cos(theta)_e[0]_costheta_-0.56_phi_-180","cos(theta)_e[0]_costheta_-0.56_phi_-135","cos(theta)_e[0]_costheta_-0.56_phi_-90","cos(theta)_e[0]_costheta_-0.56_phi_-45","cos(theta)_e[0]_costheta_-0.56_phi_0","cos(theta)_e[0]_costheta_-0.56_phi_45","cos(theta)_e[0]_costheta_-0.56_phi_90","cos(theta)_e[0]_costheta_-0.56_phi_135",
	//							  "cos(theta)_e[0]_costheta_-0.11_phi_-180","cos(theta)_e[0]_costheta_-0.11_phi_-135","cos(theta)_e[0]_costheta_-0.11_phi_-90","cos(theta)_e[0]_costheta_-0.11_phi_-45","cos(theta)_e[0]_costheta_-0.11_phi_0","cos(theta)_e[0]_costheta_-0.11_phi_45","cos(theta)_e[0]_costheta_-0.11_phi_90","cos(theta)_e[0]_costheta_-0.11_phi_135",
	//							  "cos(theta)_e[0]_costheta_-0.33_phi_-180","cos(theta)_e[0]_costheta_-0.33_phi_-135","cos(theta)_e[0]_costheta_-0.33_phi_-90","cos(theta)_e[0]_costheta_-0.33_phi_-45","cos(theta)_e[0]_costheta_-0.33_phi_0","cos(theta)_e[0]_costheta_-0.33_phi_45","cos(theta)_e[0]_costheta_-0.33_phi_90","cos(theta)_e[0]_costheta_-0.33_phi_135",
	//							  "cos(theta)_e[0]_costheta_0.11_phi_-180","cos(theta)_e[0]_costheta_0.11_phi_-135","cos(theta)_e[0]_costheta_0.11_phi_-90","cos(theta)_e[0]_costheta_0.11_phi_-45","cos(theta)_e[0]_costheta_0.11_phi_0","cos(theta)_e[0]_costheta_0.11_phi_45","cos(theta)_e[0]_costheta_0.11_phi_90","cos(theta)_e[0]_costheta_0.11_phi_135",
	//							  "cos(theta)_e[0]_costheta_0.33_phi_-180","cos(theta)_e[0]_costheta_0.33_phi_-135","cos(theta)_e[0]_costheta_0.33_phi_-90","cos(theta)_e[0]_costheta_0.33_phi_-45","cos(theta)_e[0]_costheta_0.33_phi_0","cos(theta)_e[0]_costheta_0.33_phi_45","cos(theta)_e[0]_costheta_0.33_phi_90","cos(theta)_e[0]_costheta_0.33_phi_135",
	//							  "cos(theta)_e[0]_costheta_0.55_phi_-180","cos(theta)_e[0]_costheta_0.55_phi_-135","cos(theta)_e[0]_costheta_0.55_phi_-90","cos(theta)_e[0]_costheta_0.55_phi_-45","cos(theta)_e[0]_costheta_0.55_phi_0","cos(theta)_e[0]_costheta_0.55_phi_45","cos(theta)_e[0]_costheta_0.55_phi_90","cos(theta)_e[0]_costheta_0.55_phi_135",
	//							  "cos(theta)_e[0]_costheta_0.78_phi_-180","cos(theta)_e[0]_costheta_0.78_phi_-135","cos(theta)_e[0]_costheta_0.78_phi_-90","cos(theta)_e[0]_costheta_0.78_phi_-45","cos(theta)_e[0]_costheta_0.78_phi_0","cos(theta)_e[0]_costheta_0.78_phi_45","cos(theta)_e[0]_costheta_0.78_phi_90","cos(theta)_e[0]_costheta_0.78_phi_135"};

	////4x4
	//TString histlist1D[16] = {"cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-1.00_phi_-90","cos(theta)_e[0]_costheta_-1.00_phi_0","cos(theta)_e[0]_costheta_-1.00_phi_90",
	//						  "cos(theta)_e[0]_costheta_-0.50_phi_-180","cos(theta)_e[0]_costheta_-0.50_phi_-90","cos(theta)_e[0]_costheta_-0.50_phi_0","cos(theta)_e[0]_costheta_-0.50_phi_90",
	//						  "cos(theta)_e[0]_costheta_0.00_phi_-180","cos(theta)_e[0]_costheta_0.00_phi_-90","cos(theta)_e[0]_costheta_0.00_phi_0","cos(theta)_e[0]_costheta_0.00_phi_90",
	//						  "cos(theta)_e[0]_costheta_0.50_phi_-180","cos(theta)_e[0]_costheta_0.50_phi_-90","cos(theta)_e[0]_costheta_0.50_phi_0","cos(theta)_e[0]_costheta_0.50_phi_90"};
	//TString histlist1D[9] = {"cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-1.00_phi_-60","cos(theta)_e[0]_costheta_-1.00_phi_60",
	//						  "cos(theta)_e[0]_costheta_-0.33_phi_-180","cos(theta)_e[0]_costheta_-0.33_phi_-60","cos(theta)_e[0]_costheta_-0.33_phi_60",
	//						  "cos(theta)_e[0]_costheta_0.33_phi_-180","cos(theta)_e[0]_costheta_0.33_phi_-60","cos(theta)_e[0]_costheta_0.33_phi_60"};

	//TString histlist2D[1] = {"B1_map"};
	//12x6 = 72
	TString histlist2D[] = { "MFPAD3D_engate_costheta_-1.00_phi_-180", "MFPAD3D_engate_costheta_-0.67_phi_-180","MFPAD3D_engate_costheta_-0.33_phi_-180","MFPAD3D_engate_costheta_-0.00_phi_-180", "MFPAD3D_engate_costheta_0.33_phi_-180", "MFPAD3D_engate_costheta_-0.67_phi_-180",
		"MFPAD3D_engate_costheta_-1.00_phi_-150","MFPAD3D_engate_costheta_-0.67_phi_-150","MFPAD3D_engate_costheta_-0.33_phi_-150","MFPAD3D_engate_costheta_-0.00_phi_-150","MFPAD3D_engate_costheta_0.33_phi_-150","MFPAD3D_engate_costheta_-0.67_phi_-150",
		"MFPAD3D_engate_costheta_-1.00_phi_-120","MFPAD3D_engate_costheta_-0.67_phi_-120","MFPAD3D_engate_costheta_-0.33_phi_-120","MFPAD3D_engate_costheta_-0.00_phi_-120","MFPAD3D_engate_costheta_0.33_phi_-120","MFPAD3D_engate_costheta_-0.67_phi_-120",
		"MFPAD3D_engate_costheta_-1.00_phi_-90", "MFPAD3D_engate_costheta_-0.67_phi_-90", "MFPAD3D_engate_costheta_-0.33_phi_-90", "MFPAD3D_engate_costheta_-0.00_phi_-90", "MFPAD3D_engate_costheta_0.33_phi_-90", "MFPAD3D_engate_costheta_-0.67_phi_-90",
		"MFPAD3D_engate_costheta_-1.00_phi_-60", "MFPAD3D_engate_costheta_-0.67_phi_-60", "MFPAD3D_engate_costheta_-0.33_phi_-60", "MFPAD3D_engate_costheta_-0.00_phi_-60", "MFPAD3D_engate_costheta_0.33_phi_-60", "MFPAD3D_engate_costheta_-0.67_phi_-60",
		"MFPAD3D_engate_costheta_-1.00_phi_-30", "MFPAD3D_engate_costheta_-0.67_phi_-30", "MFPAD3D_engate_costheta_-0.33_phi_-30", "MFPAD3D_engate_costheta_-0.00_phi_-30", "MFPAD3D_engate_costheta_0.33_phi_-30", "MFPAD3D_engate_costheta_-0.67_phi_-30",
		"MFPAD3D_engate_costheta_-1.00_phi_0",	 "MFPAD3D_engate_costheta_-0.67_phi_0",   "MFPAD3D_engate_costheta_-0.33_phi_0",   "MFPAD3D_engate_costheta_-0.00_phi_0",   "MFPAD3D_engate_costheta_0.33_phi_0",   "MFPAD3D_engate_costheta_-0.67_phi_0",
		"MFPAD3D_engate_costheta_-1.00_phi_30",	 "MFPAD3D_engate_costheta_-0.67_phi_30",  "MFPAD3D_engate_costheta_-0.33_phi_30",  "MFPAD3D_engate_costheta_-0.00_phi_30",  "MFPAD3D_engate_costheta_0.33_phi_30",  "MFPAD3D_engate_costheta_-0.67_phi_30",
		"MFPAD3D_engate_costheta_-1.00_phi_60",	 "MFPAD3D_engate_costheta_-0.67_phi_60",  "MFPAD3D_engate_costheta_-0.33_phi_60",  "MFPAD3D_engate_costheta_-0.00_phi_60",  "MFPAD3D_engate_costheta_0.33_phi_60",  "MFPAD3D_engate_costheta_-0.67_phi_60",
		"MFPAD3D_engate_costheta_-1.00_phi_90",	 "MFPAD3D_engate_costheta_-0.67_phi_90",  "MFPAD3D_engate_costheta_-0.33_phi_90",  "MFPAD3D_engate_costheta_-0.00_phi_90",  "MFPAD3D_engate_costheta_0.33_phi_90",  "MFPAD3D_engate_costheta_-0.67_phi_90",
		"MFPAD3D_engate_costheta_-1.00_phi_120", "MFPAD3D_engate_costheta_-0.67_phi_120", "MFPAD3D_engate_costheta_-0.33_phi_120", "MFPAD3D_engate_costheta_-0.00_phi_120", "MFPAD3D_engate_costheta_0.33_phi_120", "MFPAD3D_engate_costheta_-0.67_phi_120",
		"MFPAD3D_engate_costheta_-1.00_phi_150", "MFPAD3D_engate_costheta_-0.67_phi_150", "MFPAD3D_engate_costheta_-0.33_phi_150", "MFPAD3D_engate_costheta_-0.00_phi_150", "MFPAD3D_engate_costheta_0.33_phi_150", "MFPAD3D_engate_costheta_-0.67_phi_150"};


	//3 canvas strig name, from Till suggestion! Not necessary now.
	const char *cname[] = {"b1fitR","b1fitS","b1fitRS","b1mapR","b1mapS","b1mapRS"};
	//const char *cname[] = {"PECD_e[0]_R_CR-CL","PECD_e[0]_S_CR-CL","PECD_e[0]_R_CL-S_CR","b1_map_R_CR-CL","b1_map_S_CR-CL","b1_map_R_CL-S_CR"};

	const char *cname2D[] = {"PECD_e[0]_R_CR-CL","PECD_e[0]_S_CR-CL","PECD_e[0]_R-S_CL-CR","PECD_e[0]_R_CR-CL_smooth","PECD_e[0]_S_CR-CL_smooth","PECD_e[0]_R-S_CL-CR_smooth"};
	TCanvas *canv[6]; // dim(3 cos(theta) + 3 b1__maps) - > larger for debugging
	TCanvas *canv2D[3]; // dim(3 MFPAD) - > larger for debugging
	TCanvas *canv2D_smooth[3]; // dim(3 MFPAD) - > larger for debugging
	TCanvas *canv_fit[3];
	std::cout << "Vector of TCanvas created. \n" << endl;

	char tempj[500];
	char tempk[500];
	char tempjk[500];

	// LOAD AND NORMALIZED ALL GRAPHS
	load_TH1D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH1D_Array, S_RCP_TH1D_Array, R_LCP_TH1D_Array, R_RCP_TH1D_Array, general_Array, dirlist, histlist1D, dir);
	std::cout << "load_all TH1D : DONE! \n" << endl;

	load_TH2D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH2D_Array, S_RCP_TH2D_Array, R_LCP_TH2D_Array, R_RCP_TH2D_Array, general_Array2D, dirlist, histlist2D, dir);
	std::cout << "load_all TH2D : DONE! \n" << endl;

	//Suppress INFO messages ROOT
	//gErrorIgnoreLevel = kWarning;

	//generation of phi and cos(theta) 16 bins: be sure that the increment and the incremented variabel is the same at in histlist1D 
	// 12x6 : 12 cos(theta), 6 phi
	int counter = 0;
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 6; j++) {
			phi[counter] = i*30-180+15;
			cos_theta[counter] = j*0.333-1+0.1667; // (bin-period/2+bin/2)
			counter++;
		}
	}

	int rcounter = 0;
	//3 canvas of i*j plots each cos(theta) + 3 cavas for b1_map
	for (int j = 0; j < 6; j++) {
		if (j < 3) {
			canv[j] = new TCanvas(cname[j], cname[j]);
			canv[j]->SetCanvasSize(1920, 1080);
			canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
#if NOPRINTOUTLOAD == 1
			canv[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
									  //12x6
			canv[j]->Divide(12, 6);						  
			canv[j]->SetFillStyle(4000);

			// all cos(theta) -> counter goes for ixj entries
			//12X6
			for (int k = 0; k < 72; k++) {
				std::sprintf(tempj, "clean");
				std::sprintf(tempk, "clean");
				if (j == 0) { // R enantiomer, different ellicity
					std::sprintf(tempk, "%s", R_RCP_TH1D_Array[k]->GetTitle());
					std::sprintf(tempj, "%s", R_LCP_TH1D_Array[k]->GetTitle());
					//Working_Array_1[k] = (TH1D*)R_RCP_TH1D_Array[k]->Add(R_RCP_TH1D_Array[k+72]);
					//Working_Array_2[k] = (TH1D*)R_LCP_TH1D_Array[k]->Add(R_LCP_TH1D_Array[k+72]);
					Working_Array_1[k] = R_RCP_TH1D_Array[k];
					Working_Array_2[k] = R_LCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk R_R and it is  "<< R_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else if (j == 1) {  // S enantiomer, different ellicity
					std::sprintf(tempk, "%s", S_RCP_TH1D_Array[k]->GetTitle());
					std::sprintf(tempj, "%s", S_LCP_TH1D_Array[k]->GetTitle());
					//Working_Array_1[k] = (TH1D*)S_RCP_TH1D_Array[k]->Add(S_RCP_TH1D_Array[k+72]);
					//Working_Array_2[k] = (TH1D*)S_LCP_TH1D_Array[k]->Add(S_LCP_TH1D_Array[k+72]);
					Working_Array_1[k] = S_RCP_TH1D_Array[k];
					Working_Array_2[k] = S_LCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be S_L and it is  "<< S_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else if (j == 2){  // DIFFERENT ellicity DIFFERENT enantiomer
					std::sprintf(tempj, "%s", R_LCP_TH1D_Array[k]->GetTitle());
					std::sprintf(tempk, "%s", S_RCP_TH1D_Array[k]->GetTitle());
					//Working_Array_1[k] = (TH1D*)R_LCP_TH1D_Array[k]->Add(R_LCP_TH1D_Array[k+72]);
					//Working_Array_2[k] = (TH1D*)S_RCP_TH1D_Array[k]->Add(S_RCP_TH1D_Array[k+72]);
					Working_Array_1[k] = R_LCP_TH1D_Array[k];
					Working_Array_2[k] = S_RCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else {
					std::cout << "j >= 3 there is a problem with the counter " << rcounter << endl;
				}

				// tempj and tempk taken from Working_Array_2
				std::sprintf(tempjk, "PECD %s vs. %s - cos(theta) = %2.2f - phi = %3.0f", tempj, tempk, cos_theta[k], phi[k]);
				TPaveStats *sk = (TPaveStats*)Working_Array_1[k]->GetListOfFunctions()->FindObject("stats");
				sk->SetOptStat(11);

				ratio_genarray[rcounter] = ratio_TH1D(Working_Array_1[k], Working_Array_2[k]);

				//ratio_genarray[rcounter]->Fit("PECD","QR"); // console output could e annoying: option "Q" to suppress it
				//fit1[k] = ratio_genarray[rcounter]->GetFunction("PECD");

				ratio_genarray[rcounter]->Fit("PECD6","QR"); // console output could e annoying: option "Q" to suppress it
				fit1[k] = ratio_genarray[rcounter]->GetFunction("PECD6");

				b1[rcounter] = fit1[k]->GetParameter(0);
				b1_err[rcounter] = fit1[k]->GetParError(0);

				//fitfunction(ratio_genarray[rcounter], fit1[k], b1[rcounter], b1_err[rcounter]);

				ratio_genarray[rcounter]->SetMarkerStyle(20);
				ratio_genarray[rcounter]->GetXaxis()->SetTitle("cos(theta)");
				ratio_genarray[rcounter]->GetXaxis()->SetRangeUser(-1.,1.);
				ratio_genarray[rcounter]->GetYaxis()->SetTitle("asymm [a.u.]");
				ratio_genarray[rcounter]->GetYaxis()->SetRangeUser(0.4, 1.6); // centered at 1
				ratio_genarray[rcounter]->SetTitle("");
				//const char* title = ("B1 = " + b1[i]); // example: ratio_genarray[rcounter]->SetTitle((const char*)histlist1D[j]);
				//ratio_genarray[rcounter]->SetTitle("B1_correct_title_char");

				canv[j]->cd(k+1);
				ratio_genarray[rcounter]->Draw("AP");
				fit1[k]->Draw("SAME");

				DrawTextInPad(0.2, 0.93, tempjk, 0.05);
				gStyle->SetPalette(200,MyPalette);
				gStyle->SetOptTitle(0);
				gStyle->SetOptFit();

				std::sprintf(tempjk, "b1fit_cos(theta)=%2.2f_phi=%3.0f_counter=%d", cos_theta[k], phi[k], rcounter);
				//cout << "check tempjk = " << tempjk << endl;
				savegraphas(tempjk,ratio_genarray[rcounter],fit1[k]);

				//// Set the leggend
				//TLegend *leg = canv[i]->cd(counter)->BuildLegend();
				//leg->SetTextFont(30);
				//leg->SetY1NDC(2);

				//cout << "ratio counter before = " << rcounter << endl;
				rcounter++;
				//cout << "ratio counter after = " << rcounter << endl;

				if (fp!=NULL) {
					// 12x6
					fprintf(fp,"%4.2f \t %4.0f \t %4.8f \t %4.8f \n", cos_theta[k], phi[k], b1[k+j*72], b1_err[k+j*72]);
					if (k == 71) {
						fprintf(fp, "\n", cos_theta[k], phi[k], b1[k + j * 72], b1_err[k + j * 72]);
						printf("saving b1 and b1_err file for each cos(theta) phi: DONE! \n");
					}
					//// 4x4
					//fstd::printf(fp,"%4.2f \t %4.0f \t %4.8f \t %4.8f \n", cos_theta[i], phi[i], b1[i+j*16], b1_err[i+j*16]);
					//// 3x3
					//fstd::printf(fp,"%4.2f \t %4.0f \t %4.8f \t %4.8f \n", cos_theta[i], phi[i], b1[i+j*9], b1_err[i+j*9]);
					//output a file with b1_parameter and error in two column!
				} else {
					printf("the b1 map .dat is not correctly opened \n");
					return 0;
				}
			}// end of for k
			canv[j]->Update();
			//canv[j]->Write(0,TObject::kWriteDelete); // you need to write the canvas to file
			canv[j]->Write(); // you need to write the canvas to file
			canv[j]->Close();
			gSystem->ProcessEvents();
			//std::cout << "plotting and saving cos(theta) canvas " << cname[j] << " DONE! \n" << endl;
			printf("plotting and saving cos(theta) canvas %s DONE! \n", cname[j]);
		} else {
			std::sprintf(tempjk, "%s_12x6", cname[j]);
			canv[j] = new TCanvas(cname[j], cname[j]);

			canv[j]->SetBatch(kFALSE);  // kTRUE to suppress the graphical putput
										// here kFALSE to display the b1maps
										//#if NOPRINTOUTLOAD == 1
										//			canv[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput	
										//#endif
										//canv[j]->SetCanvasSize(1920, 1080);
										//canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[j]->SetFillStyle(4000);
			canv[j]->cd();

			// 12x6
			b1_map[j-3] = new TGraph2D();
			b1_map[j-3]->SetNpx(12);
			b1_map[j-3]->SetNpy(6);			
			b1_map[j-3]->SetTitle(tempjk); // is is necessary?
			canv[j]->Update();

			////filling b1_maps
			// 12x6
			for (int k = 0; k < 72; k++) //k = all phi and cos(theta) = 12x6
				b1_map[j-3]->SetPoint(k, phi[k], cos_theta[k], b1[k+(j-3)*72]);

			b1_map[j-3]->GetHistogram();  // By default returns a pointer to the Delaunay histogram.
										  // the ffollwing settings should be preceed from ->GetHistogram()
										  //b1_map_copy[j-3]=(TH2D*)b1_map[j-3];
										  //b1_map[j-3]->SetTitle("  ; phi [DEG]; cos(theta) [adm]"); SetTitle("global title;x axis;y axis")
										  //b1_map[j-3]->GetXaxis()->SetRangeUser(-180., 180.); 
			b1_map[j-3]->GetXaxis()->SetTitle("cos(theta) [adm]");
			b1_map[j-3]->GetXaxis()->CenterTitle(true);
			//b1_map[j-3]->GetYaxis()->SetRangeUser(-1., 1.);	    //
			b1_map[j-3]->GetYaxis()->SetTitle("b1 [adm]");
			b1_map[j-3]->GetYaxis()->CenterTitle(true);
			//b1_map[j-3]->GetZaxis()->ResetAttAxis();
			b1_map[j-3]->GetZaxis()->SetRangeUser(-0.8, 0.8);
			//b1_map[j-3]->GetZaxis()->SetLimits(-0.6, 0.6);

			b1_map[j-3]->Draw("COLZ"); // this will introduce directly the smoothing
									   //b1_map[j-3]->Draw("LINE");
									   //Save as the b1_maps
			canv[j]->Update();
			// final touches: title and legend
			//DrawTextInPad(0.28, 0.93, tempjk, 0.03);  // double work
			//TLegend *leg = canv[j]->cd()->BuildLegend();
			//leg->SetTextFont(30);
			//leg->SetY1NDC(2);

			std::sprintf(tempjk, "%s.png", tempjk);
			canv[j]->SaveAs(tempjk);
			//canv[j]->Update();
			b1_map[j-3]->GetZaxis()->SetRangeUser(-0.8, 0.8);
			canv[j]->Write(); //for the b1 maps
			printf("plotting and saving b1_maps canvas %s: DONE! \n", cname[j]);
		}
	}


	for (int j = 0; j < 3; j++) { // MFPDAS smoothed asymmetries

								  //Canvas setup
								  //gStyle->SetPadTopMargin(0.);
								  //gStyle->SetPadBottomMargin(0.);
								  //gStyle->SetPadLeftMargin(0.);
								  //gStyle->SetPadRightMargin(0.);

								  /*Double_t mlb = 0.3;
								  Double_t mrt = 0.1;
								  Double_t nx  = 12;
								  Double_t ny  = 6;
								  DividegPad(nx,ny,mlb,mrt,mrt,mlb);*/

		canv2D[j] = new TCanvas(cname2D[j], cname2D[j]);
		canv2D[j]->SetCanvasSize(1920, 1080);
		canv2D[j]->SetWindowSize(1920 + 4, 1080 + 28);
		canv2D[j]->Divide(12, 6);
		canv2D[j]->SetFillStyle(4000);
		canv2D[j]->SetFrameBorderMode(0);
		canv2D[j]->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
									 //#if NOPRINTOUTLOAD2D == 1
									 //		canv2D[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
									 //#endif

		canv2D_smooth[j] = new TCanvas(cname2D[j+3], cname2D[j+3]);
		canv2D_smooth[j]->SetCanvasSize(1920, 1080);
		canv2D_smooth[j]->SetWindowSize(1920 + 4, 1080 + 28);
		canv2D_smooth[j]->Divide(12, 6);
		canv2D_smooth[j]->SetFillStyle(4000);
		canv2D_smooth[j]->SetFrameBorderMode(0);

		canv2D_smooth[j]->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
											//#if NOPRINTOUTLOAD2D == 1
											//		canv2D_smooth[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
											//#endif

		std::sprintf(tempj, "FitSlicesY b1 - canvas %s", cname2D[j+3]);
		canv_fit[j] = new TCanvas(tempj, tempj);
		canv_fit[j]->SetCanvasSize(1920, 1080);
		canv_fit[j]->SetWindowSize(1920 + 4, 1080 + 28);
		canv_fit[j]->Divide(12, 6);
		canv_fit[j]->SetFillStyle(4000);

		canv_fit[j]->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOPRINTOUTLOAD2D == 1
		canv_fit[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
		gStyle->SetPalette(200,MyPalette);
		gStyle->SetOptTitle(0);
		gSystem->ProcessEvents();

		int rcounter_2D = 0; // ?????
		int counter_2D = 72;

		for (int k = 0; k < 72; k++) {
			std::sprintf(tempj, "clean");
			std::sprintf(tempk, "clean");
			if (j == 0) { // R enantiomer, different ellicity
				std::sprintf(tempk, "%s", R_RCP_TH2D_Array[k]->GetName());
				std::sprintf(tempj, "%s", R_LCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = R_RCP_TH2D_Array[k];
				Working_Array2D_2[k] = R_LCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk R_R and it is  "<< R_RCP_TH1D_Array[i]->GetTitle() << endl;
			}
			else if (j == 1) {  // S enantiomer, different ellicity
				std::sprintf(tempk, "%s", S_RCP_TH2D_Array[k]->GetName());
				std::sprintf(tempj, "%s", S_LCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = S_RCP_TH2D_Array[k];
				Working_Array2D_2[k] = S_LCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be S_L and it is  "<< S_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl;
			}
			else if (j == 2) {  // DIFFERENT ellicity DIFFERENT enantiomer
				std::sprintf(tempj, "%s", R_LCP_TH2D_Array[k]->GetName());
				std::sprintf(tempk, "%s", S_RCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = R_LCP_TH2D_Array[k];
				Working_Array2D_2[k] = S_RCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl;
			}else{
				std::cout << "j exceeds 3!" << endl;
				return 0;
			}

			canv2D[j]->cd(counter_2D);
			std::sprintf(tempjk, "ratio=%d_k=%d_j=%d",rcounter_2D, k, j);
			asymmTH2D[rcounter_2D] = asymm_TH2D(Working_Array2D_1[k], Working_Array2D_2[k], cname2D[j]);
			//cout << "Number of entries for ratio TH2D is " << asymmTH2D[rcounter_2D]->GetEntries() << endl;
			asymmTH2D[rcounter_2D]->SetTitle("");
			asymmTH2D[rcounter_2D]->SetStats(0);
			//TPaveStats *sk = (TPaveStats*)Working_Array2D_1[k]->GetListOfFunctions()->FindObject("stats");
			//sk->SetOptStat(11);
			asymmTH2D[rcounter_2D]->GetZaxis()->SetRangeUser(-0.2, 0.2); // doesn't work in printing
			asymmTH2D[rcounter_2D]->Draw();  //avoid COLZ; CONTZ smooths but doesn't work for this configuration
			canv2D[j]->Update();
			DrawTextInPad(0.2, 0.93, tempjk, 0.06);
			canv2D[j]->Update();
			gSystem->ProcessEvents();

			canv2D_smooth[j]->Write(canv2D_smooth[j]->GetName(),TObject::kWriteDelete);

			//gPad->Modified()
			//gPad->Update();
			//gSystem->ProcessEvents();

			// smoothing
			std::sprintf(tempjk, "ratio=%d_k=%d_j=%d_smoothed",rcounter_2D, k, j);
			canv2D_smooth[j]->cd(counter_2D);
			//asymmTH2D[rcounter_2D]->GetZaxis()->SetLimits(-0.2, 0.2);
			//asymmTH2D[rcounter_2D]->GetZaxis()->SetRangeUser(-0.2, 0.2); // doesn't work in printing
			asymmTH2D[rcounter_2D]->Smooth(1,"k3a");  //Smooth() : 1 = n times, "k5a", "k5b" (5x5) and k3a (3x3)
			asymmTH2D[rcounter_2D]->Draw();  //avoid COLZ; CONTZ smooths but doesn't work for this configuration
			canv2D_smooth[j]->Update();
			DrawTextInPad(0.2, 0.93, tempjk, 0.06);
			canv2D_smooth[j]->Update();
			gSystem->ProcessEvents();
			//gPad->Modified()
			//gPad->Update();
			//gSystem->ProcessEvents();

			//Just the smoothed version is printed out
			//cout << " readinf file " << tempjk << endl;
			savegraphas(tempjk, asymmTH2D[rcounter_2D]);

			// FITSLICEY for b1 slice
			std::sprintf(tempjk, "canvas=%d_k=%d_j=%d",rcounter_2D, k, j);
			TObjArray *aSlice = new TObjArray();
			aSlice->SetOwner(kTRUE);

			//adaptentries(asymmTH2D[rcounter_2D]);
			asymmTH2D[rcounter_2D]->FitSlicesY(f1, 0, -1, -1, "QR", aSlice);  // function PECD fit 
																			  //asymmTH2D[rcounter_2D]->FitSlicesY(f2, 0, -1, -1, "QR", aSlice);	// function PECD4 fit
			copyaSliceTH1D[k] = (TH1D*)(*aSlice)[0]->Clone();
			//adaptentries_red(copyaSliceTH1D[k]);
			canv_fit[j]->cd(counter_2D);
			//canv_fit[j]->cd(counter_2D)->SetBottomMargin(0.05);
			//canv_fit[j]->cd(counter_2D)->SetTopMargin(0.05);
			//canv_fit[j]->cd(counter_2D)->SetLeftMargin(0.05);
			//canv_fit[j]->cd(counter_2D)->SetRightMargin(0.05);
			copyaSliceTH1D[k]->SetTitle(tempjk);
			copyaSliceTH1D[k]->SetName(tempjk);
			copyaSliceTH1D[k]->SetMarkerStyle(20);
			copyaSliceTH1D[k]->SetMarkerSize(1);
			copyaSliceTH1D[k]->SetStats(0);
			copyaSliceTH1D[k]->GetYaxis()->SetTitle("b1");
			//copyaSliceTH1D[k]->GetYaxis()->SetRangeUser(-0.2, 0.2);
			copyaSliceTH1D[k]->Draw("E1");
			//(*aSlice)[0]->Draw();  // draw b1 histogram
			canv_fit[j]->Update();
			gSystem->ProcessEvents(); //is it usefull?

			savegraphas(tempjk, copyaSliceTH1D[k]); // weird passing o aSlice

			rcounter_2D++;
			counter_2D--;

		}
		canv2D[j]->Update(); 
		gSystem->ProcessEvents(); //is it usefull?
		canv2D[j]->Write();
		canv2D[j]->Close();
		printf("plotting and saving TH2D smoothed canvas %s: DONE! \n \n", cname2D[j]);

		//canv2D_amooth[j]->Update(); 
		//gSystem->ProcessEvents(); //is it usefull?
		//canv2D_smooth[j]->Write(canv2D_smooth[j]->GetName(),TObject::kWriteDelete);
		canv2D_smooth[j]->Close();
		printf("plotting and saving TH2D smoothed canvas %s: DONE! \n \n", cname2D[j+3]);

		canv_fit[j]->Update();
		gSystem->ProcessEvents(); //is it usefull?
		canv_fit[j]->Write();
		canv_fit[j]->Close();
		printf("plotting and saving FitSlicesY smoothed canvas %s: DONE! \n \n", cname2D[j+3]);
	}

	f->Close();
	std::cout << "gPad : SAVED! \n " << endl;
	fclose(fp);
	// normal exit with no delete of variable in c++: exit(0) is functionally indentical to return 0
	exit(0);
	theApp.Run();
	//terminate();

	return 0;
}
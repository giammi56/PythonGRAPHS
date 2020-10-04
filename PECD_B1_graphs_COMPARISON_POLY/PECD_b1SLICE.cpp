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
#define NOPRINT 0
#define NOPRINTLOAD 0
#define NOOUTPUT 0
#define NOOUTPUTLOAD 0

using namespace std;


void DrawTextInPad(double xCoord, double yCoord, char* text, Float_t size) {
	TText* t1 = new TText(xCoord, yCoord, text);
	t1->SetNDC();
	t1->SetTextColor(1);
	t1->SetTextSize(size);
	t1->SetTextFont(62);
	t1->Draw();
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

void savegraphas(char* tempjk, TH2D* ratioTH2D) {

#if NOPRINT == 1
	return;
#endif

	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	sprintf(tempjk, "%s.png", tempjk);
	cdummy->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOOUTPUT == 1
	cdummy->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
	cdummy->cd();
	centrelabelaxis(ratioTH2D);
	ratioTH2D->Draw();
	cdummy->SaveAs(tempjk);
	if (cdummy) { 
		cdummy->Close();
		gSystem->ProcessEvents();
	}

	//std::cout << ratioTH2D->GetName() << ": Exported as .png \n" << endl << endl;

	return;
}

void savegraphas(char* tempjk, TGraphAsymmErrors* ratio_genarray, TF1* fit1) {

#if NOPRINTLOAD == 1
	return;
#endif

	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	//TLegend *leg = cdummy->cd()->BuildLegend();
	//leg->SetTextFont(30);
	//leg->SetY1NDC(2);
	sprintf(tempjk, "%s.png", tempjk);
#if NOOUTPUTLOAD == 1
	cdummy->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
	cdummy->cd();
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

#if NOPRINT == 1
	return;
#endif
	
	//ratioTH1D->SetTitle(tempjk);
	//ratioTH1D->SetTitle("");
	TCanvas* cdummy = new TCanvas(tempjk, tempjk);
	sprintf(tempjk, "%s.png", tempjk);
	cdummy->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOOUTPUT == 1
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

	//std::cout << ratioTH2D->GetName() << ": Exported as .png \n" << endl << endl;

	return;
}

//alternative method to define the fitting function
Double_t buildfunction(Double_t* x, Double_t* par) {
	Double_t arg = 0;
	arg = (x[0] - par[1])/par[2];
	Double_t f1new = par[0] * x[0] / (1 - 0.125*par[1] * x[0] * x[0]) + 1;
	return f1new;
}

TGraphAsymmErrors* ratio_TH1D(TH1D* &Working_Array_1, TH1D* &Working_Array_2, const char* &cname) {
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

TH2D* ratio_TH2D(TH2D* &Working_Array2D_1, TH2D* &Working_Array2D_2, const char* &cname) {
	TH2D* ratioTH2D = (TH2D*)Working_Array2D_1->GetAsymmetry(Working_Array2D_2);
	//cout << "Total number of entries = " << ratioTH2D->GetEntries() << endl;
	//cout << "Size of ratioTH2D: " << (ratioTH2D->Sizeof())/sizeof(int) << endl << endl;
	//cout << "Division TH2D DONE!" << endl << endl;

	return ratioTH2D;
}

TH2D* adaptentries(TH2D* &ratioTH2D) {
	
	//higer entries to allow FitSliceY()
	int nx = ratioTH2D->GetXaxis()->GetNbins();
	int ny = ratioTH2D->GetYaxis()->GetNbins();
	for(int x=0; x<nx; x++){
		for(int y=0; y<ny; y++){
			Double_t temp = ratioTH2D->GetBinContent(x, y);
			ratioTH2D->SetBinContent(x,y,temp * 10000);
			//cout << "New bin content: old = " << temp << " new = " << ratioTH2D->GetBinContent(x, y) << endl;
		}
	}

	//cout << "Total number of entries = " << ratioTH2D->GetEntries() << endl;
	//cout << "Size of ratioTH2D: " << (ratioTH2D->Sizeof())/sizeof(int) << endl << endl;
	//cout << "Division TH2D DONE!" << endl << endl;

	return ratioTH2D;
}

TH1D* adaptentries_red(TH1D* &ratioTH1D) {

	int nx = ratioTH1D->GetXaxis()->GetNbins();
		for(int x=0; x<nx; x++){
			Double_t temp = ratioTH1D->GetBinContent(x);
			Double_t err = ratioTH1D->GetBinError(x);
			ratioTH1D->SetBinError(x,err / 100); // sqrt(10000)
			ratioTH1D->SetBinContent(x,temp / 10000);
			//cout << "New bin content: old = " << temp << " new = " << ratioTH2D->GetBinContent(x, y) << endl;
		}
	//cout << "Total number of entries = " << ratioTH2D->GetEntries() << endl;
	//cout << "Size of ratioTH2D: " << (ratioTH2D->Sizeof())/sizeof(int) << endl << endl;
	//cout << "Division TH2D DONE!" << endl << endl;

	return ratioTH1D;
}

// TAKE CARE THAT THE FEEDING PARAMETERS ARE IN TEH SAME ORDER!!!
void load_TH1D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH1D** S_LCP_TH1D_Array, TH1D** S_RCP_TH1D_Array, TH1D** R_LCP_TH1D_Array, TH1D** R_RCP_TH1D_Array, TH1D** general_Array, TString* dirlist, TString* histlist1D, const char dir[]) {
	int counter = 0;
	int i = 0;
	//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
	TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory(dir+dirlist[i]));
	TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory(dir+dirlist[i]));
	TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory(dir+dirlist[i]));
	TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory(dir+dirlist[i]));
	gDirectory->pwd();
	//for (int j = 0; j < LCPhistlist->Sizeof(); j++) {

	//12x6
	for (int k=0;k<72;k++) {
		////4x4
		//for (int k=0;k<16;k++) {
		////3x3
		//for (int k=0;k<9;k++) {
		S_RCP_TH1D_Array[k]=(TH1D*)S_Rdir->FindObjectAny(histlist1D[k]);				 
		S_RCP_TH1D_Array[k]->SetName("S_RCP");								 
		S_RCP_TH1D_Array[k]->SetTitle("S_RCP");								 
		S_RCP_TH1D_Array[k]->Sumw2();													 
		S_RCP_TH1D_Array[k]->Scale(1. / S_RCP_TH1D_Array[k]->Integral());			 
		//std::cout << "Loading and normalisation of " << S_RCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << S_RCP_TH1D_Array[k]->GetEntries() << endl << endl;

		general_Array[counter++] = S_RCP_TH1D_Array[k];			

		S_LCP_TH1D_Array[k]=(TH1D*)S_Ldir->FindObjectAny(histlist1D[k]);
		S_LCP_TH1D_Array[k]->SetName("S_LCP");
		S_LCP_TH1D_Array[k]->SetTitle("S_LCP");
		S_LCP_TH1D_Array[k]->Sumw2();
		S_LCP_TH1D_Array[k]->Scale(1. / S_LCP_TH1D_Array[k]->Integral());
		//std::cout << "Loading and normalisation of " << S_LCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << S_LCP_TH1D_Array[k]->GetEntries() << endl << endl;

		general_Array[counter++] = S_LCP_TH1D_Array[k];							 

		R_RCP_TH1D_Array[k] = (TH1D*)R_Rdir->FindObjectAny(histlist1D[k]);				 
		R_RCP_TH1D_Array[k]->SetName("R_RCP");								 
		R_RCP_TH1D_Array[k]->SetTitle("R_RCP");								 
		R_RCP_TH1D_Array[k]->Sumw2();													 
		R_RCP_TH1D_Array[k]->Scale(1. / R_RCP_TH1D_Array[k]->Integral());			 
		//std::cout << "Loading and normalisation of " << R_RCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << R_RCP_TH1D_Array[k]->GetEntries() << endl << endl;

		general_Array[counter++] = R_RCP_TH1D_Array[k];

		R_LCP_TH1D_Array[k] = (TH1D*)R_Ldir->FindObjectAny(histlist1D[k]);				 
		R_LCP_TH1D_Array[k]->SetName("R_LCP");								 
		R_LCP_TH1D_Array[k]->SetTitle("R_LCP");								 
		R_LCP_TH1D_Array[k]->Sumw2();													 
		R_LCP_TH1D_Array[k]->Scale(1. / R_LCP_TH1D_Array[k]->Integral());			 
		//std::cout << "Loading and normalisation of " << R_LCP_TH1D_Array[k]->GetName() << "/" << histlist1D[k] << " done! # entries = " << R_LCP_TH1D_Array[k]->GetEntries() << endl << endl;

		general_Array[counter++] = R_LCP_TH1D_Array[k];
	}

	std::cout << "Total number of TH1D histos = " << counter << endl << endl;

	//Autoput of all loaded graphs
	const char *cname[4] = { "S_RCP","S_LCP","R_RCP", "R_LCP" };

	//plotting all cos(theta) plots
	for (int i = 0; i < 4; i++) {
		char canvastitle[100];
		sprintf(canvastitle, "all entries = %d", counter / 4);
		char canvasname[100];
		sprintf(canvasname, "All %s plots", cname[i]);
		TCanvas* cplot = new TCanvas(canvasname, canvastitle);
		cplot->SetCanvasSize(1920, 1080);
		cplot->SetWindowSize(1920 + 4, 1080 + 28);
		cplot->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOOUTPUTLOAD == 1
		cplot->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif 
								//12x6
		cplot->Divide(12, 6);
		for (int k = 0; k < 72; k++) {
			//// 4x4
			//cplot->Divide(16,4);
			//for (int k=0; k<64; k++) {	
			//// 3x3
			//cplot->Divide(9,4);
			//for (int k=0; k<36; k++) {
			cplot->cd(k + 1);
			//general_Array[k]->SetOption("E");
			general_Array[k * 4 + i]->Draw("E1");
			//gPad->Modified();
			gPad->Update();
		}
		cplot->Write();
		cplot->Close();
		gSystem->ProcessEvents();
	}

	std::cout << endl << endl;
	return;
}

void load_TH2D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH2D** S_LCP_TH2D_Array, TH2D** S_RCP_TH2D_Array, TH2D** R_LCP_TH2D_Array, TH2D** R_RCP_TH2D_Array, TH2D** general_Array2D, TString* dirlist, TString* histlist2D, const char dir[]) {

	int counter = 0;
	int i = 0;
	//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
	TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory(dir+dirlist[i]));
	TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory(dir+dirlist[i]));
	TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory(dir+dirlist[i]));
	TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory(dir+dirlist[i]));

	gDirectory->pwd();

	for (int k=0;k<72;k++) {
		S_RCP_TH2D_Array[k]=(TH2D*)S_Rdir->FindObjectAny(histlist2D[k]);				 
		S_RCP_TH2D_Array[k]->SetName("S_RCP"+dirlist[i]);								 
		S_RCP_TH2D_Array[k]->SetTitle("S_RCP"+dirlist[i]);								 
		//S_RCP_TH2D_Array[k+i*3]->Sumw2();													 
		normalize_TH2D(S_RCP_TH2D_Array[k]);
		general_Array2D[counter++] = S_RCP_TH2D_Array[k];

		S_LCP_TH2D_Array[k]=(TH2D*)S_Ldir->FindObjectAny(histlist2D[k]);
		S_LCP_TH2D_Array[k]->SetName("S_LCP"+dirlist[i]);
		S_LCP_TH2D_Array[k]->SetTitle("S_LCP"+dirlist[i]);
		//S_LCP_TH2D_Array[k+i*3]->Sumw2();
		normalize_TH2D(S_LCP_TH2D_Array[k]);
		general_Array2D[counter++] = S_LCP_TH2D_Array[k];										 

		R_RCP_TH2D_Array[k] = (TH2D*)R_Rdir->FindObjectAny(histlist2D[k]);				 
		R_RCP_TH2D_Array[k]->SetName("R_RCP"+dirlist[i]);								 
		R_RCP_TH2D_Array[k]->SetTitle("R_RCP"+dirlist[i]);								 
		//R_RCP_TH2D_Array[k+i*3]->Sumw2();													 
		normalize_TH2D(R_RCP_TH2D_Array[k]);
		general_Array2D[counter++] = R_RCP_TH2D_Array[k];

		R_LCP_TH2D_Array[k] = (TH2D*)R_Ldir->FindObjectAny(histlist2D[k]);				 
		R_LCP_TH2D_Array[k]->SetName("R_LCP"+dirlist[i]);								 
		R_LCP_TH2D_Array[k]->SetTitle("R_LCP"+dirlist[i]);								 
		//R_LCP_TH2D_Array[k+i*3]->Sumw2();													 
		normalize_TH2D(R_LCP_TH2D_Array[k]);
		general_Array2D[counter++] = R_LCP_TH2D_Array[k];
	}

	std::cout << "Total number of TH2D histos = " << counter << endl << endl;

	const char *cname[4] = { "S_RCP","S_LCP","R_RCP", "R_LCP" };

	for (int j = 0; j < 4; j++) { //0S_CR, 1S_CL, 2R_CR, 3R_CL
		char canvastitle_2D[100];
		char canvasname[100];
		sprintf(canvasname, "All %s plots", cname[i]);
		sprintf(canvastitle_2D, "MFPADs all entries = %d", counter / 8);
		TCanvas* cplot = new TCanvas(canvasname, canvastitle_2D);
		cplot->SetCanvasSize(1920, 1080);
		cplot->SetWindowSize(1920 + 4, 1080 + 28);
		cplot->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
#if NOOUTPUTLOAD == 1
		cplot->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
		cplot->Divide(12, 6);
		for (int k = 0; k < 72; k++) {
			cplot->cd(k + 1);
			//general_Array[k]->SetOption("E");
			general_Array2D[k * 4 * (i + 1) + j]->Draw("COLZ");
			//gPad->Modified();
			gPad->Update();
		}
		cplot->Write();
		cplot->Close();
		gSystem->ProcessEvents();
	}

	std::cout << endl << endl;

	return;

}

int main(int argc, char** argv)
//int main(int argc, char* argv[])

{
	TApplication theApp("App", &argc, argv);

	// single TH1D arrays: adapt to the requested dimension i*j
	int dim = 288; //12x6x4
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

	// TH2D arrays!!!
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

	TH2D** ratioTH2D;
	ratioTH2D = new TH2D*[dim];

	TH1D** copyaSliceTH1D;
	copyaSliceTH1D = new TH1D*[dim];

	// Fititign parameters
	TF1** fit1;
	fit1 = new TF1*[];
	//Alternative description of fitting formula
	Double_t* x;
	x = new Double_t[];
	Double_t* par;
	par = new Double_t[];

	//Adapt dimensio to request
	Double_t b1[288];
	Double_t b1_err[288];

	double cos_theta[288/4];
	double phi[288/4];

	int ratiocounter;
	int counter;

	TGraph2D** b1_map;
	b1_map = new TGraph2D*[3];

	Int_t MyPalette[200];
	Double_t r[]    = {0.0, 1.0, 1.0, 1.0};
	Double_t g[]    = {0.0, 1.0, 1.0, 0.0};
	Double_t b[]    = {1.0, 1.0, 0.0, 0.0};
	Double_t stop[] = {0.0, 0.5, 0.75, 1.0};
	Int_t FI = TColor::CreateGradientColorTable(4, stop, r, g, b, 200);
	for (int i=0;i<200;i++) MyPalette[i] = FI+i;
	gStyle->SetPalette(200,MyPalette);

	// DO NOT CHANGE THE ORDER OF THESE INPUTS R-C3H3F3O_546eV_CL_10800-2350ns_multiCH11_MFPAD_30.root
	//550 eV
	TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	if (S_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_LFile->GetName());
	TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	if (S_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_RFile->GetName());
	TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	if (R_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_LFile->GetName());
	TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30_t0.root", "READ");
	if (R_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_RFile->GetName());

	//550eV
	//TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (S_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_LFile->GetName());
	//TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (S_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_RFile->GetName());
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (R_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (R_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_RFile->GetName());

	//// DO NOT CHANGE THE ORDER OF THESE INPUTS
	//TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/S-C3H3F3O_poly_30-69m_550eV_CL_reaction_POLY_3x3_test.root", "READ");
	//if (S_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_LFile->GetName());
	//TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/S-C3H3F3O_poly_30-69m_550eV_CR_reaction_POLY_3x3_test.root", "READ");
	//if (S_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", S_RFile->GetName());
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CL_reaction_POLY_3x3_test.root", "READ");
	//if (R_LFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CR_reaction_POLY_3x3_test.root", "READ");
	//if (R_RFile->IsOpen()) std::printf("%s: File opened successfully!\n", R_RFile->GetName());

	//File to save outout b1
	FILE *fp = fopen("b1_550eV_ID_ALL_CL_CH11m69_12x8_TEST.dat","w+");
	//FILE *fp = fopen("b1_ID_ALL_m69_3x3_test.dat","w");
	fprintf(fp,"cos(theta) \t phi \t p1 \t p1_err \n");

	//FILE NAME
	// POLYATOMIC - THREE CHARGED FRAGMENTS
	//TFile *f = new TFile("FMtox_ID_ALL_m69_3x3_test.root","RECREATE");
	//TFile *f = new TFile("FMtox_ID3_m41_69.root", "RECREATE");
	//POLIATOMIC - NEUTRAL
	TFile *f = new TFile("FMtox_550eV_ID_ALL_CL_CH11m69_12x6_TEST.root","RECREATE");

	if (f->IsOpen()) {
		std::printf("File CREATED successfully\n");
	}

	const char dir[] = {"angular_distr_el/ID_ALL_mol_e0_valid"};
	//TString dirlist[1] = {"/EN_gate/MFPADs_CR"};
	TString dirlist[1] = {"/EN_gate/MFPADs_CL"};
	//12x6
	TString histlist1D[72] = { "cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-0.67_phi_-180","cos(theta)_e[0]_costheta_-0.33_phi_-180","cos(theta)_e[0]_costheta_-0.00_phi_-180","cos(theta)_e[0]_costheta_0.33_phi_-180","cos(theta)_e[0]_costheta_0.67_phi_-180",
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
	//12x6
	TString histlist2D[72] = { "MFPAD3D_engate_costheta_-1.00_phi_-180", "MFPAD3D_engate_costheta_-0.67_phi_-180","MFPAD3D_engate_costheta_-0.33_phi_-180","MFPAD3D_engate_costheta_-0.00_phi_-180", "MFPAD3D_engate_costheta_0.33_phi_-180", "MFPAD3D_engate_costheta_-0.67_phi_-180",
		"MFPAD3D_engate_costheta_-1.00_phi_-150", "MFPAD3D_engate_costheta_-0.67_phi_-150","MFPAD3D_engate_costheta_-0.33_phi_-150","MFPAD3D_engate_costheta_-0.00_phi_-150", "MFPAD3D_engate_costheta_0.33_phi_-150", "MFPAD3D_engate_costheta_-0.67_phi_-150",
		"MFPAD3D_engate_costheta_-1.00_phi_-120", "MFPAD3D_engate_costheta_-0.67_phi_-120","MFPAD3D_engate_costheta_-0.33_phi_-120","MFPAD3D_engate_costheta_-0.00_phi_-120", "MFPAD3D_engate_costheta_0.33_phi_-120", "MFPAD3D_engate_costheta_-0.67_phi_-120",
		"MFPAD3D_engate_costheta_-1.00_phi_-90",	 "MFPAD3D_engate_costheta_-0.67_phi_-90", "MFPAD3D_engate_costheta_-0.33_phi_-90", "MFPAD3D_engate_costheta_-0.00_phi_-90",  "MFPAD3D_engate_costheta_0.33_phi_-90",  "MFPAD3D_engate_costheta_-0.67_phi_-90",
		"MFPAD3D_engate_costheta_-1.00_phi_-60",	 "MFPAD3D_engate_costheta_-0.67_phi_-60", "MFPAD3D_engate_costheta_-0.33_phi_-60", "MFPAD3D_engate_costheta_-0.00_phi_-60",  "MFPAD3D_engate_costheta_0.33_phi_-60",  "MFPAD3D_engate_costheta_-0.67_phi_-60",
		"MFPAD3D_engate_costheta_-1.00_phi_-30",	 "MFPAD3D_engate_costheta_-0.67_phi_-30", "MFPAD3D_engate_costheta_-0.33_phi_-30", "MFPAD3D_engate_costheta_-0.00_phi_-30",  "MFPAD3D_engate_costheta_0.33_phi_-30",  "MFPAD3D_engate_costheta_-0.67_phi_-30",
		"MFPAD3D_engate_costheta_-1.00_phi_0",	 "MFPAD3D_engate_costheta_-0.67_phi_0",   "MFPAD3D_engate_costheta_-0.33_phi_0",   "MFPAD3D_engate_costheta_-0.00_phi_0",    "MFPAD3D_engate_costheta_0.33_phi_0",    "MFPAD3D_engate_costheta_-0.67_phi_0",
		"MFPAD3D_engate_costheta_-1.00_phi_30",	 "MFPAD3D_engate_costheta_-0.67_phi_30",  "MFPAD3D_engate_costheta_-0.33_phi_30",  "MFPAD3D_engate_costheta_-0.00_phi_30",   "MFPAD3D_engate_costheta_0.33_phi_30",   "MFPAD3D_engate_costheta_-0.67_phi_30",
		"MFPAD3D_engate_costheta_-1.00_phi_60",	 "MFPAD3D_engate_costheta_-0.67_phi_60",  "MFPAD3D_engate_costheta_-0.33_phi_60",  "MFPAD3D_engate_costheta_-0.00_phi_60",   "MFPAD3D_engate_costheta_0.33_phi_60",   "MFPAD3D_engate_costheta_-0.67_phi_60",
		"MFPAD3D_engate_costheta_-1.00_phi_90",	 "MFPAD3D_engate_costheta_-0.67_phi_90",  "MFPAD3D_engate_costheta_-0.33_phi_90",  "MFPAD3D_engate_costheta_-0.00_phi_90",   "MFPAD3D_engate_costheta_0.33_phi_90",   "MFPAD3D_engate_costheta_-0.67_phi_90",
		"MFPAD3D_engate_costheta_-1.00_phi_120",	 "MFPAD3D_engate_costheta_-0.67_phi_120", "MFPAD3D_engate_costheta_-0.33_phi_120", "MFPAD3D_engate_costheta_-0.00_phi_120",  "MFPAD3D_engate_costheta_0.33_phi_120",  "MFPAD3D_engate_costheta_-0.67_phi_120",
		"MFPAD3D_engate_costheta_-1.00_phi_150",	 "MFPAD3D_engate_costheta_-0.67_phi_150", "MFPAD3D_engate_costheta_-0.33_phi_150", "MFPAD3D_engate_costheta_-0.00_phi_150",  "MFPAD3D_engate_costheta_0.33_phi_150",  "MFPAD3D_engate_costheta_-0.67_phi_150"};



	// canvas strig name, from Till suggestion! Not necessary now.
	const char *cname[] = {"PECD_e[0]_R_CR-CL","PECD_e[0]_S_CR-CL","PECD_e[0]_R_CL-S_CR", "b1_map_R_CR-CL","b1_map_S_CR-CL","b1_map_R_CL-S_CR"};
	const char *cname2D[] = {"MFPAD_e[0]_R_CR-CL","MFPAD_e[0]_S_CR-CL","MFPAD_e[0]_R-S_CL-CR"};
	TCanvas *canv[6]; // dim(3 cos(theta) + 3 b1__maps)
	TCanvas *canv2D[3]; // dim(3 MFPAD)
	TCanvas *canv_fit[3];
	std::cout << "Vector of TCanvas created. \n" << endl << endl;

	char tempj[250];
	char tempk[250];
	char tempjk[250];

	// LOAD AND NORMALIZED ALL GRAPHS
	load_TH1D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH1D_Array, S_RCP_TH1D_Array, R_LCP_TH1D_Array, R_RCP_TH1D_Array, general_Array, dirlist, histlist1D, dir);
	std::cout << "load_all TH1D : DONE! \n" << endl << endl;

	load_TH2D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH2D_Array, S_RCP_TH2D_Array, R_LCP_TH2D_Array, R_RCP_TH2D_Array, general_Array2D, dirlist, histlist2D, dir);
	std::cout << "load_all TH2D : DONE! \n" << endl << endl;

	//Suppress INFO messages ROOT
	//gErrorIgnoreLevel = kWarning;

	ratiocounter = 0;
	//generation of phi and cos(theta) 16 bins: be sure that the increment and the incremented variabel is the same at in histlist1D 
	counter = 0;
	// 12x6 : 12 cos(theta), 6 phi
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 6; j++) {
			phi[counter] = i*30-180+15;
			cos_theta[counter] = j*0.333-1+0.1667; // (bin-period/2+bin/2)
			counter++;
		}
	}

	//// 4x4
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		cos_theta[counter] = i*0.5-1+0.25; // (bin-period/2+bin/2)
	//		phi[counter] = j*90-180+45;
	//		counter++;
	//	}
	//}

	// 3x3
	//for (int i = 0; i < 3; i++) {
	//	for (int j = 0; j < 3; j++) {
	//		cos_theta[counter] = i*0.667-1+0.33; // (bin-period/2+bin/2)
	//		phi[counter] = j*120-180+60;
	//		counter++;
	//	}
	//}

	//3 canvas of i*j plots each cos(theta) + 3 cavas for b1_map
	for (int j = 0; j < 6; j++) {
		if (j < 3) {
			canv[j] = new TCanvas(cname[j], cname[j]);
			canv[j]->SetCanvasSize(1920, 1080);
			canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
#if NOOUTPUTLOAD == 1
			canv[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif
									  //12x6
			canv[j]->Divide(12, 6);						  
			// 9x8
			//canv[j]->Divide(9, 8);
			//// 4x4
			//canv[j]->Divide(4, 4);
			//// 3x3
			//canv[j]->Divide(3, 3);
			canv[j]->SetFillStyle(4000);

			// all cos(theta) -> counter goes for ixj entries
			//12X6
			for (int k = 0; k < 72; k++) {
				//4x4
				//for (int i = 0; i < 16; i++) {
				//3X3
				//for (int i = 0; i < 9; i++) {
				sprintf(tempj, "clean");
				sprintf(tempk, "clean");
				if (j == 0) { // R enantiomer, different ellicity
					sprintf(tempk, "%s", R_RCP_TH1D_Array[k]->GetTitle());
					sprintf(tempj, "%s", R_LCP_TH1D_Array[k]->GetTitle());
					Working_Array_1[k] = R_RCP_TH1D_Array[k];
					Working_Array_2[k] = R_LCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk R_R and it is  "<< R_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else if (j == 1) {  // S enantiomer, different ellicity
					sprintf(tempk, "%s", S_RCP_TH1D_Array[k]->GetTitle());
					sprintf(tempj, "%s", S_LCP_TH1D_Array[k]->GetTitle());
					Working_Array_1[k] = S_RCP_TH1D_Array[k];
					Working_Array_2[k] = S_LCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be S_L and it is  "<< S_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else if (j == 2){  // DIFFERENT ellicity DIFFERENT enantiomer
					sprintf(tempj, "%s", R_LCP_TH1D_Array[k]->GetTitle());
					sprintf(tempk, "%s", S_RCP_TH1D_Array[k]->GetTitle());
					Working_Array_1[k] = R_LCP_TH1D_Array[k];
					Working_Array_2[k] = S_RCP_TH1D_Array[k];
					//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
				}
				else {
					cout << "j >= 3 there is a problem!" << ratiocounter << endl << endl;
				}

				// tempj and tempk taken from Working_Array_2
				sprintf(tempjk, "PECD %s vs. %s - cos(theta) = %2.2f - phi = %3.0f", tempj, tempk, cos_theta[k], phi[k]);
				TPaveStats *sk = (TPaveStats*)Working_Array_1[k]->GetListOfFunctions()->FindObject("stats");
				sk->SetOptStat(11);

				//// ASYMMETY plots TH1D
				//asymm_genarray[ratiocounter] = asymm_TH1D(Working_Array_1[k], Working_Array_2[k], cname[j]);
				//asymm_genarray[ratiocounter]->Fit("PECD", "Q"); // console output could e annoying: option "Q" to suppress it
				//fit1[k] = asymm_genarray[ratiocounter]->GetFunction("PECD");
				//b1[ratiocounter] = fit1[k]->GetParameter(0);
				//b1_err[ratiocounter] = fit1[k]->GetParError(0);
				////asymm array for asymmetry plotting b1 with regression curve
				//asymm_genarray[ratiocounter]->SetMarkerStyle(20);
				//asymm_genarray[ratiocounter]->GetXaxis()->SetTitle("cos(theta)");
				//asymm_genarray[ratiocounter]->GetYaxis()->SetTitle("asymm [a.u.]");
				//asymm_genarray[ratiocounter]->GetYaxis()->SetRangeUser(-0.6, 0.6);
				//asymm_genarray[ratiocounter]->SetTitle("");
				////const char* title = ("B1 = " + b1[i]); // example: ratio_genarray[ratiocounter]->SetTitle((const char*)histlist1D[j]);
				////ratio_genarray[ratiocounter]->SetTitle("B1_correct_title_char");
				//canv[j]->cd(k+1);
				//asymm_genarray[ratiocounter]->Draw("AP");
				//fit1[k]->Draw("SAME");

				ratio_genarray[ratiocounter] = ratio_TH1D(Working_Array_1[k], Working_Array_2[k], cname[j]);

				////pol1 fitting
				//ratio_genarray[ratiocounter]->Fit("pol1"); // console output could e annoying: option "Q" to suppress it
				//fit1[k] = ratio_genarray[ratiocounter]->GetFunction("pol1");
				//b1[ratiocounter] = fit1[k]->GetParameter(1);
				//b1_err[ratiocounter] = fit1[k]->GetParError(1);

				//PECD fitting
				// TF1 *f1 = new TF1("ols_lin", "[0]*x+[1]", -1, 1);; //Ordinary least squares fitting from /fitLinearRobust.C
				// TF1 *f1 = new TF1("ols_lin", "pol1", -1, 1); //this should be equivalent at the expressio above
				TF1 *f1 = new TF1("PECD", "([0]*x)/(1-0.125*[1]*(3*x*x-1))+1", -1, 1);; //1D PECD fitting mod (+1 for the divide, 1/8 instead of 1/4 seocnd legendre polinomial)
				f1->SetParameters(0,0);
				f1->SetParLimits(0,-1,1);
				f1->SetParLimits(1,-6.5,6.5);
				f1->SetParNames("b1","b2");
				ratio_genarray[ratiocounter]->Fit("PECD", "Q"); // console output could e annoying: option "Q" to suppress it
				fit1[k] = ratio_genarray[ratiocounter]->GetFunction("PECD");
				b1[ratiocounter] = fit1[k]->GetParameter(0);
				b1_err[ratiocounter] = fit1[k]->GetParError(0);

				//fitfunction(ratio_genarray[ratiocounter], fit1[k], b1[ratiocounter], b1_err[ratiocounter]);

				ratio_genarray[ratiocounter]->SetMarkerStyle(20);
				ratio_genarray[ratiocounter]->GetXaxis()->SetTitle("cos(theta)");
				ratio_genarray[ratiocounter]->GetYaxis()->SetTitle("asymm [a.u.]");
				ratio_genarray[ratiocounter]->GetYaxis()->SetRangeUser(0.4, 1.6); // centered at 1
				ratio_genarray[ratiocounter]->SetTitle("");
				//const char* title = ("B1 = " + b1[i]); // example: ratio_genarray[ratiocounter]->SetTitle((const char*)histlist1D[j]);
				//ratio_genarray[ratiocounter]->SetTitle("B1_correct_title_char");

				canv[j]->cd(k+1);
				ratio_genarray[ratiocounter]->Draw("AP");
				fit1[k]->Draw("SAME");

				DrawTextInPad(0.2, 0.93, tempjk, 0.05);
				gStyle->SetPalette(200,MyPalette);
				gStyle->SetOptTitle(0);
				gStyle->SetOptFit();

				sprintf(tempjk, "b1fit_cos(theta)=%2.2f_phi=%3.0f_counter=%d", cos_theta[k], phi[k], ratiocounter);
				savegraphas(tempjk,ratio_genarray[ratiocounter],fit1[k]);

				//// Set the leggend
				//TLegend *leg = canv[i]->cd(counter)->BuildLegend();
				//leg->SetTextFont(30);
				//leg->SetY1NDC(2);

				//cout << "ratio counter before = " << ratiocounter << endl << endl;
				ratiocounter++;
				//cout << "ratio counter after = " << ratiocounter << endl<< endl;

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
			canv[j]->Write();
			//canv[j]->Close();
			//gSystem->ProcessEvents();
			cout << "plotting and saving cos(theta) canvas " << cname[j] << " DONE! \n" << endl << endl;	
		} else {
			sprintf(tempjk, "b1_map_of_%s", cname[j]);
			canv[j] = new TCanvas(cname[j], cname[j]);

			canv[j]->SetBatch(kFALSE);  // kTRUE to suppress the graphical putput
										// here kFALSE to display the b1maps
#if NOOUTPUT == 1
			canv[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput	
#endif
									  //canv[j]->SetCanvasSize(1920, 1080);
									  //canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[j]->SetFillStyle(4000);
			canv[j]->cd();

			// 12x6
			b1_map[j-3] = new TGraph2D();
			b1_map[j-3]->SetNpx(12);
			b1_map[j-3]->SetNpy(6);			
			b1_map[j-3]->SetTitle("; phi [DEG]; cos(theta) [adm]; b1 [adm]");
			canv[j]->Update();

			////filling b1_maps
			// 12x6
			for (int k = 0; k < 72; k++) //k = all phi and cos(theta) = 12x6
				b1_map[j-3]->SetPoint(k, phi[k], cos_theta[k], b1[k+(j-3)*72]);
			b1_map[j-3]->GetXaxis()->SetRangeUser(-180., 180.); // doesn't work with TGraph2D
			b1_map[j-3]->GetYaxis()->SetRangeUser(-1., 1.);		// doesn't work with TGraph2D
			b1_map[j-3]->GetZaxis()->SetRangeUser(-0.6, 0.6);	// doesn't work with TGraph2D
			b1_map[j-3]->Draw("COLZ"); // this will introduce directly the smoothing
									   //b1_map[j-3]->Draw("LINE");

									   //Save as the b1_maps
			sprintf(tempjk, "%s.png", tempjk);
			canv[j]->SaveAs(tempjk);
			//// 3x3
			//for (int k = 0; k < 9; k++) //k = all phi and cos(theta) = 3x3
			//	b1_map[j-3]->SetPoint(k, phi[k], cos_theta[k], b1[k+(j-3)*9]);

			//// 4x4
			//for (int k = 0; k < 16; k++) //k = all phi and cos(theta) = 4x4
			//	b1_map[j-3]->SetPoint(k, phi[k], cos_theta[k], b1[k+(j-3)*16]);

			// final touches: title and legend
			DrawTextInPad(0.32, 0.93, tempjk, 0.05);
			//TLegend *leg = canv[j]->cd()->BuildLegend();
			//leg->SetTextFont(30);
			//leg->SetY1NDC(2);

			gStyle->SetPalette(200,MyPalette);
			gStyle->SetOptTitle(0);
			canv[j]->Update();
			canv[j]->Write(); //for the b1 maps
			printf("plotting and saving b1_maps canvas %s: DONE! \n", cname[j]);
		}
	}

	for (int j = 0; j < 3; j++) { // MFPDAS smoothed asymmetries
		canv2D[j] = new TCanvas(cname2D[j], cname2D[j]);
		canv2D[j]->SetCanvasSize(1920, 1080);
		canv2D[j]->SetWindowSize(1920 + 4, 1080 + 28);
		canv2D[j]->Divide(12, 6);
		canv2D[j]->SetFillStyle(4000);
		canv2D[j]->SetBatch(kFALSE); // kTRUE to suppress the graphical putput
		
#if NOOUTPUTLOAD == 1
		canv2D[j]->SetBatch(kTRUE); // kTRUE to suppress the graphical putput
#endif

		sprintf(tempj, "FitSlicesY b1 b2 - canvas %s", cname2D[j]);
		canv_fit[j] = new TCanvas(tempj, tempj);
		canv_fit[j]->SetCanvasSize(1920, 1080);
		canv_fit[j]->SetWindowSize(1920 + 4, 1080 + 28);
		canv_fit[j]->Divide(12, 6);
		canv_fit[j]->SetFillStyle(4000);

		gStyle->SetPalette(200,MyPalette);
		gStyle->SetOptTitle(0);
		int counter_2D = 72;
		int ratiocounter_2D = 0;
		for (int k = 0; k < 72; k++) {
			sprintf(tempj, "clean");
			sprintf(tempk, "clean");
			if (j == 0) { // R enantiomer, different ellicity
				sprintf(tempk, "%s", R_RCP_TH2D_Array[k]->GetName());
				sprintf(tempj, "%s", R_LCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = R_RCP_TH2D_Array[k];
				Working_Array2D_2[k] = R_LCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk R_R and it is  "<< R_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
			}
			else if (j == 1) {  // S enantiomer, different ellicity
				sprintf(tempk, "%s", S_RCP_TH2D_Array[k]->GetName());
				sprintf(tempj, "%s", S_LCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = S_RCP_TH2D_Array[k];
				Working_Array2D_2[k] = S_LCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be S_L and it is  "<< S_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
			}
			else if (j == 2) {  // DIFFERENT ellicity DIFFERENT enantiomer
				sprintf(tempj, "%s", R_LCP_TH2D_Array[k]->GetName());
				sprintf(tempk, "%s", S_RCP_TH2D_Array[k]->GetName());
				Working_Array2D_1[k] = R_LCP_TH2D_Array[k];
				Working_Array2D_2[k] = S_RCP_TH2D_Array[k];
				//cout << " j = " << j <<" tempj should be R_L and it is  "<< R_LCP_TH1D_Array[i]->GetTitle() << ", tempk S_R and it is  "<< S_RCP_TH1D_Array[i]->GetTitle() << endl << endl;
			}
			ratioTH2D[ratiocounter_2D] = ratio_TH2D(Working_Array2D_1[k], Working_Array2D_2[k], cname2D[j]);
			//cout << "Number of entries for ratio TH2D is " << ratioTH2D[ratiocounter_2D]->GetEntries() << endl;
			//ratioTH2D[ratiocounter_2D]->SetMarkerStyle(20);

			ratioTH2D[ratiocounter_2D]->GetZaxis()->SetLimits(-0.2, 0.2);
			//ratioTH2D[ratiocounter_2D]->GetZaxis()->SetRangeUser(-0.2, 0.2);
			const char* temptitle = ratioTH2D[ratiocounter_2D]->GetTitle();
			ratioTH2D[ratiocounter_2D]->SetTitle("");
			ratioTH2D[ratiocounter_2D]->Smooth();  //Smooth() : 3 kernels are proposed k5a, k5b and k3a
			gStyle->SetPalette(200,MyPalette);
			TPaveStats *sk = (TPaveStats*)Working_Array2D_1[k]->GetListOfFunctions()->FindObject("stats");
			sk->SetOptStat(11);

			canv2D[j]->cd(counter_2D);
			ratioTH2D[ratiocounter_2D]->Draw();  //avoid COLZ
			//ratioTH2D[ratiocounter_2D]->Draw("CONTZ");
			DrawTextInPad(0.2, 0.93, tempjk, 0.06);

			sprintf(tempjk, "ratio=%d_k=%d_j=%d",ratiocounter_2D, k, j);
			//cout << " readinf file " << tempjk << endl << endl;
			savegraphas(tempjk, ratioTH2D[ratiocounter_2D]);

			canv2D[j]->Update();

			// FITSLICEY for b1 slice
			TObjArray *aSlice = new TObjArray();
			aSlice->SetOwner(kTRUE);
			//TF1 *f1 = new TF1("ols_lin", "[0]*x+[1]", -1, 1);
			TF1 *f1 = new TF1("PECD", "([0]*x)/(1-0.125*[1]*(3*x*x-1))+1", -1, 1);; //1D PECD fitting mod (+1 for the divide, 1/8 instead of 1/4 seocnd legendre polinomial)
			f1->SetParameters(0,0);
			//f1->SetParLimits(0,-1,1);
			f1->SetParLimits(1,-6.5,6.5);
			f1->SetParNames("b1","b2");
			adaptentries(ratioTH2D[ratiocounter_2D]);
			ratioTH2D[ratiocounter_2D]->FitSlicesY(f1, 0, -1, -1, "QR", aSlice);
			copyaSliceTH1D[k] = (TH1D*)(*aSlice)[0]->Clone();
			adaptentries_red(copyaSliceTH1D[k]);
			canv_fit[j]->cd(counter_2D);
			copyaSliceTH1D[k]->SetTitle(temptitle);
			copyaSliceTH1D[k]->SetMarkerStyle(20);
			copyaSliceTH1D[k]->SetMarkerSize(1);
			copyaSliceTH1D[k]->SetStats(0);
			copyaSliceTH1D[k]->GetYaxis()->SetTitle("b1");
			copyaSliceTH1D[k]->GetYaxis()->SetRangeUser(-0.3, 0.3);
			//copyaSliceTH1D[k]->Draw();
			copyaSliceTH1D[k]->Draw("E1");
			//(*aSlice)[0]->Draw();  // draw b1 histogram
			canv_fit[j]->Update();
			sprintf(tempjk, "FitSlicesY b1 %s, canvas %d",cname2D[j], k);
			savegraphas(tempjk, copyaSliceTH1D[k]); // weird passing o aSlice
			//

			ratiocounter_2D++;
			counter_2D--;
		}

		canv_fit[j]->Update();
		canv_fit[j]->Write();
		canv_fit[j]->Close();
		gSystem->ProcessEvents();
		std::printf("plotting and saving FitSlicesY smoothed canvas %s: DONE! \n \n", cname2D[j]);

		canv2D[j]->Update(); // ??? is it usefull?
		canv2D[j]->Write();
		canv2D[j]->Close();
		gSystem->ProcessEvents();
		std::printf("plotting and saving TH2D smoothed canvas %s: DONE! \n \n", cname2D[j]);
	}

	f->Close();
	std::cout << "gPad : SAVED! \n " << endl << endl;
	fclose(fp);
	// normal exit with no delete of variable in c++: exit(0) is functionally indentical to return 0
	exit(0);
	theApp.Run();
	//terminate();

	return 0;
}
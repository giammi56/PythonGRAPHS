// PECD_graphs.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"
#include <atlstr.h>

#include <TApplication.h>

#include <TDirectory.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
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

#include <iostream>
#include <vector>


void DrawTextInPad(double xCoord, double yCoord, char* text, Float_t size) {
	TText* t1 = new TText(xCoord, yCoord, text);
	t1->SetNDC();
	t1->SetTextColor(1);
	t1->SetTextSize(size);
	t1->SetTextFont(62);
	t1->Draw();
}

void load_TH1D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH1D** S_LCP_TH1D_Array, TH1D** S_RCP_TH1D_Array, TH1D** R_RCP_TH1D_Array, TH1D** R_LCP_TH1D_Array, TH1D** general_Array, TString* dirlist, TString* histlist1D) {

	int counter = 0;
	for(int i=0;i<3;i++) {
		//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
		TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		gDirectory->pwd();
		//for (int j = 0; j < LCPhistlist->Sizeof(); j++) {
		for (int j=0;j<3;j++) {
			S_LCP_TH1D_Array[j+i*3]=(TH1D*)S_Ldir->FindObjectAny(histlist1D[j]);
			S_LCP_TH1D_Array[j+i*3]->SetName("S_LCP"+dirlist[i]);
			S_LCP_TH1D_Array[j+i*3]->SetTitle("S_LCP"+dirlist[i]);
			S_LCP_TH1D_Array[j+i*3]->Sumw2();
			S_LCP_TH1D_Array[j+i*3]->Scale(1. / S_LCP_TH1D_Array[j+i*3]->Integral());
			std::cout << "Normalization of " << S_LCP_TH1D_Array[j+i*3]->GetName() << "/" << histlist1D[j] << " done! # entries = " << S_LCP_TH1D_Array[j+i*3]->GetEntries() << "\n";
																								 
			general_Array[counter++] = S_LCP_TH1D_Array[j+i*3];										
																								 
																								 
			S_RCP_TH1D_Array[j+i*3]=(TH1D*)S_Rdir->FindObjectAny(histlist1D[j]);				 
			S_RCP_TH1D_Array[j+i*3]->SetName("S_RCP"+dirlist[i]);								 
			S_RCP_TH1D_Array[j+i*3]->SetTitle("S_RCP"+dirlist[i]);								 
			S_RCP_TH1D_Array[j+i*3]->Sumw2();													 
			S_RCP_TH1D_Array[j+i*3]->Scale(1. / S_RCP_TH1D_Array[j+i*3]->Integral());			 
			std::cout << "Normalization of " << S_RCP_TH1D_Array[j+i*3]->GetName() << "/" << histlist1D[j] << " done! # entries = " << S_RCP_TH1D_Array[j+i*3]->GetEntries() << "\n";
																								 
			general_Array[counter++] = S_RCP_TH1D_Array[j+i*3];									 
																								 
			R_LCP_TH1D_Array[j+i*3] = (TH1D*)R_Ldir->FindObjectAny(histlist1D[j]);				 
			R_LCP_TH1D_Array[j+i*3]->SetName("R_LCP"+dirlist[i]);								 
			R_LCP_TH1D_Array[j+i*3]->SetTitle("R_LCP"+dirlist[i]);								 
			R_LCP_TH1D_Array[j+i*3]->Sumw2();													 
			R_LCP_TH1D_Array[j+i*3]->Scale(1. / R_LCP_TH1D_Array[j+i*3]->Integral());			 
			std::cout << "Normalization of " << R_LCP_TH1D_Array[j+i*3]->GetName() << "/" << histlist1D[j] << " done! # entries = " << R_LCP_TH1D_Array[j+i*3]->GetEntries() << "\n";
																								 
			general_Array[counter++] = R_LCP_TH1D_Array[j+i*3];									 
																								 
			R_RCP_TH1D_Array[j+i*3] = (TH1D*)R_Rdir->FindObjectAny(histlist1D[j]);				 
			R_RCP_TH1D_Array[j+i*3]->SetName("R_RCP"+dirlist[i]);								 
			R_RCP_TH1D_Array[j+i*3]->SetTitle("R_RCP"+dirlist[i]);								 
			R_RCP_TH1D_Array[j+i*3]->Sumw2();													 
			R_RCP_TH1D_Array[j+i*3]->Scale(1. / R_RCP_TH1D_Array[j+i*3]->Integral());			 
			std::cout << "Normalization of " << R_RCP_TH1D_Array[j+i*3]->GetName() << "/" << histlist1D[j] << " done! # entries = " << R_RCP_TH1D_Array[j+i*3]->GetEntries() << "\n";
			
			general_Array[counter++] = R_RCP_TH1D_Array[j+i*3];

		}
	}

	std::cout << "Total number of TH1D histos = " << counter << "\n";

	//Autoput of all loaded graphs
	TCanvas* c56 = new TCanvas("all TH1D histos", "all = 36");
	c56->SetCanvasSize(1920, 1080);
	c56->SetWindowSize(1920+4, 1080+28);
	c56->Divide(6,6);
	for (int k=0; k<36; k++) {
		c56->cd(k+1);
		//general_Array[k]->SetOption("E");
		general_Array[k]->Draw("E1");
		//gPad->Modified();
		gPad->Update();
	}
	//list->Add(c56);
	c56->Write();
	std::cout << "\n";

	//TCanvas* c55 = new TCanvas("EN_gate", "all graphs = 12");
	//c55->SetCanvasSize(1920, 1080);
	//c55->SetWindowSize(1920+4, 1080+28);
	//c55->Divide(4,3);
	//for (int k=0; k<12; k++) {
	//	c55->cd(k+1);
	//	//general_Array[k]->SetOption("E");
	//	general_Array[k]->Draw("E1");
	//	//gPad->Modified();
	//	gPad->Update();
	//}
	////list->Add(c55);
	//c55->Write();
	//std::cout << "\n";
	////!Autoput all loaded graphs

	return;
	
}

void normalize_TH2D(TH2D* TH2D_Array) {

	Double_t norm = 1. / TH2D_Array->Integral();
	TH2D_Array->Scale(norm);
	std::cout << TH2D_Array->GetName() << ": TH2D normalized" << "\n \n";

	//int binx = 0;
	//int biny = 0;
	//double integral = 0.;
	////TH2D_Array->Sumw2();
	//binx = TH2D_Array->GetXaxis()->GetNbins();
	////std::cout << "BINX done!" << "\n";
	//biny = TH2D_Array->GetYaxis()->GetNbins();
	////std::cout << "BINY done!" << "\n";
	//for (int x = 1; x <= binx; ++x) {		
	//	for (int y = 1; y <= biny; ++y) {
	//		integral += TH2D_Array->GetBinContent(x, y);
	//		TH2D_Array->SetBinContent(x, y, TH2D_Array->GetBinContent(x, y) / integral);
	//	}
	//}
	
	//TH2D_Array->Scale();
	return;
}

void load_TH2D(TFile* S_LFile, TFile* S_RFile, TFile* R_LFile, TFile* R_RFile, TH2D** S_LCP_TH2D_Array, TH2D** S_RCP_TH2D_Array, TH2D** R_LCP_TH2D_Array, TH2D** R_RCP_TH2D_Array, TH2D** general2D_Array, TString* dirlist, TString* histlist2D,  TString* namelist2D) {

	int counter = 0;
	//just EN_gate -> i=0
	for(int i=0;i<1;i++) {
		//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
		TDirectory* S_Ldir = (TDirectory*)(S_LFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* S_Rdir = (TDirectory*)(S_RFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory("angular_distr_el/ID3_mol_e0_valid"+dirlist[i]));
		gDirectory->pwd();
		
		//2 PECD graphs 
		for (int j=0;j<2;j++) {
			
			//*3 depends on the size of dirlist
			S_LCP_TH2D_Array[j+i*3]=(TH2D*)S_Ldir->FindObjectAny(histlist2D[j]);
			S_LCP_TH2D_Array[j+i*3]->SetName("S_LCP"+namelist2D[j]);
			S_LCP_TH2D_Array[j+i*3]->SetTitle("S_LCP"+namelist2D[j]);
			//S_LCP_TH2D_Array[j+i*3]->Sumw2();
			normalize_TH2D(S_LCP_TH2D_Array[j+i*3]);
			general2D_Array[counter++] = S_LCP_TH2D_Array[j+i*3];									 


			S_RCP_TH2D_Array[j+i*3]=(TH2D*)S_Rdir->FindObjectAny(histlist2D[j]);				 
			S_RCP_TH2D_Array[j+i*3]->SetName("S_RCP"+namelist2D[j]);								 
			S_RCP_TH2D_Array[j+i*3]->SetTitle("S_RCP"+namelist2D[j]);								 
			//S_RCP_TH2D_Array[j+i*3]->Sumw2();													 
			normalize_TH2D(S_RCP_TH2D_Array[j+i*3]);
			general2D_Array[counter++] = S_RCP_TH2D_Array[j+i*3];									 

			R_LCP_TH2D_Array[j+i*3] = (TH2D*)R_Ldir->FindObjectAny(histlist2D[j]);				 
			R_LCP_TH2D_Array[j+i*3]->SetName("R_LCP"+namelist2D[j]);								 
			R_LCP_TH2D_Array[j+i*3]->SetTitle("R_LCP"+namelist2D[j]);								 
			//R_LCP_TH2D_Array[j+i*3]->Sumw2();													 
			normalize_TH2D(R_LCP_TH2D_Array[j+i*3]);
			general2D_Array[counter++] = R_LCP_TH2D_Array[j+i*3];									 

			R_RCP_TH2D_Array[j+i*3] = (TH2D*)R_Rdir->FindObjectAny(histlist2D[j]);				 
			R_RCP_TH2D_Array[j+i*3]->SetName("R_RCP"+namelist2D[j]);								 
			R_RCP_TH2D_Array[j+i*3]->SetTitle("R_RCP"+namelist2D[j]);								 
			//R_RCP_TH2D_Array[j+i*3]->Sumw2();													 
			normalize_TH2D(R_RCP_TH2D_Array[j+i*3]);
			general2D_Array[counter++] = R_RCP_TH2D_Array[j+i*3];

		}
	}

	std::cout << "Total number of TH2D histos = " << counter+1 << "\n";

	//Autoput of all loaded graphs
	TCanvas* c60 = new TCanvas("all full TH2 histos", "all histos = 4");
	c60->SetCanvasSize(1920, 1080);
	c60->SetWindowSize(1920+4, 1080+28);
	c60->Divide(2,2);
	for (int k=0; k<4; k++) {
		c60->cd(k+1);
		//general_Array[k]->SetOption("E");
		general2D_Array[k]->Draw("COLZ");
		//gPad->Modified();
		gPad->Update();
	}
	//list->Add(c56);
	c60->Write();
	std::cout << "\n";

	TCanvas* c61 = new TCanvas("all reduced TH2 histos", "all histos = 4");
	c61->SetCanvasSize(1920, 1080);
	c61->SetWindowSize(1920+4, 1080+28);
	c61->Divide(2,2);
	for (int k=0; k<4; k++) {
		c61->cd(k+1);
		//general_Array[k]->SetOption("E");
		general2D_Array[k+4]->Draw("COLZ");
		//gPad->Modified();
		gPad->Update();
	}
	//list->Add(c55);
	c61->Write();
	std::cout << "\n";
	////!Autoput all loaded graphs

	return;

}

TGraphAsymmErrors* ratio_TH1D(TH1D* &LCPTH1Darray, TH1D* &RCPTH1Darray, const char* &cname) {

	int a = LCPTH1Darray->GetNbinsX();
	//std::cout << "a = " << a << "\n";
	// Asymmetry = (h1 - h2)/(h1 + h2)  where h1 = this
	auto ratioTH1D = new TGraphAsymmErrors(a);
	// https://en.wikipedia.org/wiki/Poisson_distribution#Assumptions:_When_is_the_Poisson_distribution_an_appropriate_model? //
	// second option in Divide() midp = lancaster mid-p_value, n = normal propagation
	ratioTH1D->Divide(LCPTH1Darray, RCPTH1Darray, "pois n");
	printf("Division of %s and %s in %s: DONE! Size ratio vector = %d \n",LCPTH1Darray->GetTitle(), RCPTH1Darray->GetTitle(), cname, ratioTH1D->Sizeof()/sizeof(int));

	return ratioTH1D;
}

TH2D* function_ratio_TH2D(TH2D* &LCP_TH2D_array, TH2D* &RCP_TH2D_array, const char* &cname) {

	TH2D* ratioTH2D = (TH2D*)LCP_TH2D_array->GetAsymmetry(RCP_TH2D_array);
	//std::cout << "Size of ratioTH2D: " << (ratioTH2D->Sizeof())/sizeof(int) << "\n";
	std::cout << "Division TH2D DONE!" << "\n";
	
	return ratioTH2D;
}

int main(int argc, char** argv)

{
	TApplication theApp("App",&argc,argv);

	// single TH1D arrays
	TH1D** S_LCP_TH1D_Array;
	S_LCP_TH1D_Array = new TH1D*[15];

	TH1D** R_LCP_TH1D_Array;
	R_LCP_TH1D_Array = new TH1D*[15];

	TH1D** S_RCP_TH1D_Array;
	S_RCP_TH1D_Array = new TH1D*[15];

	TH1D** R_RCP_TH1D_Array;
	R_RCP_TH1D_Array = new TH1D*[15];

	TH1D** general_Array;
	general_Array = new TH1D*[60];

	// TH2D arrays!!!
	TH2D** S_LCP_TH2D_Array;
	S_LCP_TH2D_Array = new TH2D*[15];

	TH2D** R_LCP_TH2D_Array;
	R_LCP_TH2D_Array = new TH2D*[15];

	TH2D** S_RCP_TH2D_Array;
	S_RCP_TH2D_Array = new TH2D*[15];
	
	TH2D** R_RCP_TH2D_Array;
	R_RCP_TH2D_Array = new TH2D*[15];

	TH2D** general2D_Array;
	general2D_Array = new TH2D*[60];

	//Ratios
	TH2D** ratioTH2D;
	ratioTH2D = new TH2D*[15];

	TGraphAsymmErrors** ratio_genarray;
	ratio_genarray = new TGraphAsymmErrors*[60];

	// DO NOT CHANGE THE ORDER OF THESE INPUTS
	TFile* S_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CL_11800-1400ns_user_LUT.root", "READ");
	if (S_LFile->IsOpen()) printf("%s: File opened successfully!\n", S_LFile->GetName());
	TFile* S_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/S-C3H3F3O_550eV_CR_11800-1400ns_user_LUT.root", "READ");
	if (S_RFile->IsOpen()) printf("%s: File opened successfully!\n", S_RFile->GetName());
	TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_11800-1400ns_user_LUT.root", "READ");
	if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_11800-1400ns_user_LUT.root", "READ");
	if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());
	
	TString dirlist[3] = {"/EN_gate","/mom_cm.y>0","/mom_cm.y<0"};
	TString histlist1D[3] = {"cos(theta)y_e[0]","cos(theta)_e[0]","cos(theta)0_Ppol"};
	TString namelist2D[3] = {"/EN_gate","/EN_gate_red"};
	TString histlist2D[2] = {"PECD_e[0](theta)", "PECD_e[0](theta)_redTHETA"};

	//FILE NAME
	TFile *f = new TFile("FMtox_ID3_m41_69_LUT.root","RECREATE");
	if (f->IsOpen()) {
		//TFile *f = new TFile("FMtox_ID3_m41_69.root", "RECREATE");
		printf("  File CREATED successfully\n");
	}
	// LOAD AND NORMALIZED ALL GRAPHS
	load_TH1D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH1D_Array, S_RCP_TH1D_Array, R_LCP_TH1D_Array, R_RCP_TH1D_Array, general_Array, dirlist, histlist1D);
	std::cout << "load_all TH1D : DONE! \n \n";

	load_TH2D(S_LFile, S_RFile, R_LFile, R_RFile, S_LCP_TH2D_Array, S_RCP_TH2D_Array, R_LCP_TH2D_Array, R_RCP_TH2D_Array, general2D_Array, dirlist, histlist2D, namelist2D);
	std::cout << "load_all TH2D : DONE! \n \n";

	// canvas strig name, from Till suggestion! Not necessary now.
	const char *cname[] = {"e[0]_NOgate","e[0]_y>0","e[0]_y<0", "e[0]_NOgate_ASYM","e[0]_y>0_ASYM","e[0]_y<0_ASYM","PECD_map_ASYM", "PECD_reduced-map_ASYM","PECD_map_ASYM_smooth", "PECD_reduced-map_ASYM_smooth"};
	TCanvas *canv[7]; // dim(cnam NON ASYMMETRIC + TH2D graphs)
	std::cout << "Vector of TCanvas created." << "\n \n";

	char tempj[250];
	char tempk[250];
	char tempjk[250];

	int ratio_counter = 0;
	int ratio_counter_2D = 0;

	// JUST EN_gate
	for (int i = 0; i < 7; i++) { // dim(cnam NON ASYMMETRIC + TH2D graphs)
		if (i < 3) { // TH1D
			//create normal canvas
			canv[i] = new TCanvas(cname[i], cname[i]);
			canv[i]->SetCanvasSize(1920, 1080);
			canv[i]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[i]->Divide(4, 4);
			canv[i]->SetFillStyle(4000);
			gStyle->SetOptTitle(0);

			//create the ASYMMETRIC canvas
			canv[i + 3] = new TCanvas(cname[i + 3], cname[i + 3]);
			canv[i + 3]->SetCanvasSize(1920, 1080);
			canv[i + 3]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[i + 3]->Divide(4, 4);
			canv[i + 3]->SetFillStyle(4000);
			//gStyle->SetOptTitle(0);

			// counter for canvas (4*4)
			int counter = 16;
			for (int j = 0; j < 4; j++) { //S enantiomer
				sprintf(tempj, "%s", general_Array[j + i * 4]->GetTitle());
				for (int k = 0; k < 4; k++) { //R enantiomer
					canv[i]->cd(counter);
					sprintf(tempk, "%s", general_Array[k + i * 4]->GetTitle());

					// first pass of array J
					general_Array[j + i * 4]->SetMarkerStyle(20);
					general_Array[j + i * 4]->SetMarkerSize(1);
					general_Array[j + i * 4]->GetYaxis()->SetRangeUser(0.08, 0.115);
					TPaveStats *sj = (TPaveStats*)general_Array[j + i * 4]->GetListOfFunctions()->FindObject("stats");
					sj->SetOptStat(11);

					// if equal to 0 -> S is red, R is blue
					if ((std::strcmp(general_Array[j + i * 4]->GetName(), "S_RCP/EN_gate") == 0) || (std::strcmp(general_Array[j + i * 4]->GetName(), "S_LCP/EN_gate") == 0)) {
						general_Array[j + i * 4]->SetMarkerColor(kRed);
					}
					else {
						general_Array[j + i * 4]->SetMarkerColor(kBlue);
					}
					// if equal to 0 -> LCP has star
					if ((std::strcmp(general_Array[j + i * 4]->GetName(), "S_LCP/EN_gate") == 0) || (std::strcmp(general_Array[j + i * 4]->GetName(), "R_LCP/EN_gate") == 0))
						general_Array[j + i * 4]->SetMarkerStyle(29);

					general_Array[j + i * 4]->Draw("E1");  //"AP" -> you need TGrapherrors?

					// second pass of array K
					general_Array[k + i * 4]->SetMarkerStyle(20);
					general_Array[k + i * 4]->SetMarkerSize(1);
					general_Array[k + i * 4]->GetXaxis()->SetTitle("cos(theta)");
					general_Array[k + i * 4]->GetYaxis()->SetTitle("counts [a.u.]");
					general_Array[k + i * 4]->GetYaxis()->SetRangeUser(0.08, 0.115);
					//general_Array[k+i*4]->SetOption("E SAME");

					// 0 if equal -> S is red, R is blue
					if ((std::strcmp(general_Array[k + i * 4]->GetName(), "S_RCP/EN_gate") == 0) || (std::strcmp(general_Array[k + i * 4]->GetName(), "S_LCP/EN_gate") == 0)) {
						general_Array[k + i * 4]->SetMarkerColor(kRed);
					}
					else {
						general_Array[k + i * 4]->SetMarkerColor(kBlue);
					}
					// if equal to 0 -> LCP has star
					if ((std::strcmp(general_Array[k + i * 4]->GetName(), "S_LCP/EN_gate") == 0) || (std::strcmp(general_Array[k + i * 4]->GetName(), "R_LCP/EN_gate") == 0))
						general_Array[k + i * 4]->SetMarkerStyle(29);

					// j+i*4 = k+i*4 causes nname to go crazy
					sprintf(tempjk, "%s vs. %s", tempj, tempk);
					TPaveStats *sk = (TPaveStats*)general_Array[k + i * 4]->GetListOfFunctions()->FindObject("stats");
					sk->SetOptStat(11);
					general_Array[k + i * 4]->SetBit(TH1::kNoTitle);
					general_Array[k + i * 4]->Draw("E1SAMES"); // "AP" or hs->Draw("nostack,e1p");

					DrawTextInPad(0.2, 0.93, tempjk, 0.06);
					canv[i]->Update();
					//// leggend
					//TLegend *leg = canv[i]->cd(counter)->BuildLegend();
					//leg->SetTextFont(30);
					//leg->SetY1NDC(2);

					// ASYMMETRIC plots
					canv[i + 3]->cd(counter);
					ratio_genarray[ratio_counter] = ratio_TH1D(general_Array[j + i * 4], general_Array[k + i * 4], cname[i]);
					ratio_genarray[ratio_counter]->SetMarkerStyle(20);
					ratio_genarray[ratio_counter]->GetXaxis()->SetTitle("cos(theta)");
					ratio_genarray[ratio_counter]->GetYaxis()->SetTitle("asymm [a.u.]");
					ratio_genarray[ratio_counter]->GetYaxis()->SetRangeUser(0.85, 1.15);
					//ratio_genarray[ratio_counter]->SetTitle((const char*)histlist1D[j]);
					ratio_genarray[ratio_counter]->SetTitle("");
					ratio_genarray[ratio_counter]->Draw("AP");

					DrawTextInPad(0.2, 0.93, tempjk, 0.06);
					canv[i + 3]->Update();

					ratio_counter++;
					counter--;
				}
			}

			//list->Add(canv[i]);
			canv[i]->Write();
			printf("plotting and saving TH1D convas %s: DONE! \n", cname[i]);
			canv[i + 3]->Write();
			printf("plotting and saving TH1D canvas %s: DONE! \n \n", cname[i + 3]);
		}
		else if (i > 4) { // TH2D
		//canv[i] = new TCanvas(cname[i], cname[i]);
		//canv[i]->SetCanvasSize(1920, 1080);
		//canv[i]->SetWindowSize(1920 + 4, 1080 + 28);
		//canv[i]->Divide(2, 2);
		//canv[i]->SetFillStyle(4000);
		//gStyle->SetOptTitle(0);

		//just asymmetry plot for the PECD map
			int z =i-3;  //number of iteration non TH2D
			canv[i + 3] = new TCanvas(cname[i + 3], cname[i + 3]);
			canv[i + 3]->SetCanvasSize(1920, 1080);
			canv[i + 3]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[i + 3]->Divide(4, 4);
			canv[i + 3]->SetFillStyle(4000);
			gStyle->SetOptTitle(0);

			int counter_2D = 16;  //(i-2)*4 to return back to the baginning of the arraz sicne we processed two PECD graphs
			for (int j = 0; j < 4; j++) {
				sprintf(tempj, "%s", general2D_Array[j+z*4]->GetTitle());
				for (int k = 0; k < 4; k++) {
					sprintf(tempk, "%s", general2D_Array[k+z*4]->GetTitle());
					canv[i + 3]->cd(counter_2D);
					ratioTH2D[ratio_counter_2D] = function_ratio_TH2D(general2D_Array[j+z*4], general2D_Array[k+z*4], cname[i+3]);
					//ratioTH2D[ratio_counter_2D]->SetMarkerStyle(20);
					TPaveStats *sk = (TPaveStats*)general2D_Array[k+z*4]->GetListOfFunctions()->FindObject("stats");
					//sk->SetOptStat(11);
					ratioTH2D[ratio_counter_2D]->GetZaxis()->SetRangeUser(-0.17, 0.17);
					ratioTH2D[ratio_counter_2D]->SetTitle("");
					//Smooth() : 3 kernels are proposed k5a, k5b and k3a
					ratioTH2D[ratio_counter_2D]->Smooth();
					ratioTH2D[ratio_counter_2D]->Draw("COLZ");
					//ratioTH2D[ratio_counter_2D]->Draw("CONTZ");

					sprintf(tempjk, "%s vs. %s", tempj, tempk);
					DrawTextInPad(0.2, 0.93, tempjk, 0.06);

					canv[i + 3]->Update();

					ratio_counter_2D++;
					counter_2D--;
				}
			}
			canv[i + 3]->Write();
			printf("plotting and saving TH2D canvas %s: DONE! \n \n", cname[i + 3]);
		}
		else { // TH2D
		   //canv[i] = new TCanvas(cname[i], cname[i]);
		   //canv[i]->SetCanvasSize(1920, 1080);
		   //canv[i]->SetWindowSize(1920 + 4, 1080 + 28);
		   //canv[i]->Divide(2, 2);
		   //canv[i]->SetFillStyle(4000);
		   //gStyle->SetOptTitle(0);

		   //just asymmetry plot for the PECD map
			int z =i-3;  //number of iteration non TH2D
			canv[i + 3] = new TCanvas(cname[i + 3], cname[i + 3]);
			canv[i + 3]->SetCanvasSize(1920, 1080);
			canv[i + 3]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[i + 3]->Divide(4, 4);
			canv[i + 3]->SetFillStyle(4000);
			gStyle->SetOptTitle(0);

			int counter_2D = 16;
			for (int j = 0; j < 4; j++) {
				sprintf(tempj, "%s", general2D_Array[j+z*4]->GetTitle());
				for (int k = 0; k < 4; k++) {
					sprintf(tempk, "%s", general2D_Array[k+z*4]->GetTitle());
					canv[i + 3]->cd(counter_2D);
					ratioTH2D[ratio_counter_2D] = function_ratio_TH2D(general2D_Array[j+z*4], general2D_Array[k+z*4], cname[i + 3]);
					//ratioTH2D[ratio_counter_2D]->SetMarkerStyle(20);
					TPaveStats *sk = (TPaveStats*)general2D_Array[k + i]->GetListOfFunctions()->FindObject("stats");
					//sk->SetOptStat(11);
					ratioTH2D[ratio_counter_2D]->GetZaxis()->SetRangeUser(-0.17, 0.17);
					ratioTH2D[ratio_counter_2D]->SetTitle("");
					ratioTH2D[ratio_counter_2D]->Draw("COLZ");
					//ratioTH2D[ratio_counter_2D]->Draw("CONTZ");

					sprintf(tempjk, "%s vs. %s", tempj, tempk);
					DrawTextInPad(0.2, 0.93, tempjk, 0.06);

					canv[i + 3]->Update();

					ratio_counter_2D++;
					counter_2D--;
				}
			}
			canv[i + 3]->Write();
			printf("plotting and saving TH2D canvas %s: DONE! \n \n", cname[i + 3]);
		}
	}
	f->Close();
	std::cout << "gPad : SAVED!" << "\n \n";
	
	// hard exit with no delete of variable
	exit(0);
	theApp.Run();

	return 0;
}
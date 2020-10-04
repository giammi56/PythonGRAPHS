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
using namespace std;
#include <vector>


void DrawTextInPad(double xCoord, double yCoord, char* text, Float_t size) {
	TText* t1 = new TText(xCoord, yCoord, text);
	t1->SetNDC();
	t1->SetTextColor(1);
	t1->SetTextSize(size);
	t1->SetTextFont(62);
	t1->Draw();
}

// TAKE CARE THAT THE FEEDING PARAMETERS ARE IN TEH SAME ORDER!!!
void load_TH1D(TFile* R_LFile, TFile* R_RFile, TH1D** R_LCP_TH1D_Array, TH1D** R_RCP_TH1D_Array, TH1D** general_Array, TString* dirlist, TString* histlist1D, const char dir[]) {
	int counter = 0;
	int i = 0;
	//TDirectoryFile* LCPDir = (TDirectoryFile*)(LCPFile->GetDirectory("angular_distr_el/if_twoelec"+dirlist[i]));
	TDirectory* R_Ldir = (TDirectory*)(R_LFile->GetDirectory(dir+dirlist[i]));
	TDirectory* R_Rdir = (TDirectory*)(R_RFile->GetDirectory(dir+dirlist[i]));
	gDirectory->pwd();
	//for (int j = 0; j < LCPhistlist->Sizeof(); j++) {

	//12x6
	for (int j=0;j<72;j++) {
		R_LCP_TH1D_Array[j] = (TH1D*)R_Ldir->FindObjectAny(histlist1D[j]);				 
		R_LCP_TH1D_Array[j]->SetName("R_LCP");								 
		R_LCP_TH1D_Array[j]->SetTitle("R_LCP");								 
		R_LCP_TH1D_Array[j]->Sumw2();													 
		R_LCP_TH1D_Array[j]->Scale(1. / R_LCP_TH1D_Array[j]->Integral());			 
		std::cout << "Loading and normalisation of " << R_LCP_TH1D_Array[j]->GetName() << "/" << histlist1D[j] << " done! # entries = " << R_LCP_TH1D_Array[j]->GetEntries() << "\n";

		general_Array[counter++] = R_LCP_TH1D_Array[j];									 

		R_RCP_TH1D_Array[j] = (TH1D*)R_Rdir->FindObjectAny(histlist1D[j]);				 
		R_RCP_TH1D_Array[j]->SetName("R_RCP");								 
		R_RCP_TH1D_Array[j]->SetTitle("R_RCP");								 
		R_RCP_TH1D_Array[j]->Sumw2();													 
		R_RCP_TH1D_Array[j]->Scale(1. / R_RCP_TH1D_Array[j]->Integral());			 
		std::cout << "Loading and normalisation of " << R_RCP_TH1D_Array[j]->GetName() << "/" << histlist1D[j] << " done! # entries = " << R_RCP_TH1D_Array[j]->GetEntries() << "\n";

		general_Array[counter++] = R_RCP_TH1D_Array[j];
	}

	std::cout << "Total number of TH1D histos = " << counter << "\n";

	//Autoput of all loaded graphs
	char canvastitle[100];
	sprintf(canvastitle, "all entries = %d", counter);

	TCanvas* c56 = new TCanvas("all TH1D histos", canvastitle, sizeof(general_Array));
	c56->SetCanvasSize(1920, 1080);
	c56->SetWindowSize(1920+4, 1080+28);
	c56->SetBatch(kTRUE); // to suppress graphical putput

						  //12x6
	c56->Divide(9,16);
	for (int k=0; k<144; k++) {
		c56->cd(k+1);
		general_Array[k]->Draw("E1");
		//gPad->Modified();
		gPad->Update();
	}
	//list->Add(c56);
	c56->Write();
	c56->SetBatch(kFALSE); // to suppress graphical putput
	std::cout << "\n";
	return;
}

TGraphAsymmErrors* ratio_TH1D(TH1D* &Working_Array_1, TH1D* &Working_Array_2, const char* &cname) {

	int a = Working_Array_1->GetNbinsX();
	//std::cout << "a = " << a << "\n";
	// Asymmetry = (h1 - h2)/(h1 + h2)  where h1 = this
	auto ratioTH1D = new TGraphAsymmErrors(a);
	// https://en.wikipedia.org/wiki/Poisson_distribution#Assumptions:_When_is_the_Poisson_distribution_an_appropriate_model? //
	// second option in Divide() midp = lancaster mid-p_value, n = normal propagation
	ratioTH1D->Divide(Working_Array_1, Working_Array_2, "pois n");
	printf("Division of %s and %s in %s: DONE! Size ratio vector = %d \n",Working_Array_1->GetTitle(), Working_Array_2->GetTitle(), cname, ratioTH1D->Sizeof()/sizeof(int));

	return ratioTH1D;
}

int main(int argc, char** argv)

{
	TApplication theApp("App", &argc, argv);

	// single TH1D arrays: adapt to the requested dimension i*j
	int dim = 144; //12x6x2
	TH1D** R_LCP_TH1D_Array;
	R_LCP_TH1D_Array = new TH1D*[dim];

	TH1D** R_RCP_TH1D_Array;
	R_RCP_TH1D_Array = new TH1D*[dim];

	TH1D** general_Array;
	general_Array = new TH1D*[dim];

	TGraphAsymmErrors** ratio_genarray;
	ratio_genarray = new TGraphAsymmErrors*[dim];

	// Fititign parameters
	TF1** fit1;
	fit1 = new TF1*[];
	// fit1[i] = TF1("ols_lin", "[0]*x+[1]", -1, 1);; //Ordinary least squares fitting from /fitLinearRobust.C
	// fit1[i] = TF1("ols_lin", "pol1", -1, 1); //this should be equivalent at the expressio above

	//Adapt dimensio to request
	Double_t b1[144];
	Double_t b1_err[144];

	double cos_theta[144/2];
	double phi[144/2];

	TGraph2D** b1_map;
	b1_map = new TGraph2D*[];

	// DO NOT CHANGE THE ORDER OF THESE INPUTS R-C3H3F3O_546eV_CL_10800-2350ns_multiCH11_MFPAD_30.root
	//546 eV
	TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_546eV_CL_9600-3700ns_multiCH9_MFPAD_30.root", "READ");
	if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_546eV_CR_9600-3700ns_multiCH9_MFPAD_30.root", "READ");
	if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//550eV
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CL_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/R-C3H3F3O_550eV_CR_10800-2350ns_newEfield_multiCH11_MFPAD_30.root", "READ");
	//if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//// DO NOT CHANGE THE ORDER OF THESE INPUTS
	//TFile* R_LFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CL_reaction_POLY_3x3_test.root", "READ");
	//if (R_LFile->IsOpen()) printf("%s: File opened successfully!\n", R_LFile->GetName());
	//TFile* R_RFile = new TFile("E:/Soleil_2018_09_SEX/Analysis/POLY/R-C3H3F3O_poly_30-69m_550eV_CR_reaction_POLY_3x3_test.root", "READ");
	//if (R_RFile->IsOpen()) printf("%s: File opened successfully!\n", R_RFile->GetName());

	//FILE NAME
	//POLIATOMIC - NEUTRAL
	TFile *f = new TFile("FMtox_546eV_ID_ALL_CH9m69_12x6.root","RECREATE");
	// POLYATOMIC - COMPLETE
	//TFile *f = new TFile("FMtox_ID_ALL_m69_3x3_test.root","RECREATE");
	//TFile *f = new TFile("FMtox_ID3_m41_69.root", "RECREATE");

	//File to save outout b1
	FILE *fp = fopen("b1_546eV_ID_ALL_CH9m69_12x8.dat","w");
	fprintf(fp,"cos(theta) \t phi \t p1 \t p1_err \n");

	if (f->IsOpen()) {
		printf("File CREATED successfully\n");
	}
	const char dir[] = {"angular_distr_el/ID_ALL_mol_e0_valid"};

	TString dirlist[1] = {"/EN_gate/MFPADs"};
	//12x6
	TString histlist1D[72] = { "cos(theta)_e[0]_costheta_-1.00_phi_-180","cos(theta)_e[0]_costheta_-0.67_phi_-180","cos(theta)_e[0]_costheta_-0.33_phi_-180","cos(theta)_e[0]_costheta_-0.00_phi_-180","cos(theta)_e[0]_costheta_0.33_phi_-180","cos(theta)_e[0]_costheta_0.67_phi_-180",
		"cos(theta)_e[0]_costheta_-1.00_phi_-150","cos(theta)_e[0]_costheta_-0.67_phi_-150","cos(theta)_e[0]_costheta_-0.33_phi_-150","cos(theta)_e[0]_costheta_-0.00_phi_-150","cos(theta)_e[0]_costheta_0.33_phi_-150","cos(theta)_e[0]_costheta_0.67_phi_-150",
		"cos(theta)_e[0]_costheta_-1.00_phi_-120","cos(theta)_e[0]_costheta_-0.67_phi_-120","cos(theta)_e[0]_costheta_-0.33_phi_-120","cos(theta)_e[0]_costheta_-0.00_phi_-120","cos(theta)_e[0]_costheta_0.33_phi_-120","cos(theta)_e[0]_costheta_0.67_phi_-120",
		"cos(theta)_e[0]_costheta_-1.00_phi_-90","cos(theta)_e[0]_costheta_-0.67_phi_-90","cos(theta)_e[0]_costheta_-0.33_phi_-90","cos(theta)_e[0]_costheta_-0.00_phi_-90","cos(theta)_e[0]_costheta_0.33_phi_-90","cos(theta)_e[0]_costheta_0.67_phi_-90",
		"cos(theta)_e[0]_costheta_-1.00_phi_-60","cos(theta)_e[0]_costheta_-0.67_phi_-60","cos(theta)_e[0]_costheta_-0.33_phi_-60","cos(theta)_e[0]_costheta_-0.00_phi_-60","cos(theta)_e[0]_costheta_0.33_phi_-60","cos(theta)_e[0]_costheta_0.67_phi_-60",
		"cos(theta)_e[0]_costheta_-1.00_phi_-30","cos(theta)_e[0]_costheta_-0.67_phi_-30","cos(theta)_e[0]_costheta_-0.33_phi_-30","cos(theta)_e[0]_costheta_-0.00_phi_-30","cos(theta)_e[0]_costheta_0.33_phi_-30","cos(theta)_e[0]_costheta_0.67_phi_-30",
		"cos(theta)_e[0]_costheta_-1.00_phi_0","cos(theta)_e[0]_costheta_-0.67_phi_0","cos(theta)_e[0]_costheta_-0.33_phi_0","cos(theta)_e[0]_costheta_-0.00_phi_0","cos(theta)_e[0]_costheta_0.33_phi_0","cos(theta)_e[0]_costheta_0.67_phi_0",
		"cos(theta)_e[0]_costheta_-1.00_phi_30","cos(theta)_e[0]_costheta_-0.67_phi_30","cos(theta)_e[0]_costheta_-0.33_phi_30","cos(theta)_e[0]_costheta_-0.00_phi_30","cos(theta)_e[0]_costheta_0.33_phi_30","cos(theta)_e[0]_costheta_0.67_phi_30",
		"cos(theta)_e[0]_costheta_-1.00_phi_60","cos(theta)_e[0]_costheta_-0.67_phi_60","cos(theta)_e[0]_costheta_-0.33_phi_60","cos(theta)_e[0]_costheta_-0.00_phi_60","cos(theta)_e[0]_costheta_0.33_phi_60","cos(theta)_e[0]_costheta_0.67_phi_60",
		"cos(theta)_e[0]_costheta_-1.00_phi_90","cos(theta)_e[0]_costheta_-0.67_phi_90","cos(theta)_e[0]_costheta_-0.33_phi_90","cos(theta)_e[0]_costheta_-0.00_phi_90","cos(theta)_e[0]_costheta_0.33_phi_90","cos(theta)_e[0]_costheta_0.67_phi_90",
		"cos(theta)_e[0]_costheta_-1.00_phi_120","cos(theta)_e[0]_costheta_-0.67_phi_120","cos(theta)_e[0]_costheta_-0.33_phi_120","cos(theta)_e[0]_costheta_-0.00_phi_120","cos(theta)_e[0]_costheta_0.33_phi_120","cos(theta)_e[0]_costheta_0.67_phi_120",
		"cos(theta)_e[0]_costheta_-1.00_phi_150","cos(theta)_e[0]_costheta_-0.67_phi_150","cos(theta)_e[0]_costheta_-0.33_phi_150","cos(theta)_e[0]_costheta_-0.00_phi_150","cos(theta)_e[0]_costheta_0.33_phi_150","cos(theta)_e[0]_costheta_0.67_phi_150" };

	TString histlist2D[1] = {"B1_map"};

	// canvas strig name, from Till suggestion! Not necessary now.
	const char *cname[] = {"PECD_e[0]_R_CL-CR", "R_CL-CR"};
	TCanvas *canv[10]; // dim(3 cos(theta) + 3 b1__maps)
	std::cout << "Vector of TCanvas created." << "\n \n";

	char tempj[250];
	char tempk[250];
	char tempjk[250];

	// LOAD AND NORMALIZED ALL GRAPHS
	load_TH1D(R_LFile, R_RFile, R_LCP_TH1D_Array, R_RCP_TH1D_Array, general_Array, dirlist, histlist1D, dir);
	std::cout << "load_all TH1D : DONE! \n \n";
	int ratio_counter = 0;

	//generation of phi and cos(theta) 16 bins: be sure that the increment and the incremented variabel is the same at in histlist1D 
	int counter = 0;
	// 12x6
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 6; j++) {
			phi[counter] = i*30-180+15;
			cos_theta[counter] = j*0.333-1+0.1667; // (bin-period/2+bin/2)
			counter++;
		}
	}

	//1 canvas of i*j plots each cos(theta) + 1 cavas for b1_map
	for (int j = 0; j < 2; j++) {
		if (j == 0) {
			canv[j] = new TCanvas(cname[j], cname[j]);
			canv[j]->SetCanvasSize(1920, 1080);
			canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[j]->SetBatch(kTRUE); // to suppress graphical putput
			canv[j]->Divide(12, 6);
			canv[j]->SetFillStyle(4000);

			for (int i = 0; i < 72; i++) {
				sprintf(tempj, "clean");
				sprintf(tempk, "clean");
				sprintf(tempj, "%s", R_LCP_TH1D_Array[i]->GetTitle());
				sprintf(tempk, "%s", R_RCP_TH1D_Array[i]->GetTitle());

				// tempj and tempk taken from Working_Array_2
				sprintf(tempjk, "PECD %s vs. %s - cos(theta) = %2.2f ; phi = %3.0f", tempj, tempk, cos_theta[i], phi[i]);
				TPaveStats *sk = (TPaveStats*)R_LCP_TH1D_Array[i]->GetListOfFunctions()->FindObject("stats");
				sk->SetOptStat(11);

				// ASYMMETY plots
				ratio_genarray[ratio_counter] = ratio_TH1D(R_LCP_TH1D_Array[i], R_RCP_TH1D_Array[i], cname[j]);
				ratio_genarray[ratio_counter]->Fit("pol1", "Q");
				// console output could be annoying: put here a command ti suppress it
				fit1[i] = ratio_genarray[ratio_counter]->GetFunction("pol1");
				b1[ratio_counter] = fit1[i]->GetParameter(1);
				b1_err[ratio_counter] = fit1[i]->GetParError(1);

				if (fp!=NULL) {
					fprintf(fp,"%4.2f \t %4.0f \t %4.8f \t %4.8f \n", cos_theta[i], phi[i], b1[i], b1_err[i]);
				}

				//plotting B1 with regression curve
				ratio_genarray[ratio_counter]->SetMarkerStyle(20);
				ratio_genarray[ratio_counter]->GetXaxis()->SetTitle("cos(theta)");
				ratio_genarray[ratio_counter]->GetYaxis()->SetTitle("asymm [a.u.]");
				ratio_genarray[ratio_counter]->GetYaxis()->SetRangeUser(0.4, 1.6);
				ratio_genarray[ratio_counter]->SetTitle("");
				canv[j]->cd(i+1);
				ratio_genarray[ratio_counter]->Draw("AP");
				fit1[i]->Draw("SAME");

				DrawTextInPad(0.2, 0.93, tempjk, 0.05);
				gStyle->SetOptTitle(0);
				gStyle->SetOptFit();	
				canv[j]->Update();

				//// Set the leggend
				//TLegend *leg = canv[i]->cd(counter)->BuildLegend();
				//leg->SetTextFont(30);
				//leg->SetY1NDC(2);
				ratio_counter++;
			}// end of for i
			canv[j]->Write();
			canv[j]->SetBatch(kFALSE); // to suppress graphical putput
			printf("plotting and saving cos(theta) canvas %s: DONE! \n \n", cname[j]);
			//cout << "test1 \n";
		}	
		if (j == 1) {
			cout << "test2 \n";
			sprintf(tempjk, "b1 map of %s", cname[j]);
			canv[j] = new TCanvas(cname[j], cname[j]);
			//canv[j]->SetCanvasSize(1920, 1080);
			//canv[j]->SetWindowSize(1920 + 4, 1080 + 28);
			canv[j]->SetFillStyle(4000);
			canv[j]->cd();

			b1_map[j-1] = new TGraph2D();
			b1_map[j-1]->SetNpx(12);
			b1_map[j-1]->SetNpy(6);			
			b1_map[j-1]->SetTitle("; phi [DEG]; cos(theta) [adm]; b1 [adm]");
			canv[j]->Update();
			////filling b1_maps
			//for (int k = 0; k < 16; k++) //k = all phi and cos(theta) = 4*4
			//	b1_map[j-3]->SetPoint(k, phi[k], cos_theta[k], b1[k+(j-3)*16]);

			for (int k = 0; k < 72; k++) //k = all phi and cos(theta) = 9x8
				b1_map[j-1]->SetPoint(k, phi[k], cos_theta[k], b1[k]);
			b1_map[j-1]->GetXaxis()->SetRangeUser(-180., 180.); // doesn't work with TGraph2D
			b1_map[j-1]->GetYaxis()->SetRangeUser(-1., 1.);		// doesn't work with TGraph2D
			b1_map[j-1]->GetZaxis()->SetRangeUser(-0.15, 0.15); // doesn't work with TGraph2D
			b1_map[j-1]->Draw("COLZ"); // this will introduce directly the smoothing
			DrawTextInPad(0.32, 0.93, tempjk, 0.05); // final touches: title and legend
			//TLegend *leg = canv[j]->cd()->BuildLegend();
			//leg->SetTextFont(30);
			//leg->SetY1NDC(2);
			gStyle->SetOptTitle(0);
			canv[j]->Update();

			canv[j]->Write();
			printf("plotting and saving Beta1_maps canvas %s: DONE! \n", cname[j]);
		}
	}

	f->Close();
	std::cout << "gPad : SAVED!" << "\n \n";

	// hard exit with no delete of variable
	exit(0);
	theApp.Run();

	return 0;
}
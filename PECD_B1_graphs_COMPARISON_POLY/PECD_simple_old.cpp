// PECD_graphs.cpp : Defines the entry point for the console application.
//
#define _CRT_SECURE_NO_WARNINGS
#include "stdafx.h"

#include <TApplication.h>

#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TPaveLabel.h>
#include <THistPainter.h>
#include <TText.h>
#include <TPad.h>
#include <TStyle.h>

#include <TRandom.h>

#include <TF1.h>



int main(int argc, char** argv)

{

	int binx = 0;
	int biny = 0;
	int a = 0;
	double integral = 0.;
	double asym_hist1_counts = 0.;
	double asym_hist2_counts = 0.;


	TApplication theApp("App",&argc,argv);

	//LCP loading
	TFile *HFile1 = new TFile("E:/Soleil_2018_09_SEX/Analysis/R_MetOx_33eV_LCP_m30-28_invalid_mod.root", "READ");
	//HFile->ls();
	HFile1->cd("angular_distr_el/if_twoelec/NO_gate");
	TH1D *h1a = (TH1D*)HFile1->FindObjectAny("cos(theta)_e[e1]");
	TH1D *h1b = (TH1D*)HFile1->FindObjectAny("cos(theta)_e[e2]");
	TH1D *h1c = (TH1D*)HFile1->FindObjectAny("cos(theta)_pesum");
	TH2D *h2a = (TH2D*)HFile1->FindObjectAny("PECD_e[e1](theta)");
	TH2D *h2b = (TH2D*)HFile1->FindObjectAny("PECD_e[e2](theta)");
	TH2D *h2c = (TH2D*)HFile1->FindObjectAny("PECD_pesum(theta)");
	gDirectory->pwd();
	//h1a->ls();
	//h1b->ls();
	//h1c->ls();
	//h2a->ls();
	//h2b->ls();
	//h2c->ls();

	//RCP loading
	TFile *HFile2 = new TFile("E:/Soleil_2018_09_SEX/Analysis/R_MetOx_33eV_RCP_m30-28_invalid_mod.root");
	HFile2->cd("angular_distr_el/if_twoelec/NO_gate");

	TH1D *h3a = (TH1D*)HFile2->FindObjectAny("cos(theta)_e[e1]");
	TH1D *h3b = (TH1D*)HFile2->FindObjectAny("cos(theta)_e[e2]");
	TH1D *h3c = (TH1D*)HFile2->FindObjectAny("cos(theta)_pesum");
	TH2D *h4a = (TH2D*)HFile2->FindObjectAny("PECD_e[e1](theta)");
	TH2D *h4b = (TH2D*)HFile2->FindObjectAny("PECD_e[e2](theta)");
	TH2D *h4c = (TH2D*)HFile2->FindObjectAny("PECD_pesum(theta)");
	gDirectory->pwd();
	//h3a->ls();
	//h3b->ls();
	//h3c->ls();
	//h4a->ls();
	//h4b->ls();
	//h4c->ls();

	//TH2D *h7= PECD_e[e1](theta);
	//TH2D *h8= PECD_e[e1](theta);

	h1a->Sumw2();
	h1b->Sumw2();
	h1c->Sumw2();
	//h2a->Sumw2();
	//h2b->Sumw2();
	//h2c->Sumw2();
	h3a->Sumw2();
	h3b->Sumw2();
	h3c->Sumw2();
	//h4a->Sumw2();
	//h4b->Sumw2();
	//h4c->Sumw2();

	//Normalize TH1D graphs
	h1a->Scale(1./h1a->Integral()); 
	h1b->Scale(1./h1b->Integral());
	h1c->Scale(1./h1c->Integral());
	h3a->Scale(1./h3a->Integral());
	h3b->Scale(1./h3b->Integral());
	h3c->Scale(1./h3c->Integral());

	//normalizing TH2D graphs
	//h2a
	binx = h2a->GetXaxis()->GetNbins();
	biny = h2a->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h2a->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h2a->SetBinContent(x,y,h2a->GetBinContent(x,y)/integral);
	}

	//h2b
	binx = h2b->GetXaxis()->GetNbins();
	biny = h2b->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		double integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h2b->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h2b->SetBinContent(x,y,h2b->GetBinContent(x,y)/integral);
	}

	//h2c
	binx = h2c->GetXaxis()->GetNbins();
	biny = h2c->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h2c->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h2c->SetBinContent(x,y,h2c->GetBinContent(x,y)/integral);
	}

	//h4a
	binx = h4a->GetXaxis()->GetNbins();
	biny = h4a->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h4a->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h4a->SetBinContent(x,y,h4a->GetBinContent(x,y)/integral);
	}

	//h4b
	binx = h4b->GetXaxis()->GetNbins();
	biny = h4b->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h4b->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h4b->SetBinContent(x,y,h4b->GetBinContent(x,y)/integral);
	}

	//h4c
	binx = h4c->GetXaxis()->GetNbins();
	biny = h4c->GetYaxis()->GetNbins();

	for (int x = 1;x<=binx;++x) {
		integral = 0.;
		for (int y = 1;y<=biny;++y) integral += h4c->GetBinContent(x,y);
		if (fabs(integral)> 1.e-100) for (int y = 1;y<=biny;++y) h4c->SetBinContent(x,y,h4c->GetBinContent(x,y)/integral);
	}
	//h1->ResetStats();
	//h2->ResetStats();
	//h5->ResetStats();
	//h6->ResetStats();

	//create canvas
	TCanvas* c1 = new TCanvas(gDirectory->GetName(),HFile1->GetName());
	c1->Divide(2,3);
	c1->cd();
	c1->SetFillStyle(4000);
	c1->Draw();
	//Get text from the HFile1
	TText *t = new TText(.45,.95,gDirectory->GetName());
	t->SetTextAlign();
	t->Draw();


	//TCanvas *c1 = new TCanvas("c1",HFile1->GetName()); 
	//gStyle->SetOptStat(1);
	//c1->Divide(2,3);
	//gStyle->SetOptTitle(0);
	//TPaveLabel *title = new TPaveLabel(11,95,35,99,"aosjdasjdas","brndc");
	//title->Draw();
	//c1->Update();

	// Asymmetry = (h1 - h2)/(h1 + h2)  where h1 = this
	c1->cd(1);
	c1->cd(1)->SetFillStyle(4000);
	c1->cd(1)->SetFrameFillColor(4000);
	//gPad->SetFillStyle(4000);
	//c1->SetFrameFillColor(0);
	a = h1a->GetNbinsX();
	auto ratio1 = new TGraphAsymmErrors(a);
	ratio1->Divide(h1a,h3a,"pois midp");
	ratio1->SetMarkerStyle(20);
	ratio1->GetXaxis()->SetTitle("cos(theta)");
	ratio1->GetYaxis()->SetTitle("counts [a.u.]");
	ratio1->SetTitle("Histogram title");
	ratio1->Draw("AP");
	//gPad->Update();
	//ratio1->Draw("HIST E");
	//ratio1->SetMarkerStyle(20); 
	//ratio1->Draw("SAME E"); 

	c1->cd(3);
	a = h1b->GetNbinsX();
	auto ratio2 = new TGraphAsymmErrors(a);
	ratio2->Divide(h1b,h3b,"pois midp"); 
	ratio2->SetMarkerStyle(20);
	ratio2->GetXaxis()->SetTitle("cos(theta)");
	ratio2->GetYaxis()->SetTitle("counts [a.u.]");
	ratio2->SetTitle("Histogram title");
	ratio2->Draw("AP");
	//ratio2->Draw("HIST E"); 
	//ratio2->SetMarkerStyle(20); 
	//ratio2->Draw("SAME E"); 

	c1->cd(5);
	a = h1c->GetNbinsX();
	auto ratio3 = new TGraphAsymmErrors(a);
	ratio3->Divide(h1c,h3c,"pois midp");
	ratio3->SetMarkerStyle(20);
	ratio3->GetXaxis()->SetTitle("cos(theta)");
	ratio3->GetYaxis()->SetTitle("counts [a.u.]");
	ratio3->SetTitle("Histogram title");
	ratio3->Draw("AP");
	//gPad->Update();
	//ratio1->Draw("HIST E");
	//ratio1->SetMarkerStyle(20); 
	//ratio1->Draw("SAME E"); 

	// alternative merhod
	//c1->cd(2,1);
	//TH2D *hsum = new TH2D("Sum",h2b->GetName(),h2a->GetNbinsX(),h2a->GetXaxis()->GetXmin(),h2a->GetXaxis()->GetXmax(),h2a->GetNbinsY(),h2a->GetYaxis()->GetXmin(),h2a->GetYaxis()->GetXmax());
	//TH2D *hdiff= new TH2D("Diff",h2b->GetName(),h2a->GetNbinsX(),h2a->GetXaxis()->GetXmin(),h2a->GetXaxis()->GetXmax(),h2a->GetNbinsY(),h2a->GetYaxis()->GetXmin(),h2a->GetYaxis()->GetXmax());
	////hsum->SetOption("colz");
	//hsum->Add(h2a,h2b,1,1);
	//hdiff->Add(h2a,h2b,1,-1);
	//hdiff->Divide(hsum);
	////this draw option makes it identical to cd(4)
	////hdiff->Draw("COLZ");
	//hdiff->Draw();

	c1->cd(2);
	TH2D* hdiff_n1 = (TH2D*)h2a->GetAsymmetry(h4a);
	hdiff_n1->Draw("COLZ2");

	c1->cd(4);
	TH2D* hdiff_n2 = (TH2D*)h2b->GetAsymmetry(h4b);
	hdiff_n2->Draw("COLZ2");

	// // to be cleaned
	// asym_hist1_counts = h2b->GetEntries();
	// asym_hist2_counts = h4b->GetEntries();
	// TH2D* hdiff_n1 = (TH2D*)h2b->GetAsymmetry(h4b, asym_hist1_counts/asym_hist2_counts);
	// hdiff_n1->Draw("COLZ2");

	c1->cd(6);
	TH2D* hdiff_n3 = (TH2D*)h2c->GetAsymmetry(h4c);
	hdiff_n3->Draw("COLZ2");

	//// to be cleaned
	//asym_hist1_counts = h2c->GetEntries();
	//asym_hist2_counts = h4c->GetEntries();
	//TH2D* hdiff_n2 = (TH2D*)h2c->GetAsymmetry(h4c, asym_hist1_counts/asym_hist2_counts);
	//hdiff_n2->Draw("COLZ2");

	theApp.Run();

	return 0;

}
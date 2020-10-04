// Author: Joshua Williams

#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include <vector>


//Takes at set of tof and momentum points as imput and returns a TF1 fuctions that is fitted to these points.
//This TF1 can then be used to convert tof to momentum.
//It uses a polynomial function A*x^4+B*x^3+C*x^2+D*x+C, where A,B,C,D,E are fitted parameters.
TF1 *  Tof_to_P(char * name, int num_data_points, double *tofs, double *moms)
{
	TCanvas *c1 = new TCanvas(name,name,200,10,700,500);
    c1->SetGrid();
	gPad->SetRightMargin(0.12);
	gPad->SetLeftMargin(0.11);
	gPad->SetTopMargin(0.0825);
	gPad->SetBottomMargin(0.087);
	
	TGraph * gr1 = new TGraph(num_data_points, tofs, moms);
	gr1->SetMarkerColor(kBlue);
	gr1->SetMarkerStyle(21);
	gr1->SetLineWidth(0);
	gr1->SetLineColor(0);
	gr1->Draw();	
	

	double uxmin = gPad->GetUxmin();
	double uxmax = gPad->GetUxmax();
	
	char fitname[256];
	sprintf(fitname,"%s-myfit",name);
	
	TF1 * myfit = new TF1(fitname,"[4]*x**4+[3]*x**3+[2]*x**2+[1]*x+[0]",uxmin, uxmax);
	myfit->SetParName(0,"E");
	myfit->SetParName(1,"D");
	myfit->SetParName(2,"C");
	myfit->SetParName(3,"B");
	myfit->SetParName(4,"A");
	myfit->SetParameter(0, 0);
	myfit->SetParameter(1, 0);
	myfit->SetParameter(2, 0);
	myfit->SetParameter(3, 0);
	myfit->SetParameter(4, 0);
	gr1->Fit(fitname,"M");


	gPad->Modified();
	gPad->Update();
	gPad->Modified();
	gPad->Update();
	
	return myfit;
	
}











// rootstuff.cpp: implementation of the rootstuff class.
//////////////////////////////////////////////////////////////////////

#pragma warning(disable : 4996)
#pragma warning(disable : 4800)

/////////////////////////////////////////////////////////////////////////

// Project/Settings/C++/Precompiled Headers/Not using precompiled Headers !!!!

/////////////////////////////////////////////////////////////////////////
#include "OS_Version.h"

#include "stdlib.h"
#include "rootstuff.h"
//#include "TApplication.h"
//#include "TCanvas.h"
//#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
//#include "TMinuit.h"

//#include "TBrowser.h"
#include "TNtuple.h"
#include "TColor.h"

//#include "afx.h"
//#include "TKey.h" //Must be after afx include since they both have a GetClassName member!

#ifdef _DEBUG
	#include "assert.h"
#endif


////////////////////////////////////////////////////////////////
//For ROOT//////////////////////////////////////////////////////
//#include "TROOT.h"
//extern void InitGui();
//VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
//TROOT root("rint", "The ROOT Interactive Interface", initfuncs);
////////////////////////////////////////////////////////////////



rootstuff::rootstuff()
{
}



rootstuff::~rootstuff()
{
}



TTree * rootstuff::OpenRootFileGetTree(const char *name,const char *TreeName)
{
    TFile * m_RootFile = new TFile(name,"READ");

	if (m_RootFile->IsZombie()) {
		return 0;
	}

    TTree * m_RootTree = (TTree*)m_RootFile->Get(TreeName);
	return m_RootTree;
}


TFile * rootstuff::RecreateRootFile(const char *name,const char *comment)
{
    TFile * m_RootFile = new TFile(name,"RECREATE",comment);
	return m_RootFile;
}

TCanvas * rootstuff::newCanvas(char *name,char *titel,__int32 xposition,__int32 yposition,__int32 pixelsx,__int32 pixelsy)
{
	TCanvas * canvaspointer;
	canvaspointer = new TCanvas(name,titel,xposition,yposition,pixelsx,pixelsy);
	return canvaspointer;
}

TNtuple * rootstuff::newNTuple(char *name, char * title, char *varlist, __int32 buffersize)
{
   TNtuple * nTuple = new TNtuple(name,title,varlist,buffersize);
   

   return nTuple;
}




TH1D * rootstuff::newTH1D(__int32 number,char *comment,__int32 bins,double xmin,double xmax,char *titelx,char *option)
{
   char name[200];
   TH1D * hist;
   sprintf(name,"%i",number);
   hist = new TH1D(name,comment,bins,xmin,xmax);
   hist->SetOption(option);
   hist->GetXaxis()->SetTitle(titelx);
   return hist;
}



TH1D * rootstuff::newTH1D(char *name,char *comment,__int32 bins,double xmin,double xmax,char *titelx,char *option)
{
   TH1D * hist;
   hist = new TH1D(name,comment,bins,xmin,xmax);
   hist->SetOption(option);
   hist->GetXaxis()->SetTitle(titelx);
   return hist;
}


TH2D * rootstuff::newTH2D(__int32 number,char *comment,__int32 xbins,double xmin,double xmax,char *titelx,__int32 ybins,double ymin,double ymax,char *titely,char *option)
{
   char name[200];
   TH2D * hist;
   sprintf(name,"%i",number);
   hist = new TH2D(name,comment,xbins,xmin,xmax,ybins,ymin,ymax);
   hist->SetOption(option);
   hist->GetXaxis()->SetTitle(titelx);
   hist->GetYaxis()->SetTitle(titely);
   return hist;
}


TH2D * rootstuff::newTH2D(char *name,char *comment,__int32 xbins,double xmin,double xmax,char *titelx,__int32 ybins,double ymin,double ymax,char *titely,char *option)
{
   TH2D * hist;
   hist = new TH2D(name,comment,xbins,xmin,xmax,ybins,ymin,ymax);
   hist->SetOption(option);
   hist->GetXaxis()->SetTitle(titelx);
   hist->GetYaxis()->SetTitle(titely);
   return hist;
}

/////////////////////////////////////////ColorScales/////////////////////////////////////////////////////////////////////////
void rootstuff::color()
{ 
   
	//tills palette
	UInt_t Number = 4;
	Double_t Red[4]   = { 0.50, 0.00, 0.99, 0.99};
	Double_t Green[4] = { 0.99, 0.00, 0.00, 0.99};
	Double_t Blue[4]  = { 0.99, 0.99, 0.00, 0.00};
	Double_t Stops[4] = { 0.00, 0.30, 0.70, 1.00};
	TColor * c = new TColor();
	c->CreateGradientColorTable(Number,Stops,Red,Green,Blue,255);
	if (c) delete c;
	c=0;
	
	
	/*//tills std. palette
	UInt_t Number = 5;
	Double_t Red[5]   = { 0.00, 0.23, 0.99, 0.99, 0.99};
	Double_t Green[5] = { 0.00, 0.00, 0.00, 0.99, 0.99};
	Double_t Blue[5]  = { 0.30, 0.70, 0.00, 0.00, 0.99};
	Double_t Stops[5] = { 0.00, 0.10, 0.33, 0.70, 1.00};
	TColor * c = new TColor();
	c->CreateGradientColorTable(Number,Stops,Red,Green,Blue,255);
	delete c;
	//gStyle->SetFrameFillColor(229); 
	*/
	
}

TDirectory* rootstuff::getDir(TFile *rootfile, TString dirName)
{
	//first find out whether directory exists
#ifdef _DEBUG
	assert(rootfile);
#endif
	//printf("here1: %d \n",rootfile);
	if (!rootfile) return 0;
	rootfile->cd("/"); 
	TDirectory * direc = rootfile->GetDirectory(dirName.Data()); 
	if (!direc)
	{
		//direc->SetWritable(true);
		//if not create it//
		TString lhs;
		TString rhs;
		TString tmp = dirName;
		while (1)
		{
			//if there is no / then this is the last subdir
			if (tmp.Index("/") == -1)
			{
				lhs = tmp;
			}
			else //otherwise split the string to lefthandside and righthandside of "/"
			{
				lhs = tmp(0,tmp.Index("/"));
				rhs = tmp(tmp.Index("/")+1,tmp.Length());
			}

			//check wether subdir exits//
			direc = gDirectory->GetDirectory(lhs.Data());
			if (direc){
				gDirectory->cd(lhs.Data());//cd into it
			}
			else
			{
				direc = gDirectory->mkdir(lhs.Data()); //create it
				gDirectory->cd(lhs.Data()); //and cd into it
			}

			//when there is no "/" anymore break here//
			if (tmp.Index("/") == -1)
				break;

			//the new temp is all that is on the right hand side
			tmp = rhs;
		}
	}
	//return to root Path//
	rootfile->cd("/");
	return direc;
}

void rootstuff::add_hist(CH::H1d * hist, TFile * RootFile){
	//change/make directory 
	//TString dir = hist->get_dir();
	//cout <<hist->get_name()<< ", " << hist->get_dir()<<endl;
	Histograms_dir = getDir(RootFile, hist->get_dir() );
	Histograms_dir->cd();

	// create root hist
//	TH1d * root_hist =new TH1d(hist->get_name().c_str(), hist->get_title().c_str(), hist->get_X_n_bins(), hist->get_X_min(), hist->get_X_max() );
	TH1D * root_hist =new TH1D(hist->get_name(), hist->get_title(), hist->get_X_n_bins(), hist->get_X_min(), hist->get_X_max() );

	// copy contents
	
	root_hist->SetBinContent( 0 , hist->get_X_underflow() );
	for(int i=0; i < hist->get_X_n_bins(); ++i){
		root_hist->SetBinContent( i+1, hist->bins[i]);
	}
	root_hist->SetBinContent( hist->get_X_n_bins()+1 , hist->get_X_overflow() );
	//set axis title
	root_hist->SetXTitle( hist->get_X_title() );
	root_hist->SetEntries(hist->get_Entries());

	//write to root file
	// No need to write it, obviously new TH1D already creates the hist in the file.. 
	//root_hist->Write(hist->get_name() );
}



void rootstuff::add_hist(CH::H2d * hist, TFile * RootFile){
	//change/make directory 
	Histograms_dir = getDir(RootFile, hist->get_dir() );
	Histograms_dir->cd();

	// create root hist
	TH2D * root_hist =new TH2D(hist->get_name(), hist->get_title(), hist->get_X_n_bins(), hist->get_X_min(), hist->get_X_max(), hist->get_Y_n_bins(), hist->get_Y_min(), hist->get_Y_max() );
	// copy contents

	for(int i=0; i < hist->get_X_n_bins(); ++i){
		for(int j=0; j < hist->get_Y_n_bins(); ++j){
			root_hist->SetBinContent( i+1, j+1, hist->bins[i][j]);
		}
	}
	root_hist->SetBinContent( hist->get_X_n_bins()+1, 1 , hist->get_X_overflow() );
	root_hist->SetBinContent( 1, hist->get_Y_n_bins()+1 , hist->get_Y_overflow() );
	
	root_hist->SetBinContent( hist->get_X_n_bins()+1, hist->get_Y_n_bins()+1 , hist->get_overflow_xy() );
	root_hist->SetBinContent( 0, 0 , hist->get_underflow_xy() );

	root_hist->SetBinContent( 0, hist->get_Y_n_bins()+1 , hist->get_flow_under_x_over_y() );
	root_hist->SetBinContent( hist->get_X_n_bins()+1, 0, hist->get_flow_over_x_under_y() );

	root_hist->SetBinContent( 0, 1, hist->get_X_underflow() );
	root_hist->SetBinContent( 1, 0, hist->get_Y_underflow() );


	//set axis title
	root_hist->SetXTitle( hist->get_X_title() );
	root_hist->SetYTitle( hist->get_Y_title() );
	root_hist->SetEntries(hist->get_Entries());

	root_hist->SetOption("colz");
	//write to root file
	//root_hist->Write(hist->get_name() );
}

void rootstuff::add_hist(CH::H3d * hist, TFile * RootFile){
	//change/make directory 
	Histograms_dir = getDir(RootFile, hist->get_dir() );
	Histograms_dir->cd();
	
	// create root hist
	TH3D* root_hist = new TH3D(hist->get_name(), hist->get_title(), hist->get_X_n_bins(), hist->get_X_min(), hist->get_X_max(), hist->get_Y_n_bins(), hist->get_Y_min(), hist->get_Y_max(), hist->get_Z_n_bins(), hist->get_Z_min(), hist->get_Z_max() );
	
	// copy contents
	//root_hist->SetBinContent( root_hist->GetBin(0,1,1), hist->get_X_underflow() );
	//root_hist->SetBinContent( root_hist->GetBin(1,0,1), hist->get_Y_underflow() );
	//root_hist->SetBinContent( root_hist->GetBin(1,1,0), hist->get_Z_underflow() );

	for(int i=0; i < hist->get_X_n_bins(); ++i){
		for(int j=0; j < hist->get_Y_n_bins(); ++j){
			for(int k=0; k < hist->get_Z_n_bins(); ++k){
				root_hist->SetBinContent( root_hist->GetBin(i+1, j+1, k+1), hist->bins[i][j][k]);
			}
		}
	}
	//root_hist->SetBinContent( root_hist->GetBin(hist->get_X_n_bins()+1, 1, 1), hist->get_X_overflow() );
	//root_hist->SetBinContent( root_hist->GetBin(1, hist->get_Y_n_bins()+1, 1), hist->get_Y_overflow() );
	//root_hist->SetBinContent( root_hist->GetBin(1, 1, hist->get_Z_n_bins()+1), hist->get_Z_overflow() );

	//set axis title
	root_hist->SetXTitle( hist->get_X_title() );
	root_hist->SetYTitle( hist->get_Y_title() );
	root_hist->SetZTitle( hist->get_Z_title() );
	
	root_hist->SetEntries(hist->get_Entries());


	root_hist->SetOption("glbox");
	root_hist->SetFillColor(kBlue);

	//write to root file
	//root_hist->Write(hist->get_name() );
} 
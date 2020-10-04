//////////////////////////////////////////////////////////////////////
// rootstuff.h: interface for the rootstuff class.
//////////////////////////////////////////////////////////////////////
#pragma once

#include "OS_Version.h"

#include "TFile.h"
#include "TTree.h"

#include "../ColAHelL/src/CH_Hist.h"

#ifndef _ROOTSTUFF_
	#define _ROOTSTUFF_

	class TCanvas;
	class TFile;
	//class TMinuit;
	class TH1D;
	class TH2D;
	class TNtuple;

	class TTree;
	class TMarker;

	class rootstuff
	{

	private:

		TDirectory* getDir(TFile *rootfile, TString dirName);
		TDirectory * Histograms_dir;

	public:
		TTree * OpenRootFileGetTree(const char*,const char*);
		TFile * RecreateRootFile(const char*,const char*);

		TH1D * newTH1D(char*,char*,__int32,double,double,char*,char*);
		TCanvas * newCanvas(char*,char*,__int32,int,__int32,__int32);
		TH1D * newTH1D(int,char*,__int32,double,double,char*,char*);
		TH2D * newTH2D(char*,char*,__int32,double,double,char*,__int32,double,double,char*,char*);
		TH2D * newTH2D(int,char*,__int32,double,double,char*,__int32,double,double,char*,char*);
		TNtuple * newNTuple(char*,char*,char*,__int32);

		void color();

		rootstuff();
		~rootstuff();

		void add_hist(CH::H1d * hist, TFile * RootFile);
		void add_hist(CH::H2d * hist, TFile * RootFile);
		void add_hist(CH::H3d * hist, TFile * RootFile);

	};

	#define root1(name,posx,weight,explanation,bins,min,max,titelx)                                {static TH1D * name; if (!name || Ueber->start_new_root_file) name=rt->newTH1D(#name,explanation,__int32(bins+0.0000001),double(min),double(max),titelx,"");                                                               name->Fill(posx,weight);}
	#define root2(name,posx,posy,weight,explanation,binsx,minx,maxx,titelx,binsy,miny,maxy,titely) {static TH2D * name; if (!name || Ueber->start_new_root_file) name=rt->newTH2D(#name,explanation,__int32(binsx+0.0000001),double(minx),double(maxx),titelx,__int32(binsy+0.0000001),double(miny),double(maxy),titely,"COLZ");  name->Fill(posx,posy,weight);}

	#define rootNTuple(name,title,varlist,buffersize,data) static TNtuple * name; if (!name || Ueber->start_new_root_file) name=rt->newNTuple(#name,title,varlist,__int32(buffersize+0.00001)); name->Fill(data)
#endif

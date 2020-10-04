

#include "OS_Version.h"
#pragma warning(disable : 4996)
#include <TFile.h>

#include <TH3.h>
#include <TObject.h>
#include <TDirectory.h>
#include <iostream>
#include "TNtupleD.h"
#include "Histo.h"




//_______________________________________________________________________________________________________
Histo::Histo(__int32 size)
{
	arsize = size;
	tuple_size = 100;

	hiarray = new TObject*[arsize];
	hi_info_array = new hist_info*[arsize];
	for (__int32 i =0;i<arsize;++i) {
		hiarray[i]=0;
		hi_info_array[i] = 0;
	}

	tuple_array = new TObject*[tuple_size];
	for (__int32 i =0;i<tuple_size;++i) tuple_array[i] = 0;

	rootfile = 0;

	verbose = false;
}

//_______________________________________________________________________________________________________
Histo::~Histo()
{
	//--delete all histos--//
	Reset();

	delete [] hiarray;
	delete [] tuple_array;
	delete [] hi_info_array;

	//--Close root file--//
	if (rootfile) rootfile->Close();
	rootfile = 0;
}


//########################ROOT File Stuff################################################################
//_______________________________________________________________________________________________________
void Histo::openRootFile(const char *name)
{
	rootfile = TFile::Open(name,"RECREATE");
	std::cout << "Histograms will be written to: "<<name<<std::endl;
}

//_______________________________________________________________________________________________________
void Histo::flushRootFile()
{
	//--write histos to directory--//
	for (__int32 i=0;i<arsize;++i)
	{
		TObject * obj = hiarray[i];
		if (obj)
		{
			//--First change to the right directory--//
			TH1* h = (TH1*)obj;
			h->GetDirectory()->cd();
			//std::cout<<"Writing Histogram:"<<h->GetName()<<" Dir:"<<h->GetDirectory()->GetName()<<std::endl;
			//--write object to current directory--//
			obj->Write(0,TObject::kOverwrite);
			if (rootfile) rootfile->cd("/");
		}
	}

	//--save directory to root file--//
	if (rootfile) rootfile->SaveSelf();
}



//########################ROOT Histograms Stuff################################################################



void Histo::Reset()
{
	//--delete all histograms--//
	for (__int32 i=0;i<arsize;++i)
	{
		TObject * obj = hiarray[i];
		if (obj) delete obj;
		obj = 0;
		hiarray[i] = 0;
		if (hi_info_array[i]) {delete hi_info_array[i]; hi_info_array[i] = 0;}
	}

	for (__int32 i=0;i<tuple_size;++i)
	{
		TObject * obj = tuple_array[i];
		if (obj) delete obj;
		obj = 0;
		tuple_array[i] = 0;
	}
}

TDirectory* getDir(TFile *rootfile, TString dirName)
{
	//first find out whether directory exists
#ifdef _DEBUG
	assert(rootfile);
#endif

	if (!rootfile) return 0;
	rootfile->cd("/");
	TDirectory * direc = rootfile->GetDirectory(dirName.Data());
	if (!direc)
	{
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
			else //otherwise split the string to left hand side and right hand side of "/"
			{
				lhs = tmp(0,tmp.Index("/"));
				rhs = tmp(tmp.Index("/")+1,tmp.Length());
			}

			//check whether subdir exits//
			direc = gDirectory->GetDirectory(lhs.Data());
			if (direc)
				gDirectory->cd(lhs.Data());//cd into it
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




void Histo::SetBinContentAt(__int32 id, double x, double content)
{
	#ifdef _DEBUG
		assert(id >= 0);
		assert(id < arsize);
		assert(hiarray[id]);
	#endif
	TH1D * hist = (TH1D*)hiarray[id];
#ifdef _DEBUG
	assert(hist->GetDimension() == 1);
#endif
	__int32 binx = hist->GetXaxis()->FindBin(x);
	hist->SetBinContent(binx,content);
}





void Histo::SetBinContentAt(__int32 id, double x, double y, double content)
{
	#ifdef _DEBUG
		assert(id >= 0);
		assert(id < arsize);
		assert(hiarray[id]);
	#endif
	TH2D * hist = (TH2D*)hiarray[id];
#ifdef _DEBUG
	assert(hist->GetDimension() == 2);
#endif
	__int32 binx = hist->GetXaxis()->FindBin(x);
	__int32 biny = hist->GetXaxis()->FindBin(y);
	hist->SetBinContent(binx,biny,content);
}




double Histo::GetBinContentAt(__int32 id, double x)
{
	#ifdef _DEBUG
		assert(id >= 0);
		assert(id < arsize);
		assert(hiarray[id]);
	#endif
	TH1D * hist = (TH1D*)hiarray[id];
#ifdef _DEBUG
	assert(hist->GetDimension() == 1);
#endif
	__int32 binx = hist->GetXaxis()->FindBin(x);
	return hist->GetBinContent(binx);
}





double Histo::GetBinContentAt(__int32 id, double x, double y)
{
	#ifdef _DEBUG
		assert(id >= 0);
		assert(id < arsize);
		assert(hiarray[id]);
	#endif
	TH2D * hist = (TH2D*)hiarray[id];
#ifdef _DEBUG
	assert(hist->GetDimension() == 2);
#endif
	__int32 binx = hist->GetXaxis()->FindBin(x);
	__int32 biny = hist->GetXaxis()->FindBin(y);
	return hist->GetBinContent(binx,biny);
}






//___________________________________3D____________________________________________________________________
void Histo::fill(__int32 id, const char *name, double fillX, double fillY, double fillZ,
					const char *titleX, const char *titleY, const char *titleZ,
					__int32 nXbins, double xLow, double xUp,
					__int32 nYbins, double yLow, double yUp,
					__int32 nZbins, double zLow, double zUp, const char * dir,
					double weight)
{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH3D * hist = (TH3D *) hiarray[id];

	//--if the histo does not exist create it first--//
	if (!hist)
	{
		//--create a 3D histogram--//
		hist = new TH3D(name, name, nXbins, xLow, xUp, nYbins, yLow, yUp, nZbins, zLow, zUp);
		hist->SetXTitle(titleX);
		hist->GetXaxis()->CenterTitle(true);
		hist->GetXaxis()->SetTitleOffset(1.5);
		hist->SetYTitle(titleY);
		hist->GetYaxis()->CenterTitle(true);
		hist->GetYaxis()->SetTitleOffset(1.5);
		hist->SetZTitle(titleZ);
		hist->GetZaxis()->CenterTitle(true);
		hist->GetZaxis()->SetTitleOffset(1.5);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory

		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;

		if(verbose)
			std::cout <<"create 3D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 3;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
		hi_info_array[id]->nbinsy	= hist->GetNbinsY();
		hi_info_array[id]->low		= hist->GetYaxis()->GetXmin();
		hi_info_array[id]->up		= hist->GetYaxis()->GetXmax();
	}

	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;

	if (hi_info_array[id]->dim != 3) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	hist->Fill(fillX,fillY,fillZ,weight);
}




//Achim style

void Histo::fill(__int32 id, const char * name, double fillX, double fillY, double fillZ, double weight,
		  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
		  __int32 nYbins, double yLow, double yUp, const char * titleY, 
		  __int32 nZbins, double zLow, double zUp, const char * titleZ, const char * dir)

{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH3D * hist = (TH3D *) hiarray[id];

	//--if the histo does not exist create it first--//
	if (!hist)
	{
		//--create a 3D histogram--//
		hist = new TH3D(name, name, nXbins, xLow, xUp, nYbins, yLow, yUp, nZbins, zLow, zUp);
		hist->SetXTitle(titleX);
//		hist->GetXaxis()->CenterTitle(true);
		hist->GetXaxis()->SetTitleOffset(1.5);
		hist->SetYTitle(titleY);
//		hist->GetYaxis()->CenterTitle(true);
		hist->GetYaxis()->SetTitleOffset(1.5);
		hist->SetZTitle(titleZ);
//		hist->GetZaxis()->CenterTitle(true);
		hist->GetZaxis()->SetTitleOffset(1.5);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory
		hist->SetOption("glbox");
		hist->SetFillColor(kBlue);

		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;
		if(verbose)
			std::cout <<"create 3D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 3;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
		hi_info_array[id]->nbinsy	= hist->GetNbinsY();
		hi_info_array[id]->low		= hist->GetYaxis()->GetXmin();
		hi_info_array[id]->up		= hist->GetYaxis()->GetXmax();
	}

	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;

	if (hi_info_array[id]->dim != 3) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	hist->Fill(fillX,fillY,fillZ,weight);


}



void Histo::fill3(__int32 id, const char * name, double fillX, double fillY, double fillZ, double weight,
		  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
		  __int32 nYbins, double yLow, double yUp, const char * titleY, 
		  __int32 nZbins, double zLow, double zUp, const char * titleZ, const char * dir)
{
		  fill(id,name,fillX,fillY,fillZ,weight,title,nXbins,xLow,xUp,titleX,nYbins,yLow,yUp,titleY,nZbins,zLow,zUp,titleZ,dir);
}



//__________________________________2D_____________________________________________________________________
void Histo::fill(__int32 id, const char *name, double fillX, double fillY, 
				    const char *titleX, const char *titleY, 
				    __int32 nXbins, double xLow, double xUp, 
				    __int32 nYbins, double yLow, double yUp, const char * dir,
				    double weight)
{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH2D * hist = (TH2D *) hiarray[id];

	//--if the histo does not exist create it first--//
	if (!hist)
	{
		hist = new TH2D(name, name, nXbins, xLow, xUp, nYbins, yLow, yUp);
		hist->SetOption("colz");
		hist->SetXTitle(titleX);
		hist->GetXaxis()->CenterTitle(true);
		hist->SetYTitle(titleY);
		hist->GetYaxis()->CenterTitle(true);
		hist->GetYaxis()->SetTitleOffset(1.5);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory

		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;
		//--write it into rootfile--//
		if(verbose)
			std::cout <<"create 2D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 2;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
		hi_info_array[id]->nbinsy	= hist->GetNbinsY();
		hi_info_array[id]->low		= hist->GetYaxis()->GetXmin();
		hi_info_array[id]->up		= hist->GetYaxis()->GetXmax();
	}
	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
	{
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;
	}

	if (hi_info_array[id]->dim != 2) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	if (hist) Fill(id,fillX,fillY,weight);
}


// Achim style

void Histo::fill(__int32 id, const char * name, double fillX, double fillY, double weight,
		  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
		  __int32 nYbins, double yLow, double yUp, const char * titleY, const char * dir)

{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH2D * hist = (TH2D *) hiarray[id];

	//--if the histo does not exist create it first--//
	if (!hist)
	{
		hist = new TH2D(name, title, nXbins, xLow, xUp, nYbins, yLow, yUp);
		hist->SetOption("colz");
		hist->SetXTitle(titleX);
//		hist->GetXaxis()->CenterTitle(true);
		hist->SetYTitle(titleY);
//		hist->GetYaxis()->CenterTitle(true);
		hist->GetYaxis()->SetTitleOffset(1.5);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory

		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;
		//--write it into rootfile--//
		if(verbose)
			std::cout <<"create 2D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 2;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
		hi_info_array[id]->nbinsy	= hist->GetNbinsY();
		hi_info_array[id]->low		= hist->GetYaxis()->GetXmin();
		hi_info_array[id]->up		= hist->GetYaxis()->GetXmax();
	}
	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
	{
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;
	}

	if (hi_info_array[id]->dim != 2) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	if (hist) Fill(id,fillX,fillY,weight);
}


void Histo::fill2(__int32 id, std::string name, double fillX, double fillY, double weight,
		  std::string title, __int32 nXbins, double xLow, double xUp, std::string titleX, 
		  __int32 nYbins, double yLow, double yUp, std::string titleY, const char * dir)

{
	fill(id,name.c_str(),fillX,fillY,weight,title.c_str(),nXbins,xLow,xUp,titleX.c_str(),nYbins,yLow,yUp,titleY.c_str(),dir);
}


//___________________________________1D____________________________________________________________________
void Histo::fill(__int32 id, const char * name, double fillX, const char * titleX,
				    __int32 nXbins, double xLow, double xUp, const char * dir,
				    double weight)
{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH1D * hist = (TH1D *) hiarray[id];

	//--if the histo does not exist create it first--//
	if (!hist)
	{
		hist = new TH1D(name, name, nXbins, xLow, xUp);
		hist->SetXTitle(titleX);
		hist->GetXaxis()->CenterTitle(true);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory

		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;
		//--write it into root file--//
		if(verbose)
			std::cout <<"create 1D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 1;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
	}
	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
	{
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;
	}

	if (hi_info_array[id]->dim != 1) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	if (hist) Fill(id,fillX,weight);
}


// Achim style


void Histo::fill(__int32 id, const char * name, double fillX, double weight, const char * title,
			  __int32 nXbins, double xLow, double xUp, const char * titleX, const char * dir)

{
	if (id < 0) {
		std::cout << "error: histogram ID "<<id<<"<0"<<std::endl;
		return;
	}
	if (id >= arsize) {
		std::cout << "error: histogram ID "<<id<<">"<<arsize<<std::endl;
		return;
	}
	if (!rootfile) {
		std::cout << "error: rootfile handle = 0 "<<id<<">"<<arsize<<std::endl;
		return;
	}
	TH1D * hist = (TH1D *) hiarray[id];

	
	//--if the histo does not exist create it first--//
	if (!hist)
	{
		hist = new TH1D(name, title, nXbins, xLow, xUp);
		hist->SetXTitle(titleX);
//		hist->GetXaxis()->CenterTitle(true);
		hist->SetDirectory(getDir(rootfile,dir));				//put it in the wanted directory
		//--now add it to the list--//
		hiarray[id]= (TObject *) hist;
		//--write it into root file--//
		if(verbose)
			std::cout <<"create 1D: "<<dir<<"/"<<hist->GetName()<<std::endl;

		if (!hi_info_array[id]) hi_info_array[id] = new hist_info();
		hi_info_array[id]->dim = 1;
		hi_info_array[id]->buffer	= hist->fArray;
		hi_info_array[id]->pointer_to_hist = (TObject*)hist;
		hi_info_array[id]->nbinsx	= hist->GetNbinsX();
		hi_info_array[id]->left		= hist->GetXaxis()->GetXmin();
		hi_info_array[id]->right	= hist->GetXaxis()->GetXmax();
	}
	//--check also if name is the same one--//
	else if (strcmp(hist->GetName(),name))
	{
		std::cout <<"name doesn't match("<<id<<") is:"<<hist->GetName()<<" should be:"<<name<<std::endl;
	}

	if (hi_info_array[id]->dim != 1) {
		std::cout << "wrong dimension ("<< id <<") " << hist->GetName();
		return;
	}

	//--now fill it--//
	if (hist) Fill(id,fillX,weight);
}




/*
void Histo::fill1(__int32 id, const char * name, double fillX, double weight, const char * title,
			  __int32 nXbins, double xLow, double xUp, const char * titleX, const char * dir)
{
	fill(id, name,fillX,weight,title,nXbins,xLow,xUp,titleX,dir);
}
*/

void Histo::fill1(__int32 id, std::string name, double fillX, double weight, std::string title,
			  __int32 nXbins, double xLow, double xUp, std::string titleX, const char * dir)
{
	fill(id, name.c_str(),fillX,weight,title.c_str(),nXbins,xLow,xUp,titleX.c_str(),dir);
}


void Histo::NTupleD(__int32 id, const char *name, const char * title, const char *varlist, __int32 buffersize, double *data, const char * dir)
{
	MyNTuple = 0;
	if (id <0 || id >= tuple_size) {
		printf("\nError: NTuple-Id is smaller than 0 or greater than %i\n",tuple_size-1);
		return;
	}
	MyNTuple = (TNtupleD *) tuple_array[id];

	//--if the histo does not exist create it first--//
	if (!MyNTuple)
	{
		MyNTuple = new TNtupleD(name,title,varlist,buffersize);
		MyNTuple->SetDirectory(getDir(rootfile,dir)); //put it in the wanted directory

		//--now add it to the list--//
		tuple_array[id]= (TObject *) MyNTuple;
		//--write it into root file--//
		if(verbose)
			std::cout <<"create NTuple: "<<dir<<"/"<<MyNTuple->GetName()<<std::endl;
	}
	//--check also if name is the same one--//
	else if (strcmp(MyNTuple->GetName(),name))
	{
		std::cout <<"name doesn't match("<<id<<") is:"<<MyNTuple->GetName()<<" should be:"<<name<<std::endl;
	}

	//--now fill it--//
	if (MyNTuple) MyNTuple->Fill(data);
}





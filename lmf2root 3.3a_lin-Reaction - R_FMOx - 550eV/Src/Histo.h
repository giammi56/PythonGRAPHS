#pragma once

#ifndef __Histo_H_
	#include "OS_Version.h"
	#define __Histo_H_

	#include <TFile.h>
	//#include "stdlib.h"

	#include <TH1.h>
	#include <TH2.h>

	class TH3D;
	class TObject;
	class TH1;
	class TNtupleD;

	#ifdef _DEBUG
		#include "assert.h"
	#endif

	struct hist_info {
		public:
			__int32 nbinsx;
			__int32 nbinsy;
			__int32 dim;
			double left,right,up,low;
			double * buffer;
			TObject * pointer_to_hist;
	};



	class Histo  
	{
	public:
		Histo(__int32 NbrMaxHistos);
		~Histo();
	
		void Reset();

		bool is_defined(__int32 id) {
			if (id < 0 || id >= arsize) return false;
			return hiarray[id] ? true:false;
		}

	// Lutz fill functions:

		void fill(__int32 id, const char *name, double fillX, double fillY, double fillZ,
						const char *titleX, const char *titleY, const char *titleZ,
						__int32 nXbins, double xLow, double xUp,
						__int32 nYbins, double yLow, double yUp,
						__int32 nZbins, double zLow, double zUp, const char * dir = "",
						double weight=1.);		//fill function for 3d histos

		void fill(__int32 id, const char * name, double fillX, double fillY, 
				  const char * titleX, const char * titleY, 
				  __int32 nXbins, double xLow, double xUp, 
				  __int32 nYbins, double yLow, double yUp, const char * dir = "",
				  double weight=1.);			//fill function for 2d histos

		void fill(__int32 id, const char * name, double fillX, const char * titleX,
				  __int32 nXbins, double xLow, double xUp, const char * dir = "",
				  double weight=1.);			//fill function for 1d histos


	// Achim fill functions:

		void SetBinContentAt(__int32 id, double x, double content);
		void SetBinContentAt(__int32 id, double x, double y, double content);


		double GetBinContentAt(__int32 id, double x);
		double GetBinContentAt(__int32 id, double x, double y);

		void fill(__int32 id, const char * name, double fillX, double fillY, double fillZ, double weight,
				  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
				  __int32 nYbins, double yLow, double yUp, const char * titleY, 
				  __int32 nZbins, double zLow, double zUp, const char * titleZ,const char * dir = "");			//fill function for 3d histos

		void fill3(__int32 id, const char * name, double fillX, double fillY, double fillZ, double weight,
				  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
				  __int32 nYbins, double yLow, double yUp, const char * titleY, 
				  __int32 nZbins, double zLow, double zUp, const char * titleZ,const char * dir = "");			//fill function for 3d histos


		void fill(__int32 id, const char * name, double fillX, double fillY, double weight,
				  const char * title, __int32 nXbins, double xLow, double xUp, const char * titleX, 
				  __int32 nYbins, double yLow, double yUp, const char * titleY, const char * dir = "");			//fill function for 2d histos

		void fill2(__int32 id, std::string name, double fillX, double fillY, double weight,
				  std::string title, __int32 nXbins, double xLow, double xUp, std::string titleX, 
				  __int32 nYbins, double yLow, double yUp, std::string titleY, const char * dir = "");			//fill function for 2d histos


		void fill(__int32 id, const char * name, double fillX, double weight, const char * title,
				  __int32 nXbins, double xLow, double xUp, const char * titleX, const char * dir = "");			//fill function for 1d histos

		void fill1(__int32 id, std::string name, double fillX, double weight, std::string title,
				  __int32 nXbins, double xLow, double xUp, std::string titleX, const char * dir = "");			//fill function for 1d histos
	
		inline void fill1(__int32 id, double fillX, double weight){
			#ifdef _DEBUG
				assert(id >= 0);
				assert(id < arsize);
				assert(hiarray[id]);
			#endif
			Fill(id,fillX,weight);
			return;
		}

		inline void fill2(__int32 id, double fillX, double filly, double weight){
			#ifdef _DEBUG
				assert(id >= 0);
				assert(id < arsize);
				assert(hiarray[id]);
			#endif
			Fill(id,fillX,filly,weight);
			return;
		}
		inline void fill1(__int32 id, double fillX) {
			#ifdef _DEBUG
				assert(id >= 0);
				assert(id < arsize);
				assert(hiarray[id]);
			#endif
			Fill(id,fillX,1.);
			return;
		}

		inline void fill2(__int32 id, double fillX, double filly) {
			#ifdef _DEBUG
				assert(id >= 0);
				assert(id < arsize);
				assert(hiarray[id]);
			#endif
			Fill(id,fillX,filly,1.);
			return;
		}





	inline void Fill(__int32 id, double x, double w){
		#ifdef _DEBUG
			assert(id>=0);
			assert(id < arsize);
			assert(hiarray[id]);
			assert(hi_info_array[id]);
			assert(hi_info_array[id]->pointer_to_hist);
		#endif
		
			double * buf = hi_info_array[id]->buffer;
			double left  = hi_info_array[id]->left;

			TH1D * hist = (TH1D*)hiarray[id];
			hist->SetEntries(hist->GetEntries()+1.);
			if (x < left)  {++buf[0]; return;}

			__int32 nbins = hi_info_array[id]->nbinsx;

			double right = hi_info_array[id]->right;
			if (x >= right) {++buf[nbins+1]; return;}

			__int32 nn = __int32((x-left)*nbins/(right-left)) + 1;

			buf[nn] += w;

			return;
		}






	inline void Fill(__int32 id, double x, double y, double w){
		#ifdef _DEBUG
			assert(id>=0);
			assert(id < arsize);
			assert(hiarray[id]);
			assert(hi_info_array[id]);
		#endif

			double * buf = hi_info_array[id]->buffer;

			TH2D * hist = (TH2D*)hiarray[id];
			hist->SetEntries(hist->GetEntries()+1.);

			__int32 nbinsx = hi_info_array[id]->nbinsx;
			double leftx  = hi_info_array[id]->left;
			double rightx = hi_info_array[id]->right;
			__int32 nnx;
			if (x < leftx) {
				nnx = 0;
			} else {
				if (x >= rightx) {
					nnx = nbinsx+1;
				} else {
					nnx = __int32((x-leftx)*nbinsx/(rightx-leftx)) + 1;
				}
			}

			__int32 nbinsy = hi_info_array[id]->nbinsy;
			double upy  = hi_info_array[id]->up;
			double lowy = hi_info_array[id]->low;
			__int32 nny;
			if (y < lowy) {
				nny = 0;
			} else {
				if (y >= upy) {
					nny = nbinsy+1;
				} else {
					nny = __int32((y-lowy)*nbinsy/(upy-lowy)) + 1;
				}
			}

			buf[nnx + nny*(2+nbinsx)] += w;

			return;
		}




		void NTupleD(__int32 id, const char *name, const char * title, const char *varlist, __int32 buffersize, double *data, const char *  dir = "");
	
		void flushRootFile();
		void openRootFile(const char * RootFileName);
		TH1 * getHist(__int32 id)const				{return ((TH1*)hiarray[id]);}

		bool verbose;
		TFile * rootfile;

		TObject ** hiarray;
		TObject ** tuple_array;

		hist_info ** hi_info_array;

	private:

		TNtupleD * MyNTuple;
		__int32 arsize;
		__int32 tuple_size;
	};

#endif
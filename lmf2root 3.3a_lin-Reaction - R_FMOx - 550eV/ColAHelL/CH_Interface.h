#include <vector>
#include "CH_event.h"
#include "../ColAHelL/src/CH_Hist.h"

#define COLAHELLVERSION 2.1

// some forward declarations
class ColAHelL;
namespace CH {
	#define EL 0
	#define IO 1
	#define PR 2
	#define US 3
};

// Interface to ColAHelL
class CH_Int {

	private:
		ColAHelL *myCH;

	public:
		CH_Int();
		int ReadConfig();
		int ReadConfig(const char* FName);
		int ReadConfig(const char* FName, bool CrLogFile);
		int ReplaceConfig(const char* FName);
		~CH_Int();

/*
Readconfig return values:
>0 ok
-------------------- strange
-15 no momentum pars
-20 no dets
-------------------- no histogramming
-30 no hist
-------------------- no momenta
-40 no reac
-------------------- no presorters
-50 no presorters
-80 no spec 
-------------------- error
-98 no TOF
-99 no XML
*/

	void CalculateTOFs();

	std::vector<CH::CH_event_struct> PresortEvent();
	bool CheckEvent(CH::CH_event_struct * evt, int chan);

	void ProcessEvent(CH::CH_event_struct * evt, int ReactionDefNum=0);
	void ProcessEvent(int ReactionDefNum);
	void FillHistograms();

	double GetBunchSpacing();
	int GetBMChannel();
	bool GetUseBM();

	ColAHelL* CH_Int::GetColAHelL();
	double GetPDetSize();
	double GetIDetSize();
	double GetEDetSize();
	int GetMaxPrsHists();
	int GetMaxHists();
	int GetMaxUserHists();
	CH::histo_handler* GetPrsHists();
	
	void ResetEvent();
	void ResetEvent(__int64 event_number);
	void SetEvent_xyt(int det, double x, double y, double mcp_time, double method=0); 
	void SetEvent_xyttof(int det, double x, double y, double mcp_time, double tof, double method=0); 
	void SetEvent_scanvals(double scanval); 
	void SetEvent_bm(double bm); 
	void SetEvent_reaction(double reaction); 
	void SetBMTime(double* tdc_ns, int* cnt);
	void SetMasterFolder(const char* folder);
	CH::CH_event_struct* GetEvent();
	CH::histo_handler* GetHists(int Type);

	int EvtBelongsToReactions(int channel);
};


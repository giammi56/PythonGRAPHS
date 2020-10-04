#include "ColAHelL.h"
#include "../CH_Interface.h"

CH_Int::CH_Int() 
{
	myCH = new ColAHelL();
	myCH->SetVersion(COLAHELLVERSION);
}

int CH_Int::ReadConfig() 
{
	return myCH->InitAll("ColAHelL.cfg",false);
}
int CH_Int::ReadConfig(const char* FName) 
{
	return myCH->InitAll(FName, false);
}

int CH_Int::ReadConfig(const char* FName, bool CrLogFile) 
{
	return myCH->InitAll(FName, CrLogFile);
}

int CH_Int::ReplaceConfig(const char* FName) 
{
	return myCH->ReplaceCHCfg(FName);
}

CH_Int::~CH_Int() 
{
	delete myCH; myCH=0; 
}

void CH_Int::CalculateTOFs() {
	myCH->CalculateTOFs();
}

bool CH_Int::CheckEvent(CH::CH_event_struct * evt, int chan) {
	return myCH->CheckEvent(evt, chan);
}

std::vector<CH_event_struct> CH_Int::PresortEvent() {
	return myCH->PresortEvent();
}

void CH_Int::ProcessEvent(CH::CH_event_struct * evt, int ReactionDefNum) {
	myCH->ProcessEvent(evt, ReactionDefNum);
}
void CH_Int::ProcessEvent(int ReactionDefNum) {
	myCH->ProcessEvent(ReactionDefNum);
}

void CH_Int::FillHistograms() {
	myCH->FillCHHists();
}

double CH_Int::GetBunchSpacing() {
	return myCH->CTof->GetBMspacing();
}
int CH_Int::GetBMChannel() {
	return myCH->CTof->GetBMchannel();
}
bool CH_Int::GetUseBM() {
	return myCH->CTof->usesBM();
}

double CH_Int::GetPDetSize() {
	return myCH->GetPDetSize()/2.0;
}
double CH_Int::GetIDetSize() {
	return myCH->GetIDetSize()/2.0;
}
double CH_Int::GetEDetSize() {
	return myCH->GetEDetSize()/2.0;
}
int CH_Int::GetMaxPrsHists() {
	return myCH->GetMaxPrsHists();
}
histo_handler *CH_Int::GetPrsHists() {
	return myCH->GetPrsHists();
}

int CH_Int::GetMaxHists() {
	return myCH->GetMaxHists();
}
int CH_Int::GetMaxUserHists() {
	return myCH->GetMaxUserHists();
}
histo_handler *CH_Int::GetHists(int Type) {
	return myCH->GetHists(Type);
}

ColAHelL* CH_Int::GetColAHelL() {
	return myCH;
}

void CH_Int::ResetEvent() {
	myCH->ResetEventStruct(myCH->evt,0);
}
void CH_Int::ResetEvent(__int64 event_number) {
	myCH->ResetEventStruct(myCH->evt, event_number);
}
void CH_Int::SetEvent_xyt(int det, double x, double y, double mcp_time, double method) {
	switch(det) {
		case EL:
			myCH->evt->e.x.push_back(x);
			myCH->evt->e.y.push_back(y);
			myCH->evt->e.time.push_back(mcp_time);
			myCH->evt->e.tof.push_back(0.0);
			myCH->evt->e.method.push_back(method);
			myCH->evt->e.num_hits++;
			break;
		case IO:
			myCH->evt->r.x.push_back(x);
			myCH->evt->r.y.push_back(y);
			myCH->evt->r.time.push_back(mcp_time);
			myCH->evt->r.tof.push_back(0.0);
			myCH->evt->r.method.push_back(method);
			myCH->evt->r.num_hits++;
			break;
		case PR:
			myCH->evt->p.x.push_back(x);
			myCH->evt->p.y.push_back(y);
			myCH->evt->p.time.push_back(mcp_time);
			myCH->evt->p.tof.push_back(0.0);
			myCH->evt->p.method.push_back(method);
			myCH->evt->p.num_hits++;
			break;
	}
}
void CH_Int::SetEvent_xyttof(int det, double x, double y, double mcp_time, double tof, double method) {
	switch(det) {
		case EL:
			myCH->evt->e.x.push_back(x);
			myCH->evt->e.y.push_back(y);
			myCH->evt->e.time.push_back(mcp_time);
			myCH->evt->e.tof.push_back(tof);
			myCH->evt->e.method.push_back(method);
			myCH->evt->e.num_hits++;
			break;
		case IO:
			myCH->evt->r.x.push_back(x);
			myCH->evt->r.y.push_back(y);
			myCH->evt->r.time.push_back(mcp_time);
			myCH->evt->r.tof.push_back(tof);
			myCH->evt->r.method.push_back(method);
			myCH->evt->r.num_hits++;
			break;
		case PR:
			myCH->evt->p.x.push_back(x);
			myCH->evt->p.y.push_back(y);
			myCH->evt->p.time.push_back(mcp_time);
			myCH->evt->p.tof.push_back(tof);
			myCH->evt->p.method.push_back(method);
			myCH->evt->p.num_hits++;
			break;
	}
}
void CH_Int::SetEvent_scanvals(double scanval) {
	myCH->evt->scanval.push_back(scanval);
}
void CH_Int::SetEvent_bm(double bm) {
	myCH->evt->bunchmarker = bm;
}
void CH_Int::SetEvent_reaction(double r) {
	myCH->evt->reaction = r;
}
void CH_Int::SetBMTime(double* tdc_ns, int* cnt) {
	myCH->evt->bunchmarker = myCH->CTof->GetBMtime(tdc_ns, cnt);
}

void CH_Int::SetMasterFolder(const char* folder) {
	myCH->SetMasterFolder(folder);
}

CH_event_struct* CH_Int::GetEvent() {
	return myCH->evt;	
}

// returns how many different reactions were defined with channel number CHANNEL
int CH_Int::EvtBelongsToReactions(int channel) {
	return myCH->EvtBelongsToReactions(channel);
}
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "TStyle.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TProfile.h"
#include "TCutG.h"
#include "TPolyLine.h"
#include "TLine.h"
#include "TMarker.h"
#include "TObjArray.h"
#include "TF1.h"
#include "ipa.h"

void ipa_class::enable_fishtool(int bins, double min, double max) {
	this->NumPars[PPCOPAGE]= 5;
	this->fishtool_enabled = true;
	this->fishtool_min=min;
	this->fishtool_max=max;

	this->H_fishtool_x = root->newTH2D("fish_x","The Fish Filet, x-direction",bins,min,max,"time-of-flight [ns]",100,-80.0,80.0,"x [mm]","COLZ");
	this->H_fishtool_y = root->newTH2D("fish_y","The Fish Filet, y-direction",bins,min,max,"time-of-flight [ns]",100,-80.0,80.0,"y [mm]","COLZ");
	this->H_fishtool_r = root->newTH2D("wiggle","The Wiggle",bins,min,max,"time-of-flight [ns]",100,0.0,80.0,"r [mm]","COLZ");

}

void ipa_class::init_fishtool() {

	this->myCanvas2 = root->newCanvas("IPA_PPCO","IPA - the fish tool (press space to update)",40,20,330*2,330*3);	
	this->myCanvas2->SetFixedAspectRatio(false);
	this->myCanvas2->Divide(1,3);

	for(int j=0;j<3;j++) { 
		for(int i=0;i<36;i++) {
			mrk[0][j][i]= new TMarker();
			mrk[1][j][i]= new TMarker();
			mrk[2][j][i]= new TMarker();
		}
	}

	for(int i=0;i<128;i++) {
		knots[i] = new bpoint();
		firstControlPoints[i] = new bpoint();
		secondControlPoints[i] = new bpoint();
	}

	this->fishtool_t_0_elec = IPA_Proc[0]->e_fac->dt;
	this->fishtool_ee1 = 1.0;
	this->fishtool_ee2 = 5.0;
	this->fishtool_ee3 = 10.0;
	add_par(this->fishtool_ee1,0.0,1000.0,1.0,"Energy e1",70);
	add_par(this->fishtool_ee2,0.0,1000.0,1.0,"Energy e2",71);
	add_par(this->fishtool_ee3,0.0,1000.0,1.0,"Energy e3",72);
	add_par(0.0,0.0,10.0,0.1,"Yaw factor",73);
	add_par(0.0,-360.0,360.0,5.0,"Yaw angle",74);

}

void ipa_class::Fill_fishtool(CH_event_struct *evnt) 
{
	if(evnt->e.num_hits>0) {

		e_raw[0]->method = evnt->e.method[0];
		e_raw[0]->data.x = evnt->e.x[0];
		e_raw[0]->data.y = evnt->e.y[0];
		e_raw[0]->data.time = evnt->e.time[0];
		e_raw[0]->data.tof = evnt->e.tof[0];
	
		e[0][0]->set_raw(*e_raw[0]);
		e[0][0]->shift_stretch_raw(IPA_Proc[0]->e_fac);
		e[0][0]->rotate_raw(IPA_Proc[0]->e_fac->rot_ang);
		e[0][0]->calc_phi_pos();		

		//if(fabs(e[0][0]->raw.phi-90.0)<15.0 || fabs(e[0][0]->raw.phi+90.0)<15.0) 		
		if(fabs(e[0][0]->raw.data.y)<5.0)
			this->H_fishtool_x->Fill(e[0][0]->raw.data.tof,e[0][0]->raw.data.x);

	//	if(fabs(e[0][0]->raw.phi)<15.0 || fabs(e[0][0]->raw.phi-180.0)<15.0) 
		if(fabs(e[0][0]->raw.data.x)<5.0)
			this->H_fishtool_y->Fill(e[0][0]->raw.data.tof,e[0][0]->raw.data.y);
		this->H_fishtool_r->Fill(e[0][0]->raw.data.tof,sqrt(e[0][0]->raw.data.x*e[0][0]->raw.data.x+e[0][0]->raw.data.y*e[0][0]->raw.data.y));
		num_events_read++;
	}
}

void ipa_class::draw_fishtool() {
	double U[4];
	U[0]=0.0;

	for(int j=0;j<IPA_Proc[0]->Spect->electron_side->number_of_regions;j++) {
		U[j+1]=U[j]+IPA_Proc[0]->Spect->electron_side->Efields[j]*IPA_Proc[0]->Spect->ion_side->lengths[j]/10.0;
	}

	this->myCanvas2->cd(1);
	this->H_fishtool_x->Draw("COLZ");
	draw_fishtool_markers(0,U);

	this->myCanvas2->cd(2);
	this->H_fishtool_y->Draw("COLZ");
	draw_fishtool_markers(1,U);

	this->myCanvas2->cd(3);
	this->H_fishtool_r->Draw("COLZ");

	draw_fishtool_markers(2,U);

	this->myCanvas2->Modified();
	this->myCanvas2->Update();
}

void ipa_class::draw_fishtool_markers(int type, double *U) {
	double EM = 1.0;
	double Omega = 2.618e8/15. * IPA_Proc[0]->Spect->Bfield_G;
	double pau = 5.0185113E23;

	double yf[3];
	yf[0] = par[73].value * cos(par[74].value/180.0*PI);
	yf[1] = par[73].value * sin(par[74].value/180.0*PI);
	yf[2] = 0.0;

	for(int e_num=0;e_num<3;e_num++) {
		for(int i = 0 ; i < (type==2?18:36); i++) {

			double te = ((double)(i*10) + 0.5);

			double pe = sqrt(2.*par[70+e_num].value/27.6/EM);				// total momentum
			double pz = pe * cos(te/180.*3.14152);			// momentum in field
			double py = pe * sin(te/180.*3.14152);	// momentum in pos direction	

			double mt = TTOT(pz, U);
			double my = (py*2.*sin(Omega * mt*1.e-9 / 2.)/(Omega*MEKG*pau))*1000.; //mm

			// add yaw
			mt += my*yf[type]/10.0;

			mrk[type][e_num][i]->SetMarkerStyle(7);
			mrk[type][e_num][i]->SetX(mt);
			mrk[type][e_num][i]->SetY(my);
			mrk[type][e_num][i]->SetMarkerColor(e_num%2 + 1);
			mrk[type][e_num][i]->SetMarkerSize(1.5f);
			mrk[type][e_num][i]->Draw();
		}
	}

	// Bezier curve
	for(int i=0; i<36; i++) {
		knots[i]->x = (float)mrk[0][2][i]->GetX();
		knots[i]->y = (float)mrk[0][2][i]->GetY();
	}
	CalcControlPoints(36);

	for(int i=0;i<(type==2?17:35);i++) {
		bpoint p1 = BezierCube(knots[i], firstControlPoints[i], secondControlPoints[i], knots[i+1], 0.33f);
		bpoint p2 = BezierCube(knots[i], firstControlPoints[i], secondControlPoints[i], knots[i+1], 0.66f);
		TLine *l= new TLine(knots[i]->x,knots[i]->y,p1.x,p1.y);
		l->Draw();
		l= new TLine(p1.x,p1.y,p2.x,p2.y);
		l->Draw();
		l= new TLine(p2.x,p2.y,knots[i+1]->x,knots[i+1]->y);
		l->Draw();
	}
}

double ipa_class::TTOT(double pz, double *U) 
{                                    
	
	double TGES = 0.;
	double Q = 1.0;
	double EM = 1.0;
	double V0 = 0.;
	double V1 = 0.;
	double vs = 0.;

	double EStart =pz*pz/2.*EM*27.6;
	if(pz<0.0) {
		double dzturn = EStart / IPA_Proc[0]->Spect->electron_side->Efields[0];
		vs = 2.188e-1*sqrt(EStart/EM/27.21*2.)+0.0000000001;
		TGES = 2. * (dzturn*2./vs);  // 2 x turnaround
	}
	for(int i = 0; i < IPA_Proc[0]->Spect->electron_side->number_of_regions ; i++) {  

		V0=2.188e-1*sqrt((Q*U[i]+EStart)/EM/27.21*2.); 
		V1=2.188e-1*sqrt((Q*U[i+1]+EStart)/EM/27.21*2.); 

		TGES = TGES+(IPA_Proc[0]->Spect->electron_side->lengths[i]/10.0)*2./(V0+V1);
	}	
    
	return TGES; 
}

// Calculate Bezier-Curve
// --------------------------------------------------------------------------------
float ipa_class::getPt(float n1, float n2 , float perc)
{
    float diff = n2 - n1;
    return n1 + ( diff * perc );
}  

// P0: start
// P2: end
// P1: Control point
bpoint ipa_class::BezierQuad(bpoint *P0, bpoint *P1, bpoint *P2, float perc)  
{
    bpoint out;
	
    float xa = getPt( P0->x , P1->x , perc );
    float ya = getPt( P0->y , P1->y , perc );
    float xb = getPt( P1->x , P2->x , perc );
    float yb = getPt( P1->y , P2->y , perc );

    out.x = getPt( xa , xb , perc );
    out.y = getPt( ya , yb , perc );

    return out;
}

// P0: start
// P3: end
// P1: Control point 1 (close to P0)
// P2: Control point 2 (close to P3)
bpoint ipa_class::BezierCube(bpoint *P0, bpoint *P1, bpoint *P2, bpoint *P3, float perc)  
{
    bpoint out;
	
    float xa = getPt( P0->x , P1->x , perc );
    float ya = getPt( P0->y , P1->y , perc );
    float xb = getPt( P1->x , P2->x , perc );
    float yb = getPt( P1->y , P2->y , perc );
    float xc = getPt( P2->x , P3->x , perc );
    float yc = getPt( P2->y , P3->y , perc );

    float xm = getPt( xa , xb , perc );
    float ym = getPt( ya , yb , perc );
    float xn = getPt( xb , xc , perc );
    float yn = getPt( yb , yc , perc );

    out.x = getPt( xm , xn , perc );
    out.y = getPt( ym , yn , perc );

    return out;
}

void ipa_class::CalcControlPoints(int num_knots)
{
	double x[128];
	double y[128];
	double rhs[128];

	int n = num_knots - 1;

	// Set right hand side X values
	for (int i = 1; i < n - 1; ++i)
		rhs[i] = 4 * knots[i]->x + 2 * knots[i + 1]->x;
	rhs[0] = knots[0]->x + 2 * knots[1]->x;
	rhs[n - 1] = (8 * knots[n - 1]->x + knots[n]->x) / 2.0;
	// Get first control points X-values
	GetFirstControlPoints(rhs, n, x);

	// Set right hand side Y values
	for (int i = 1; i < n - 1; ++i)
		rhs[i] = 4 * knots[i]->y + 2 * knots[i + 1]->y;
	rhs[0] = knots[0]->y + 2 * knots[1]->y;
	rhs[n - 1] = (8 * knots[n - 1]->y + knots[n]->y) / 2.0;
	// Get first control points Y-values
	GetFirstControlPoints(rhs, n, y);

	for (int i = 0; i < n; ++i)
	{
		// First control point
		firstControlPoints[i]->x = float(x[i]);
		firstControlPoints[i]->y = float(y[i]);

		// Second control point
		if (i < n - 1) {
			secondControlPoints[i]->y =  2.0f * knots[i + 1]->y - float(y[i + 1]);
			secondControlPoints[i]->x =  2.0f * knots[i + 1]->x - float(x[i + 1]);
		}
		else {
			secondControlPoints[i]->x = (knots[n]->x + float(x[n - 1])) / 2.0f; 
			secondControlPoints[i]->y = (knots[n]->x + float(y[n - 1])) / 2.0f;
		}
	}
}

 void ipa_class::GetFirstControlPoints(double *rhs, int n, double res[])
{
	double tmp[128];

	double b = 2.0;
	res[0] = rhs[0] / b;
	for (int i = 1; i < n; i++) // Decomposition and forward substitution.
	{
		tmp[i] = 1 / b;
		b = (i < n - 1 ? 4.0 : 3.5) - tmp[i];
		res[i] = (rhs[i] - res[i - 1]) / b;
	}
	for (int i = 1; i < n; i++)
		res[n - i - 1] -= tmp[n - i] * res[n - i]; // Backsubstitution.
}
// --------------------------------------------------------------------------------
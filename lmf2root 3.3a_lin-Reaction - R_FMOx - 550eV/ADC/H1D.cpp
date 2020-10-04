#include "H1D.h"
#pragma warning(disable : 4996)

AXIS::AXIS(int nbins, double low, double up, const char* name) {
	m_nbr_bins = nbins;
	m_low = low;
	m_up = up;
	m_name = name;
	m_range = up - low;
	m_bin_size = m_range / nbins;
	m_bin_size_half = m_bin_size / 2.;
}

int AXIS::FindBin(double x) {
	int bin = (int)((x - m_low) / m_range * m_nbr_bins + 1);
	if (bin < 1)// noch grenzen drinne machen..
		return 0;
	if (bin > m_nbr_bins)
		return m_nbr_bins + 1;
	return bin;
}

double AXIS::GetXlow() {
	return m_low;
}

double AXIS::GetXup() {
	return m_up;
}

double AXIS::GetXRange() {
	return m_range;
}

double AXIS::GetBinSize() {
	return m_bin_size;
}

double AXIS::GetBinSizeHalf() {
	return m_bin_size_half;
}

int AXIS::GetNbins() {
	return m_nbr_bins;
}

double AXIS::GetBinCenter(int bin) {
	if (bin > 0 && bin <= m_nbr_bins) {
		return m_low + (bin-1) * m_bin_size + m_bin_size_half;
	}
	return 0;
}

H1D::H1D() {
	m_name			= nullptr;
	m_title			= nullptr;
}

H1D::H1D(const char *name,const char *title, int nbinsx, double xlow, double xup) {

	m_name = new char[strlen(name)+1];
	sprintf(m_name,"%s",name);

	m_title = new char[strlen(title)+1];
	sprintf(m_title, "%s", title);

	if (nbinsx > 0) {
		m_data = new double[nbinsx+2];
		for (int i = 0; i < nbinsx + 2; i++)
			m_data[i] = 0;
	}

	xaxis = new AXIS(nbinsx, xlow, xup);
}

H1D::~H1D() {
	delete[] m_data;
	return;
}

void H1D::Fill(double x) {
	m_data[this->GetXaxis()->FindBin(x)] += 1;
	return;
}

void H1D::Fill(double x, double w) {
	m_data[this->GetXaxis()->FindBin(x)] += w;
	return;
}

bool H1D::SetBinContent(int bin, double val) {
	if (bin > 0 && bin <= this->GetXaxis()->GetNbins()) {
		m_data[bin] = val;
		return true;
	}
	else {
		return false;
	}

}

double H1D::GetBinContent(int bin) {
	if (bin >= 0 && bin <= this->GetXaxis()->GetNbins() + 1)
		return m_data[bin];
	return 0;
}

double H1D::GetBinCenter(int bin) {
	return xaxis->GetBinCenter(bin);
}

bool H1D::PlotData() {
	std::cout << "\n----------------------------------------";
	std::cout << "\n Name  :  " << m_name;
	std::cout << "\n Title :  " << m_title;
	std::cout << "\n From  :" << this->GetXaxis()->GetXlow() << " to " << this->GetXaxis()->GetXup() << " (" << this->GetXaxis()->GetNbins() << " bins)";
	std::cout << "\n----------------------------------------";
	std::cout << "\n#Bin\t\tx-Val\t\tData";
	for (int bin = 1; bin <= this->GetXaxis()->GetNbins(); bin++)
		std::cout << "\n" << bin << "\t\t" << this->GetXaxis()->GetBinCenter(bin) << "\t\t" << m_data[bin];
	std::cout << "\nUnderflow : " << m_data[0];
	std::cout << "\nOverflow  : " << m_data[this->GetXaxis()->GetNbins() + 1];
	return true;
}

int H1D::GetNbins() {
	return this->GetXaxis()->GetNbins();
}

double H1D::GetXmin() {
	return this->GetXaxis()->GetXlow();
}

double H1D::GetXmax() {
	return this->GetXaxis()->GetXup();
}

const char*	H1D::GetName() {
	return m_name;
}

void H1D::SetName(const char* name) {
	if (m_name)
		delete m_name;
	m_name = 0;

	m_name = new char[strlen(name) + 1];
	sprintf(m_name, "%s", name);
}

const char* H1D::GetTitle() {
	return m_title;
}
void H1D::SetTitle(const char* title){
	if (m_title)
		delete m_title;
	m_title = 0;

	m_title = new char[strlen(title) + 1];
	sprintf(m_title, "%s", title);
}

AXIS* H1D::GetXaxis() {
	return xaxis;
}

int H1D::GetBinMaximum() {
	int bin_max = 1;
	double max = this->m_data[1];
	for (int i = 2; i <= this->GetXaxis()->GetNbins(); i++) {
		if (this->m_data[i] > max) {
			bin_max = i;
			max		= this->m_data[i];
		}
	}
	return bin_max;
}

double H1D::GetMaximum(){
	return m_data[this->GetBinMaximum()];
}

void H1D::NormalizeTo(double value) {
	double norm = value / this->GetMaximum();
	for (int i = 1; i <= this->GetXaxis()->GetNbins(); i++) {
		m_data[i] = m_data[i] * norm;
	}
}

double H1D::GetIntegral() {
	return GetIntegral(1, this->GetXaxis()->GetNbins());
}

double H1D::GetIntegral(int bin_low, int bin_up) {
	if (bin_low > 0 && bin_up <= this->GetXaxis()->GetNbins()) {
		double integral = 0; 
		for (int i = bin_low; i <= bin_up; i++) {
			integral += this->GetBinContent(i);
		}
		return integral;
	}
	return 0.;
}

double H1D::InterpolateLinear(double x) {
	// out of range
	double a = this->GetXaxis()->GetXlow();
	double b = this->GetXaxis()->GetXup();

	if (x < this->GetXaxis()->GetXlow() || x > this->GetXaxis()->GetXup())
		return 0;

	int xbin = this->GetXaxis()->FindBin(x);

	if (x < this->GetXaxis()->GetBinCenter(xbin) || xbin == this->GetXaxis()->GetNbins()) {
	// take left neighbor
		return this->GetBinContent(xbin-1) + (this->GetBinContent(xbin) - this->GetBinContent(xbin-1)) * (x - this->GetXaxis()->GetBinCenter(xbin-1)) / this->GetXaxis()->GetBinSize();
	}
	// take right neighbor
	return this->GetBinContent(xbin) + (this->GetBinContent(xbin + 1) - this->GetBinContent(xbin)) * (x - this->GetXaxis()->GetBinCenter(xbin)) / this->GetXaxis()->GetBinSize();
}


double H1D::parabola_fit(double y1, double y2, double y3, double xpos)
{
	//    y1---y2---y3
	// x:  0--- 1--- 2 (xpos)
	double a = 0.5*(y1 + y3) - y2;
	double b = 2 * y2 - 1.5*y1 - 0.5*y3;
	return (a*xpos + b)*xpos + y1;
}

double H1D::InterpolateParabola(double x) {
	if (x < this->GetXmin() || x > this->GetXmax())
		return 0;

	int xbin = this->GetXaxis()->FindBin(x);
	
	if (xbin < 3) {
		double dx = (x - this->GetXaxis()->GetBinCenter(2)) / this->GetXaxis()->GetBinSize()+1;
		return this->parabola_fit(this->GetBinContent(1), this->GetBinContent(2), this->GetBinContent(3), dx);
	}

	if (xbin > this->GetXaxis()->GetNbins() - 2) {
		double dx = (x - this->GetXaxis()->GetBinCenter(this->GetNbins()-1)) / this->GetXaxis()->GetBinSize() + 1;
		return this->parabola_fit(this->GetBinContent(this->GetNbins()-2), this->GetBinContent(this->GetNbins()-1), this->GetBinContent(this->GetNbins()), dx);
	}

	if (x < this->GetXaxis()->GetBinCenter(xbin)) {
		double dx = (x - this->GetXaxis()->GetBinCenter(xbin)) / this->GetXaxis()->GetBinSize();
		
		double y1 = this->parabola_fit(this->GetBinContent(xbin - 1), this->GetBinContent(xbin), this->GetBinContent(xbin + 1), 1 + dx);
		double y2 = this->parabola_fit(this->GetBinContent(xbin - 2), this->GetBinContent(xbin - 1), this->GetBinContent(xbin), 2 + dx);

		return (1 + dx)*y1 - dx*y2;
	}
	else {
		double dx = (x - this->GetXaxis()->GetBinCenter(xbin)) / this->GetXaxis()->GetBinSize();

		double y1 = this->parabola_fit(this->GetBinContent(xbin - 1), this->GetBinContent(xbin), this->GetBinContent(xbin + 1), 1 + dx);
		double y2 = this->parabola_fit(this->GetBinContent(xbin), this->GetBinContent(xbin + 1), this->GetBinContent(xbin + 2), dx);

		return (1 - dx)*y1 + dx*y2;
	}
	return 0;
}




/* This functin will find the borders symmetrical around the center of mass to hold the given percentage of mass (plus a safty percentage) */
void H1D::GetRangeHisto(double *range_low, double *range_up, double percentage, double safty_percentage) {
	int xbins = this->GetXaxis()->GetNbins();
	double integral = this->GetIntegral();

	int sum_center = 0;
	int bin_center = 0;
	for (int i = 1; sum_center <= integral / 2; i++) {
		sum_center += (int)this->GetBinContent(i);
		bin_center = i;
	}

	double per = percentage / 100 / 2;

	// go left 
	double sum_temp = this->GetBinContent(bin_center) / 2;
	for (int i = bin_center - 1; i > 0; i--) {
		sum_temp += this->GetBinContent(i);
		if (sum_temp / integral >= per) {
			*range_low = this->GetXaxis()->GetBinCenter(i);
			break;
		}
	}

		//go right
	sum_temp = this->GetBinContent(bin_center) / 2;
	for (int i = bin_center + 1; i < xbins; i++) {
		sum_temp += this->GetBinContent(i);
		if (sum_temp / integral >= per) {
			*range_up = this->GetXaxis()->GetBinCenter(i);
			break;
		}
	}

	// give safty 
	*range_low	= *range_low - (safty_percentage / 100) * abs(*range_up - *range_low);
	*range_up	= *range_up  + (safty_percentage / 100) * abs(*range_up - *range_low);
}

void H1D::Smooth(int ntimes) {
	for (int n = 0; n < ntimes; n++) {
		// go right
		for (int i = 2; i < this->GetXaxis()->GetNbins(); i++) {
			this->SetBinContent(i, 0.25*this->GetBinContent(i - 1) + 0.5*this->GetBinContent(i) + 0.25*this->GetBinContent(i + 1));
		}
		//go left
		for (int i = this->GetXaxis()->GetNbins()-1; i > 1; i--) {
			this->SetBinContent(i, 0.25*this->GetBinContent(i - 1) + 0.5*this->GetBinContent(i) + 0.25*this->GetBinContent(i + 1));
		}
	}
}
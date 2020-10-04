#include "H1D.h"
#pragma warning(disable : 4996)

void H2D::Fill(double x, double y) {
	m_data[this->GetXaxis()->FindBin(x)][this->GetYaxis()->FindBin(y)] += 1;
	return;
}

void H2D::Fill(double x, double y, double w) {
	m_data[this->GetXaxis()->FindBin(x)][this->GetYaxis()->FindBin(y)] += w;
	return;
}

H2D::H2D(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup) {
	m_name = new char[strlen(name) + 1];
	sprintf(m_name, "%s", name);

	m_title = new char[strlen(title) + 1];
	sprintf(m_title, "%s", title);
	
	xaxis = new AXIS(nbinsx, xlow, xup);
	yaxis = new AXIS(nbinsy, ylow, yup);

	m_data = new double*[xaxis->GetNbins()+2];
	for (int x = 0; x < xaxis->GetNbins()+2; x++) {
		m_data[x] = new double[yaxis->GetNbins()+2];
		for (int y = 0; y < yaxis->GetNbins()+2; y++)
			m_data[x][y] = 0;
	}
}

AXIS* H2D::GetYaxis() {
	return yaxis;
}

bool H2D::SetBinContent(int xbin, int ybin, double val) {
	if (xbin > 0 && xbin <= this->GetXaxis()->GetNbins() && ybin > 0 && ybin <= this->GetYaxis()->GetNbins()) {
		m_data[xbin][ybin] = val;
		return true;
	}
	else {
		return false;
	}
}

double H2D::GetBinContent(int xbin, int ybin) {
	if (xbin >= 0 && xbin <= this->GetXaxis()->GetNbins()+1 && ybin >= 0 && ybin <= this->GetYaxis()->GetNbins()+1) {
		return m_data[xbin][ybin];
	}
	else {
		return 0.;
	}
}

bool H2D::GetDimensions(int *xbins, int *ybins) {
	*xbins = xaxis->GetNbins();
	*ybins = yaxis->GetNbins();
	return true;
}

int H2D::GetBinMaximum(int &x_bin_max, int &y_bin_max) {
	int x_bin_max_out = 1;
	int y_bin_max_out = 1;
	double max = this->m_data[1][1];
	for (int i = 1; i <= this->GetXaxis()->GetNbins(); i++) {
		for (int j = 1; j <= this->GetYaxis()->GetNbins(); j++) {
			if (this->m_data[i][j] > max) {
				x_bin_max_out = i;
				y_bin_max_out = j;
				max = this->m_data[i][j];
			}
		}
	}

	x_bin_max = x_bin_max_out; 
	y_bin_max = y_bin_max_out;

	return 1;
}

double H2D::GetMaximum() {
	int xmax, ymax;
	this->GetBinMaximum(xmax, ymax);
	return this->m_data[xmax][ymax];
}

void H2D::NormalizeTo(double value) {
	double norm = value / this->GetMaximum();
	for (int x = 1; x <= this->GetXaxis()->GetNbins(); x++) {
		for (int y = 1; y <= this->GetYaxis()->GetNbins(); y++) {
			m_data[x][y] = m_data[x][y] * norm;
		}
	}

}


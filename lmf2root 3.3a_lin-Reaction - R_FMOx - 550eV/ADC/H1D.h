#ifndef __H1D_H__
#define __H1D_H__

#include <iostream>
#include <stdio.h>

class AXIS {
public:
	AXIS(int nbins, double low, double up, const char* name = "");

	int		FindBin(double x);
	int		GetNbins();
	double	GetXlow();
	double	GetXup();
	double	GetXRange();
	double	GetBinCenter(int bin);
	double	GetBinSize();
	double	GetBinSizeHalf();

private:
	const char* m_name;
	int m_nbr_bins;
	double m_low, m_up;
	double m_range, m_bin_size, m_bin_size_half;
};

class H1D {
	friend AXIS;
public: 
	H1D();
	H1D(const char* name,const char* title, int nbinsx, double xlow, double xup);
	~H1D();

	void		Fill(double x);
	void		Fill(double x, double w);

	AXIS*		GetXaxis();
	bool		SetBinContent(int bin, double val);
	double		GetBinContent(int bin);
	double		GetBinCenter(int bin);

	int			GetNbins();
	double		GetXmin();
	double		GetXmax();

	const char*	GetName();
	void		SetName(const char* name);
	const char* GetTitle();
	void		SetTitle(const char* title);

	int			GetBinMaximum();
	double		GetMaximum();

	double		GetIntegral();
	double		GetIntegral(int bin_low, int bin_up);

	void		GetRangeHisto(double *range_low, double *range_up, double percentage, double safty_percentage = 0.);

	void		NormalizeTo(double value = 1);
	void		Smooth(int ntimes);

	bool		PlotData();
	
	double		InterpolateLinear(double x);
	double		InterpolateParabola(double x);

	
protected:
	char		*m_name;
	char		*m_title;
	
	double		*m_data;
	AXIS		*xaxis;	

private: 
	double		parabola_fit(double y1, double y2, double y3, double xpos);
};

class H2D : public H1D {
public:
	H2D(const char* name, const char* title, int nbinsx, double xlow, double xup, int nbinsy, double ylow, double yup);
	AXIS*		GetYaxis();

	void		Fill(double x, double y);
	void		Fill(double x, double y, double w);

	bool		SetBinContent(int xbin, int ybin, double val);
	double		GetBinContent(int xbin, int ybin);
	bool		GetDimensions(int *xbins, int *ybins);

	int			GetBinMaximum(int &x_bin_max, int &y_bin_max);
	void		NormalizeTo(double value = 1);
	double		GetMaximum();

private:
	double		**m_data;
	AXIS		*yaxis;
};

#endif

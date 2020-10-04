#include "CH_particle.h"
#include "CH_Spectrometer.h"
#include "CH_Basics.h"
#include "CH_FUN_Lowlevel.h"

#include "Math.h"

namespace CH
{

	particle_class::particle_class()
	{
		this->raw.m = 0.0;
		this->raw.q = 0.0;
		this->channel = 0;
		this->valid = false;

		this->mom.x = -1000.0;
		this->mom.y = -1000.0;
		this->mom.z = -1000.0;

		this->cor.dt = 0.0;
		this->cor.dx = 0.0;
		this->cor.dy = 0.0;
		this->cor.overall_stretch = 1.0;
		this->cor.x_stretch = 1.0;
		this->cor.y_stretch = 1.0;

		this->t_mean = 0.0;
		this->t_width = -1.0;
	}

	particle_class::~particle_class()
	{
	}

	void particle_class::set_raw(particle_raw *rawdata) {
		this->raw.data = rawdata->data;
		this->raw.m = rawdata->m;
		this->raw.q = rawdata->q;
		this->raw.method = rawdata->method;
		this->raw.phi = rawdata->phi;
	}

	void particle_class::set_raw(double x, double y, double t, int method) {
		this->raw.data.x = x;
		this->raw.data.y = y;
		this->raw.data.tof = t;

		this->raw.phi = 0;
		this->raw.method = method;
	}

	void particle_class::calc_phi_pos() {
		this->raw.phi = atan2(this->raw.data.x, this->raw.data.y)/PI*180.0;
	}

	void particle_class::set_channel(int channel) {
		this->channel = channel;
	}

	void particle_class::set_channel(int channel, bool valid) {
		this->channel = channel;
		this->valid = valid;
	}
	void particle_class::set_valid(bool valid) {
		this->valid = valid;
	}

	void particle_class::set_t_mean(double t) {
		this->t_mean = t;
	}

	void particle_class::set_t_width(double t) {
		this->t_width = t;
	}

	void particle_class::shift_stretch_raw(cor_param * fac) {
		this->raw.data.x += fac->dx;
		this->raw.data.x *= fac->x_stretch;
		this->raw.data.x *= fac->overall_stretch;

		this->raw.data.y += fac->dy;
		this->raw.data.y *= fac->y_stretch;
		this->raw.data.y *= fac->overall_stretch;

		this->raw.data.tof += fac->dt;

		if(fac->mir_x)
			this->raw.data.x *= -1.0;
		if(fac->mir_y)
			this->raw.data.y *= -1.0;
	}

	void particle_class::shift_stretch_mom(mom_cor_param * fac) {
		this->mom.x += fac->dx;
		this->mom.y += fac->dy;
		this->mom.z += fac->dz;

		this->mom.x *= fac->x_stretch * fac->overall_stretch;
		this->mom.y *= fac->y_stretch * fac->overall_stretch;
		this->mom.z *= fac->z_stretch * fac->overall_stretch;

		if(fac->mir_x)
			this->mom.x *= -1.0;
		if(fac->mir_y)
			this->mom.y *= -1.0;
		if(fac->mir_z)
			this->mom.z *= -1.0;
	}

	void particle_class::rotate_raw(double ang) {
		ang = ang/180.*PI;

		double x_temp = cos(ang) * this->raw.data.x - sin(ang) * this->raw.data.y;
		double y_temp = sin(ang) * this->raw.data.x + cos(ang) * this->raw.data.y;
		this->raw.data.x = x_temp;
		this->raw.data.y = y_temp;
	}

	void particle_class::autotune() {

		int bins_phi = 36;
		int bins_ctheta = 24;
		double corr_tab[36][24] = { { 0.993299 ,0.990797 ,0.978642 ,0.97613 ,0.973376 ,0.972084 ,0.978078 ,0.977131 ,0.979863 ,0.982942 ,0.981677 ,0.984789 ,0.986699 ,0.980286 ,0.977187 ,0.97533 ,0.9713 ,0.970842 ,0.970154 ,0.976847 ,0.98733 ,1.01012 ,1.0285 ,1.05724 },
		{ 0.997405 ,0.992249 ,0.98551 ,0.971333 ,0.969571 ,0.969663 ,0.973168 ,0.97353 ,0.974551 ,0.978601 ,0.977588 ,0.977526 ,0.977846 ,0.975413 ,0.97378 ,0.967058 ,0.963132 ,0.957198 ,0.956288 ,0.962171 ,0.973124 ,1.00048 ,1.02669 ,1.05863 },
		{ 0.994815 ,0.987579 ,0.982722 ,0.974433 ,0.973731 ,0.972871 ,0.969789 ,0.975821 ,0.976639 ,0.977371 ,0.979257 ,0.978924 ,0.977928 ,0.974584 ,0.968946 ,0.964776 ,0.958559 ,0.951005 ,0.951459 ,0.96026 ,0.966964 ,0.994927 ,1.02619 ,1.06205 },
		{ 0.98964 ,0.983762 ,0.981169 ,0.976265 ,0.975698 ,0.978073 ,0.982245 ,0.979292 ,0.976992 ,0.978345 ,0.981765 ,0.977914 ,0.975118 ,0.973282 ,0.964747 ,0.962422 ,0.953554 ,0.943159 ,0.939487 ,0.942747 ,0.963149 ,0.986612 ,1.01652 ,1.04285 },
		{ 0.991294 ,0.981341 ,0.976669 ,0.976155 ,0.971811 ,0.974996 ,0.97181 ,0.976238 ,0.978257 ,0.982039 ,0.978855 ,0.979994 ,0.975242 ,0.973087 ,0.966553 ,0.96328 ,0.955646 ,0.947295 ,0.939286 ,0.94544 ,0.953126 ,0.979654 ,1.01391 ,1.04035 },
		{ 0.988699 ,0.985171 ,0.979434 ,0.969437 ,0.971025 ,0.969241 ,0.970899 ,0.974446 ,0.981349 ,0.977808 ,0.978405 ,0.977154 ,0.976288 ,0.975389 ,0.968909 ,0.958577 ,0.955534 ,0.943439 ,0.936975 ,0.940121 ,0.952686 ,0.974882 ,1.00724 ,1.0498 },
		{ 0.992599 ,0.98278 ,0.985314 ,0.977184 ,0.983489 ,0.980932 ,0.976804 ,0.985058 ,0.981803 ,0.982282 ,0.976519 ,0.976908 ,0.977368 ,0.970119 ,0.96296 ,0.958135 ,0.946804 ,0.936517 ,0.929131 ,0.925508 ,0.939013 ,0.969422 ,1.00902 ,1.0349 },
		{ 0.980682 ,0.989454 ,0.983906 ,0.979144 ,0.98323 ,0.982808 ,0.982349 ,0.986599 ,0.983608 ,0.985297 ,0.981735 ,0.977149 ,0.974539 ,0.96842 ,0.963742 ,0.951732 ,0.942246 ,0.9293 ,0.921389 ,0.921584 ,0.944367 ,0.979528 ,1.01148 ,1.05859 },
		{ 0.962586 ,0.97375 ,0.978824 ,0.976039 ,0.973616 ,0.977466 ,0.969095 ,0.982043 ,0.980939 ,0.984774 ,0.98291 ,0.97696 ,0.97441 ,0.975013 ,0.96606 ,0.959589 ,0.948258 ,0.933358 ,0.930426 ,0.928371 ,0.945628 ,0.977431 ,1.01141 ,1.06851 },
		{ 0.970167 ,0.969782 ,0.965777 ,0.968444 ,0.968585 ,0.962822 ,0.971121 ,0.973055 ,0.973305 ,0.978359 ,0.980142 ,0.978424 ,0.975008 ,0.969588 ,0.968106 ,0.963918 ,0.956193 ,0.946939 ,0.941395 ,0.940769 ,0.96174 ,0.997452 ,1.04767 ,1.08973 },
		{ 0.969083 ,0.966315 ,0.967272 ,0.964508 ,0.9712 ,0.96553 ,0.970662 ,0.970993 ,0.976583 ,0.978717 ,0.98625 ,0.980948 ,0.97988 ,0.976692 ,0.971788 ,0.964455 ,0.957278 ,0.951076 ,0.944978 ,0.949904 ,0.968954 ,1.00693 ,1.05198 ,1.07673 },
		{ 0.982804 ,0.978357 ,0.970216 ,0.96157 ,0.962076 ,0.967893 ,0.975259 ,0.97523 ,0.973137 ,0.977608 ,0.982427 ,0.985277 ,0.984619 ,0.988002 ,0.975676 ,0.972035 ,0.960627 ,0.960799 ,0.953454 ,0.955712 ,0.968217 ,0.989266 ,1.03618 ,1.07899 },
		{ 0.984329 ,0.979032 ,0.972908 ,0.963286 ,0.965157 ,0.972431 ,0.973907 ,0.978667 ,0.976402 ,0.979853 ,0.980115 ,0.992664 ,0.990575 ,0.985894 ,0.986049 ,0.982318 ,0.974736 ,0.969791 ,0.963541 ,0.965463 ,0.981748 ,0.991875 ,1.03103 ,1.06733 },
		{ 0.99304 ,0.975182 ,0.968496 ,0.960099 ,0.963577 ,0.96193 ,0.975247 ,0.974381 ,0.979061 ,0.986147 ,0.986861 ,0.987487 ,0.9943 ,0.99376 ,0.992214 ,0.983973 ,0.981484 ,0.976992 ,0.970078 ,0.973678 ,0.986794 ,0.999207 ,1.02059 ,1.07428 },
		{ 0.983272 ,0.977947 ,0.975087 ,0.964998 ,0.961232 ,0.97188 ,0.966942 ,0.969954 ,0.978587 ,0.989158 ,0.994775 ,0.996366 ,0.998757 ,0.998174 ,0.991684 ,0.997729 ,0.989186 ,0.984641 ,0.978577 ,0.98513 ,0.995555 ,1.00722 ,1.03962 ,1.06433 },
		{ 0.987462 ,0.977607 ,0.971703 ,0.972301 ,0.962686 ,0.961227 ,0.97098 ,0.980649 ,0.984982 ,0.991906 ,0.997156 ,1.00195 ,1.00466 ,1.00464 ,1.0077 ,1.00109 ,0.997109 ,0.991815 ,0.990239 ,0.998871 ,1.00817 ,1.01455 ,1.03382 ,1.07071 },
		{ 0.982182 ,0.972869 ,0.972723 ,0.965276 ,0.968042 ,0.969771 ,0.98308 ,0.985351 ,0.99607 ,1.00284 ,1.00733 ,1.01064 ,1.01547 ,1.01466 ,1.01267 ,1.00907 ,1.00437 ,1.00357 ,0.998742 ,0.993876 ,1.01404 ,1.02519 ,1.03504 ,1.07209 },
		{ 0.983914 ,0.982093 ,0.970492 ,0.972619 ,0.968742 ,0.980429 ,0.984735 ,0.995658 ,1.00433 ,1.01286 ,1.01654 ,1.02432 ,1.02467 ,1.02651 ,1.02485 ,1.02047 ,1.01638 ,1.01327 ,1.00549 ,1.00168 ,1.01928 ,1.02844 ,1.04067 ,1.07949 },
		{ 0.989 ,0.979226 ,0.974968 ,0.973066 ,0.975886 ,0.980976 ,0.992099 ,1.00366 ,1.0126 ,1.02201 ,1.02997 ,1.03397 ,1.03578 ,1.03497 ,1.0341 ,1.03099 ,1.02762 ,1.01799 ,1.01358 ,1.01624 ,1.02648 ,1.03594 ,1.05524 ,1.08354 },
		{ 0.992498 ,0.981951 ,0.975581 ,0.976662 ,0.979067 ,0.99013 ,1.00083 ,1.01378 ,1.01975 ,1.03029 ,1.03806 ,1.04444 ,1.04704 ,1.04621 ,1.04604 ,1.0451 ,1.03546 ,1.03253 ,1.0279 ,1.02603 ,1.03818 ,1.04426 ,1.0623 ,1.09111 },
		{ 0.980287 ,0.979761 ,0.981388 ,0.979615 ,0.988177 ,0.997178 ,1.01259 ,1.02286 ,1.03141 ,1.03889 ,1.0464 ,1.05323 ,1.05715 ,1.05605 ,1.05331 ,1.05192 ,1.04555 ,1.04454 ,1.03634 ,1.03447 ,1.04328 ,1.0563 ,1.06265 ,1.0831 },
		{ 0.989191 ,0.984749 ,0.980598 ,0.992321 ,0.994191 ,1.00301 ,1.01333 ,1.02702 ,1.03615 ,1.04536 ,1.05472 ,1.06067 ,1.06332 ,1.06641 ,1.0646 ,1.06151 ,1.05578 ,1.05039 ,1.05029 ,1.04369 ,1.05087 ,1.06032 ,1.07427 ,1.08586 },
		{ 0.988192 ,0.982686 ,0.980531 ,0.983856 ,0.995175 ,1.00751 ,1.02257 ,1.0316 ,1.04482 ,1.05199 ,1.06159 ,1.06736 ,1.07116 ,1.07474 ,1.07327 ,1.07054 ,1.0646 ,1.05881 ,1.05675 ,1.05284 ,1.04972 ,1.06534 ,1.07588 ,1.08927 },
		{ 0.988256 ,0.981071 ,0.984533 ,0.993433 ,1.00157 ,1.01869 ,1.02621 ,1.03674 ,1.04766 ,1.05813 ,1.0647 ,1.07657 ,1.0783 ,1.07975 ,1.08154 ,1.07779 ,1.0734 ,1.06726 ,1.06983 ,1.05869 ,1.06255 ,1.06877 ,1.08093 ,1.09761 },
		{ 0.987592 ,0.987137 ,0.988462 ,0.99758 ,1.00712 ,1.01679 ,1.0264 ,1.04066 ,1.05141 ,1.0606 ,1.07235 ,1.07794 ,1.0784 ,1.08254 ,1.0859 ,1.08098 ,1.07341 ,1.07227 ,1.06341 ,1.06317 ,1.0702 ,1.07081 ,1.08534 ,1.09864 },
		{ 0.982334 ,0.985034 ,0.985885 ,0.993895 ,1.00408 ,1.01251 ,1.03018 ,1.03639 ,1.05166 ,1.06113 ,1.06452 ,1.06944 ,1.07836 ,1.07917 ,1.08117 ,1.07354 ,1.07761 ,1.07033 ,1.06788 ,1.06808 ,1.07709 ,1.09373 ,1.09849 ,1.10869 },
		{ 0.980137 ,0.983206 ,0.994643 ,0.994006 ,1.00233 ,1.01713 ,1.02659 ,1.03778 ,1.04879 ,1.05328 ,1.06146 ,1.06763 ,1.0731 ,1.07111 ,1.07682 ,1.07017 ,1.06667 ,1.06462 ,1.06924 ,1.07001 ,1.08471 ,1.08598 ,1.11571 ,1.12233 },
		{ 0.982608 ,0.996479 ,0.999176 ,1.00339 ,1.0099 ,1.01776 ,1.03029 ,1.03535 ,1.04555 ,1.04967 ,1.05944 ,1.06144 ,1.06374 ,1.0697 ,1.0645 ,1.06654 ,1.05646 ,1.06081 ,1.05743 ,1.06261 ,1.06185 ,1.0753 ,1.08249 ,1.08909 },
		{ 0.992601 ,0.997085 ,1.00314 ,1.00493 ,1.00636 ,1.01454 ,1.0218 ,1.03263 ,1.03792 ,1.0397 ,1.04707 ,1.05798 ,1.05094 ,1.06264 ,1.0555 ,1.0621 ,1.05199 ,1.05375 ,1.05388 ,1.05695 ,1.05525 ,1.06294 ,1.06821 ,1.08819 },
		{ 1.00079 ,0.997935 ,1.00136 ,0.994838 ,1.00571 ,1.01707 ,1.0103 ,1.01914 ,1.03083 ,1.04212 ,1.03284 ,1.03221 ,1.04291 ,1.05283 ,1.0429 ,1.04545 ,1.04566 ,1.04218 ,1.04478 ,1.05698 ,1.05052 ,1.06126 ,1.07659 ,1.07924 },
		{ 0.991295 ,0.995319 ,0.995369 ,0.996285 ,0.997468 ,1.00458 ,1.01068 ,1.01241 ,1.02353 ,1.02896 ,1.02915 ,1.03398 ,1.0357 ,1.03052 ,1.03795 ,1.03529 ,1.03538 ,1.02934 ,1.03983 ,1.03863 ,1.04068 ,1.05238 ,1.06832 ,1.07863 },
		{ 1.00091 ,0.998129 ,0.988076 ,0.990407 ,0.99787 ,1.00261 ,1.00813 ,1.00947 ,1.01295 ,1.02113 ,1.02633 ,1.02384 ,1.03044 ,1.02855 ,1.02895 ,1.02491 ,1.03 ,1.02479 ,1.02322 ,1.03032 ,1.03832 ,1.0418 ,1.05483 ,1.06327 },
		{ 0.997697 ,1.00313 ,0.98927 ,0.992368 ,0.990535 ,0.989642 ,1.00504 ,1.00253 ,1.00991 ,1.01419 ,1.01778 ,1.0173 ,1.01805 ,1.02015 ,1.01674 ,1.01585 ,1.01067 ,1.01149 ,1.00816 ,1.01437 ,1.02417 ,1.03766 ,1.05119 ,1.06273 },
		{ 0.999226 ,1.00074 ,0.985401 ,0.986661 ,0.984454 ,0.986371 ,0.989734 ,0.996035 ,0.9966 ,1.00252 ,1.00399 ,1.00194 ,1.00231 ,1.00732 ,1.00458 ,0.999048 ,0.997224 ,0.999706 ,0.997466 ,1.00717 ,1.01156 ,1.02777 ,1.04767 ,1.07252 },
		{ 0.999711 ,0.987442 ,0.986139 ,0.983055 ,0.983977 ,0.982746 ,0.985777 ,0.991024 ,0.98937 ,0.994431 ,0.994506 ,0.999707 ,0.993179 ,0.994636 ,0.991404 ,0.989373 ,0.98693 ,0.984244 ,0.985022 ,0.997963 ,1.00229 ,1.019 ,1.04215 ,1.06123 },
		{ 1.00096 ,0.995955 ,0.981038 ,0.983865 ,0.978107 ,0.977254 ,0.976812 ,0.987334 ,0.987478 ,0.986842 ,0.991762 ,0.991512 ,0.993827 ,0.991213 ,0.989929 ,0.984125 ,0.980334 ,0.973955 ,0.97467 ,0.981495 ,0.997121 ,1.01501 ,1.03593 ,1.05608 }
		};

		int pbin = (int)(((this->mom.Phi()+3.14152)/3.14152/2.0)*(bins_phi));
		if(pbin==bins_phi)
			pbin-=1;
		int ctbin = (int)((cos(this->mom.Theta())+1.0)/2.0*(bins_ctheta));
		
		// parabola-spline-interpolation (1D, in ctheta)
		// ---------------------------------------
		/*
		double ctbin_double = ((cos(this->mom.Theta())+1.0)/2.0*(bins_ctheta));
		double cfac = ps_get_y_at_x(corr_tab[pbin], bins_ctheta, ctbin_double);
		*/
		// bilinear-interpolation
		// ---------------------------------------
		/*
		double ctbin_double = ((cos(this->mom.Theta())+1.0)/2.0*(bins_ctheta));
		double pbin_double = (((this->mom.Phi()+3.14152)/3.14152/2.0)*(bins_phi));
		int x1 = pbin;
		int x2 = pbin + 1;
		int y1 = ctbin;
		int y2 = ctbin + 1;
		
		double cfac_lin1 =  ((double)x2 - pbin_double) * corr_tab[x1][y1] + (pbin_double - (double)x1) * corr_tab[x2 % bins_phi][y1];
		double cfac_lin2 =  ((double)x2 - pbin_double) * corr_tab[x1][y2] + (pbin_double - (double)x1) * corr_tab[x2 % bins_phi][y2];
		double cfac = ((double)y2 - ctbin_double) * cfac_lin1 + (ctbin_double - (double)y1) * cfac_lin2;
		*/
		// no interpolation
		// ---------------------------------------
		double cfac = corr_tab[pbin][ctbin]; 
		
		this->mom = this->mom * cfac;
	}

	void particle_class::process(spectrometer_class *spect) { //, const std::vector<std::vector<double>> &tofs_and_moms, const std::vector<double> & adjustment_coeff) {
		double x = this->raw.data.x;
		double y = this->raw.data.y;
		double t = this->raw.data.tof;

		// first correct raw data with individual correction parameters
		x += cor.dx;
		x *= cor.x_stretch;
		x *= cor.overall_stretch;

		y += cor.dy;
		y *= cor.y_stretch;
		y *= cor.overall_stretch;

		t += cor.dt;
		
		// include jet velocity for ions (angle is with respect to y-axis)
		if(this->raw.m >= 1.0) {
			double sh_vjet = spect->VJet * 1000.0 * t * 1e-9; 
			x -= sh_vjet * sin(spect->AngJet/180.0*PI);
			y -= sh_vjet * cos(spect->AngJet/180.0*PI);
		}

		// calculate momenta
		this->mom.x = calc_px(t, x, y, this->raw.m, this->raw.q, spect->Bfield_ns, spect->Bfield_clockwise);
		this->mom.y = calc_py(t, x, y, this->raw.m, this->raw.q, spect->Bfield_ns, spect->Bfield_clockwise);

		if(this->raw.m < 1.0) { 
			if (spect->electron_side->linear_approximation) {
				//GN
				this->mom.z = spect->electron_side->Efields[0] * (t - spect->MeanTOFe) / 124.38;
				//if(this->mom.z>0)
				//	this->mom.z = sqrt(11.66*2/27.211-this->mom.x*this->mom.x-this->mom.y*this->mom.y);
				//else
				//	this->mom.z = -sqrt(11.66*2/27.211-this->mom.x*this->mom.x-this->mom.y*this->mom.y);
				//GN_!
			} else
				this->mom.z = - tof2mom_3accel(t, spect->electron_side->lengths[0], spect->electron_side->lengths[1], spect->electron_side->lengths[2], spect->electron_side->Efields[0], spect->electron_side->Efields[1], spect->electron_side->Efields[2], fabs(this->raw.q), this->raw.m);
		} else {
			if(spect->ion_side->linear_approximation) {	
				this->mom.z = - spect->ion_side->Efields[0] * this->raw.q * (t - this->t_mean) / 124.38;  // according to Markus' diploma thesis
			} 
/*			else if(spect->ion_side->interpolation_approximation ){
				if (tofs_and_moms.size() > 1) {
					this->mom.z = interpolation2mom(t, tofs_and_moms);
					this->mom.z = this->mom.z + adjustment_coeff[0] * x + adjustment_coeff[1] * (x*x) + adjustment_coeff[2] * y + adjustment_coeff[3] * (y*y) + adjustment_coeff[4] * sqrt((x*x) + (y*y));
				}
				else {
					printf("\nERROR: Interpolation point array is empty!!!\n");
				}

			}*/
			else {
				this->mom.z = tof2mom_3accel(t, spect->ion_side->lengths[0], spect->ion_side->lengths[1], spect->ion_side->lengths[2], spect->ion_side->Efields[0], spect->ion_side->Efields[1], spect->ion_side->Efields[2], this->raw.q, this->raw.m);			
			}
		}
	}

	void particle_class::process(spectrometer_class *spect, double x, double y, double t) { //, const std::vector<std::vector<double>> &tofs_and_moms, const std::vector<double> & adjustment_coeff) { //for ion_matrix (polyatomic class)
		
		//this->raw.data.x / y / tof MUST NOT BE CHANGED inside this function, otherwise polyatomic::sort_ion_matrix will explode!
		// first correct raw data with individual correction parameters
		x += cor.dx;
		x *= cor.x_stretch;
		x *= cor.overall_stretch;

		y += cor.dy;
		y *= cor.y_stretch;
		y *= cor.overall_stretch;

		t += cor.dt;
		
		// include jet velocity for ions (angle is with respect to y-axis)
		if(this->raw.m >= 1.0) {
			double sh_vjet = spect->VJet * 1000.0 * t * 1e-9; 
			x -= sh_vjet * sin(spect->AngJet/180.0*PI);
			y -= sh_vjet * cos(spect->AngJet/180.0*PI);
		}

		// calculate momenta
		this->mom.x = calc_px(t, x, y, this->raw.m, this->raw.q, spect->Bfield_ns, spect->Bfield_clockwise);
		this->mom.y = calc_py(t, x, y, this->raw.m, this->raw.q, spect->Bfield_ns, spect->Bfield_clockwise);

		if(this->raw.m < 1.0) { 
			if(spect->electron_side->linear_approximation) 
				this->mom.z = spect->electron_side->Efields[0] * (t - spect->MeanTOFe) / 124.38;  
			else
				this->mom.z = - tof2mom_3accel(t, spect->electron_side->lengths[0], spect->electron_side->lengths[1], spect->electron_side->lengths[2], spect->electron_side->Efields[0], spect->electron_side->Efields[1], spect->electron_side->Efields[2], fabs(this->raw.q), this->raw.m);
		} else {
			if(spect->ion_side->linear_approximation) {	
				this->mom.z = - spect->ion_side->Efields[0] * this->raw.q * (t - this->t_mean) / 124.38;  // according to Markus' diploma thesis
			} 
			/*else if(spect->ion_side->interpolation_approximation ){
				if (tofs_and_moms.size() != 0) {
					this->mom.z = interpolation2mom(t, tofs_and_moms);
					this->mom.z = this->mom.z + adjustment_coeff[0] * x + adjustment_coeff[1] * (x*x) + adjustment_coeff[2] * y + adjustment_coeff[3] * (y*y) + adjustment_coeff[4] * sqrt((x*x) + (y*y));
				}
				else {
					printf("\nERROR: Interpolation point array is empty!!!\n");
				}

			}	*/		
			else {
				this->mom.z = tof2mom_3accel(t, spect->ion_side->lengths[0], spect->ion_side->lengths[1], spect->ion_side->lengths[2], spect->ion_side->Efields[0], spect->ion_side->Efields[1], spect->ion_side->Efields[2], this->raw.q, this->raw.m);			
			}
		}
	}

	double particle_class::energy() 
	{
		return this->mom.Mag() * this->mom.Mag() / (2.0 * this->raw.m * MASSAU) * EVAU;
	}

	bool particle_class::check_tof(double tof_from, double tof_to, int channel, bool invalidate)
	{
		if((this->raw.data.tof < tof_from) || (raw.data.tof > tof_to))
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(channel>=0)
			this->channel = channel;	
		return true;
	}

	bool particle_class::check_mom(double px_width, double py_width, double pz_width, int channel, bool invalidate)
	{
		if(fabs(mom.x)>px_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(fabs(mom.y)>py_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(fabs(mom.z)>pz_width)
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(channel>=0)
			this->channel = channel;
		return true;
	}


	bool particle_class::check_energy(double energy_from, double energy_to, int channel, bool invalidate)
	{
		if((this->energy() < energy_from) || (this->energy() > energy_to))
		{
			if(invalidate)
				this->valid = false;
			return false;
		}
		if(channel>=0)
			this->channel = channel;
		return true;
	}

	bool particle_class::check_validity(double var, double min, double max, bool negate)
	{
		if(negate) {
			if(var < min || var >= max)
				return false;
			else 
				return this->valid;
		} else {
			if(var > min && var <= max)
				return false;
			else 
				return this->valid;
		}
	}

	void particle_class::invalidate(rtag_struct* rtag) 
	{
		switch (rtag->tag_type) {
		
		case RTAG_X:
			this->valid = check_validity(this->raw.data.x,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_Y:
			this->valid = check_validity(this->raw.data.y,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_TOF:
			this->valid = check_validity(this->raw.data.tof,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PX:
			this->valid = check_validity(this->mom.x,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PY:
			this->valid = check_validity(this->mom.y,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PZ:
			this->valid = check_validity(this->mom.z,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PXY:
			this->valid = check_validity(sqrt(this->mom.x*this->mom.x + this->mom.y*this->mom.y),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PXZ:
			this->valid = check_validity(sqrt(this->mom.x*this->mom.x + this->mom.z*this->mom.z),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PYZ:
			this->valid = check_validity(sqrt(this->mom.y*this->mom.y + this->mom.z*this->mom.z),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_P:
			this->valid = check_validity(this->mom.Mag(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PHI_DEG:
			this->valid = check_validity(this->mom.Phi_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_PHIPOS:
			this->valid = check_validity(this->raw.phi,rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_THETA_DEG:
			this->valid = check_validity(this->mom.Theta_deg(),rtag->min,rtag->max,rtag->negate);
			break;
		case RTAG_ENERGY:
			this->valid = check_validity(this->energy(),rtag->min,rtag->max,rtag->negate);
			break;

		default:
			break;
		}

	}

}
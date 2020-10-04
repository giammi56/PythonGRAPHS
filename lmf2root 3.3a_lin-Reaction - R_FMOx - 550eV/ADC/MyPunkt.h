#pragma once

#ifndef _MyPunkt_h_
#define _MyPunkt_h_

//#include <TObject.h>

class MyPunkt
{
public: 
	MyPunkt():fX(0),fY(0){}
	MyPunkt(double x,double y):fX(x),fY(y){}
public:
	double		&x()				{return fX;}
	double		 x()const			{return fX;}
	double		&y()				{return fY;}
	double		 y()const			{return fY;}

private:
	double fX; 
	double fY;
};

#endif
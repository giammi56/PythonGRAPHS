#include "OS_Version.h"

#ifndef dont_use_MFC

#include <iostream>

	#include "afxwin.h"
	//#include "windows.h"
	#pragma warning(disable : 4996)
	
	/////////////////////////////////////////////////////////////////////////////
	bool init_MFC()
	/////////////////////////////////////////////////////////////////////////////
	{
		// initialize MFC and print and error on failure
		if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
		{
			// TODO: change error code to suit your needs
			std::cerr << _T("Fatal Error: MFC initialization failed") << std::endl;
			return false;
		}
		return true;
	}

#endif
//#include "stdafx.h"
#include "CH_console.h"

namespace CH
{
#ifndef LINUX
	using namespace std;

	void cls()
	{
		gotoXY(0,0);
		for(int i=0 ; i<25 ; i++)
			printf("                                                                                                             ");
		
		gotoXY(0,0);
	}

	void gotoXY(int x, int y)
	{
			//Initialize the coordinates
			COORD coord = {short(x), short(y)};
			//Set the position
			if(!SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord)){

			}
			return;
	}

	void gotoX(int x)
	{
			//Initialize the coordinates
			COORD coord = getXY();
			coord.X = x;

			//Set the position
			if(!SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), coord)){

			}
			return;
	}

	COORD getXY()
	{ 
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		COORD coord = {0, 0};
		if(GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi))
			coord = csbi.dwCursorPosition;
		return coord;
	}


	void Green(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_GREEN);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_GREEN);
	}
	void Red(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED);
	}
	void Blue(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_BLUE);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_BLUE);
	}
	void White(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),15);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),7);
	}
	void Yellow(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_RED|FOREGROUND_GREEN);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_RED|FOREGROUND_GREEN);
	}
	void Cyan(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_BLUE|FOREGROUND_GREEN);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_BLUE|FOREGROUND_GREEN);
	}
	void Magenta(bool highlite) {
		if (highlite) 
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_INTENSITY|FOREGROUND_BLUE|FOREGROUND_RED);
		else
			SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),FOREGROUND_BLUE|FOREGROUND_RED);
	}
	__int32 keyhit()
	{
		if (!_kbhit()) return 0;
		return _getch();
	}

	int getwidth() {
		int columns;
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
		columns = csbi.srWindow.Right - csbi.srWindow.Left + 1;
		return columns;
	}

	int getheight() {
		int rows;
		CONSOLE_SCREEN_BUFFER_INFO csbi;
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
		rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;
		return rows;
	}

	static int current_text_colors = 07;  // presumably...
	int gettextcolor()
	{
	  return current_text_colors & 0x0F;
	}

	int gettextbackground()
	  {
	  return (current_text_colors >> 4) & 0x0F;
	}

	void settextcolor( int color )
	  {
	  current_text_colors = (current_text_colors & 0xF0) | (color & 0x0F);
	  SetConsoleTextAttribute(
		GetStdHandle( STD_OUTPUT_HANDLE ),
		current_text_colors
		);
	}

	void settextbackground( int color )
	  {
	  current_text_colors = (current_text_colors & 0x0F) | ((color & 0x0F) << 4);
	  SetConsoleTextAttribute(
		GetStdHandle( STD_OUTPUT_HANDLE ),
		current_text_colors
		);
	}

#endif


#ifdef LINUX
	void gotoXY(int x, int y) {}

	COORD getXY()
	{
		return 0;
	}

	void Green(bool highlite = false) {
	}

	void Red(bool highlite = false) {
	}

	void White(bool highlite = false) {
	}

	__int32 keyhit(void)
	{
		struct termios term, oterm;
		__int32 fd = 0;
		__int32 c = 0;
		tcgetattr(fd, &oterm);
		memcpy(&term, &oterm, sizeof(term));
		term.c_lflag = term.c_lflag & (!ICANON);
		term.c_cc[VMIN] = 0;
		term.c_cc[VTIME] = 1;
		tcsetattr(fd, TCSANOW, &term);
		c = getchar();
		tcsetattr(fd, TCSANOW, &oterm);
		return ((c != -1) ? c : 0);
	}

#endif
}
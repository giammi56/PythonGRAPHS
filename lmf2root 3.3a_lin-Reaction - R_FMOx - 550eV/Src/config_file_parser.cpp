#include "OS_Version.h"

#include "config_file_parser.h"
#include "stdlib.h"
#ifdef LINUX
#include "string.h"
#endif
#include <sys/stat.h>

#pragma warning(disable : 4996)








#define MAX_NUMBER_OF_BYTES_IN_CONFIG_FILE (200000)




config_file_parser::config_file_parser(std::string name)
{
	init_after_construction((char*)name.c_str());
}




config_file_parser::config_file_parser(const char * input_file)
{
	init_after_construction(input_file);
}



void config_file_parser::init_after_construction(const char * input_file)
{
	file_name[0] = 0;
	if (!FileExists(input_file)) throw;
	number_of_read_bytes = 0;
	all_read_bytes = new char[MAX_NUMBER_OF_BYTES_IN_CONFIG_FILE];
	all_read_bytes[0] = 13; number_of_read_bytes++;
	all_read_bytes[1] = 10; number_of_read_bytes++;
	all_read_bytes[2] = 13; number_of_read_bytes++;
	all_read_bytes[3] = 10; number_of_read_bytes++;
	all_read_bytes[number_of_read_bytes] = 0; // no number_of_read_bytes++; here is ok and needed. 
	number_of_open_files = 0;
	filep = 0;
	error_flag = 0;
	error_text = "";
	external_command = "";
	number_of_external_arguments = 0;
	open_new_file(input_file);
}





config_file_parser::~config_file_parser()
{
	if (all_read_bytes) delete[] all_read_bytes;
	all_read_bytes = 0;
	while(close_last_file());
}






bool config_file_parser::FileExists(const std::string strFilename) {
	return FileExists((const char*)strFilename.c_str());
}

bool config_file_parser::FileExists(const char * strFilename) {
#ifdef LINUX

	if ( access(strFilename,F_OK) )

		return(false);

	else if ( access(strFilename,R_OK) )

	{

		return(false);

	}

	else

		return(true);

#else
	FILE * fi;
	fi = fopen(strFilename,"rb");
	if (!fi) return false;
	if (ferror(fi)) {
		clearerr(fi); fclose(fi);	fi = 0;
		return false;
	}
	char dummy_char;
	fread(&dummy_char,1,1,fi);
	if (ferror(fi)) {
		clearerr(fi); fclose(fi);	fi = 0;
		return false;
	}
	if (fi) {
		fclose(fi);	fi = 0;
	}
	return true;
#endif
}





bool config_file_parser::open_new_file(std::string input_file)
{
	this->remove_trailing_whitespaces(input_file);
	this->strip_quotation_marks(input_file);
	return open_new_file((char*)input_file.c_str());
}


bool config_file_parser::open_new_file(char * input_file)
{
	if (number_of_open_files < MAX_NUMBER_OF_FILES) {
		if (!FileExists(input_file)) return false;
		FILE * f = fopen(input_file,"rt");
		if (!f) return false;

		file_pointer[number_of_open_files] = f;
		file_names[number_of_open_files] = input_file;
		++number_of_open_files;
		filep = file_pointer[number_of_open_files-1];
		line_number[number_of_open_files-1] = 0;
		sprintf(file_name,input_file);
		return true;
	} else return false;
}


bool config_file_parser::close_last_file()
{
	if (number_of_open_files > 0) {
		fclose(file_pointer[number_of_open_files-1]);
		file_pointer[number_of_open_files-1] = 0;
		--number_of_open_files;
		if (number_of_open_files > 0) {
			sprintf(file_name,file_names[number_of_open_files-1].c_str());
			filep = file_pointer[number_of_open_files-1]; 
		} else {filep = 0; file_name[0] = 0;}
		return true;
	}
	return false;
}


void config_file_parser::write_byte_to_config(char c)
{
	if (number_of_read_bytes < MAX_NUMBER_OF_BYTES_IN_CONFIG_FILE-6) {
		all_read_bytes[number_of_read_bytes++] = c;
		all_read_bytes[number_of_read_bytes] = 13; // makes sure that the string is always correctly 0-terminated.
		all_read_bytes[number_of_read_bytes+1] = 10; // makes sure that the string is always correctly 0-terminated.
		all_read_bytes[number_of_read_bytes+2] = 13; // makes sure that the string is always correctly 0-terminated.
		all_read_bytes[number_of_read_bytes+3] = 10; // makes sure that the string is always correctly 0-terminated.
		all_read_bytes[number_of_read_bytes+4] = 0; // makes sure that the string is always correctly 0-terminated.
	}
}


__int32 config_file_parser::getc()
{
	if (!filep) return -1;

	L100:

	__int32 i = fgetc(filep);
	if (feof(filep)) {
		close_last_file();
		if (number_of_open_files > 0) goto L100;
		return -1;
	}
	write_byte_to_config(char(i));

	if (i == '\r') goto L100;
	if (i == '\t') i = ' ';
	if (i == '\n') ++line_number[number_of_open_files-1];

	return i;
}



void config_file_parser::read_line_until_end()
{
	while (true) {
		char c = fgetc(filep);
		if (feof(filep)) break;
		write_byte_to_config(c);
		if (c == '\r') continue;
		if (c == '\n') {++line_number[number_of_open_files-1]; break;}
	}
}



void config_file_parser::read_comment_block()
{
	while (true) {
		char c = fgetc(filep); 
		if (feof(filep)) break;
		write_byte_to_config(c);
		if (c == '\r') continue;
		if (c =='\n') ++line_number[number_of_open_files-1];
		if (c == '*') {
			char c = fgetc(filep);
			if (feof(filep)) break;
			write_byte_to_config(c);
			if (c == '\r') {
				c = fgetc(filep);
				if (feof(filep)) break;
				write_byte_to_config(c);
			}
			if (c == '\n') ++line_number[number_of_open_files-1];
			if (c == '/') break;
		}
	}
}




bool config_file_parser::read_white_spaces(char &c)
{
	c = 0;
	while (true) {
		c = this->getc();
		if (c == -1) return false;
		if (c == '\t') continue;
		if (c == ' ') continue;
		break;
	}
	return true;
}





void config_file_parser::remove_trailing_whitespaces(std::string &arg)
{
	char * st = new char[arg.length()+1];
	sprintf(st,"%s",arg.c_str());
	while (true) {
		__int32 i = strlen(st);
		if (i < 1) break;
		if (st[i-1] == ' ')  {
			st[i-1] = 0;
			continue;
		}
		if (st[i-1] == '\t') {st[i-1] = 0; continue;}
		break;
	}
	
	arg = st;
	sprintf(st,"%s",arg.c_str());
	__int32 i = 0;
	while (true) {
		if (i >= __int32(strlen(st))) break;
		if (st[i] == ' ' || st[i] == '\t') {++i; arg = st+i; continue;}
		break;
	}

	delete[] st; st = 0;
}




//////////////////////////////////////////////////////////
__int32 config_file_parser::GetCurrentLineNumber()
//////////////////////////////////////////////////////////
{
	if (number_of_open_files > 0) return this->line_number[number_of_open_files-1]; else return -1;
}





//////////////////////////////////////////////////////////
bool config_file_parser::read_next_command(std::string &command,std::string &arg)
//////////////////////////////////////////////////////////
{
	char c;
	std::string line = "";
	command = "";
	arg = "";

L100:
	while (true) {
		c = this->getc();
		if (c == -1) goto L200;
		if (c == ';') {
			read_line_until_end();
			continue;
		}
		if (c!='\n' && c!=' ' && c!='\t' && c!=';') break;
	}

	if (c == '/') {
		c = this->getc();
		if (c != '/' && c != '*') {
			error_flag = 1;
			char errort[200];
			sprintf(errort,"error in line %i in file %s",this->line_number[number_of_open_files-1]+1,this->file_name);
			error_text = errort;
			return false;
		}
		if (c == '/') {read_line_until_end(); goto L100;}
		if (c == '*') {read_comment_block();  goto L100;}
	}

	line += c;
	while(true) {
		c = this->getc();
		if (c == -1) break;
		if (c==';') break;
		if (c=='/') {
			error_flag = 1;
			char errort[200];
			sprintf(errort,"error in line %i in file %s",this->line_number[number_of_open_files-1]+1,this->file_name);
			error_text = errort;
			return false;
		}
		line += c;
	}

L200:
	remove_trailing_whitespaces(line);
	command = this->extract_first_string(line);
	remove_trailing_whitespaces(line);
	arg = line;

	return true;
}



std::string config_file_parser::convert_backslash_to_slash(std::string arg)
{
	char * st = new char[arg.length()+1];
	sprintf(st,"%s",arg.c_str());
	for (__int32 i=0;i<__int32(arg.length());++i) {
		if (st[i] == '\\') {
			st[i] = '/';
		}
	}
	std::string return_string = st;
	delete[] st; st = 0;
	return return_string;
}


double config_file_parser::read_double_from_string(std::string &arg)
{
	std::string first = extract_first_string(arg);
	double a = atof(first.c_str());
	return a;
}



void config_file_parser::strip_quotation_marks(std::string & arg)
{
  char * st = new char[arg.length()+1];
  sprintf(st,"%s",arg.c_str());
  while (true) {
	__int32 i = strlen(st);
	if (i<1) break;
	if (st[i-1] == '\'' || st[i-1] == '\"') {st[i-1] = 0; continue;}
	break;
  }
  __int32 i = 0;
  while (true) {
	if (i >= __int32(strlen(st))) break;
	if (st[i] == '\'' || st[i] == '\"') {++i; continue;}
	break;
  }
  arg = st+i;
  delete[] st; st = 0;
}









std::string config_file_parser::extract_first_string(std::string &arg)
{
	remove_trailing_whitespaces(arg);
	char *st = new char[arg.length()+1];
	sprintf(st,"%s",arg.c_str());
	std::string result = "";

	if (strlen(st) > 0) {
		if (st[0] == '\"') {
			__int32 i = 1;
			while (true) {
				if (st[i] == '\"') {arg = st+i+1; break;}
				result += st[i];
				++i;
				if (i >= __int32(strlen(st))) {arg = ""; break;}
			}
		} else {
			__int32 i = -1;
			while (true) {
				++i;
				if (i >= __int32(strlen(st))) {arg = ""; break;}
				if (st[i] == ' ' || st[i] == '\t') {arg = st+i+1; break;}
				result += st[i];
			}
		}
	}

	delete[] st;
	st = 0;
	remove_trailing_whitespaces(arg);
	remove_trailing_whitespaces(result);
	return result;
}








void config_file_parser::run()
{
	std::string command;
	std::string arg;
	external_command = "";
	external_arguments[0] = "";

	while (read_next_command(command,arg)) {
		external_command = "";
		external_arguments[0] = "";
		number_of_external_arguments = 0;

		
		if (!command.compare("execute")) {
#ifdef LINUX
			std::string carg = arg;
			convert_backslash_to_slash(carg);
			if (!this->open_new_file(carg)) printf("WARNING: Could not open file \"%s\"!\n",carg.c_str()); 
#else
			if (!this->open_new_file(arg)) printf("WARNING: Could not open file \"%s\"!\n",arg.c_str()); 
#endif
			continue;
		}

		// must be an external command:
		external_command = command;
		while (true) {
			std::string first =	this->extract_first_string(arg);
			this->remove_trailing_whitespaces(first);
			this->strip_quotation_marks(first);
			if (first.length() == 0) break;
			this->external_arguments[number_of_external_arguments] = first;
			++number_of_external_arguments;
			if (number_of_external_arguments > 100) break;
		}
		break;
	}
}





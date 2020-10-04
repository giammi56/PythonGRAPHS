#pragma once

#ifndef CONFIG_FILE_PARSER_ALREADY_INCLUDED
	#define CONFIG_FILE_PARSER_ALREADY_INCLUDED

	#include "OS_Version.h"

	#define MAX_NUMBER_OF_FILES 400

	#include "fstream"

	class config_file_parser  
	{
	public:
		// constructors:
		config_file_parser(const std::string input_file_name);
		config_file_parser(const char * input_file_name);

		~config_file_parser();

		void init_after_construction(const char * input_file);
		static bool FileExists(const char * strFilename);
		static bool FileExists(const std::string strFilename);
		__int32 GetCurrentLineNumber();
		bool open_new_file(char * input_file_name);
		bool open_new_file(std::string input_file_name);
		static std::string convert_backslash_to_slash(std::string arg);
		bool close_last_file();
		void run();
		__int32 getc();
		bool read_next_command(std::string &command, std::string &arg);
		void read_line_until_end();
		void read_comment_block();
		static void remove_trailing_whitespaces(std::string &arg);
		bool read_white_spaces(char &c);
		void write_byte_to_config(char c);
		static double read_double_from_string(std::string &arg);
		static void strip_quotation_marks(std::string &arg);
		static std::string extract_first_string(std::string &arg);

		__int32 number_of_open_files;
		FILE * file_pointer[MAX_NUMBER_OF_FILES];
		std::string file_names[MAX_NUMBER_OF_FILES];
		FILE * filep;
		__int32 error_flag;
		std::string error_text;
		__int32 line_number[MAX_NUMBER_OF_FILES];
		char file_name[512];
		std::string external_command;
		std::string external_arguments[100];
		__int32 number_of_external_arguments;

		unsigned __int64 number_of_read_bytes;
		char * all_read_bytes;
	};

#endif
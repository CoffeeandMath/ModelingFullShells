#ifndef READ_FILE_H
#define READ_FILE_H

#include <fstream>
#include <iostream>
#include <stdio.h>
using namespace std;
#define MAXLINE 1024

class some_data {

public:
	some_data() {
	}
	;
	~some_data() {
	}
	;


	int refinelevel = 5;
	double L = 1.;


};

class read_file {
public:
	read_file() {
	}
	;
	~read_file() {
	}
	;

	void readInputFile(char *filename, some_data &dat);

private:
	void getNextDataLine(FILE *const filePtr, char *nextLinePtr,
			int const maxSize, int *const endOfFileFlag);
};

#endif

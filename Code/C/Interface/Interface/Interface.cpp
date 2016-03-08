// FileReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>

void loadMatrix(char fname[]);

int main()
{
	char filename[FILENAME_MAX];	// The string we store filenames in

	printf("Please enter the filename of the matrix to load:\n");
	scanf("%s", &filename);

	loadMatrix(filename);

	return 0;
}

void loadMatrix(char fname[])
{
	FILE *myFile;	// Pointer to the data file

	myFile = fopen(fname, "r");	// Opens the file (Creates a handle)

	if (myFile == NULL) perror("Error opening file.");
	else
	{
		char buff[255];	// Creates a buffer to read each line into

		// Read in the file line by line
		while (fgets(buff, 255, myFile) != NULL) {

			printf(buff);
		}
	}
}
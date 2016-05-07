// FileReader.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

//#############################################################################
// Declare Structures:

// A list of available commands to display to the user.
struct commandList {
	int noCommands;		// Number of available command options
	char** comNames;	// List of the names of each command
};

// A rasterized rectangular array 
// to store a floating point matrix row by row
struct Matrix {
	int size,		// The total number of values in the matrix
		stride;		// The number of values per row of the matrix
	float* data;	// The values of the matrix stored row by row
};

// A collection of all the available matrices in memory and their names
struct matrixList {
	int noMatrices;				// Number of matrices in memory
	struct Matrix* matrices;	// The matrices in memory
	char** names;				// The name of each matrix
};

//=============================================================================
// Declare Functions:



struct Matrix *loadCSV(char* fname);
struct matrixList 

//#############################################################################

int main()
{
	//#########################################################################
	// Initialize Variables and Memory:

	char filename[FILENAME_MAX];	// The string we store the filename to load

	//=========================================================================

	// Prepare the lists of Matrices in memory and available commands

	// Request the filename to load the csv data from.
	printf("Please enter the filename of the matrix to load:\n");
	scanf("%s", &filename);

	// Prepare the places in memory to store the float matrix data.
	int stride,	// The number of tokens per row
		size;	// The number of total tokens

	// Load the float matrix from the csv file with its dimensions.
	float* matrix = fCSV(filename,&size,&stride);

	//=========================================================================

	// Print the loaded matrix back to the user.
	if (matrix != NULL) {
		// For each row
		for (int i = 0;i < (size/stride); i++) {
			// Print each column in the row token by token as floats
			for (int j = 0; j < stride; j++) {
				printf("%f	", matrix[(i*stride) + j]);
			}
			// End each column with a newlone so that it displays row by row
			printf("\n");
		}
	}

	return 0;
	//#########################################################################
}

float* fCSV(char* fname, int* sizeHandle, int* strideHandle)
{
	//#########################################################################

	FILE *file;	// Variable for handle to the data file

	file = fopen(fname, "r");	// Opens the file (Creates the handle)

	// If openning the file name gives a NULL handle,
	// we announce it, and skip parsing the matrix.
	if (file == NULL) perror("Error opening file (NULL Handle).");

	//#########################################################################

	// If the handle exists, we proceed to interpret it as a csv file.
	else
	{
		// Initialize our variables:

		char c;					// Current character in the file.

		int tokenCounter = 0,	// Counts the current token in the row.
			charCounter = 0,	// Counts the current character in the token.
			avgTokenLength = 0;	// Counts the average token length to preallocate memory efficiently.

		// The array for the set of tokens as converted to floats
		float* matrixBuffer = (float*)malloc(sizeof(float));
		int matrixBufferSize = 1;	// The number of tokens that can fit in the current matrixBuffer 

		// The array for each token as it is read in as a string
		char* tokenBuffer = (char*)malloc(sizeof(char));
		int tokenBufferSize = 1;	// The number of characters that can fit in the current tokenBuffer

		// Initialize the flags to mark delimeters as actual characters
		bool subDelimiterFlag = false,
			sSubDelimiterFlag = false,

			strideFlag = true;	// Initialize the flag that decides when to count the stride

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		// For each character in the file
		//=====================================================================
		do
		{
			// Get the current character in the file
			c = getc(file);

			//=================================================================
			// If there's an error reading the character, 
			// inform the user and break the loop.
			if (c == ferror(file))
			{
				fprintf(stderr, "Error reading character %d from token %d.\n", charCounter, tokenCounter);
				break;
			}

			//-----------------------------------------------------------------
			// If getc returned EOF,
			// finish the current token, the current row, and close the file:
			if (c == EOF)
			{
				// If this isn't the first character of a token
				if (charCounter != 0) {
					tokenCounter++;		// Increase the token counter.

					// If the new token would overload the matrixBuffer, 
					// we have to resize the matrixBuffer before we do anything.
					if (tokenCounter > matrixBufferSize) {

						// Double the size of the matrixBuffer
						float* temp = (float*)realloc(matrixBuffer, 2 * matrixBufferSize * sizeof(float));

						// If the reassignment worked, copy the temp variable over
						// and double the size variable
						if (temp) {
							matrixBuffer = temp;
							matrixBufferSize *= 2;
						}
						// Else there was an error doubling the size: inform the use and break.
						else {
							fprintf(stderr, "Error doubling the matrixBuffer at token %d.\n", tokenCounter);
							break;
						}
					}

					// Once we have space, 
					// read the token as a float, and copy it into the matrixBuffer
					matrixBuffer[tokenCounter - 1] = strtof(tokenBuffer, NULL);
				}

				// Now we have to end the array

				tokenCounter++;		// Increase the token counter.

				// If the ending token would overload the matrixBuffer, 
				// we have to resize the matrixBuffer before we do anything.
				if (tokenCounter > matrixBufferSize) {

					// Double the size of the matrixBuffer
					float* temp = (float*)realloc(matrixBuffer, 2 * matrixBufferSize * sizeof(float));

					// If the reassignment worked, copy the temp variable over
					// and double the size variable
					if (temp) {
						matrixBuffer = temp;
						matrixBufferSize *= 2;
					}
					// Else there was an error doubling the size: inform the use and break.
					else {
						fprintf(stderr, "Error doubling the size of the matrixBuffer at token %d.\n", tokenCounter);
						break;
					}
				}

				// Set a final NULL character at the end of the array
				matrixBuffer[tokenCounter - 1] = NULL;

				free(tokenBuffer);		// destroy the tokenBuffer,
				fclose(file);			// close the file
				*sizeHandle = tokenCounter-1;	// set the size of the resulting matrix
				return matrixBuffer;	// return the final set of tokens in the matrixBuffer
			}

			//-----------------------------------------------------------------
			// If we hit a newline character outside of subdelimiters
			// We've hit the end of the line:
			if (c == '\n' && !subDelimiterFlag) {

				// To deal with the end of a line:

				// Increase the token counter.
				tokenCounter++;
				// If the stride hasn't been set, do so.
				if (strideFlag) {
					*strideHandle = tokenCounter;
					strideFlag = false;
				};

				// If the new token would overload the matrixBuffer, 
				// we have to resize the matrixBuffer before we do anything.
				if (tokenCounter > matrixBufferSize) {

					// Double the size of the matrixBuffer
					float* temp = (float*)realloc(matrixBuffer, 2 * matrixBufferSize * sizeof(float));

					// If the reassignment worked, copy the temp variable over
					if (temp) {
						matrixBuffer = temp;
						matrixBufferSize *= 2;
					}
					// Else there was an error doubling the size: inform the use and break.
					else {
						fprintf(stderr, "Error doubling the size of the matrixBuffer at token %d.\n", tokenCounter);
						break;
					}
				}

				// Once we have space, 
				// read the token as a float, and copy it into the matrixBuffer
				matrixBuffer[tokenCounter - 1] = strtof(tokenBuffer, NULL);

				// Clear the token buffer
				memset(tokenBuffer, NULL, tokenBufferSize);
				charCounter = 0;
				continue;
			}

			//---------------------------------------------------------------------
			// If we've hit a delimiter outside of subDelimiters,
			// then we've hit the end of the token.
			if (c == ',' && !subDelimiterFlag) {

				tokenCounter++;		// Increase the tokenCounter.

				// make sure the matrixBuffer can handle another token
				
				// If the new token would overload the matrixBuffer, 
				// we have to resize the matrixBuffer before we do anything.
				if (tokenCounter > matrixBufferSize) {

					// Double the size of the matrixBuffer
					float* temp = (float*)realloc(matrixBuffer, 2 * matrixBufferSize * sizeof(float));

					// If the reassignment worked, copy the temp variable over
					// and double the size variable
					if (temp) {
						matrixBuffer = temp;
						matrixBufferSize *= 2;
					}
					// Else there was an error doubling the size: inform the use and break.
					else {
						fprintf(stderr, "Error doubling the size of the matrixBuffer at token %d.\n", tokenCounter);
						break;
					}
				}

				// read the token as a float and copy it into the matrixBuffer
				matrixBuffer[tokenCounter - 1] = strtof(tokenBuffer, NULL);

				// clear the tokenBuffer, reset the character counter, and restart
				memset(tokenBuffer, NULL, tokenBufferSize);
				charCounter = 0;
				continue;
			}

			//-----------------------------------------------------------------
			// If we're at a subDelimiter,
			if (c == '"') {

				// if it's the first of a set:
				if (!subDelimiterFlag) {

					// and if we havn't immediately gone through a set of subDelimiters,
					// turn on the subDelimiterFlag to allow literal commas, and move to the next character.
					if (!sSubDelimiterFlag) {
						subDelimiterFlag = true;
						continue;
					}

					// else this is a subDelimiter right after another subDelimiter,
					// we recognize it as an actual character this cycle,
					// but turn off the flag for allowing literal quotations on the next run.
					else { sSubDelimiterFlag = false; }
				}

				// else if we're at the second of a set of subDelimiters;
				// We reset the subDelimiterFlag,
				// set the potential flag for a quotation, and restart.
				else {
					subDelimiterFlag = false;
					sSubDelimiterFlag = true;
					continue;
				}
			}

			//-----------------------------------------------------------------
			// If we are at anything after a second subDelimiter that is not another subDelimiter,
			// we ignore the potential of literal quotation marks by resetting the flag
			if (sSubDelimiterFlag && (c != '"')) { sSubDelimiterFlag = false; }

			//=================================================================
			// If we did not get an error, hit the end of a line, the end of the file, 
			// or have to deal with delimiters or subDelimiters, then we are at a valid character.

			charCounter++;		// Increase the characterCounter

			// Make sure we can fit another character onto the token

			//If the new character overloads the tokenBuffer, we must resize it
			if (charCounter > tokenBufferSize) {

				// Double the size of the tokenBuffer
				char* temp = (char*)realloc(tokenBuffer, 2 * tokenBufferSize * sizeof(char));

				// If the assignment works, then we can copy the temp variable over.
				if (temp) {
					tokenBuffer = temp;
					tokenBufferSize *= 2;
				}
				// otherwise there was an error, so we inform the user and break.
				else {
					fprintf(stderr, "Error doubling the size of the token buffer for token %d.\n", tokenCounter);
				}
			}

			// Store the character in the tokenBuffer
			tokenBuffer[charCounter - 1] = c;

			//=================================================================

		} while (true);

		// Throw an error and close the file if we somehow end up outside the loop
		fprintf(stderr, "Outside of loop.\n");

		fclose(file);

		return NULL;

		//#####################################################################
	}
}
/* readCSV ********************************************************************

	Opens the specified csv file and parses the result 
	character by character into rows of float values.

parameter input:
	string filePath = directory path of the csv file
output:

******************************************************************************/
/* Pseudocode *****************************************************************
####
Initialization:
	====
	Cactus Arguments:
	====
	Main Arguments:
####
Main Code:
	++++
	Initialize Reader Variables:
	++++
	Read and Format Character by Character:
		====
		EOF Procedure:
			~~~~
			Min/Max Check:
			~~~~
			~~~~
			Update Grid:
			~~~~
		----
		Newline Procedure:
			~~~~
			Min/Max Check:
			~~~~
		----
		Delimiter Procedure:
			~~~~
			Min/Max Check:
			~~~~
		----
		Subdelimiter Procedure:
		----
		Sub-Flag Reset:
		====
		New Character Procedure:
	++++
####
******************************************************************************/
// Initialization: ############################################################

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_WarnLevel.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// Cactus Arguments: ==========================================================

void readCSV(CCTK_ARGUMENTS)
{
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

// Main Arguments: ============================================================

	FILE *file;	// Variable for handle to the data file

// Main Code: #################################################################

	file = fopen(filePath, "r");	// Opens the file (Creates the handle)

	int check = 0;

	// If openning the path given returns a NULL handle,
	// we announce the error, and skip the rest of the code.
	if (file == NULL) CCTK_VWarn(CCTK_WARN_ABORT,__LINE__,__FILE__,
		CCTK_THORNSTRING,"Error opening file '%s' (NULL Handle).",filePath);

	// If it works, we read the file through its handle as a csv file
	else
	{

// Initialize Reader Variables: +++++++++++++++++++++++++++++++++++++++++++++++

		char c;					// Current character in the file.

		int tokenCounter = 0,	// Counts the current token in the row.
			charCounter = 0,	// Counts the current character in the token.
			avgTokenLength = 0,	// Counts the average token length to 
								// preallocate memory efficiently.
			columns=0,			// Counts the number of columns, or components,
								// which should be constant for all the data.
			rows=0;				// Counts the number of entries in the matrix.

		// Minumum, maximum, and spacing values for each coordinate.
		float* localMin = (float*)malloc(6*sizeof(float));
		float* localMax = (float*)malloc(6*sizeof(float));
		float* localDelta = (float*)malloc(6*sizeof(float));
		int lSize = 6;

		// The array for each set of tokens as converted to floats.
		float* matrixBuffer = (float*)malloc(sizeof(float));
		// The number of tokens that fit in the matrixBuffer.
		int matrixBufferSize = 1; 

		// The array for each token as it is read in as a string.
		char* tokenBuffer = (char*)malloc(sizeof(char));
		// The number of characters that fit in the tokenBuffer.
		int tokenBufferSize = 1;

		// No such thing as booleans in C,
		// so we define them as a type of integer.
		typedef int bool;
		#define true 1
		#define false 0
		
		// Initialize the flags to mark delimeters as actual characters.
		bool subDelimiterFlag = false,
			sSubDelimiterFlag = false,

			strideFlag = true;	// Flag that decides when to count the stride.

// Read and Format Character by Character: ++++++++++++++++++++++++++++++++++++
		
		do
		{
			// Update the current character in the file.
			c = getc(file);

			// If there's an error reading the character, 
			// inform the user and break the loop.
			if (c == ferror(file))
			{
				CCTK_VWarn(CCTK_WARN_ABORT,__LINE__,__FILE__,CCTK_THORNSTRING,
					"Error reading character %d from token %d.\n",
					charCounter,tokenCounter);
				break;
			}

// EOF Procedure: =============================================================

			// Finish the current token, the current row, and close the file:
			if (c == EOF)
			{
				// If this isn't the first character of a token
				if (charCounter != 0) {

					// If the new token would overload the matrixBuffer, 
					// we have to resize the matrixBuffer before anything else.
					if (tokenCounter + 1 > matrixBufferSize) {

						// Double the size of the matrixBuffer.
						float* temp = (float*)realloc(matrixBuffer,
							2 * matrixBufferSize * sizeof(float));

						// If the reassignment worked, 
						// copy the temp variable over
						// and double the size variable.
						if (temp) {
							matrixBuffer = temp;
							matrixBufferSize *= 2;
						}
						// Else there was an error doubling the size:
						// inform the use and break.
						else {
							CCTK_VWarn(CCTK_WARN_ABORT,
								__LINE__, __FILE__, CCTK_THORNSTRING,
								"Error doubling matrixBuffer at token %d.\n",
								tokenCounter + 1);
							break;
						}
					}

					// Once we have space, 
					// read the token as a float, 
					// and copy it into the matrixBuffer.
					matrixBuffer[tokenCounter] = strtof(tokenBuffer, NULL);
				}					

				tokenCounter++;		// Increase the token counter.

				// If the stride hasn't been set, do so.
				if (strideFlag) {
					columns = tokenCounter;
					strideFlag = false;

					// Update the local min and max 
					// to have the correct number of components.
					float* temp = 
						(float*)realloc(localMin, columns*sizeof(float));
					if (temp) { localMin = temp; }
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,
							__LINE__, __FILE__, CCTK_THORNSTRING,
							"Error resising localMin.\n");
						break;
					}
					float* temp2 = 
						(float*)realloc(localMax, columns*sizeof(float));
					if (temp2) { localMax = temp2; }
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,
							__LINE__, __FILE__, CCTK_THORNSTRING,
							"Error resizing localMax.\n");
						break;
					}
					float* temp3 =
						(float*)realloc(localDelta, columns*sizeof(float));
					if (temp3) { localDelta = temp3; }
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,
							__LINE__, __FILE__, CCTK_THORNSTRING,
							"Error resizing localDelta.\n");
						break;
					}

// Min/Max Check: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

					// If the token is smaller/larger than the local Min/Max,
					// update the Min/Max according to column number

					// The first row is used for the initial values
					localMin[columns - 1] = matrixBuffer[columns - 1];
					localMax[columns - 1] = matrixBuffer[columns - 1];

				}
				else {
					if (localMin[columns - 1]
								> matrixBuffer[columns - 1]) {
						localMin[columns - 1]
							= matrixBuffer[columns - 1];
					}

					if (localMax[columns - 1]
						< matrixBuffer[columns - 1]) {
						localMax[columns - 1]
							= matrixBuffer[columns - 1];
					}
				}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				// Set the number of rows of data
				rows = tokenCounter / columns;

				// Display final results for Min/Max of each column
				for (int i = 0; i < columns; i++) {
					CCTK_VInfo(CCTK_THORNSTRING,
						"Column(%d): Min = %f, Max = %f\n",
						i, localMin[i], localMax[i]);
				}

// Update Grid: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				
				// set up the spacing between each coordinate
				// to be half the minimum distance between each 
				// of the points' coordinates
	
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				free(tokenBuffer);		// Destroy the tokenBuffer,
				free(matrixBuffer);		// the row buffer,
				fclose(file);			// and close the file.
				break;
			}

// Newline Procedure: ---------------------------------------------------------

			// If we hit a newline character outside of subdelimiters
			// We've hit the end of the line:
			if (c == '\n' && !subDelimiterFlag) {

				// If the new token would overload the matrixBuffer, 
				// we have to resize the matrixBuffer before we do anything.
				if (tokenCounter + 1 > matrixBufferSize) {

					// Double the size of the matrixBuffer
					float* temp = (float*)realloc(matrixBuffer, 
						2 * matrixBufferSize * sizeof(float));

					// If the reassignment worked,
					// copy the temp variable over.
					if (temp) {
						matrixBuffer = temp;
						matrixBufferSize *= 2;
					}
					// Else there was an error doubling the size:
					// inform the use and break.
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,__LINE__,__FILE__,
							CCTK_THORNSTRING,
							"Error doubling matrixBuffer at token '%d'.\n",
							tokenCounter+1);
						break;
					}
				}

				// Once we have space, 
				// read the token as a float, 
				// and copy it into the matrixBuffer.
				matrixBuffer[tokenCounter] = strtof(tokenBuffer, NULL);

				// If the stride hasn't been set, do so.
				if (strideFlag) {
					columns = tokenCounter+1;
					strideFlag = false;

					// Update the local min and max 
					// to have the correct number of components.
					float* temp =
						(float*)realloc(localMin, columns*sizeof(float));
					if (temp) localMin = temp;
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,
							__LINE__, __FILE__, CCTK_THORNSTRING,
							"Error resising localMin.\n");
						break;
					}
					float* temp2 =
						(float*)realloc(localMax, columns*sizeof(float));
					if (temp2) localMax = temp2;
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,
							__LINE__, __FILE__, CCTK_THORNSTRING,
							"Error resizing localMax.\n");
						break;
					}

// Min/Max Check: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

					// Assign them their initial values
					// for the first row of data.
					localMin[columns-1]	= matrixBuffer[columns - 1];
					localMax[columns - 1] = matrixBuffer[columns - 1];

				}else {
					if (((tokenCounter+1)%columns) != 0) {
						CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__,
							CCTK_THORNSTRING,
							"Stride mismatch in line '%d'.\n",
							(tokenCounter+1)/ columns);
					}
					if (localMin[columns - 1]
							> matrixBuffer[tokenCounter]) {
						localMin[columns - 1]
							= matrixBuffer[tokenCounter];
					}
					if (localMax[columns - 1]
						< matrixBuffer[tokenCounter]) {
						localMax[columns - 1]
							= matrixBuffer[tokenCounter];
					}
				}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				tokenCounter++;		// Increase the token counter.

				// Clear the tokenBuffer.
				memset(tokenBuffer, 0, tokenBufferSize);
				charCounter = 0;
				continue;	// Restart.
			}

// Delimeter Procedure: -------------------------------------------------------

			// If we've hit a delimiter outside of subDelimiters,
			// then we've hit the end of the token.
			if (c == ',' && !subDelimiterFlag) {

				// If the new token would overload the matrixBuffer, 
				// we have to resize the matrixBuffer before we do anything.
				if ((tokenCounter+1)> matrixBufferSize) {

					// Double the size of the matrixBuffer
					float* temp = (float*)realloc(matrixBuffer,
						2 * matrixBufferSize * sizeof(float));

					// If the reassignment worked, copy the temp variable over
					// and double the size variable.
					if (temp) {
						matrixBuffer = temp;
						matrixBufferSize *= 2;
					}
					// Else there was an error doubling the size:
					// inform the use and break.
					else {
						CCTK_VWarn(CCTK_WARN_ABORT,__LINE__,__FILE__,
							CCTK_THORNSTRING,
							"Error doubling the matrixBuffer at token %d.\n",
							tokenCounter+1);
						break;
					}
				}

				// Read the token as a float and copy it into the matrixBuffer.
				matrixBuffer[tokenCounter] = strtof(tokenBuffer, NULL);

// Min/Max Check: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				// If the token is smaller/larger than the local Min/Max,
				// update the Min/Max according to column number.

				// The first row is auto assigned as starting values.
				if (strideFlag) {

					// If the new coordinate would overload the locaMin/Max, 
					// resize their size before we do anything else.
					if ((tokenCounter+1)> lSize) {

						// Double the size of the localMin/Max vectors.
						float* temp = (float*)realloc(localMin,
							2 * lSize * sizeof(float));
						float* temp2 = (float*)realloc(localMin,
							2 * lSize * sizeof(float));

						// If the reassignment worked, copy the temp variables
						// and double the size variable.
						if (temp) {
							localMin = temp;
							localMax = temp2;
							lSize *= 2;
						}
						// Else there was an error doubling the size:
						// inform the use and break.
						else {
							CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__,
								CCTK_THORNSTRING,
								"Error doubling the localMin/Max at %d.\n",
								tokenCounter + 1);
							break;
						}
					}

					// Assign the first row of values for localMin/Max
					localMin[tokenCounter] = matrixBuffer[tokenCounter];
					localMax[tokenCounter] = matrixBuffer[tokenCounter];

				}else {
					if (localMin[tokenCounter%columns]
							> matrixBuffer[tokenCounter]) {
						localMin[tokenCounter%columns]
							= matrixBuffer[tokenCounter];
					}
					if (localMax[tokenCounter%columns]
							< matrixBuffer[tokenCounter]) {
						localMax[tokenCounter%columns]
							= matrixBuffer[tokenCounter];
					}
				}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				tokenCounter++;		// Increase the tokenCounter.

				// Clear the tokenBuffer, 
				// reset the character counter, 
				// and restart.
				memset(tokenBuffer, 0, tokenBufferSize);
				charCounter = 0;
				continue;
			}

// Subdelimiter Procedure: ----------------------------------------------------

			if (c == '"') {

				// If it's the first of a set,
				if (!subDelimiterFlag) {

					// and if we havn't immediately gone through
					// a set of subDelimiters, then
					// turn on the subDelimiterFlag to allow literal commas,
					// and move to the next character.
					if (!sSubDelimiterFlag) {
						subDelimiterFlag = true;
						continue;
					}

					// Otherwise this is a subDelimiter 
					// right after another subDelimiter.
					// We recognize it as an actual character this cycle,
					// but turn off the flag for allowing 
					// literal quotations on the next run.
					else { sSubDelimiterFlag = false; }
				}

				// If we're at the second of a set of subDelimiters;
				// We reset the subDelimiterFlag,
				// set the potential flag for a quotation, and restart.
				else {
					subDelimiterFlag = false;
					sSubDelimiterFlag = true;
					continue;
				}
			}

// Sub-Flag Reset: ------------------------------------------------------------

			// If we are at anything after a second subDelimiter
			// that is not another subDelimiter,
			// we ignore the potential of literal quotation marks
			// by resetting the flag.
			if (sSubDelimiterFlag && (c != '"')) { sSubDelimiterFlag = false; }

// New Character Procedure: ===================================================

			// If we did not get an error, hit the end of a line,
			// the end of the file, or have to deal with delimiters 
			// or subDelimiters, then we are at a valid character.

			// If the new character overloads the tokenBuffer,
			// we must resize it.
			if (charCounter + 1 > tokenBufferSize) {

				// Double the size of the tokenBuffer
				char* temp = (char*)realloc(tokenBuffer,
					2 * tokenBufferSize * sizeof(char));

				// If the assignment works, 
				// then we copy the temp variable over.
				if (temp) {
					tokenBuffer = temp;
					tokenBufferSize *= 2;
				}
				// Otherwise there was an error,
				// so we inform the user and break.
				else {
					CCTK_VWarn(CCTK_WARN_ABORT,__LINE__,__FILE__,
						CCTK_THORNSTRING,
						"Error doubling the matrixBuffer at token %d.\n",
						tokenCounter);
					break;
				}
			}

			// Store the character in the tokenBuffer.
			tokenBuffer[charCounter] = c;

			charCounter++;		// Increase the characterCounter.

		} while (true);

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	}

}

//#############################################################################
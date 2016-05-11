/* readCSV ********************************************************************

	Opens the specified csv file and parses the result character by character
	as rows of float values with column components.

parameter input:
	string filePath = directory path of the csv file
output:

******************************************************************************/
/* mergeSort ******************************************************************

	A recursive sorting algorithm that splits an array into halves
	and then reconstructs them. This algorithm is destructive; that is
	it sorts the array passed to it directly rather than returning a
	separate one.

input:
	a = array to sort
	b = a working buffer of the same size needed for storage
	pa = integer array to be sorted with "a" for tracking permutations
	pb = a working buffer for sorting pa along side a
	i = inclusive (zero index) starting place for sort
	f = exclusive number of elements from each array to sort (last index+1)
output:

******************************************************************************/
/* readCSV Pseudocode *********************************************************
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
			Sort & Delta Check:
			~~~~
			Display Stats:
			~~~~
			Cleanup:
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
/* mergeSort Pseudocode *******************************************************
####
Main Code:
	====
	Terminating Case:
	====
	Non-Terminating Case:
	====
	Reconstruction:
####
******************************************************************************/

/* readCSV *******************************************************************/
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
					// update the Min/Max according to column number.

					// The first row is used for the initial values.
					localMin[columns - 1] = matrixBuffer[columns - 1];
					localMax[columns - 1] = matrixBuffer[columns - 1];

				}
				else {
					// Check if the value is the new minimum for the column.
					if (localMin[columns - 1]
								> matrixBuffer[columns - 1]) {
						localMin[columns - 1]
							= matrixBuffer[columns - 1];
					}
					// Check if the value is the new maximum for the column.
					else if (localMax[columns - 1]
						< matrixBuffer[columns - 1]) {
						localMax[columns - 1]
							= matrixBuffer[columns - 1];
					}
				}

				// Set the number of rows of data.
				rows = tokenCounter / columns;

// Sort & Delta Check:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				
				// We need transpose integer arrays to track the permutations
				// of each column individually and serve as output space.
				int* seq = (int*)malloc(rows*sizeof(int));
				int* perm = (int*)malloc(rows*sizeof(int));
				
				for (int i = 0; i < rows; i++) {
					// The extra row of the array is the indexing sequence
					// used as a baseline for the sorting algorithm.
					seq[i] = i;
				}

				// Also we need float arrays to store each column to be sorted.
				float* col = (float*)malloc(rows*sizeof(float));
				float* sorted = (float*)malloc(rows*sizeof(float));

				// Copy each column into a working array and sort.
				for (int j = (columns-1); j >= 0; j--) {

					for (int i = 0; i < rows; i++) {
						col[i] = matrixBuffer[(columns*i) + j];
					}

					mergeSort(col, sorted, seq, perm, 0, rows);
					
					// After each column is sorted, 
					// resort the entire matrixBuffer column by column
					// using the resulting permutation.
					for (int j2 = 0; j2 < columns; j2++) {

						// Copy each column into the working matrix
						for (int i = 0; i < rows; i++) {
							sorted[i] = matrixBuffer[(columns*i) + j2];
						}

						// Search for the smallest difference 
						localDelta[j] = col[1] - col[0];
						for (int i = 1; i < rows; i++) {
							if (localDelta[j] > col[i] - col[i - 1]) {
								localDelta[j] = col[i] - col[i - 1];
							}
						}

						// Replace the entry with the permuted one
						for (int i = 0; i < rows; i++) {
							matrixBuffer[(columns*i) + j2] = sorted[seq[i]];
							seq[i] = i;
						}
					}
				}

// Diplay Stats: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				// Display final results of sort for first & last columns.
				
				CCTK_INFO("First & Last Column:\n");
				for (int i = 0; i < rows; i++) {
					CCTK_VInfo(CCTK_THORNSTRING,
						"Row %d: %f, %f",
						i+1,matrixBuffer[(columns*i)], matrixBuffer[(columns*i) + columns - 1]);
				}

				// Display number of rows and columns.
				CCTK_VInfo(CCTK_THORNSTRING,
					"%d Rows, %d Columns\n", rows, columns);

				// Display final results for Min/Max of each column.
				for (int i = 0; i < columns; i++) {
					CCTK_VInfo(CCTK_THORNSTRING,
						"Column %d: Min = %f, Max = %f, minDelta = %f\n",
						i+1, localMin[i], localMax[i],localDelta[i]);
				}

// Cleanup: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

				free(sorted);			// Destroy the sorting arrays.
				free(col);
				free(perm);
				free(localMin);			// Destroy the local Stats
				free(localMax);	
				free(localDelta);
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

					// The first row is used for the initial values.
					localMin[columns-1]	= matrixBuffer[columns - 1];
					localMax[columns - 1] = matrixBuffer[columns - 1];

				}else {
					// Check to make sure each row has
					// the same length as the first one (constant #columns).
					if (((tokenCounter+1)%columns) != 0) {
						CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__,
							CCTK_THORNSTRING,
							"Stride mismatch in line '%d'.\n",
							(tokenCounter+1)/ columns);
					}
					// Check if the value is the new minimum for the column.
					if (localMin[columns - 1]
							> matrixBuffer[tokenCounter]) {
						localMin[columns - 1]
							= matrixBuffer[tokenCounter];
					}
					// Check if the value is the new maximum for the column.
					else if (localMax[columns - 1]
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

					// The first row is used for the initial values.
					localMin[tokenCounter] = matrixBuffer[tokenCounter];
					localMax[tokenCounter] = matrixBuffer[tokenCounter];

				}else {
					// Check if the value is the new minimum for the column.
					if (localMin[tokenCounter%columns]
							> matrixBuffer[tokenCounter]) {
						localMin[tokenCounter%columns]
							= matrixBuffer[tokenCounter];
					}
					// Check if the value is the new maximum for the column.
					else if (localMax[tokenCounter%columns]
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

/* mergeSort *****************************************************************/
// Main Code: #################################################################

void mergeSort(float* a, float* b, int* pa,int* pb,int i,int f) {
	
// Terminating Case: ==========================================================

	// If the matrix is one element long, then return.
	if ((f - i) < 2) { return; }

// Non-Terminating Case: ======================================================

	// If the matrix is longer than one element,
	// call it in two segments from the middle recursively.
	int m = (i + f) / 2;
	mergeSort(a, b, pa, pb, i, m);	// First half
	mergeSort(a, b, pa, pb, m, f);	// Second half

// Reconstruction: ============================================================
	
	// Once the pieces are broken up
	// then we repiece them back together in sections
	int it = i, mt = m;
	
	// Run through the output matrix and populate it 
	// with the smallest of each section's leading elements in order.
	for (int j = i; j < f; j++) {
		if (it < m && (mt >= f || a[it] < a[mt])) {
			b[j] = a[it];
			pb[j] = pa[it];
			it = it++;
		}
		else {
			b[j] = a[mt];
			pb[j] = pa[mt];
			mt = mt++;
		}
	}
	for (int j = i; j < f; j++) {
		a[j] = b[j];
		pa[j] = pb[j];
	}
}
//#############################################################################
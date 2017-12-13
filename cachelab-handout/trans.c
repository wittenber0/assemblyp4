/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

char transpose_skip2_desc[] = "Transpose skip2";
void transpose_skip2(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;

	for (j = 0; j < M; j++) {
		for (i = 0; i < N; i+= 2) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 1; i < N; i+= 2) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}
}

char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]){


	int block_r, block_c, r, c;

	if (N == 32){
		for (block_r = 0; block_r < 4; block_r++){
			for (block_c = 0; block_c < 4; block_c++){
				for (r = block_r * 8; r < block_r * 8 + 8; r++){
					if (block_r == block_c){
						B[r][r] = A[r][r];
					}
					for (c = block_c * 8; c < block_c * 8 + 8; c++){
						if (r != c){
							B[r][c] = A[c][r];
						}
					}
				}
			}
		}
	}else if (N == 64){
		int r0, r1, r2, r3;
		for (block_r = 0; block_r < 8; block_r++){
			for (block_c = 0; block_c < 8; block_c++){
				for (c = block_c * 8; c < block_c * 8 + 4; c++){
					for (r = block_r * 8; r < block_r * 8 + 4; r++){
						if (r != c){
							B[r][c] = A[c][r];
						}
					}
					if (block_r == block_c){
						B[c][c] = A[c][c];
					}
				}
				for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
					for (r = block_r * 8; r < block_r * 8 + 4; r++){
						if (r != c - 4){
							B[r][c] = A[c - 4][r + 4];
						}
					}
					if (block_c == block_r){
						B[c - 4][c] = A[c - 4][c];
					}
				}
				r = block_r * 8;
				c = block_c * 8 + 4;
				r0 = B[r][c];
				r1 = B[r][c + 1];
				r2 = B[r][c + 2];
				r3 = B[r][c + 3];
				for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
					B[r][c] = A[c][r];
				}
				r = block_r * 8 + 4;
				c = block_c * 8;
				B[r][c] = r0;
				B[r][c + 1] = r1;
				B[r][c + 2] = r2;
				B[r][c + 3] = r3;

				r = block_r * 8 + 1;
				c = block_c * 8 + 4;
				r0 = B[r][c];
				r1 = B[r][c + 1];
				r2 = B[r][c + 2];
				r3 = B[r][c + 3];
				for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
					B[r][c] = A[c][r];
				}
				r = block_r * 8 + 5;
				c = block_c * 8;
				B[r][c] = r0;
				B[r][c + 1] = r1;
				B[r][c + 2] = r2;
				B[r][c + 3] = r3;

				r = block_r * 8 + 2;
				c = block_c * 8 + 4;
				r0 = B[r][c];
				r1 = B[r][c + 1];
				r2 = B[r][c + 2];
				r3 = B[r][c + 3];
				for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
					B[r][c] = A[c][r];
				}
				r = block_r * 8 + 6;
				c = block_c * 8;
				B[r][c] = r0;
				B[r][c + 1] = r1;
				B[r][c + 2] = r2;
				B[r][c + 3] = r3;

				r = block_r * 8 + 3;
				c = block_c * 8 + 4;
				r0 = B[r][c];
				r1 = B[r][c + 1];
				r2 = B[r][c + 2];
				r3 = B[r][c + 3];
				for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
					B[r][c] = A[c][r];
				}
				r = block_r * 8 + 7;
				c = block_c * 8;
				B[r][c] = r0;
				B[r][c + 1] = r1;
				B[r][c + 2] = r2;
				B[r][c + 3] = r3;

				for (r = block_r * 8 + 4; r < block_r * 8 + 8; r++){
					if (block_c == block_r){
						B[r][r] = A[r][r];
					}
					for (c = block_c * 8 + 4; c < block_c * 8 + 8; c++){
						if (r != c){
							B[r][c] = A[c][r];
						}
					}
				}
			}
		}
	}else{
		int temp;
		int d;
		int blockSize = 16;
		for (int blockForCol = 0; blockForCol < M; blockForCol += blockSize)	{
			for (int blockForRow = 0; blockForRow < N; blockForRow += blockSize)	{
				for(r = blockForRow; (r < N) && (r < blockForRow + blockSize); r++){
					for(c = blockForCol; (c < M) && (c < blockForCol + blockSize); c++){
						if (r != c){
							B[c][r] = A[r][c];
						}else{
							temp = A[r][c];
							d = r;
						}
					}
					if(blockForRow == blockForCol){
						B[d][d] = temp;
					}
				}
			}
		}
	}
}

char transpose_submit2_desc[] = "Transpose submission2";
void transpose_submit2(int M, int N, int A[N][M], int B[M][N])
{
	int blockSize; //variable for size of block, used in each of the iterations, N ==32, N ==63 and the else
		int blockForRow, blockForCol; //to iterate over blocks, user in outer loops
		int r, c; //to iterate through each block, used in inner loops
		int temp = 0, d = 0; //d stands for diagonal, temp is just a temporary variable
		int v0,v1,v2,v3,v4; //Variables to be used in the N==64 case for various assignments within it
		/*
		Using blockSize = 8 in this case. Only N == 32 is used in the condition since matrix transpose can
		occur for any a*b and c*a where only a needs to be same and b and c can vary.
		Blocking is used here.
		4 levels of loop sare used here. 2 outer loops iterate accross blocks (in column major iteration) while the 2 inner loops iterate through each block.
		 */
		if (N == 32)
		{
			blockSize = 8;
			for(blockForCol = 0; blockForCol < N; blockForCol += 8)
			{
				for(blockForRow = 0; blockForRow < N; blockForRow += 8)
				{
					for(r = blockForRow; r < blockForRow + 8; r++)
					{
						for(c = blockForCol; c < blockForCol + 8; c++)
						{
							if(r != c)
							{
								B[c][r] = A[r][c];
							}

							else
							{
								temp = A[r][c];
								d = r;
							}
						}
						//We don't move elements on diagonals since we are transposing a square matrix
						if (blockForRow == blockForCol)
						{
							B[d][d] = temp;
						}
					}
				}
			}
		}

		/* Using blockSize = 4 here.
		2 levels of loops are used
		We assign elements in each row individually. Causes reduced missess. */
		else if (N == 64)
		{
			blockSize = 4;
			for(r = 0; r < N; r += blockSize)
			{
				for(c = 0; c < M; c += blockSize)
				{
					/*Elements in A[r][], A[r+1][], A[r+2][] are assigned to the variables for use throughout this loop
					This is becuase we are only allowed to modify the second matrix B but not the matrix A */
					v0 = A[r][c];
					v1 = A[r+1][c];
					v2 = A[r+2][c];
					v3 = A[r+2][c+1];
					v4 = A[r+2][c+2];
					//Elements in B[c+3][] are assigned
					B[c+3][r] = A[r][c+3];
					B[c+3][r+1] = A[r+1][c+3];
					B[c+3][r+2] = A[r+2][c+3];
					//Elements in B[c+2][] are assigned
					B[c+2][r] = A[r][c+2];
					B[c+2][r+1] = A[r+1][c+2];
					B[c+2][r+2] = v4;
					v4 = A[r+1][c+1];
					//Elements in B[c+1][] are assigned
					B[c+1][r] = A[r][c+1];
					B[c+1][r+1] = v4;
					B[c+1][r+2] = v3;
					//Elements in B[c][] are assigned
					B[c][r] = v0;
					B[c][r+1] = v1;
					B[c][r+2] = v2;
					//Elements in row A[r+3][] are assigned to the left out elements in B (where B has r+3)
					B[c][r+3] = A[r+3][c];
					B[c+1][r+3] = A[r+3][c+1];
					B[c+2][r+3] = A[r+3][c+2];
					v0 = A[r+3][c+3];
					//Finally, elements in row B[c+3][] are assigned
					B[c+3][r+3] = v0;
				}
			}
		}

		/* This is the case for a random matrix size. We use blockSize = 16
		2 levels of loops are used to iterate over blocks in column major iteration and 2 levels are used to go through the blocks	*/
		else

		{
			blockSize = 16;

			for (blockForCol = 0; blockForCol < M; blockForCol += blockSize)
			{
				for (blockForRow = 0; blockForRow < N; blockForRow += blockSize)
				{
					/*Since our sizes can be odd, not all blocks will be square. Special case: if (blockForRow + 16 > N), we get an invalid access.
					We also do regular check for i<N and j<M */
					for(r = blockForRow; (r < N) && (r < blockForRow + blockSize); r++)
					{
						for(c = blockForCol; (c < M) && (c < blockForCol + blockSize); c++)
						{
							//row and column are not same
							if (r != c)
							{
								B[c][r] = A[r][c];
							}

							//row and column same
							else
							{
								temp = A[r][c];
								d = r;
							}
						}

						//Row and column number are same in the blocks, diagonal element assigned
						if(blockForRow == blockForCol)
						{
							B[d][d] = temp;
						}
					}
				}
			}
		}
	}







/*
int i;
int j;
int row_block;
int col_block;
int diag = 0;
int temp = 0;


if (N == 32)
{
	//Finding optimal block sizes included some guess-work; increased B-size for larger M & N when M == N proved less 		 //efficient. However, when M != N, i.e. rectangular matrix, larger B was the way to go.
	for (col_block = 0; col_block < N; col_block += 8)
	{
		for (row_block = 0; row_block < N; row_block += 8)
		{

			for (i = row_block; i < row_block + 8; i ++)
			{

				for (j = col_block; j < col_block + 8; j ++)
				{
					if (i != j)
					{
						B[j][i] = A[i][j];
					}
					 else {
						//Reduce misses m < i*j in B by storing in temp instead of missing in B[j][i]
						temp = A[i][j];
						diag = i;
					}
				}

				//Transpose of a square-matrix has a unique property; no need to move elements on the diagonal.

				if (row_block == col_block)
				{
					//Misses in B reduced to m < i
					B[diag][diag] = temp;
				}
			}

		}
	}

}

else if (N == 64)
{

	//Iterate through matrix using column-major iteration over blocks
	for (col_block = 0; col_block < N; col_block += 4)
	{
		for (row_block = 0; row_block < N; row_block += 4)
		{
			//Iterate over each row using row-major iteration
			for (i = row_block; i < row_block + 4; i ++)
			{
				for (j = col_block; j < col_block + 4; j ++)
				{
					if (i != j)
					{
						B[j][i] = A[i][j];
					}
					else {
						//On the diagonal
						temp = A[i][j];
						diag = i;
					}
				}

				if (row_block == col_block)
				{
					B[diag][diag] = temp;
				}
			}

		}
	}


}
else {

	//Iterate through matrix using column-major iteration over blocks
	for (col_block = 0; col_block < M; col_block += 16)
	{
		for (row_block = 0; row_block < N; row_block += 16)
		{
			//Since our sizes are prime, not all blocks will be square sub-matrices
			//Consider corner-case when (row_block + 16 > N) => invalid access. Explicit check for i, j < n, m
			for (i = row_block; (i < row_block + 16) && (i < N); i ++)
			{
				for (j = col_block; (j < col_block + 16) && (j < M); j ++)
				{

					if (i != j)
					{
						B[j][i] = A[i][j];
					}
					else
					{
						temp = A[i][j];
						diag = i;
					}
				}

				if (row_block == col_block) {
					B[diag][diag] = temp;
				}

			}

 		}
	}

}



}


	int k;
	int h;
	int diagonalLine = 0;
	int holderBlock = 0;
	int increment1 = 8;
	int horizBlock;
	int vertiBlock;



	if (N == 32) {

		for (vertiBlock = 0; vertiBlock < N; vertiBlock += increment1) {

			for (horizBlock = 0; horizBlock < N; horizBlock += increment1) {

				for (k = horizBlock; k < (horizBlock + increment1); k++) {
					if (horizBlock == vertiBlock) {
						B[diagonalLine][diagonalLine] = holderBlock;
					}

					for (h = vertiBlock; h < (vertiBlock + increment1); h++) {
						if (k != h) {
							B[h][k] = A[k][h];
						}

						else {
							diagonalLine = k;
							holderBlock = A[k][h];
						}
					}
				}
			}
		}
	}

	if (N == 64) {

		for (vertiBlock = 0; vertiBlock < N; vertiBlock += (increment1 - 4)) {

			for (horizBlock = 0; horizBlock < N; horizBlock += (increment1 - 4)) {

				for (k = horizBlock; k < horizBlock + (increment1 -4); k ++) {

					for (h = vertiBlock; h < vertiBlock + (increment1 -4); h ++) {
						if (k != h) {
							B[h][k] = A[k][h];
						}
						else {

							diagonalLine = k;
							holderBlock = A[k][h];
						}
					}

					if (horizBlock == vertiBlock) {

						B[diagonalLine][diagonalLine] = holderBlock;
					}
				}
			}
		}
	}

	else {
		for (vertiBlock = 0; vertiBlock < M; vertiBlock += (2 * increment1)) {

			for (horizBlock = 0; horizBlock < N; horizBlock += (2 * increment1)) {

				for (k = horizBlock; (k < horizBlock + (increment1 * 2)) && (k < N); k++) {
					if (horizBlock == vertiBlock) {
						B[diagonalLine][diagonalLine] = holderBlock;
					}

					for (h = vertiBlock; (h < vertiBlock + (increment1 * 2)) && (h < M); h++) {
						if (k != h) {
							B[h][k] = A[k][h];
						}

						else {

							diagonalLine = k;
							holderBlock = A[k][h];
						}
					}
				}
			}
		}
	}
}
*/

char transpose_internet_desc[] = "Transpose internet";
void transpose_internet(int M, int N, int A[N][M], int B[M][N])
{
	/*
   we will navigate three cases in this function:
   1) square matrix of size 32
   2) square matrix of size 64
   3) when matrix is anything else (for example: 61*67)
	 */

	int blockSize; //variable for size of block, used in each of the iterations, N ==32, N ==63 and the else
	int blockForRow, blockForCol; //to iterate over blocks, user in outer loops
	int r, c; //to iterate through each block, used in inner loops
	int temp = 0, d = 0; //d stands for diagonal, temp is just a temporary variable
	int v0,v1,v2,v3,v4; //Variables to be used in the N==64 case for various assignments within it
	/*
	Using blockSize = 8 in this case. Only N == 32 is used in the condition since matrix transpose can
	occur for any a*b and c*a where only a needs to be same and b and c can vary.
	Blocking is used here.
	4 levels of loop sare used here. 2 outer loops iterate accross blocks (in column major iteration) while the 2 inner loops iterate through each block.
	 */
	if (N == 32)
	{
		blockSize = 8;
		for(blockForCol = 0; blockForCol < N; blockForCol += 8)
		{
			for(blockForRow = 0; blockForRow < N; blockForRow += 8)
			{
				for(r = blockForRow; r < blockForRow + 8; r++)
				{
					for(c = blockForCol; c < blockForCol + 8; c++)
					{
						if(r != c)
						{
							B[c][r] = A[r][c];
						}

						else
						{
							temp = A[r][c];
							d = r;
						}
					}
					//We don't move elements on diagonals since we are transposing a square matrix
					if (blockForRow == blockForCol)
					{
						B[d][d] = temp;
					}
				}
			}
		}
	}

	/* Using blockSize = 4 here.
	2 levels of loops are used
	We assign elements in each row individually. Causes reduced missess. */
	else if (N == 64)
	{
		blockSize = 4;
		for(r = 0; r < N; r += blockSize)
		{
			for(c = 0; c < M; c += blockSize)
			{
				/*Elements in A[r][], A[r+1][], A[r+2][] are assigned to the variables for use throughout this loop
				This is becuase we are only allowed to modify the second matrix B but not the matrix A */
				v0 = A[r][c];
				v1 = A[r+1][c];
				v2 = A[r+2][c];
				v3 = A[r+2][c+1];
				v4 = A[r+2][c+2];
				//Elements in B[c+3][] are assigned
				B[c+3][r] = A[r][c+3];
				B[c+3][r+1] = A[r+1][c+3];
				B[c+3][r+2] = A[r+2][c+3];
				//Elements in B[c+2][] are assigned
				B[c+2][r] = A[r][c+2];
				B[c+2][r+1] = A[r+1][c+2];
				B[c+2][r+2] = v4;
				v4 = A[r+1][c+1];
				//Elements in B[c+1][] are assigned
				B[c+1][r] = A[r][c+1];
				B[c+1][r+1] = v4;
				B[c+1][r+2] = v3;
				//Elements in B[c][] are assigned
				B[c][r] = v0;
				B[c][r+1] = v1;
				B[c][r+2] = v2;
				//Elements in row A[r+3][] are assigned to the left out elements in B (where B has r+3)
				B[c][r+3] = A[r+3][c];
				B[c+1][r+3] = A[r+3][c+1];
				B[c+2][r+3] = A[r+3][c+2];
				v0 = A[r+3][c+3];
				//Finally, elements in row B[c+3][] are assigned
				B[c+3][r+3] = v0;
			}
		}
	}

	/* This is the case for a random matrix size. We use blockSize = 16
	2 levels of loops are used to iterate over blocks in column major iteration and 2 levels are used to go through the blocks	*/
	else

	{
		blockSize = 16;

		for (blockForCol = 0; blockForCol < M; blockForCol += blockSize)
		{
			for (blockForRow = 0; blockForRow < N; blockForRow += blockSize)
			{
				/*Since our sizes can be odd, not all blocks will be square. Special case: if (blockForRow + 16 > N), we get an invalid access.
				We also do regular check for i<N and j<M */
				for(r = blockForRow; (r < N) && (r < blockForRow + blockSize); r++)
				{
					for(c = blockForCol; (c < M) && (c < blockForCol + blockSize); c++)
					{
						//row and column are not same
						if (r != c)
						{
							B[c][r] = A[r][c];
						}

						//row and column same
						else
						{
							temp = A[r][c];
							d = r;
						}
					}

					//Row and column number are same in the blocks, diagonal element assigned
					if(blockForRow == blockForCol)
					{
						B[d][d] = temp;
					}
				}
			}
		}
	}
}

char transpose_skip3_desc[] = "Transpose skip3";
void transpose_skip3(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;

	for (j = 0; j < M; j++) {
		for (i = 0; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 1; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 2; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

}

char transpose_skip4_desc[] = "Transpose skip4";
void transpose_skip4(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;

	for (j = 0; j < M; j++) {
		for (i = 0; i < N; i+= 4) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 1; i < N; i+= 4) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 2; i < N; i+= 4) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 3; i < N; i+= 4) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

}

/*
char blocking_desc[] = "blocking";
void blocking(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;
	int bsize = 4;
	int en = bsize * (n/bsize);
	int kk, jj;

	for (kk = 0; kk < en; kk += bsize) {
		for (jj = 0; jj < en; jj += bsize) {
			for (i = 0; i < N; i++) {
				for (j = jj; j < jj + bsize; j++) {
					tmp = A[i][j];
					B[j][i] = tmp;
				}
			}
		}
	}

}

 */

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 */

/*
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;

	for (j = 0; j < M; j++) {
		for (i = 0; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 1; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

	for (j = 0; j < M; j++) {
		for (i = 2; i < N; i+= 3) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

}

*/

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
	int i, j, tmp;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++) {
			tmp = A[i][j];
			B[j][i] = tmp;
		}
	}

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
	/* Register your solution function */
	registerTransFunction(transpose_submit, transpose_submit_desc);

	/* Register any additional transpose functions */
	registerTransFunction(trans, trans_desc);

	//registerTransFunction(transpose_skip2, transpose_skip2_desc);

	registerTransFunction(transpose_internet, transpose_internet_desc);

	//registerTransFunction(transpose_internet2, transpose_internet2_desc);

	registerTransFunction(transpose_skip3, transpose_skip3_desc);

	//registerTransFunction(transpose_skip4, transpose_skip4_desc);



}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < M; ++j) {
			if (A[i][j] != B[j][i]) {
				return 0;
			}
		}
	}
	return 1;
}




/* 
 * trans.c - Matrix transpose B = A^T


 *  Claire LASSERRE 115148


 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{

  int blocksize = 0; // size of the block used in the blocking strategy
  int blockrow = 0;
  int blockcolumn = 0;
  
  if((N == 32 && M == 32)){ 
     blocksize = 8; // optimized value, here are some test (blocksize, misses)
                     // (1,1183),(2,695),(4,439),(8,287),(16,1141),(32,1148)
                     // it is logical to find a blocksize of 8 because s=5 E=1 and b=5 which means. We have 32 sets direct mapped and 32 bytes per block. With a blocksize=8, the cache can hold a entirely 8*8 block of data.
     int i = 0;
     int j = 0;
     int tmp=0;
     int diag_value =0;
     int diag_index =0;
     int diag_to_consider =0; 
     
     for(blockrow = 0; blockrow < N; blockrow += blocksize){ //blocking strategy for the lines as in the doc for the lab
        for(blockcolumn = 0; blockcolumn < M; blockcolumn += blocksize){ //blocking strategy for the columns as in the doc for the lab
           for(i = blockrow; i < blockrow + blocksize; i++){
              for(j = blockcolumn; j < blockcolumn + blocksize; j++){
                 if(i > 31 || j > 31){ //this row and column are already transposed
                    continue; //break out of this loop iteration, it is a strategy to to take blocksize value that it is not a divisor of 32
	               }
	               else {
	                  if (i!=j){
		                  tmp = A[i][j];
		                  B[j][i] = tmp;
	                  }
	                  else{ // strategy to overcome the conflict misses along the diagonale
		                  diag_value = A[i][j]; // memorize the value 
		                  diag_index =j; // memorize the index
		                  diag_to_consider=1; // activation of the flag for diagonal update at the end of the 4th loop
	                  }
	               }       
	            }
	            if (diag_to_consider==1){
	               B[diag_index][diag_index]=diag_value;
	               diag_to_consider =0;
	            }
	         }
         }
      }
    }
   

    
    //methode 1 for 64, similar to the method for N==32 but it is not effective as the best block size leade to 17895 misses
    
   if (N == 64 && M == 64) {
     blocksize = 4; // optimized value, here are some test (blocksize, misses)
                     // (1,4723),(2,2835),(4,1795),(8,4635),(16,4651),(32,4659), (64,4663)
     int i = 0;
     int j = 0;
     int tmp=0;
     int diag_value =0;
     int diag_index =0;
     int diag_to_consider =0; 
       for(blockrow = 0; blockrow < N; blockrow += blocksize){
        for(blockcolumn = 0; blockcolumn < M; blockcolumn += blocksize){
           for(i = blockrow; i < blockrow + blocksize; i++){
              for(j = blockcolumn; j < blockcolumn + blocksize; j++){
                 if(i > 63 || j > 63){ //this row and column are already transposed
                    continue; //break out of this loop iteration, it is a strategy to to take blocksize value that it is not a divisor of 32
	               }
	               else {
	                  if (i!=j){
		                  tmp = A[i][j];
		                  B[j][i] = tmp;
	                  }
	                  else { // strategy to overcome the conflict misses along the diagonale
		                  diag_value = A[i][j]; // memorize the value 
		                  diag_index =j; // memorize the index
		                  diag_to_consider=1; // activation of the flag for diagonal update at the end of the 4th loop
	                  }
	               }       
	            }
	            if (diag_to_consider==1){
	               B[diag_index][diag_index]=diag_value;
	               diag_to_consider =0;
	            }
	         }
         }
      }
   }

  

  if(M == 61 && N == 67){
    blocksize = 23; //optimized value, here are some test (blocksize, misses): (4,2421),(5,2290),(8,2112), (12,2051), (14,1989), (15, 2013), (17, 1943), (18, 1953), (19, 1971), (20, 1949), (22, 1953), (23, 1921), (25, 2098), (26, 2192), (27,2286)
   int i = 0;
   int j = 0;
   int tmp=0;
   int diag_value =0;
   int diag_index =0;
   int diag_to_consider =0; 
   for(blockrow = 0; blockrow < N; blockrow += blocksize){ //blocking strategy for the lines as in the doc for the lab
	  for(blockcolumn = 0; blockcolumn < M; blockcolumn += blocksize){ //blocking strategy for the columns as in the doc for the lab
	    for(i = blockrow; i < blockrow + blocksize; i++){
	      for(j = blockcolumn; j < blockcolumn + blocksize; j++){
		      if(i > 66 || j > 60){ //this row and column are already transposed
		         continue; //break out of this loop iteration
		      }
		      else{
		         if (i!=j){
		          tmp = A[i][j];
		          B[j][i] = tmp;
		         }
		         else {// strategy to overcome the conflict misses along the diagonale
		            diag_value = A[i][j];  // memorize the value 
		            diag_index =j;  // memorize the index
		            diag_to_consider=1; // activation of the flag for diagonal update at the end of the 4th loop
		         }
		      }
	      }
	      if (diag_to_consider==1){
		      B[diag_index][diag_index]=diag_value;
		      diag_to_consider =0;
	     }
	   }
	 }
   }	
  }
  
  
}
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


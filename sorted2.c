//gcc -Wall sorted2.c -o sorted2

#include <stdio.h>
#include <stdlib.h>
#define N 4

typedef struct{
  int a;
  int b;
  int c;
  int d;
} my_data_t;

int less (const void* a, const void* b){
  const my_data_t* ptr_a = a;
  const my_data_t* ptr_b = b;
  return (int)(ptr_a->data[0] - ptr_b->data[0]);
}

/*----------------------------------------------------------------------------------*/
int main(int argc, char *argv[]){

  int i,j,k,y;//c,q,w,y;
  int row=10, count=10;//, counted, n_max, freq[row];//, sortsort[4];
  //int gid[][10] = {{2,6,7,1}, {4,2,6,1}, {7,0,7,1}, {8,8,7,1}, {4,4,8,2}, {5,3,7,2}, {2,2,7,3}, {4,2,7,3}, {8,3,7,9}, {9,5,7,10}};
  //int gid_ordsort[10][4] = {{0}};//, snap[10] = {0},
  //int gid_only[10] = {0};
  //int gid_out[10][4] = {{0}};

    
/*----------------------------------------------------------------------------------*/
/*                                            0                                     */
  my_data_t gid [count];
  gid ={2,6,7,1, 4,2,6,1, 7,0,7,1, 8,8,7,1, 4,4,8,2, 5,3,7,2, 2,2,7,3, 4,2,7,3, 8,3,7,9, 9,5,7,10};


  printf ("The given matrix is \n"); //incident matrix
  for (i=0;i<row;++i) {
    //    for (j=0;j<col;++j) {
//printf (" %d",gid[i]);//[j]);
    //}
    printf ("\n");
  }

  /////
  /*NEED FOR PROJECT*/
  y=0;
  for (j = 0; j < row; j++){ //make array out of a column from incident matrix
// gid_only[i] = gid[i+y];
    //    snap[i] = gid[i][3];
    y=y+4;
  }
  printf("\n last column \n ");
  for (k = 0; k < row; k++){
    //printf(" %d\n ", gid_only[i]);//, snap[i]);
  }
  //qsort (gid_only, count, sizeof(*gid_only), less);

  /////
  /*NEED FOR PROJECT*/
    //for (i = 0; i < row; i++){    //sort snap array in ascending order 
    //for ( j = 0; j < row; j++){ //Loop for comparing other values
    //if (gid_only[j] > gid_only[i]){ //Comparing other array elements
    //tmp = gid_only[i];            //Use temporary variable for storing last value
    //gid_only[i] = gid_only[j];    //replacing value
    //    gid_only[j] = tmp;
    //}
    //}
    //}


  ////  printf("\n\nAscending : ");
  ////for (i = 0; i < row; i++){
  //// printf(" %d ", gid_only[i]);
  ////}
  ////printf ("\n\n");

  /////
  /*NEED FOR PROJECT*/
  ////for(i=0; i<count; i++){ //make dummy array for tracking repeats                
  ////freq[i] = -1; //Initially initialize frequencies to -1                          
    //printf("%d ",freq[i]);
  ////}
  // printf ("\n");

  /////
  /*NEED FOR PROJECT BUT CAN COMMENT OUT AFTER 1ST RUN*/
  ////counted = 0;
  ////for(i=0; i<count; i++){/* If duplicate element is found */                      
  ////if(gid_only[i+1]==gid_only[i]){/*Make sure not to count freq of same element again*/
  ////counted = counted + 1;                                                          
  ////freq[i] = 0;                                                                    
  ////}                                                                              
  ////if(freq[i] != 0){/* If frequency of current element is not counted */           
  ////freq[i] = counted + 1;                                                          
  ////counted = 0;                                                                    
  ////}  
  ////printf("%d %d \n",gid_only[i],freq[i]);                                        
  //// }                                                                              
  ////printf ("\n");

  /////                                                                   

  ////for(i=0; i<count; i++){ //frequency of each element in array           
  ////if(freq[i] != 0){
  ////printf("%d occurs %d times\n", gid_only[i], freq[i]);
  ////}
  ////}
  ////printf ("\n");

  /////
  /*NEED FOR PROJECT BUT CAN COMMENT OUT AFTER 1ST RUN*/
  ////n_max = freq[0];
  ////for (c = 1; c < count; c++){
  ////if (freq[c] > n_max){ //find max value an element was ever repeated            
  ////n_max  = freq[c];
      //location = c+1;                                                               
  ////}
  ////}
  ////printf("maximum repeats of an element is %d.\n", n_max);

  /////
  /*NEED FOR PROJECT*/
  // k=0;
  //loop:do { //building new array sorted by column
  //for (i=0;i<row;i++) {
  //  if (gid_only[k] == gid[i][0]) {
  //	for (j=0;j<col;j++) {
  //	  gid_ordsort[k][j] = gid[i][j];
  //	  gid[i][j] = 0;
  //	}
  //	k=k+1;
  //	goto loop;
  //  }
  //}  
  //}while(k<row);
  
  //  printf ("\n");
  //for (i=0;i<row;++i) {
  //for (j=0;j<col;++j) {
  //  printf (" %d",gid_ordsort[i][j]);
  //}
  //printf ("\n");
  //}


  /*new addition*/  
  //for (j=1;j<23;j++) {
  //for (i=0;i<count;i++) {
  //  if (gid_ordsort[i][3] == j) {
  //	for (k=0;k<col;k++) { 
  //	  gid_out[y][k]=gid_ordsort[i][k];
  //	}
  //	y=y+1;
  //  }
  //}
  //}

  // for (i=0;i<row;i++) {                                                                                                                                             
  // loop2:                                                                                                                                                             
  //if (freq[i] !=0 && freq[i] > 1){                                                                                                                                        
  //q = freq[i];                                                                                                                                                         
  //w = i+2-q;                                                                                                                                                          
  //printf (" %d",q);                                                                                                                                     
  //if (gid_ordsort[w-1][3] > gid_ordsort[w][3]){                                                                                                                         
  //for (j=0;j<col;j++) {                                                                                                                                                
  // sortsort[j] = gid_ordsort[w][j];                                                                                                                                   
  // gid_ordsort[w][j] = gid_ordsort[w-1][j];                                                                                                                             
  // gid_ordsort[w-1][j] = sortsort[j];                                                                                                                                  
  //}                                                                                                                                                                     
  //       goto loop2;                                                                                                                                                  
  //}                                                                          
  //}                                                                            
  //}                                                                              

  //  w=0; 
    //    if (freq[i] !=0 && freq[i] > 1){
    //q = freq[i];
    //y=0;
      //      do {
      //for (w = (i+2-q); w<i+1; w++){  
    //if (gid_ordsort[w][3] < gid_ordsort[w-1][3]){
    //	  for (j=0;j<col;j++) {
    //	    sortsort[j] = gid_ordsort[w][j];
    //	    gid_ordsort[w][j] = gid_ordsort[w-1][j];
    //	    gid_ordsort[w-1][j] = sortsort[j];
    //	  }
    //	  y=y+1;
    //	  printf ("y= %d \n",y);
    //	  printf ("w= %d \n",w);	  
    //	}
    //}
      //      }while(y<w);
    //}
  /*new addition*/

  //  printf ("\n");
  //for (i=0;i<row;++i) {
  //for (j=0;j<col;++j) {
  //  printf (" %d",gid_out[i][j]);
  //}
  //printf ("\n");
  //}
  
  return 0; 
}

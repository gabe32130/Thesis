//gcc -Wall sorted.c -o sorted

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

  int i,j,k,c;//,y;//,q,w,y;
  int row=20, col=4, tmp, count=20, counted, n_max, freq[row];//, sortsort[4];
  int gid[][20] = {{1,6,7,1}, {4,2,6,1}, {5,0,7,1}, {8,8,7,1}, {9,4,8,1}, {1,3,7,2}, {3,2,7,2}, {4,2,7,2}, {6,3,7,2}, {1,5,7,3}, {1,6,7,4}, {4,2,6,4}, {5,0,7,4}, {8,8,7,4}, {9,4,8,5}, {1,3,7,6}, {3,2,7,6}, {4,2,7,6}, {6,3,7,6}, {1,5,7,8}};
  int gid_ordsort[20][4] = {{0}};//, snap[10] = {0};
  int gid_only[20] = {0};
  //int gid_out[20][4] = {{0}};

/*---------------------------------------------------------------------------------------------------------------------------------*/
/*                                                                  0                                                              */
  printf ("The given matrix is \n"); //incident matrix
  printf ("id eu grpid snap \n");
  for (i=0;i<row;++i) {
    for (j=0;j<col;++j) {
      printf (" %d",gid[i][j]);
    }
    printf ("\n");
  }

  /////
  /*NEED FOR PROJECT*/
  for (i = 0; i < row; i++){ //make array out of a column from incident matrix
    gid_only[i] = gid[i][0];
    //snap[i] = gid[i][3];
  }
  printf("\n first column \n ");
  for (i = 0; i < row; i++){
    printf(" %d \n ", gid_only[i]);//, snap[i]);
  }

  /////
  /*NEED FOR PROJECT*/
  for (i = 0; i < row; i++){    //sort snap array in ascending order 
    for ( j = 0; j < row; j++){ //Loop for comparing other values
      if (gid_only[j] > gid_only[i]){ //Comparing other array elements
        tmp = gid_only[i];            //Using temporary variable for storing last value
        gid_only[i] = gid_only[j];    //replacing value
        gid_only[j] = tmp;	

      }
    }
  }
  printf("\n\nAscending : ");
  for (i = 0; i < row; i++){
    printf(" %d ", gid_only[i]);
  }
  printf ("\n\n");

  /////
  /*NEED FOR PROJECT*/
  for(i=0; i<count; i++){ //make dummy array for tracking repeats                    
    freq[i] = -1; //Initially initialize frequencies to -1                           
    //printf("%d ",freq[i]);
  }
  // printf ("\n");

  /////
  /*NEED FOR PROJECT BUT CAN COMMENT OUT AFTER 1ST RUN*/
  counted = 0;
  for(i=0; i<count; i++){/* If duplicate element is found */                      
    if(gid_only[i+1]==gid_only[i]){/*Make sure not to count freq of same element again*/
      counted = counted + 1;                                                          
      freq[i] = 0;                                                                    
    }                                                                                 
    if(freq[i] != 0){/* If frequency of current element is not counted */             
      freq[i] = counted + 1;                                                          
      counted = 0;                                                                    
    }  
    printf("%d %d \n",gid_only[i],freq[i]);                                             
  }                                                                               
  printf ("\n");

  /////                                                                   

  for(i=0; i<count; i++){ //frequency of each element in array           
    if(freq[i] != 0){
      printf("%d occurs %d times\n", gid_only[i], freq[i]);
    }
  }
  printf ("\n");

  /////
  /*NEED FOR PROJECT BUT CAN COMMENT OUT AFTER 1ST RUN*/
  n_max = freq[0];
  for (c = 1; c < count; c++){
    if (freq[c] > n_max){ //find max value an element was ever repeated               
      n_max  = freq[c];
      //location = c+1;                                                               
    }
  }
  printf("maximum repeats of an element is %d.\n", n_max);

  /////
  /*NEED FOR PROJECT*/
  k=0;
 loop:do { //building new array sorted by column
    for (i=0;i<row;i++) {
      if (gid_only[k] == gid[i][0]) {
	for (j=0;j<col;j++) {
	  gid_ordsort[k][j] = gid[i][j];
	  gid[i][j] = 0;
	}
	k=k+1;
	goto loop;
      }
    }  
  }while(k<row);
  
  printf ("\n");
  printf ("id eu grpid snap \n");
  for (i=0;i<row;++i) {
    for (j=0;j<col;++j) {
      printf (" %d",gid_ordsort[i][j]);
    }
    printf ("\n");
  }


  /*new addition*/  
  //for (j=1;j<23;j++) {
  //for (i=0;i<count;i++) {
  //  if (gid_ordsort[i][0] == j) {
  //	for (k=0;k<col;k++) { 
  //	  gid_out[y][k]=gid_ordsort[i][k];
  //	}
  //	y=y+1;
  //  }
  //}
  //}

  // for (i=0;i<row;i++){
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

    //w=0;
    //if (freq[i] !=0 && freq[i] > 1){
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

  //  printf ("\n\n\n\n");
  //for (i=0;i<row;++i) {
  //for (j=0;j<col;++j) {
  //  printf (" %d",gid_out[i][j]);
  //}
  //printf ("\n");
  //}
  
  return 0; 
}

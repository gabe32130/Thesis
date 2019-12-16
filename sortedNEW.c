//gcc -Wall sortedNEW.c -o sortedNEW

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

  int i,j,k,c,y;
  int row=80, tmp, count=20, counted, n_max, freq[row];
  int gid[80] = {1,6,17,1,4,2,6,1,5,0,7,1,8,8,5,1,9,4,9,1,1,3,11,2,3,2,7,2,4,2,7,2,6,3,7,2,1,5,7,3,1,6,7,4,4,2,6,4,5,0,7,4,8,8,7,4,9,4,8,5,1,3,12,6,3,2,7,6,4,2,2,6,6,3,7,6,1,5,1,8};
  int gid_ordsort[80] = {0};
  int gid_only[20] = {0};

/*---------------------------------------------------------------------------------------------------------------------------*/
/*                                                            0                                                              */
  printf ("The given matrix is \n"); //incident matrix
  printf ("id eu grpid snap \n");
  y=0;
  for (i=0;i<row;++i) {
      printf (" %d",gid[i]);
      
      if (i==(3+y)){
	printf ("\n");
	y=y+4;
      }
  }

    
  /*NEED FOR PROJECT*/
  y=0;
  for (i = 0; i < count; i++){ //make array out of a column from incident matrix
    gid_only[i] = gid[y];
    y=y+4;
  }
  //printf("\n first column \n ");
  for (i = 0; i < count; i++){
    //printf(" %d \n ", gid_only[i]);
  }


    
  /*NEED FOR PROJECT*/
  for (i = 0; i < count; i++){    //sort snap array in ascending order
    for ( j = 0; j < count; j++){ //Loop for comparing other values
      if (gid_only[j] > gid_only[i]){ //Comparing other array elements
	tmp = gid_only[i];            //Using temporary variable for storing last value
	gid_only[i] = gid_only[j];    //replacing value
	gid_only[j] = tmp;

      }
    }
  }
  printf("\n\nAscending : ");
  for (i = 0; i < count; i++){
    printf(" %d ", gid_only[i]);
  }
  printf ("\n\n");


    
  /*NEED FOR PROJECT*/
  for(i=0; i<count; i++){ //make dummy array for tracking repeats
    freq[i] = -1; //Initially initialize frequencies to -1
    //printf("%d ",freq[i]);
  }
  // printf ("\n");


    
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
    //printf("%d %d \n",gid_only[i],freq[i]);
  }
  printf ("\n");



  for(i=0; i<count; i++){ //frequency of each element in array
    if(freq[i] != 0){
      printf("%d occurs %d times\n", gid_only[i], freq[i]);
    }
  }
  printf ("\n");


    
  /*NEED FOR PROJECT BUT CAN COMMENT OUT AFTER 1ST RUN*/
  n_max = freq[0];
  for (c = 1; c < count; c++){
    if (freq[c] > n_max){ //find max value an element was ever repeated
      n_max  = freq[c];
      //location = c+1;
    }
  }
  printf("maximum repeats of an element is %d.\n", n_max);


    
  /*NEED FOR PROJECT*/
  k=0;
  y=0;
   loop:do { //building new array sorted by column
    for (i=0;i<row;i=i+4) {
      if (gid_only[k] == gid[i]) {
	  gid_ordsort[y] = gid[i];
	  gid_ordsort[y+1] = gid[i+1];
	  gid_ordsort[y+2] = gid[i+2];
	  gid_ordsort[y+3] = gid[i+3];
	  gid[i] = 0;
	  k=k+1;
	  y=y+4;
	  goto loop;
      }
    }
  }while(k<count);

  printf ("\n");
  printf ("id eu grpid snap \n");
  y=0;
  for (i=0;i<row;++i) {
    printf (" %d",gid_ordsort[i]);
    if (i==(3+y)){
      printf ("\n");
      y=y+4;
    }
  }


  //  for (i=0;i<row;++i) {
  //for (j=0;j<col;++j) {
  //  printf (" %d",gid_ordsort[i][j]);
  //}
  //printf ("\n");
  //}

  return 0;
}

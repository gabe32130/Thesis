//gcc -Wall counted.c -o counted

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
  int gid_new[20], freq[20];
  int count=20, i;
  int counted, c, maximum;
    
  gid_new[0]=1;
  gid_new[1]=1;
  gid_new[2]=1;
  gid_new[3]=1;
  gid_new[4]=2;
  gid_new[5]=2;
  gid_new[6]=3;
  gid_new[7]=3;
  gid_new[8]=11;
  gid_new[9]=11;
  gid_new[10]=11;
  gid_new[11]=11;
  gid_new[12]=11;
  gid_new[13]=11;
  gid_new[14]=7;
  gid_new[15]=7;
  gid_new[16]=7;
  gid_new[17]=8;
  gid_new[18]=9;
  gid_new[19]=9;

  /* Input elements in array */
  for(i=0; i<count; i++){/* Initially initialize frequencies to -1 */

    freq[i] = -1;
    printf("%d %d \n",gid_new[i],freq[i]);
  }
  counted = 0;
  for(i=0; i<count; i++){/* If duplicate element is found */
    if(gid_new[i+1]==gid_new[i]){/* Make sure not to count frequency of same element again */
      counted = counted +1;
      freq[i] = 0;
    }
    if(freq[i] != 0){/* If frequency of current element is not counted */
      freq[i] = counted + 1;
      counted = 0;
    }
    printf("%d %d \n",gid_new[i],freq[i]);
  }
  //Print frequency of each element
   printf("\nFrequency of all elements of array : \n");
   for(i=0; i<count; i++){
     if(freq[i] != 0){
       printf("%d occurs %d times\n", gid_new[i], freq[i]);
     }
   }
  maximum = freq[0];
  for (c = 1; c < count; c++){
    if (freq[c] > maximum){
      maximum  = freq[c];
      //location = c+1;
    }
  }
  printf("Maximum repeats of an element is %d.\n", maximum);
  
  return 0; 
}

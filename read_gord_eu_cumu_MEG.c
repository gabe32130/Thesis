//gcc -o read_gord_eu_cumu_MEG read_gord_eu_cumu_MEG.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm
/* read the gord structure and extract the gas particles that ever had received Eu injection at certain time, should be a smaller number than the total gas particles */
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "Tipsy.h"
#include "assert.h"
#include "rpc/xdr.h"
#include "tipsmall_defs3.h"
#include <malloc.h>
#include <rpc/types.h>

int compare_integer(const void *a, const void *b){
  return (*(int*)a - *(int*)b); 
}

int compare_gassmall(const void *a, const void *b){
  struct gasout *ia = (struct gasout *)a;
  struct gasout *ib = (struct gasout *)b;
  return (int)(ia->iord - ib->iord);
}

int unique(int *array, int ntot, int *uniq, int nuniq){
  /* array -- sorted array with ascending values */
  /* ntot -- total number of elements in the array */
  /* uniq -- return value */ 
  int i, j, count;
  if (nuniq == -1){ //number of unique elements
    count = 1;  /* the first element is always unique */
    for(i=0;i<ntot-1;i++){
      if (array[i+1]!=array[i]){
	count = count+1; 
      }
    }
    return count; 
  } else 
    {
      j=0; 
      uniq[j]=array[0]; //makes shorter array with only unique (no repeating elements)
      for(i=0;i<ntot-1; i++){
	if (array[i+1]!=array[i]){
	  j=j+1; 
	  uniq[j]=array[i+1];
	}
      }
    }
  return 0; 
}

/*--------------------------------------------------------------------*/
/*                                BLOCK1                              */
int main(int argc, char *argv[]){
  FILE *fpord, *fpeu, *fpout, *fpprev; 
  //  struct gassmall *gord;
  struct gasout *gout, *gord, *gprev, *gid, *gid_new;
  struct starlog *slog;
  int i,j,k,nord, nfiles=400, norigin, nslog, n_max;
  char filename[200], filename1[200], filename2[200], filename3[200];
  int *gidcount, *gidcount_new, *uniq, nuniq, count, count1, counted, nread, npeu, result, id, n, nwrite, *ns_hist, *freq;
  float eumass, mass, euj, *eu_dwsum, eu_sum, eu_cumu, eu_prev, massloss_prev;
  int ngd, ndark;

  printf("                 1 \n");
  count = 0;  
  for(i=0;i<=nfiles-1;i++){
    //input files here are the output files for step 1, from paint_rp.c, only the path and the extension needs to be changed, please do not change L90Mpc8000_hithres.%05d. because that's the main name and step number.
    //in the extension I used idx2p0 or idx2 below, because it was for a case with dalay time distribution of index n=-2, but for the fidual case the index is -1
    sprintf(filename1, "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/L90Mpc8000_hithres.%05d.painteu.gidlist_idx1p0",i+1);
    fpeu = fopen(filename1, "r"); //opens outputs from step1
    if(fpeu == NULL){
      printf("file %d doesn't exist \n", i+1); //if step1 output for step i dne then print so
      continue; 
    }
    //assert(fpeu !=NULL); 
    fscanf(fpeu,"%d", &npeu);
    printf("number of gas particles painted %d, step %d \n", npeu, i+1);
    count = count + npeu;
    fclose(fpeu);
  }
  nord = count; //total number of painted gas particles
  printf("total number of painted gas particles %d \n", nord);
  gid = (struct gasout *)malloc(nord*sizeof(struct gasout));
  gidcount = (int*)malloc(nord*sizeof(int));

/*--------------------------------------------------------------------*/
/*                                BLOCK2                              */
  printf("                 2 \n");
  count = 0; 
  for(i=0;i<=nfiles-1;i++){
    count1 = 0; 
    // input files are output files from step 2, i.e. from make_ord_sml.c, again only path and extension needs to be changed 
    sprintf(filename, "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/L90Mpc8000_hithres.%05d.gord_eu_idx1",i+1);
    fpord = fopen(filename, "r"); //opens outputs from step2
    if(fpord == NULL){
      printf("file %d doesn't exist \n", i+1); //if step2 output for step i dne then print so (files 0-12)
      continue; 
    }
    fread(&norigin,sizeof(int),1, fpord);
    gord = (struct gasout *)malloc(norigin*sizeof(struct gasout)); //gassmall[id, eu, mass, fmassloss, grpid, snap]
    nread=fread(gord, sizeof(struct gasout), norigin, fpord); //reads in the order of the gassmall structure
    //printf("%d    %g    \n\n\n", gord[3778].id, gord[3778].eu);
    assert(nread == norigin); //this loop is only looking at the gasout eu to see which particles are noticable painted with eu over all 400 snapshots
    for(j=0;j<=norigin-1;j++){ //norigin = 12936000
      if(gord[j].eu_curr > 1e-37){  //if the eu accumulation is greater than 1e-37
	gid[count].iord = j;                      //absolute id of gas particle
	gid[count].id = gord[j].id;               //gas id for the eu that is being saved
	gid[count].eu_curr = gord[j].eu_curr;     //amount of eu a gas particle has at snap i
	gid[count].mass = gord[j].mass;           //amount of mass a gas particle has
	gid[count].fmassloss = gord[j].fmassloss; //amount of mass a gas particle has loss
	gid[count].grpid = gord[j].grpid;         //the group id of a gas particle
	gid[count].snap = gord[j].snap;           //the snapshot

	gidcount[count]=j;
	count=count+1;   //total number of gas particles painted so far in step i
	count1=count1+1; //number of particles painted in this step i
      }
    }
    printf("number of painted particles in this step %d, total number painted %d, step %d \n", count1, count, i+1);
    free(gord);
    fclose(fpord);
  }
  //for (k=0;k<=150;k++){
  //printf("%d    %d    %g    %g    %d \n", gid[k].iord, gid[k].id, gid[k].eu, gid[k].mass, gid[k].snap);
  //}

  //printf("                 2.2 \n");
  printf("number gas particles painted in gord file %d %d \n", nord, count); //out of the total number of painted gas particles (nord=12021984) the total number of gas particles painted so far in step i with eu greater than 1e-37 (count=7491095)
  assert(nord >=count);
  /* trim off */
  gid_new = (struct gasout *)malloc(count*sizeof(struct gasout));
  gidcount_new = (int*)malloc(count*sizeof(int)); //new array of gas id's
  for(i=0;i<=count-1;i++){ //count = number of gas particles of interest
    gid_new[i].iord = gid[i].iord; //makes new shorter array of the particle id's over all 400 snapshots
    gid_new[i].id = gid[i].id;
    gid_new[i].eu_cumu = -1.0;
    gid_new[i].eu_curr = gid[i].eu_curr;
    gid_new[i].eu_sum = -1.0;
    gid_new[i].mass = gid[i].mass;
    gid_new[i].fmassloss = gid[i].fmassloss;
    gid_new[i].grpid = gid[i].grpid;
    gid_new[i].snap = gid[i].snap;

    gidcount_new[i]=gidcount[i];
  }
  //printf("size of gid %d\n", sizeof(gid));
  //printf("size of gid_new %d\n", sizeof(gid_new));
  //for (k=0;k<=150;k++){
  //printf("%d    %g    %d    %d \n", gid_new[k].id, gid_new[k].eu, gid_new[k].grpid, gid_new[k].snap);
  //}
  //printf("\n\n");
  //for (k=count-150;k<=count;k++){
  //printf("%d %g %d\n", gid_new[k].id, gid_new[k].eu, gid_new[k].grpid);
  //printf("%d    %g    %d    %d \n", gid_new[k].id, gid_new[k].eu, gid_new[k].grpid, gid_new[k].snap);
  //}
  free(gid);
  free(gidcount);
  //printf("                 2.3 \n");
  qsort(gid_new,count,sizeof(struct gasout), compare_gassmall); //compares particle id and sorts
  qsort(gidcount_new,count,sizeof(int), compare_integer);
  
  //for (k=0;k<=200;k++){
    //printf("%d %g %d\n", gid_new[k].id, gid_new[k].eu, gid_new[k].grpid);
    //printf("%d     %d     %g     %d     %d\n", gid_new[k].iord, gid_new[k].id, gid_new[k].eu, gid_new[k].grpid, gid_new[k].snap);
  //}
  
  n = -1;
  nuniq = unique(gidcount_new,count, uniq, n); //total number of unique elements
  uniq = (int*)malloc(nuniq*sizeof(int)); //makes shorter array with only unique (no repeating elements)
  result=unique(gidcount_new, count, uniq, nuniq);
  printf("number of unique gas particles gets painted %d \n", nuniq);

/*--------------------------------------------------------------------*/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^SUB-BLOCK2-A^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*find out how many times which gas id repeats and how many times*/
  //printf("                  sub2a \n");
  freq = (int*)malloc(count*sizeof(int));
  for(i=0; i<count; i++){ //make dummy array for tracking repeats
    freq[i] = -1; //Initially initialize frequencies to -1
  }
  //printf("                  sub2a.1 \n"); 
  counted = 0;
  for(i=0; i<count; i++){/* If duplicate element is found */
    if(gidcount_new[i+1]==gidcount_new[i]){/*Make sure not to count freq of same element again*/
      counted = counted + 1;
      freq[i] = 0;
    }
    if(freq[i] != 0){/* If frequency of current element is not counted */
      freq[i] = counted + 1;
      counted = 0;
    }
  }
  //printf("                   sub2a.2 \n");
  for(i=0; i<count; i++){ //frequency of each element in array
    if(freq[i] != 0){
      //printf("%d occurs %d times \n", gidcount_new[i], freq[i]);
    }
  }
  //printf("                   sub2a.3 \n");
  n_max = freq[0];
  for (i = 1; i < count; i++){
    if (freq[i] > n_max){ //find max value an element was ever repeated
      n_max  = freq[i];
      //location = i+1;
    }
  }
  printf("maximum repeats of an element is %d.\n", n_max);
  free(freq);

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^SUB-BLOCK2-B^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*sum the amount of eu a particle with grpid > 1 has for all particles with grpid > 1*/
  //printf("                  sub2b \n");
  eu_dwsum = (float*)malloc(count*sizeof(float));
  for(i=0; i<count; i++){ //make dummy array for tracking repeats
    eu_dwsum[i] = -1.0;
  }
  //printf("                  sub2b.1 \n");
  euj = 0.0;
  for(i=0; i<count; i++){/* If duplicate element is found */
    if(gid_new[i+1].iord==gid_new[i].iord){/*Make sure not to count the same element again*/
      eu_dwsum[i] = 0.0;
      if(gid_new[i].grpid > 1){
	euj = gid_new[i].eu_curr + euj;
	eu_dwsum[i] = euj;
      }
      //printf("%d  %g  %d  %g  %g\n", gid_new[i].iord, gid_new[i].eu_curr, gid_new[i].grpid, eu_dwsum[i], euj);
    }
    if(eu_dwsum[i] == -1.0){/* If frequency of current element is not counted */
      if(gid_new[i].grpid > 1){
	eu_dwsum[i] = gid_new[i].eu_curr + euj;
      }else{
	eu_dwsum[i]=0;
      }
      //printf("%d  %g  %d  %g  %g\n", gid_new[i].iord, gid_new[i].eu_curr, gid_new[i].grpid, eu_dwsum[i],euj);
      euj = 0.0;
    }
  }
  //for (k=80;k<=260;k++){
  //printf("%d  %g  %d  %g\n", gid_new[k].iord, gid_new[k].eu_curr, gid_new[k].grpid, eu_dwsum[k]);
  //printf("%g  %g\n", gid_new[k].eu_sum, gid_new[k].eu_cumu);
  //}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^SUB-BLOCK2-C^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*sum the amount of eu for a particle over each snapshot i*/
/*the eu in this step + the eu from previous step*/
  //printf("                  sub2c \n");
  eu_sum = 0.0;
  for(i=0; i<count; i++){                    /* If duplicate element is found */
    //for (i=80;i<=260;i++){
    if(gid_new[i+1].iord == gid_new[i].iord){/*Make sure not to count the same element again*/
      eu_sum = gid_new[i].eu_curr + eu_sum;
      gid_new[i].eu_sum = eu_sum;
      //printf("if\n");
      //printf("%d  %g  %g  %d\n", gid_new[i].iord, gid_new[i].eu_curr, gid_new[i].eu_sum, gid_new[i].snap);
    }else{                                   /* If frequency of current element is not counted */
      gid_new[i].eu_sum = gid_new[i].eu_curr + eu_sum;
      eu_sum = 0.0;
      //printf("else\n");
      //printf("%d  %g  %g  %d\n", gid_new[i].iord, gid_new[i].eu_curr, gid_new[i].eu_sum, gid_new[i].snap);
    }
  }
  //for (k=80;k<=260;k++){
  //printf("%d  %g  %g  %d\n", gid_new[k].iord, gid_new[k].eu_curr, gid_new[k].eu_sum, gid_new[k].snap);
  //}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^SUB-BLOCK2-D^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*accumulate the amount of eu for a particle over each snapshot i while accounting for mass loss*/
/*the eu in this step + the eu from previous step - the eu lost in previous star formation*/
  //printf("                  sub2d \n");
  eu_cumu = 0.0;
  eu_prev = 0.0;
  massloss_prev = 0.0;
  for(i=0; i<count; i++){                    /* If duplicate element is found */
    //for (i=80;i<=260;i++){
    if(gid_new[i+1].iord == gid_new[i].iord){/*Make sure not to count the same element again*/
      gid_new[i].eu_cumu = (gid_new[i].eu_curr + eu_prev) - (eu_prev * massloss_prev);
      //printf("%d  %g  %g  %g  %g %g\n",gid_new[i].iord, gid_new[i].eu_sum, gid_new[i].eu_cumu, gid_new[i].eu_curr, eu_prev, (eu_prev * massloss_prev));
      eu_prev = gid_new[i].eu_cumu;
      massloss_prev = gid_new[i].fmassloss; /*!!!!!this is getting reset to zero before ELSE does this matter!!!!!*/
    }else{                                   /* If frequency of current element is not counted */
      gid_new[i].eu_cumu = (gid_new[i].eu_curr + eu_prev) - (eu_prev * massloss_prev);
      //printf("%d  %g  %g  %g  %g %g\n",gid_new[i].iord, gid_new[i].eu_sum, gid_new[i].eu_cumu, gid_new[i].eu_curr, eu_prev, (eu_prev * massloss_prev));
      eu_cumu = 0.0;
      eu_prev= 0.0;
      massloss_prev = 0.0;
    }
  }
  //for (k=80;k<=260;k++){
  //printf("%d  %g  %g  %g  %d\n", gid_new[k].iord, gid_new[k].eu_curr, gid_new[k].eu_sum, gid_new[k].eu_cumu, gid_new[k].snap);
  //printf("%d  %08g  %08g  %08g\n",gid_new[k].iord, gid_new[k].eu_sum, gid_new[k].eu_cumu, gid_new[k].fmassloss); 
  //}

/*--------------------------------------------------------------------*/
/*                              BLOCK3.1                              */
  gout = (struct gasout *)malloc(count*sizeof(struct gasout)); //gasout[iord, id, eu, mass, fmassloss, grpid, snap]
  printf("go through each file again... \n");
  printf("                 3 \n");

  for(i=0;i<=nfiles-1;i++){
    //  return 0;
    // this is again the input file from step 2.
    /*So far the gord_euonly_cumu files have been read and show number of gas particles painted and number of painted gas particles in each snapshot. the total number of gas particles recieving eu along with which ones and how many in each snapshot, if a gas particle recieved eu more than once its id will only be saved once in the array (non-repeating numbers). if the gord_euonly_cumu file dne baseline values are set for the gasout structure but if there is a gord_euonly_cumu file the gasout structure will be read in*/
    //printf("                 3.1 \n");
    sprintf(filename, "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/L90Mpc8000_hithres.%05d.gord_eu_idx1",i+1);
    fpord = fopen(filename, "r"); 
    if(fpord == NULL){
      printf("current file %d doesn't exist \n", i+1); //if file dne like 0-12 assign these baseline values
      continue;
    }
    fread(&norigin,sizeof(int),1, fpord);
    gord = (struct gasout *)malloc(norigin*sizeof(struct gasout)); //gassmall[id, eu, mass, fmassloss, *grpid]
    nread=fread(gord, sizeof(struct gasout), norigin, fpord); 
    assert(nread == norigin);
/*BUILD OUTPUT FILE*/
    printf("iord      id          eu_cumu       eu_curr        eu_sum        mass         fmloss   grpid  dwarf  snap\n");
    for(k=0;k<=count-1;k++){ //out array gasout structure
      //for(k=80;k<=260-1;k++){ 
      if(gid_new[k].snap == i+1){
	gout[k].iord = gid_new[k].iord;                     /*absolute ID original snapshot(i) id name for gas particle*/
	gout[k].id = gid_new[k].id;                         /*tipsy ID gas particles id*/
	gout[k].eu_cumu = gid_new[k].eu_cumu;               /*accumulated eu of a particle accounting for massloss over each snapshot*/
	gout[k].eu_curr = gid_new[k].eu_curr;               /*eu of the current time step*/
	gout[k].eu_sum = gid_new[k].eu_sum;                 /*sum of eu for each particle up to snap i*/
	gout[k].mass = gid_new[k].mass;                     /*gas particle mass in new master out array*/
	gout[k].fmassloss = gid_new[k].fmassloss;           /*mass losed due to star formation*/
	gout[k].grpid = gid_new[k].grpid;                   /*gas particles satellite number*/
	gout[k].dwarfrac = eu_dwsum[k] / gid_new[k].eu_sum; /*eu contributions from dwarf galaxies*/
	gout[k].snap = i+1;                                 /*snapshot number of the particle*/

	//if(i == 398){
	//printf("%d   %d   %06g   %06g   %06g   %06g   %06g    %d   %06g   %d\n", gout[k].iord, gout[k].id, gout[k].eu_cumu, gout[k].eu_curr, gout[k].eu_sum, gout[k].mass, gout[k].fmassloss, gout[k].grpid, gout[k].dwarfrac, gout[k].snap);
	//}
      }
    }

    //if(i==99){
    //for(k=80;k<=260-1;k++){
    //printf("iord     id      eu_cumu      eu_curr       eu_sum       mass     fmloss   grpid   dwarf  snap\n");
    //for(k=0;k<=count-1;k++){
    //	if(gout[k].grpid > 1){
    //	  printf("%d  %d  %06g  %06g  %06g  %06g  %06g     %d    %06g    %d\n", gout[k].iord, gout[k].id, gout[k].eu_cumu, gout[k].eu_curr, gout[k].eu_sum, gout[k].mass, gout[k].fmassloss, gout[k].grpid, gout[k].dwarfrac, gout[k].snap);	
    //	}
    //}
    //}

/*--------------------------------------------------------------------*/
/*                              BLOCK3.2                              */
    // this is output file from the current time step, change extensions but just to be consistent with the previous ones   
    //printf("                 3.1 \n");
    sprintf(filename2, "/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/L90Mpc8000_hithres.%05d.gord_euonly_cumu_idx1",i+1);
    fpout = fopen(filename2, "w"); 
    fwrite(&count, sizeof(int), 1, fpout);
    /*fprintf(fpout, "%d \n", count);    
    for (k=0;k<=count-1;k++){
      printf("%d  %d  %g  %g  %g  %g  %g     %d    %g    %d\n", gout[k].iord, gout[k].id, gout[k].eu_cumu, gout[k].eu_curr, gout[k].eu_sum, gout[k].mass, gout[k].fmassloss, gout[k].grpid, gout[k].dwarfrac, gout[k].snap);
    }
    fclose(fpout);*/
    nwrite = fwrite(gout, sizeof(struct gasout), count, fpout); 
    assert(nwrite == count);
    fclose(fpout); 
    printf("Done file %d \n\n", i+1);
    free(gord);
  }

  return 0; 
}

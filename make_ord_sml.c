// gcc -o make_ord_sml make_ord_sml.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm
// Mac version, changes: add #include <stdlib.h>, change %d to %ld. 
//use ./runmake_ord_sml.sh to run all 400 jobs
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "Tipsy.h"
#include "assert.h"
#include "rpc/xdr.h"
#include "tipsmall_defs3.h"
#include <malloc.h>
#include <rpc/types.h>
int main(int argc, char *argv[]){

  FILE *fptipsy, *fpord, *fpout, *fpcool, *fpgrp, *fpeu, *fphis, *fpslog; 
  struct Tipsy parts; 
  struct gassmall *gout; 
  struct starlog *slog;
  XDR xdrs; 
  int i, j, ng, nfiles=400, nd, norigin, nwrite, ntotgrp, ntot, ntoteu, nstart, nend, ns, nslog, nread, n; 
  int *iord, *igrp, *ideu, *sid, *ns_hist, *grpid;
  float *ctime, *masseu, *dist, *massloss; 

  fphis= fopen("/mn/stornext/u3/shens/scratch/rprocess/L90Mpc8000_hithres.ns_hist", "r"); //reads in some file with star info
  assert(fphis != NULL); 
  ns_hist = (int *)malloc(nfiles*sizeof(int));
  fscanf(fphis,"%d", &n); 
  for(i=0;i<=nfiles-1;i++){
    fscanf(fphis, "%d", &ns_hist[i]); 
  }
  fclose(fphis); 

  fpslog=fopen("/mn/stornext/u3/shens/scratch/rprocess/L90Mpc8000_hithres.starlog.bin", "r"); //reads in some star log file
  assert(fpslog !=NULL); 
  fread(&nslog,sizeof(int),1, fpslog);
  printf("total particles in the starlog file: %d \n", nslog);
  slog = (struct starlog *)malloc(nslog*sizeof(struct starlog)); 
  nread=fread(slog, sizeof(struct starlog), nslog, fpslog); 
  assert(nread == nslog);
  fclose(fpslog);
  for (j=0;j<=10;j++){//prints first 10 lines of starlog
    printf("ordstar %d, ordgas %d, x %g, y %g, z %g, vx%g, vy %g, vz %g, mass %g, rho %g T %g \n", slog[j].iordstar, slog[j].iordgas, slog[j].x, slog[j].y, slog[j].z, slog[j].vx, slog[j].vy, slog[j].vz, slog[j].massform,slog[j].rhoform,slog[j].tempform); 
  }
  // the main output with no extensions after the time step, like xxx.00001 
  fptipsy = fopen(argv[1], "r"); //reads in 
  assert(fptipsy != NULL); 
  Tipsy_in(fptipsy, &parts); 
  ng = parts.h.nsph; //total number of gas particles
  nd = parts.h.ndark;//total number of dark matter particles
  ns = parts.h.nstar;//total number of star particles
  fclose(fptipsy); 
  printf("Number of Gas particles:  %d \n", ng);
  printf("total number of particles: %d \n", ng+nd+parts.h.nstar);

  /* determine the start and end ID formed between this step and next step, note that the formation will occur in the future, correpondingly, the gas particle will lose certain Eu mass in next time step, not this one */ 
  for (i=0;i<=nfiles-1;i++){
    if(ns == ns_hist[i]){
      nstart=ns_hist[i];     //start point of some distribution
      nend = ns_hist[i+1]-1; //end point of some distribution 
    }
  }
  // the file with extension of .iord 
  fpord = fopen(argv[2], "rb"); //reads in .iord files
  assert(fpord != NULL); 
  xdrstdio_create(&xdrs, fpord, XDR_DECODE); 
  xdr_int(&xdrs, &ntot); 
  printf("total particles in the iord file: %d \n", ntot);
  
  //return 0;
  iord = (int *)malloc(ntot*sizeof(int)); 
  for (i=0; i<=ntot-1; i++){
    xdr_int(&xdrs, &iord[i]); 
  }  
  fclose(fpord); 
 

  // this is the output file from step 1
  fpeu = fopen(argv[3], "r"); //reads in outputs from step1 .painteu.gidlist_idx1p0
  assert(fpeu != NULL); 
  fscanf(fpeu, "%d \n", &ntoteu);
  printf("total particles in the r-process file: %d \n", ntoteu);/*prints total number of partiples enriched with eu*/
  if (ntoteu == 0){
    printf("no NSNS merger in this timestep \n"); 
    return 0; 
  }
  ideu = (int *)malloc(ntoteu*sizeof(int)); 
  sid = (int *)malloc(ntoteu*sizeof(int));
  dist = (float *)malloc(ntoteu*sizeof(float)); 
  masseu = (float *)malloc(ntoteu*sizeof(float)); 
  
  for (i=0; i<=ntoteu-1; i++){
    fscanf(fpeu, "%d %d %f %f", &ideu[i],&sid[i], &dist[i],&masseu[i]); 
  }
  fclose(fpeu);
  
  for (i=ntoteu-1; i>=ntoteu-10; i--){
    printf("%d, %d, %g, %g \n", ideu[i], sid[i], dist[i],masseu[i]); 
  } 

  fpgrp = fopen(argv[4], "r"); //4th argument added to track halo id
  assert(fpgrp != NULL);
  fscanf(fpgrp, "%d \n", &ntotgrp);
  printf("total particles in the group id file: %d \n", ntotgrp);/*prints total number of particles in group file*/
  assert(ntotgrp == ntot);
  grpid=(int *)malloc(ntotgrp*sizeof(int));

  for (i=0; i<=ntotgrp-1; i++){
    fscanf(fpgrp, "%d", &grpid[i]);
  }
  fclose(fpgrp);

  printf("the last 50 group ID \n");
  for (i=ntotgrp-1; i>=ntotgrp-50; i--){
    printf("%d \n", grpid[i]); //Print the last 50 group ID
  }

  //   return 0;
  /* fpcool = fopen(argv[3], "rb"); 
  assert(fpcool != NULL); 
  xdrstdio_create(&xdrs, fpcool, XDR_DECODE); 
  xdr_int(&xdrs, &ntot); 
  printf("total particles in the cooltime file: %d \n", ntot); 
  
  ctime = (float *)malloc(ntot*sizeof(float)); 
  for (i=0; i<=ntot-1; i++){
    xdr_float(&xdrs, &ctime[i]); 
  }
  
  fclose(fpcool); */

  /*  for (i=ng-100; i<=ng-1; i++){
    printf("%g \n",ctime[i]);
    } */

  // return 0; 
  
  //norigin = iord[ng-1]+1;
  // norigin = 16777216;  
  //  norigin = 2097152;
    norigin = 12936000;  //Eris
  // norigin = 1583603;
  // norigin = 6076728;  //dwarf
 
  gout = (struct gassmall *)malloc(norigin *sizeof(struct gassmall)); 
  for (i=0; i<=norigin-1; i++){
    gout[i].id = -1; 
    gout[i].eu = 0;  //mass of Europium received in the current timestep
    //gout[i].eu_curr = 0;
    gout[i].mass = 0; 
    gout[i].fmassloss = 0;
    gout[i].grpid = -1;    //halo id
    //gout[i].snap = 0;
  }

  for (i=0; i<=ng-1; i++){
    gout[iord[i]].id = i;
    gout[iord[i]].grpid = grpid[i]; 
    /* for (j=0;j<=ntoteu-1;j++){
      if (i==ideu[j]){
	gout[iord[i]].eu = gout[iord[i]].eu + masseu[j];	
      }
    }*/
    gout[iord[i]].mass = parts.gp[i].mass;
  }



  // printf("########## TEST ########## \n");
  //printf("          last 20 \n");
  //for (i=ntotgrp-1; i>=ntotgrp-20; i--){//ntotgrp=60461129
  //  printf("%d \n", grpid[i]);
  //}
  //printf("          middle 20 \n");
  //for (i=30230564; i<=30230584; i++){
  //  printf("%d \n", grpid[i]);
  //}
  //printf("          first 20 \n");
  //for (i=0; i<=20; i++){
  //  printf("%d \n", grpid[i]);
  //}

  //printf("\n");
  //printf("          first 20 \n");
  //for (i=0; i<=20; i++){
  //  printf("%d, %d \n", gout[i].grpid, grpid[i]);
  //}
  //printf("          middle 20 \n");
  //for (i=5899534; i<=5899554; i++){
  //  printf("%d, %d \n", gout[i].grpid, grpid[i]);
  //}
  //printf("          last 20 \n");
  //for(i=norigin-1; i>=norigin-21; i--){//norigin=11799068
  //  printf("%d, %d \n", gout[i].grpid, grpid[i]);
  //}

  //printf("############ END TEST ############ \n");



  for(j=0;j<=ntoteu-1;j++){
    gout[iord[ideu[j]]].eu = gout[iord[ideu[j]]].eu + masseu[j]; /*in case the same gas particle receive twice of r-process within a single timestep*/
  }
  massloss = (float *)malloc(norigin*sizeof(float)); 

  for(i=0;i<=norigin-1;i++){
    massloss[i] = 0.; 
  }

  /* in case a gas particle sprout > 1 stellar particles between the two output */
  for(j=nstart; j<=nend-1; j++){
    massloss[slog[j].iordgas] = massloss[slog[j].iordgas]+slog[j].massform; 
  }

  /* predicting the FUTURE step mass loss, however there is still uncertainties  because gas mass may grow during the two snapshots, thus if the ratio is larger than 1, set it to 1*/
  for(i=0;i<=norigin-1;i++){
    
    if(gout[i].mass > 0.) gout[i].fmassloss = massloss[i]/gout[i].mass; 
    if (gout[i].fmassloss > 1.0) gout[i].fmassloss = 1.0;
    //if(gout[i].fmassloss !=0.) printf("mass loss fraction %g \n",gout[i].fmassloss);
  }
 
    for(i=0; i<=20; i++){
    printf("first 20 data \n");
    printf("%d, %g, %g, %g, %d \n", gout[i].id, gout[i].eu, gout[i].mass, gout[i].fmassloss, gout[i].grpid); //added gout[i].grpid as 5th output 
  }
  
  for(i=norigin-1; i>=norigin-21; i--){
    printf("last 20 data \n");
    printf("%d, %g, %g, %g, %d \n", gout[i].id, gout[i].eu, gout[i].mass, gout[i].fmassloss, gout[i].grpid); //added gout[i].grpid as 5th output 
  }
  printf("write to file \n");
  // this is output 
  fpout = fopen(argv[5], "w"); //Output argument
  assert(fpout != NULL); 
  fwrite(&norigin, sizeof(int), 1, fpout);
  //  fprintf(fpout, "%d \n", norigin); 
  nwrite = fwrite(gout, sizeof(struct gassmall), norigin, fpout); 
  assert(nwrite == norigin); 
  // for (i=0; i<=norigin-1; i++){
  //  fwrite(gout[i], sizeof(struct gassmall), 1, fpout); 
  //  fprintf(fpout, "%d %g %g \n",gout[i].id, gout[i].eu, gout[i].mass); 
  //  } 
  fclose(fpout); 

  return 0; 
}

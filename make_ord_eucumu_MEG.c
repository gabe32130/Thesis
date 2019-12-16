//gcc -o make_ord_eucumu_MEG make_ord_eucumu_MEG.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm 
//Mac version, changes: add #include <stdlib.h>, change %d to %ld. 
#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "Tipsy.h"
#include "assert.h"
#include "rpc/xdr.h"
#include "tipsmall_defs3.h"
#include <malloc.h>
#include <rpc/types.h>

#define pi 3.14159265358979323846
#define msolunit 1.078e17 
#define kpcunit 90000.

int compare_function(const void *a, const void *b){
  //float x = *(float *)a; 
  //float y = *(float *)b; 
  float x =((const struct gassort*)a)->dist; 
  float y =((const struct gassort*)b)->dist; 
  if (x < y) return -1; 
  else if (x > y) return 1; 
  return 0; 
}

int main(int argc, char *argv[]){
  FILE *fpord, *fpout,*fpslog, *fpeu, *fptipsy, *fphis, *fpiord, *fpparam;
  struct Tipsy parts;
  //struct gassmall *gord;
  struct gasout *geu, *gasout, *gord;
  struct starlog *slog;
  struct gassort *gas, *gas1;
  XDR xdrs;
  int i,j,k,l,ng,nd,ns,nslog,norigin,nwrite,ntot,ntoteu,nread,count,nfiles=400,n,nstart,nend,ns_step,nsm=128,count2,count1,giord,gid,neu,nuniq;
  int *iord, *igrp, *ideu, *sid, *ns_hist;
  float *nseu_out, *nseu_sm, dist, search_r, xc, yc, zc, z, mtot, meu_tot; 
  float *dwarf, meu_totdwarf, eu_dwarf, *dwarfdirect;

/*--------------------------------------------------------------------*/
/*                                BLOCK0                              */
/*read in the starformation history and star logs*/
  printf("Reading in stellar info \n");
  fphis= fopen("/mn/stornext/u3/shens/scratch/rprocess/L90Mpc8000_hithres.ns_hist", "r"); 
  assert(fphis != NULL); 
  ns_hist = (int *)malloc(nfiles*sizeof(int));
  fscanf(fphis,"%d", &n); 
  for(i=0;i<=nfiles-1;i++){
    fscanf(fphis, "%d", &ns_hist[i]); 
  }
  fclose(fphis); 

  fpslog=fopen("/mn/stornext/u3/shens/scratch/rprocess/L90Mpc8000_hithres.starlog.bin", "r"); 
  assert(fpslog !=NULL); 
  fread(&nslog,sizeof(int),1, fpslog);
  printf("total particles in the starlog file: %d \n", nslog);
  slog = (struct starlog *)malloc(nslog*sizeof(struct starlog)); 
  nread=fread(slog, sizeof(struct starlog), nslog, fpslog); 
  assert(nread ==nslog); 
  fclose(fpslog);
  for (j=nslog-3;j<=nslog-1;j++){
    printf("ordstar %d, ordgas %d, x %g, y %g, z %g, vx%g, vy %g, vz %g, mass %g, rho %g T %g \n", slog[j].iordstar, slog[j].iordgas, slog[j].x, slog[j].y, slog[j].z, slog[j].vx, slog[j].vy, slog[j].vz, slog[j].massform,slog[j].rhoform,slog[j].tempform);
  }
  //return 0; 
  //fplog = fopen(argv[6], "w"); 
  //assert(fplog != NULL);
  //fpparam = fopen(argv[1], "r");
  //assert(fpparam !=NULL); 
  //fscanf(fpparam, "%d", &nsm); 
  //printf("nsm = %d \n", nsm);
  //return 0;

/*--------------------------------------------------------------------*/
/*                                BLOCK1                              */
  printf("                 1 \n");
  fptipsy = fopen(argv[1], "r");//reads in simulation files L90Mpc8000_hithres
  assert(fptipsy != NULL); 
  Tipsy_in(fptipsy, &parts); 
  ng = parts.h.nsph;
  nd = parts.h.ndark;
  ns = parts.h.nstar;
  fclose(fptipsy);
  printf("Number of Gas particles:  %d \n", ng);
  printf("total number of particles %d \n", ng+nd+ns);
  z = (1-parts.h.time)/parts.h.time; 
  printf("redshift %f \n", z); 

  /* determine the start and end ID formed between this step and next step, note that the formation will occur in the future */ 
  for (i=0;i<=nfiles-1;i++){
    if(ns == ns_hist[i]){
      nstart=ns_hist[i]; 
      nend = ns_hist[i+1]-1; 
    }
  }
  if (nend < nstart){
    printf("no new SF"); 
    fprintf(stderr, "no new SF"); 
    return 0; 
  }
  ns_step = nend-nstart+1; 
  assert(ns_step >0); 
  printf("number of stars formed in this step %d \n", ns_step); 
  fprintf(stderr, "number of stars formed in this step %d \n", ns_step); 

  // this is the file from step 2;
  printf("                 1.2 \n");
  fpord = fopen(argv[2], "r");//reads in outputs .gord_eu_idx1
  if(fpord == NULL){
    printf("gord file %d doesn't exist \n", i); 
    fprintf(stderr, "gord file %d doesn't exist \n", i); 
      return 0; 
  }
  fread(&norigin,sizeof(int),1, fpord);
  printf("norigin = %d \n", norigin); 
  fprintf(stderr, "norigin = %d \n", norigin); 
  gord = (struct gasout *)malloc(norigin*sizeof(struct gasout));
  nread=fread(gord, sizeof(struct gasout), norigin, fpord); 
  assert(nread == norigin);
  fclose(fpord);

  // this is the file from step 3; 
  //changed all nuniq to count
  printf("                 1.3 \n");
  fpeu = fopen(argv[3], "r"); //reads in outputs .gord_euonly_cumu_idx1
  assert(fpeu != NULL);
  fread(&count,sizeof(int),1, fpeu);
  printf("number of particles in the eu-only file: %d \n", count);
  fprintf(stderr, "number of particles in the eu-only file: %d \n", count);
  geu = (struct gasout *)malloc(count*sizeof(struct gasout));
  nread = fread(geu, sizeof(struct gasout), count, fpeu); 
  assert(nread == count);
  fclose(fpeu);

  ////
  //printf("count                           %d \n",count);
  //for (i=0;i<=count-1;i++){
  //if(geu[i].iord > 0){
  //printf("geu[i].iord, geu[i].eu_cumu, geu[i].snap, geu[i].dwarfrac, gord[geu[i].iord].dwarfrac \n");
  //printf("%d %g %d %g \n",geu[i].iord, geu[i].eu_cumu, geu[i].snap, geu[i].dwarfrac);
  //printf("geu[i].snap                %d \n", geu[i].snap);
  //printf("geu[i].dwarfrac            %g \n", geu[i].dwarfrac);
  //printf("gord[geu[i].iord].dwarfrac %g \n", gord[geu[i].iord].dwarfrac);
  //}
  //}
  ////

  for (i=0;i<=count-1;i++){
    //gord[geu[i].iord].eu_sum = geu[i].eu_sum; /* updated to cumulative Eu */ //same as eu in step3 without massloss uses tag:_sum
    gord[geu[i].iord].eu_cumu = geu[i].eu_cumu; /* updated to cumulative Eu */ //same as eu in step3 with massloss tag:_cumu
    gord[geu[i].iord].dwarfrac = geu[i].dwarfrac;
  }
  neu = 0;
  for (i=0;i<=norigin-1;i++){
    if(gord[i].eu_cumu >1e-37){
      neu++; /* updated to cumulative Eu */
      //printf("eu_cumu %g, dwarf frac %g, id %d, mass %g, fmassloss %g, grpid %d \n", gord[i].eu_cumu, gord[i].dwarfrac, i, gord[i].mass, gord[i].fmassloss, gord[i].grpid); 
    }
  }
  printf("total number of particles with eumass non zero, %d \n",neu);
  fprintf(stderr, "total number of particles with eumass non zero, %d \n",neu);

  //return; 
  /* .iord file still necessary for smoothing */
  // this is the file with extension of .iord in the data directory
  printf("                 1.4 \n");
  fpiord = fopen(argv[4], "rb");//reads in iord files from simulations
  assert(fpiord != NULL); 
  xdrstdio_create(&xdrs, fpiord, XDR_DECODE); 
  xdr_int(&xdrs, &ntot); 
  printf("total particles in the iord file: %d \n", ntot);
  fprintf(stderr, "total particles in the iord file: %d \n", ntot);
  //return 0;
  iord = (int *)malloc(ntot*sizeof(int)); 
  for (i=0; i<=ntot-1; i++){
    xdr_int(&xdrs, &iord[i]); 
  }
  fclose(fpiord);
  /* for(i=0; i<=20; i++){
    printf("first 20 data \n");
    printf("%d, %g, %g \n", gout[i].id, gout[i].eu_sum, gout[i].mass); 
  }
  for(i=norigin-1; i>=norigin-21; i--){
    printf("last 20 data \n");
    printf("%d, %g, %g \n", gout[i].id, gout[i].eu_sum, gout[i].mass); 
    } */

/*--------------------------------------------------------------------*/
/*                                BLOCK2                              */
  printf("                 2 \n");
  /* paint stars */ 
  nseu_out = (float *)malloc(ns_step*sizeof(float));
  nseu_sm = (float *)malloc(ns_step*sizeof(float));
  dwarf = (float *)malloc(ns_step*sizeof(float)); //added dwarf arry
  dwarfdirect = (float *)malloc(ns_step*sizeof(float)); //added dwarf direct arry
  for(i=0;i<=ns_step-1;i++){
    nseu_out[i]=0.0; 
    nseu_sm[i]=0.0;
    dwarf[i] = 0.0; 
    dwarfdirect[i]=0.0;
  }
  count2 = 0;
  for (i=nstart;i<=nend;i++){
    if(i/100*100 == i) {
      printf("i, nend, %d, %d \n", i, nend); 
      fprintf(stderr, "i, nend, %d, %d \n", i, nend); 
    }
    giord = slog[i].iordgas; 
    gid = gord[giord].id; 
    if(gord[giord].eu_cumu >1e-37){
      count2++; 
      if(gid == -1){
	printf("error: gas doesn't exist but expect to form a star %d\n", giord);
	fprintf(stderr, "error: gas doesn't exist but expect to form a star %d\n", giord);
	fclose(stderr); 
	return 1;
      }
      /*direct painting*/
      nseu_out[i-nstart]=gord[giord].eu_cumu/gord[giord].mass;
      dwarfdirect[i-nstart] = gord[giord].dwarfrac;
      //printf("eucumu %g, frac %g \n", nseu_out[i-nstart], dwarfdirect[i-nstart]);
      ////if(dwarfdirect[i-nstart] > 0){
      ////printf("dwarfrac       %g \n",dwarfdirect[i-nstart]);////
      ////}

      /* smoothed painting */
      if (1 == 1){
	xc =parts.gp[gid].pos[0]; 
	yc =parts.gp[gid].pos[1]; 
	zc =parts.gp[gid].pos[2]; 
	gas = (struct gassort *)malloc(ng*sizeof(struct gassort)); //added dwarf to gassort
	search_r = 4.0; 
	for(l=0;l<=6;l++){
	  count1 = 0;
	  k = 0;
	  for (j=0;j<=ng-1;j++){
	    dist = sqrt((parts.gp[j].pos[0]-xc)*(parts.gp[j].pos[0]-xc) + (parts.gp[j].pos[1]-yc)*(parts.gp[j].pos[1]-yc) + (parts.gp[j].pos[2]-zc)*(parts.gp[j].pos[2]-zc)); 
	    if(dist*kpcunit/(1+z) <=search_r){ 
	      count1++;
	      gas[k].dist=dist;
	      gas[k].id = j;    /* id of gas */ 
	      gas[k].sid = gid; /* ID of the star */
	      gas[k].eu = gord[iord[j]].eu_cumu;
	      gas[k].mass = parts.gp[j].mass;
	      gas[k].dwarfrac = gord[iord[j]].dwarfrac; /*ADDED TO TRACK THE DWARF FRACTION*/
	      k=k+1;
	    }
	  }
	  if (count1 >=nsm) break; 
	  search_r = search_r + 2.0; 
	}
	assert(count1 >=nsm &&"search-r too small, not enough neighbors"); 
	
	gas1 = (struct gassort *)malloc(count1*sizeof(struct gassort));
	for (j=0; j<=count1-1; j++){
	  gas1[j].dist = gas[j].dist; 
	  gas1[j].id = gas[j].id;
	  gas1[j].sid = gas[j].sid;
	  gas1[j].eu = gas[j].eu;
	  gas1[j].mass = gas[j].mass;
	  gas1[j].dwarfrac = gas[j].dwarfrac; /*ADDED TO TRACK THE DWARF FRACTION*/
	}
	/*for(j=0;j<=ng-1;j++){
	  dist = sqrt((parts.gp[j].pos[0]-xc)*(parts.gp[j].pos[0]-xc) + (parts.gp[j].pos[1]-yc)*(parts.gp[j].pos[1]-yc) + (parts.gp[j].pos[2]-zc)*(parts.gp[j].pos[2]-zc));
	  if(dist*kpcunit/(1+z) <= search_r){
	  gas[k].dist=dist;
	  gas[k].id = j;
	  gas[k].sid = gid; 
	  gas[k].eu_sum = gord[iord[j]].eu_sum;
	  gas[k].mass = parts.gp[j].mass;
	  k=k+1;
	  }
	  }*/
	qsort(gas1, count1, sizeof(struct gassort), compare_function); //sorted by distance
	mtot = 0.0; 
	meu_tot = 0.0;
	meu_totdwarf = 0.0;
	for(k=0;k<=nsm-1;k++){
	  mtot = mtot+gas1[k].mass; 
	  meu_tot = meu_tot + gas1[k].eu;
	  meu_totdwarf = meu_totdwarf + (gas1[k].dwarfrac * gas1[k].eu); //make sure .eu is eu_cumu /*ADDED TO TRACK THE DWARF FRACTION*/
	}
	nseu_sm[i-nstart]=meu_tot/mtot;
	dwarf[i-nstart] = meu_totdwarf/meu_tot; /*ADDED TO TRACK THE DWARF FRACTION*/
	////if(dwarf[i-nstart] > 0){
	////printf("dwarf        %g \n", dwarf[i-nstart]);
	////}
	free(gas);
	free(gas1);
	//printf("smoothed eucumu %g, frac %g \n", nseu_sm[i-nstart], dwarf[i-nstart]);
      }
    }
  }
  printf("number of contaminated particles %d \n", count2);
  fprintf(stderr, "number of contaminated particles %d \n", count2);
  
  printf("%d %d\n", nstart, nend);
  
  //for(k=nstart;k<=nend;k++){
  //printf("dwarfdirect   %g \n", dwarfdirect[k]);
  //printf("dwarf         %g \n", dwarf[k]);
  //}


/*--------------------------------------------------------------------*/
/*                                BLOCK3                              */
  printf("                 3 \n");
  // return 0;
  printf("writing to file \n");
  /* this is output file */
  fpout = fopen(argv[5], "w"); // output file for step4 .stellar_rp
  assert(fpout != NULL);
  fwrite(&nstart, sizeof(int), 1, fpout);
  fwrite(&nend, sizeof(int), 1, fpout);
  //write smoothed eu
  nwrite = fwrite(nseu_out, sizeof(float), ns_step, fpout); 
  assert(nwrite == ns_step);
  //write direct eu
  nwrite = fwrite(nseu_sm, sizeof(float), ns_step, fpout); 
  assert(nwrite == ns_step);
  //write smoothed dwarf fraction
  nwrite = fwrite(dwarf, sizeof(float), ns_step, fpout);
  assert(nwrite == ns_step);
  //write direct dwarf fraction
  nwrite = fwrite(dwarfdirect, sizeof(float), ns_step, fpout);
  assert(nwrite == ns_step);
  // for (i=0; i<=norigin-1; i++){
  // fwrite(gout[i], sizeof(struct gassmall), 1, fpout); 
  // fprintf(fpout, "%d %g %g \n",gout[i].id, gout[i].eu_sum, gout[i].mass); 
  // } 
  fclose(fpout);
  //fclose(fplog);

  //for(k=0;k<=ns_step-1;k++){
  //  if(nseu_sm[k] > 0){
  //printf("smooth %g, direct %g,  dwarf smoothed %g, dwarf direct %g \n", nseu_sm[k],nseu_out[k], dwarf[k], dwarfdirect[k]);
  // }else{
  //printf("no negative values\n");
  //  }
  //}
  
  printf("step4 finished");

  for (j=0;j<=count-1;j++){
    printf("%d %g %g %g\n", slog[j].iordstar, slog[j].x, slog[j].y, slog[j].z);
  }

  return 0;
}

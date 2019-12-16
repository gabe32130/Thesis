// gcc -o read_gord_eu_pos read_gord_eu_pos.c tipsyxdr.c Tipsy.c tipsmall_defs3.h -lm
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

int main(int argc, char *argv[]){
  FILE *fpord, *fpeu, *fpout,*fptipsy, *fpcool, *fpgrp, *fphis, *fpslog;
  struct Tipsy parts;
  struct gassmall *gord;
  struct gasout *gout;
  struct starlog *slog;
  int i,j,k,nord, nfiles=400, norigin, nslog, ng, nd, ntotgrp, ntot, ntoteu, nstart, nend, ns, n;
  char filename[200], filename1[200], filename2[200], filename3[200]; 
  int *gid, *gid_new, *uniq, nuniq, count, count1, nread,npeu,result,*id, *eu, nwrite;
  int *iord, *igrp, *ideu, *sid, *grpid;
  float *ctime, *masseu, *dist, *massloss, *fmassloss;
  float *eumass, *mass;
  int ngd, ndark;

//reads in starlog
  fpslog=fopen("/mn/stornext/u3/shens/scratch/rprocess/L90Mpc8000_hithres.starlog.bin", "r");
  assert(fpslog !=NULL);
  fread(&nslog,sizeof(int),1, fpslog);
  printf("total particles in the starlog file: %d \n", nslog);
  slog = (struct starlog *)malloc(nslog*sizeof(struct starlog));
  nread=fread(slog, sizeof(struct starlog), nslog, fpslog);
  assert(nread == nslog);
  fclose(fpslog);


//reads in step1 paint_rp.c
  fpeu = fopen(argv[1], "r");
  assert(fpeu != NULL);
  fscanf(fpeu, "%d \n", &ntoteu);
  printf("total particles in the r-process file: %d \n", ntoteu);
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





//input files here are the output files from step 2 from make_ord_sml.c
  fpord = fopen(argv[2], "r"); 
  assert(fpord != NULL);
  fscanf(fpord, "%d \n", norigin);
  
  id = (int *)malloc(ntoteu*sizeof(int));
  eu = (int *)malloc(ntoteu*sizeof(int));
  mass = (float *)malloc(ntoteu*sizeof(float));
  fmassloss = (float *)malloc(ntoteu*sizeof(float));
  grpid = (int *)malloc(ntoteu*sizeof(int));

  for (i=0; i<=ntoteu-1; i++){
    fscanf(fpord, "%d, %g, %g, %g, %d \n", &id[i], &eu[i], &mass[i], &fmassloss[i], &grpid[i]);
  }
  fclose(fpeu);

  for (i=ntoteu-1; i>=ntoteu-10; i--){
    printf("%d, %g, %g, %g, %d \n",&id[i], &eu[i],&mass[i], &fmassloss[i], &grpid[i]);
  }




// this is  output file
// want this to print id, eu, grpid, x, y, z
// but to txt file so to be used for plotting




  for (j=0; j<=norigin; j++){
    //printf("%g %d %g %g %g \n", gout[j].eu, gout[j].grpid, slog[j].x, slog[j].y, slog[j].z);
  }

  return 0; 
}

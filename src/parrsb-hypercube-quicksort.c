#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <stdlib.h>

GenmapScalar g_min,g_max,g_delta;

void parRSBHypercubeQuickSortLocalSort(GenmapHandle h,int field,buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    // sort locally according to Fiedler vector
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt, fiedler,
                TYPE_DOUBLE, buf0);
  } else if(GENMAP_GLOBALID) {
    // sort locally according to globalId
    sarray_sort(struct GenmapElement_private, elements, (GenmapUInt)lelt,
                globalId0, TYPE_LONG, buf0);
  }

}

void parRSBFiedlerMinMax(GenmapHandle h,GenmapComm c,GenmapScalar *min,GenmapScalar *max) {
  *min = 1; *max = -1;

  GenmapElements e = GenmapGetElements(h);
  GenmapInt i;
  for(i = 0; i < GenmapGetNLocalElements(h); i++) {
    if(e[i].fiedler < *min) {
      *min = e[i].fiedler;
    }
    if(e[i].fiedler > *max) {
      *max = e[i].fiedler;
    }
  }

  GenmapGop(c,min,1,GENMAP_SCALAR,GENMAP_MIN);
  GenmapGop(c,max,1,GENMAP_SCALAR,GENMAP_MAX);
}

void parRSBHypercubeQuickSortInitProbes(GenmapHandle h,GenmapComm c,int field) {
  GenmapElements elements=GenmapGetElements(h);
  GenmapInt lelt=GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  int nsplitters=3;
  int pos=size/2;

  parRSBFiedlerMinMax(h,c,&g_min,&g_max);

  if(field==GENMAP_FIEDLER){
    g_delta=(g_max-g_min)/size;
    h->histogram->probes[0]=g_min;
    h->histogram->probes[1]=g_min+pos*g_delta;
    h->histogram->probes[2]=g_max;
  }
}

void parRSBHypercubeQuickSortUpdateCounts(GenmapHandle h,int nsplitters,int field) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  for(int i=0; i<nsplitters; i++) {
    h->histogram->count[i]=0;
  }

  GenmapElements p,e;
  if(field == GENMAP_FIEDLER){
    for(p=elements,e=p+lelt; p!=e; p++){
      for(int i=0; i<nsplitters; i++){
        if(p->fiedler<h->histogram->probes[i]){
          h->histogram->count[i]++;
        }
      }
    }
  }
}

int parRSBHypercubeQuickSortReachedThreshold(GenmapHandle h,GenmapComm c,GenmapLong *count,
                                             GenmapInt threshold,int field)
{
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  GenmapLong lelgt=GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;
  GenmapInt nrem=lelgt-partition_size*size;
  int pos=size/2;

  GenmapInt converged=1;

  if(rank==0) {
    GenmapLong expected=pos*partition_size+((pos<nrem)?pos:nrem);
    if(abs(count[1]-expected)>threshold) {
      converged=0;
    }
  }

  GenmapBcast(c,&converged,1,GENMAP_INT);
  return converged;
}

void parRSBHistoSortUpdateSplitter(GenmapHandle h,GenmapComm c,int i,GenmapInt threshold) {
  int size=GenmapCommSize(c);

  GenmapLong lelgt=GenmapGetNGlobalElements(h);
  GenmapInt partition_size=lelgt/size;
  GenmapInt nrem=lelgt-partition_size*size;
  int pos=size/2;

  GenmapLong expected=pos*partition_size+((pos<nrem)?pos:nrem);

  int indx=3*i-2;

  GenmapLong *count=h->histogram->count;
  GenmapScalar *probes=h->histogram->probes;

  if(abs(count[indx]-expected)<threshold) {
    //splitter is achieved
    return;
  }
  if(count[indx]<expected) { // we are going to dump the smaller probe
    count [indx-1]=count[indx];
    probes[indx-1]=probes[indx];
    count [indx]=count[indx]+(count[indx+1]-count[indx])/2;
    probes[indx]=probes[indx]+(probes[indx+1]-probes[indx])/2;
  } else { // we are going to dump the larger probe
    count [indx+1]=count[indx];
    probes[indx+1]=probes[indx];
    count [indx]=count[indx-1]+(count[indx]-count[indx-1])/2;
    probes[indx]=probes[indx-1]+(probes[indx]-probes[indx-1])/2;
  }
}

void parRSBHypercubeQuickSortUpdateProbes(GenmapHandle h,GenmapComm c,
                                          GenmapLong *count,GenmapInt threshold,int field)
{
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int rank=GenmapCommRank(c);
  int size=GenmapCommSize(c);
  int nsplitters=3;

  GenmapLong lelgt = GenmapGetNGlobalElements(h);

  if(rank==0) {
    for(int i=0; i<nsplitters; i++) {
      h->histogram->count[i]=count[i];
    }

    if(field == GENMAP_FIEDLER) {
      parRSBHistoSortUpdateSplitter(h,c,1,threshold);
    }
  }
}

int parRSBHypercubeQuickSortSetProc(GenmapHandle h,GenmapComm c,int nsplitters,
                                    int field,buffer *buf0)
{
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);
  
  // do a scane of all the elements less than the probe;
  // interesting stuff
  GenmapLong out[2][2],buf[2][2],in[2];

  parRSBHypercubeQuickSortUpdateCounts(h,nsplitters,field);

  //calculate destination procs
  GenmapLong n1=h->histogram->count[1]; 
  GenmapLong n2=lelt-n1;
  in[0]=n1; in[1]=n2;
  comm_scan(out,&(c->gsComm),genmap_gs_long,gs_add,in,2,buf);

  GenmapLong start1=out[0][0];
  GenmapLong start2=out[0][1];
  GenmapLong nelg1=out[1][0];
  GenmapLong nelg2=out[1][1];

  //let's do the first half first
  int np1=size/2;
  GenmapInt pNel=nelg1/np1;
  GenmapInt nrem=nelg1-pNel*np1;
  int id=start1/pNel;

  GenmapLong up=(id+1)*pNel+((id<nrem)?(id+1):nrem);
  printf("rank=%d np1=%d nelg1=%lld start1=%lld n1=%lld up=%lld lelt=%d\n",rank,
      np1,nelg1,start1,n1,up,lelt);
  GenmapInt i=0;
  while(i<n1) {
    if((i+start1)<up){
      elements[i].proc=id;
      //printf("fiedler=%lf,id=%d\n",elements[i].fiedler,elements[i].proc);
      i++;
    }else{
      id++;
      up=(id+1)*pNel+((id<nrem)?(id+1):nrem);
    }
  }

  //now the second half
  int np2=size-np1;
  pNel=nelg2/np2;
  nrem=nelg2-pNel*np2;
  id=start2/pNel;

  up=(id+1)*pNel+((id<nrem)?(id+1):nrem);
  printf("rank=%d np2=%d nelg2=%lld start2=%lld n2=%lld up=%lld lelt=%d\n",rank,
      np2,nelg2,start2,n2,up,lelt);
  i=0;
  while(i<n2) {
    if((i+start2)<up){
      elements[n1+i].proc=np1+id;
      //printf("fiedler=%lf,id=%d\n",elements[n1+i].fiedler,elements[n1+i].proc);
      i++;
    }else{
      id++;
      up=(id+1)*pNel+((id<nrem)?(id+1):nrem);
    }
  }

  return 0;
}

void parRSBHypercubeQuickSortTransferToProc(GenmapHandle h, int field, buffer *buf0) {
  GenmapElements elements = GenmapGetElements(h);
  GenmapInt lelt = GenmapGetNLocalElements(h);

  if(field == GENMAP_FIEDLER) { // Fiedler
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    GenmapScan(h, GenmapGetLocalComm(h));
  } else if(field == GENMAP_GLOBALID) {
    sarray_transfer(struct GenmapElement_private, &(h->elementArray), proc, 0,
                    &(h->cr));
    GenmapScan(h, GenmapGetLocalComm(h));
  }
}

void parRSBHypercubeQuickSort(GenmapHandle h,GenmapComm c,int field,buffer *buf0,int level) {
  GenmapScan(h,c);
  int size=GenmapCommSize(c);
  int rank=GenmapCommRank(c);

  // 10% of load balanced partition size 
  GenmapInt threshold=(GenmapGetNGlobalElements(h)/(10*size));
  if(threshold<2) threshold=1;

  // sort locally.
  parRSBHypercubeQuickSortLocalSort(h,field,buf0);

  // We are done if size==1
  if(size==1) return;

  // Else we continue
  int nsplitters=3;

  // Allocate space for temp counts
  GenmapLong *count=NULL;
  // Allocate space for probes and counts
  if(level==0) {
    GenmapMalloc(nsplitters,&h->histogram->probes);
    GenmapMalloc(nsplitters,&h->histogram->count);
  }
  if(rank==0) {
    GenmapMalloc(nsplitters,&count);
  }

  // init probes values
  parRSBHypercubeQuickSortInitProbes(h,c,field);
  if(rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("level: %d iter: 0 probe[%d]= " GenmapScalarFormat "\n",level,i,
          h->histogram->probes[i]);
    }

  // update counts locally 
  parRSBHypercubeQuickSortUpdateCounts(h,nsplitters,field);

  // global reduction
  GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
  if(rank==0)
    for(int i=0; i<nsplitters; i++) {
      printf("level: %d iter: 0 count[%d]= " GenmapLongFormat "\n",level,i,count[i]);
    }

  int iter=0;
  while(!parRSBHypercubeQuickSortReachedThreshold(h,c,count,threshold,field)){
    parRSBHypercubeQuickSortUpdateProbes(h,c,count,threshold,field);
#if 0
    if(rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("level: %d iter: %d probe[%d]= " GenmapScalarFormat "\n",level,iter+1,
            i,h->histogram->probes[i]);
      }
    if(rank==0) printf("Update probes: done\n");
#endif

    GenmapBcast(c,h->histogram->probes,nsplitters,GENMAP_LONG);
#if 0
    if(rank==1)
      for(int i=0; i<nsplitters; i++){
        printf("%d: %d probe[%d]= " GenmapScalarFormat "\n",rank,iter,i,h->histogram->probes[i]);
      }
    if(rank==0) printf("Bcast: done\n");
#endif

    parRSBHypercubeQuickSortUpdateCounts(h,nsplitters,field);

    // global reduction
    GenmapReduce(c,count,h->histogram->count,nsplitters,GENMAP_LONG,GENMAP_SUM);
#if 0
    if(rank==0)
      for(int i=0; i<nsplitters; i++){
        printf("level: %d iter: %d count[%d]= " GenmapLongFormat "\n",level,iter+1,
            i,count[i]);
      }
#endif
    iter++;
  }

  // set destination processor id for each element
  parRSBHypercubeQuickSortSetProc(h,c,nsplitters,field,buf0);

  // send elements to right processor
  GenmapCrystalInit(h,c);
  parRSBHypercubeQuickSortTransferToProc(h,field,buf0);
  GenmapCrystalFinalize(h);

  //split comms here
  int bin;
  if(rank<size/2) bin=0;
  else bin=1;

  GenmapCommExternal local;
  MPI_Comm_split(c->gsComm.c,bin,rank,&local);
  GenmapComm new; GenmapCreateComm(&new,local);
  MPI_Comm_free(&local);

  parRSBHypercubeQuickSort(h,new,field,buf0,1);

  // sort locally.
  parRSBHypercubeQuickSortLocalSort(h,field,buf0);

  if(level==0) {
    GenmapFree(h->histogram->probes);
    GenmapFree(h->histogram->count);
  }
  if(rank==0) {
    GenmapFree(count);
  }
}

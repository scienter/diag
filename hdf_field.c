// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "math.h"
#include "mpi.h"
#include "string.h"
#define FIELD 0
#define DENSITY 1
#define FIELDNDFX 2
#define pi 3.14159265359

void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt);
void saveIntMeta(char *fileName,char *dataName,int *data);
void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz);
void saveFile(double ***EzR,double ***EzI,double ***ErR,double ***ErI,double ***EpR,double ***EpI,double ***BzR,double ***BzI,double ***BrR,double ***BrI,double ***BpR,double ***BpI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode,char *name,int stX,int stY,double minY,double maxY);
void saveEnergy(double *energy,double ***ErR,double ***ErI,double ***EpR,double ***EpI,int angleDivision,int nx,int ny,int step,int stX,int stY,int cellNum);
void saveRhoFile(double ***fR,double ***fI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode,char *name,int stX,int stY,double minY,double maxY,int shift);
void saveRhoSumFile(double ***dataLow,double ***dataUp,double *dataX,double *dataY,double angle,int nx,int ny,int stX, int stY,char *name,double minY,double maxY);
void saveRho(double ***dataLow,double ***dataUp,double ***fR,double ***fI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode);
void deleteField(double ***field,int nx,int ny,int nz);
double ***memoryAsign(int nx, int ny, int nz);
int whatFunctionMode(char *str);
void sumRho(double ***sumLow,double ***sumUp,double ***dataLow,double ***dataUp,int nx,int ny);

void main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],fileType[100],**name;
   int nx,ny,mode,numMode,stX,stY,initial,final,timeStep,step,dataType;
   int i,j,s,numS,angleDivision,index,cellNum,shift;
   double angle,minY,maxY;
   double ***EzR,***EzI,***ErR,***ErI,***EpR,***EpI,*dataX,*dataY;
   double ***BzR,***BzI,***BrR,***BrI,***BpR,***BpI,***fR,***fI;
   double ***dataLow,***dataUp,***sumLow,***sumUp,*energy;

    int myrank, nTasks;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   if(argc < 6)
   {
      printf("hdf_field mode initial final timestep stX stY\n");
      printf("mode(0): fileType angle minY maxY shift\n");
      printf("mode(1): fileType angleDivision cellNum\n");
      printf("mode(2): minY maxY angle numS species1 species2 ... \n");
      printf("mode(3 : summing): minY maxY angle numS species1 species2 ... \n");
      exit(0);
   }
   mode=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);
   stX=atoi(argv[5]);
   stY=atoi(argv[6]);

   switch (mode) {
   case 3 :
     minY=atof(argv[7]);
     maxY=atof(argv[8]);
     angle=atof(argv[9]);
     numS=atoi(argv[10]);

     name=(char **)malloc(numS*sizeof(char *));
     for(s=0; s<numS; s++)
       name[s]=(char *)malloc(100*sizeof(char ));

     for(s=0; s<numS; s++)
       sprintf(name[s],"%sdensity",argv[s+11]);

     sprintf(fileName,"%s%d.h5",name[0],initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"numMode",&numMode);
    
     dataX=(double *)malloc(nx*sizeof(double));
     dataY=(double *)malloc(ny*sizeof(double)); 
     fR=memoryAsign(nx,ny,numMode);
     fI=memoryAsign(nx,ny,numMode);
     dataLow=memoryAsign(nx,ny,1);
     dataUp=memoryAsign(nx,ny,1);
     sumLow=memoryAsign(nx,ny,1);
     sumUp=memoryAsign(nx,ny,1);

     for(step=initial; step<=final; step+=timeStep)
     {
       for(i=0; i<nx; i++)
         for(j=0; j<ny; j++) {
           sumLow[i][j][0]=0.0;
           sumUp[i][j][0]=0.0;
         }

       for(s=0; s<numS; s++)
       {
         for(i=0; i<nx; i++)
           for(j=0; j<ny; j++) {
             dataLow[i][j][0]=0.0;
             dataUp[i][j][0]=0.0;
           }
         sprintf(fileName,"%s%d.h5",name[s],step);
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;

         restoreFloatArray(fileName,"X",dataX,nx);
         restoreFloatArray(fileName,"Y",dataY,ny);
         restoreFieldComp(fR,fileName,"R",nx,ny,numMode);
         restoreFieldComp(fI,fileName,"I",nx,ny,numMode);
         printf("%s is done\n",fileName);

         sprintf(fileName,"%s%d","denSum",step);
         saveRho(dataLow,dataUp,fR,fI,dataX,dataY,angle,nx,ny,numMode);      
         sumRho(sumLow,sumUp,dataLow,dataUp,nx,ny);      
       }
       sprintf(fileName,"sumDensity%d",step);
       saveRhoSumFile(sumLow,sumUp,dataX,dataY,angle,nx,ny,stX,stY,fileName,minY,maxY);       
       printf("%s_%g is made.\n",fileName,angle);
     }

     deleteField(fR,nx,ny,numMode);
     deleteField(fI,nx,ny,numMode);
     deleteField(dataLow,nx,ny,1);
     deleteField(dataUp,nx,ny,1);
     deleteField(sumLow,nx,ny,1);
     deleteField(sumUp,nx,ny,1);
     free(dataX); free(dataY);
     for(s=0; s<numS; s++) free(name[s]); free(name);

     break;
   case 2 :
     minY=atof(argv[7]);
     maxY=atof(argv[8]);
     angle=atof(argv[9]);
     numS=atoi(argv[10]);

     name=(char **)malloc(numS*sizeof(char *));
     for(s=0; s<numS; s++)
       name[s]=(char *)malloc(100*sizeof(char ));

     for(s=0; s<numS; s++)
       sprintf(name[s],"%sdensity",argv[s+11]);

     sprintf(fileName,"%s%d.h5",name[0],initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;
     saveIntMeta(fileName,"nx",&nx);
     saveIntMeta(fileName,"ny",&ny);
     saveIntMeta(fileName,"numMode",&numMode);
    
     dataX=(double *)malloc(nx*sizeof(double));
     dataY=(double *)malloc(ny*sizeof(double)); 
     fR=memoryAsign(nx,ny,numMode);
     fI=memoryAsign(nx,ny,numMode);
     dataLow=memoryAsign(nx,ny,1);
     dataUp=memoryAsign(nx,ny,1);

     for(step=initial; step<=final; step+=timeStep)
     {
       for(i=0; i<nx; i++)
         for(j=0; j<ny; j++) {
           dataLow[i][j][0]=0.0;
           dataUp[i][j][0]=0.0;
         }

       for(s=0; s<numS; s++)  {
         sprintf(fileName,"%s%d.h5",name[s],step);
         if(fopen(fileName,"r")==NULL)  {
           printf("%s is not exited.\n",fileName);
           exit(0);
         } else ;

         restoreFloatArray(fileName,"X",dataX,nx);
         restoreFloatArray(fileName,"Y",dataY,ny);
         restoreFieldComp(fR,fileName,"R",nx,ny,numMode);
         restoreFieldComp(fI,fileName,"I",nx,ny,numMode);
         printf("%s is done\n",fileName);

         sprintf(fileName,"%s%d","denSum",step);
         saveRho(dataLow,dataUp,fR,fI,dataX,dataY,angle,nx,ny,numMode);       
         saveRhoSumFile(dataLow,dataUp,dataX,dataY,angle,nx,ny,stX,stY,fileName,minY,maxY);       
       }
       printf("%s is made.\n",fileName);
     }

     deleteField(fR,nx,ny,numMode);
     deleteField(fI,nx,ny,numMode);
     deleteField(dataLow,nx,ny,1);
     deleteField(dataUp,nx,ny,1);
     free(dataX); free(dataY);
     for(s=0; s<numS; s++) free(name[s]); free(name);

     break;
   case 1 :
     angleDivision=atoi(argv[8]);
     cellNum=atoi(argv[9]);
     sprintf(fileType,"%s",argv[7]);     
     dataType=whatFunctionMode(fileType);
     sprintf(fileName,"%s%d.h5",fileType,initial);
     if(fopen(fileName,"r")==NULL)  {
       printf("%s is not exited.\n",fileName);
       exit(0);
     } else ;

     step=(final-initial)/timeStep;
     energy=(double *)malloc((step+1)*sizeof(double)); 

     for(step=initial; step<=final; step+=timeStep)
     {
       sprintf(fileName,"%s%d.h5",fileType,step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;
       saveIntMeta(fileName,"nx",&nx);
       saveIntMeta(fileName,"ny",&ny);
       saveIntMeta(fileName,"numMode",&numMode);
    
       dataX=(double *)malloc(nx*sizeof(double));
       dataY=(double *)malloc(ny*sizeof(double)); 
    
       ErR=memoryAsign(nx,ny,numMode);
       ErI=memoryAsign(nx,ny,numMode);
       EpR=memoryAsign(nx,ny,numMode);
       EpI=memoryAsign(nx,ny,numMode);

       restoreFloatArray(fileName,"X",dataX,nx);
       restoreFloatArray(fileName,"Y",dataY,ny);
       restoreFieldComp(ErR,fileName,"ErR",nx,ny,numMode);
       restoreFieldComp(ErI,fileName,"ErI",nx,ny,numMode);
       restoreFieldComp(EpR,fileName,"EpR",nx,ny,numMode);
       restoreFieldComp(EpI,fileName,"EpI",nx,ny,numMode);
       printf("%s is done\n",fileName);

       sprintf(fileName,"%s%d",fileType,step);
       index=(step-initial)/timeStep;
       saveEnergy(&energy[index],ErR,ErI,EpR,EpI,angleDivision,nx,ny,step,stX,stY,cellNum);

       free(dataX); free(dataY);
     }

     out=fopen("laserEnergy","w");
     for(step=initial; step<final; step+=timeStep) {
       index=(step-initial)/timeStep;
       fprintf(out,"%d %g %g\n",step,energy[index],energy[index]-energy[index+1]);
     }
     fclose(out);
     printf("laserEnergy is made\n");

     deleteField(ErR,nx,ny,numMode);
     deleteField(ErI,nx,ny,numMode);
     deleteField(EpR,nx,ny,numMode);
     deleteField(EpI,nx,ny,numMode);
     free(energy); 

     break;
   case 0 :
     angle=atof(argv[8]);
     minY=atof(argv[9]);
     maxY=atof(argv[10]);
     shift=atoi(argv[11]);
     sprintf(fileType,"%s",argv[7]);     
     dataType=whatFunctionMode(fileType);

     for(step=initial; step<=final; step+=timeStep)
     {
       sprintf(fileName,"%s%d.h5",fileType,step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;
       saveIntMeta(fileName,"nx",&nx);
       saveIntMeta(fileName,"ny",&ny);
       saveIntMeta(fileName,"numMode",&numMode);
       dataX=(double *)malloc(nx*sizeof(double));
       dataY=(double *)malloc(ny*sizeof(double)); 
       if(dataType==FIELD || dataType==FIELDNDFX) {
         EzR=memoryAsign(nx,ny,numMode);
         EzI=memoryAsign(nx,ny,numMode);
         ErR=memoryAsign(nx,ny,numMode);
         ErI=memoryAsign(nx,ny,numMode);
         EpR=memoryAsign(nx,ny,numMode);
         EpI=memoryAsign(nx,ny,numMode);
         BzR=memoryAsign(nx,ny,numMode);
         BzI=memoryAsign(nx,ny,numMode);
         BrR=memoryAsign(nx,ny,numMode);
         BrI=memoryAsign(nx,ny,numMode);
         BpR=memoryAsign(nx,ny,numMode);
         BpI=memoryAsign(nx,ny,numMode);
       } else {
         fR=memoryAsign(nx,ny,numMode);
         fI=memoryAsign(nx,ny,numMode);
       }

       restoreFloatArray(fileName,"X",dataX,nx);
       restoreFloatArray(fileName,"Y",dataY,ny);
       if(dataType==DENSITY) {
         restoreFieldComp(fR,fileName,"R",nx,ny,numMode);
         restoreFieldComp(fI,fileName,"I",nx,ny,numMode);
         printf("%s is done\n",fileName);

         sprintf(fileName,"%s%d",fileType,step);
         saveRhoFile(fR,fI,dataX,dataY,angle,nx,ny,numMode,fileName,stX,stY,minY,maxY,shift);
       } else if(dataType==FIELD) {
         restoreFieldComp(EzR,fileName,"EzR",nx,ny,numMode);
         restoreFieldComp(EzI,fileName,"EzI",nx,ny,numMode);
         restoreFieldComp(ErR,fileName,"ErR",nx,ny,numMode);
         restoreFieldComp(ErI,fileName,"ErI",nx,ny,numMode);
         restoreFieldComp(EpR,fileName,"EpR",nx,ny,numMode);
         restoreFieldComp(EpI,fileName,"EpI",nx,ny,numMode);
         restoreFieldComp(BzR,fileName,"BzR",nx,ny,numMode);
         restoreFieldComp(BzI,fileName,"BzI",nx,ny,numMode);
         restoreFieldComp(BrR,fileName,"BrR",nx,ny,numMode);
         restoreFieldComp(BrI,fileName,"BrI",nx,ny,numMode);
//         restoreFieldComp(BpR,fileName,"BpR",nx,ny,numMode);
//         restoreFieldComp(BpI,fileName,"BpI",nx,ny,numMode);
         restoreFieldComp(BpR,fileName,"FR",nx,ny,numMode);
         restoreFieldComp(BpI,fileName,"FI",nx,ny,numMode);
         printf("%s is done\n",fileName);

         sprintf(fileName,"%s%d",fileType,step);
         saveFile(EzR,EzI,ErR,ErI,EpR,EpI,BzR,BzI,BrR,BrI,BpR,BpI,dataX,dataY,angle,nx,ny,numMode,fileName,stX,stY,minY,maxY);       
       } else if(dataType==FIELDNDFX) {
         restoreFieldComp(EzR,fileName,"EzR",nx,ny,numMode);
         restoreFieldComp(EzI,fileName,"EzI",nx,ny,numMode);
         restoreFieldComp(ErR,fileName,"PrR",nx,ny,numMode);
         restoreFieldComp(ErI,fileName,"PrI",nx,ny,numMode);
         restoreFieldComp(EpR,fileName,"PlR",nx,ny,numMode);
         restoreFieldComp(EpI,fileName,"PlI",nx,ny,numMode);
         restoreFieldComp(BzR,fileName,"BzR",nx,ny,numMode);
         restoreFieldComp(BzI,fileName,"BzI",nx,ny,numMode);
         restoreFieldComp(BrR,fileName,"SrR",nx,ny,numMode);
         restoreFieldComp(BrI,fileName,"SrI",nx,ny,numMode);
//         restoreFieldComp(BpR,fileName,"SlR",nx,ny,numMode);
//         restoreFieldComp(BpI,fileName,"SlI",nx,ny,numMode);
         restoreFieldComp(BpR,fileName,"FR",nx,ny,numMode);
         restoreFieldComp(BpI,fileName,"FI",nx,ny,numMode);
         printf("%s is done\n",fileName);

         sprintf(fileName,"%s%d",fileType,step);
         saveFile(EzR,EzI,ErR,ErI,EpR,EpI,BzR,BzI,BrR,BrI,BpR,BpI,dataX,dataY,angle,nx,ny,numMode,fileName,stX,stY,minY,maxY);       
       }

       if(dataType==DENSITY) {
         deleteField(fR,nx,ny,numMode);
         deleteField(fI,nx,ny,numMode);
       } else if(dataType==FIELD || dataType==FIELDNDFX) {
         deleteField(EzR,nx,ny,numMode);
         deleteField(EzI,nx,ny,numMode);
         deleteField(ErR,nx,ny,numMode);
         deleteField(ErI,nx,ny,numMode);
         deleteField(EpR,nx,ny,numMode);
         deleteField(EpI,nx,ny,numMode);
         deleteField(BzR,nx,ny,numMode);
         deleteField(BzI,nx,ny,numMode);
         deleteField(BrR,nx,ny,numMode);
         deleteField(BrI,nx,ny,numMode);
         deleteField(BpR,nx,ny,numMode);
         deleteField(BpI,nx,ny,numMode);
       }
       free(dataX); free(dataY);  

     }


     break;
   }
    MPI_Finalize();
}

void saveEnergy(double *energy,double ***ErR,double ***ErI,double ***EpR,double ***EpI,int angleDivision,int nx,int ny,int step,int stX,int stY,int cellNum)
{
  int i,j,m,n;
  double Er,Ep,sinP,cosP,reangle,delPhi,sum;
  char fileName[100];
  FILE *out;
  
  delPhi=2*3.14159265359/(double)angleDivision;
//  out=fopen("laserEnergy","a");
  m=1;

  sum=0;
  for(n=0; n<angleDivision; n++) {
    reangle=n*delPhi;
    cosP=cos(reangle);
    sinP=sin(reangle);
    for(i=nx-cellNum; i<nx; i+=stX) 
      for(j=0; j<ny; j+=stY) {        
        Er=ErR[i][j][m]*cosP-ErI[i][j][m]*sinP;
        Ep=EpR[i][j][m]*cosP-EpI[i][j][m]*sinP;
        sum+=(Er*Er+Ep*Ep)*j*delPhi*stX*stY;
      }
  }
  *energy=sum;
//  fprintf(out,"%d %g\n",step,sum);
//  fclose(out);
//  printf("step%d is done.\n",step);
}

//lala
void saveRhoSumFile(double ***dataLow,double ***dataUp,double *dataX,double *dataY,double angle,int nx,int ny,int stX, int stY,char *name,double minY,double maxY)
{
  int i,j,m;
  double f,sinP,cosP,sinP1,cosP1,reangle;
  char fileName[100];
  FILE *out;
  
  sprintf(fileName,"%s_%g",name,angle);
  out=fopen(fileName,"w");
  for(i=0; i<nx; i+=stX)  {
    for(j=ny-1; j>=0; j-=stY) {        
      if(-dataY[j]>minY && -dataY[j]<maxY)
        fprintf(out,"%g %g %g\n",dataX[i],-dataY[j],dataLow[i][j][0]);
      else ;
    }
    for(j=1; j<ny; j+=stY) {        
      if(dataY[j]>minY && dataY[j]<maxY)
        fprintf(out,"%g %g %g\n",dataX[i],dataY[j],dataUp[i][j][0]);
      else;
    }
    fprintf(out,"\n");
  }
  fclose(out);
}

void sumRho(double ***sumLow,double ***sumUp,double ***dataLow,double ***dataUp,int nx,int ny)
{
  int i,j;
  
  for(i=0; i<nx; i++)  {
    for(j=ny-1; j>=0; j--)       
      sumLow[i][j][0]+=dataLow[i][j][0];
    for(j=1; j<ny; j++)         
      sumUp[i][j][0]+=dataUp[i][j][0];
  }
}


void saveRho(double ***dataLow,double ***dataUp,double ***fR,double ***fI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode)
{
  int i,j,m;
  double f,sinP,cosP,sinP1,cosP1,reangle;
  char fileName[100];
  FILE *out;
  
  reangle=angle*pi/180.0;
  for(i=0; i<nx; i++)  {
    for(j=ny-1; j>=0; j--) {        
      for(m=0; m<numMode; m++) {
        cosP1=cos(m*(reangle+pi));
        sinP1=sin(m*(reangle+pi));
        dataLow[i][j][0]+=fR[i][j][m]*cosP1-fI[i][j][m]*sinP1;
      }
    }
    for(j=1; j<ny; j++) {        
      f=0.0;
      for(m=0; m<numMode; m++) {
        cosP=cos(m*reangle);
        sinP=sin(m*reangle);
        dataUp[i][j][0]+=fR[i][j][m]*cosP-fI[i][j][m]*sinP;
      }
    }
  }
}

void saveRhoFile(double ***fR,double ***fI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode,char *name,int stX,int stY,double minY,double maxY,int shift)
{
  int i,j,m;
  double f,sinP,cosP,sinP1,cosP1,reangle,dY;
  char fileName[100];
  FILE *out;
  
  reangle=angle*pi/180.0;
//  dY=dataY[1]-dataY[0];
/*
  for(m=0; m<numMode; m++) {
    cosP=cos(m*reangle);
    sinP=sin(m*reangle);
    cosP1=cos(m*(reangle+pi));
    sinP1=sin(m*(reangle+pi));
    sprintf(fileName,"%s_%d_%g",name,m,angle);
    out=fopen(fileName,"w");
    for(i=0; i<nx; i+=stX)  {
      for(j=ny-1; j>0; j-=stY) {        
        f=fR[i][j][m]*cosP1-fI[i][j][m]*sinP1;
        fprintf(out,"%g %g %g\n",dataX[i],-dataY[j],fabs(f));
      }

      f=fR[i][0][0];
      fprintf(out,"%g %g %g\n",dataX[i],dataY[0],2.0*fabs(f));

      for(j=1; j<ny; j+=stY) {        
        f=fR[i][j][m]*cosP-fI[i][j][m]*sinP;
        fprintf(out,"%g %g %g\n",dataX[i],dataY[j],fabs(f));
      }
      fprintf(out,"\n");
    }
    fclose(out);
    printf("%s is made.\n",fileName);
  }
*/
  //summation
  sprintf(fileName,"%s_sum_%g",name,angle);
  out=fopen(fileName,"w");
  for(i=0; i<nx; i+=stX)  {
    for(j=ny-1; j>=0; j-=stY) {        
      f=0.0;
      for(m=0; m<numMode; m++) {
        cosP1=cos(m*(reangle+pi));
        sinP1=sin(m*(reangle+pi));
        f+=fR[i][j][m]*cosP1-fI[i][j][m]*sinP1;
      }
      if(-dataY[j]>minY && -dataY[j]<maxY)
        fprintf(out,"%g %g %g\n",dataX[i],-dataY[j],fabs(f));
      else ;
    }
//    if(shift==0) {
//      f=fR[i][0][0];
//      fprintf(out,"%g %g %g\n",dataX[i],dataY[0],2.0*fabs(f));
//    } else ;

    for(j=shift; j<ny-shift; j+=stY) {        
      f=0.0;
      for(m=0; m<numMode; m++) {
        cosP=cos(m*reangle);
        sinP=sin(m*reangle);
        f+=fR[i][j+shift][m]*cosP-fI[i][j+shift][m]*sinP;
      }
      if(dataY[j]>minY && dataY[j]<maxY)
        fprintf(out,"%g %g %g\n",dataX[i],dataY[j],fabs(f));
      else ;
    }
    fprintf(out,"\n");
  }
  fclose(out);
  printf("%s is made.\n",fileName);
}

void saveFile(double ***EzR,double ***EzI,double ***ErR,double ***ErI,double ***EpR,double ***EpI,double ***BzR,double ***BzI,double ***BrR,double ***BrI,double ***BpR,double ***BpI,double *dataX,double *dataY,double angle,int nx,int ny,int numMode,char *name,int stX,int stY,double minY,double maxY)
{
  int i,j,m;
  double Ez,Er,Ep,Bz,Br,Bp,sinP,cosP,sinP1,cosP1,reangle;
  char fileName[100];
  FILE *out;
 
  reangle=angle*pi/180.0;
/*
  for(m=0; m<numMode; m++) {
    cosP=cos(m*reangle);
    sinP=sin(m*reangle);
    cosP1=cos(m*(reangle+pi));
    sinP1=sin(m*(reangle+pi));
    sprintf(fileName,"%s_%d_%g",name,m,angle);
    out=fopen(fileName,"w");
    for(i=0; i<nx; i+=stX)  {
      for(j=ny-1; j>0; j-=stY) {        
        Ez=EzR[i][j][m]*cosP1-EzI[i][j][m]*sinP1;
        Er=ErR[i][j][m]*cosP1-ErI[i][j][m]*sinP1;
        Ep=EpR[i][j][m]*cosP1-EpI[i][j][m]*sinP1;
        Bz=BzR[i][j][m]*cosP1-BzI[i][j][m]*sinP1;
        Br=BrR[i][j][m]*cosP1-BrI[i][j][m]*sinP1;
        Bp=BpR[i][j][m]*cosP1-BpI[i][j][m]*sinP1;
        fprintf(out,"%g %g %g %g %g %g %g %g\n",dataX[i],-dataY[j],Ez,Er,Ep,Bz,Br,Bp);
      }
      for(j=0; j<ny; j+=stY) {        
        Ez=EzR[i][j][m]*cosP-EzI[i][j][m]*sinP;
        Er=ErR[i][j][m]*cosP-ErI[i][j][m]*sinP;
        Ep=EpR[i][j][m]*cosP-EpI[i][j][m]*sinP;
        Bz=BzR[i][j][m]*cosP-BzI[i][j][m]*sinP;
        Br=BrR[i][j][m]*cosP-BrI[i][j][m]*sinP;
        Bp=BpR[i][j][m]*cosP-BpI[i][j][m]*sinP;
        fprintf(out,"%g %g %g %g %g %g %g %g\n",dataX[i],dataY[j],Ez,Er,Ep,Bz,Br,Bp);
      }
      fprintf(out,"\n");
    }
    fclose(out);
    printf("%s is made.\n",fileName);
  }
*/
  //summation
  sprintf(fileName,"%s_sum_%g",name,angle);
  out=fopen(fileName,"w");
  for(i=0; i<nx; i+=stX)  {
    for(j=ny-1; j>=1; j-=stY) {        
      Ez=Er=Ep=Bz=Br=Bp=0.0;
      for(m=0; m<numMode; m++) {
        cosP1=cos(m*(reangle+pi));
        sinP1=sin(m*(reangle+pi));
        Ez+=EzR[i][j][m]*cosP1-EzI[i][j][m]*sinP1;
        Er+=ErR[i][j][m]*cosP1-ErI[i][j][m]*sinP1;
        Ep+=EpR[i][j][m]*cosP1-EpI[i][j][m]*sinP1;
        Bz+=BzR[i][j][m]*cosP1-BzI[i][j][m]*sinP1;
        Br+=BrR[i][j][m]*cosP1-BrI[i][j][m]*sinP1;
        Bp+=BpR[i][j][m]*cosP1-BpI[i][j][m]*sinP1;
      }
      if(-dataY[j]>minY && -dataY[j]<maxY)
        fprintf(out,"%g %g %g %g %g %g %g %g\n",dataX[i],-dataY[j],Ez,Er,Ep,Bz,Br,Bp);
    }
    for(j=0; j<ny; j+=stY) {        
      Ez=Er=Ep=Bz=Br=Bp=0.0;
      for(m=0; m<numMode; m++) {
        cosP=cos(m*reangle);
        sinP=sin(m*reangle);
        Ez+=EzR[i][j][m]*cosP-EzI[i][j][m]*sinP;
        Er+=ErR[i][j][m]*cosP-ErI[i][j][m]*sinP;
        Ep+=EpR[i][j][m]*cosP-EpI[i][j][m]*sinP;
        Bz+=BzR[i][j][m]*cosP-BzI[i][j][m]*sinP;
        Br+=BrR[i][j][m]*cosP-BrI[i][j][m]*sinP;
        Bp+=BpR[i][j][m]*cosP-BpI[i][j][m]*sinP;
      }
      if(dataY[j]>minY && dataY[j]<maxY)
        fprintf(out,"%g %g %g %g %g %g %g %g\n",dataX[i],dataY[j],Ez,Er,Ep,Bz,Br,Bp);
    }
    fprintf(out,"\n");
  }
  fclose(out);

  printf("%s is made.\n",fileName);
}

void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz)
{
  int i,j,k,start;
  double *field;
  char name[100];
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny;
  dimsf[1]=nx;
  dimsf[2]=nz;
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=ny;
  count[1]=nx;
  count[2]=nz;
  offset[0]=0;
  offset[1]=0;
  offset[2]=0;
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nx*ny*nz*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(j=0; j<ny; j++)
    for(i=0; i<nx; i++)
    {
      for(k=0; k<nz; k++)
        data[i][j][k]=field[start+k];
      start+=nz;
    }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreFloatArray(char *fileName,char *dataName,double *data,int totalCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=totalCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

  filespace=H5Screate_simple(1,metaDim,NULL);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}


void saveIntMeta(char *fileName,char *dataName,int *data)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=1;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         field[i][j][k]=0.0;

   return field;
}

void deleteField(double ***field,int nx,int ny,int nz)
{
   int i,j,k;
   for(i=0; i<nx; i++)
   {
     for(j=0; j<ny; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


int whatFunctionMode(char *str)
{
   int result=DENSITY;
   if(strstr(str,"field")) { result=FIELD; } else ;
   if(strstr(str,"fieldSplit")) { result=FIELDNDFX; } else;
   if(strstr(str,"0density")) { result=DENSITY; } else ;
   return result;
}


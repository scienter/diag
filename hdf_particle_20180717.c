#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "mpi.h"
#include "math.h"

void restoreData(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC,int columns);
void restore1Data(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC);
void calSub(int totalCnt,int *cntSub,int *start);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveParticleFile(double *data,double minPx,double minX,int cntSub,int step,int species,double minE,double maxE,int atomNum);
void saveDensityFile(double *data,double minX,double maxX,int numX,double minY,double maxY,int numY,double minZ,double maxZ,int step,int species,int cntSub);
double whatMass(int atomNum);

void main(int argc, char *argv[])
{
   int mode,species,initial,final,timeStep,atomNum;
   int i,n,step,totalCnt,cntSub,start,column,sharePNum,index;
   int *recv,numX,numY;
   double minPx,minE,maxE,minX,maxX,minY,maxY,minZ,maxZ;
   double slope,rangeX;
   double *data;
   FILE *out;
   char fileName[100],dataName[100],dataSet[100];
   int myrank, nTasks;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


   recv=(int *)malloc(nTasks*sizeof(int ));

   if(argc < 6)   {
     if(myrank==0) {
       printf("mpirun -np [ ] hdf_particle mode species initial final step\n");
       printf("mode(1) : minPx minX atomNum\n");
       printf("mode(2) : minE maxE atomNum\n");
       printf("mode(3) : slope rangeX numX minY maxY numY minZ maxZ\n");
     } else ;
     MPI_Finalize();  
     exit(0);
   }  else;

   mode=atoi(argv[1]);
   species=atoi(argv[2]);
   initial=atoi(argv[3]);
   final=atoi(argv[4]);
   timeStep=atoi(argv[5]);

   if(mode==1) {
     minPx=atof(argv[6]); minE=-1e9; maxE=1e9;
     minX=atof(argv[7]);
     atomNum=atoi(argv[8]);
   } else if(mode==2) {
     minPx=-1e5; minE=atof(argv[6]); maxE=atof(argv[7]);
     atomNum=atoi(argv[8]);
   } else if(mode==3) {
     slope=atof(argv[6]); rangeX=atof(argv[7]); numX=atoi(argv[8]);
     minY=atof(argv[9]); maxY=atof(argv[10]); numY=atoi(argv[11]);
     minZ=atof(argv[12]); maxZ=atof(argv[13]);
   } else ;

     for(step=initial; step<=final; step+=timeStep) 
     {
       sharePNum=0;

       sprintf(fileName,"%dParticle%d.h5",species,step);
       if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
       } else ;

       sprintf(dataName,"totalCnt");
       if(myrank==0)  
         restoreIntMeta(fileName,dataName,&totalCnt,1);
       else ;
       MPI_Barrier(MPI_COMM_WORLD);
       MPI_Bcast(&totalCnt,nTasks,MPI_INT,0,MPI_COMM_WORLD);

       calSub(totalCnt,&cntSub,&start);
       data = (double *)malloc(cntSub*9*sizeof(double ));
       sprintf(dataName,"%d",species);
       restoreData(fileName,dataName,totalCnt,cntSub,start,data,0,9);

       if(mode==1 || mode==2) {
         saveParticleFile(data,minPx,minX,cntSub,step,species,minE,maxE,atomNum); 
       } else if(mode==3) {
         maxX=slope*step; minX=maxX-rangeX;
         saveDensityFile(data,minX,maxX,numX,minY,maxY,numY,minZ,maxZ,step,species,cntSub);
       }

       free(data);
       MPI_Barrier(MPI_COMM_WORLD);
     }

   free(recv);
   MPI_Finalize();
}

void saveDensityFile(double *data,double minX,double maxX,int numX,double minY,double maxY,int numY,double minZ,double maxZ,int step,int species,int cntSub)
{
   int i,j,n,ii,jj,shareNum,start,dataNum=9;
   double x,y,z,weight,density,dx,dy,dz,wx[2],wy[2];
   double *recvData,*sendData,**nn,lambda3;
   char fileName[100];
   FILE *out;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   nn=(double **)malloc((numX+1)*sizeof(double *));
   for(i=0; i<numX+1; i++)
     nn[i]=(double *)malloc((numY+1)*sizeof(double ));
   for(i=0; i<numX+1; i++)
     for(j=0; j<numY+1; j++)
       nn[i][j]=0.0;
   minX=0; maxX=4e-5;
   dx=(maxX-minX)/(double)(numX);
   dy=(maxY-minY)/(double)(numY);
   dz=maxZ-minZ;
//   lambda3=0.8e-6*0.8e-6*0.8e-6;
   for(n=0; n<cntSub; n++)  {
     x=data[n*dataNum+0]; y=data[n*dataNum+1]; z=data[n*dataNum+2];
     weight=data[n*dataNum+8];
     i=(int)((x-minX)/dx);
     j=(int)((y-minY)/dy);
     wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
     wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
     if(i>=0 && i<numX && j>=0 && j<numY && z>=minZ && z<maxZ) {
       for(ii=0; ii<2; ii++)
         for(jj=0; jj<2; jj++) 
           nn[i+ii][j+jj]+=wx[ii]*wy[jj]*weight/dx/dy/dz;
     } else ;
   }
  
   shareNum=(numX+1)*(numY+1);
   sendData=(double *)malloc(shareNum*sizeof(double ));
   recvData=(double *)malloc(shareNum*sizeof(double ));

   start=0;
   for(i=0; i<numX+1; i++) {
     for(j=0; j<numY+1; j++) sendData[start+j]=nn[i][j];
     start+=numY+1;
   }

   if(myrank!=0)  {
     MPI_Send(sendData,shareNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
   } else {
     for(n=1; n<nTasks; n++)  {
       MPI_Recv(recvData,shareNum,MPI_DOUBLE,n,n,MPI_COMM_WORLD,&status);
       start=0;
       for(i=0; i<numX+1; i++) {
         for(j=0; j<numY+1; j++) nn[i][j]+=recvData[start+j];
         start+=numY+1;
       }
     }
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if(myrank==0) {
     sprintf(fileName,"%dDenXY%d",species,step);
     out=fopen(fileName,"w");
     for(i=0; i<numX+1; i++) { 
       for(j=0; j<numY+1; j++) {              
         x=minX+i*dx;
         y=minY+j*dy;
         density=nn[i][j];
         fprintf(out,"%g %g %g\n",x,y,density);
       }
       fprintf(out,"\n");
     }
     fclose(out);
     printf("%s is made\n",fileName);
   }
   free(sendData); free(recvData);
   for(i=0; i<numX+1; i++) free(nn[i]); free(nn);
}

void saveParticleFile(double *data,double minPx,double minX,int cntSub,int step,int species,double minE,double maxE,int atomNum)
{
   int i,n,index,sharePNum,dataNum=9,intId;
   int *recvDataCnt;
   double x,y,z,px,py,pz,id,core,weight,energy,mass;
   double *recvData,*sendData;
   char fileName[100];
   FILE *out;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   mass=whatMass(atomNum);
   recvDataCnt=(int *)malloc(nTasks*sizeof(int ));

   sharePNum=0;
   for(i=0; i<cntSub; i++)  {
     if(data[i*dataNum+3]>minPx) sharePNum+=1; else;
   }
   for(i=1; i<nTasks; i++)   {
     if(myrank==i)
       MPI_Send(&sharePNum,1,MPI_INT,0,myrank,MPI_COMM_WORLD);
   }
   if(myrank==0)  {
     for(i=1; i<nTasks; i++)
       MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&status);
     recvDataCnt[0]=sharePNum;
   }
   MPI_Barrier(MPI_COMM_WORLD);
if(myrank==0)  {
  for(i=0; i<nTasks; i++)
    printf("i=%d, recvDataCnt=%d\n",i,recvDataCnt[i]);
}  
   for(i=1; i<nTasks; i++)   {
     if(myrank==i)  {
       sendData=(double *)malloc(sharePNum*dataNum*sizeof(double ));
       index=0;
       for(n=0; n<cntSub; n++)  {
         if(data[n*dataNum+3]>minPx) {
           sendData[index*dataNum+0]=data[n*dataNum+0];
           sendData[index*dataNum+1]=data[n*dataNum+1];
           sendData[index*dataNum+2]=data[n*dataNum+2];
           sendData[index*dataNum+3]=data[n*dataNum+3];
           sendData[index*dataNum+4]=data[n*dataNum+4];
           sendData[index*dataNum+5]=data[n*dataNum+5];
           sendData[index*dataNum+6]=data[n*dataNum+6];
           sendData[index*dataNum+7]=data[n*dataNum+7];
           sendData[index*dataNum+8]=data[n*dataNum+8];
           index++;
         }  else	;
       }
       MPI_Send(sendData,sharePNum*dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
       free(sendData);
     } else;
   }

   if(myrank==0)  {
     sprintf(fileName,"%dParticle%d",species,step);
     out=fopen(fileName,"w");
     for(n=0; n<cntSub; n++)  {
       if(data[n*dataNum+3]>minPx) {
         x=data[n*dataNum+0];
         y=data[n*dataNum+1];
         z=data[n*dataNum+2];
         px=data[n*dataNum+3];
         py=data[n*dataNum+4];
         pz=data[n*dataNum+5];
         id=data[n*dataNum+6];
         core=data[n*dataNum+7];
         weight=data[n*dataNum+8];
         energy=0.511*mass*(sqrt(1.0+px*px+py*py+pz*pz)-1.0);
         intId=(int)id;
         if(energy>=minE && energy<maxE && x>minX)
           fprintf(out,"%.10g %g %g %g %g %g %.15g %g %g\n",x,y,z,px,py,pz,id,core,weight);
         else ;
       } else	;
     }
     fclose(out);
printf("minX=%g\n",minX);
     for(i=1; i<nTasks; i++)  {
       recvData=(double *)malloc(recvDataCnt[i]*dataNum*sizeof(double ));
       MPI_Recv(recvData,recvDataCnt[i]*dataNum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
       out=fopen(fileName,"a");
       for(n=0; n<recvDataCnt[i]; n++)  {
         x=recvData[n*dataNum+0];
         y=recvData[n*dataNum+1];
         z=recvData[n*dataNum+2];
         px=recvData[n*dataNum+3];
         py=recvData[n*dataNum+4];
         pz=recvData[n*dataNum+5];
         id=recvData[n*dataNum+6];
         core=recvData[n*dataNum+7];
         weight=recvData[n*dataNum+8];
         energy=0.511*mass*(sqrt(1.0+px*px+py*py+pz*pz)-1.0);
         if(energy>=minE && energy<maxE && x>minX)
           fprintf(out,"%g %g %g %g %g %g %.15g %g %g\n",x,y,z,px,py,pz,id,core,weight);
         else ;
       }
       fclose(out);
       free(recvData);
     }
     printf("%s is made\n",fileName);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   free(recvDataCnt);
}

void restore1Data(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC)
{
   hid_t file_id,dset_id,plist_id;
   hid_t filespace,memspace;
   hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
   herr_t ierr;

   //open file
   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   //set dataset
   dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
   //file space
   dimsf[0]=totalCnt;
   dimsf[1]=9;
   filespace=H5Screate_simple(2,dimsf,NULL);
   //memory space
   dimsf[0]=cntSub;
   dimsf[1]=1;
   memspace=H5Screate_simple(2,dimsf,NULL);

   stride[0]=1;   stride[1]=1;
   count[0]=1;    count[1]=1;

   //hyperslab in file space
   block[0]=cntSub;  block[1]=1;
   offSet[0]=start;  offSet[1]=startC;
   H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);
   //hyperslab in memory space
   offSet[0]=0;  offSet[1]=0;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);
   //read data
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
   H5Dclose(dset_id);
   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Fclose(file_id);
}
//lala
void restoreData(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC,int columns)
{
   hid_t file_id,dset_id,plist_id;
   hid_t filespace,memspace;
   hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
   herr_t ierr;

   //open file
   plist_id=H5Pcreate(H5P_FILE_ACCESS);
   H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
   file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
   H5Pclose(plist_id);

   //set dataset
   dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

   //file space
   dimsf[0]=totalCnt;
   dimsf[1]=9;
   filespace=H5Screate_simple(2,dimsf,NULL);
   //memory space
   dimsf[0]=cntSub;
   dimsf[1]=columns;
   memspace=H5Screate_simple(2,dimsf,NULL);

   stride[0]=1;   stride[1]=1;
   count[0]=1;    count[1]=1;

   //hyperslab in file space
   block[0]=cntSub;  block[1]=columns;
   offSet[0]=start;  offSet[1]=startC;
   H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);
   //hyperslab in memory space
   offSet[0]=0;      offSet[1]=0;
   H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);
   //read data
   plist_id=H5Pcreate(H5P_DATASET_XFER);
   H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

   H5Pclose(plist_id);
   H5Sclose(memspace);
   H5Sclose(filespace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}


void calSub(int totalCnt,int *cntSub,int *start)
{
   int i,sub,remain,rank,tmp,*recv,subCnt;;
   int myrank, nTasks;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   recv=(int *)malloc(nTasks*sizeof(int ));

   sub=totalCnt/nTasks;
   remain=totalCnt%nTasks;
   for(rank=0; rank<nTasks; rank++) {
     if(rank<remain)  tmp=sub+1;
     else             tmp=sub;
     if(myrank==rank)  subCnt=tmp; else;
   }
   *cntSub=subCnt;
   MPI_Gather(&subCnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

   tmp=0;
   for(i=0; i<myrank; i++) tmp+=recv[i];
   *start=tmp;
  
   free(recv);
}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

double whatMass(int atomNum)
{
  float eMassU;

  eMassU=5.485799e-4;
  if(atomNum==0)        return 1.0;
  else if(atomNum==1)   return 1.0/eMassU;
  else if(atomNum==2)   return 4.0/eMassU;
  else if(atomNum==6)   return 12.0/eMassU;
  else { printf("no list in atomNum data\n"); exit(0); }
}


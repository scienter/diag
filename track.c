#include "stdio.h"
#include "stdlib.h"
#include <string.h>
#include "malloc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "mpi.h"
#include "math.h"

struct Head_Track {
   struct Track *tr;
};

struct Track {
	double id,core,weight;
	int index;
	struct Track *next;
};

struct Head_type {
   struct Particle *pt;
};

struct Particle {
	int step;
	double r[3];
	double u[3];
	struct Particle *next;
};



void restoreData(char *fileName,char *dataName,int totalCnt,int cntSub,int start,double *data,int startC,int columns);
void calSub(int totalCnt,int *cntSub,int *start);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveFieldComp(double *data,char *fileName,char *dataName,int nx,int ny);

void main(int argc, char *argv[])
{
   int species,initial,final,timeStep,dataNum=9,FLAG;
   int i,j,n,step,cnt,reCnt,totalCnt,cntSub,start,index,linkIdx,stepCnt,stepIdx;
   double minPx,id,core,px,dt,dz,weight,velocityC;
   double *data,*idcore,*idcore_0,*recv,*idcore2,*recv2;
	int *coreCnt,tmpCnt,recvCnt,sendCnt,rank,*trackCore;
   FILE *out;
   char fileName[100],dataName[100];
	double *trackX,*trackY,*trackZ,*trackUx,*trackUy,*trackUz,*trackWe,*trackId;
	struct Track *tr,*NewTr,*tmpTr,*prevTr;
	struct Head_Track *HTr;
	

   int myrank, nTasks;
   MPI_Status status;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   HTr = malloc(sizeof(struct Head_Track));
	HTr->tr=NULL;

   if(argc < 6)   {
     if(myrank==0) {
       printf("mpirun -np [ ] track species initial final step minPx dz\n");
     } else ;
     MPI_Finalize();  
     exit(0);
   }  else;

   species=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   timeStep=atoi(argv[4]);
   minPx=atof(argv[5]);
   dz=atof(argv[6]);
	velocityC=299792458;
	dt=dz/velocityC;

	for(step=initial; step<=final; step+=timeStep) 
	{
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
      MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

		calSub(totalCnt,&cntSub,&start);
     	data = (double *)malloc(cntSub*9*sizeof(double ));
     	sprintf(dataName,"%d",species);
     	restoreData(fileName,dataName,totalCnt,cntSub,start,data,0,dataNum);

		for(i=0; i<cntSub; i++) {
			px=data[i*dataNum+3];
			if(px>minPx) {
				id=data[i*dataNum+6];
				core=data[i*dataNum+7];
				weight=data[i*dataNum+8];

				FLAG=0;
				tr=HTr->tr;
				while(tr) {
					if(id==tr->id && core==tr->core) {
						FLAG=1;
						break;
					}	else ;
					tr=tr->next;
				}

				if(FLAG==0) {
					NewTr=malloc(sizeof(struct Track));
					NewTr->next = HTr->tr;
					HTr->tr = NewTr;
					NewTr->id = id;
					NewTr->core = core;
					NewTr->weight = weight;
				} else ;
			} else ;
		}

   	free(data);
     	MPI_Barrier(MPI_COMM_WORLD);

		cnt=0;
		tr=HTr->tr;
		while(tr) {
			cnt++;
			tr=tr->next;
		}
//printf("step=%d, myrank=%d, cnt=%d\n",step,myrank,cnt);

		if(myrank==0) {
			printf("Constructing id and core, step%d is done.\n",step); 
		}	else ;
  	}

	//constructing track memory
	cnt=0;
	tr=HTr->tr;
	while(tr) {
		cnt++;
		tr=tr->next;
	}
	sendCnt=cnt;
	coreCnt=(int *)malloc(nTasks*sizeof(int ));
	for(i=0; i<nTasks; i++)		coreCnt[i]=0;
	if(myrank!=0) 
		MPI_Send(&sendCnt,1,MPI_INT,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(&recvCnt,1,MPI_INT,rank,rank,MPI_COMM_WORLD,&status);
			coreCnt[rank]=recvCnt;
		}
	}
   MPI_Barrier(MPI_COMM_WORLD);

	idcore=(double *)malloc(3*cnt*sizeof(double ));
	index=0;
	tr=HTr->tr;
	while(tr) {
		idcore[index*3+0]=tr->id;
		idcore[index*3+1]=tr->core;
		idcore[index*3+2]=tr->weight;
		index+=1;
		tr=tr->next;
	}

	if(myrank!=0) 
		MPI_Send(idcore,cnt*3,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			tmpCnt=coreCnt[rank];
			recv=(double *)malloc(3*tmpCnt*sizeof(double ));
			MPI_Recv(recv,tmpCnt*3,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);

			for(j=0; j<tmpCnt; j++) {
				id=recv[j*3+0];
				core=recv[j*3+1];
				weight=recv[j*3+2];
				FLAG=0;
				for(i=0; i<cnt; i++) {
					if(id==idcore[i*3+0] && core==idcore[i*3+1]) {
						FLAG=1;
					} else ;
				}
				if(FLAG==0) {
					NewTr=malloc(sizeof(struct Track));
					NewTr->next = HTr->tr;
					HTr->tr = NewTr;
					NewTr->id = id;
					NewTr->core = core;
					NewTr->weight = weight;
				} else ;
			}
			free(recv);
		}
	}
	free(idcore);
   MPI_Barrier(MPI_COMM_WORLD);


	if(myrank==0) {
		cnt=0;
		tr=HTr->tr;
		while(tr) {
			cnt++;
			tr=tr->next;
		}
	} else ;
	MPI_Bcast(&cnt,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
		

	// test duplication
	idcore_0=(double *)malloc(3*cnt*sizeof(double ));
	idcore2=(double *)malloc(3*cnt*sizeof(double ));
	
	if(myrank==0) {
		index=0;
		tr=HTr->tr;
		while(tr) {
			idcore_0[index*3+0]=tr->id;
			idcore_0[index*3+1]=tr->core;
			idcore_0[index*3+2]=tr->weight;
			index+=1;
			tr=tr->next;
		}	
		memcpy(idcore2,idcore_0,sizeof(double)*3*cnt);

		for(i=0; i<cnt; i++) {
			id=idcore2[i*3];
			core=idcore2[i*3+1];
			for(j=0; j<cnt; j++) {
				if(i!=j && id==idcore_0[j*3] && core==idcore_0[j*3+1] && id>0) {
					idcore2[j*3]=0;
				} else ;
			}
		}		//End of for(cnt)

		reCnt=0;
		for(i=0; i<cnt; i++) {
			if(idcore2[i*3+0]>0) reCnt++; else;
		}
		printf("\n total tracking count is %d.\n\n",reCnt);

	}	else ;
   MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(idcore2,cnt*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&reCnt,1,MPI_INT,0,MPI_COMM_WORLD);

	tr=HTr->tr;
	while(tr) {
		HTr->tr=tr->next;
		tr->next=NULL;
		free(tr);
		tr=HTr->tr;
	}

//-------------------------- tracking particle -----------------------------------//
	stepCnt=(final-initial)/timeStep+1;
	trackX=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackY=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackZ=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackUx=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackUy=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackUz=(double *)malloc(reCnt*stepCnt*sizeof(double ));
	trackId=(double *)malloc(reCnt*sizeof(double ));
	trackCore=(int *)malloc(reCnt*sizeof(int ));
	trackWe=(double *)malloc(reCnt*sizeof(double ));
	for(i=0; i<reCnt*stepCnt; i++) { 
		trackX[i]=0.0;
		trackY[i]=0.0;
		trackZ[i]=0.0;
		trackUx[i]=0.0;
		trackUy[i]=0.0;
		trackUz[i]=0.0;
	}
	for(i=0; i<reCnt; i++) { 
		trackId[i]=0.0;
		trackCore[i]=0;
		trackWe[i]=0.0;
	}
	index=0;
	for(i=0; i<cnt; i++) {
		if(idcore2[i*3]>0) {
			trackId[index]=idcore2[i*3+0];
			trackCore[index]=idcore2[i*3+1];
			trackWe[index]=idcore2[i*3+2];
			index++;
		} else ;
	}

	// save track memory in rinked list
	for(step=initial; step<=final; step+=timeStep) 
	{
		stepIdx=(step-initial)/timeStep;
		
		cnt=0;
		tr=HTr->tr;
		while(tr) {
			tr=tr->next;
			cnt++;
		}
		if(cnt>0) {
			printf("myrank=%d, step=%d, cnt=%d\n",myrank,step,cnt); 
			exit(0); 
		} else;

		// save track memory in rinked list 
		for(i=0; i<reCnt; i++) {
			NewTr=malloc(sizeof(struct Track));
			NewTr->next = HTr->tr;
			HTr->tr = NewTr;
			NewTr->id = trackId[i];
			NewTr->core = trackCore[i];
			NewTr->weight = trackWe[i];
			NewTr->index = i;
		}

		sprintf(fileName,"%dParticle%d.h5",species,step);
      sprintf(dataName,"totalCnt");
     	if(myrank==0)  
      	restoreIntMeta(fileName,dataName,&totalCnt,1);
     	else ;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&totalCnt,1,MPI_INT,0,MPI_COMM_WORLD);

		calSub(totalCnt,&cntSub,&start);
     	data = (double *)malloc(cntSub*9*sizeof(double ));
     	sprintf(dataName,"%d",species);
     	restoreData(fileName,dataName,totalCnt,cntSub,start,data,0,dataNum);

		for(n=0; n<cntSub; n++) {
			id=data[n*dataNum+6];
			core=data[n*dataNum+7];
			px=data[n*dataNum+3];
			if(px>minPx) {
				linkIdx=1;
				tr=HTr->tr;
				while(tr) {
					if(linkIdx==1)	
						prevTr=tr;	
					else ;
	
					FLAG=0;
					if(id==tr->id && core==tr->core) {
						index=tr->index;
						trackX[index*stepCnt+stepIdx]=data[n*dataNum+0];
						trackY[index*stepCnt+stepIdx]=data[n*dataNum+1];
						trackZ[index*stepCnt+stepIdx]=data[n*dataNum+2];
						trackUx[index*stepCnt+stepIdx]=data[n*dataNum+3];
						trackUy[index*stepCnt+stepIdx]=data[n*dataNum+4];
						trackUz[index*stepCnt+stepIdx]=data[n*dataNum+5];
						FLAG=1;
					}	else ;
				
					if(FLAG==1) {
						if(linkIdx==1) {
							tmpTr=tr->next;
							HTr->tr=tmpTr;
							tr->next=NULL;
							free(tr);
							tr=HTr->tr;
							linkIdx=1;
						} else {
							prevTr->next = tr->next;
							tr->next=NULL;
							free(tr);
							tr=prevTr->next;
						}
						break;
					} else {
						prevTr = tr;
						tr = tr->next;
						linkIdx++;
					}
				
				}		// End of while(tr)
			} else ; 	// End of for(cntSub)
			
		}
		
		tr=HTr->tr;
		while(tr) {
			HTr->tr=tr->next;
			tr->next=NULL;
			free(tr);
			tr=HTr->tr;
		}

		if(myrank==0) {
			printf("Tracking, step%d is done.\n",step); 
		}	else ;

	}
   MPI_Barrier(MPI_COMM_WORLD);

  	recv2 = (double *)malloc(reCnt*stepCnt*sizeof(double ));

	if(myrank!=0) 
		MPI_Send(trackX,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackX[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) 
		MPI_Send(trackY,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackY[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) 
		MPI_Send(trackZ,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackZ[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) 
		MPI_Send(trackUx,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackUx[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) 
		MPI_Send(trackUy,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackUy[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0) 
		MPI_Send(trackUz,reCnt*stepCnt,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv2,reCnt*stepCnt,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			for(i=0; i<reCnt*stepCnt; i++)	trackUz[i]+=recv2[i];
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);


	hid_t file_id,group_id,dset_id,filespace;
	herr_t hdf_status;
	int *steps;

	steps=(int *)malloc(stepCnt*sizeof(int ));

	if(myrank==0) {
		for(i=0; i<stepCnt; i++)
			steps[i]=initial+i*timeStep;


		sprintf(fileName,"track.h5");
		file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
		H5Fclose(file_id);
		saveIntMeta(fileName,"/particleCnt",&reCnt,1);
		saveIntMeta(fileName,"/stepCnt",&stepCnt,1);
		saveDoubleMeta(fileName,"/dt",&dt,1);
		saveIntMeta(fileName,"/timeSteps",steps,stepCnt);

		saveDoubleMeta(fileName,"trackId",trackId,reCnt);
		saveDoubleMeta(fileName,"trackWeight",trackWe,reCnt);
		saveIntMeta(fileName,"trackCore",trackCore,reCnt);
      saveFieldComp(trackX,fileName,"trackX",reCnt,stepCnt);
      saveFieldComp(trackY,fileName,"trackY",reCnt,stepCnt);
      saveFieldComp(trackZ,fileName,"trackZ",reCnt,stepCnt);
      saveFieldComp(trackUx,fileName,"trackUx",reCnt,stepCnt);
      saveFieldComp(trackUy,fileName,"trackUy",reCnt,stepCnt);
      saveFieldComp(trackUz,fileName,"trackUz",reCnt,stepCnt);
		printf("%s is made.\n",fileName);


	} else ;
	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank==0) {
		out=fopen("trackFile","w");
		for(i=0; i<reCnt; i++)
			fprintf(out,"%.20g %d %g\n",trackId[i],trackCore[i],trackWe[i]);
		fclose(out);
		printf("trackFile is made.\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);



	free(idcore2);
	free(idcore_0);
	free(recv2);
	free(trackX);
	free(trackY);
	free(trackZ);
	free(trackUx);
	free(trackUy);
	free(trackUz);
	free(trackId);
	free(trackCore);
	free(trackWe);
	free(steps);

   MPI_Finalize();
}

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


void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
	hid_t file_id,dset_id,filespace;
	hsize_t metaDim[1];
	herr_t status;

	metaDim[0]=dataCnt;

	file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
	filespace=H5Screate_simple(1,metaDim,NULL);
	dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
	H5Dclose(dset_id);
	H5Sclose(filespace);
	H5Fclose(file_id);
}

void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
	hid_t file_id,dset_id,filespace;
	hsize_t metaDim[1];
	herr_t status;

	metaDim[0]=dataCnt;

	file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
	filespace=H5Screate_simple(1,metaDim,NULL);
	dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
	H5Dclose(dset_id);
	H5Sclose(filespace);
	H5Fclose(file_id);
}

void saveFieldComp(double *data,char *fileName,char *dataName,int nx,int ny)
{
    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[2],count[2],offset[2];

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);

    dimsf[0]=nx;
    dimsf[1]=ny;
    filespace=H5Screate_simple(2,dimsf,NULL);

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

	 H5Dclose(dset_id);
	 H5Sclose(filespace);
	 H5Fclose(file_id);
}


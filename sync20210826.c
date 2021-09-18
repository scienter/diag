#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define m 9.10938291e-31
#define e 1.602e-19
#define c 299792458
#define eps0 8.854187817e-12
#define mu0 1.256637062e-6			// [H/m]
//#define hbar 4.135667696e-15			// [eV .s]
#define hbar 6.582119569e-16		// [eV .s]

struct Particle {
	double E[3],B[3];
	double rNow[3],rNext[3];
	double uOld[3],uNext[3];
	struct Particle *next;
};

struct Domain_type {
	int mode,correction,numeric,testMode;
	int numW;
	int maxN;
	int maxStep,jump;
	int scrN;
	double scrD,scrAngle,dr,screenX0;
	double minW,dW;
	double dt;
	double B0,undL,lambdaU;
};

struct Head_type {
	struct Particle *pt;
};

typedef struct Head_type *Head;
typedef struct Domain_type *Domain;
typedef struct FCOMPLEX {double r,i;} fcomplex;
fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex Complex(double re, double im);
fcomplex Cdiv(fcomplex a, fcomplex b); 
double Cabs(fcomplex z) ;
fcomplex RCmul(double x, fcomplex a) ;


static void terminate(const char *message);
void particlePush(Head head,Domain D);
void calRadiation(Domain D,double ***spectra1,double ***spectra2,double *trackParticle,double weight,int index,int end);
void saveSpectra(Domain D,double ***spectra1,double ***spectra2,int pCnt);
void interpolation(Head head,Domain D,double *E,double *B);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void restoreAttr(char *fileName,char *dataName,char *attrName,double *data);
void restoreData(char *fileName,char *dataName,double *data);
void calSub(int totalCnt,int *cntSub,int *start);
void fresnel(double x,double *C,double *S);


void main(int argc, char *argv[])
{
	int i,j,k,step,cnt,nn,rank;
	int pCnt,start,end,cntSub,dataNum;
	double energy,maxX,gamma0,dz,E0,pickId,pickCore;
	double maxPhotonE,minPhotonE,wc,gy_radius;
	double x,y,z,ux,uy,uz,weight,id,core;
	double ux0,uy0,uz0,t,E[3],B[3],*send,*recv;
	double ***spectra1,***spectra2,*trackParticle;
	struct Particle *New,*p;
	FILE *out,*in;
	char fileName[100],dataName[100];
	Head head;
	Domain D;

	 int myrank, nTasks;	
	 MPI_Status status;

	 MPI_Init(&argc,&argv);
	 MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


	D = malloc(sizeof(struct Domain_type));

   if(argc < 10)   { 
		if(myrank==0) {
	      printf("sync mode correction scrAngle scrD scrN minPhotonE maxPhotonE numPhoton screenX0 numeric\n");
   	   printf("mode 0 : energy B0 Ex undL lambdaU dz maxStep\n");
      	printf("mode 1 (PIC) : testMode(on:0 / off:1) jump\n");
	      printf("mode 2 (PIC) : id core \n");
		} else ;
		MPI_Finalize();
      exit(0);
   } else ;

	D->mode=atoi(argv[1]);
	D->correction=atoi(argv[2]);
   D->scrAngle=atof(argv[3]);
   D->scrD=atof(argv[4]);
   D->scrN=atoi(argv[5]);
	minPhotonE=atof(argv[6]);
	maxPhotonE=atof(argv[7]);
	D->numW=atoi(argv[8]);
	D->screenX0=atof(argv[9]);
	D->numeric=atoi(argv[10]);

	D->jump=1;
	D->testMode=0;

	// set spectra memory
	D->dr=D->scrD/(1.0*D->scrN);
	D->dW=(maxPhotonE-minPhotonE)/(1.0*D->numW)/hbar;
	D->minW=minPhotonE/hbar;

	spectra1 = (double ***)malloc(D->numW*sizeof(double **));
	spectra2 = (double ***)malloc(D->numW*sizeof(double **));
	for(i=0; i<D->numW; i++)  {
		spectra1[i] = (double **)malloc(D->scrN*sizeof(double *));
		spectra2[i] = (double **)malloc(D->scrN*sizeof(double *));
		for(j=0; j<D->scrN; j++)  {
			spectra1[i][j] = (double *)malloc(3*sizeof(double ));
			spectra2[i][j] = (double *)malloc(3*sizeof(double ));
		}
	}
	for(i=0; i<D->numW; i++)  
		for(j=0; j<D->scrN; j++)  
			for(k=0; k<3; k++)  {
				spectra1[i][j][k]=0.0;
				spectra2[i][j][k]=0.0;
			}

	 switch (D->mode) {
	 case 0 :
		energy=atof(argv[11]);
		D->B0=atof(argv[12]);
		E0=atof(argv[13]);
		D->undL=atof(argv[14]);
		D->lambdaU=atof(argv[15]);
		dz=atof(argv[16]);
		D->maxStep=atoi(argv[17]);

		E[0]=E0;		  E[1]=0.0;		  E[2]=0.0;
		B[0]=0.0;		  B[1]=0.0;		  B[2]=D->B0;
		gamma0=energy/0.511;
		D->undL=10000e-6;  //2;
		D->lambdaU=2500e-6; //3.5e-2;

		ux0=sqrt(gamma0*gamma0-1.0);	uy0=uz0=0.0;
		gy_radius=ux0*m*c/e/D->B0;
		wc=3.0/2.0*gamma0*gamma0*gamma0*c/gy_radius;

		D->dt=dz/c;
		printf("energy=%g, gamma0=%g, ux0=%g, B0=%g, gy_radius=%g, critical_energy=%g[eV], maxEnergy=%g\n",energy,gamma0,ux0,D->B0,gy_radius,wc*hbar,1.0/D->dt*hbar/6.28);

		// initialize particle
		head = malloc(sizeof(struct Head_type));
		if(head==NULL) 	terminate("Error in creating head pointer.\n");
		else head->pt = NULL;

		New = malloc(sizeof(struct Particle ));
		New->next = head->pt;
		head->pt = New;
		New->rNow[0] = 0.0;		New->rNow[1] = 0.0; 	New->rNow[2] = 0.0;
		New->rNext[0] = 0.0;		New->rNext[1] = 0.0; 	New->rNext[2] = 0.0;
		New->uOld[0] = ux0;	   New->uOld[1] = uy0;		New->uOld[2] = uz0;
		New->uNext[0] = ux0;	   New->uNext[1] = uy0;		New->uNext[2] = uz0;
		for(i=0; i<3; i++) { 
			New->E[i]=E[i];
			New->B[i]=B[i];
		}

		trackParticle=(double *)malloc(D->maxStep*6*sizeof(double));

		// Main loop
		step=0;
		while(step<D->maxStep) {	 
			p=head->pt;
			while(p) {
				for(i=0; i<3; i++)	trackParticle[step*6+i]=p->rNext[i];
				for(i=0; i<3; i++)  trackParticle[step*6+i+3]=p->uNext[i];
				p=p->next;
			}
				
			interpolation(head,D,E,B);

			particlePush(head,D);

			step++;
		}

		sprintf(fileName,"particle");
		out=fopen(fileName,"w");
		for(i=0; i<D->maxStep; i++) {
			x=trackParticle[i*6+0];
			y=trackParticle[i*6+1];
			z=trackParticle[i*6+2];
			ux=trackParticle[i*6+3];
			uy=trackParticle[i*6+4];
			uz=trackParticle[i*6+5];
			fprintf(out,"%g %g %g %g %g %g\n",x,y,z,ux,uy,uz);
		}
		fclose(out);
		printf("%s is made.\n",fileName);

		p=head->pt;
		while(p) {
			head->pt=p->next;
			p->next=NULL;
			free(p);
			p=head->pt;
		}
		free(head);

		calRadiation(D,spectra1,spectra2,trackParticle,1,0,1);

		break;

	case 1 :
	case 2 :
		if(D->mode==2) {
			pickId=atof(argv[11]);
			pickCore=atof(argv[12]);
		} else if(D->mode==1)  {
			D->testMode=atoi(argv[11]);
			D->jump=atoi(argv[12]);
		} else ;

		sprintf(fileName,"track.h5");
		if(fopen("track.h5","r")==NULL)  {
			printf("%s does not exit.\n",fileName);
			exit(0);
		} else ;

		if(myrank==0) {
			restoreIntMeta(fileName,"particleCnt",&pCnt,1);
			restoreIntMeta(fileName,"stepCnt",&D->maxStep,1);
			restoreDoubleMeta(fileName,"dt",&D->dt,1);
		}  else ;
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&pCnt,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&D->maxStep,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&D->dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
		start=0;
		calSub(pCnt,&cntSub,&start);
		cntSub=pCnt/nTasks;
		end=start+cntSub;
		if(myrank==0) {
			printf("missing particles are %d in total %d\n",pCnt-cntSub*nTasks,pCnt);
		} else ;

		printf("start=%d, end=%d, cntSub=%d, dt=%g, maxStep=%d\n",start,end,cntSub,D->dt,D->maxStep);

		trackParticle=(double *)malloc(D->maxStep*6*sizeof(double));

		if(D->mode==1) {
			if(D->testMode==1) end=start+1;
			else ;
			for(i=start; i<end; i++) 	{
				sprintf(dataName,"%d",i);
				restoreData(fileName,dataName,trackParticle);
				restoreAttr(fileName,dataName,"weight",&weight);
				restoreAttr(fileName,dataName,"id",&id);
				restoreAttr(fileName,dataName,"core",&core);
			
				printf("myrank=%d, index=%d,id=%g, core=%g, weight=%g\n",myrank,i,id,core,weight);
				calRadiation(D,spectra1,spectra2,trackParticle,weight,i,end);

				if(myrank==0) printf("%d/%d is done.\n",i,end-1); else ;
			} 
		} else if(D->mode==2) {
			for(i=start; i<end; i+=D->jump) 	{
				sprintf(dataName,"%d",i);
				restoreAttr(fileName,dataName,"weight",&weight);
				restoreAttr(fileName,dataName,"id",&id);
				restoreAttr(fileName,dataName,"core",&core);

				if(pickId==id && pickCore==core) {
					restoreData(fileName,dataName,trackParticle);
					printf("myrank=%d, index=%d,id=%g, core=%g, weight=%g\n",myrank,i,id,core,weight);
					calRadiation(D,spectra1,spectra2,trackParticle,weight*D->jump,i,end);
				}	else ;
			}
		}	else ;

		break;
	}			//End of switch (mode)

   dataNum=D->numW*D->scrN*3*2;
   send=(double *)malloc(dataNum*sizeof(double ));
   recv=(double *)malloc(dataNum*sizeof(double ));
   start=0;
   for(i=0; i<D->numW; i++)
		for(j=0; j<D->scrN; j++)
			for(k=0; k<3; k++)  {
				send[start+0]=spectra1[i][j][k];
				send[start+1]=spectra2[i][j][k];
				start+=2;
			}
	if(myrank!=0)
		MPI_Send(send,dataNum,MPI_DOUBLE,0,myrank,MPI_COMM_WORLD);
	else {
		for(rank=1; rank<nTasks; rank++) {
			MPI_Recv(recv,dataNum,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,&status);
			start=0;
			for(i=0; i<D->numW; i++)
				for(j=0; j<D->scrN; j++)
					for(k=0; k<3; k++)  {
						spectra1[i][j][k]+=recv[start+0];
						spectra2[i][j][k]+=recv[start+1];
						start+=2;
					}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(myrank==0) 	saveSpectra(D,spectra1,spectra2,pCnt);	else ;

	free(recv);
	free(send);

	for(i=0; i<D->numW; i++) {
		for(j=0; j<D->scrN; j++) {
			free(spectra1[i][j]);
			free(spectra2[i][j]);
		}
		free(spectra1[i]);
		free(spectra2[i]);
	}
	free(spectra1);
	free(spectra2);
	free(trackParticle);
	free(D);

	MPI_Finalize();
	 

}


void calRadiation(Domain D,double ***spectra1,double ***spectra2,double *trackParticle,double weight,int index,int end)
{
	int i,ii,numW,idx,startI,nn;
	double coef,R2,invR,minW,dW,w,invGamma,gamma,prevG,nextG;
	double t,dt,Phi_p,Phi_m,absChi2,sign,sgn,abs;
	double u[3],n[3],delR[3],r[3],Ireal[3],Iimag[3],SCREEN[3],nextU[3],prevU[3];
	double r1[3],r2[3],u1[3],chi0,chi1,chi2,Th_p,Th_m,Psi_p[3],Psi_m[3];
	double fresC_p,fresC_m,fresS_p,fresS_m,z,I0,I1,I2;
	double sv1[3],sv2[3],cosX,sinX,nDotI,absI,tmp;
	double below,nDotR,nDotV,phase1,phase2,cosP1,sinP1,cosP2,sinP2;
	struct Particle *p;

	int myrank, nTasks;	

	MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	numW = D->numW;
	minW = D->minW;
	dW = D->dW;
	dt=D->dt;

	if(minW==0.0) 	startI=1;
	else				startI=0;

	for(idx=0; idx<D->scrN; idx++)  
	{
		SCREEN[0]=D->screenX0;
		SCREEN[1]=idx*D->dr*cos(D->scrAngle*M_PI/180.0);
		SCREEN[2]=idx*D->dr*sin(D->scrAngle*M_PI/180.0);

		for(ii=startI; ii<numW; ii++) 
		{
			w=minW+ii*dW;
				
			for(nn=1; nn<D->maxStep-1; nn++) 
			{
				t=nn*dt;

				for(i=0; i<3; i++)  nextU[i]=trackParticle[(nn+1)*6+i+3];
				if(nextU[1]!=0) {
					for(i=0; i<3; i++)  prevU[i]=trackParticle[(nn+0)*6+i+3];
					for(i=0; i<3; i++)  u[i]=0.5*(nextU[i]+prevU[i]);
					prevG=1.0; 	for(i=0; i<3; i++) prevG+=prevU[i]*prevU[i];
					nextG=1.0; 	for(i=0; i<3; i++) nextG+=nextU[i]*nextU[i];
					prevG=sqrt(prevG);	nextG=sqrt(nextG);
					gamma=0.5*(prevG+nextG);
					invGamma = 1.0/gamma;

					for(i=0; i<3; i++) {
						u1[i]=(nextU[i]-prevU[i])/(1.0*dt);
						r1[i]=u[i]*c*invGamma;
						r2[i]=u1[i]*0.5*c*invGamma;
						r[i]=trackParticle[(nn)*6+i];
					}

					for(i=0; i<3; i++) delR[i]=SCREEN[i]-r[i];
					R2=0.0;					for(i=0; i<3; i++)	R2+=delR[i]*delR[i];
					invR=1.0/sqrt(R2);	for(i=0; i<3; i++)	n[i]=delR[i]*invR;

					nDotV=0.0;     for(i=0; i<3; i++) { nDotV+=n[i]*u[i]; }

					if(D->numeric==0) {
						chi0=w*t;		
						chi1=w;
						chi2=0.0;	
						for(i=0; i<3; i++) {
							chi0-=w/c*n[i]*r[i];
							chi1-=w/c*n[i]*r1[i];
							chi2-=w/c*n[i]*r2[i];
						}
						for(i=0; i<3; i++) {
							Ireal[i]=0.0;
							Iimag[i]=0.0;
						}

						if(chi2!=0.0) {			
							absChi2=fabs(chi2);
	
							if(absChi2*dt*dt<1) {
								z=chi1*dt*0.5;
								I0=sin(z)/z*dt;
								I1=dt/chi1*(sin(z)/z-cos(z));
								I2=dt*dt*dt*0.25*sin(z)/z-2.0/chi1*I1;					
								for(i=0; i<3; i++) {
									Ireal[i]=u[i]*I0*invGamma;
									Iimag[i]=(u1[i]*I1+chi2*u[i]*I2)*invGamma;
								}
							} 
							else {					
								sign=copysign(1.0,chi2);
								tmp=1.0/sqrt(2*M_PI*absChi2);
								Th_p=(chi1+absChi2*dt)*tmp;
								Th_m=(chi1-absChi2*dt)*tmp;

								tmp=sqrt(2*M_PI/absChi2);
								cosX=cos(chi1*chi1*0.25/absChi2);
								sinX=sin(chi1*chi1*0.25/absChi2);
								for(i=0; i<3; i++) {
									Psi_p[i]=tmp*(2*absChi2*u[i]-sign*chi1*u1[i])*cosX;
									Psi_m[i]=tmp*(2*absChi2*u[i]-sign*chi1*u1[i])*sinX;
								}

								z=Th_p;
								sgn=copysign(1.0,z);
								abs=fabs(z);
								fresC_p=( sin(M_PI*0.5*z*z)/M_PI/(abs+20*M_PI*exp(-200*M_PI*sqrt(abs)))+8.0/25.0*(1-exp(-0.69*M_PI*abs*abs*abs))+2.0/25.0*(1-exp(-4.5*M_PI*abs*abs))+0.1*(1-exp(-1.55294068198794*M_PI*abs)) )*sgn;
								fresS_p=( -cos(M_PI*0.5*z*z)/M_PI/(abs+16.7312774552827*M_PI*exp(-1.57638860756614*M_PI*sqrt(abs)))+8.0/25.0*(1-exp(-0.608707749430681*M_PI*abs*abs*abs))+2.0/25.0*(1-exp(-1.71402838165388*M_PI*abs*abs))+0.1*(1-exp(-0.9*M_PI*abs)) )*sgn;
								z=Th_m;
								sgn=copysign(1.0,z);
								abs=fabs(z);
								fresC_m=( sin(M_PI*0.5*z*z)/M_PI/(abs+20*M_PI*exp(-200*M_PI*sqrt(abs)))+8.0/25.0*(1-exp(-0.69*M_PI*abs*abs*abs))+2.0/25.0*(1-exp(-4.5*M_PI*abs*abs))+0.1*(1-exp(-1.55294068198794*M_PI*abs)) )*sgn;
								fresS_m=( -cos(M_PI*0.5*z*z)/M_PI/(abs+16.7312774552827*M_PI*exp(-1.57638860756614*M_PI*sqrt(abs)))+8.0/25.0*(1-exp(-0.608707749430681*M_PI*abs*abs*abs))+2.0/25.0*(1-exp(-1.71402838165388*M_PI*abs*abs))+0.1*(1-exp(-0.9*M_PI*abs)) )*sgn;

//						fresnel(Th_p,&fresC_p,&fresS_p);
//					fresnel(Th_m,&fresC_m,&fresS_m);

								Phi_p=0.25*dt*dt*absChi2+sign*0.5*dt*chi1;
								Phi_m=0.25*dt*dt*absChi2-sign*0.5*dt*chi1;
								for(i=0; i<3; i++) {
									Ireal[i]=0.25*invGamma/absChi2*(Psi_p[i]*(fresC_p-fresC_m)+Psi_m[i]*(fresS_p-fresS_m)+2*u1[i]*(sin(Phi_p)-sin(Phi_m)));
									Iimag[i]=0.25*invGamma/absChi2*(Psi_p[i]*(fresS_p-fresS_m)-Psi_m[i]*(fresC_p-fresC_m)-2*u1[i]*(cos(Phi_p)-cos(Phi_m)))*sign;
								}
							}		// Enf of Fresnel integral
						}
						else {
							for(i=0; i<3; i++) {
								Ireal[i]=2*u[i]/chi1*invGamma*sin(chi1*0.5*dt);
								Iimag[i]=0.0;
							}
						}

						// cross product
						cosX=cos(chi0);
						sinX=sin(chi0);

						for(i=0; i<3; i++) { 
							sv1[i]=Iimag[i]*cosX+Ireal[i]*sinX;
							sv2[i]=Iimag[i]*sinX-Ireal[i]*cosX;
						}

						nDotI=0.0;	for(i=0; i<3; i++) { nDotI+=n[i]*sv1[i]; }		
						for(i=0; i<3; i++)  
							spectra1[ii][idx][i]+=(n[i]*nDotI-sv1[i])*w;
					
						nDotI=0.0;	for(i=0; i<3; i++) { nDotI+=n[i]*sv2[i]; }		
						for(i=0; i<3; i++) 
							spectra2[ii][idx][i]+=(n[i]*nDotI-sv2[i])*w;

					}		  //End of numeric==0
					else if(D->numeric==1) {
						phase1=w*t;
						for(i=0; i<3; i++) phase1-=w/c*n[i]*r[i];
						cosP1=cos(phase1);
						sinP1=sin(phase1);
	
						for(i=0; i<3; i++) {
							spectra1[ii][idx][i]+=w*(n[i]*nDotV-u[i])*sinP1*weight*invGamma*dt;
							spectra2[ii][idx][i]-=w*(n[i]*nDotV-u[i])*cosP1*weight*invGamma*dt;
						}
					}


					// First correction
					if(D->correction==1) {
//						for(i=0; i<3; i++) { u[i]=0.5*(trackParticle[(nn+2)*6+i+3]+trackParticle[(nn+1)*6+i+3]); }
						nDotV=0.0;		for(i=0; i<3; i++) { nDotV+=n[i]*nextU[i]; }
						nDotR=0.0;		for(i=0; i<3; i++) { nDotR+=n[i]*(r[i]+r1[i]*0.5*dt+r2[i]*0.25*dt*dt); }		
						phase1 = w*(t+dt*0.5-nDotR/c);
						cosP1=cos(phase1); sinP1=sin(phase1);
						for(i=0; i<3; i++) {
							spectra1[ii][idx][i]+=(n[i]*nDotV-nextU[i])/(nextG-nDotV)*cosP1;
							spectra2[ii][idx][i]+=(n[i]*nDotV-nextU[i])/(nextG-nDotV)*sinP1;
						}

//					for(i=0; i<3; i++) { u[i]=0.5*(trackParticle[(nn+1)*6+i+3]+trackParticle[(nn+0)*6+i+3]); }
						nDotV=0.0;		for(i=0; i<3; i++) { nDotV+=n[i]*prevU[i]; }
						nDotR=0.0;		for(i=0; i<3; i++) { nDotR+=n[i]*(r[i]-r1[i]*0.5*dt+r2[i]*0.25*dt*dt); }		
						phase2 = w*(t-dt*0.5-nDotR/c);
						cosP2=cos(phase2); sinP2=sin(phase2);
						for(i=0; i<3; i++) {
							spectra1[ii][idx][i]-=(n[i]*nDotV-prevU[i])/(prevG-nDotV)*cosP2;
							spectra2[ii][idx][i]-=(n[i]*nDotV-prevU[i])/(prevG-nDotV)*sinP2;
						}
					} else ;

				} else ;	//End of if(nextU[1]!=0)
		  
			}		// End of for(nn)
		}			// End of for(W)
		if(myrank==0) { printf("screen %d/%d is done. for index=%d/%d\n",idx,D->scrN,index,end); } else ;
	}				// end of (idx)
	


}


void saveSpectra(Domain D,double ***spectra1,double ***spectra2,int pCnt)
{
	int i,ii,numW,idx;
	double minW,dW,w,sp1,sp2,coef,coef2,r,angle,sum[3],I0[3],dA;
	FILE *out1,*out2;
	char name1[100],name2[100];

	numW = D->numW;
	minW = D->minW;
	dW = D->dW;

	coef = mu0*e*e*c/16.0/(M_PI*M_PI*M_PI);
	coef2 = sqrt(coef);
	
	if(D->mode==0 || D->mode==2 || D->testMode==1) pCnt=1; 
	else { pCnt=pCnt/D->jump; }
   sprintf(name1,"%gspectra%d_%d",D->scrAngle,D->correction,pCnt);
   sprintf(name2,"%gspectraSum%d_%d",D->scrAngle,D->correction,pCnt);
	out1=fopen(name1,"w");
	out2=fopen(name2,"w");
	dA=D->dr/D->screenX0;
	for(ii=0; ii<numW; ii++) {
		w=minW+ii*dW;
		for(i=0; i<3; i++) sum[i]=0.0;
		for(idx=0; idx<D->scrN; idx++) {
			r=idx*D->dr;
			angle=r/D->screenX0;
			fprintf(out1,"%g %g %g ",w*hbar,angle,r);
			for(i=0; i<3; i++) {
				sp1=spectra1[ii][idx][i];
				sp2=spectra2[ii][idx][i];
				I0[i]=(sp1*sp1+sp2*sp2)*coef;
				sum[i]+=I0[i];
			}
			fprintf(out1,"%g %g %g",I0[0],I0[1],I0[2]);
			if(D->scrN>1) fprintf(out1,"\n"); else ;
		}
		fprintf(out1,"\n");
		fprintf(out2,"%g %g %g %g\n",w*hbar,sum[0]*dA,sum[1]*dA,sum[2]*dA);
	}
	fclose(out1);
	fclose(out2);
	printf("%s is made.\n",name1);
	printf("%s is made.\n",name2);


}


void interpolation(Head head,Domain D,double *E,double *B)
{
	struct Particle *p;
	double lambdaU,ku,x,B0,L,maxL;
	int num,i;

	lambdaU=D->lambdaU;
	L=D->undL;
	ku=2*M_PI/lambdaU;
	num=(int)(L/lambdaU);
	maxL=num*lambdaU;


	p=head->pt;
	while(p) {
		 x=p->rNext[0];
		 B0=D->B0*sin(ku*x);
		 if(x<lambdaU) p->B[2]=B0*x/lambdaU;
		 else if(x>=maxL+lambdaU) p->B[2]=0.0;
		 else if(x>=maxL) p->B[2]=B0*(1-(x-maxL)/lambdaU);
		 else					p->B[2]=B0;
		 for(i=0; i<3; i++) p->E[i]=E[i];
		 p=p->next;
	}
}

void particlePush(Head head,Domain D)
{
	int i,j;
	struct Particle *p;
	double uMinus[3],uPlus[3],T[3],S[3],M[3][3],u[3];
	double T2,invGamma,gamma2,coefE,coefB,delR,dt;

	coefE=-e/m/c*D->dt*0.5;
	coefB=-e/m*D->dt*0.5;
	dt=D->dt;

	p=head->pt;
	while(p) {
		for(i=0; i<3; i++)	{ 
			p->rNow[i]=p->rNext[i];
			p->uOld[i]=p->uNext[i];
		}
		for(i=0; i<3; i++)	uMinus[i]=p->uOld[i]+coefE*p->E[i];
		
		gamma2=1.0; for(i=0; i<3; i++) gamma2+=uMinus[i]*uMinus[i];
		invGamma = 1.0/sqrt(gamma2);

		for(i=0; i<3; i++) T[i]=coefB*invGamma*p->B[i];

		T2=0;
		for(i=0; i<3; i++) T2+=T[i]*T[i];
		for(i=0; i<3; i++) S[i]=2.0*T[i]/(1.0+T2);

		M[0][0]=1.0-S[2]*T[2]-S[1]*T[1];
		M[0][1]=S[1]*T[0]+S[2];
		M[0][2]=S[2]*T[0]-S[1];
		M[1][0]=S[1]*T[1]-S[2];
		M[1][1]=1.0-S[0]*T[0]-S[2]*T[2];
		M[1][2]=S[2]*T[1]+S[0];
		M[2][0]=S[0]*T[2]+S[1];
		M[2][1]=S[1]*T[2]-S[0];
		M[2][2]=1.0-S[0]*T[0]-S[1]*T[1];

		for(i=0; i<3; i++) {
      	uPlus[i]=0;
			for(j=0; j<3; j++)
				uPlus[i]+=M[i][j]*uMinus[j];
		}

		for(i=0; i<3; i++) u[i]=uPlus[i]+coefE*p->E[i];

		// update position
		gamma2=1.0; for(i=0; i<3; i++) gamma2+=u[i]*u[i];
		invGamma = 1.0/sqrt(gamma2);

		for(i=0; i<3; i++) { 
        p->uNext[i]=u[i];
        delR=u[i]*invGamma*c*dt;
        p->rNext[i]=p->rNow[i]+delR;
		}

		p=p->next;
	}

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

void restoreDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

void restoreAttr(char *fileName,char *dataName,char *attrName,double *data)
{
  hid_t file_id,dset_id,attr_id;
  hsize_t metaDim[1];
  herr_t status;

  file_id=H5Fopen(fileName,H5F_ACC_RDONLY,H5P_DEFAULT);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  attr_id=H5Aopen(dset_id,attrName,H5P_DEFAULT);

  status=H5Aread(attr_id,H5T_NATIVE_DOUBLE,data);
  H5Aclose(attr_id);
  H5Dclose(dset_id);
  H5Fclose(file_id);
}

void restoreData(char *fileName,char *dataName,double *data)
{
   hid_t file_id,dset_id,plist_id;
   hid_t dataspace;
   hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
   herr_t ierr;

   //open file
	plist_id=H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
	file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
	H5Pclose(plist_id);

	//set dataset
	dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

	dataspace=H5Dget_space(dset_id);	

   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
//   ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,dataspace,H5P_DEFAULT,data);

   H5Sclose(dataspace);
   H5Dclose(dset_id);
   H5Fclose(file_id);
}



static void terminate(const char *message)
{
	printf("%s\n",message);
	exit(EXIT_FAILURE);
}
          
void fresnel(double x,double *C,double *S)
{
	int k,n,odd;
	double a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
	fcomplex b,cc,d,h,del,cs,ONE;
	double EPS,FPMIN=1.0e-30,XMIN=1.5,PI=3.1415927,PIBY2;
	int MAXIT=100,TRUE=1;
	PIBY2=(PI/2.0);
	ONE=Complex(1.0,0.0);
	EPS=6e-8;

	ax=fabs(x);
	if (ax < sqrt(FPMIN)) {
		*S=0.0;
		*C=ax;
	}	else if (ax <= XMIN) { 
		sum=sums=0.0;
		sumc=ax;
		sign=1.0;
		fact=PIBY2*ax*ax;
		odd=TRUE;
		term=ax;
		n=3;
		for (k=1;k<=MAXIT;k++) {
			term *= fact/k;
			sum += sign*term/n;
			test=fabs(sum)*EPS;
			if (odd) {
				sign = -sign;
				sums=sum;
				sum=sumc;
			} else {
				sumc=sum;
				sum=sums;
			}
			if (term < test) break;
			odd=!odd;
			n += 2;
		}
		if (k > MAXIT) { printf("series failed in frenel"); exit(0); } else ;
		*S=sums;
		*C=sumc;
	}	else { 
		pix2=PI*ax*ax;
		b=Complex(1.0,-pix2);
		cc=Complex(1.0/FPMIN,0.0);
		d=h=Cdiv(ONE,b);
		n = -1;
		for (k=2;k<=MAXIT;k++) {
			n += 2;
			a = -n*(n+1);
			b=Cadd(b,Complex(4.0,0.0));
			d=Cdiv(ONE,Cadd(RCmul(a,d),b));
			cc=Cadd(b,Cdiv(Complex(a,0.0),cc));
			del=Cmul(cc,d);
			h=Cmul(h,del);
			if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;			
		}
		if (k > MAXIT) { printf("cf failed in frenel, x=%g",x); exit(0); }
		h=Cmul(Complex(ax,-ax),h);
		cs=Cmul(Complex(0.5,0.5),Csub(ONE,Cmul(Complex(cos(0.5*pix2),sin(0.5*pix2)),h)));
		*C=cs.r;
		*S=cs.i;
	}
	if (x < 0.0) { 
		*C = -(*C);
		*S = -(*S);
	}
}

fcomplex Cadd(fcomplex a, fcomplex b)
{
	fcomplex C;
	C.r=a.r+b.r;
	C.i=a.i+b.i;
	return C;
}
fcomplex Csub(fcomplex a, fcomplex b)
{
	fcomplex C;
	C.r=a.r-b.r;
	C.i=a.i-b.i;
	return C;
}
fcomplex Cmul(fcomplex a, fcomplex b)
{
	fcomplex C;
	C.r=a.r*b.r-a.i*b.i;
	C.i=a.i*b.r+a.r*b.i;
	return C;
}
fcomplex Complex(double re, double im)
{
	fcomplex C;
	C.r=re;
	C.i=im;
	return C;
}

fcomplex Cdiv(fcomplex a, fcomplex b) {
	fcomplex C;
	double r,den;
	if (fabs(b.r) >= fabs(b.i)) {
		r=b.i/b.r;
		den=b.r+r*b.i;
		C.r=(a.r+r*a.i)/den;
		C.i=(a.i-r*a.r)/den;
	} else {
		r=b.r/b.i;
		den=b.i+r*b.r;
		C.r=(a.r*r+a.i)/den;
		C.i=(a.i*r-a.r)/den;
	}
	return C;
}
double Cabs(fcomplex z) {
	double x,y,ans,temp;
	x=fabs(z.r);
	y=fabs(z.i);
	if (x == 0.0)
		ans=y;
		else if (y == 0.0)
			ans=x;
			else if (x > y) {
				temp=y/x;
				ans=x*sqrt(1.0+temp*temp);
			} else {
				temp=x/y;
				ans=y*sqrt(1.0+temp*temp);
			}
			return ans;
}
fcomplex RCmul(double x, fcomplex a) {
	fcomplex C;
	C.r=x*a.r;
	C.i=x*a.i;
	return C;
}


	

double FresnelS(double x,int maxN)
{
	int n;
	double b,deno,prev,sum,tmp;
	
	b=M_PI*0.5; 

	//n = 1
	deno=-b*b*x*x*x*x;
	prev=deno/2.0;
	sum=0.5+1.0/(6.0)*prev;

	for(n=2; n<maxN; n++) {
		prev=prev*deno/(2.0*n)/(2.0*(n-1));
		tmp=1.0/(4*n+2.0)*prev;
		sum+=tmp;
	}

	return sum*b*x*x*x;
}



/*
				z=Th_p;
				fresC_p=0.5+(1.0+0.926*z)/(2.0+1.792*z+3.104*z*z)*sin(M_PI*z*z*0.5)
		  								-1.0/(2.0+4.142*z+3.492*z*z+6.67*z*z*z)*cos(M_PI*z*z*0.5);
				fresS_p=0.5-(1.0+0.926*z)/(2.0+1.792*z+3.104*z*z)*cos(M_PI*z*z*0.5)
		  								-1.0/(2.0+4.142*z+3.492*z*z+6.67*z*z*z)*sin(M_PI*z*z*0.5);
				z=Th_m;
				fresC_m=0.5+(1.0+0.926*z)/(2.0+1.792*z+3.104*z*z)*sin(M_PI*z*z*0.5)
		  								-1.0/(2.0+4.142*z+3.492*z*z+6.67*z*z*z)*cos(M_PI*z*z*0.5);
				fresS_m=0.5-(1.0+0.926*z)/(2.0+1.792*z+3.104*z*z)*cos(M_PI*z*z*0.5)
		  								-1.0/(2.0+4.142*z+3.492*z*z+6.67*z*z*z)*sin(M_PI*z*z*0.5);
*/

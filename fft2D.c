#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/*
 this computes an in-place complex to complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir = 1 gives forward transform
 dir = -1 gives reverse transform
*/

int main(int argc, char *argv[])
{
   FILE *in,*out;
   int cenX,cenY,startI[2],endI[2],startJ[2],endJ[2];
   float norData1,norData2,norData3,norData4,norData5,norData6,fx,fy,lowFreX,lowFreY,rangeX,rangeY,center;
   float cData1,cData2,cData3,cData4,cData5,cData6;
   float *xold,*yold,**oldData1,**oldData2,**oldData3,**oldData4,**oldData5,**oldData6;
   float *xnew,*ynew,**newData1,**newData2,**newData3,**newData4,**newData5,**newData6;
   float w00,w01,w10,w11,xx,yy,x,y;
   float xtmpold,xtmp,tmp,dxOld,dyOld,dxNew,dyNew, xrange,x0,y0,minX;
   float lambda,lambdaDivision,rU,c,dt,period,maxX;
   char name[100],outFile[100];
   int binOrderX,binOrderY,dataNumX,dataNumY,centerMode,initial,rnk,rank,step,nx,ny,dataNum;
   int i,j,ii,jj,column;
   void four1();   
   void fourn();   

   if (argc < 3) {
      printf("2Dfft [file] [minX] [maxX]\n");
      exit(0);
   }

   minX = atof(argv[2]);
   maxX = atof(argv[3]);

     in = fopen(argv[1],"r");
     nx=ny=1;
     fscanf(in,"%g %g %g %g %g %g %g %g",&xtmpold,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp);
     while(fscanf(in,"%g %g %g %g %g %g %g %g",&xtmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp)!=EOF)
     {
     if(xtmp!=xtmpold) { nx++; ny=0; }
       ny++;
       xtmpold = xtmp;
     }
     fclose(in);

   binOrderX=((int)(log(nx)/log(2)))+1;
   binOrderY=((int)(log(ny)/log(2)))+1;
   dataNumX=dataNumY=1;
   for (i=1; i<=binOrderX; i++) {
      dataNumX=dataNumX*2;
   }
   for (i=1; i<=binOrderY; i++) {
      dataNumY=dataNumY*2;
   }

   xold =(float *)malloc((nx+1)*sizeof(float));
   yold =(float *)malloc((ny+1)*sizeof(float));
     oldData1=(float **)malloc((nx+1)*sizeof(float *));
     oldData2=(float **)malloc((nx+1)*sizeof(float *));
     oldData3=(float **)malloc((nx+1)*sizeof(float *));
     oldData4=(float **)malloc((nx+1)*sizeof(float *));
     oldData5=(float **)malloc((nx+1)*sizeof(float *));
     oldData6=(float **)malloc((nx+1)*sizeof(float *));
   for(i=0; i<=nx; i++)
   {
       oldData1[i]=(float *)malloc((ny+1)*sizeof(float ));
       oldData2[i]=(float *)malloc((ny+1)*sizeof(float ));
       oldData3[i]=(float *)malloc((ny+1)*sizeof(float ));
       oldData4[i]=(float *)malloc((ny+1)*sizeof(float ));
       oldData5[i]=(float *)malloc((ny+1)*sizeof(float ));
       oldData6[i]=(float *)malloc((ny+1)*sizeof(float ));
   }
   xnew =(float *)malloc((dataNumX+1)*sizeof(float));
   ynew =(float *)malloc((dataNumY+1)*sizeof(float));
     newData1=(float **)malloc((dataNumX+1)*sizeof(float *));
     newData2=(float **)malloc((dataNumX+1)*sizeof(float *));
     newData3=(float **)malloc((dataNumX+1)*sizeof(float *));
     newData4=(float **)malloc((dataNumX+1)*sizeof(float *));
     newData5=(float **)malloc((dataNumX+1)*sizeof(float *));
     newData6=(float **)malloc((dataNumX+1)*sizeof(float *));
   for(i=0; i<=dataNumX; i++)
   {
       newData1[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
       newData2[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
       newData3[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
       newData4[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
       newData5[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
       newData6[i]=(float *)malloc((dataNumY+1)*sizeof(float ));
   }
   cData3=cData4=cData1=cData2=cData5=cData6=0.0;

   in = fopen(argv[1],"r");
   for (i=1; i<=nx; i++)
     for (j=1; j<=ny; j++)
     {
       fscanf(in,"%g %g %g %g %g %g %g %g",&xold[i],&yold[j],&oldData1[i][j],&oldData2[i][j],&oldData3[i][j],&oldData4[i][j],&oldData5[i][j],&oldData6[i][j]);
       if(xold[i]<minX || xold[i]>maxX)  {
           oldData1[i][j]=0.0;
           oldData2[i][j]=0.0;
           oldData3[i][j]=0.0;
           oldData4[i][j]=0.0;
           oldData5[i][j]=0.0;
           oldData6[i][j]=0.0;
       }
         cData3+=oldData3[i][j];
         cData4+=oldData4[i][j];
         cData1+=oldData1[i][j];
         cData2+=oldData2[i][j];
         cData5+=oldData5[i][j];
         cData6+=oldData6[i][j];
     }
   fclose(in);

   cData1=cData1/(nx*ny);
   cData2=cData2/(nx*ny);
   cData3=cData3/(nx*ny);
   cData4=cData4/(nx*ny);
   cData5=cData5/(nx*ny);
   cData6=cData6/(nx*ny);

   dxOld=(xold[nx]-xold[1])/nx;
   dyOld=(yold[ny]-yold[1])/ny;
   dxNew=(xold[nx]-xold[1])/dataNumX;
   dyNew=(yold[ny]-yold[1])/dataNumY;

   x0=xnew[1]=xold[1];
   y0=ynew[1]=yold[1];
   xnew[dataNumX]=xold[nx];
   ynew[dataNumY]=yold[ny];

   for (i=2; i<=dataNumX; i++)	
     xnew[i]=x0+(i-1)*dxNew;
   for (j=2; j<=dataNumY; j++)	
     ynew[j]=y0+(j-1)*dyNew;


   for (i=1; i<dataNumX; i++) 
     for (j=1; j<dataNumY; j++)
     {
       x=(i-1)*dxNew;
       y=(j-1)*dyNew;
       ii=(int)(x/dxOld);
       jj=(int)(y/dyOld);
       xx=x/dxOld-ii;
       yy=y/dyOld-jj;
       w00=xx*yy;
       w10=(1-xx)*yy;
       w01=xx*(1-yy);
       w11=(1-xx)*(1-yy);
     
         newData1[i][j] = w00*oldData1[ii+2][jj+2]+w10*oldData1[ii+1][jj+2]+w01*oldData1[ii+2][jj+1]+w11*oldData1[ii+1][jj+1];
         newData2[i][j] = w00*oldData2[ii+2][jj+2]+w10*oldData2[ii+1][jj+2]+w01*oldData2[ii+2][jj+1]+w11*oldData2[ii+1][jj+1];
         newData3[i][j] = w00*oldData3[ii+2][jj+2]+w10*oldData3[ii+1][jj+2]+w01*oldData3[ii+2][jj+1]+w11*oldData3[ii+1][jj+1];
         newData4[i][j] = w00*oldData4[ii+2][jj+2]+w10*oldData4[ii+1][jj+2]+w01*oldData4[ii+2][jj+1]+w11*oldData4[ii+1][jj+1];
         newData5[i][j] = w00*oldData5[ii+2][jj+2]+w10*oldData5[ii+1][jj+2]+w01*oldData5[ii+2][jj+1]+w11*oldData5[ii+1][jj+1];
         newData6[i][j] = w00*oldData6[ii+2][jj+2]+w10*oldData6[ii+1][jj+2]+w01*oldData6[ii+2][jj+1]+w11*oldData6[ii+1][jj+1];
     }
/*
   for(i=0; i<nx; i++)
   {
     for(j=0; j<ny; j++)
     {
       printf("%g %g %g %g %g %g %g %g\n",xold[i],yold[j],oldData1[i][j],oldData2[i][j],oldData3[i][j],oldData4[i][j],oldData5[i][j],oldData6[i][j]);
     }
     printf("\n");
   }
       
   for(i=0; i<dataNumX; i++)
   {
     for(j=0; j<dataNumY; j++)
     {
       printf("%g %g %g %g %g %g %g %g\n",xnew[i],ynew[j],newData1[i][j],newData2[i][j],newData3[i][j],newData4[i][j],newData5[i][j],newData6[i][j]);
     }
     printf("\n");
   }
*/
   dataNum=dataNumX*2*dataNumY;
   float *resultX,*Data1,*Data2,*Data3,*Data4,*Data5,*Data6;
   resultX=(float *)malloc((dataNum+1)*sizeof(float));
     Data1=(float *)malloc((dataNum+1)*sizeof(float ));
     Data2=(float *)malloc((dataNum+1)*sizeof(float ));
     Data3=(float *)malloc((dataNum+1)*sizeof(float ));
     Data4=(float *)malloc((dataNum+1)*sizeof(float ));
     Data5=(float *)malloc((dataNum+1)*sizeof(float ));
     Data6=(float *)malloc((dataNum+1)*sizeof(float ));

   for (j=1; j<=dataNumY ; j++) 
   {
     for (i=1; i<=dataNumX ; i++) 
     {
       resultX[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data1[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data2[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data3[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data4[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data5[(i-1)*2+(j-1)*dataNumX*2]=0.0;
         Data6[(i-1)*2+(j-1)*dataNumX*2]=0.0;
       resultX[(i-1)*2+1+(j-1)*dataNumX*2]=xnew[i];
         Data1[(i-1)*2+1+(j-1)*dataNumX*2]=newData1[i][j];
         Data2[(i-1)*2+1+(j-1)*dataNumX*2]=newData2[i][j];
         Data3[(i-1)*2+1+(j-1)*dataNumX*2]=newData3[i][j];
         Data4[(i-1)*2+1+(j-1)*dataNumX*2]=newData4[i][j];
         Data5[(i-1)*2+1+(j-1)*dataNumX*2]=newData5[i][j];
         Data6[(i-1)*2+1+(j-1)*dataNumX*2]=newData6[i][j];
     }
   }   

   int n[3];
//   n=(int *)malloc(3*sizeof(int));
   n[0]=0;
   n[2]=dataNumX;
   n[1]=dataNumY;

     fourn(Data1, n, 2, 1);
     fourn(Data2, n, 2, 1);
     fourn(Data3, n, 2, 1);
     fourn(Data4, n, 2, 1);
     fourn(Data5, n, 2, 1);
     fourn(Data6, n, 2, 1);

   sprintf(outFile,"fft%s",argv[1]);
   out=fopen(outFile,"w");
   cenX=dataNumX/2; cenY=dataNumY/2;
   lowFreX=-cenX;   lowFreY=-cenY;
   rangeX=xold[nx]-xold[1];
   rangeY=yold[ny]-yold[1];

   startI[0]=cenX; endI[0]=dataNumX;
   startI[1]=1;    endI[1]=cenX;
   startJ[0]=cenY; endJ[0]=dataNumY;
   startJ[1]=1;    endJ[1]=cenY;
   for(ii=0; ii<2; ii++)
   {
     for(i=startI[ii]; i<endI[ii]; i++) 
     {
       fx=lowFreX/rangeX;
       for(jj=0; jj<2; jj++)
       {
         for(j=startJ[jj]; j<endJ[jj]; j++) 
         {
           fy=lowFreY/rangeY;
           norData1=(Data1[(i-1)*2+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+(j-1)*dataNumX*2]+Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]);
           norData2=(Data2[(i-1)*2+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+(j-1)*dataNumX*2]+Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]);
           norData3=(Data3[(i-1)*2+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+(j-1)*dataNumX*2]+Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]);
           norData4=(Data4[(i-1)*2+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+(j-1)*dataNumX*2]+Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]);
           norData5=sqrt(Data5[(i-1)*2+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+(j-1)*dataNumX*2]+Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]);
           norData6=sqrt(Data6[(i-1)*2+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+(j-1)*dataNumX*2]+Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]);

           fprintf(out,"%g %g %g %g %g %g %g %g\n",fx,fy,norData1,norData2,norData3,norData4,norData5,norData6);
           lowFreY=lowFreY+1.0;
         }
       }
       fprintf(out,"\n");
       lowFreY=-cenY;
       lowFreX=lowFreX+1.0;
     }
   }

/*
   for (i=1; i<cenX; i++) 
   {
     fx=lowFreX/rangeX;
     for (j=cenY; j<dataNumY; j++) 
     {
       fy=lowFreY/rangeY;
         norData1=(Data1[(i-1)*2+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+(j-1)*dataNumX*2]+Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData2=(Data2[(i-1)*2+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+(j-1)*dataNumX*2]+Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData3=(Data3[(i-1)*2+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+(j-1)*dataNumX*2]+Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData4=(Data4[(i-1)*2+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+(j-1)*dataNumX*2]+Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData5=sqrt(Data5[(i-1)*2+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+(j-1)*dataNumX*2]+Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData6=sqrt(Data6[(i-1)*2+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+(j-1)*dataNumX*2]+Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]);

       lowFreY=lowFreY+1.0;
//         printf("%g %g %g %g %g %g %g %g\n",fx,fy,norData1,norData2,norData3,norData4,norData5,norData6);
       fprintf(out,"%g %g %g %g %g %g %g %g\n",fx,fy,norData1,norData2,norData3,norData4,norData5,norData6);
     }
     for (j=1; j<cenY; j++) 
     {
       fy=lowFreY/rangeY;
         norData1=(Data1[(i-1)*2+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+(j-1)*dataNumX*2]+Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data1[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData2=(Data2[(i-1)*2+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+(j-1)*dataNumX*2]+Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data2[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData3=(Data3[(i-1)*2+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+(j-1)*dataNumX*2]+Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data3[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData4=(Data4[(i-1)*2+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+(j-1)*dataNumX*2]+Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data4[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData5=sqrt(Data5[(i-1)*2+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+(j-1)*dataNumX*2]+Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data5[(i-1)*2+1+1+(j-1)*dataNumX*2]);
         norData6=sqrt(Data6[(i-1)*2+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+(j-1)*dataNumX*2]+Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]*Data6[(i-1)*2+1+1+(j-1)*dataNumX*2]);

       lowFreY=lowFreY+1.0;
//         printf("%g %g %g %g %g %g %g %g\n",fx,fy,norData1,norData2,norData3,norData4,norData5,norData6);
       fprintf(out,"%g %g %g %g %g %g %g %g\n",fx,fy,norData1,norData2,norData3,norData4,norData5,norData6);
     }
     lowFreY=-cenY;
     lowFreX=lowFreX+1.0;
     fprintf(out,"\n");
   }
*/
   fclose(out);
   printf("%s is made.\n",outFile);

   free(xold);
   free(yold);
   for(i=0; i<nx; i++)
   {  
     free(oldData1[i]);
     free(oldData2[i]);
     free(oldData3[i]);
     free(oldData4[i]);
     free(oldData5[i]);
     free(oldData6[i]);
   }

   free(oldData1);
   free(oldData2);
   free(oldData3);
   free(oldData4);
   free(oldData5);
   free(oldData6);
   free(xnew);
   free(ynew);
   for(i=0; i<dataNumX; i++)
   {
     free(newData1[i]);
     free(newData2[i]);
     free(newData3[i]);
     free(newData4[i]);
     free(newData5[i]);
     free(newData6[i]);
   }
   free(newData1);
   free(newData2);
   free(newData3);
   free(newData4);
   free(newData5);
   free(newData6);
   free(resultX);
   free(Data1);
   free(Data2);
   free(Data3);
   free(Data4);
   free(Data5);
   free(Data6);
}

void fourn(float *data, int *nn, int ndim, int isign)
{
   int idim;
   unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
   unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
   float tempi,tempr;
   double theta,wi,wpi,wpr,wr,wtemp;

   for(ntot=1,idim=1;idim<=ndim;idim++)
     ntot*=nn[idim];
   nprev=1;
   for(idim=ndim;idim>=1;idim--)
   {
     n=nn[idim];
     nrem=ntot/(n*nprev);
     ip1=nprev<<1;
     ip2=ip1*n;
     ip3=ip2*nrem;
     i2rev=1;
     for(i2=1;i2<=ip2;i2+=ip1)
     {
       if(i2<i2rev)
       {
         for(i1=i2;i1<=i2+ip1-2;i1+=2)
	 {
           for(i3=i1;i3<=ip3;i3+=ip2)
	   {
	     i3rev=i2rev+i3-i2;
	     SWAP(data[i3],data[i3rev]);
	     SWAP(data[i3+1],data[i3rev+1]);
	   }
	 }
       }
       ibit=ip2>>1;
       while(ibit>=ip1 && i2rev>ibit)
       {
         i2rev -= ibit;
         ibit >>= 1;
       }
       i2rev += ibit;
     }
   
     ifp1=ip1;
     while(ifp1<ip2)
     {
//printf("ifp1=%ld, 1p2=%ld\n",ifp1,ip2);
       ifp2=ifp1<<1;
       theta=isign*6.28318530717959/(ifp2/ip1);
       wtemp=sin(0.5*theta);
       wpr=-2.0*wtemp*wtemp;
       wpi=sin(theta);
       wr=1.0;
       wi=0.0;
       for(i3=1;i3<=ifp1;i3+=ip1)
       {
         for(i1=i3;i1<=i3+ip1-2;i1+=2)
	 {
	   for(i2=i1;i2<=ip3;i2+=ifp2)
	   {
	     k1=i2;
	     k2=k1+ifp1;
	     tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
	     tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
	     data[k2]=data[k1]-tempr;
	     data[k2+1]=data[k1+1]-tempi;
	     data[k1]+=tempr;
	     data[k1+1]+=tempi;
	   }
	 }
	 wr=(wtemp=wr)*wpr-wi*wpi+wr;
	 wi=wi*wpr+wtemp*wpi+wi;
       }
       ifp1=ifp2;
     }
     nprev*=n;
   }
}

void four1(float *data, int n, int isign)
{
  int nn, mmax,m,j,istep,i;
  float wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

 /* This is the bit-reversal section of the routine.
    Exchange the two complex numbers. */
  nn = n << 1;
  j = 1;
  for (i=1; i<nn; i+=2) {
    if (j > i) {
      SWAP(data[j-1], data[i-1]);
      SWAP(data[j], data[i]);
    }
    m = n;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }  
    j += m; 
  }

/* Here begins the Danielson-Lanczos section of the routine.*/
  mmax = 2;
  while (nn > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1; m<mmax; m+=2) {
      for (i=m; i<=nn; i+=istep) {
        j=i+mmax;
        tempr = wr*data[j-1]-wi*data[j];
        tempi = wr*data[j] + wi*data[j-1];
        data[j-1] = data[i-1] -tempr;
        data[j] = data[i] - tempi;
        data[i-1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr+wtemp*wpi+wi;
    }
    mmax = istep;
  }
} 
        

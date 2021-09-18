// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "malloc.h"
#include "math.h"

float whatMass(int atomNum);

void main(int argc, char *argv[])
{
   int i,j,l,n,numS,s,mode,step,initial,final,saveStep,core;
   int cnt,atomNum,cores,ii,jj,numX,numY,numE,species,steps,cntStep;
   float minE,maxE,dE,slope,rangeX,dX,dY,wE,wX,**nE,mass,minP,tmp;
   float minX,maxX,minY,maxY,minZ,maxZ,slopeX,dx,dy,dz,wx[2],wy[2];
   float minPx,maxPx,shiftX;
   float x,y,z,px,py,pz,id,weight,energy,rangeE,den,emitY,emitZ;
   char fileName[100],outFile[100],**name,dataType[100];
   float **density,**dataX,**dataY,**dataZ,**dataPx;
   float *dataID,*trackID;
   int *dataCore,*trackCore,**dataT;
   FILE *out,*in;

   if(argc < 5)
   {
      printf("particle mode initial final saveStep\n");
      printf("mode(1) : dataType minE maxE numE numS atomNum spc1 spc2 ...\n");
      printf("mode(2) : minE maxE numE minX maxX numX shiftX numS spc1 spc2 ...\n");
      printf("mode(3) : minE maxE numE slope rangeX numX shiftX numS spc1 spc2 ...\n");
      printf("mode(4) : minX maxX numX minY maxY numY minZ maxZ species\n");
      printf("mode(5) : slopeX rangeX numX minY maxY numY minZ maxZ species\n");
      printf("mode(6) : inFile species\n");
      printf("mode(7) : species  /*with 'track' file*/\n");
      printf("mode(8) : inFile minX maxX minY maxY minZ maxZ minPx maxPx\n");
      printf("mode(9) : inFile species /*extracting id filed*/\n");
      printf("mode(10) : dataType minEmitY maxEmitY minEmitZ maxEmitZ numEmit\n");
      printf("mode(11) : dataType minX maxX numX minY maxY numY numS spc1 spc2...\n");
      exit(0);
   }
   mode=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   saveStep=atoi(argv[4]);

   switch (mode) {
   case 11 :
     sprintf(dataType,"%s",argv[5]);
     minX=atof(argv[6]);
     maxX=atof(argv[7]);
     numX=atoi(argv[8]);
     minY=atof(argv[9]);
     maxY=atof(argv[10]);
     numY=atoi(argv[11]);
     numS=atoi(argv[12]);

     name=(char **)malloc(numS*sizeof(char *));
     for(s=0; s<numS; s++)
       name[s]=(char *)malloc(100*sizeof(char ));

     for(s=0; s<numS; s++) 
       sprintf(name[s],"%s",argv[s+13]);

     nE=(float **)malloc((numX+1)*sizeof(float *) );
     for(i=0; i<=numX; i++)
       nE[i]=(float *)malloc((numY+1)*sizeof(float ) );

     dX=(maxX-minX)/((float)(numX));
     dY=(maxY-minY)/((float)(numY));

     mass=1.0;
     for(step=initial; step<=final; step+=saveStep)
     {
       for(i=0; i<=numX; i++)
         for(j=0; j<=numE; j++)
           nE[i][j]=0.0;

       for(s=0; s<numS; s++)       {
         sprintf(fileName,"%s%s%d",name[s],dataType,step);
         in=fopen(fileName,"r");
         while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&px,&py,&pz,&tmp,&tmp,&weight)!=EOF)  {
           i=(int)((x-minX)/dX);
           wx[1]=(x-minX)/dX-i; wx[0]=1.0-wx[1];
           j=(int)((y-minY)/dY);
           wy[1]=(y-minY)/dY-j; wy[0]=1.0-wy[1];
           if(i>=0 && i<numX && j>=0 && j<numY) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++)
                 nE[i+ii][j+jj]+=wx[ii]*wy[jj]*weight;
           }
         }
         fclose(in);
         printf("%s is done.\n",fileName);
       }

       sprintf(fileName,"denXY%d",step);
       out=fopen(fileName,"w");
       for(i=0; i<=numX; i++) {
         for(j=0; j<=numY; j++)
           fprintf(out,"%g %g %g\n",minX+i*dX,minY+j*dY,fabs(nE[i][j]/dX/dE));
         fprintf(out,"\n");
       }
       printf("%s is made\n",fileName);
     }

     for(i=0; i<numX+1; i++) free(nE[i]); free(nE);
     free(name);
    
     break;

   case 10 :
     sprintf(dataType,"%s",argv[5]);
     minX=atof(argv[6]);
     maxX=atof(argv[7]);
     minY=atof(argv[8]);
     maxY=atof(argv[9]);
     numX=atoi(argv[10]);
     numY=numX;

     nE=(float **)malloc((numX+1)*sizeof(float *) );
     for(i=0; i<=numX; i++)
       nE[i]=(float *)malloc((numY+1)*sizeof(float ) );

     dX=(maxX-minX)/((float)(numX));
     dY=(maxY-minY)/((float)(numY));

     for(step=initial; step<=final; step+=saveStep)
     {
       for(i=0; i<=numX; i++)
         for(j=0; j<=numY; j++)
           nE[i][j]=0.0;

         sprintf(fileName,"%s%d",dataType,step);
         in=fopen(fileName,"r");
         while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&px,&py,&pz,&tmp,&tmp,&weight)!=EOF)  {
           emitY=py/px;
           emitZ=pz/px;
           i=(int)((emitY-minX)/dX);
           wx[1]=(emitY-minX)/dX-i; wx[0]=1.0-wx[1];
           j=(int)((emitZ-minY)/dY);
           wy[1]=(emitZ-minY)/dY-j; wy[0]=1.0-wy[1];
           if(i>=0 && i<numX && j>=0 && j<numY) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++)
                 nE[i+ii][j+jj]+=wx[ii]*wy[jj]*weight;
           }
         }
         fclose(in);
         printf("%s is done.\n",fileName);

       sprintf(fileName,"emit%d",step);
       out=fopen(fileName,"w");
       for(i=0; i<=numX; i++) {
         for(j=0; j<=numY; j++)
           fprintf(out,"%g %g %g\n",minX+i*dX,minY+j*dY,fabs(nE[i][j]/dX/dY));
         fprintf(out,"\n");
       }
       printf("%s is made\n",fileName);
     }

     for(i=0; i<numX+1; i++) free(nE[i]); free(nE);
     break;

   case 9 :
     species=atoi(argv[6]);

     sprintf(fileName,"%s",argv[5]);
     in=fopen(fileName,"r");
     cnt=0;
     while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  
       cnt++;
     fclose(in);

     dataID=(float *)malloc(cnt*sizeof(float ) );
     dataCore=(int *)malloc(cnt*sizeof(int ) );
     steps=(final-initial)/saveStep+1;
     dataX=(float **)malloc(steps*sizeof(float *) );
     dataY=(float **)malloc(steps*sizeof(float *) );
     dataZ=(float **)malloc(steps*sizeof(float *) );
     dataPx=(float **)malloc(steps*sizeof(float *) );
     dataT=(int **)malloc(steps*sizeof(int *) );
     for(i=0; i<steps; i++) {
       dataX[i]=(float *)malloc(cnt*sizeof(float ) );
       dataY[i]=(float *)malloc(cnt*sizeof(float ) );
       dataZ[i]=(float *)malloc(cnt*sizeof(float ) );
       dataPx[i]=(float *)malloc(cnt*sizeof(float ) );
       dataT[i]=(int *)malloc(cnt*sizeof(int ) );
     }     
     
     in=fopen(fileName,"r");
     for(i=0; i<cnt; i++)
       fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&dataID[i],&dataCore[i],&den);
     fclose(in);

     for(step=initial; step<=final; step+=saveStep)
     {
       n=(step-initial)/saveStep;

       sprintf(fileName,"%dParticle%d",species,step);
       in=fopen(fileName,"r");

       while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  {
         for(i=0; i<cnt; i++) 
           if(dataID[i]==id && dataCore[i]==core) {
             dataX[n][i]=x;
             dataY[n][i]=y;
             dataZ[n][i]=z;
             dataPx[n][i]=px;
             dataT[n][i]=initial+n*saveStep;
           }
       }
       fclose(in);
       printf("%s is done.\n",fileName);
     }

     for(i=0; i<cnt; i++) {
       id=dataID[i];
       sprintf(outFile,"%d_%g",species,id);
       out=fopen(outFile,"w");
       for(n=0; n<steps; n++) {
         x=dataX[n][i];
         y=dataY[n][i];
         z=dataZ[n][i];
         px=dataPx[n][i];
         step=dataT[n][i];
         fprintf(out,"%g %g %g %g %d\n",x,y,z,px,step);
       }
       fclose(out);
       printf("%s is made.\n",outFile);
     }
  
     free(dataID); free(dataCore);
     for(n=0; n<steps; n++) {
       free(dataX[n]); free(dataY[n]); free(dataZ[n]);
       free(dataPx[n]); free(dataT[n]);
     }
     free(dataX); free(dataY); free(dataZ); free(dataPx); free(dataT);
     break;

   case 7 :
     species=atoi(argv[5]);

     sprintf(fileName,"track");
     in=fopen(fileName,"r");
     cnt=0;
     while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  
       cnt++;
     fclose(in);

     trackID=(float *)malloc(cnt*sizeof(float ) );
     trackCore=(int *)malloc(cnt*sizeof(int ) );
     in=fopen(fileName,"r");
     for(i=0; i<cnt; i++)
       fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&trackID[i],&trackCore[i],&den);
     fclose(in);

     for(i=0; i<cnt; i++) 
     {
       sprintf(outFile,"track%g_%d",trackID[i],trackCore[i]);
       out=fopen(outFile,"w");

       for(step=initial; step<=final; step+=saveStep) {
         sprintf(fileName,"%dParticle%d",species,step);
         in=fopen(fileName,"r");

         while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF) {
           if(trackID[i]==id && trackCore[i]==core) 
             fprintf(out,"%g %g %g %.12g %g %g %d\n",x,y,z,px,py,pz,step);
           else ;
         }
         fclose(in);
         printf("For %s, %s is done\n",outFile,fileName);
       }
       fclose(out);
       printf("%s is made\n",outFile);
     }

    
     free(trackID); free(trackCore);
     break;

   case 6 :
     species=atoi(argv[6]);

     sprintf(fileName,"%s",argv[5]);
     in=fopen(fileName,"r");
     cnt=0;
     while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  
       cnt++;
     fclose(in);

     dataID=(float *)malloc(cnt*sizeof(float ) );
     dataCore=(int *)malloc(cnt*sizeof(int ) );
     
     in=fopen(fileName,"r");
     for(i=0; i<cnt; i++)
       fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&dataID[i],&dataCore[i],&den);
     fclose(in);

     for(step=initial; step<=final; step+=saveStep)
     {
       sprintf(outFile,"%did%d",species,step);

       out=fopen(outFile,"w");

       sprintf(fileName,"%dParticle%d",species,step);
       in=fopen(fileName,"r");

       while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  {
         for(i=0; i<cnt; i++) 
           if(dataID[i]==id && dataCore[i]==core) 
             fprintf(out,"%g %g %g %g %g %g %.15g %d %g\n",x,y,z,px,py,pz,id,core,den);
       }

       fclose(in);

       fclose(out);

       printf("%s is made\n",outFile);
     }
  
     free(dataID); free(dataCore);
     break;

   case 4 :
   case 5 :
     if(mode==4) {
       minX=atof(argv[5]);
       maxX=atof(argv[6]);
     } else {
       slope=atof(argv[5]);
       rangeX=atof(argv[6]);
     }
     numX=atoi(argv[7]);
     minY=atof(argv[8]);
     maxY=atof(argv[9]);
     numY=atoi(argv[10]);
     minZ=atof(argv[11]);
     maxZ=atof(argv[12]);
     species=atoi(argv[13]);

     nE=(float **)malloc((numX+1)*sizeof(float *) );
     for(i=0; i<numX+1; i++) 
       nE[i]=(float *)malloc((numY+1)*sizeof(float ) );
     
     dy=(maxY-minY)/(float)(numY);
     dz=maxZ-minZ; if(dz==0.0) dz=1.0; else ;
     if(mode==4) dx=(maxX-minX)/(float)(numX);
     else        dx=rangeX/(float)numX;
     
     for(step=initial; step<=final; step+=saveStep)
     {
       if(mode==5) { maxX=slope*step; minX=maxX-rangeX; } else ;
       for(i=0; i<numX+1; i++) 
         for(j=0; j<numY+1; j++) 
           nE[i][j]=0.0;
       
       sprintf(fileName,"%dParticle%d",species,step);
       in=fopen(fileName,"r");
       while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&weight)!=EOF)  {
           i=(int)((x-minX)/dx);
           j=(int)((y-minY)/dy);
//printf("i=%d,minX=%g,x=%g,dx=%g\n",i,minX,x,dx);
//printf("j=%d,minY=%g,y=%g,dy=%g\n",j,minY,y,dy);
           wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
           wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
           if(i>=0 && i<numX && j>=0 && j<numY && z>=minZ && z<=maxZ) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 nE[i+ii][j+jj]+=wx[ii]*wy[jj]/dx/dy/dz;
//                 printf("nE=%g\n",nE[i+ii][j+jj]);
               }
           } else ;
       }
       fclose(in);
       printf("%s is done.\n",fileName);
       
       sprintf(fileName,"%dDenXY%d",species,step);
       out=fopen(fileName,"w");
       for(i=0; i<numX; i++)  {
         for(j=0; j<numY; j++) {
           x=minX+i*dx;
           y=minY+j*dy;
           fprintf(out,"%g %g %g\n",x,y,nE[i][j]);
         }
         fprintf(out,"\n");
       }
       fclose(out);
       printf("%s is made\n",fileName);
     }
     
     for(i=0; i<numX+1; i++) free(nE[i]);
     free(nE);
/*
     slope=atof(argv[5]);
     rangeX=atof(argv[6]);
     numX=atoi(argv[7]);
     minY=atof(argv[8]);
     maxY=atof(argv[9]);
     numY=atoi(argv[10]);
     minZ=atof(argv[11]);
     maxZ=atof(argv[12]);
     cores=atof(argv[13]);

     nE=(float **)malloc((numX+1)*sizeof(float *) );
     for(i=0; i<numX+1; i++) 
       nE[i]=(float *)malloc((numY+1)*sizeof(float ) );
     
     dy=(maxY-minY)/(float)(numY);
     dz=maxZ-minZ;
     for(step=initial; step<=final; step+=saveStep)
     {
       maxX=slope*step; minX=maxX-rangeX;
       dx=(maxX-minX)/(float)(numX);
       for(i=0; i<numX+1; i++) 
         for(j=0; j<numY+1; j++) 
           nE[i][j]=0.0;
       
       for(n=0; n<cores; n++) {
         sprintf(fileName,"0Particle%d_%d",step,n);
         in=fopen(fileName,"r");
         while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&y,&z,&px,&py,&pz,&id,&core)!=EOF)  {
           i=(int)((x-minX)/dx);
           j=(int)((y-minY)/dy);
           wx[1]=(x-minX)/dx-i; wx[0]=1.0-wx[1];
           wy[1]=(y-minY)/dy-j; wy[0]=1.0-wy[1];
           if(i>=0 && i<numX && j>=0 && j<numY && z>=minZ && z<maxZ) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++)
                 nE[i+ii][j+jj]+=wx[ii]*wy[jj]/dx/dy/dz;
           } else ;
         }
         fclose(in);
         printf("%s is done.\n",fileName);
       }
       
       sprintf(fileName,"0DenXY%d",step);
       out=fopen(fileName,"w");
       for(i=0; i<numX; i++)  {
         for(j=0; j<numY; j++) {
           x=minX+i*dx;
           y=minY+j*dy;
           fprintf(out,"%g %g %g\n",x,y,nE[i][j]);
         }
         fprintf(out,"\n");
       }
       fclose(out);
       printf("%s is made\n",fileName);
     }
     for(i=0; i<numX+1; i++) free(nE[i]);
     free(nE);
*/
     break;

   case 8 :
     minX=atof(argv[6]);
     maxX=atof(argv[7]);
     minY=atof(argv[8]);
     maxY=atof(argv[9]);
     minZ=atof(argv[10]);
     maxZ=atof(argv[11]);
     minPx=atof(argv[12]);
     maxPx=atof(argv[13]);
//lala
     sprintf(fileName,"%s",argv[5]);
     in=fopen(fileName,"r");
     sprintf(outFile,"re%s",argv[5]);
     out=fopen(outFile,"w");
     while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&den)!=EOF)  {
       if(x>minX && x<maxX && y>minY && y<maxY && z>minZ && z<maxZ && px>minPx && px<maxPx) {
         fprintf(out,"%g %g %g %g %g %g %.15g %d %g\n",x,y,z,px,py,pz,id,core,den);
       } else ;
     }
     fclose(in);
     fclose(out);
     printf("%s is made\n",outFile);
     break;

   case 1 :
     sprintf(dataType,"%s",argv[5]);
     minE=atof(argv[6]);
     maxE=atof(argv[7]);
     numE=atoi(argv[8]);
     numS=atoi(argv[9]);
     atomNum=atoi(argv[10]);
     mass=whatMass(atomNum);

     name=(char **)malloc(numS*sizeof(char *));
     for(s=0; s<numS; s++)
       name[s]=(char *)malloc(100*sizeof(char ));

     for(s=0; s<numS; s++) 
       sprintf(name[s],"%s",argv[s+11]);

     nE=(float **)malloc(1*sizeof(float *) );
     for(i=0; i<1; i++)
       nE[i]=(float *)malloc((numE+1)*sizeof(float ) );
     dE=(maxE-minE)/((float)(numE));
//lala
     
     for(step=initial; step<=final; step+=saveStep)
     {
       for(i=0; i<numE+1; i++) nE[0][i]=0.0;
       
       for(s=0; s<numS; s++)  {
         sprintf(fileName,"%s%s%d",name[s],dataType,step);
         in=fopen(fileName,"r");
         while(fscanf(in,"%g %g %g %g %g %g %g %d %g",&x,&y,&z,&px,&py,&pz,&id,&core,&weight)!=EOF)  {
           energy=0.511*mass*(sqrt(1.0+px*px+py*py+pz*pz)-1.0);
           i=(int)((energy-minE)/dE);
           wE=(energy-minE)/dE-i;
           if(i>=0 && i<numE) {
             nE[0][i]+=(1.0-wE)*weight;
             nE[0][i+1]+=wE*weight;
           } else ;
         }
         fclose(in);
       }
       sprintf(fileName,"spectrum%d",step);
       out=fopen(fileName,"w");
       for(i=0; i<numE; i++) 
         fprintf(out,"%g %g\n",minE+i*dE,fabs(nE[0][i]/dE));
       fclose(out);
       printf("%s is made\n",fileName);
     }

     for(i=0; i<1; i++) free(nE[i]);
     free(nE);
     for(s=0; s<numS; s++) free(name[s]);
     free(name);
     break;

   case 2 :
   case 3 :
     minE=atof(argv[5]);
     maxE=atof(argv[6]);
     numE=atoi(argv[7]);
     if(mode==2) {
       minX=atof(argv[8]);
       maxX=atof(argv[9]);
     } else {
       slope=atof(argv[8]);
       rangeX=atof(argv[9]);
     }
     numX=atoi(argv[10]);
     shiftX=atof(argv[11]);
     numS=atoi(argv[12]);

     name=(char **)malloc(numS*sizeof(char *));
     for(s=0; s<numS; s++)
       name[s]=(char *)malloc(100*sizeof(char ));

     for(s=0; s<numS; s++) 
       sprintf(name[s],"%s",argv[s+13]);

     nE=(float **)malloc((numX+1)*sizeof(float *) );
     for(i=0; i<=numX; i++)
       nE[i]=(float *)malloc((numE+1)*sizeof(float ) );

     dE=(maxE-minE)/((float)(numE));
     if(mode==2) dX=(maxX-minX)/((float)(numX));
     else        dX=rangeX/((float)(numX));
//lala
     mass=1.0;
     for(step=initial; step<=final; step+=saveStep)
     {
       if(mode==3) {
         maxX=step*slope+shiftX;
         minX=maxX-rangeX;
       } else ;
       for(i=0; i<=numX; i++)
         for(j=0; j<=numE; j++)
           nE[i][j]=0.0;

       for(s=0; s<numS; s++)       {
         sprintf(fileName,"%sParticle%d",name[s],step);
         in=fopen(fileName,"r");
         while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&px,&py,&pz,&tmp,&tmp,&weight)!=EOF)  {
           energy=0.511*mass*(sqrt(1.0+px*px+py*py+pz*pz)-1.0);
           i=(int)((x-minX)/dX);
           wx[1]=(x-minX)/dX-i; wx[0]=1.0-wx[1];
           j=(int)((energy-minE)/dE);
           wy[1]=(energy-minE)/dE-j; wy[0]=1.0-wy[1];
           if(i>=0 && i<numX && j>=0 && j<numE) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++)
                 nE[i+ii][j+jj]+=wx[ii]*wy[jj]*weight;
           }
         }
         fclose(in);
         printf("%s is done.\n",fileName);
       }

       sprintf(fileName,"denSpec%d",step);
       out=fopen(fileName,"w");
       for(i=0; i<=numX; i++) {
         for(j=0; j<=numE; j++)
           fprintf(out,"%g %g %g\n",minX+i*dX,minE+j*dE,fabs(nE[i][j]/dX/dE));
         fprintf(out,"\n");
       }
       printf("%s is made\n",fileName);
     }

     for(i=0; i<numX+1; i++) free(nE[i]); free(nE);
     free(name);
    
     break;
   }

}

float whatMass(int atomNum)
{
  float eMassU;

  eMassU=5.485799e-4;
  if(atomNum==0)        return 1.0;
  else if(atomNum==1)   return 1.0/eMassU;
  else if(atomNum==2)   return 4.0/eMassU;
  else if(atomNum==6)   return 12.0/eMassU;
  else { printf("no list in atomNum data\n"); exit(0); }
}


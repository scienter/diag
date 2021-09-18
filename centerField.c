#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int main(int argc, char *argv[])
{
   double x,y,Ex,Ey,Ez,Bx,By,Bz,den,D,F,dx,dt,minX,weight,rangeX,oldEx,shiftX;
   double x1,x2,val,pickX;
   double *dataX,*dataT,***dataF,**dataD,*data1,*data2,*data3;
   FILE *in,*out;
   char fileName[100],fileName1[100],fileName2[100],outFile[100],dataType[100],comp[100];
   int dimension,mode,core,initial,final,saveStep;
   int cnt,cntT,i,j,k,rangeI,step,index,cenFlag;

   if(argc < 6)   { 
      printf("centerfield mode dimension initial final saveStep dataType\n");
      printf("mode 0(x-T) : rangeX shiftX dt dx core \n");
      printf("mode 1(x-T,cylind) : rangeX shiftX dt dx core \n");
      printf("mode 2(pickX) : component(Ex,Ey,0,1..) pickX dx \n");
      exit(0);
   }

   mode = atoi(argv[1]);
   dimension = atoi(argv[2]);
   initial = atoi(argv[3]);
   final = atoi(argv[4]);
   saveStep = atoi(argv[5]);
   cntT=(final-initial)/saveStep+1;

   sprintf(dataType,"%s",argv[6]);

   switch (mode) {
   case 0 :
   case 1 :
     rangeX = atof(argv[7]);
     shiftX = atof(argv[8]);
     dt = atof(argv[9]);
     dx = atof(argv[10]);
     core = atoi(argv[11]);

     step=initial;
     sprintf(fileName,"%s%d_%d",dataType,step,core); 

     if(fopen(fileName,"r")==NULL) {
       printf("%s is not exited.\n",fileName);
       cnt=1;
       step=final+saveStep;
     } else {
       in = fopen(fileName,"r");
       cnt=0;
       if(dimension==1) { 
         while(fscanf(in,"%lf %lf",&x,&Ex)!=EOF) cnt++;
       } else {
         while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF) 
           cnt++;
       }
       fclose(in);
       step=initial;
     }

     minX=-rangeX+shiftX;
     rangeI=(int)(rangeX/dx)+1;

     dataX=(double *)malloc(cnt*sizeof(double));
     data1=(double *)malloc(cnt*sizeof(double));
     data2=(double *)malloc(cnt*sizeof(double));
     data3=(double *)malloc(cnt*sizeof(double));
     dataT=(double *)malloc(cntT*sizeof(double));
     dataF=(double ***)malloc(cntT*sizeof(double **));
     dataD=(double **)malloc(cntT*sizeof(double *));
     for(i=0; i<cntT; i++) {
       dataF[i]=(double **)malloc(rangeI*sizeof(double *));
       for(j=0; j<rangeI; j++) 
         dataF[i][j]=(double *)malloc(2*sizeof(double ));
       dataD[i]=(double *)malloc(rangeI*sizeof(double));
     }

     for(j=0; j<cntT; j++) 
       for(i=0; i<rangeI; i++) {
         for(k=0; k<2; k++) dataF[j][i][k]=0.0;
         dataD[j][i]=0.0;
       }

     while(step<=final) {
       j=(step-initial)/saveStep;
       dataT[j]=step;
       if(dimension==1) {
         sprintf(fileName,"cenfieldE%d_Ex",step);
         sprintf(fileName1,"cenfieldE%d_Ey",step);
         sprintf(fileName2,"cen0density%d_0",step);
         if(fopen(fileName,"r")==NULL) printf("%s is not exited.\n",fileName);
         if(fopen(fileName1,"r")==NULL) printf("%s is not exited.\n",fileName1);
         if(fopen(fileName2,"r")==NULL) printf("%s is not exited.\n",fileName2);
       }
       else { 
//         if(mode==0)   sprintf(fileName,"cenPukhov%d_%d",step,core);
//         if(mode==0)   sprintf(fileName,"cenYee%d_%d",step,core);
//         else          sprintf(fileName,"%s%d_%d",dataType,step,core);
//         if(fopen(fileName,"r")==NULL) printf("%s is not exited.\n",fileName);
          sprintf(fileName,"%s%d_%d",dataType,step,core);
       }
  
       
       if(dimension==1) {
         in = fopen(fileName,"r");
         for(i=0; i<cnt; i++) {
           fscanf(in,"%lf %lf",&x,&data1[i]); dataX[i]=x-dt*step; }
         fclose(in);

         in = fopen(fileName1,"r");
         for(i=0; i<cnt; i++) fscanf(in,"%lf %lf",&x,&data2[i]); 
         fclose(in);

         in = fopen(fileName2,"r");
         for(i=0; i<cnt; i++) fscanf(in,"%lf %lf",&x,&data3[i]); 
         fclose(in);
       } else {
         in = fopen(fileName,"r");
         for(i=0; i<cnt; i++) { 
           fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&data1[i],&data2[i],&Ez,&Bx,&By,&Bz,&data3[i]); 
           dataX[i]=x-dt*step;
         }
         fclose(in);
       }
       printf("%s is done\n",fileName);
      
       for(i=0; i<cnt; i++) { 
         x=minX+i*dx;
         index=(int)((dataX[i]-minX)/dx);
         weight=(dataX[i]-minX)/dx-index;
         if(index>=0 && index<rangeI) {
           dataF[j][index][0]=weight*data1[i+1]+(1.0-weight)*data1[i];
           dataF[j][index][1]=weight*data2[i+1]+(1.0-weight)*data2[i];
           dataD[j][index]=weight*data3[i+1]+(1.0-weight)*data3[i];
         } else ;
       }

       step+=saveStep;
     }

     sprintf(fileName,"cenXT");
     out = fopen(fileName,"w");
printf("cntT=%d, rangeI=%d\n",cntT,rangeI);
     for(j=0; j<cntT; j++) {
       for(i=0; i<rangeI; i++) {
         x=minX+i*dx;
         fprintf(out,"%.10g %lf %lf %lf %lf\n",x,dataT[j],dataF[j][i][0],dataF[j][i][1],dataD[j][i]);
       }
       fprintf(out,"\n");
     }
     fclose(out);    
     printf("%s is made.\n",fileName);
     
     for(i=0; i<cntT; i++)  {
       for(j=0; j<rangeI; j++) free(dataF[i][j]);
       free(dataF[i]);
       free(dataD[i]);
     }
     free(dataF); free(dataD); 
     free(dataX); free(dataT);
     free(data1); free(data2); free(data3);
     break;
/*
   case 1 :
     rangeX = atof(argv[6]);
     dt = atof(argv[7]);
     dx = atof(argv[8]);
     core = atoi(argv[9]);

     step=initial;
     sprintf(fileName,"cen%d_0",step); 

     if(fopen(fileName,"r")==NULL) {
       printf("%s is not exited.\n",fileName);
       cnt=1;
       step=final+saveStep;
     } else {
       in = fopen(fileName,"r");
       cnt=0;
       while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF) 
           cnt++;
       fclose(in);
       step=initial;
     }

     minX=-rangeX;
     rangeI=(int)(rangeX/dx)+1;

     dataX=(double *)malloc(cnt*sizeof(double));
     data1=(double *)malloc(cnt*sizeof(double));
     data2=(double *)malloc(cnt*sizeof(double));
     data3=(double *)malloc(cnt*sizeof(double));
     dataT=(double *)malloc(cntT*sizeof(double));
     dataF=(double ***)malloc(cntT*sizeof(double **));
     dataD=(double **)malloc(cntT*sizeof(double *));
     for(i=0; i<cntT; i++) {
       dataF[i]=(double **)malloc(rangeI*sizeof(double *));
       for(j=0; j<rangeI; j++) 
         dataF[i][j]=(double *)malloc(2*sizeof(double ));
       dataD[i]=(double *)malloc(rangeI*sizeof(double));
     }

     for(j=0; j<cntT; j++) 
       for(i=0; i<rangeI; i++) {
         for(k=0; k<2; k++) dataF[j][i][k]=0.0;
         dataD[j][i]=0.0;
       }

     while(step<=final) {
       j=(step-initial)/saveStep;
       dataT[j]=step;
       sprintf(fileName,"cen%d_0",step); 
       if(fopen(fileName,"r")==NULL) printf("%s is not exited.\n",fileName);
  
       in = fopen(fileName,"r");
       for(i=0; i<cnt; i++) {
         fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&data1[i],&data2[i],&Ez,&Bx,&By,&Bz,&data3[i]); 
         dataX[i]=x-dt*step;
       }
       fclose(in);
       printf("%s is done\n",fileName);
      
       for(i=0; i<cnt; i++) { 
         x=minX+i*dx;
         index=(int)((dataX[i]-minX)/dx);
         weight=(dataX[i]-minX)/dx-index;
         if(index>=0 && index<rangeI) {
           dataF[j][index][0]=weight*data1[i+1]+(1.0-weight)*data1[i];
           dataF[j][index][1]=weight*data2[i+1]+(1.0-weight)*data2[i];
           dataD[j][index]=weight*data3[i+1]+(1.0-weight)*data3[i];
         } else ;
       }

       step+=saveStep;
     }

     sprintf(fileName,"cenXT");
     out = fopen(fileName,"w");
     for(j=0; j<cntT; j++) {
       for(i=0; i<rangeI; i++) {
         x=minX+i*dx;
         fprintf(out,"%.10g %lf %lf %lf %lf\n",x,dataT[j],dataF[j][i][0],dataF[j][i][1],dataD[j][i]);
       }
       fprintf(out,"\n");
     }
     fclose(out);    
     printf("%s is made.\n",fileName);
     
     for(i=0; i<cntT; i++)  {
       for(j=0; j<rangeI; j++) free(dataF[i][j]);
       free(dataF[i]);
       free(dataD[i]);
     }
     free(dataF); free(dataD); 
     free(dataX); free(dataT);
     free(data1); free(data2); free(data3);
     break;
*/
   case 2 :
     sprintf(comp,"%s",argv[7]);
     pickX=atof(argv[8]);
     dx=atof(argv[9]);
 
     step=initial;
//     sprintf(fileName,"%s%d%s",dataType,step,comp);
     sprintf(fileName,"%s%d_%s",dataType,step,comp);
     in = fopen(fileName,"r");
     fscanf(in,"%lf %lf %lf",&x1,&y,&val);
     fscanf(in,"%lf %lf %lf",&x2,&y,&val);
     fclose(in);

     while(step<=final) {
       sprintf(fileName,"%s%d_%s",dataType,step,comp);
       sprintf(outFile,"pick%s",fileName);
       out = fopen(outFile,"w");

       in = fopen(fileName,"r");
       while(fscanf(in,"%lf %lf %lf",&x,&y,&val)!=EOF) {
         if(x>=pickX && x<pickX+dx) fprintf(out,"%g %g %g\n",x,y,val); else;
       }
       fclose(in);
       fclose(out);
       printf("%s is made.\n",outFile);
       step+=saveStep;
     }
     break;
/*
   case 1 :
     core = atoi(argv[5]);
     dx = atof(argv[6]);
     minX = atof(argv[7]);
     
     data1=(double *)malloc(cntT*sizeof(double));
     for(i=0; i<cntT; i++) data1[i]=0.0;

     for(step=initial; step<=final; step+=saveStep) {
       i=(step-initial)/saveStep;
       sprintf(fileName,"field%d_%d",step,core);

       if(fopen(fileName,"r")==NULL) {
         printf("%s is not exited.\n",fileName);
         step=final+saveStep;
       } else {
         in = fopen(fileName,"r");
         while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz)!=EOF)        {
           x=x-dx*step;
           if(x>minX) data1[i]+=Ey*Ey*dx; else ;
         }
         fclose(in);
       }
       printf("%s is done\n",fileName);
     }

     sprintf(outFile,"sumLaser");
     out = fopen(outFile,"w");
     for(i=0; i<cntT-1; i++)
       fprintf(out,"%d %lf %lf\n",initial+i*saveStep,data1[i],(data1[i]-data1[i+1])/((double)(saveStep)));
     fclose(out);
     printf("%s is done\n",outFile);

     free(data1);
     break;
*/
/*
   case 0 :
     core = atoi(argv[5]);
     dt = atof(argv[6]);

     step=initial;
     sprintf(fileName,"field%d_%d",step,core);
//     sprintf(fileName,"cen%d_%d",step,core);
     if(fopen(fileName,"r")==NULL) {
       printf("%s is not exited.\n",fileName);
       cnt=1;
       step=final+saveStep;
     } else {
       in = fopen(fileName,"r");
       cnt=0;
//       while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF) 
       while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz)!=EOF) 
         cnt++;
       fclose(in);
       step=initial;
     }

     dataX=(double *)malloc(cnt*sizeof(double));
     dataT=(double *)malloc(cntT*sizeof(double));
     data1=(double *)malloc(cnt*sizeof(double));
     data2=(double *)malloc(cnt*sizeof(double));
     dataF=(double **)malloc(cntT*sizeof(double *));
     dataD=(double **)malloc(cntT*sizeof(double *));
     for(i=0; i<cntT; i++) {
       dataF[i]=(double *)malloc(cnt*sizeof(double));
       dataD[i]=(double *)malloc(cnt*sizeof(double));
     }
     for(j=0; j<cntT; j++) 
       for(i=0; i<cnt; i++) { dataF[j][i]=0.0; dataD[j][i]=0.0; }

     while(step<=final) {
       j=(step-initial)/saveStep;
       dataT[j]=step;
//       sprintf(fileName,"cen%d_%d",step,core);
       sprintf(fileName,"field%d_%d",step,core);
       in = fopen(fileName,"r");
       for(i=0; i<cnt; i++) { 
//         fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&data1[i],&Ey,&Ez,&Bx,&By,&Bz,&data2[i]); 
         fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",&x,&data1[i],&Ey,&Ez,&Bx,&By,&Bz); 
         dataX[i]=x-dt*step;
       }
       fclose(in);
       printf("%s is done\n",fileName);
      
       dx=dataX[1]-dataX[0];
       minX=-(cnt-1)*dx;

       for(i=0; i<cnt; i++) { 
         x=minX+i*dx;
         index=(int)((x-minX)/dx);
         weight=(x-minX)/dx-index;
         if(index>=0 && index<cnt-1) {
           dataF[j][i]=weight*data1[index+1]+(1.0-weight)*data1[index];
           dataD[j][i]=weight*data2[index+1]+(1.0-weight)*data2[index];
         } else ;
       }
       
       step+=saveStep;
     }

     sprintf(fileName,"cenXT");
     out = fopen(fileName,"w");
     for(j=0; j<cntT; j++) {
       step=dataT[j];
       for(i=0; i<cnt; i++) {
         x=minX+i*dx;
         fprintf(out,"%lf %d %lf %lf\n",x,step,dataF[j][i],dataD[j][i]);
//         printf("%lf %d %lf %lf\n",x,step,dataF[j][i],dataD[j][i]);
       }
       fprintf(out,"\n");
     }
     fclose(out);    
     printf("%s is made.\n",fileName);
     
     free(dataX); free(dataT);
     for(i=0; i<cntT; i++) { free(dataF[i]); free(dataD[i]); }
     free(dataF); free(dataD);
     free(data1); free(data2);

     break;
*/
   }

   return 0;

}

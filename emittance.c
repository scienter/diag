#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void main(int argc, char *argv[])
{
   int initial,final,timeStep,step,species,nStep,i,j,indexI,mode;
   int ii,jj,nx,ny;
   double x,y,z,px,py,pz,id,core,ptcls,gamma,sliceI;
   double cnt,aveX2,aveXPrime2,aveCrsX,aveGam,emittanceX,xPrime;
   double aveY2,aveYPrime2,aveCrsY,emittanceY,yPrime;
   double *emitX,*betaX,*gammaX,*alphaX;
   double NemitX,NbetaX,NgammaX,NalphaX;
   double NemitY,NbetaY,NgammaY,NalphaY;
   double dx,dy,rangeX,rangeY,minX,minY,wx,wy,**den;
   FILE *in,*out;
   char fileName[100],strTemp[255];

   if(argc < 2)
   {  printf("emittance mode \n");
      printf("mode 1 : start end step species\n");
      printf("mode 2 (simplex) twiss : fileName\n");
      printf("mode 3 (0Particle) density : fileName rangeX nx rangeY ny\n");
      printf("mode 4 (simplex) density : fileName rangeX nx rangeY ny\n");
      exit(0);
   } else ;


   mode=atoi(argv[1]);
   switch (mode) {
   case 1 :
     initial=atoi(argv[2]);
     final=atoi(argv[3]);
     timeStep=atoi(argv[4]);
     species=atoi(argv[5]);
     nStep=(final-initial)/timeStep+1;

     emitX=(double *)malloc(nStep*sizeof(double ));
     betaX=(double *)malloc(nStep*sizeof(double ));
     gammaX=(double *)malloc(nStep*sizeof(double ));
     alphaX=(double *)malloc(nStep*sizeof(double ));
     for(i=0; i<nStep; i++) {
        emitX[i]=0.0;
        betaX[i]=0.0;
        alphaX[i]=0.0;
        gammaX[i]=0.0;
     }

     out=fopen("twiss","w");
   
     for(step=initial; step<=final; step+=timeStep) {
       sprintf(fileName,"%dParticle%d",species,step);
       in=fopen(fileName,"r");

        cnt=0.0;
        aveX2=0.0; aveXPrime2=0.0; aveCrsX=0.0; aveGam=0.0;
        while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&z,&x,&y,&pz,&px,&py,&id,&core,&ptcls)!=EOF) {
//        while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&z,&x,&y,&gamma,&px,&py,&ptcls,&id,&sliceI,&core)!=EOF) {
          gamma=sqrt(1+px*px+py*py+pz*pz);
          pz=sqrt(gamma*gamma-1+px*px+py*py);
  	  xPrime=px/pz;
	  if(isnan(xPrime)) printf("pz=%g\n",pz);

	  aveX2+=x*x;
 	  aveXPrime2+=xPrime*xPrime;
	  aveCrsX+=x*xPrime;
	  aveGam+=gamma;

	  cnt+=1.0;
        }
        emittanceX=sqrt((aveX2*aveXPrime2-aveCrsX*aveCrsX)/cnt/cnt);
	printf("emittanceX=%g\n",emittanceX);
        aveGam/=cnt;

        indexI=(step-initial)/timeStep;
        emitX[indexI]=emittanceX;
        betaX[indexI]=aveX2/cnt/emittanceX;
        gammaX[indexI]=aveXPrime2/cnt/emittanceX;
        alphaX[indexI]=-aveCrsX/cnt/emittanceX;

        fprintf(out,"%d %g %g %g %g %g\n",step,emitX[indexI],betaX[indexI],gammaX[indexI],alphaX[indexI],aveGam);
        fclose(in);
        printf("%s is done\n",fileName);
     }
     fclose(out);
     printf("twiss is made\n");

     free(emitX);
     free(betaX);
     free(gammaX);
     free(alphaX);

     break ;

   case 2 :
     NemitX=NemitY=0.0;
     NbetaX=NbetaY=0.0;
     NgammaX=NgammaY=0.0;
     NalphaX=NalphaY=0.0;

     out=fopen("twiss","w");
     sprintf(fileName,"%s",argv[2]);
     in=fopen(fileName,"r");

     cnt=0.0;
     aveX2=0.0; aveXPrime2=0.0; aveCrsX=0.0;
     aveY2=0.0; aveYPrime2=0.0; aveCrsY=0.0;
     fgets(strTemp,sizeof(strTemp),in);
     while(fscanf(in,"%lf %lf %lf %lf %lf %lf",&z,&x,&y,&px,&py,&gamma)!=EOF) {
       aveX2+=x*x;
       aveXPrime2+=px*px;
       aveCrsX+=x*px;
       aveY2+=y*y;
       aveYPrime2+=py*py;
       aveCrsY+=y*py;
       aveGam+=gamma;

       cnt+=1.0;
     }
     emittanceX=sqrt((aveX2*aveXPrime2-aveCrsX*aveCrsX)/cnt/cnt);
     emittanceY=sqrt((aveX2*aveXPrime2-aveCrsX*aveCrsX)/cnt/cnt);
     aveGam/=cnt;

     NemitX=emittanceX;
     NbetaX=aveX2/cnt/emittanceX;
     NgammaX=aveXPrime2/cnt/emittanceX;
     NalphaX=-aveCrsX/cnt/emittanceX;

     NemitY=emittanceY;
     NbetaY=aveY2/cnt/emittanceY;
     NgammaY=aveYPrime2/cnt/emittanceY;
     NalphaY=-aveCrsY/cnt/emittanceY;

     printf("emitX=%g, betaX=%g, gammaX=%g, alphaX=%g\n",NemitX,NbetaX,NgammaX,NalphaX);
     printf("emitY=%g, betaY=%g, gammaY=%g, alphaY=%g, aveGam=%g\n",NemitY,NbetaY,NgammaY,NalphaY,aveGam);

     break ;

   case 3 :
   case 4 :
     rangeX=atof(argv[3]);
     nx=atoi(argv[4]);
     rangeY=atof(argv[5]);
     ny=atoi(argv[6]);

     minX=-1.0*rangeX*0.5;
     minY=-1.0*rangeY*0.5;
     dx=rangeX/(1.0*nx);
     dy=rangeY/(1.0*ny);

     den=(double **)malloc((nx+1)*sizeof(double *));
     for(i=0; i<=nx; i++)
       den[i]=(double *)malloc((ny+1)*sizeof(double ));
     for(i=0; i<=nx; i++) 
       for(j=0; j<=ny; j++) 
         den[i][j]=0.0;

     sprintf(fileName,"%s",argv[2]);
     in=fopen(fileName,"r");
     if(mode==3) {
       while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&z,&x,&y,&gamma,&px,&py,&ptcls,&id,&sliceI,&core)!=EOF) {
         ii=(x-minX)/dx;
	 jj=(y-minY)/dy;
	 wx=(x-minX)/dx-ii;
	 wy=(y-minY)/dy-jj;
         if(ii>=0 && ii<nx && jj>=0 && jj<ny) {
           den[ii][jj]+=(1.0-wx)*(1.0-wy);	
           den[ii][jj+1]+=(1.0-wx)*wy;	
           den[ii+1][jj]+=wx*(1.0-wy);	
           den[ii+1][jj+1]+=wx*wy;	
	 }
       }
     } else {
       fgets(strTemp, sizeof(strTemp), in );
       while(fscanf(in,"%lf %lf %lf %lf %lf %lf",&z,&x,&y,&px,&py,&gamma)!=EOF) {
 	 ii=(x-minX)/dx;
	 jj=(y-minY)/dy;
	 wx=(x-minX)/dx-ii;
	 wy=(y-minY)/dy-jj;
         if(ii>=0 && ii<nx && jj>=0 && jj<ny) {
           den[ii][jj]+=(1.0-wx)*(1.0-wy);	
           den[ii][jj+1]+=(1.0-wx)*wy;	
           den[ii+1][jj]+=wx*(1.0-wy);	
           den[ii+1][jj+1]+=wx*wy;	
	 }
       }
     }
     fclose(in);

     out=fopen("density","w");
     for(i=0; i<=nx; i++) {
       x=minX+i*dx;
       for(j=0; j<=ny; j++) {
         y=minY+j*dy;
	 fprintf(out,"%g %g %g\n",x,y,den[i][j]);
       }
       fprintf(out,"\n");
     }     
     fclose(out);
     printf("density is made\n");

     for(i=0; i<=nx; i++) free(den[i]); free(den);

     break ;

   }

}

// From x-y ascii data from xgrafix, reconstruct the density profile
#include "stdio.h"
#include "stdlib.h"
#include "math.h"


void main(int argc, char *argv[])
{
   FILE *in,*out;
   char fileName[100],dataName[100],fileType[100];
   int mode,initial,final,saveStep,step,polarity,index,core,species,cnt;
   float x,y,Ex,Ey,Ez,Bx,By,Bz,minEy,maxEy,testEy,sum,position,kp,prevX,prevEx,low,up;
   float px,py,pz,maxPx,rangeX,minX,vg,gamma0,dx,backX,den;
   kp=2.66e5;
   dx=4e-8;
   gamma0=29.509;
   vg=1.0-0.5/gamma0/gamma0;

   if(argc < 4)
   {
      printf("findMax mode initial final saveStep\n");
      printf("mode (0) : fieldType core rangeX polarity\n");
      printf("mode (1) : initial final saveStep species\n");
      exit(0);
   }

   mode=atoi(argv[1]);
   initial=atoi(argv[2]);
   final=atoi(argv[3]);
   saveStep=atoi(argv[4]);


   if(mode==0) 
   {
      core=atoi(argv[6]);
      rangeX=atof(argv[7]);
      polarity=atof(argv[8]);
		sprintf(fileType,"%s",argv[5]);
      sprintf(fileName,"%s%d_%d",fileType,initial,core);
      if(fopen(fileName,"r")==NULL)  {
         printf("%s is not exited.\n",fileName);
         exit(0);
      } else ;		
     
      out = fopen("laserMax","w");
      for(step=initial; step<=final; step+=saveStep)
      {
//       sprintf(fileName,"cenField%d_%d",step,core);
        sprintf(fileName,"%s%d_%d",fileType,step,core);
        in = fopen(fileName,"r");
//     fgets(str,100,in);
        maxEy=-1000;
        minEy=10000;
        sum=0.0;
        minX=step*dx-rangeX;
        prevEx=-1;
        cnt=1;
//       while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&y,&Ex,&Ey,&Ez,&Bx,&By,&Bz)!=EOF)
        while(fscanf(in,"%g %g %g %g %g %g %g %g",&x,&Ex,&Ey,&Ez,&Bx,&By,&Bz,&den)!=EOF)
        {
          if(polarity==2)
          {
            if(Ey>maxEy) maxEy=Ey; else;
            if(Ey<minEy) minEy=Ey; else;
            sum+=Ey*Ey;

          }
          else if(polarity==3)
          {
            if(Ey>maxEy) maxEy=Ey; else;
            if(Ey<minEy) minEy=Ey; else;
            sum+=Ez*Ez;
          }
          else if(polarity==1)
          {
            testEy=Ex;
            if(testEy<maxEy && x>minX)  {
              maxEy=testEy;
              position=x*kp;
              if(prevEx>0 && Ex<0 && cnt==1)  {
                backX=x*kp;
                cnt++;
              }
            }
            prevEx=Ex;
          }
        }
        fclose(in);
        printf("%s is made\n",fileName);

        low=fabs(minEy);
        up=fabs(maxEy);
        if(low>up) maxEy=up;
        else 	   maxEy=low;

        if(polarity==1) 
          printf("%d, maxE=%g, position=%g, backX=%g\n",step,maxEy,position,backX);
        else 
          fprintf(out,"%d %g %g\n",step,maxEy,sum);
      }
      fclose(out);
      printf("laserMax is made\n");
   }

}



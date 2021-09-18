#include "stdio.h"
#include "stdlib.h"
#include "math.h"
main(int argc, char *argv[])
{
   float x,y,z,ux,uy,uz,id,core,gamma,min,max,coefY,coefZ;
   float emitY,emitZ,dy,dz,minY,minZ,maxY,maxZ,den;
   int num,i,j,indexY,indexZ,cnt;
   FILE *in;
   char str[100];

   if(argc < 2)
   {  printf("emittance file numData \n");  exit(0);  }

   num=atoi(argv[2]);

   in=fopen(argv[1],"r");
   cnt=0;
   minY=minZ=100;
   maxY=maxZ=-100;
   while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&id,&core,&den)!=EOF) 
   {
     emitY=uy/ux;
     emitZ=uz/ux;
     if(maxY<emitY)  maxY=emitY;
     if(maxZ<emitZ)  maxZ=emitZ;
     if(minY>emitY)  minY=emitY;
     if(minZ>emitZ)  minZ=emitZ;
     cnt=cnt+1;
   }
   fclose(in);
   if(minY>minZ) min=minZ;
   else 	 min=minY;
   if(maxY>maxZ) max=maxY;
   else 	 max=maxZ;
   dz=dy=(max-min)/(float)num;
printf ("max=%g, min=%g,dy=%g,dz=%g\n",max,min,dy,dz);

   float nEmit[num+1][num+1],yy[num+1],zz[num+1];
   for(i=0; i<num+1; i++)
     for(j=0; j<num+1; j++)
       nEmit[i][j]=0.0;
   for(i=0; i<num+1; i++)
     yy[i]=i*dy+min;
   for(j=0; j<num+1; j++)
     zz[j]=j*dz+min;

   in = fopen(argv[1],"r");
//   fgets(str,100,in);
   while(fscanf(in,"%g %g %g %g %g %g %g %g %g",&x,&y,&z,&ux,&uy,&uz,&id,&core,&den)!=EOF)
   {      
     emitY=uy/ux;
     emitZ=uz/ux;
     indexY=(int)((emitY-min)/dy);
     indexZ=(int)((emitZ-min)/dz);
     coefY=(emitY-min)/dy-indexY;
     coefZ=(emitZ-min)/dz-indexZ;
     nEmit[indexY][indexZ]+=(1.0-coefY)*(1.0-coefZ)*den;
     nEmit[indexY][indexZ+1]+=(1.0-coefY)*coefZ*den;
     nEmit[indexY+1][indexZ]+=coefY*(1.0-coefZ)*den;
     nEmit[indexY+1][indexZ+1]+=coefY*coefZ*den;
   }
   fclose(in);

   for(i=0; i<num+1; i++)
   {
     for(j=0; j<num+1; j++)
     {
       printf("%g %g %g\n",yy[i],zz[j],fabs(nEmit[i][j]));
     }
     printf("\n");
   }
}

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

void main(int argc, char *argv[])
{
   int mode;
   double x,y,z,px,py,pz,macro,tmp,id;
   double energy0,sum,unitP,eCharge;
   FILE *in,*out;
   char fileName[100],outFile[100];

   unitP=0.511e6;
   eCharge=1.602e-19;

   if(argc == 0 || argc < 2 ) { 
     printf("conAtoP mode\n");
     printf("mode 1 : fileName mean_energy[MeV]\n");
     exit(0);
   } else ;

   mode=atoi(argv[1]);

   switch (mode) {
   case 1 :
     sprintf(fileName,"%s",argv[2]);
     energy0=atof(argv[3])*1e6;	 //eV
     in=fopen(fileName,"r");
     sum=0.0;
     sprintf(outFile,"beam");
     out=fopen(outFile,"w");
     while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&px,&py,&pz,&tmp,&macro,&id,&tmp)!=EOF) {
        pz=(energy0+pz)/unitP;
        px=px/unitP;
        py=py/unitP;
	macro*=1e-9/eCharge;
	sum+=macro;
        if(z<1) fprintf(out,"%g %g %g %g %g %g %g 0 %g\n",z,x,y,pz,px,py,id,-macro); else ;
     }
     fclose(in);

     fclose(out);
     printf("%s is made\n",outFile);

     printf("total Charge=%g\n",sum*eCharge);


     
     break;
   }
}

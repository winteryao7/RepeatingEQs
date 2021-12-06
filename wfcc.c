/*********************************************************************
*	wfcc.c:
*
*	Usage:
*		wfcc -Dlength_of_correlation_window/max_time_shift [-W] [-N]
*		then input name of file from stdin, in following format:
*			name arr
*		where arr is the approx. arrival time for the trace.
*		The first line is for the master trace.
*		
*		The outputs are in the same format so it can be used
*		as input for iteration
*
*	Author:  Zhigang Peng
*
*	Revision History
*		June 1997	Initial coding by Zhu
*		06/23/97	input component names from stdio with
*				a shift0 by Zhu.
*		09/04/97	change shift0 to arrival time by Zhu.
*		02/11/02	output cross-correlation value by Peng.
*		11/08/03	loop through all the possible pairs by Peng.
*		02/08/07	update the distance part by Peng.
*		02/08/07	change the stdin input to be arrival, lon, lat, and depth
*		06/11/07	fixed the bug for calculating distance
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"
#include "Complex.h"

#define MAXSEISMOGRAMS 40000
#define DEGREETOKM 111.195

#define TAPER	0.05		/* portion of window being tapered */

void usage(char *program);

int main(int argc, char **argv) {
  int 		i, j, k, nn, mm, mt8, t8, max_shift, error;
  int		shift, start, end, ntrace, overWrite, normalize;
  char		line[128], input[128];
  size_t	wndw_size;
  float		tshift;	/*max. time shift in sec.*/
  float		tBefore, tAfter, arr, max, min, moddata, delta, tmp1, tmp2, tmp3;
  float 	norm, dt, *data, *dpt, *crl;
  SACHEAD	hd;
  void		taper(float *, int);
  double         distance(double, double, double, double);

/* added variables */
  float         delay, cccf;
  double        ihd, ihdcf;
  float *seismogram[MAXSEISMOGRAMS], modseis[MAXSEISMOGRAMS], tarr[MAXSEISMOGRAMS];
  double lat[MAXSEISMOGRAMS], lon[MAXSEISMOGRAMS], dep[MAXSEISMOGRAMS];
  char filenames[MAXSEISMOGRAMS][128];
  int verbose;
  
  error = 0;
  tBefore = -5;
  tAfter = 10;
  tshift = 1;
  overWrite = 0;
  normalize = 1;

  ihdcf = 100;
  cccf = 0.4;

  verbose = 0;
/* inter-hypocentral distance cutoff value 
 * default is 10 km 
 * cross-correlation cutoff is 0.4 
 */

  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
  
       case 'C':
         sscanf(&argv[i][2],"%f/%f",&ihdcf,&cccf);
         break;
  
       case 'D':
         sscanf(&argv[i][2],"%f/%f/%f",&tBefore,&tAfter,&tshift);
         break;
  
       case 'V':
	 verbose = 1;
	 break;

       default:
         error = 1;
         break;
  
       }
    } else {
       error  = 1;
    }
  }

  if (argc == 1 || error) {
     fprintf(stderr,"saclst a evlo evla evdp f *.SAC | usage: %s -Dt1/t2/max_shift [-Cihdcf/cccf] [-V]\n",argv[0]);
     return -1;
  }

  /* input all traces and store it into memories */
  ntrace = 0;
  while (fgets(line,128,stdin)) {
     if (ntrace == MAXSEISMOGRAMS) {
	fprintf(stderr,"Error. this program only can load maximum %d seismograms, skipping the rest\n", ntrace);
	break;
     }
     sscanf(line,"%s %f %f %f %f",input,&arr,&tmp1,&tmp2,&tmp3);
     strcpy(filenames[ntrace],input);
     tarr[ntrace] = arr;
     lon[ntrace] = tmp1;
     lat[ntrace] = tmp2;
     dep[ntrace] = tmp3;
     if (verbose) {
        printf("Input file name is %s, arrival is %f, lat %f lon %f dep %f km\n",filenames[ntrace],tarr[ntrace],lon[ntrace],lat[ntrace],dep[ntrace]);
     }
     if ( (data=read_sac(filenames[ntrace],&hd)) == NULL ) return -1;
     if (ntrace == 0) {
        dt = hd.delta;
        max_shift = 2*rint(tshift/dt);
        mm = rint((tAfter-tBefore)/dt);
        wndw_size = mm*sizeof(float);
     }
     mt8 = rint((arr+tBefore-hd.b)/dt);
     if (mt8 < 0) {
       fprintf(stderr,"%s time before arr. is not long enough\n",filenames[ntrace]);
       return -1;
     }
     if ((dpt = (float *) malloc(wndw_size)) == NULL ) {
       fprintf(stderr, "unable to allocate memory\n");
       exit(-1);
     }

/* read the headers */
//     lat[ntrace] = hd.evla;
//     lon[ntrace] = hd.evlo;
//     dep[ntrace] = hd.evdp;

     if (verbose) {
       printf("arr = %f, tBefore = %f, hd.b = %f, dt = %f, mt8 = %d, mm = %d\n",arr,tBefore,hd.b,dt,mt8,mm);
     }
     for (k = mt8, j = 0; j < mm; j++,k++) {
       dpt[j] = data[k];
     }

     taper(dpt, mm);
     for(moddata=0.,i=0;i<mm;i++) moddata += dpt[i]*dpt[i];
     moddata = sqrt(moddata);
     modseis[ntrace] = moddata;
     

     if (verbose) {
        printf ("%s mod is %f\n",filenames[ntrace],modseis[ntrace]);
     }

     seismogram[ntrace] = dpt;
     ntrace++;
     if (verbose) {
        printf("Finish reading data for %s, No. %d, data length is %d\n",filenames[ntrace],ntrace, mm);
     }
  }
     
  for (j = 0; j < ntrace; j++) {
     for (k = j+1; k < ntrace; k++) {
       /* calculate the inter-hypo distance */
        delta = distance(lat[j],lon[j],lat[k],lon[k]);
        ihd = sqrt((dep[j]-dep[k])*(dep[j]-dep[k])+DEGREETOKM*DEGREETOKM*delta*delta);
	if (verbose) {
   		printf("%s %.4f %.4f %.3f km, %s %.4f %.4f %.3f km, inter-hypocentral distance is %.4f deg, %f km\n",filenames[j],lon[j],lat[j],dep[j],filenames[k],lon[k],lat[k],dep[k],delta,ihd);
	}
//        ihd = sqrt((dep[j]-dep[k])*(dep[j]-dep[k])+(lat[j]-lat[k])*DEGREETOKM*(lat[j]-lat[k])*DEGREETOKM + (lon[j]-lon[k])*DEGREETOKM*(lon[j]-lon[k])*DEGREETOKM);
	if (ihd <= ihdcf) {
           crl = crscrl(mm, seismogram[j], seismogram[k], max_shift);
           shift = 0;
           norm = 0;
           for(i=0;i<=max_shift;i++) {
             if (crl[i]>norm) {
	       shift = i;
	       norm = crl[i];
             }
           }
          free(crl);
          norm = norm/modseis[j]/modseis[k];
          shift -= max_shift/2;
          delay = shift*dt;
	  if (norm >= cccf) {
	    printf("%s %9.4f %s %9.4f %7.4f %7.4f %6.4f\n",filenames[j],tarr[j],filenames[k],tarr[k],ihd,delay,norm);
          fflush(stdout);
	  }else {
	    if (verbose) {	      
	      printf("%s %s %7.4f %7.4f %6.4f < cutoff cc value %6.4f\n",filenames[j],filenames[k],ihd,delay,norm,cccf);
	    }
	  }
       } else {
          if (verbose) {
             printf("%s %s inter-hypocentral distance is %f km larger than the cut-off value %f km \n",filenames[j],filenames[k],ihd,ihdcf);
          }
       }
    }
  }
  return 0;
}

void	taper(float *aa, int n)
{
  int i, m;
  float	tt, pi1;
  m = TAPER*n;
  pi1 = 3.1415926/m;
  for (i=0; i<m; i++) {
    tt = 0.5*(1.-cos(i*pi1));
    aa[i] *= tt;
    aa[n-i-1] *= tt;
  }
}

double distance(double stalat, double stalon, double evtlat, double evtlon) {
   double delta, az, baz;
   double scolat, slon, ecolat, elon;
   double a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk;
   double rhs1,rhs2,sph,rad,del,daz,dbaz,pi,piby2;
   pi=3.14159;
   piby2=pi/2.0; 
   rad=2.*pi/360.0;
/* 
c  
c scolat and ecolat are the geocentric colatitudes
c as defined by Richter (pg. 318)
c
c Earth Flattening of 1/298.257 take from Bott (pg. 3)
c
*/
   sph=1.0/298.257;

   scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(stalat*rad));
   ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(evtlat*rad));
   slon=stalon*rad;
   elon=evtlon*rad;
/*
c
c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
c     These are defined for the pt. 1
c
*/ 
   a=sin(scolat)*cos(slon);
   b=sin(scolat)*sin(slon);
   c=cos(scolat);
   d=sin(slon);
   e=-cos(slon);
   g=-c*e;
   h=c*d;
   k=-sin(scolat);
/*
c
c  aa - ee are the same as a - e, except for pt. 2
c
*/ 
   aa=sin(ecolat)*cos(elon);
   bb=sin(ecolat)*sin(elon);
   cc=cos(ecolat);
   dd=sin(elon);
   ee=-cos(elon);
   gg=-cc*ee;
   hh=cc*dd;
   kk=-sin(ecolat);
/*
c
c  Bullen, Sec 10.2, eqn. 4
c
*/
   del=acos(a*aa + b*bb + c*cc);
   delta=del/rad;
/*
c
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 1 is unprimed, so this is technically the baz
c
c  Calculate baz this way to avoid quadrant problems
c
*/ 
   rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.;
   rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.;
   dbaz=atan2(rhs1,rhs2);
   if (dbaz<0.0) {
      dbaz=dbaz+2*pi;
   } 
   baz=dbaz/rad;
/*
c
c  Bullen, Sec 10.2, eqn 7 / eqn 8
c
c    pt. 2 is unprimed, so this is technically the az
c
*/
   rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.;
   rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.;
   daz=atan2(rhs1,rhs2);
   if(daz<0.0) {
      daz=daz+2*pi;
   }
   az=daz/rad;
/*
c
c   Make sure 0.0 is always 0.0, not 360.
c
*/
   if(abs(baz-360.) < .00001) baz=0.0;
   if(abs(az-360.) < .00001) az=0.0;
   return delta;
}

/* Copyright etc. */

/*  bet.c - Brain Extraction Tool

    Stephen Smith, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    WWW:      http://www.fmrib.ox.ac.uk/fsl
    Email:    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or (at
    your option) any later version.
    
    This program is distributed in the hope that it will be useful, but in
    order that the University as a charitable foundation protects its
    assets for the benefit of its educational and research purposes, the
    University makes clear that no condition is made or to be implied, nor
    is any warranty given or to be implied, as to the accuracy of FSL, or
    that it will be suitable for any particular purpose or for use under
    any specific conditions.  Furthermore, the University disclaims all
    responsibility for the use which is made of FSL.  See the GNU General
    Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
    USA */

/* defines, includes and typedefs */

#include "libss/libss.h"
#include "libss/libavw.h"
#include "libss/libminc.h"
#include "libss/libtessa.h"

void usage(void);

/* usage */

void usage(void)
{
  printf("\nUsage: bet <input fileroot> <output fileroot> [options]\n\n");

  printf("-o : generate brain surface outline overlaid onto original image\n");
  printf("-m : generate binary brain mask\n");
  printf("-s : generate approximate skull image\n");
  printf("-n : don't generate segmented brain image output\n");
  printf("-b : generate surface mask in BIC brain-view format\n");
  printf("-f <fractional threshold> : fractional intensity threshold (0->1); default=0.5; smaller values give larger brain outline estimates\n");
  printf("-g <fractional threshold> : vertical gradient in fractional intensity threshold (-1->1); default=0; positive values give larger brain outline at bottom, smaller at top\n");
  printf("-h <threshold> : ratio for hyperintense voxels (>1) (best for t1);\n");
  printf("                 smaller values remove more brain; default=none\n");
  printf("-r : reversed intensities, like t2 and pd\n");
  printf("-t : apply thresholding to segmented brain image and mask\n\n");

  /*printf("-x : generate xtopol format output; file extensions=coo,dat\n");
  printf("-c : generate approximate image cost function output; fileroot=<output filroot>_cost\n");
  printf("-S : as above but \"colour code\" brain-skull distance\n");*/

  exit(1);
}

/* main */

#define TESSELATE_ORDER             6       /* 5 */
#define ITERATIONS_CHECK            500     /* 500 */
#define ITERATIONS                  2000    /* 2000 */
#define PASSES                      10      /* 10 */
#define BRAIN_THRESHOLD_DEFAULT     0.5     /* 0.5 */
#define COST_SEARCH                 7.0     /* 20 mm */
#define COST_UPPER_SEARCH           1.0     /* 2 mm */

// Weight for surface smoothness - tangential component. Higher value
// will favor equal edge lengths.
#define LAMBDA_TANGENT              0.5     /* 0.5 */

// Weight for image (volume) fit. May also be viewed as an integration
// step: a large value takes larger steps, requires fewer iterations.
// However, it's better to have small steps, more iterations.
#define LAMBDA_FIT                  0.05    /* 0.05 */

// For surface smoothness, put more weight on smoothing (normal component) when
// radius of curvature is less than RADIUSMIN, less weight when radius of curvature
// is more than RADIUSMAX.
#define RADIUSMAX                   10.0    /* 10.0 mm */
#define RADIUSMIN                   3.33    /* 3.33 mm */

#define SELF_INTERSECTION_THRESHOLD       1000    /* 4000 */
#define MAX_SELF_INTERSECTION_THRESHOLD  50000

// Used only with -s and -S options
#define SKULL_SEARCH                30      /* 30 mm */
#define SKULL_START                 -3      /* -3 mm */

/*#define DEBUG_NORMALS*/
/*#define DEBUG_EVOLVE*/

int main(argc, argv)
  int   argc;
  char  *argv [];
{
  /* vars */

FDT *in, *mask=NULL, *raw=NULL, threshold, thresh2, thresh98,
  hist_min=0, hist_max=0, medianval;
int x_size, y_size, z_size, x, y, z, i, pc=0, iters, pass=1,
  output_brain=1, output_bic=0, output_xtopol=0, output_cost=0, output_mask=0,
  output_overlay=0, output_skull=0, apply_thresholding=0, code_skull=0,
  fix_upper=0, reversed_intensities=0;
double cgx, cgy, cgz, radius, scale, ml=0, ml0=0, tmpf,
  brain_threshold=BRAIN_THRESHOLD_DEFAULT, 
  threshold_upper=0.0, ratio_upper=0.0,
  gradthresh=0, incfactor=0,
  rE = 0.5 * (1/RADIUSMIN + 1/RADIUSMAX), rF = 6 / (1/RADIUSMIN - 1/RADIUSMAX);
char filename[1000];
image_struct im;
points_struc *v;


  /* process arguments */

if (argc<3)
     usage();

minc_read(argv[1], &im);

in=im.i;
x_size=im.x; y_size=im.y; z_size=im.z;
scale=MIN(im.xv,MIN(im.yv,im.zv)); /* minimum voxel side length */

scale /= 2.0;   /* do oversampling along the normal line - much better! */

int nevals = (int)( COST_SEARCH / scale );
if( COST_SEARCH > nevals * scale ) nevals++;

for (i = 3; i < argc; i++) {
  if (!strcmp(argv[i], "-x"))
    output_xtopol=1;
  else if (!strcmp(argv[i], "-c"))
    output_cost=1;
  else if (!strcmp(argv[i], "-s"))
    output_skull=1;
  else if (!strcmp(argv[i], "-S"))
    /* colour-coded skull image */

    {
      output_skull=1;
      code_skull=1;
    }
  else if (!strcmp(argv[i], "-m"))
    output_mask=1;
  else if (!strcmp(argv[i], "-o"))
    output_overlay=1;
  else if (!strcmp(argv[i], "-n"))
    output_brain=0;
  else if (!strcmp(argv[i], "-b"))
    output_bic=1;
  else if (!strcmp(argv[i], "-t"))
    apply_thresholding=1;
  else if (!strcmp(argv[i], "-f"))
    /* fractional brain threshold */

    {
      i++;
      if (argc<i+1) /* option following f hasn't been given */
      {
	printf("Error: no value given following -f\n");
	usage();
      }
      brain_threshold=atof(argv[i]);
      if ( (brain_threshold<=0) || (brain_threshold>=1) )
      {
	printf("Error: value following -f must lie between 0 and 1\n");
	usage();
      }
    }
  else if (!strcmp(argv[i], "-h")) {
      i++;
      if (argc<i+1) { /* option following h hasn't been given */
	printf("Error: no value given following -h\n");
	usage();
      }
      ratio_upper=atof(argv[i]);
      fix_upper = 1;
    }
  else if (!strcmp(argv[i], "-r")) {
      reversed_intensities=1;
    }
  else if (!strcmp(argv[i], "-g"))
    /* gradient fractional brain threshold */

    {
      i++;
      if (argc<i+1) /* option following g hasn't been given */
      {
	printf("Error: no value given following -g\n");
	usage();
      }
      gradthresh=atof(argv[i]);
      if ( (gradthresh<-1) || (gradthresh>1) )
      {
	printf("Error: value following -g must lie between -1 and 1\n");
	usage();
      }
    }

  else
    usage();
}

brain_threshold=pow(brain_threshold,0.275);

if ( !output_xtopol && !output_cost && !output_skull && !output_mask && !output_overlay &&
     !output_brain && !output_bic ) {
  printf("No outputs requested!\n");
  usage();
}

  /* image preprocessing */

// this is not the best thing to do: remove it later.
im.min=im.max=0;

// Note: We must have threshold > thresh2 in the min/max below
//       to avoid possible division by zero.
// - The value of threshold is used only for the min value of lmax,
//   so it's not very important.
// - A larger value of thresh2 (similar to threshold) has the effect
//   of removing more brain, in particular bright meninges, but can
//   also forget tips of temporal lobes and cerebellum. A smaller
//   value does not hurt.
// - A smaller value of thresh98 causes a smaller value for medianval,
//   thus will have the effect of including more brain tissue. However,
//   a smaller value of thresh98 is more logical to use as it will 
//   exclude more very bright voxels in the calculation of medianval.
// - A smaller value of medianval will include more brain than a
//   larger value of medianval. medianval should be a dividing 
//   threshold between non-tissue and tissue.

find_thresholds(&im);
hist_min=im.min; 
hist_max=im.max; 
thresh2=im.thresh2; 
thresh98=im.thresh98; 
threshold=im.thresh;
medianval = im.medianval;

printf("hist_min=%f thresh2=%f thresh=%f thresh98=%f hist_max=%f\n",
       im.min,(double)thresh2,(double)threshold,(double)thresh98,im.max);
printf("THRESHOLD %f\n",(double)threshold);
printf("MEDIANVAL %f\n",(double)medianval);

c_of_g (im,&cgx,&cgy,&cgz);   // in voxel coords
cgx*=im.xv; cgy*=im.yv; cgz*=im.zv;     // now converted to real coords
printf("CofG (%f,%f,%f) mm\n",cgx,cgy,cgz);

radius = find_radius (im,im.xv*im.yv*im.zv);
printf("RADIUS %f\n",radius);

if( fix_upper ) {
  /* Upper threshold: 
     For t1 image, used to clip hyperintense white voxels in bone marrow.
     Should be about equal to white_mean + 3 * white_stddev as obtained
     from the preliminary classified image. For t2/pd image, do nothing. */
  threshold_upper = thresh98;
  printf("UPPER THRESHOLD %f\n", threshold_upper );
}

if (output_cost) /* prepare cost function image for writing into */
{
  raw = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);
  memset(raw,(unsigned char)0,sizeof(FDT)*x_size*y_size*z_size);
}

object *old = &ico;                /* start with icosohedral tessellation */
tessa((int)TESSELATE_ORDER,&old);  /* create tessellated triangles; level 5 gives 2562 points */
pc=points_list(old,&v);            /* convert triangles to vertex list */
free(old);                         /* don't need triangular tessellation any more */

/*printf("VERTICES %d\n",pc);*/

/* measure initial spacing for use later in self-intersection calculations */
ml0 = sqrt( (v[0].xorig-v[v[0].n[0]].xorig)*(v[0].xorig-v[v[0].n[0]].xorig) +
            (v[0].yorig-v[v[0].n[0]].yorig)*(v[0].yorig-v[v[0].n[0]].yorig) +
            (v[0].zorig-v[v[0].n[0]].zorig)*(v[0].zorig-v[v[0].n[0]].zorig) );
/*printf("ml0=%f\n",ml0);*/

#ifdef DEBUG_NORMALS
v = (points_struc *)realloc((void *)v,sizeof(points_struc)*2*(pc+10)); /* make space for storing normals */
#endif

double intersection=0;

while (pass>0) {

  /* correct centre of mass based on current mask (replaces neck-cropping) */
  if( pass > 1 ) {

    mask = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);
    im.i=mask;
    fill_surface(&im,v,pc);

    cgx = 0.0;
    cgy = 0.0;
    cgz = 0.0;
    double total = 0.0;
    double value;

    for(z=0; z<z_size; z++) {
      for(y=0; y<y_size; y++) {
        for(x=0; x<x_size; x++) {
          if( IA(mask,x,y,z) > 0.50 ) {
            value = IA(in,x,y,z);
            if( value >= threshold && value <= thresh98 ) {
              total += 1.0;
              cgx += (double) x;
              cgy += (double) y;
              cgz += (double) z;
            }
          }
        }
      }
    }
    cgx /= total; cgy /= total; cgz /= total;
    cgx *= im.xv; cgy *= im.yv; cgz *= im.zv;     // now converted to real coords
    printf("CofG (%f,%f,%f) mm\n",cgx,cgy,cgz);
    /* assume that radius and medianval are suitable */

    free( mask );
    mask = NULL;
    im.i = in;
  }

  /* initialize tessellation */

  /* scale vertex positions for this image, set surface small and allow to grow */
  for(i=0; i<pc; i++) {
    v[i].x = v[i].xorig * radius * 0.5 + cgx;
    v[i].y = v[i].yorig * radius * 0.5 + cgy;
    v[i].z = v[i].zorig * radius * 0.5 + cgz;
  }

  /* find brain surface */

  for(iters=0; iters<ITERATIONS; iters++) {
    /* find local surface normals */

    for(i=0; i<pc; i++) {
      double nx, ny, nz, tmpf;
      int k, l;

      nx=ny=nz=0.0;

      for(k=0; v[i].n[k]>-1; k++); /* find number of connections */

      for(l=0; l<k; l++) { /* for each pair of consecutive neighbours form a vector product to get normal */
        double adx = v[v[i].n[l]].x - v[i].x,
               ady = v[v[i].n[l]].y - v[i].y,
               adz = v[v[i].n[l]].z - v[i].z,
               bdx = v[v[i].n[(l+1)%k]].x - v[i].x,
               bdy = v[v[i].n[(l+1)%k]].y - v[i].y,
               bdz = v[v[i].n[(l+1)%k]].z - v[i].z;
        nx += ady*bdz - adz*bdy;
        ny += adz*bdx - adx*bdz;
        nz += adx*bdy - ady*bdx;
      }

      /* make the normal vector of length 1 */
      tmpf = sqrt(nx*nx+ny*ny+nz*nz);
      nx/=(double)tmpf; ny/=(double)tmpf; nz/=(double)tmpf;

      /* debug normals */

#ifdef DEBUG_NORMALS
      if (iters==ITERATIONS-1) { /* final iteration */
        v[pc+i].x=v[i].x+nx*20.0;
        v[pc+i].y=v[i].y+ny*20.0;
        v[pc+i].z=v[i].z+nz*20.0;
        v[pc+i].n[0]=i;
        v[pc+i].n[1]=-1;
      }
#endif

      v[i].nx=nx;
      v[i].ny=ny;
      v[i].nz=nz;
    }

    /* find mean connection length every now and then */
    /* add the 50 as the rate of change is highest at start; thus do 0,50,100,200,..... */

    if ( (iters==50) || (iters%100==0) ) {
      int l;

      ml=0;

      for(i=0; i<pc; i++) {
        double mml=0;
        for(l=0; v[i].n[l]>-1; l++) {
          mml += sqrt( (v[i].x-v[v[i].n[l]].x)*(v[i].x-v[v[i].n[l]].x) +
                       (v[i].y-v[v[i].n[l]].y)*(v[i].y-v[v[i].n[l]].y) +
                       (v[i].z-v[v[i].n[l]].z)*(v[i].z-v[v[i].n[l]].z) );
        }
        ml += mml/l;
      }
      ml /= pc;   // average edge length over all surface
    }

    /* increased smoothing for pass>1 */

    if (pass>1) {
      incfactor=pow((double)10.0,(double)pass-1.0);

      if (iters>ITERATIONS*0.75) {
        incfactor=4.0*(1.0-((double)iters)/ITERATIONS)*(incfactor-1.0) + 1.0;
      }
    }

    for(i=0; i<pc; i++) {       /* calculate tessellation update */
      /* variables, and setup k and normal */

      FDT lmin, lmax, lavg;
      int k, l;
      double d;
      double nx=v[i].nx, ny=v[i].ny, nz=v[i].nz, sx, sy, sz, mml, fit, sn, stx, sty, stz;

      /* average position of neighbours: smoothing vector */
      for(k=0; v[i].n[k]>-1; k++); /* find number of connections */

      /* s is vector from current vertex to the mean position of its neighbours */

      sx=sy=sz=mml=0;

      for(l=0; l<k; l++) {
        sx += v[v[i].n[l]].x;
        sy += v[v[i].n[l]].y;
        sz += v[v[i].n[l]].z;
        mml += sqrt( (v[i].x-v[v[i].n[l]].x)*(v[i].x-v[v[i].n[l]].x) +
                     (v[i].y-v[v[i].n[l]].y)*(v[i].y-v[v[i].n[l]].y) +
                     (v[i].z-v[v[i].n[l]].z)*(v[i].z-v[v[i].n[l]].z) );
      }

      sx = sx/k - v[i].x;
      sy = sy/k - v[i].y;
      sz = sz/k - v[i].z;
      mml /= k;

      /* part of s normal to surface, sn = n * (s.n)
         part of s tangent to surface, st = s - sn */

      sn = sx*nx + sy*ny + sz*nz; /* this is just the s.n part - will multiply by n later */

      stx = sx - nx * sn;
      sty = sy - ny * sn;
      stz = sz - nz * sn;

      /* find intensity-based part of cost function - local max */

      FDT lnew;
      int all_inside_image=1;
      double lthresh, local_brain_threshold=brain_threshold;

      fit=0;
      lmin=medianval;
      lmax=threshold;
      lavg=0.0;

      /* change local threshold if gradient threshold used */

      if (gradthresh!=0) {
        tmpf = brain_threshold + gradthresh * ( (v[i].z-cgz) / radius );
        local_brain_threshold = MIN( 1.0 , MAX( tmpf , 0.0 ) );
      }
   
      d=scale;   /* start d at 1 not 0 so that boundary is just outside brain not on the actual edge */
      x=FTOI((v[i].x-d*nx)/im.xv);   // convert from real to voxel coords
      y=FTOI((v[i].y-d*ny)/im.yv); 
      z=FTOI((v[i].z-d*nz)/im.zv);
      if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) ) {
        lnew=IA(in,x,y,z);
        lmin = MIN(lmin,lnew);
        lmax = MAX(lmax,lnew);
        lavg += lnew;
      } else {
        all_inside_image=0;
      }

      d=COST_SEARCH;  /* furthest point inside mask, towards center of brain */
      x=FTOI((v[i].x-d*nx)/im.xv);   // convert from real to voxel coords
      y=FTOI((v[i].y-d*ny)/im.yv); 
      z=FTOI((v[i].z-d*nz)/im.zv);
      if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) ) {
        lnew=IA(in,x,y,z);
        if( reversed_intensities ) {
          lmax = MAX(lmax,lnew);   // like t2, pd
        } else {
          lmin = MIN(lmin,lnew);   // like t1
        }
        lavg += lnew;
      } else {
        all_inside_image=0;
      }

      if (all_inside_image) {
        for(d=2.0*scale;d<COST_SEARCH;d+=scale) {
          x=FTOI((v[i].x-d*nx)/im.xv);   // convert from real to voxel coords
          y=FTOI((v[i].y-d*ny)/im.yv); 
          z=FTOI((v[i].z-d*nz)/im.zv);
          lnew=IA(in,x,y,z);
          if( reversed_intensities ) {
            lmax = MAX(lmax,lnew);   // like t2, pd
          } else {
            lmin = MIN(lmin,lnew);   // like t1
          }
          lavg += lnew;
          if (d<0.5*COST_SEARCH) { /* only look relatively locally for maximum intensity */
            if( reversed_intensities ) {
              lmin = MIN(lmin,lnew);   // like t2, pd
            } else {
              lmax = MAX(lmax,lnew);   // like t1
            }
          }
	}

        /* This works fine for all image types t1, t2, pd to stop at a minimum. */ 
        lmin = MAX(lmin,thresh2);   /* surely these two lines are now redundant? */
        lmax = MIN(lmax,medianval); /* this effectively limits the growth in the */
                                    /* middle of the brain while lmax>medianval */
        lthresh = (lmax - thresh2)*local_brain_threshold + thresh2;
        fit = (lmin - lthresh) / ((lmax - thresh2)*0.5); /* scale range to around -1:1 */

        // This is a test for inside large ventricles: if the maximum value is 
        // small enough, like all csf and no real tissue, then assume we are in
        // a large ventricle and don't allow inward motion. This vertex will be
        // moved by its neighbours using the smoothness constraint. We apply this
        // correction in the first half of the growth of the surface mask, which
        // should be enough to grow beyond the ventricles.

        if( iters <= ITERATIONS/2 &&
            !reversed_intensities && lthresh <= threshold ) fit = 0.0;

        /* Stop just before a local max for hyperintense voxels. */

        if( fix_upper && fit > 0 ) {
          FDT lmax_upper=medianval;
          /* search outside for the local max */
          for(d=-COST_UPPER_SEARCH;d<=COST_UPPER_SEARCH;d+=scale) {
            x=FTOI((v[i].x+d*nx)/im.xv);   // convert from real to voxel coords
            y=FTOI((v[i].y+d*ny)/im.yv); 
            z=FTOI((v[i].z+d*nz)/im.zv);
            if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) && (z>=0) && (z<z_size) ) {
              lmax_upper = MAX(lmax_upper,IA(in,x,y,z));
            }
          }
          lavg /= (FDT)nevals;
          if( lmax_upper > threshold_upper && lmax_upper > ratio_upper * lavg ) {
            fit = -( lmax_upper - threshold_upper ) / ( threshold_upper - medianval );
          }
        }

        if (output_cost) {
	  x=(v[i].x/im.xv); y=(v[i].y/im.yv); z=(v[i].z/im.zv);
	  if ( (x>=0) && (x<x_size) && (y>=0) && (y<y_size) &&(z>=0) && (z<z_size) )
	    IA(raw,x,y,z)=thresh98/2 + thresh2/2 + (0.5*fit*((double)thresh98 - (double)thresh2));
        }  
      }

      /* estimate the update */

      /* normal component of smoothing */
      tmpf = 2 * ABS(sn) / (mml*mml);            /* 1/r ; r=local radius of curvature */
      tmpf = 0.5 * ( 1 + tanh((tmpf-rE)*rF) );   /* f(1/r) ; this makes big r have low correction and vice versa */

      if ( (pass>1) && (sn > 0) ) { /* if need increased smoothing, to avoid self-intersection; */
                                    /* only apply to concave curvature (as seen from outside surface) */
        tmpf *= incfactor;
        tmpf = MIN(tmpf,1);
      }

      tmpf = sn*tmpf;                            /* multiply normal component by correction factor */

      tmpf += ml * LAMBDA_FIT * fit ;            /* combine normal smoothing and image fit */

      v[i].xnew = v[i].x + stx*LAMBDA_TANGENT + tmpf*nx;
      v[i].ynew = v[i].y + sty*LAMBDA_TANGENT + tmpf*ny;
      v[i].znew = v[i].z + stz*LAMBDA_TANGENT + tmpf*nz;

    }

    /* update tessellation */

    for(i=0; i<pc; i++) {
      v[i].x = v[i].xnew;
      v[i].y = v[i].ynew;
      v[i].z = v[i].znew;
    }

    /* test surface for non-spherecity */

    if( ( (iters+1)%ITERATIONS_CHECK == 0 ) || ( iters == ITERATIONS-1 ) ) {
      int j, l;

      intersection=0;
      for(i=0; i<pc; i++) { /* loop round all points */
        for(j=0; j<pc; j++) { /* inner loop round all points */
          int do_it=1;
	  
          if (j==i) /* other point is same as current one - don't use */
            do_it=0;
      
          for(l=0; v[i].n[l]>-1; l++) /* other point is connected to current one - don't use */
            if (j==v[i].n[l])
              do_it=0;

          if ( (do_it) &&
             ( (v[i].x-v[j].x)*(v[i].x-v[j].x) + (v[i].y-v[j].y)*(v[i].y-v[j].y) + (v[i].z-v[j].z)*(v[i].z-v[j].z) < ml*ml ) ) {
            double dist = sqrt ( (v[i].x-v[j].x)*(v[i].x-v[j].x) +
                                 (v[i].y-v[j].y)*(v[i].y-v[j].y) +
                                 (v[i].z-v[j].z)*(v[i].z-v[j].z) ),
            distorig = sqrt ( (v[i].xorig-v[j].xorig)*(v[i].xorig-v[j].xorig) +
                              (v[i].yorig-v[j].yorig)*(v[i].yorig-v[j].yorig) +
                              (v[i].zorig-v[j].zorig)*(v[i].zorig-v[j].zorig) );

            tmpf = (distorig/ml0) - (dist/ml);    /* orig distance (in units of mean length) - current distance */
            tmpf *= tmpf;                         /* weight higher values more highly */
            intersection += tmpf;
          }
        }
      }
      printf("self-intersection total after %d iterations = %f (threshold=%f)\n",
             iters+1, intersection,(double)SELF_INTERSECTION_THRESHOLD);
      if( intersection > MAX_SELF_INTERSECTION_THRESHOLD ) break;
    }
  }

  if (pass<PASSES) {
    if (intersection>SELF_INTERSECTION_THRESHOLD) {
      printf("thus will rerun with higher smoothness constraint\n");
      pass++;
    } else {
      pass = 0;
    }
  } else {
    printf( "Couldn't converge mincbet. This is the best you can have.\n" );
    pass=0;
  }
}

/* write brain image and mask */

if ( (output_mask) || (output_brain) ) {
  mask = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);
  im.i=mask;
  fill_surface(&im,v,pc);

  if (apply_thresholding) {
      for(i=0; i<z_size*y_size*x_size; i++)
	if (in[i]<threshold)
	  mask[i]=0;
  }
}

if (output_mask) {
  im.min=0;
  im.max=1;
  sprintf(filename,"%s_mask",argv[2]);
  minc_write(argv[1],MINCBET_BYTE,filename,im);
}

if (output_brain) {
  int goesneg=0;

  im.min=thresh2;
  im.max=thresh98;

  if ( (im.dtmin<0) && ((double)hist_min<0) )
    goesneg=1;

  for(i=0; i<z_size*y_size*x_size; i++) {
      if ( goesneg ) {
	  if (mask[i]<0.5)
	    mask[i]=0;
	  else
	    mask[i]=(FDT)( ((double)(MIN(MAX(in[i],thresh2),thresh98))-(double)thresh2) * (im.dtmax-1) / ((double)thresh98-(double)thresh2) + 1 );

	  im.min=0;
	  im.max=(FDT)im.dtmax;
      } else
	mask[i]=mask[i]*in[i];
  }

  im.i=mask;
  minc_write(argv[1],MINCBET_FLOAT,argv[2],im);
}

/* output cost */

if (output_cost) {
  im.i=raw;
  im.min=raw[0];
  im.max=raw[0];

  for(i=0; i<z_size*y_size*x_size; i++) {
      if (raw[i]>im.max) im.max=raw[i];
      if (raw[i]<im.min) im.min=raw[i];
  }

  sprintf(filename,"%s_cost",argv[2]);
  minc_write(argv[1],MINCBET_FLOAT,filename,im);
  free(raw);
}

/* find skull and output */

if (output_skull) {
  FDT *skull = (FDT *) malloc(sizeof(FDT)*x_size*y_size*z_size);

  im.i=skull;
  memset(skull,(unsigned char)0,sizeof(FDT)*x_size*y_size*z_size);
  draw_surface(&im,1,1,v,pc); /* tell skull searching where to start */

  /* create skull image */

  im.i=in;

  for(z=0; z<z_size; z++)
    for(y=0; y<y_size; y++)
      for(x=0; x<x_size; x++)
        if (IA(skull,x,y,z)==1) {
	  FDT val, val2, minval, maxval;
	  double nx, ny, nz, d, d_max=0, d_min, grad, maxgrad, X, Y, Z, maxJ, lastJ;
	  int xx, yy, zz, j=0;

	  /* zero this point */

          IA(skull,x,y,z)=0;

	  /* find nearest node and setup normal */

          double mind=10000000;

          for(i=0; i<pc; i++) { 
            tmpf = (x*im.xv-v[i].x)*(x*im.xv-v[i].x) +
                   (y*im.yv-v[i].y)*(y*im.yv-v[i].y) +
                   (z*im.zv-v[i].z)*(z*im.zv-v[i].z);
            if (tmpf<mind) {
              j=i;
	      mind=tmpf;
            }
          }

          nx=v[j].nx;
          ny=v[j].ny;
          nz=v[j].nz;

          /* find minval, d_max and maxval up to SKULL_SEARCH distance */

          maxval=threshold;
          minval=IA(in,x,y,z);

	  for(d=0; d<SKULL_SEARCH; d+=scale*0.5) {
	      xx=x+FTOI(d*nx/im.xv); yy=y+FTOI(d*ny/im.yv); zz=z+FTOI(d*nz/im.zv);

	      if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) ) {
		  val=IA(in,xx,yy,zz);

		  if (val>maxval) {
		      maxval=val;
		      d_max=d;
		  }

		  if (val<minval)
		      minval=val;
              }
	  }

	  if (maxval>threshold) { /* can we see an intensity peak on the far side of the skull? */
	      /* find furthest out point that has small intensity */

            maxgrad=1;
            d_min=SKULL_START;
            maxJ =-1000000;
            lastJ=-2000000;

	    for(d=SKULL_START; d<d_max; d+=scale*0.5) {
	      xx=x+FTOI(d*nx/im.xv); yy=y+FTOI(d*ny/im.yv); zz=z+FTOI(d*nz/im.zv);

	      if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) ) {
		  tmpf=d/30 - IA(in,xx,yy,zz)/((double)(thresh98-thresh2));

		  /*		  if ( (tmpf>maxJ) && (lastJ+1.0/30>tmpf*0.5) )*/
		  if (tmpf>maxJ) {
		    maxJ=tmpf;
		    d_min=d;
		  }
		  lastJ=tmpf;
              }
            }

	    /* find _first_ max gradient out from d_min */

            maxgrad=0;

            X=x+d_min*nx/im.xv; Y=y+d_min*ny/im.yv; Z=z+d_min*nz/im.zv;

            if ( (X>0) && (X<x_size-1) && (Y>0) && (Y<y_size-1) && (Z>0) && (Z<z_size-1) ) {
              val2 = TLI(im,X,Y,Z);

              for(d=d_min+scale; d<d_max; d+=scale*0.5) {
	        X=x+d*nx/im.xv; Y=y+d*ny/im.yv; Z=z+d*nz/im.zv;

	        if ( (X>0) && (X<x_size-1) && (Y>0) && (Y<y_size-1) && (Z>0) && (Z<z_size-1) ) {
		  val = TLI(im,X,Y,Z);
		  /*printf("%d %d %d   %f %f %f   %d %d\n",x,y,z,X,Y,Z,(int)val,(int)val2);*/

		  grad=val-val2;
		  val2=val;		 

		  if (grad>0) { /* this so that we don't do anything if we're still in the same voxel */
		      if (grad > maxgrad) {
			  maxgrad=grad;
			  d_min=d;
                      } else {
			d=d_max;
                      }
		  }
		}
	      }
            }

	    /* mark this point as skull */

            if (maxgrad>0) {
              xx=x+FTOI(d_min*nx/im.xv);
              yy=y+FTOI(d_min*ny/im.yv);
              zz=z+FTOI(d_min*nz/im.zv);

              if ( (xx>=0) && (xx<x_size) && (yy>=0) && (yy<y_size) && (zz>=0) && (zz<z_size) ) {
                if (code_skull)
	          IA(skull,xx,yy,zz)=(FDT)(d_min*100.0);
                else
	          IA(skull,xx,yy,zz)=100;
              }
            }

	  }
	}

  im.i=skull;

  /* output */

  im.min=0;
  im.max=100;
  sprintf(filename,"%s_skull",argv[2]);
  minc_write(argv[1],MINCBET_BYTE,filename,im);

  free(skull);
}

/* output overlay (corrupts input image) */

if (output_overlay) {
  im.min=thresh2;
  im.max=hist_max;
  im.i=in;
  draw_surface(&im,im.max,im.max,v,pc);
  sprintf(filename,"%s_overlay",argv[2]);
  minc_write(argv[1],MINCBET_SAME,filename,im);
}

/* output xtopol surface (corrupts tessellation) */

if (output_xtopol) {
  int xtopol_pc=pc;
  double ax=0, ay=0, az=0, tmpf=0;

#ifdef DEBUG_NORMALS
  xtopol_pc *= 2;
#endif

  for(i=0; i<xtopol_pc; i++) {
      ax += v[i].x;
      ay += v[i].y;
      az += v[i].z;
  }
  ax/=(double)xtopol_pc; ay/=(double)xtopol_pc; az/=(double)xtopol_pc; 

  for(i=0; i<xtopol_pc; i++)
    tmpf += sqrt( (v[i].x-ax)*(v[i].x-ax) + (v[i].y-ay)*(v[i].y-ay) + (v[i].z-az)*(v[i].z-az) );
  tmpf/=(double)xtopol_pc;

  for(i=0; i<xtopol_pc; i++) {
      v[i].x=(v[i].x-ax)/tmpf;;
      v[i].y=(v[i].y-ay)/tmpf;;
      v[i].z=(v[i].z-az)/tmpf;;
  }
  xtopol_output(v,xtopol_pc,argv[2]);
}

// output in brain-view .obj format
if (output_bic) {
  sprintf(filename,"%s.obj",argv[2]);
  FILE * fp = fopen( filename, "w" );
  fprintf( fp, "P 0.3 0.3 0.4 10 1 %d\n", pc );
  for(i=0;i<pc;i++){
    fprintf( fp, "%f %f %f\n", v[i].x, v[i].y, v[i].z );
  }
  fprintf( fp, "\n" );

  int ntri = 0;
  for(i=0;i<pc;i++){
    double nx, ny, nz, tmpf;
    int k, l;
    nx=ny=nz=0.0;
    for(k=0; v[i].n[k]>-1; k++);
    ntri += k;
    for(l=0; l<k; l++) {
      double adx = v[v[i].n[l]].x - v[i].x,
             ady = v[v[i].n[l]].y - v[i].y,
             adz = v[v[i].n[l]].z - v[i].z,
             bdx = v[v[i].n[(l+1)%k]].x - v[i].x,
             bdy = v[v[i].n[(l+1)%k]].y - v[i].y,
             bdz = v[v[i].n[(l+1)%k]].z - v[i].z;

      nx += ady*bdz - adz*bdy;
      ny += adz*bdx - adx*bdz;
      nz += adx*bdy - ady*bdx;
    }
    tmpf = sqrt(nx*nx+ny*ny+nz*nz);
    nx/=(double)tmpf; ny/=(double)tmpf; nz/=(double)tmpf;
    fprintf( fp, "%f %f %f\n", nx, ny, nz );
  }

  ntri /= 3;
  fprintf( fp, "%d\n", ntri );
  fprintf( fp, "0 1 1 1 1\n" );
  for( i = 1; i <= ntri; i++ ) {
    fprintf( fp, " %d", 3*i );
    if( i%8 == 0 ) fprintf( fp, "\n" );
  }

  int count = 0;
  for(i=0;i<pc;i++){
    int k, l;
    for(k=0; v[i].n[k]>-1; k++);
    for(l=0; l<k; l++) {
      int n1 = v[i].n[l];
      int n2 = v[i].n[(l+1)%k];
      if( i < n1 && i < n2 ) {
        fprintf( fp, " %d", i ); count++; if( count%8 == 0 ) fprintf( fp, "\n" );
        fprintf( fp, " %d", n1 ); count++; if( count%8 == 0 ) fprintf( fp, "\n" );
        fprintf( fp, " %d", n2 ); count++; if( count%8 == 0 ) fprintf( fp, "\n" );
      }
    }

  }

  fclose( fp );
}

  return(0);
}



/*  libss.h - basic C include stuff

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

#include "libssbase.h"

typedef struct {
  FDT *i;                                    /* image data */
  int x, y, z, t,                            /* image dimensions */
    dt, bpv;                                 /* datatype and bytes-per-voxel */
  double xv, yv, zv,                         /* voxel size (always all +ve) */
    xv0, yv0, zv0,                           /* voxel size as set in original header */
    xo, yo, zo,                              /* origin */
    tr,                                      /* seconds between volumes */
    min, max, thresh2, thresh98, thresh,     /* various brightness thresholds */
    medianval,
    lthresh, uthresh,                        /* lower and upper thresholds - some procs will ignore data outside these */
    dtmin, dtmax;                            /* min and max values for data type */
  char orient,                               /* orientation */
    lut[24];                                 /* lut file */
  float info;                                /* slice order information */
} image_struct;  

FDT median(double, FDT*, int);
double dmedian(double, double*, int);
FDT mean(FDT*, int);
FDT TLI(image_struct, double, double, double);
int make_isometric (image_struct, image_struct*);
int p2c(int, int, int, double, double, double, double, double, double, double*, double*, double*);
FDT getp2c(image_struct, double, double, double, double, double, double);
int find_histogram (image_struct*, int*, int, float);
/*void find_roi_histogram (image_struct*, int, int, int, int, int, int, int*, int);*/
void find_thresholds (image_struct*);
void c_of_g (image_struct, double*, double*, double*);
void invert_y (image_struct);
double find_radius (image_struct, double);
void find_ellipsoid( image_struct, double *, double *, double * );
void print_image_struct(image_struct);
void init_image_struct(image_struct*);
void tempfilt(image_struct, double, double);
void spatgauss(FDT*, int, int, int, double, double, double);
void spatfilt(image_struct, double);
void intnorm(image_struct, int, double);
void fmrimask(image_struct, image_struct);
void resample(image_struct, image_struct*, float);
void sample(image_struct, image_struct, float);

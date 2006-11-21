
/*  libssbase.h - basic C include stuff

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>

#ifndef MAXSHORT
#define MAXSHORT SHRT_MAX
#endif

#ifndef MAXINT
#define MAXINT INT_MAX
#endif

#ifndef FDT
#define FDT unsigned char
#endif

/* some weird stuff to expand the data type into a string */
#define FDTSSS(s) #s
#define FDTSS(s) FDTSSS(s)
#define FDTS FDTSS(FDT)

#ifndef MAX
#define MAX(a,b)   ( (a) > (b) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b)   ( (a) < (b) ? (a) : (b) )
#endif

#define SGN(a)  ( (a) > 0   ? (1) : (-1) )
#define FTOI(a) ( (a) < 0 ? ((int)(a-0.5)) : ((int)(a+0.5)) )
#define ABS(a)  ( (a) > 0   ? (a) : (-(a)) )

/* a cheap and grotty way of indexing 3/4D data - sometimes useful */
#define IA(i,x,y,z)   ( i[((z)*y_size+(y))*x_size+(x)] )
#define IT(i,x,y,z,t) ( i[(((t)*z_size+(z))*y_size+(y))*x_size+(x)] )

#if !defined(M_PI)
# define M_PI 3.14159265358979323846
#endif

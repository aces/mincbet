/* {{{ Copyright */

/*  libpicio - collection of 2D image io routines

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

/* }}} */
/* {{{ includes */

#include "libss.h"

/* }}} */
/* {{{ write_pgm */

int write_pgm ( char *filename, int x_size, int y_size,	unsigned char *i )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"w"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P5\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      fwrite(&i[y*x_size+x],1,1,ofp);

  fclose(ofp);
  return(0);
}

/* }}} */
/* {{{ write_ppm */

int write_ppm ( char *filename, int x_size, int y_size,
		unsigned char *r, unsigned char *g, unsigned char *b )
{
  FILE *ofp;
  int  x, y;

  if ((ofp=fopen(filename,"w"))==NULL)
    {
      printf("Can't open %s for writing\n",filename);
      return(1);
    }

  fprintf(ofp,"P6\n");
  fprintf(ofp,"%d %d\n",x_size,y_size);
  fprintf(ofp,"255\n");

  for(y=0; y<y_size; y++)
    for(x=0; x<x_size; x++)
      {
	fwrite(&r[y*x_size+x],1,1,ofp);
	fwrite(&g[y*x_size+x],1,1,ofp);
	fwrite(&b[y*x_size+x],1,1,ofp);
      }

  fclose(ofp);
  return(0);
}

/* }}} */

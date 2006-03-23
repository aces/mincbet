/*
 * Analyze IO stuff, originally based on code by Mayo:
 *
 * (c) Copyright, 1986-1995
 * Biomedical Imaging Resource
 * Mayo Foundation
 *
 */

/* {{{ includes and typedefs */

#include "libss.h"
#include "libavw.h"

/* }}} */
/* {{{ swap byte ordering */

/* {{{ swap_4 */

void swap_4(unsigned char *pntr)
{
  unsigned char tmp;

  tmp=*pntr;
  *pntr=*(pntr+3);
  *(pntr+3)=tmp;

  tmp=*(pntr+1);
  *(pntr+1)=*(pntr+2);
  *(pntr+2)=tmp;
}

/* }}} */
/* {{{ swap_2 */

void swap_2(unsigned char *pntr)
{
  unsigned char tmp;
  
  tmp=*pntr;
  *pntr=*(pntr+1);
  *(pntr+1)=tmp;
}

/* }}} */
/* {{{ swap_hdr */

void swap_hdr(struct dsr *pntr)
{
  swap_2((unsigned char *)&pntr->hist.originator[0]) ;
  swap_2((unsigned char *)&pntr->hist.originator[2]) ;
  swap_2((unsigned char *)&pntr->hist.originator[4]) ;
  swap_2((unsigned char *)&pntr->hist.originator[6]) ;
  swap_2((unsigned char *)&pntr->hist.originator[8]) ;
  swap_4((unsigned char *)&pntr->hk.sizeof_hdr) ;
  swap_4((unsigned char *)&pntr->hk.extents) ;
  swap_2((unsigned char *)&pntr->hk.session_error) ;
  swap_2((unsigned char *)&pntr->dime.dim[0]) ;
  swap_2((unsigned char *)&pntr->dime.dim[1]) ;
  swap_2((unsigned char *)&pntr->dime.dim[2]) ;
  swap_2((unsigned char *)&pntr->dime.dim[3]) ;
  swap_2((unsigned char *)&pntr->dime.dim[4]) ;
  swap_2((unsigned char *)&pntr->dime.dim[5]) ;
  swap_2((unsigned char *)&pntr->dime.dim[6]) ;
  swap_2((unsigned char *)&pntr->dime.dim[7]) ;
  swap_2((unsigned char *)&pntr->dime.unused1) ;
  swap_2((unsigned char *)&pntr->dime.datatype) ;
  swap_2((unsigned char *)&pntr->dime.bitpix) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[0]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[1]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[2]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[3]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[4]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[5]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[6]) ;
  swap_4((unsigned char *)&pntr->dime.pixdim[7]) ;
  swap_4((unsigned char *)&pntr->dime.vox_offset) ;
  swap_4((unsigned char *)&pntr->dime.funused1) ;
  swap_4((unsigned char *)&pntr->dime.funused2) ;
  swap_4((unsigned char *)&pntr->dime.cal_max) ;
  swap_4((unsigned char *)&pntr->dime.cal_min) ;
  swap_4((unsigned char *)&pntr->dime.compressed) ;
  swap_4((unsigned char *)&pntr->dime.verified) ;
  swap_2((unsigned char *)&pntr->dime.dim_un0) ;
  swap_4((unsigned char *)&pntr->dime.glmax) ;
  swap_4((unsigned char *)&pntr->dime.glmin) ;
}

/* }}} */

/* }}} */
/* {{{ avw_read */

void avw_read(char *filename, image_struct *image)
{
/* {{{ vars */

char filename2[1000];

FILE *fd;
struct dsr hdr;
int swap=0, i, size;
static int DataTypeSizes[32] = {0,1,8,0,16,0,0,0,32,0,0,0,0,0,0,0,32};
void *tmpbuf;

/* }}} */

  /* {{{ read header */

if ( (strcmp(".hdr",filename+strlen(filename)-4)==0) ||
     (strcmp(".img",filename+strlen(filename)-4)==0) )
     filename[strlen(filename)-4]=0;

sprintf(filename2,"%s.hdr",filename);

if((fd=fopen(filename2,"r"))==NULL)
{
  fprintf(stderr,"Can't open %s\n",filename2);
  exit(1);
}

fread(&hdr,1,sizeof(struct dsr),fd);
fclose(fd);

/*if(hdr.dime.dim[0] < 0 || hdr.dime.dim[0] > 15)*/ /* old test for byte order */
if(hdr.hk.sizeof_hdr!=348) /* better test for byte order */
{
  /*fprintf(stderr,"Byte swapping\n");*/
  swap=1;
  swap_hdr(&hdr);
}

/* }}} */
  /* {{{ dimensions */

image->x=MAX(1,hdr.dime.dim[1]);
image->y=MAX(1,hdr.dime.dim[2]);
image->z=MAX(1,hdr.dime.dim[3]);
image->t=MAX(1,hdr.dime.dim[4]);

size=image->x *image->y *image->z *image->t;

/* }}} */
  /* {{{ voxel size */

image->xv0=hdr.dime.pixdim[1];
image->yv0=hdr.dime.pixdim[2];
image->zv0=hdr.dime.pixdim[3];

image->xv=fabs(image->xv0);
image->yv=fabs(image->yv0);
image->zv=fabs(image->zv0);

if (image->xv<0.00001) { image->xv=1; image->xv0=1; }
if (image->yv<0.00001) { image->yv=1; image->yv0=1; }
if (image->zv<0.00001) { image->zv=1; image->zv0=1; }

/* }}} */
  /* {{{ origin */

{
  short short_array[100];
  memcpy(short_array,hdr.hist.originator,5*sizeof(short));
  image->xo=short_array[0];
  image->yo=short_array[1];
  image->zo=short_array[2];
}

/* }}} */
  /* {{{ info, orient, min+max, lut */

image->info=hdr.dime.funused1;

image->orient=hdr.hist.orient;

/* image->min=hdr.dime.glmin; */ /* old method */
/* image->max=hdr.dime.glmax; */
image->min=hdr.dime.cal_min;
image->max=hdr.dime.cal_max;

strncpy(image->lut,hdr.hist.aux_file,24);
image->lut[23]=0;

/* }}} */
  /* {{{ datatype and image mallocs */

image->dt=hdr.dime.datatype;
image->bpv=DataTypeSizes[image->dt]/8;

tmpbuf=malloc(image->bpv * size);

if (sizeof(FDT)>image->bpv)
  image->i=(FDT*)malloc(sizeof(FDT) * size); /* will need more memory than the file size */
else
  image->i=(FDT*)tmpbuf;

/* }}} */
  /* {{{ read data */

sprintf(filename2,"%s.img",filename);

if((fd=fopen(filename2,"r"))==NULL)
{
  fprintf(stderr,"Can't open %s\n",filename2);
  exit(1);
}

if ( fread(tmpbuf,image->bpv,size,fd) < size )
{
  fprintf(stderr,"Not enough data in %s\n",filename2);
  exit(1);
}

fclose(fd);

/* }}} */
  /* {{{ swap byte order if required */

if ( swap && ( (image->bpv==2) || (image->bpv==4) ) )
  for (i=0; i<size; i++)
  {
    if (image->bpv==2)
      swap_2((unsigned char *)tmpbuf+i*2);

    if (image->bpv==4)
      swap_4((unsigned char *)tmpbuf+i*4);
  }

/* }}} */
  /* {{{ setup dtmin and dtmax and convert type if required */

switch (image->dt) {
 case DT_UNSIGNED_CHAR:
   if (strcmp(FDTS,"unsigned char")!=0) { for(i=0;i<size;i++) image->i[i] = ((unsigned char*)tmpbuf)[i]; }
   break;
 case DT_SIGNED_SHORT:
   if (strcmp(FDTS,"signed short")!=0)  { for(i=0;i<size;i++) image->i[i] = ((short*)tmpbuf)[i]; }
   break;
 case DT_SIGNED_INT:
   if (strcmp(FDTS,"signed int")!=0)    { for(i=0;i<size;i++) image->i[i] = ((int*)tmpbuf)[i]; }
   break;
 case DT_FLOAT:
   if (strcmp(FDTS,"float")!=0)         { for(i=0;i<size;i++) image->i[i] = ((float*)tmpbuf)[i]; }
   break;
 default:
   { fprintf(stderr,"Wrong input data type for program\n"); exit(1); }
   break;
}

/* set correct data type - needed in case type was changed */
if (strcmp(FDTS,"unsigned char")==0) { image->dt=DT_UNSIGNED_CHAR; image->dtmin=0;           image->dtmax=255;      }
if (strcmp(FDTS,"signed short")==0)  { image->dt=DT_SIGNED_SHORT;  image->dtmin=-MAXSHORT-1; image->dtmax=MAXSHORT; }
if (strcmp(FDTS,"signed int")==0)    { image->dt=DT_SIGNED_INT;    image->dtmin=-MAXINT-1;   image->dtmax=MAXINT;   }
if (strcmp(FDTS,"float")==0)         { image->dt=DT_FLOAT;         image->dtmin=-1e10;       image->dtmax=1e10;     }
image->bpv=DataTypeSizes[image->dt]/8;

/* }}} */
  /* {{{ setup image_struct stuff that's not set by AVW */

image->tr=3;                     /* well it's possible.... */
image->thresh2=image->min;       /* better than nothing until it gets set properly using find_thresholds */
image->thresh98=image->max;      /* ditto */
image->thresh=image->min;        /* can use as a flag for whether find_thresholds has been run */
image->lthresh=image->dtmin;     /* relevant procs know not to use these unless changed from these defaults */
image->uthresh=image->dtmax;

/* }}} */
}

/* }}} */
/* {{{ avw_write */

void avw_write(char *filename, image_struct image)
{
/* {{{ vars */

FILE *fd;
char filename2[1000];
struct dsr hdr;
static int DataTypeSizes[32] = {0,1,8,0,16,0,0,0,32,0,0,0,0,0,0,0,32};

/* }}} */

  /* {{{ copy image info */

/* {{{ set uninteresting defaults */

hdr.hk.sizeof_hdr=348;
hdr.hk.data_type[0]=0;
hdr.hk.db_name[0]=0;
hdr.hk.extents=0;
hdr.hk.session_error=0;
hdr.hk.regular='r';
hdr.hk.hkey_un0=0;

sprintf(hdr.dime.vox_units,"mm");
hdr.dime.cal_units[0]=0;
hdr.dime.unused1=0;
hdr.dime.dim_un0=0;
hdr.dime.vox_offset=0;
hdr.dime.funused2=0;
hdr.dime.funused3=0;
hdr.dime.compressed=0;
hdr.dime.verified=0;

hdr.hist.descrip[0]=0;
hdr.hist.generated[0]=0;
hdr.hist.scannum[0]=0;
hdr.hist.patient_id[0]=0;
hdr.hist.exp_date[0]=0;
hdr.hist.exp_time[0]=0;
hdr.hist.hist_un0[0]=0;
hdr.hist.views=0;
hdr.hist.vols_added=0;
hdr.hist.start_field=0;
hdr.hist.field_skip=0;
hdr.hist.omax=0;
hdr.hist.omin=0;
hdr.hist.smax=0;
hdr.hist.smin=0;

/* }}} */

hdr.dime.dim[0]=4;
hdr.dime.dim[1]=image.x;
hdr.dime.dim[2]=image.y;
hdr.dime.dim[3]=image.z;
hdr.dime.dim[4]=image.t;
hdr.dime.dim[5]=0;
hdr.dime.dim[6]=0;
hdr.dime.dim[7]=0;

hdr.dime.pixdim[0]=0;
hdr.dime.pixdim[1]=image.xv0;
hdr.dime.pixdim[2]=image.yv0;
hdr.dime.pixdim[3]=image.zv0;
hdr.dime.pixdim[4]=0;
hdr.dime.pixdim[5]=0;
hdr.dime.pixdim[6]=0;
hdr.dime.pixdim[7]=0;

hdr.dime.funused1=image.info;

{
  short short_array[100];
  short_array[0]=image.xo;
  short_array[1]=image.yo;
  short_array[2]=image.zo;
  memcpy(hdr.hist.originator,short_array,5*sizeof(short));
}

hdr.hist.orient=image.orient;

strncpy(hdr.hist.aux_file,image.lut,24);
hdr.hist.aux_file[23]=0;

hdr.dime.glmin=image.min;   /* it is normally these that are used for display range */
hdr.dime.glmax=image.max;
hdr.dime.cal_min=image.min; /* but these are float so are needed for some images */
hdr.dime.cal_max=image.max;

hdr.dime.datatype=image.dt;
hdr.dime.bitpix=DataTypeSizes[image.dt];

/* }}} */
  /* {{{ write header */

if ( (strcmp(".hdr",filename+strlen(filename)-4)==0) ||
     (strcmp(".img",filename+strlen(filename)-4)==0) )
     filename[strlen(filename)-4]=0;

sprintf(filename2,"%s.hdr",filename);

if((fd=fopen(filename2,"w"))==NULL)
{
  fprintf(stderr,"Can't open %s\n",filename2);
  exit(1);
}

if ( fwrite(&hdr,1,sizeof(struct dsr),fd) < 1 )
{
  fprintf(stderr,"Couldn't write %s\n",filename2);
  exit(1);
}

fclose(fd);

/* }}} */
  /* {{{ write image */

sprintf(filename2,"%s.img",filename);

if((fd=fopen(filename2,"w"))==NULL)
{
  fprintf(stderr,"Can't open %s\n",filename2);
  exit(1);
}

if ( fwrite(image.i,hdr.dime.bitpix/8,image.x * image.y * image.z * image.t,fd) < image.x * image.y * image.z * image.t )
{
  fprintf(stderr,"Couldn't write %s\n",filename2);
  exit(1);
}

fclose(fd);

/* }}} */
}

/* }}} */

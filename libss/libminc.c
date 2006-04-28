/*
 * Andrew Janke - a.janke@gmail.com
 * Copy and munge as you choose, heck even cut this 
 * line out if you feel that way inclined
 *
 */

#include "libss.h"
#include "libavw.h"
#include <volume_io.h>

void minc_read(char *filename, image_struct * image)
{
   int      i, size, x, y, z;
   static int DataTypeSizes[32] = { 0, 1, 8, 0, 16, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 32 };

   /* minc vars */
   Volume   volume;
   char    *axis_order[3] = { MIzspace, MIyspace, MIxspace };
   Real     real_min, real_max;
   Real     voxel_min, voxel_max;
   Real     steps[MAX_VAR_DIMS];
   Real     starts[MAX_VAR_DIMS];
   int      sizes[MAX_VAR_DIMS];

   /* read in the volume */
   input_volume(filename, MAX_VAR_DIMS, axis_order,
                NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, &volume, NULL);

   get_volume_real_range(volume, &real_min, &real_max);
   get_volume_voxel_range(volume, &voxel_min, &voxel_max);
   get_volume_sizes(volume, sizes);
   get_volume_starts(volume, starts);
   get_volume_separations(volume, steps);

   image->x = MAX(1, sizes[2]);
   image->y = MAX(1, sizes[1]);
   image->z = MAX(1, sizes[0]);
   image->t = 1;

   size = image->x * image->y * image->z * image->t;

   /* voxel size */

   image->xv0 = steps[2];
   image->yv0 = steps[1];
   image->zv0 = steps[0];

   image->xv = fabs(image->xv0);
   image->yv = fabs(image->yv0);
   image->zv = fabs(image->zv0);

   if(image->xv < 0.00001){
      image->xv = 1;
      image->xv0 = 1;
      }
   if(image->yv < 0.00001){
      image->yv = 1;
      image->yv0 = 1;
      }
   if(image->zv < 0.00001){
      image->zv = 1;
      image->zv0 = 1;
      }

   /* origin */
   image->xo = starts[2];
   image->yo = starts[1];
   image->zo = starts[0];

   /* info, orient, min+max, lut */
   image->info = 0;

   /* this is set by the axis_order so it will be this (always.) */
   image->orient = 0;

   image->min = real_min;
   image->max = real_max;
   
   /* in minc a lut is never associated with a file */
   strncpy(image->lut, "                        ", 24);
   image->lut[23] = 0;

   /* datatype and image mallocs */
   /* let's do everything in float */
   image->dt = 16;
   image->bpv = DataTypeSizes[image->dt] / 8;

   /* malloc memory */
   image->i = (FDT *) malloc(sizeof(FDT) * size);

   /* read data */
   i = 0;
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){

            image->i[i] = (FDT)get_volume_real_value(volume, z, y, x, 0, 0);
            i++;
            }
         }
      }
   

   /* setup dtmin and dtmax and convert type if required */
  /* set correct data type - needed in case type was changed */
   if(strcmp(FDTS, "unsigned char") == 0){
      image->dt = DT_UNSIGNED_CHAR;
      image->dtmin = 0;
      image->dtmax = 255;
      }
   if(strcmp(FDTS, "signed short") == 0){
      image->dt = DT_SIGNED_SHORT;
      image->dtmin = -MAXSHORT - 1;
      image->dtmax = MAXSHORT;
      }
   if(strcmp(FDTS, "signed int") == 0){
      image->dt = DT_SIGNED_INT;
      image->dtmin = -MAXINT - 1;
      image->dtmax = MAXINT;
      }
   if(strcmp(FDTS, "float") == 0){
      image->dt = DT_FLOAT;
      image->dtmin = -1e10;
      image->dtmax = 1e10;
      }
   image->bpv = DataTypeSizes[image->dt] / 8;

   /* setup some other image_struct stuff */
   image->tr = 3;                      /* well it's possible.... */
   image->thresh2 = image->min;        /* better than nothing until it gets set properly using find_thresholds */
   image->thresh98 = image->max;       /* ditto */
   image->thresh = image->min;         /* can use as a flag for whether find_thresholds has been run */
   image->lthresh = image->dtmin;      /* relevant procs know not to use these unless changed from these defaults */
   image->uthresh = image->dtmax;

   }

/* minc_write */
/* subtle dig at FSL developers to switch to minc                                        */
/* ----------------------------------------------                                        */
/* this is a bit snarfling as we loose most of the good volume info darn internal format */
/* that doesn't hold goodies like direction cosines, history, etc etc etc                */
/* besides it's good enough for what it was written for. BET                             */
void minc_write(char *filename, image_struct image)
{
   char filename2[1000];
   Volume volume;
   int x,y,z,i;
   char    *axis_order[3] = { MIzspace, MIyspace, MIxspace };
   
   nc_type datatype;
   BOOLEAN data_signed;
   Real     steps[MAX_VAR_DIMS];
   Real     starts[MAX_VAR_DIMS];
   int      sizes[MAX_VAR_DIMS];
   
   sprintf(filename2, "%s.mnc", filename);
   
   /* figger out the datatype */
   switch(image.dt){
      default: 
      case DT_UNSIGNED_CHAR:
         datatype = NC_BYTE;
         data_signed = FALSE;
         break;
      case DT_SIGNED_SHORT:
         datatype = NC_SHORT;
         data_signed = TRUE;
         break;
      case DT_FLOAT:
         datatype = NC_FLOAT;
         data_signed = TRUE;
         break;
      }
   
   volume = create_volume(3, axis_order, datatype, data_signed, image.dtmin, image.dtmax);
   set_volume_real_range(volume, image.min, image.max);
   
   sizes[2] = image.x;
   sizes[1] = image.y;
   sizes[0] = image.z;
 
   steps[2] = image.xv0;
   steps[1] = image.yv0;
   steps[0] = image.zv0;
   
   starts[2] = image.xo;
   starts[1] = image.yo;
   starts[0] = image.zo;
   
   set_volume_sizes(volume, sizes);
   set_volume_starts(volume, starts);
   set_volume_separations(volume, steps);
   alloc_volume_data(volume);
   
   /* write data */
   i = 0;
   for(z = 0; z < sizes[0]; z++){
      for(y = 0; y < sizes[1]; y++){
         for(x = 0; x < sizes[2]; x++){

            set_volume_real_value(volume, z, y, x, 0, 0, (Real)image.i[i]);
            i++;
            }
         }
      }
   
   output_volume(filename2, NC_UNSPECIFIED, TRUE, 0, 0, volume, NULL, (minc_output_options*)NULL);
   
   }


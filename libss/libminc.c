/*
 * Andrew Janke - rotor@cmr.uq.edu.au
 * Copy and munge as you choose, heck even cut this 
 * line out if you feel that way inclined
 *
 */

#include "libss.h"
#include "libavw.h"
#include "libminc.h"
#include <volume_io.h>

void minc_read(char *filename, image_struct * image)
{
   int      i, size, x, y, z;
   static int DataTypeSizes[32] = { 0, 1, 8, 0, 16, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 32 };

   enum     {XX, YY, ZZ};

   /* minc vars */
   Volume   volume;
   Real     real_min, real_max;
   Real     steps[MAX_VAR_DIMS];
   Real     starts[MAX_VAR_DIMS];
   int      sizes[MAX_VAR_DIMS];

   /* read in the volume */
   input_volume(filename, MAX_VAR_DIMS, NULL, 
                NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, &volume, NULL);

   /* obtain order of axes */
   char ** dimnames = create_output_dim_names( volume, filename, 
                                               NULL, sizes );
   delete_volume( volume );

   /* re-open volume with original order of axes,
      not the default (z,y,x) one. */
   input_volume( filename, MAX_VAR_DIMS, dimnames,
                 NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, &volume, NULL);

   get_volume_real_range(volume, &real_min, &real_max);

   get_volume_sizes(volume, sizes);
   get_volume_starts(volume, starts);
   get_volume_separations(volume, steps);

   image->x = MAX(1, sizes[XX]);
   image->y = MAX(1, sizes[YY]);
   image->z = MAX(1, sizes[ZZ]);
   image->t = 1;

   size = image->x * image->y * image->z * image->t;

   /* voxel size */

   image->xv0 = steps[XX];
   image->yv0 = steps[YY];
   image->zv0 = steps[ZZ];

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
   image->xo = starts[XX];
   image->yo = starts[YY];
   image->zo = starts[ZZ];

   /* info, orient, min+max, lut */
   image->info = 0;
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
   for(z = 0; z < image->z; z++){
      for(y = 0; y < image->y; y++){
         for(x = 0; x < image->x; x++){
            image->i[i] = (FDT)get_volume_real_value(volume, x, y, z, 0, 0);
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

   /* Delete volume after usage */
   delete_volume( volume );

   }

/* minc_write */
/* subtle dig at FSL developers to switch to minc                                        */
/* ----------------------------------------------                                        */
/* this is a bit snarfling as we loose most of the good volume info darn internal format */
/* that doesn't hold goodies like direction cosines, history, etc etc etc                */
/* besides it's good enough for what it was written for. BET                             */

void minc_write(char *input_filename, int dt, char *output_filename, image_struct image) {

   char filename2[1000];
   Volume input_vol;
   Volume volume;
   int x,y,z,i;
   int      sizes[MAX_VAR_DIMS];
   
   nc_type datatype;
   BOOLEAN data_signed;
   Real    dtmin, dtmax;
   
   /* read in the volume */
   input_volume(input_filename, MAX_VAR_DIMS, NULL,
                NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, &input_vol, NULL);

   /* obtain order of axes */
   char ** dimnames = create_output_dim_names( input_vol, input_filename, 
                                               NULL, sizes );
   delete_volume( input_vol );

   /* re-open volume with original order of axes, not the
      default (z,y,x) one. Now the output volume will be
      written in the same order as the input volume. */
   input_volume(input_filename, MAX_VAR_DIMS, dimnames,
                NC_UNSPECIFIED, TRUE, 0.0, 0.0, TRUE, &input_vol, NULL);

   sprintf(filename2, "%s.mnc", output_filename);

   switch( dt ) {
     // case default:
     case MINCBET_SAME:
       datatype = NC_UNSPECIFIED;
       data_signed = FALSE;
       dtmin = 0;
       dtmax = 0;
       break;
     case MINCBET_BYTE:
       datatype = NC_BYTE;
       data_signed = FALSE;
       dtmin = 0;
       dtmax = 255;
       break;
     case MINCBET_SHORT:
       datatype = NC_SHORT;
       data_signed = TRUE;
       dtmin = -MAXSHORT - 1;
       dtmax = MAXSHORT;
       break;
     case MINCBET_FLOAT:
       datatype = NC_FLOAT;
       data_signed = TRUE;
       dtmin = -1e10;
       dtmax = 1e10;
       break;
   }

   // Copy the header definition from the original volume.

   volume = copy_volume_definition( input_vol, datatype, data_signed, 
                                    dtmin, dtmax );
   delete_volume( input_vol );
   
   set_volume_real_range(volume, image.min, image.max);

   /* write data */
   i = 0;
   for(z = 0; z < image.z; z++){
     for(y = 0; y < image.y; y++){
       for(x = 0; x < image.x; x++){
         set_volume_real_value(volume, x, y, z, 0, 0, (Real)image.i[i]);
         i++;
       }
     }
   }
   
   output_volume(filename2, NC_UNSPECIFIED, TRUE, 0, 0, volume, 
                 NULL, (minc_output_options*)NULL);
   
}


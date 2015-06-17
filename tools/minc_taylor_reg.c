/* ----------------------------- MNI Header -----------------------------------
@NAME       :  minc_taylor_reg
@DESCRIPTION:  Taylor regularization of longitudinal deformation field
@COPYRIGHT  :
              Copyright 2015 Nicolas Guizard, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include <math.h>
#include <volume_io.h>
#include <time_stamp.h> 

int  main(
    int   argc,
    char  *argv[] )
{
    VIO_Volume     vec0,vecdx,vecdy,vecdz;
    VIO_Status     status;
    int        x, y, z, sizes[VIO_N_DIMENSIONS];
    char       *vec0_filename,*vecdx_filename,*vecdy_filename,*vecdz_filename, *output_filename, *history;
    float     value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,value11,value12,value13,new_voxel,min,max;
    VIO_BOOL    thresholding;

    if( argc < 3 )
    {
        print( "Usage: %s  vec0.mnc vecx.mnc vecy.mnc vecz.mnc output.mnc\n",
               argv[0] );
        return( 1 );
    }

    vec0_filename = argv[1];
    vecdx_filename = argv[2];
    vecdy_filename = argv[3];
    vecdz_filename = argv[4];
    output_filename = argv[5];

    status = input_volume( vec0_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &vec0, (minc_input_options *) NULL ) ;

    if( status != VIO_OK )
        return( 1 );
    
    status = input_volume( vecdx_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &vecdx, (minc_input_options *) NULL ) ;

    if( status != VIO_OK )
        return( 1 );

     status = input_volume( vecdy_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &vecdy, (minc_input_options *) NULL ) ;

    if( status != VIO_OK )
        return( 1 );
    
     status = input_volume( vecdz_filename, 3, File_order_dimension_names,
                      NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                      TRUE, &vecdz, (minc_input_options *) NULL ) ;

    if( status != VIO_OK )
        return( 1 );
    
    /* --- convert new_value to voxel */    
	
    get_volume_sizes( vec0, sizes );
    min=256;
    max=-256;

        for_less( x, 1, sizes[VIO_X]-1 )
        {
            for_less( y, 1, sizes[VIO_Y]-1 )
            {
                for_less( z, 1, sizes[VIO_Z]-1 )
                {
                    GET_VALUE_3D( value1, vec0, x, y, z );
                    GET_VALUE_3D( value2, vec0, x-1, y, z );
                    GET_VALUE_3D( value3, vecdx, x-1, y, z );
                    GET_VALUE_3D( value4, vec0, x+1, y, z );
                    GET_VALUE_3D( value5, vecdx, x+1, y, z );
                    GET_VALUE_3D( value6, vec0, x, y-1, z );
                    GET_VALUE_3D( value7, vecdy, x, y-1, z );
                    GET_VALUE_3D( value8, vec0, x, y+1, z );
                    GET_VALUE_3D( value9, vecdy, x, y+1, z );
                    GET_VALUE_3D( value10, vec0, x, y, z-1 );
                    GET_VALUE_3D( value11, vecdz, x, y, z-1 );
                    GET_VALUE_3D( value12, vec0, x, y, z+1 );		    
                    GET_VALUE_3D( value13, vecdz, x, y, z+1 );		    
                    
                    new_voxel =  (value2+value4+value6+value8+value10+value12+value3-value5+value7-value9+value11-value13)/6; 
                    
                    if(new_voxel<min)
                      min = new_voxel;
                    if(new_voxel>max)
                      max = new_voxel;
                    
                    SET_VOXEL_3D( vec0, x, y, z, new_voxel );                    		    
                    //SET_VOXEL_3D( vec0, x, y, z, value1 );                    		    
                }
            }
        }
    

    history = "Taylor regularization of longitudinal deformation field";

    status = output_volume( output_filename, NC_FLOAT, FALSE, min, max,
                            vec0, history, (minc_output_options *) NULL );

    return( 0 );
}


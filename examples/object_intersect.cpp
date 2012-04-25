/* ----------------------------- MNI Header -----------------------------------
@NAME       : 
@DESCRIPTION: 
@COPYRIGHT  :
              Copyright 2006 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include "minc_wrappers.h"
//#include "minc_io.h"
#include <itkBSplineInterpolateImageFunction.h>


#include <getopt.h>
#include  <bicpl.h>
#include <fstream>

using namespace  std;
using namespace  minc;

// typedef minc::mincVectorBSplineInterpolate<minc::def3d,double> Vector_Interpolator;
// typedef minc::VectorMagnitudeFilter<minc::def3d,minc::image3d> _MagFilterType;


typedef itk::BSplineInterpolateImageFunction< minc::image3d, double, double >  InterpolatorType;


void show_usage (const char * prog)
{
  std::cerr<<"Usage:"<<prog<<" in_volume.mnc in.obj out.csv "<<std::endl
      <<"[ "<<std::endl
      <<"  --clobber clobber output file(s)" <<std::endl
      <<"  --order <n> interpolation order" <<std::endl
      <<"]"<<std::endl; 
}

int main (int argc, char **argv)
{
  int verbose=0,clobber=0;
  int order=2;
  
  static struct option long_options[] = {
    {"verbose", no_argument,       &verbose, 1},
    {"quiet",   no_argument,       &verbose, 0},
    {"clobber", no_argument,       &clobber, 1},
    {"order", required_argument,       0, 'o'},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vo:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
      case 0:
        break;
      case 'o':
        order=atoi(optarg);break;
      case 'v':
        cout << "Version: 0.1" << endl;
        return 0;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }

  if((argc - optind) < 3) {
    show_usage (argv[0]);
    return 1;
  }
  std::string in_volume=argv[optind];
  std::string in_obj=argv[optind+1];
  std::string out_csv=argv[optind+2];
  
  
  if (!clobber && !access(out_csv.c_str(), F_OK))
  {
    cerr << out_csv.c_str() << " Exists!" << endl;
    return 1;
  }
  
  try
  {
    
    itk::ObjectFactoryBase::RegisterFactory(itk::MincImageIOFactory::New());
    itk::ImageFileReader<minc::image3d >::Pointer reader = itk::ImageFileReader<minc::image3d >::New();
    
    //initializing the reader
    reader->SetFileName(in_volume.c_str());
    reader->Update();
    
    minc::image3d::Pointer in=reader->GetOutput();
    
    std::cout<<"Building BSpline interpolator , order="<<order<<std::endl;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetSplineOrder(order);
    interpolator->SetInputImage(in);
    
    std::ofstream output(out_csv.c_str());
    
    File_formats         format;
    object_struct        **object_list;
    int n_objects=0;
    
    //bicpl sucks!
    if( input_graphics_file( (char*)in_obj.c_str(), &format, &n_objects,
        &object_list ) != OK )
    {
      std::cerr << " Error reading "<<in_obj.c_str() << std::endl;
      return 1;
    }
    
    //const int segment_length=10;
    
    std::vector<minc::tag_point> lines;
    std::vector<int> line_index;
    std::cout<<"Processing "<<n_objects<<" objects"<<std::endl;
    output<<"x,y,z,v"<<std::endl;
    for(int i=0;i< n_objects;i++ )
    {
      if( get_object_type( object_list[i] ) == POLYGONS )
      {
        polygons_struct      *polygons;
        polygons = get_polygons_ptr(object_list[i]);
        
        /*object_struct * lobj= create_object( LINES );
        lines_struct * lines = get_lines_ptr(lobj);
        
        initialize_lines( lines, WHITE );
        lines->n_points = polygons->n_items*segment_length;
        ALLOC( lines->points, lines->n_points );*/
        
        for(int pnt=0;pnt< polygons->n_points;pnt++ )
        {
          Point p_=polygons->points[pnt];//seeding point 
          minc::tag_point p_orig,p;
          p_orig[0]=p_.coords[0];p_orig[1]=p_.coords[1];p_orig[2]=p_.coords[2];
          //line_index.push_back(lines.size());
          //lines.push_back(p_orig);
          p=p_orig;
          
          double _in=interpolator->Evaluate(p);
          
/*          polygons->points[pnt].coords[0]=p[0];
          polygons->points[pnt].coords[1]=p[1];
          polygons->points[pnt].coords[2]=p[2];*/
          //int size = GET_OBJECT_SIZE( *polygons, poly );
          //if(size<3) continue; //?
          
/*          p1=polygons->points[POINT_INDEX( polygons->end_indices, poly, 0 )];
          p2=polygons->points[POINT_INDEX( polygons->end_indices, poly, 1 )];
          p3=polygons->points[POINT_INDEX( polygons->end_indices, poly, 2 )];*/
          /*
          std::cout<<poly<<"\t"<<POINT_INDEX( polygons->end_indices, poly, 0 )<<" "
              <<POINT_INDEX( polygons->end_indices, poly, 1 )<<" "
              <<POINT_INDEX( polygons->end_indices, poly, 2 )<<std::endl;*/
          //std::cout<<poly<<"\t"<<p1.coords[0]<<" "<<p1.coords[1]<<" "<<p1.coords[2]<<std::endl;
          //center
          output<<p[0]<<","<<p[1]<<","<<p[2]<<","<<
                _in<<std::endl;
        }
      }
    } 
    //int status = output_graphics_file( (char*)out_objf.c_str(), format,n_objects, object_list );
    //free up memory
    delete_object_list( n_objects, object_list );
    //return( status != OK );
    /*
    //outputting lines
    FILE *fp=fopen(out_objf.c_str(),"w");

  // IL want to write in binary

    fprintf(fp,"L 0.1 %d\n",lines.size());

  //WRITING cordinates

    for (int i=0;i<lines.size();i++) {
      fprintf(fp,"%f %f %f\n",lines[i][0],lines[i][1],lines[i][2]);
    }

    fprintf(fp,"%d\n1 ",line_index.size()); //1 color per line

    for (int i=0;i<line_index.size();i++) { //random colors
      fprintf(fp,"%f %f %f 1.0\n",rand()/(double)RAND_MAX,rand()/(double)RAND_MAX,rand()/(double)RAND_MAX);
    }

    for (int i=0;i<line_index.size();i++) {
      
      if(i<(line_index.size()-1))
        fprintf(fp,"%d ",line_index[i+1]);
      else
        fprintf(fp,"%d ",lines.size());//last line lasts untill the end
    }

    fprintf(fp,"\n");

    for (int i=0;i<line_index.size();i++) 
    {
      int end=i<(line_index.size()-1)?line_index[i+1]:lines.size();
      for(int j=line_index[i];j<end;j++)
        fprintf(fp," %d",j);
      
      fprintf(fp,"\n");
    }
    fclose(fp);*/
    
  } catch (const minc::generic_error & err) {
    cerr << "Got an error at:" << err.file () << ":" << err.line () << endl;
    return 1;
  }
  return 0;
};

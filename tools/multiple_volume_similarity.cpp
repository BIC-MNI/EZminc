/* ----------------------------- MNI Header -----------------------------------
@NAME       :  multiple_volume_similarity
@DESCRIPTION:  tool to calculate multiple volume discrete similarity
@COPYRIGHT  :
              Copyright 2011 Vladimir Fonov, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */

#include <iostream>
#include <fstream>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>

#include <getopt.h>
#include <map>
#include <set>

#if ( ITK_VERSION_MAJOR < 4 )
#include <time_stamp.h>    // for creating minc style history entry
#include <itkMincHelpers.h>
#include <itkMincImageIOFactory.h>
using namespace minc;
#else
//! allocate volume of the same dimension,spacing and origin
template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::ConstPointer &sample)
{
  image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
  image->SetBufferedRegion(sample->GetLargestPossibleRegion());
  image->SetRequestedRegion(sample->GetLargestPossibleRegion());
  image->SetSpacing( sample->GetSpacing() );
  image->SetOrigin ( sample->GetOrigin() );
  image->SetDirection(sample->GetDirection());
  image->Allocate();
}
template<class T,class S> void allocate_same(typename T::Pointer &image,const typename S::Pointer &sample)
{
  image->SetLargestPossibleRegion(sample->GetLargestPossibleRegion());
  image->SetBufferedRegion(sample->GetLargestPossibleRegion());
  image->SetRequestedRegion(sample->GetLargestPossibleRegion());
  image->SetSpacing( sample->GetSpacing() );
  image->SetOrigin ( sample->GetOrigin() );
  image->SetDirection(sample->GetDirection());
  image->Allocate();
}
#endif

typedef unsigned short  IOPixelType;
typedef unsigned char LabelPixelType;
typedef itk::Image<IOPixelType, 3>  IOImageType;
typedef itk::Image<LabelPixelType, 3>  LabelImageType;
typedef itk::Image<float, 3>  RealImageType;
typedef itk::VectorImage<float, 3>  VectorImageType;


typedef itk::ImageFileReader<LabelImageType>    LabelReaderType;
typedef itk::ImageFileReader<IOImageType>    ReaderType;
typedef itk::ImageFileWriter<IOImageType>    WriterType;
typedef itk::ImageFileWriter<RealImageType>  RealWriterType;

typedef itk::ImageRegionConstIterator<IOImageType>    IOImageConstIterator;
typedef itk::ImageRegionConstIterator<LabelImageType> LabelImageConstIterator;
typedef itk::ImageRegionConstIterator<VectorImageType> VectorImageConstIterator;

typedef itk::ImageRegionIterator<IOImageType> IOImageIterator;
typedef itk::ImageRegionIterator<LabelImageType> LabelImageIterator;
typedef itk::ImageRegionIterator<VectorImageType> VectorImageIterator;
typedef itk::ImageRegionIterator<RealImageType> RealImageIterator;



void show_usage (const char * prog)
{
  std::cout<<"Program calculates multiple volume similarity metrics for discrete labels "<<std::endl
           <<"or Generalized Tanimoto coefficient (GTC)" <<std::endl
           <<"based on :  William R. Crum, Oscar Camara, and Derek L. G. Hill"
           <<"\"Generalized Overlap Measures for Evaluation and Validation in Medical Image Analysis \""
           <<" IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 25, NO. 11, NOVEMBER 2006"<<std::endl
           <<"http://dx.doi.org/10.1109/TMI.2006.880587"<<std::endl<<std::endl
           <<"WARNING: Program will use a lot of memory , proportional to the number of classes"<<std::endl
           <<"Usage: "<<prog<<" <input1.mnc> <input2.mnc> ... <inputN.mnc>  "<<std::endl
           <<"[--verbose "<<std::endl
           <<" --bg include background in overlap "<<std::endl
           <<" --mask <mask.mnc>"<<std::endl
           <<" --classes <n>"<<std::endl
           <<" --list <file.list>"<<std::endl
           <<" --majority <output>"<<std::endl
           <<" --overlap <output>"<<std::endl
           <<" --relabel map.txt ]"<<std::endl;
}

int main(int argc,char **argv)
{
#if ( ITK_VERSION_MAJOR < 4 )
	char *history = time_stamp(argc, argv); 
	std::string minc_history=history;
	free(history);
#else
  std::string minc_history="";
#endif
  
  int verbose=0;
  std::string mask_f,list_f,majority_f,overlap_f;
	int classes=3;
  std::string map_f;
  int with_bg=0;
  int max_label=0;
	
  static struct option long_options[] = {
    {"verbose", no_argument,            &verbose, 1},
    {"quiet",   no_argument,            &verbose, 0},
    {"bg",      no_argument,            &with_bg, 1},
    {"mask",    required_argument,      0,'m'},
    {"classes", required_argument,      0,'c'},
    {"list",    required_argument,      0,'l'},
    {"majority", required_argument,     0,'a'},
    {"overlap", required_argument,      0,'o'},
    {"relabel", required_argument,      0,'r'},
    {0, 0, 0, 0}
  };
  
  for (;;) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "vqm:c:l:a:r:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;

    switch (c)
    {
      case 0:
        break;
      case 'v':
        std::cout << "Version: 0.1" << std::endl;
        return 0;
      case 'r':
        map_f=optarg;
        break;
      case 'm':
        mask_f=optarg;
        break;
      case 'l':
        list_f=optarg;
        break;
      case 'c':
        classes=atoi(optarg);
        break;
      case 'a':
        majority_f=optarg;
        break;
      case 'o':
        overlap_f=optarg;
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage (argv[0]);
        return 1;
    }
  }

  if(list_f.empty() && (argc - optind) < 2) {
    show_usage (argv[0]);
    return 1;
  }
  
  try
  {
    std::map<int,int> label_map;
    std::map<int,int> label_map2;
    
    if(!map_f.empty())
    {
      if(verbose)
        std::cout<<"Going to relabel input according to:"<<map_f.c_str()<<std::endl;
      std::ifstream in_map(map_f.c_str());
      
      while(!in_map.eof() && in_map.good())
      {
        int l1,l2;
        in_map>>l1>>l2;
        std::map<int,int>::iterator pos=label_map.find(l1);
        if(pos==label_map.end())
          label_map[l1]=l2;
        else if(label_map[l1]!=l2)
        {
          std::cerr<<"Warning: label "<<l1<<" mapped more then once!"<<std::endl;
          std::cerr<<"Previous map:"<<(*pos).second<<" New map:"<<l2<<std::endl;
          return 1;
        }
      }
    }
    
#if ( ITK_VERSION_MAJOR < 4 )
    itk::RegisterMincIO();
#endif
    
    
    VectorImageType::Pointer distribution=VectorImageType::New();
    ReaderType::Pointer reader=ReaderType::New();
    LabelReaderType::Pointer label_reader=LabelReaderType::New();
    
    LabelImageType::Pointer mask=0;
    LabelImageConstIterator* itm=0;
    
    if(!mask_f.empty())
    {
      
      label_reader->SetFileName(mask_f);
      label_reader->Update();
      mask=label_reader->GetOutput();
      mask->DisconnectPipeline();
      
      itm=new LabelImageConstIterator(mask,mask->GetLargestPossibleRegion());
      itm->GoToBegin();
    }

    std::set<IOPixelType> labels_set;
    
    int nfiles=argc-optind;
    
    std::vector<std::string> input_files;
    
    for(int i=0;i<nfiles;i++)
    {
      input_files.push_back(argv[i+optind]);
    }
    
    if(!list_f.empty())
    {
      std::ifstream inf(list_f.c_str());
      while(inf.good() && !inf.eof())
      {
        char tmp[1024];
        inf>>tmp;
        if(!strlen(tmp)) break;
        input_files.push_back(tmp);
      }
    }
    
    nfiles=input_files.size();
    
    if(!nfiles)
    {
      std::cerr<<"No input files provided!, aborting!"<<std::endl;
      return 1 ;
    }
    
    if(verbose)
      std::cout<<"Processing:"<<nfiles<< " files"<<std::endl;
  
    for(int i=0;i<nfiles;i++)
    {
      const char* fname=input_files[i].c_str();
      if(verbose)
        std::cout<<fname<<"\t"<<std::flush;
      
      reader->SetFileName(fname);
      reader->Update();
      IOImageType::Pointer img=reader->GetOutput();
      //store information from the image in the db
      
      if(!i)
      {
        distribution->SetLargestPossibleRegion(img->GetLargestPossibleRegion());
        distribution->SetBufferedRegion(img->GetLargestPossibleRegion());
        distribution->SetRequestedRegion(img->GetLargestPossibleRegion());
        distribution->SetSpacing( img->GetSpacing() );
        distribution->SetOrigin ( img->GetOrigin() );
        distribution->SetDirection(img->GetDirection());

        for(IOImageConstIterator it(img, img->GetBufferedRegion()); !it.IsAtEnd(); ++it)
        {
          IOPixelType val = it.Get();
          
          if(!label_map.empty())
          {
            std::map<int,int>::const_iterator f=label_map.find(val);
            if(f!=label_map.end())
              val=(*f).second;
            //TODO: if not found unchanged?
          }

          if(val>0 || with_bg) //we are only doing non-zero labels
            labels_set.insert(val);

          if(val>max_label) max_label=val;
        }
        if(verbose)
        {
          std::cout<<"Found:"<<labels_set.size()<<" labels"<<std::endl;
          for(std::set<IOPixelType>::iterator ls=labels_set.begin();ls!=labels_set.end() ;++ls)
            std::cout<<(*ls)<<"\t";
          std::cout<<std::endl;
          std::cout<<"Max label="<<max_label<<std::endl;
        }
        classes=labels_set.size();
        int j=0;
        
        for(std::set<IOPixelType>::iterator it=labels_set.begin();it!=labels_set.end();++it)
        {
          label_map2[*it]=j;
          j++;
        }
        
        distribution->SetNumberOfComponentsPerPixel(classes);
        itk::VariableLengthVector<float> zero(classes);
        zero.Fill(0);
        distribution->Allocate();
        distribution->FillBuffer(zero);
      } else {  
      	if(distribution->GetSpacing()!=img->GetSpacing())
      	{
      		std::cerr<<fname<<" spacing mismatch!"<<std::endl;return 1;
      	}
      	if(distribution->GetOrigin()!=img->GetOrigin())
      	{
      		std::cerr<<fname<<" origin mismatch!"<<std::endl;return 1;
      	}
      	if(distribution->GetDirection()!=img->GetDirection())
      	{
      		std::cerr<<fname<<" direction cosines mismatch!"<<std::endl;return 1;
      	}
      }
			
			IOImageConstIterator it1(img,img->GetLargestPossibleRegion());
			VectorImageIterator  it2(distribution,distribution->GetLargestPossibleRegion());
			
			if(mask) {
				itm->GoToBegin();
			}
			
			for(it1.GoToBegin(),it2.GoToBegin();!it1.IsAtEnd()&&!it2.IsAtEnd();++it1,++it2)
			{
				if(mask && !itm->Get()) {
					++(*itm);
					continue;
				}
				if(mask) ++(*itm);
				IOPixelType label=it1.Get();
        
				if(!label_map.empty())
        {
          std::map<int,int>::const_iterator f=label_map.find(label);
          if(f!=label_map.end())
            label=(*f).second;
          //TODO: if not found unchanged?
        }
        
        if(label>0 || with_bg )
        {
          label=label_map2[label];
          
          itk::VariableLengthVector<float> val=it2.Get();
          
          val[label]+=1;
          it2.Set(val);
        }
      }
    }
    
    IOImageType::Pointer majority=IOImageType::New();
    
    allocate_same<IOImageType,VectorImageType>(majority,distribution);
    majority->FillBuffer(0);
    
    RealImageType::Pointer overlap=RealImageType::New();
    allocate_same<RealImageType,VectorImageType>(overlap,distribution);
    overlap->FillBuffer(0.0);
    
    //now, let's see what we collected
    IOImageIterator it1(majority,majority->GetLargestPossibleRegion());
    VectorImageConstIterator it2(distribution,distribution->GetLargestPossibleRegion());
    
    RealImageIterator it3(overlap,overlap->GetLargestPossibleRegion());
    
    if(mask) 
      itm->GoToBegin();
    
    double a=0;
    double b=0;
    
    for(it1.GoToBegin(),it2.GoToBegin(),it3.GoToBegin();!it1.IsAtEnd()&&!it2.IsAtEnd();++it1,++it2,++it3)
    {
      if(mask && !itm->Get()) {
        ++(*itm);
        continue;
      }
      
      if(mask) ++(*itm);
      
      itk::VariableLengthVector<float> val=it2.Get();
      int mclass=-1;
      float mcount=0;
      double I=0.0;
      double U=0.0;
      for(int j=0;j<classes;j++)
      {
        if(val[j]>mcount)
        {
          mcount=val[j];
          mclass=j;
        }
        //discrete version of the formula from http://dx.doi.org/10.1109/TMI.2006.880587
        
        //group-wise union (nfiles-1)+ ...+(nfiles-val[j])
        U+=nfiles*(nfiles-1)/2-(nfiles-val[j])*(nfiles-val[j]-1)/2;
        //group-wise intersection: val[j]+val[j]-1+.....1 
        I+=val[j]*(val[j]-1)/2;
        
      }
      b+=U;//TODO add label weights
      a+=I;
      
      std::set<IOPixelType>::iterator ls;
      
      if(mclass>=0)
      {
        for(ls=labels_set.begin();ls!=labels_set.end() && mclass>0;++ls)
        {
          mclass--;
        }
        
        if(ls!=labels_set.end())
          it1.Set(*ls);
        else
          it1.Set(0); //ERROR?
      } else {
        it1.Set(0); //background!
      }
      if(U>0.0)
        it3.Set(I/U);
    }
    
    std::cout.precision(10);
    if(verbose)
      std::cout<<"Global overlap:";
    std::cout<<a/b<<std::endl;

#if ( ITK_VERSION_MAJOR < 4 )
    if( max_label<=255 )
      minc::set_minc_storage_type ( majority, NC_BYTE, false );
    else //if( max_label<=65535 )
      minc::set_minc_storage_type ( majority, NC_SHORT, false );
//     else
//       minc::set_minc_storage_type ( majority, NC_INT, false );

    minc::append_history ( majority, minc_history );

    minc::set_minc_storage_type ( overlap, NC_SHORT, false );
    minc::append_history ( overlap, minc_history );
#endif

    if(!majority_f.empty())
    {
      WriterType::Pointer writer=WriterType::New();
      writer->SetFileName(majority_f.c_str());
      writer->SetInput(majority);
      writer->Update();
    }
    
    if(!overlap_f.empty())
    {
      RealWriterType::Pointer writer=RealWriterType::New();
      writer->SetFileName(overlap_f.c_str());
      writer->SetInput(overlap);
      writer->Update();
    }
    
  } 
#if ( ITK_VERSION_MAJOR < 4 )
  catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
#endif  
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return 2;
  } 
  return 0;
}

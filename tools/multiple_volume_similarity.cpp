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
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageIOFactory.h>
#include <itkMincHelpers.h>
#include <getopt.h>
#include <map>

#include <time_stamp.h>    // for creating minc style history entry
#include <itkMincImageIOFactory.h>

typedef itk::Image<unsigned char, 3>  LabelImageType;
typedef itk::Image<float, 3>  RealImageType;
typedef itk::VectorImage<float, 3>  VectorImageType;

typedef itk::ImageFileReader<LabelImageType> ReaderType;
typedef itk::ImageFileWriter<LabelImageType> WriterType;
typedef itk::ImageFileWriter<RealImageType>  RealWriterType;

typedef itk::ImageRegionConstIterator<LabelImageType> LabelImageConstIterator;
typedef itk::ImageRegionConstIterator<VectorImageType> VectorImageConstIterator;
typedef itk::ImageRegionIterator<LabelImageType> LabelImageIterator;
typedef itk::ImageRegionIterator<VectorImageType> VectorImageIterator;
typedef itk::ImageRegionIterator<RealImageType> RealImageIterator;

using namespace minc;

void show_usage (const char * prog)
{
  std::cout<<"Program calculates multiple volume similarity metrics for discrete labels "<<std::endl
           <<"or Generalized Tanimoto coefficient (GTC)" <<std::endl
           <<"based on :  William R. Crum, Oscar Camara, and Derek L. G. Hill"
           <<"\"Generalized Overlap Measures for Evaluation and Validation in Medical Image Analysis \""
           <<" IEEE TRANSACTIONS ON MEDICAL IMAGING, VOL. 25, NO. 11, NOVEMBER 2006"<<std::endl
           <<"http://dx.doi.org/10.1109/TMI.2006.880587"<<std::endl<<std::endl
           <<"Usage: "<<prog<<" <input1.mnc> <input2.mnc> ... <inputN.mnc>  "<<std::endl
           <<"[--verbose --mask <mask.mnc> --classes <n> --list <file.list> --majority <output> --overlap <output> --relabel map.txt]"<<std::endl;
}

int main(int argc,char **argv)
{
	
	char *history = time_stamp(argc, argv); 
	std::string minc_history=history;
	free(history);

  int verbose=0;
  std::string mask_f,list_f,majority_f,overlap_f;
	int classes=3;
  std::string map_f;
	
  static struct option long_options[] = {
    {"verbose", no_argument,             &verbose, 1},
    {"quiet",   no_argument,             &verbose, 0},
    {"mask",    required_argument,       0,'m'},
		{"classes", required_argument,       0,'c'},
		{"list",    required_argument,       0,'l'},
		{"majority", required_argument,      0,'a'},
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
				std::cerr<<"Sorry, not implemented yet!"<<std::endl;
				return 1;
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

	//TODO: add reading from file list
  if((argc - optind) < 2) {
    show_usage (argv[0]);
    return 1;
  }
  
  try
  {
    std::map<int,int> label_map;
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
    
    itk::RegisterMincIO();
		int nfiles=argc-optind;
		
		if(verbose)
			std::cout<<"Processing:"<<nfiles<< " files"<<std::endl;
		
		VectorImageType::Pointer distribution=VectorImageType::New();
		ReaderType::Pointer reader=ReaderType::New();
		LabelImageType::Pointer mask=0;
		LabelImageConstIterator* itm=0;
		
		if(!mask_f.empty())
		{
			
			reader->SetFileName(mask_f);
			reader->Update();
			mask=reader->GetOutput();
			mask->DisconnectPipeline();
			
			itm=new LabelImageConstIterator(mask,mask->GetLargestPossibleRegion());
			itm->GoToBegin();
		}
		
		for(int i=0;i<nfiles;i++)
		{
			const char* fname=argv[i+optind];
			if(verbose)
				std::cout<<fname<<"\t"<<std::flush;
			
			reader->SetFileName(fname);
			reader->Update();
			LabelImageType::Pointer img=reader->GetOutput();
			//store information from the image in the db
			
			if(!i)
			{
				distribution->SetLargestPossibleRegion(img->GetLargestPossibleRegion());
				distribution->SetBufferedRegion(img->GetLargestPossibleRegion());
				distribution->SetRequestedRegion(img->GetLargestPossibleRegion());
				distribution->SetSpacing( img->GetSpacing() );
				distribution->SetOrigin ( img->GetOrigin() );
				distribution->SetDirection(img->GetDirection());
				distribution->SetNumberOfComponentsPerPixel(classes);
				distribution->Allocate();
				itk::VariableLengthVector<float> zero(classes);
				zero.Fill(0);
				distribution->FillBuffer(zero);
				
			}else {
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
			LabelImageConstIterator it1(img,img->GetLargestPossibleRegion());
			VectorImageIterator it2(distribution,distribution->GetLargestPossibleRegion());
			
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
				unsigned char label=it1.Get();
        
				
				if(!label_map.empty())
        {
          std::map<int,int>::const_iterator f=label_map.find(label);
          if(f!=label_map.end())
            label=(*f).second;
          //TODO: if not found unchanged?
        }
        if(!label) { //ignore zero label
          continue; 
        }
        label--;
        
        
				if(label>=classes) { //ignore label  with value that is too high
					continue; 
				}
				itk::VariableLengthVector<float> val=it2.Get();
        
				val[label]+=1;
				it2.Set(val);
			}
		}
		
		LabelImageType::Pointer majority=LabelImageType::New();
		allocate_same<LabelImageType,VectorImageType>(majority,distribution);
		majority->FillBuffer(0);
		
		RealImageType::Pointer overlap=RealImageType::New();
		allocate_same<RealImageType,VectorImageType>(overlap,distribution);
		overlap->FillBuffer(0.0);
		
		//now, let's see what we collected
		LabelImageIterator it1(majority,majority->GetLargestPossibleRegion());
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
			it1.Set(mclass+1);
			
			if(U>0.0)
				it3.Set(I/U);
			
		}
		
		std::cout.precision(10);
		if(verbose)
			std::cout<<"Global overlap:";
		std::cout<<a/b<<std::endl;
		
		minc::set_minc_storage_type ( majority, NC_BYTE, false );
		minc::append_history ( majority, minc_history );
		
		minc::set_minc_storage_type ( overlap, NC_SHORT, false );
		minc::append_history ( overlap, minc_history );
		
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
    
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
  return 0;
}

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include <iostream>
#include <getopt.h>
#include <minc_io_simple_volume.h>
#include <minc_1_simple.h>
#include <minc_1_simple_rw.h>
#include "minc_histograms.h"
#include <time_stamp.h>    // for creating minc style history entry
#include <math.h>
#include <pthread.h>
#include <vector>

using namespace minc;

class anlm_proc
{
  protected:
    
    struct thread_parameters{
      int ini;
      int fin;
      anlm_proc* _this; 
    };
    
  //thread running stuff
  static void* worker_thread(void *arg);
  void* thread_process(int ini,int fin);
  
  static double bessi0(double x);
  static double bessi1(double x);
  static double Epsi(double snr);
  
  static void Regularize(const simple_volume<double>& in,
                  simple_volume<double>& out,int r);
  
  static double distance(const simple_volume<double> & ima,
                          int x,int y,int z,
                          int nx,int ny,int nz,
                          int f);
  
  static double distance2(const simple_volume<double> & ima,
                   const simple_volume<double> & medias,
                   int x,int y,int z,
                   int nx,int ny,int nz,
                   int f);

  static void Average_block(const simple_volume<double>& ima,
                     int x,int y,int z,
                     int neighborhoodsize,
                     std::vector<double> &average, 
                     double weight,bool rician);
  
  static void Value_block(simple_volume<double> &Estimate, 
                   simple_volume<int> &Label,
                   int x,int y,int z,
                   int neighborhoodsize,
                   const std::vector<double>  &average, 
                   double global_sum);
  
  
  public:

  const simple_volume<double>& ima;
  simple_volume<double> fima;
  simple_volume<double> means;
  simple_volume<double> variances;
  
  simple_volume<int> Label;
  simple_volume<double> bias;
  simple_volume<double> Estimate;
  simple_volume<double> distances;
  
  int search_radius;
  int patch_radius;
  double imax;
  
  bool rician;
  double beta;
  bool debug;
  
  anlm_proc(const simple_volume<double>&img,int search,int patch,bool rician,double beta,bool _debug=false);
  
  void exec(int Nthreads);
};


void show_usage(const char *name)
{
  std::cerr 
      << "This program implements adaptative non-local denoising algorithm published in "<<std::endl
      << "Jose V. Manjon, Pierrick Coupe, Luis Marti-Bonmati, D. Louis Collins, Montserrat Robles \"Adaptive non-local means denoising of MR images with spatially varying noise levels\""
      << " Journal of Magnetic Resonance Imaging Volume 31, Issue 1, pages 192–203, January 2010"<<std::endl
      << " DOI: 10.1002/jmri.22003"<<std::endl
      << std::endl
      << "Usage: "<<name<<" <source> <output_prefix>" << std::endl
      << "\t--rician correct for rician noise (remove bias)"<<std::endl
      << "\t--search <n> search radius"<<std::endl
      << "\t--patch <n> patch radius"<<std::endl
      << "\t--beta <f> adjust weight of smoothing, default  1, <1.0 - less >1.0 - more"<<std::endl
      << "\t--mt <n> use N threads, default is 1, WARNING: currently mutlithreaded version produces different results"<<std::endl
      << "\t--verbose be verbose" << std::endl
      << "\t--double store output as double"<<std::endl
      << "\t--float store output as float"<<std::endl
      << "\t--short store output as short"<<std::endl
      << "\t--byte store output as byte"<<std::endl;
}

int main(int argc,char **argv)
{
  int verbose=0;
  int debug=0;
  int clobber=0;
  int patch_radius=1;
  int search_radius=2;
  int rician=0;
  int store_float=0;
  int store_short=0;
  int store_double=0;
  int store_byte=0;
  int threads=1;
  double beta=1.0;
  
  char *history = time_stamp(argc, argv); //maybe we should free it afterwards
  
  static struct option long_options[] =
  {
    {"verbose", no_argument , &verbose,      1},
    {"clobber", no_argument , &clobber,      1},
    {"debug",   no_argument , &debug,        1},
    {"quiet",   no_argument , &verbose,      0},
    {"float",   no_argument , &store_float,  1},
    {"short",   no_argument , &store_short,  1},
    {"byte",    no_argument , &store_byte,   1},
    {"double",  no_argument , &store_double, 1},
    {"rician",  no_argument , &rician,       1},
    {"mt",      required_argument,  0, 't'},
    {"threads", required_argument,  0, 't'},
    {"beta",    required_argument,  0, 'b'},
    {"search",  required_argument,  0, 's'},
    {"patch",    required_argument, 0, 'p'},
    {0, 0, 0, 0}
  };
  
  int c;
  for (;;)
  {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long (argc, argv, "t:b:s:p:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
      case 0:
        break;
      case 't':
        threads=atoi(optarg);
        if(threads<1) 
        {
          std::cerr<<"Warning! Number of threads should be >= 1!"<<std::endl;
          threads=1;
        }
        break;
      case 'b':
        beta=fabs(atof(optarg));
        break;
      case 's':
        search_radius=atoi(optarg);
        break;
      case 'p':
        patch_radius=atoi(optarg);
        break;
      case '?':
        /* getopt_long already printed an error message. */
      default:
        show_usage(argv[0]);
        return 1;
    }
  }

  if ((argc - optind) < 2)
  {
    show_usage(argv[0]);
    return 1;
  }
  std::string input_src_f=argv[optind];
  std::string output_prefix=argv[optind+1];
  
  //std::string output_f=output_prefix+"_denoised.mnc";
  std::string output_f=output_prefix;
  
  std::string output_distance=output_prefix+"_distance.mnc";
  std::string output_variances=output_prefix+"_variances.mnc";
  std::string output_means=output_prefix+"_means.mnc";
  std::string output_counts=output_prefix+"_counts.mnc";
  
  
  if (!clobber && !access(output_f.c_str(), F_OK))
  {
    std::cerr << output_f.c_str () << " Exists!" << std::endl;
    return 1;
  }  
 
  try
  {
    simple_volume<double> src;
    
    minc_1_reader rdr1;
    rdr1.open(input_src_f.c_str());
    load_simple_volume<double>(rdr1,src);
    
    nc_type store_datatype= store_double?NC_DOUBLE:store_float?NC_FLOAT:store_short?NC_SHORT:store_byte?NC_BYTE:rdr1.datatype();
		
    if(debug)
      std::cout<<"Image dimensions:"<<src.dim(0)<<","<<src.dim(1)<<","<<src.dim(2)<<std::endl;
		
    anlm_proc anlm(src,search_radius,patch_radius,rician,beta,debug>0);
    
    anlm.exec(threads);
    std::cout<<"Done..."<<std::endl;
    
    minc_1_writer wrt;
    wrt.open(output_f.c_str(),rdr1.info(),2,store_datatype);
    wrt.append_history(history);
    save_simple_volume<double>(wrt,anlm.fima);
    
    if(debug)
    {
      std::cerr<<"Outputting debug information..."<<std::endl;
      
      minc_1_writer wrt2;
      wrt2.open(output_distance.c_str(),rdr1.info(),2,store_datatype);
      wrt2.append_history(history);
      save_simple_volume<double>(wrt2,anlm.distances);
      
      
      minc_1_writer wrt3;
      wrt3.open(output_counts.c_str(),rdr1.info(),2,store_datatype);
      wrt3.append_history(history);
      save_simple_volume<int>(wrt3,anlm.Label);
    }
    
    return 0;
  } catch (const minc::generic_error & err) {
    std::cerr << "Got an error at:" << err.file () << ":" << err.line () << std::endl;
    std::cerr << err.msg()<<std::endl;
    return 1;
  }
}


/*Returns the modified Bessel function I0(x) for any real x.*/
double anlm_proc::bessi0(double x)
{
  double ax,ans,a;
  double y; 
  if ((ax=fabs(x)) < 3.75) 
  { 
    y=x/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } 
  else 
  {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax));
    a=y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))));
    ans=ans*(0.39894228 + y*(0.1328592e-1 +y*(0.225319e-2+y*(-0.157565e-2+a))));    
  }
  return ans;
}

/*Returns the modified Bessel function I1(x) for any real x.*/
double anlm_proc::bessi1(double x)
{
  double ax,ans;
  double y; 
  if ((ax=fabs(x)) < 3.75)
  { 
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  } 
  else 
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2+y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans *= (exp(ax)/sqrt(ax));
  }
  return x < 0.0 ? -ans : ans;
}

double anlm_proc::Epsi(double snr)
{
  double val;    
  val=2 + snr*snr - (M_PI/8)*exp(-(snr*snr)/2)*((2+snr*snr)*bessi0((snr*snr)/4) + (snr*snr)*bessi1((snr*snr)/4))*((2+snr*snr)*bessi0((snr*snr)/4) + (snr*snr)*bessi1((snr*snr)/4));     
  if(val<0.001) val=1;
  if(val>10) val=1;    
  return val;
}

/* Function which compute the weighted average for one block */
void  anlm_proc::Average_block(const simple_volume<double>& ima,
                   int x,int y,int z,
                   int neighborhoodsize,
                   std::vector<double> &average, 
                   double weight,
                   bool rician)
{
  int x_pos,y_pos,z_pos;
  bool is_outside; 
  int a,b,c,count;

  count = 0;

  for (c = -neighborhoodsize; c<=neighborhoodsize;c++)
  {
    for (b = -neighborhoodsize; b<=neighborhoodsize;b++)
    {
      for (a = -neighborhoodsize; a<=neighborhoodsize;a++)
      {	
        x_pos = x+a;
        y_pos = y+b;
        z_pos = z+c;
	
        is_outside=!ima.hit(x_pos,y_pos,z_pos);
		
        if(rician)
        {
          if (is_outside)
            average[count] += ima(x,y,z)*ima(x,y,z)*weight;
          else	
            average[count] += ima(x_pos,y_pos,z_pos)*ima(x_pos,y_pos,z_pos)*weight;
        }
        else
        {
          if (is_outside)
            average[count] += ima(x,y,z)*weight;
          else	
            average[count] += ima(x_pos,y_pos,z_pos)*weight;
        }
        count++;
      }
    }
  }
}

/* Function which computes the value assigned to each voxel */
void  anlm_proc::Value_block(simple_volume<double> &Estimate, 
                 simple_volume<int> &Label,
                 int x,int y,int z,
                 int neighborhoodsize,
                 const std::vector<double>  &average, 
                 double global_sum)
{
  int x_pos,y_pos,z_pos;
  double value = 0.0;
  double label = 0.0;
  int count=0 ;
  int a,b,c;


  for (c = -neighborhoodsize; c<=neighborhoodsize;c++)
  {
    for (b = -neighborhoodsize; b<=neighborhoodsize;b++)
    {
      for (a = -neighborhoodsize; a<=neighborhoodsize;a++)
      {	
        x_pos = x+a;
        y_pos = y+b;
        z_pos = z+c;
	
        if (Estimate.hit(x_pos,y_pos,z_pos))
        {
          Estimate(x_pos,y_pos,z_pos) += average[count]/global_sum;
          Label(x_pos,y_pos,z_pos)    += 1;
        }
        count++;
      }
    }
  }
}

double  anlm_proc::distance(const simple_volume<double> & ima,int x,int y,int z,int nx,int ny,int nz,int f)
{
  double d,acu,distancetotal;
  int i,j,k,ni1,nj1,ni2,nj2,nk1,nk2;

  distancetotal=0.0;
  
  for(k=-f;k<=f;k++)
  {
    nk1=z+k;
    nk2=nz+k;  

    for(j=-f;j<=f;j++)
    {
      nj1=y+j;
      nj2=ny+j;    

      for(i=-f;i<=f;i++)
      {
        ni1=x+i;
        ni2=nx+i;
        double d=ima.safe_get(ni1,nj1,nk1)-ima.safe_get(ni2,nj2,nk2);
        distancetotal+= d*d; //L2 norm
      }
    }
  }

  acu=(2*f+1)*(2*f+1)*(2*f+1);
  d=distancetotal/acu;

  return d;
}

double  anlm_proc::distance2(const simple_volume<double> & ima,
                 const simple_volume<double> & medias,
                 int x,int y,int z,
                 int nx,int ny,int nz,
                 int f)
{
  double d,acu,distancetotal;
  int i,j,k,ni1,nj1,ni2,nj2,nk1,nk2;

  acu=0;
  distancetotal=0.0;
	
  for(k=-f;k<=f;k++)
  {
    nk1=z+k;
    nk2=nz+k;  
    
    for(j=-f;j<=f;j++)
    {
      nj1=y+j;
      nj2=ny+j;    
      for(i=-f;i<=f;i++)
      {
        ni1=x+i;
        ni2=nx+i;    
			
        double d1=ima.safe_get(ni1,nj1,nk1)-medias.safe_get(ni1,nj1,nk1);
        double d2=ima.safe_get(ni2,nj2,nk2)-medias.safe_get(ni2,nj2,nk2);
        
        distancetotal += (d1-d2)*(d1-d2);//L2 norm on high band 
      }
    }
  }

  acu=(2*f+1)*(2*f+1)*(2*f+1);
  d=distancetotal/acu;

  return d;
}

void  anlm_proc::Regularize(const simple_volume<double>& in,
                simple_volume<double>& out,int r)
{
  double acu;
  int ind,i,j,k,ni,nj,nk,ii,jj,kk;

  //double * temp=new double[sx*sy*sz];
  simple_volume<double> temp(in.size());

  for(k=0;k<in.dim(2);k++)
    for(j=0;j<in.dim(1);j++)
      for(i=0;i<in.dim(0);i++)
  {
    if(in(i,j,k)==0.0) continue; //VF: check this
  
    acu=0;
    ind=0;
    for(ii=-r;ii<=r;ii++)
    {
      ni=i+ii;
      double v=in.safe_get(ni,j,k);
      if(v>0)
      {
        acu+=v;
        ind++;		
      }
    }
    if(ind==0) ind=1; 
    out(i,j,k)=acu/ind;
  }
  
  for(k=0;k<in.dim(2);k++)
    for(j=0;j<in.dim(1);j++)
      for(i=0;i<in.dim(0);i++)
  {
    if(out(i,j,k)==0.0) continue;
  
    acu=0;
    ind=0;
    for(jj=-r;jj<=r;jj++)
    {	
      nj=j+jj;		
      double v=out.safe_get(i,nj,k);
      if(v>0)
      {
        acu+=v;
        ind++;		
      }
    }
    if(ind==0) ind=1;
    temp(i,j,k)=acu/ind;
  }
  
  for(k=0;k<in.dim(2);k++)
    for(j=0;j<in.dim(1);j++)
      for(i=0;i<in.dim(0);i++)
  {
    if(temp(i,j,k)==0.0) continue;
  
    acu=0;
    ind=0;
    for(kk=-r;kk<=r;kk++)
    {
      nk=k+kk;			
      double v=temp.safe_get(i,j,nk);
      if(v>0)
      {
        acu+=v;
        ind++;		
      }
    }
    if(ind==0) ind=1;
    out(i,j,k)=acu/ind;
  }
}

void * anlm_proc::thread_process(int ini,int fin)
{
  int i,j,k,ii,jj,kk,ni,nj,nk;
  if(debug)
  {
    std::cerr<<"Processing:"<<ini<<" - "<<fin<<std::endl;
    std::cerr<<"Search radius:"<<search_radius<<" Patch radius: "<<patch_radius<<std::endl;
  }
  
  const double epsilon = 1e-5;
  const double mu1  = 0.95;
  const double var1 = 0.5;
  const double undef_distance=1e10;
  int   patch_size=(2*patch_radius+1)*(2*patch_radius+1)*(2*patch_radius+1);

  std::vector<double> average(patch_size);


  for(k=ini;k<fin;k+=2)
    for(j=0;j<ima.dim(1);j+=2)
      for(i=0;i<ima.dim(0);i+=2)
  { 

    double wmax=0.0;
    double totalweight=0.0;
    double distanciaminima=undef_distance;
    
    for (int init=0 ; init < patch_size; init++) 
      average[init]=0.0;	 

       
    if(ima(i,j,k)>0.0       && //VF: why is that?
       means(i,j,k)>epsilon && 
       variances(i,j,k)>epsilon)
    {
      for(kk=-search_radius; kk<=search_radius; kk++)
      {
        nk=k+kk;
        for(jj=-search_radius; jj<=search_radius; jj++)
        {
          nj=j+jj;
          for(ii=-search_radius; ii<=search_radius; ii++)
          {
            if(!ii && !jj && !kk) continue;
            //if(ii==0 || jj==0 || kk==0) continue;
            ni=i+ii;
            if(ima.hit(ni,nj,nk))
            {									
              if (ima      (ni,nj,nk)>0.0 &&  //VF: why is that?
                  means    (ni,nj,nk)>epsilon && 
                  variances(ni,nj,nk)>epsilon)
              {				
                double t1 =  means(i,j,k)/means(ni,nj,nk);  
                double t1i= (imax-means(i,j,k))/(imax-means(ni,nj,nk));
                double t2 = (variances(i,j,k))/(variances(ni,nj,nk));
	
                if( (t1>mu1  && t1<(1.0/mu1) || 
                     t1i>mu1 && t1i<(1.0/mu1)) && 
                     t2>var1 && 
                     t2<(1.0/var1)) //VF: check this
                {
                  double d=distance2(ima,means,i,j,k, ni,nj,nk, patch_radius);

                  if(d<distanciaminima) distanciaminima=d;
                }
              }
            }
          }
        }
      }
      
      
      distances(i,j,k)=(distanciaminima==undef_distance)?0.0:distanciaminima;
      
      if( distanciaminima < epsilon)
        distanciaminima=1.0;
      
          
      if(rician)
      {
        for(kk=-patch_radius;kk<=patch_radius;kk++)
        {
          nk=k+kk;
          for(ii=-patch_radius;ii<=patch_radius;ii++)
          {
            ni=i+ii;
            for(jj=-patch_radius;jj<=patch_radius;jj++)
            {
              nj=j+jj;							
              if(bias.hit(ni,nj,nk))
              {   
                if(distanciaminima==1e10) 
                  bias(ni,nj,nk)=0.0; 
                else 
                  bias(ni,nj,nk)=distanciaminima;
              }
            }
          }
        }
      }
      
      //block design 
      for(kk=-search_radius;kk<=search_radius;kk++)
      {
        nk=k+kk;
        for(jj=-search_radius;jj<=search_radius;jj++)
        {
          nj=j+jj;
          for(ii=-search_radius;ii<=search_radius;ii++)
          {
            //if(!ii || !jj || !kk) continue;  
            ni=i+ii;														
            if(!ii && !jj && !kk) continue;
				
            if(ima.hit(ni,nj,nk))
            {									
              if (ima(ni,nj,nk)>0.0 &&  //Why 0.0
                  means(ni,nj,nk)> epsilon && 
                  variances(ni,nj,nk)>epsilon)
              {				
                double t1 = means(i,j,k)/means(ni,nj,nk);  
                double t1i= (imax-means(i,j,k))/(imax-means(ni,nj,nk));  
                double t2 = variances(i,j,k)/variances(ni,nj,nk);
	
                if( (t1>mu1  && t1<(1.0/mu1) || t1i>mu1 && t1i<(1.0/mu1)) && 
                     t2>var1 && t2<(1.0/var1) ) //VF: check this
                {
                  double d=distance(ima,i,j,k,ni,nj,nk,patch_radius);
                  double w=0.0;
                  
                  if(d <= 0.0 || std::isnan(d) || std::isinf(d)) d=epsilon;
                  
                  if( d > 3.0*distanciaminima ) 
                    w=0.0;
                  else
                  {
                    w = exp(-d/(distanciaminima*beta)); //
                    
                    if(std::isnan(w) || std::isinf(w)) 
                      w=0.0;
                  }
                  
                  if(w>wmax) 
                    wmax = w;
										
                  if(w>0.0)
                  {
                    Average_block(ima,ni,nj,nk,patch_radius,average,w,rician);
                    totalweight += w;
                  }
                }
              }
            }
          }
        }
      }
      
      if(wmax==0.0) 
        wmax=1.0;
      
      Average_block(ima,i,j,k,patch_radius,average,wmax,rician);
      totalweight += wmax;
      
      Value_block(Estimate,Label,i,j,k,patch_radius,average,totalweight);
    }
    else  //local means and variances are too small?
    {
      wmax=1.0;
      Average_block(ima,i,j,k,patch_radius,average,wmax,rician);	
      totalweight = totalweight + wmax;
      Value_block(Estimate,Label,i,j,k,patch_radius,average,totalweight);
    }
  }

  //std::cout<<"Finishing thread "<<ini<<" - "<<fin<<std::endl;
  
	return 0;
}


anlm_proc::anlm_proc(
          const simple_volume<double>& _ima,
          int search,int patch,
          bool _rician,
          double _beta,
          bool _debug):
    ima(_ima),
    fima(ima.size()),
    means(ima.size()),
    variances(ima.size()),
    Label(ima.size()),
    bias(ima.size()),
    Estimate(ima.size()),
    distances(ima.size()),
    search_radius(search),
    patch_radius(patch),
    imax(-1e10),
    rician(_rician),beta(_beta),debug(_debug)
{
}

void anlm_proc::exec(int Nthreads)
{
  double SNR,h;
  int i,j,k,ii,jj,kk,ni,nj,nk,ini,fin;
  
  imax=-1e10;
  
  //Ndims = pow((2*patch_radius+1),3);
  /*Allocate memory and assign output pointer*/

  for (i = 0; i < ima.c_buf_size();i++)
  {
    Estimate.c_buf()[i] = 0.0;
    Label.c_buf()[i] = 0;
    fima.c_buf()[i] = 0.0;
    if(rician) bias.c_buf()[i]=0.0;
    distances.c_buf()[i]=0.0;
    means.c_buf()[i]=0.0;
  }

  std::cout<<"Calculating means..."<<std::endl;
  
  for(k=0;k<ima.dim(2);k++)
  {
    for(j=0;j<ima.dim(1);j++)
    {
      for(i=0;i<ima.dim(0);i++)
      {
        if(ima(i,j,k)>imax) 
          imax=ima(i,j,k);

        double m=0;
        int c=0;
        
        for(ii=-1;ii<=1;ii++)
        {
          for(jj=-1;jj<=1;jj++)
          {
            for(kk=-1;kk<=1;kk++)
            {
              m += ima.safe_get(i+ii,j+jj,k+kk);
              c++;
            }
          }
        }
        means.set(i,j,k,m/c);
      }
    }
  }
  std::cout<<"Image maximum:"<<imax<<std::endl;
  std::cout<<"Calculating variances..."<<std::endl;
  for(k=0;k<ima.dim(2);k++)
  {
    for(j=0;j<ima.dim(1);j++)
    {
      for(i=0;i<ima.dim(0);i++)
      {
        double var=0;
        int c=0;
        for(ii=-1;ii<=1;ii++)
        {
          for(jj=-1;jj<=1;jj++)
          {
            for(kk=-1;kk<=1;kk++)
            {
              ni=i+ii;
              nj=j+jj;
              nk=k+kk;
              if(ima.hit(ni,nj,nk))
              {
                double d=ima(ni,nj,nk)-means(i,j,k);
                var += d*d;
                c++;
              }
            }
          }
        }
        variances.set(i,j,k,var/(c-1));
      }
    }
  }


  std::cout<<"Launching: "<<Nthreads<<" thread(s)..."<<std::endl;
  
  std::vector<anlm_proc::thread_parameters> ThreadArgs(Nthreads);
  std::vector<pthread_t> ThreadList(Nthreads);
  
  for (i=0; i<Nthreads; i++)
  {         
    ThreadArgs[i].ini=(i*ima.dim(2))/Nthreads;
    ThreadArgs[i].fin=(i+1)*ima.dim(2)/Nthreads;
    ThreadArgs[i]._this=this;
  }
  ThreadArgs[Nthreads-1].fin=ima.dim(2);
    
  for (i=1; i<Nthreads; i++)
  { 
    if(pthread_create(&ThreadList[i], NULL, worker_thread,&ThreadArgs[i]))
    {
      std::cerr<<"Threads cannot be created"<<std::endl;
      exit(1);
    }
  }
  //VF: launch one thread here
  thread_process(ThreadArgs[0].ini,ThreadArgs[0].fin);
    
  for (i=1; i<Nthreads; i++)
  {
    pthread_join(ThreadList[i],NULL);
  }

  if(rician)
  {
    int r=5;  
    Regularize(bias,variances,r);
    for(i=0;i<ima.c_buf_size();i++)
    {
      if(variances.c_buf()[i]>0) 
      {
        SNR=means.c_buf()[i]/sqrt(variances.c_buf()[i]);
        bias.c_buf()[i]=2*(variances.c_buf()[i]/Epsi(SNR));
        if(std::isnan(bias.c_buf()[i])) bias.c_buf()[i]=0;
      }
    }
  }
  /* Aggregation of the estimators (i.e. means computation) */
  for(i=0;i<ima.c_buf_size();i++)
  {
    int label = Label.c_buf()[i];
    double estimate=0.0;
    if (label < 1 ) 
      fima.c_buf()[i] = ima.c_buf()[i];
    else
    {
      estimate = Estimate.c_buf()[i]/label;
      if(rician)
      {
        double d=estimate-bias.c_buf()[i];
        estimate = sqrt(d<0.0 ? 0.0:d);
      }

      if(!std::isnan(estimate) && !std::isinf(estimate))
        fima.c_buf()[i]=estimate;	
      else
        fima[i]=0.0;
    }
  }
}

void* anlm_proc::worker_thread(void *arg)
{
  anlm_proc::thread_parameters *par=(anlm_proc::thread_parameters *)arg;
  return par->_this->thread_process(par->ini,par->fin);
}

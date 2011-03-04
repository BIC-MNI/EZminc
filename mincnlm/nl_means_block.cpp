#include "nl_means_utils.h"
#include "nl_means.h"
#include "nl_means_block.h"

extern int      verbose;
extern int      debug;

// Function which compute the weighted average for one block
void Average_block ( float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float weight, int* vol_size )
{
  int x_pos,y_pos,z_pos;
  bool is_outside;


  int count = 0;

  for ( int c = 0; c< ( 2*neighborhoodsize[2]+1 );c++ )
  {
    for ( int b = 0; b< ( 2*neighborhoodsize[1]+1 );b++ )
    {
      for ( int a = 0;a< ( 2*neighborhoodsize[0]+1 );a++ )
      {

        is_outside = false;
        x_pos = x+a-neighborhoodsize[0];
        y_pos = y+b-neighborhoodsize[1];
        z_pos = z+c-neighborhoodsize[2];

        if ( ( z_pos < 0 ) || ( z_pos > vol_size[2]-1 ) ) is_outside = true;
        if ( ( y_pos < 0 ) || ( y_pos > vol_size[1]-1 ) ) is_outside = true;
        if ( ( x_pos < 0 ) || ( x_pos > vol_size[0]-1 ) ) is_outside = true;

        if ( is_outside )
          average[count] = average[count] + ima_in[z* ( vol_size[0]*vol_size[1] ) + ( y*vol_size[0] ) +x]*weight;
        else
          average[count] = average[count] + ima_in[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos]*weight;

        count++;
      }
    }
  }
}

// Function which compute the weighted average for one block
void Average_block_Rician ( float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float weight, int* vol_size )
{
  int x_pos,y_pos,z_pos;
  bool is_outside;


  int count = 0;

  for ( int c = 0; c< ( 2*neighborhoodsize[2]+1 );c++ )
  {
    for ( int b = 0; b< ( 2*neighborhoodsize[1]+1 );b++ )
    {
      for ( int a = 0;a< ( 2*neighborhoodsize[0]+1 );a++ )
      {

        is_outside = false;
        x_pos = x+a-neighborhoodsize[0];
        y_pos = y+b-neighborhoodsize[1];
        z_pos = z+c-neighborhoodsize[2];

        if ( ( z_pos < 0 ) || ( z_pos > vol_size[2]-1 ) ) is_outside = true;
        if ( ( y_pos < 0 ) || ( y_pos > vol_size[1]-1 ) ) is_outside = true;
        if ( ( x_pos < 0 ) || ( x_pos > vol_size[0]-1 ) ) is_outside = true;
        if ( is_outside )
          average[count] = average[count] + ima_in[z* ( vol_size[0]*vol_size[1] ) + ( y*vol_size[0] ) +x]*ima_in[z* ( vol_size[0]*vol_size[1] ) + ( y*vol_size[0] ) +x]*weight;
        else
          average[count] = average[count] + ima_in[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos]*ima_in[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos]*weight;

        count++;
      }
    }
  }
}


// Function which computes the value assigned to each voxel
void Value_block ( pthread_mutex_t *mp, volatile float *Estimate, volatile float* Label, float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float global_sum, int* vol_size )
{
  int x_pos,y_pos,z_pos;
  int ret;
  bool is_outside;
  float value = 0.0;
  float label = 0;
  int count=0 ;

  for ( int c = 0; c< ( 2*neighborhoodsize[2]+1 );c++ )
  {
    for ( int b = 0; b< ( 2*neighborhoodsize[1]+1 );b++ )
    {
      for ( int a = 0;a< ( 2*neighborhoodsize[0]+1 );a++ )
      {


        is_outside = false;
        x_pos = x+a-neighborhoodsize[0];
        y_pos = y+b-neighborhoodsize[1];
        z_pos = z+c-neighborhoodsize[2];

        if ( ( z_pos < 0 ) || ( z_pos > vol_size[2]-1 ) ) is_outside = true;
        if ( ( y_pos < 0 ) || ( y_pos > vol_size[1]-1 ) ) is_outside = true;
        if ( ( x_pos < 0 ) || ( x_pos > vol_size[0]-1 ) ) is_outside = true;
        if ( !is_outside )
        {

          value = Estimate[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos];
          value = value + ( average[count]/global_sum );
          label = Label[ ( x_pos + y_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1] ) ];

          // Lock to avoid multiple access by thread
          ret = pthread_mutex_lock ( mp );
          Estimate[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos] = value;
          Label[ ( x_pos + y_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1] ) ] = label +1;
          ret = pthread_mutex_unlock ( mp );

        }
        count++;
      }
    }
  }
}

// Function which computes the value assigned to each voxel
void Value_block_Rician ( pthread_mutex_t *mp, volatile float *Estimate, volatile float* Label, float *ima_in,int x,int y,int z,int *neighborhoodsize,float *average, float global_sum, int* vol_size, float filtering_param )
{
  int x_pos,y_pos,z_pos;
  int ret;
  bool is_outside;
  float value = 0.0;
  float label = 0;
  int count=0 ;
  float denoised_value= 0.0;

  for ( int c = 0; c< ( 2*neighborhoodsize[2]+1 );c++ )
  {
    for ( int b = 0; b< ( 2*neighborhoodsize[1]+1 );b++ )
    {
      for ( int a = 0;a< ( 2*neighborhoodsize[0]+1 );a++ )
      {


        is_outside = false;
        x_pos = x+a-neighborhoodsize[0];
        y_pos = y+b-neighborhoodsize[1];
        z_pos = z+c-neighborhoodsize[2];

        if ( ( z_pos < 0 ) || ( z_pos > vol_size[2]-1 ) ) is_outside = true;
        if ( ( y_pos < 0 ) || ( y_pos > vol_size[1]-1 ) ) is_outside = true;
        if ( ( x_pos < 0 ) || ( x_pos > vol_size[0]-1 ) ) is_outside = true;
        if ( !is_outside )
        {
          value = Estimate[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos];

          denoised_value  = ( average[count]/global_sum ) - ( 2*filtering_param*filtering_param );

          if ( denoised_value > 0 )
            denoised_value = sqrt ( denoised_value );
          else denoised_value = 0.0;

          value = value + denoised_value;

          label = Label[ ( x_pos + y_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1] ) ];

          // Lock to avoid multiple access by thread
          ret = pthread_mutex_lock ( mp );
          Estimate[z_pos* ( vol_size[0]*vol_size[1] ) + ( y_pos*vol_size[0] ) +x_pos] = value;
          Label[ ( x_pos + y_pos*vol_size[0] + z_pos *vol_size[0] * vol_size[1] ) ] = label +1;
          ret = pthread_mutex_unlock ( mp );

        }
        count++;
      }
    }
  }
}

//Main function

void *Sub_denoise_block_mt ( void *arguments )
{
  nl_mean_block_mt arg;
  arg=* ( nl_mean_block_mt * ) arguments;

  int PID=0;
  PID=arg.thread_num;
  int ret;

  int x_min,x_max,y_min,y_max,z_min,z_max;
  double epsilon=0.0001;

  int NbElement2 = ( 2*arg.neighborhoodsize[0]+1 ) * ( 2*arg.neighborhoodsize[1]+1 ) * ( 2*arg.neighborhoodsize[2]+1 );

  typedef float Voisinages[NbElement2];

  Voisinages V1,V2;

  float w_max=0;

  int offset_k = 0;
  int offset_j = 0;
  int offset_kk = 0;
  int offset_jj = 0;

  for ( int k=arg.debut; k < arg.fin ; k+=arg.b_space )
  {
    double count=0;
    offset_k = k* ( arg.vol_size[0]*arg.vol_size[1] );

    for ( int j=0; j < arg.vol_size[1]; j+=arg.b_space )
    {

      offset_j = j*arg.vol_size[0];

      for ( int i=0; i < arg.vol_size[0]; i+=arg.b_space )
      {

        float global_sum=0.;
        float average[NbElement2];
        for ( int init=0 ; init < NbElement2; init++ )
        {
          average[init]=0.;
        }
        int addr=offset_k + offset_j + i;


        x_min = MAX ( 0,i-arg.search[0] );      x_max = MIN ( arg.vol_size[0]-1,i+arg.search[0] );
        y_min = MAX ( 0,j-arg.search[1] );      y_max = MIN ( arg.vol_size[1]-1,j+arg.search[1] );
        z_min = MAX ( 0,k-arg.search[2] );      z_max = MIN ( arg.vol_size[2]-1,k+arg.search[2] );

        // Only the voxel with a local mean and variance are processed.
        // By this way the computational time is reduced since the background is not taken into account (or the B-scan ;ask for US image)

        if ( ( arg.mean_map[addr] > epsilon ) && ( arg.var_map[addr] > epsilon ) )
        {

          if ( ( arg.weight_method == 0 ) || ( arg.weight_method == 2 ) ) Neiborghood ( arg.in,i,j,k,arg.neighborhoodsize,V1,arg.vol_size,arg.weight_method );
          if ( arg.weight_method == 1 ) Neiborghood ( arg.in,i,j,k,arg.neighborhoodsize,V1,arg.vol_size,arg.weight_method );

          w_max = 0;

          for ( int kk = z_min; kk <= z_max; kk++ )
          {

            offset_kk = kk* ( arg.vol_size[0]*arg.vol_size[1] );

            for ( int jj = y_min; jj <= y_max; jj++ )
            {

              offset_jj = jj*arg.vol_size[0];

              for ( int ii = x_min; ii <= x_max; ii++ )
              {

                int addr2=offset_kk + offset_jj + ii;
                // To avoid the compute the weigh with null patch and the division by zero in ratio and ratio 2

                if ( ( arg.mean_map[addr2] > epsilon ) && ( arg.var_map[addr2] > epsilon ) )
                {

                  float ratio  = arg.mean_map[addr]/arg.mean_map[addr2];
                  float ratio2 = arg.var_map[addr]/arg.var_map[addr2];

                  float weight = 0.;

                  if ( ( arg.weight_method == 0 ) || ( arg.weight_method == 2 ) )
                  {
                    if ( ( arg.m_min <= ratio ) && ( ratio <= ( 1/arg.m_min ) ) && ( ( arg.v_min <= ratio2 ) && ( ratio2 <= ( 1/arg.v_min ) ) ) )
                    {
                      if ( ( ii != i ) || ( jj != j ) || ( kk != k ) )
                      {
                        count = count + 1;
                        Neiborghood ( arg.in,ii,jj,kk,arg.neighborhoodsize,V2,arg.vol_size,arg.weight_method );
                        weight = Weight ( L2_norm ( V1,V2,arg.neighborhoodsize ), arg.beta/2.0, arg.filtering_param );
                        global_sum = global_sum + weight;

                        if ( arg.weight_method == 0 )
                          Average_block ( arg.hallucinate?arg.hallucinate:arg.in,
                            ii,jj,kk,arg.neighborhoodsize,average,weight,arg.vol_size );

                        else
                          Average_block_Rician ( arg.hallucinate?arg.hallucinate:arg.in,
                            ii,jj,kk,arg.neighborhoodsize,average,weight,arg.vol_size );

                        if ( weight > w_max ) w_max = weight;
                      }
                    }
                  } else if ( arg.weight_method == 1 ) {

                    if ( ( arg.m_min <= ratio ) && ( ratio <= ( 1/arg.m_min ) ) )
                    {
                      if ( ( ii != i ) || ( jj != j ) || ( kk != k ) )
                      {
                        count = count + 1;
                        //no normalization of thr neighborhood
                        Neiborghood ( arg.in,ii,jj,kk,arg.neighborhoodsize,V2,arg.vol_size,arg.weight_method );
                        // Pearson distance computation
                        weight = Weight ( Pearson_distance ( V1,V2,arg.neighborhoodsize ),arg.beta/2.0, arg.filtering_param );
                        global_sum = global_sum + weight;
                        Average_block ( arg.hallucinate?arg.hallucinate:arg.in,
                           ii,jj,kk,arg.neighborhoodsize,average,weight,arg.vol_size );
                        if ( weight > w_max ) w_max = weight;
                      }
                    }
                  }
                }
              }
            }
          }


          //  w_max =1;
          global_sum = global_sum + w_max;

          if ( arg.weight_method == 2 )
            Average_block_Rician ( arg.hallucinate?arg.hallucinate:arg.in,
              i,j,k,arg.neighborhoodsize,average,w_max,arg.vol_size );
          else
            Average_block ( arg.hallucinate?arg.hallucinate:arg.in,
              i,j,k,arg.neighborhoodsize,average,w_max,arg.vol_size );

          if ( global_sum != 0 )
          {
            if ( arg.weight_method == 2 )
              Value_block_Rician ( arg.mp,arg.Estimate,arg.Label,arg.hallucinate?arg.hallucinate:arg.in,
                i,j,k,arg.neighborhoodsize,average,global_sum,arg.vol_size,arg.filtering_param );
            else
              Value_block ( arg.mp,arg.Estimate,arg.Label,arg.hallucinate?arg.hallucinate:arg.in,
                i,j,k,arg.neighborhoodsize,average,global_sum,arg.vol_size );
          }
        }

      }
    }


    if ( debug )
      std::cout<<"Thread ( "<< setw ( 2 ) <<PID+1<< " ) has finished the slice :"<<setw ( 3 ) <<k<<std::endl;

    ret = pthread_mutex_lock ( arg.mp );
    *arg.mean_neighboors = *arg.mean_neighboors + count/ ( ( arg.vol_size[0]*arg.vol_size[1] ) / ( arg.b_space*arg.b_space ) );
    ret = pthread_mutex_unlock ( arg.mp );

  }

  if ( debug )
    cout<<"End of the thread : "<<setw ( 2 ) <<PID+1<<endl;

  pthread_exit ( 0 );


}


//Creation of the threads
void denoise_block_mt ( float *in, float *out, float *mean_map, float *var_map, double filtering_param,double beta, int *neighborhoodsize, int *search,int testmean,int testvar,double m_min, double v_min,int weight_method, int b_space, int *vol_size,float *hallucinate )
{

  int ii;
  pthread_t *thread;
  nl_mean_block_mt *arguments;
  void *retval;
  thread = ( pthread_t * ) calloc ( nb_thread, sizeof ( pthread_t ) );
  arguments= ( nl_mean_block_mt* ) calloc ( nb_thread,sizeof ( nl_mean_block_mt ) );

  float *Estimate, *Label;
  Estimate = new float[vol_size[0]*vol_size[1]*vol_size[2]];
  Label = new float[vol_size[0]*vol_size[1]*vol_size[2]];
  for ( int i = 0; i < vol_size[0] *  vol_size[1] * vol_size[2]; i++ )
  {
    Estimate[i] = 0.0;
    Label[i] = 0.0;
  }


  int ret;
  float mean_neighboors=0;
  pthread_mutex_t mp = PTHREAD_MUTEX_INITIALIZER;
  /* initialize a mutex to its default value */
  ret = pthread_mutex_init ( &mp, NULL );


  for ( ii=0;ii<nb_thread;ii++ )
  {

    arguments[ii].in=in;
    arguments[ii].Estimate=Estimate;
    arguments[ii].Label=Label;
    arguments[ii].filtering_param=filtering_param;
    arguments[ii].beta=beta;
    arguments[ii].search=search;
    arguments[ii].neighborhoodsize=neighborhoodsize;
    arguments[ii].mean_map=mean_map;
    arguments[ii].var_map=var_map;
    arguments[ii].testmean=testmean;
    arguments[ii].testvar=testvar;
    arguments[ii].m_min=m_min;
    arguments[ii].v_min=v_min;
    arguments[ii].vol_size=vol_size;
    arguments[ii].weight_method = weight_method;
    arguments[ii].thread_num=ii;
    arguments[ii].debut= ( int ) ( ii*vol_size[2]/nb_thread );
    arguments[ii].fin= ( int ) ( ( ii+1 ) *vol_size[2]/nb_thread );
    arguments[ii].mp=&mp;
    arguments[ii].b_space=b_space;
    arguments[ii].mean_neighboors=&mean_neighboors;
    arguments[ii].hallucinate=hallucinate;
  }


  if ( debug )
  {
    cout<<"\n Creation of threads\n"<<endl;
    cout<<" Number of slices : "<<vol_size[2]<<"\n"<<endl;
  }
  for ( ii=0;ii<nb_thread;ii++ )
    if ( pthread_create ( &thread[ii],NULL,Sub_denoise_block_mt,&arguments[ii] ) )
    {
      std::cerr<<"Creation de thread impossible"<<std::endl;
      exit ( 1 );
    }

  for ( ii=0;ii<nb_thread;ii++ )
    if ( pthread_join ( thread[ii],&retval ) )
    {
      std::cerr<<"Synchronisation de thread impossible"<<std::endl;
      exit ( 1 );
    }
  free ( thread );
  free ( arguments );

  float label = 0;
  float estimate = 0.0;


  // Aggregation of the estimators (i.e. means computation)
  for ( int kk = 0; kk < vol_size[2]; kk++ )
  {
    for ( int jj = 0; jj < vol_size[1]; jj++ )
    {
      for ( int ii = 0; ii < vol_size[0]; ii++ )
      {
        int addr=kk* ( vol_size[0]*vol_size[1] ) + ( jj*vol_size[0] ) +ii;
        label = Label[addr];
        if ( label == 0 )
        {
          out[addr] = hallucinate?hallucinate[addr]:in[addr];
        }
        else
        {
          estimate = Estimate[addr];
          estimate = ( estimate/label );
          out[addr]=estimate;

        }
      }
    }
  }


  delete [] Estimate;
  delete [] Label;

  if ( verbose )
  {
    cout << "\nMean number of neighboors per voxel: " << mean_neighboors*b_space/vol_size[2] << endl;

    cout << "\nVolume denoising is finished"<< endl;
  }
}


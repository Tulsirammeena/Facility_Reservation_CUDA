#include <iostream>
#include <stdio.h>
#include <cuda.h>

#define max_N 100000
#define max_P 30
#define BLOCKSIZE 1024

using namespace std;

//*******************************************

// Write down the kernels here

/*Initially, I calculate the prefix sum of the facility capacities using the reduction method to facilitate the mapping of center numbers and facility IDs for processing requests. Then, I create threads that correspond to the total number of facilities available. Each thread is assigned to represent a specific facility.*/

/*Subsequently, I create a worklist to accommodate all requests that will utilize a particular facility. I also initialize a slot array that has 25 hours and set all its elements to zero. After that, I check if a facility is available for a specific request, and if it is, I proceed with processing that request. If not, I leave the facility as it is without updating anything.*/

__global__ void  exclusivesum(int n,int *garr , int* ccopy,int *csum) // to calculate exclusive sum of array, calling as kernel
{
    int id=threadIdx.x+ blockDim.x*blockIdx.x;
    int th = threadIdx.x;
    int blc = blockDim.x*blockIdx.x;
   
    if((th + blc)<n)
    {
      int aa = ccopy[th + blc];
      garr[id]= garr[id] - aa;
    }
    if(id < 0)
      return;

     if((threadIdx.x+ blockDim.x*blockIdx.x) == n-1)
      {
        int temp = ccopy[n-1] + garr[n-1];
        csum[0] = temp;
      }
}
__global__ void  dkernel(int n,int *a,int start, int end)
{     
     int th = threadIdx.x;
     
     int id=start+th;
     
     int N,M;
     M=end-start;
     N = M+1;
     if(id >= n)
     {
      return;
     }
  
     if((start+th)<n)
     {
          int  tmp;
          
          
        int off = 1;
        while(off < (M+1))
        {
                
            if(th>=off)
            { 
              int tp = id - off;
                tmp=a[tp];
            }
    
            __syncthreads();
            if(start < 0 && N == th)
            return;
 
           if(th>=off || th < 0)
              a[start+threadIdx.x]= a[id] + tmp;
   
           __syncthreads();
           off = off * 2;
      
        }
          
      
    }
  
}


__global__ void  dke(int *sum, int *gcentre,int *gfacility,int *gcapacity,int *gfac_ids,int *gsucc_reqs,int *gtot_reqs,int *greq_id,int *greq_cen,int *greq_fac,int *greq_start,int *greq_slots,int *gfacps,int R,int N)
{
    if(!(threadIdx.x+1))
    {
      return;
    }
    int thd = threadIdx.x;
    int id=threadIdx.x + blockIdx.x * blockDim.x;
    int blk = blockIdx.x;
    if(id >= sum[0])
    {
      return;
    }
    
    if((blk * blockDim.x + thd)< sum[0])
    {
     
    
    int slot[25];
    int worklist[BLOCKSIZE];
    
    
   
    for(int ii = 0; ii < 25; ii++)
    slot[ii] = 0;
    
    
    
    int centerno=-1,facilityno=-1;
    for(int i = 0; i < N; i++)
    {
      int facid = gfacps[i];
      if( i < slot[0])
      return;
      if(id < facid)
        {
          if(id < sum[0])
          {centerno = i - 1;
          int gg = gfacps[i - 1];
          facilityno =  (blk * blockDim.x + thd) - gg;
          break;}
        }
    }

    if(id < sum[0] && centerno==-1)
     {
        int lastInd = N - 1;
        int aai = gfacps[lastInd];
        centerno=N-1;
        facilityno = id - aai;
     }
     int counter  = 0;
    for(int ii = 0; ii < BLOCKSIZE; ii++)
    worklist[ii] = BLOCKSIZE + 1;
    int j = 0;
    while(j < R)
    {
      int cen_req = greq_cen[j];
      int fac_req = greq_fac[j];
      if(cen_req == centerno && fac_req == facilityno)
      {
        if(greq_fac[j] == facilityno)
          {
            worklist[counter] = greq_id[j];
            counter++;
          }
      }
      j = j+1;
    }
    int nor = counter;
    

    for(int i=0;i<nor;i++)
    {
        int var = 0;
        int check= -1;
        for(int p=greq_start[worklist[i]];p<greq_start[worklist[i]] + greq_slots[worklist[i]];p++)
        {
           if(slot[p]>=gcapacity[id])
            { 
              check=1;
              break;
            } 
        }

        if(check == -1)
        {
          int k = greq_start[worklist[i]];
          while(k < greq_start[worklist[i]] + greq_slots[worklist[i]])
            {
               if(check == -1)
               {
                slot[k]++;
               }
               k++;
            }
           
           atomicAdd(&gsucc_reqs[centerno],1);

        }
        
        var = var + check;

    }
    }

    
}
//***********************************************
int main(int argc,char **argv)
{
	// variable declarations...
    int N,*centre,*facility,*capacity,*fac_ids, *succ_reqs, *tot_reqs;
    

    FILE *inputfilepointer;
    
    //File Opening for read
    char *inputfilename = argv[1];
    inputfilepointer    = fopen( inputfilename , "r");

    if ( inputfilepointer == NULL )  {
        printf( "input.txt file failed to open." );
        return 0; 
    }

    fscanf( inputfilepointer, "%d", &N ); // N is number of centres
	
    // Allocate memory on cpu
    centre=(int*)malloc(N * sizeof (int));  // Computer  centre numbers
    facility=(int*)malloc(N * sizeof (int));  // Number of facilities in each computer centre
    fac_ids=(int*)malloc(max_P * N  * sizeof (int));  // Facility room numbers of each computer centre
    capacity=(int*)malloc(max_P * N * sizeof (int));  // stores capacities of each facility for every computer centre 


    int success=0;  // total successful requests
    int fail = 0;   // total failed requests
    tot_reqs = (int *)malloc(N*sizeof(int));   // total requests for each centre
    succ_reqs = (int *)malloc(N*sizeof(int)); // total successful requests for each centre

    // Input the computer centres data
    int k1=0 , k2 = 0;
    for(int i=0;i<N;i++)
    {
      fscanf( inputfilepointer, "%d", &centre[i] );
      fscanf( inputfilepointer, "%d", &facility[i] );
      
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &fac_ids[k1] );
        k1++;
      }
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &capacity[k2]);
        k2++;     
      }
    }

    // variable declarations
    int *req_id, *req_cen, *req_fac, *req_start, *req_slots;   // Number of slots requested for every request
    
    // Allocate memory on CPU 
	int R;
	fscanf( inputfilepointer, "%d", &R); // Total requests
    req_id = (int *) malloc ( (R) * sizeof (int) );  // Request ids
    req_cen = (int *) malloc ( (R) * sizeof (int) );  // Requested computer centre
    req_fac = (int *) malloc ( (R) * sizeof (int) );  // Requested facility
    req_start = (int *) malloc ( (R) * sizeof (int) );  // Start slot of every request
    req_slots = (int *) malloc ( (R) * sizeof (int) );   // Number of slots requested for every request
    
    // Input the user request data
    for(int j = 0; j < R; j++)
    {
       fscanf( inputfilepointer, "%d", &req_id[j]);
       fscanf( inputfilepointer, "%d", &req_cen[j]);
       fscanf( inputfilepointer, "%d", &req_fac[j]);
       fscanf( inputfilepointer, "%d", &req_start[j]);
       fscanf( inputfilepointer, "%d", &req_slots[j]);
       tot_reqs[req_cen[j]]+=1;  
    }
		


    //*********************************
    // Call the kernels here
    
    int *csum;
    int *greq_id, *greq_cen,*gcapacity;
    int *gfac_ids, *greq_fac, *greq_start, *greq_slots;
    cudaMalloc(&greq_fac,sizeof(int)*R);

    
    int *gcentre,*gfacility, *gsucc_reqs, *gtot_reqs;
    
    cudaMalloc(&csum,sizeof(int));
    cudaMalloc(&greq_id,sizeof(int)*R);
    cudaMalloc(&greq_cen,sizeof(int)*R);
    cudaMalloc(&greq_start,sizeof(int)*R);
    cudaMemcpy(greq_fac, req_fac, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_start, req_start, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc(&gsucc_reqs,sizeof(int)*N);
    cudaMemcpy(greq_id, req_id, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(greq_cen, req_cen, R * sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMalloc(&greq_slots,sizeof(int)*R);
    cudaMalloc(&gcentre,sizeof(int)*N);
    
    
    
    cudaMemcpy(greq_slots, req_slots, R * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gcentre,centre , N* sizeof(int), cudaMemcpyHostToDevice);
    
   
    
    
    cudaMemcpy(gsucc_reqs,succ_reqs ,  N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc(&gfacility,sizeof(int)*N);
    cudaMalloc(&gcapacity,sizeof(int)*max_P * N );
    
    cudaMalloc(&gfac_ids,sizeof(int)*max_P * N );
    
    int *gfacps,*ccopy;
    cudaMalloc(&gfacps,sizeof(int)*(N));

    cudaMalloc(&gtot_reqs,sizeof(int)*N);
    cudaMemcpy(gfacility, facility, N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(gcapacity,capacity , max_P * N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemset(gsucc_reqs, 0, N * sizeof(int));
    cudaMemcpy(gfac_ids, fac_ids, max_P * N * sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc(&ccopy,sizeof(int)*(N));
    
    cudaMemcpy(gtot_reqs,tot_reqs, N * sizeof(int), cudaMemcpyHostToDevice);

    
    
    
    
    cudaMemcpy(ccopy,facility, N * sizeof(int), cudaMemcpyHostToDevice);
    
    int noofblocks= ceil((float)N/1023);
    cudaMemcpy(gfacps,facility, N * sizeof(int), cudaMemcpyHostToDevice);

    int sum[1];
    int i = 0;
    while(i < noofblocks)
    {  
        int start=i*1023;
        int end;
        
        if(i < noofblocks)
        {
        int nn = noofblocks -1;
        if(i==nn)
         end=N-1;
        else
         end= (i+1)*1022;
        }
       
       if(i==0)
       dkernel<<<1,1023>>>(N,gfacps,start,end);
       else
       dkernel<<<1,1023+1>>>(N,gfacps,start-1,end);
       cudaDeviceSynchronize();
       i++;
    }
    if(N)
    {
      exclusivesum<<<noofblocks,BLOCKSIZE>>>(N,gfacps,ccopy,csum);
      cudaDeviceSynchronize();
    }
    cudaMemcpy(sum,csum,sizeof(int),cudaMemcpyDeviceToHost);

    
    int numBlocks = (sum[0] + BLOCKSIZE - 1);
    int final = numBlocks/BLOCKSIZE;
    
    dke<<<final, BLOCKSIZE>>>(csum, gcentre, gfacility, gcapacity, gfac_ids, gsucc_reqs, gtot_reqs, greq_id, greq_cen, greq_fac, greq_start, greq_slots, gfacps, R, N);
    cudaFree(gcentre);
    cudaFree(gfacility);
    cudaFree(gcapacity);
    cudaMemcpy(sum,csum,sizeof(int),cudaMemcpyDeviceToHost);

    if(sum[0] < 0)
    {
      sum[0] = success;
    }
    cudaMemcpy(succ_reqs ,gsucc_reqs,  N * sizeof(int), cudaMemcpyDeviceToHost);

    
    int iii = 0;
    while(iii<N)
     {
        int a = iii++;
        success = success+ succ_reqs[a];

     }
    int ff;
    if(success + 1)
    
    {
      ff = R - success + sum[0];
    } 
    fail = ff - sum[0];

    //********************************
    // Output
    char *outputfilename = argv[2]; 
    FILE *outputfilepointer;
    outputfilepointer = fopen(outputfilename,"w");

    fprintf( outputfilepointer, "%d %d\n", success, fail);
    //printf("%d %d\n", success, fail);
    for(int j = 0; j < N; j++)
    {
        fprintf( outputfilepointer, "%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
        //printf("%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
    }
    fclose( inputfilepointer );
    fclose( outputfilepointer );
    cudaDeviceSynchronize();
	return 0;
}

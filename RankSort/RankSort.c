#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int isActive( int rank, int p, int l );
int isReciever( int rank, int p, int l );
int isSender( int rank, int p, int l );

double * merge_array(int n, double * a, int m, double * b);
void merge_sort(int n, double * a);
void swap (double * a, double * b);

int main (int argc, char *argv[])
{

	int rank, size;

	int n = 160000, q, l, i, j, k, x, *nr;
	double m = 10.0;
	double *a, *b, time1;
	int *ranking, *overallRanking;

	MPI_Status status;

	MPI_Init(&argc, &argv);

    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);

	a = (double *) calloc(n,sizeof(double));
	b = (double *) calloc(n,sizeof(double));
	ranking = (int *) calloc(n/size,sizeof(int)); // the array processor rank generates, n/size processors
	overallRanking = (int *) calloc(n,sizeof(int));

	if( rank == 0 )
	{

	   //initialise the array with random values, then scatter to all processors
	   srand( ((unsigned)time(NULL)+rank) );

	   for( i = 0; i < n; i++ )
	   {
	      a[i]=((double)rand()/RAND_MAX)*m;
	      //printf( "Initial: %f\n", a[i] );
	   }

	}

	time1 = MPI_Wtime();
	
	// Broadcast the array to the processor
	MPI_Bcast(&a[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// P rank generates an array ranking with ranking[i] is the rank of a[i+rank*n/size] in the array
	for(i = 0; i < n/size; i++) {
	  ranking[i] = 0;
	  
	  for(j = 0; j < n; j++) {
	    
	    if( a[i+rank*n/size] > a[j] ) { // scan all elements of a
	      ranking[i]++;
	    }
	    
	  }
	  
	}
	
	// Gather the array ranking to finalRanking
	MPI_Gather(ranking, n/size, MPI_INT, overallRanking, n/size, MPI_INT, 0, MPI_COMM_WORLD);
	
	//if processor 0 then restore the order in the array b
	if(rank == 0) {
	  
	  for(i=0; i<n; i++) {
	      // what is the place of a[i]
	      b[overallRanking[i]] = a[i];
	  }
	  
	}
	
	time1 = MPI_Wtime() - time1;
	
	printf("Processor %d takes %lf sec \n", ran, time1);
	
	
	if( rank == 0 )
	{
	   for( i = 0; i < n; i++ )
	   {
	      //printf( "Output : %f\n", b[i] );
	   }
	}
	MPI_Finalize();

}




// ------------------------------------------------------------
//
// these functions deal with sorting and merging
//

double * merge_array(int n, double * a, int m, double * b){

   int i,j,k;
   double * c = (double *) calloc(n+m, sizeof(double));

   for(i=j=k=0;(i<n)&&(j<m);)
      if(a[i]<=b[j])c[k++]=a[i++];
      else c[k++]=b[j++];
      if(i==n)for(;j<m;)c[k++]=b[j++];
      else for(;i<n;)c[k++]=a[i++];
   return c;
}

void merge_sort(int n, double * a){

   double * c;
   int i;

   if (n<=1) return;

   if(n==2) {

      if(a[0]>a[1])swap(&a[0],&a[1]);
      return;
   }



   merge_sort(n/2,a);merge_sort(n-n/2,a+n/2);

   c=merge_array(n/2,a,n-n/2,a+n/2);

   for(i=0;i<n;i++)a[i]=c[i];

}

void swap (double * a, double * b){

   double temp;

   temp=*a;*a=*b;*b=temp;

}

	

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <mpi.h>
         
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm );
int MPI_Sort_odd_even( int n, double * a, int root, MPI_Comm comm );
int MPI_Is_Sorted(int n, double * a, int root, MPI_Comm comm, int * answer);
double * merge( int n, double * array, int m, double * b );
void merge_sort( int n, double * array );
void swap ( double * array, double * b );
     

int main( int argc, char ** argv ) {
    int size, rank, result, i, *answer;
    int n = 1000000;
    double m = 10.0;
    double *array, processorTime;
    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    array = ( double * ) calloc( n, sizeof( double ) );
       
    if( rank == 0 ){
      srand( ( ( unsigned ) time( NULL ) + rank ) );
      for( i = 0; i < n; i++ ){
        array[ i ] = ( ( double ) rand( ) / RAND_MAX ) * m;
        printf( "Initial: %f\n", array[ i ] );
        }
    }
    processorTime = MPI_Wtime( );
    result = MPI_Sort_odd_even( n, array, 0, MPI_COMM_WORLD );
    if( result != MPI_SUCCESS ){
      return result;
    }
     
    processorTime = MPI_Wtime( ) - processorTime;
    printf( "Processor %d takes %lf sec\n", rank, processorTime );
     if( rank == 0 ){
       for( i = 0; i < n; i++ ){
	printf( "Output : %f\n", array[ i ] );
        }
     }
     MPI_Finalize( );
    }
     
int MPI_Sort_odd_even( int n, double *a, int root, MPI_Comm comm ) {
    int rank, size, i, sorted_result;
     double *local_a;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

       local_a = (double *) calloc(n / size, sizeof(double));
       MPI_Scatter(a, n / size, MPI_DOUBLE, local_a, n / size, MPI_DOUBLE,
		root, comm);
       merge_sort(n / size, local_a);
	 for (i = 1; i <= size; i++) {
           if ((i + rank) % 2 == 0) { 
		if (rank < size - 1) {
		     MPI_Exchange(n / size, local_a, rank, rank + 1, comm);
		    }
	      } else if (rank > 0) {
		    MPI_Exchange(n / size, local_a, rank - 1, rank, comm);
            }
		MPI_Is_Sorted(n/size,local_a,root,comm,&sorted_result);
		if(sorted_result==0){break;}
	 }
	  MPI_Gather(local_a, n / size, MPI_DOUBLE, a, n / size, MPI_DOUBLE,
		   root, comm);

	return MPI_SUCCESS
      
 }
     
int MPI_Exchange( int n, double * array, int rank1, int rank2, MPI_Comm comm ) {
      int rank, size, result, i, tag1 = 0, tag2 = 1;
      double * b = ( double * ) calloc( n, sizeof( double ) );
      double * c;
       
      MPI_Status status;
      MPI_Comm_rank( comm, &rank );
      MPI_Comm_size( comm, &size );
      if( rank == rank1 )
      {
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank2, tag1, comm );
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank2, tag2, comm, &status );
        c = merge( n, array, n, b );
        for( i = 0; i < n; i++ )
        {
          array[ i ] = c[ i ];
        }
      }
      else if( rank == rank2 )
      {
        result = MPI_Recv( &b[ 0 ], n, MPI_DOUBLE, rank1, tag1, comm, &status );
        result = MPI_Send( &array[ 0 ], n, MPI_DOUBLE, rank1, tag2, comm) ;
        c = merge( n, array, n, b );
        for( i =0; i < n; i++ )
        {
          array[ i ] = c[ i + n ];
        }
      }
      return MPI_SUCCESS;
    }
     
double * merge( int n, double * a, int m, double * b ) {
       int i, j, k;
       double * c = ( double * ) calloc( n + m, sizeof( double ) );
     
       for( i=j=k=0; ( i < n ) && ( j < m ); )
       {
          if( a[ i ] <= b[ j ] )
          {
            c[ k++ ] = a[ i++ ];
          }
          else
          {
            c[ k++ ] = b[ j++ ];
          }
       }
      if( i == n )
      {
        for( ; j < m; )
        {
          c[ k++ ] = b[ j++ ];
        }
      }
      else
      {
        for( ; i < n; )
        {
          c[ k++ ] = a[ i++ ];
        }
      }
      return c;
    }
void merge_sort( int n, double * a ) {
      double * c;
      int i;
     
      if ( n <= 1 )
      {
        return;
      }
      if( n == 2 )
      {
        if( a[ 0 ] > a[ 1 ] )
        {
          swap( &a[ 0 ], &a[ 1 ] );
        }
        return;
      }
     
      merge_sort( n / 2, a );
      merge_sort( n - n / 2, a + n / 2 );
      c = merge( n / 2, a, n - n / 2, a + n / 2);
      for( i = 0; i < n; i++ )
      {
        a[ i ] = c[ i ];
      }
    }
     

void swap ( double * a, double * b ) {
       double temp;
       temp = *a;
       *a = *b;
       *b = temp;
    }


int MPI_Is_Sorted(int n, double * a, int root, MPI_Comm comm, int * answer) {
    
	int rank, size, i, result;
	double *local_b, *first, *last;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);
	
	first = (double *) calloc(size,sizeof(double));
	last = (double *) calloc(size,sizeof(double));
	result = MPI_Gather( &a[0], 1, MPI_DOUBLE, first, 1, MPI_DOUBLE, root, comm );
	if(result != MPI_SUCCESS) { return result; }
	result = MPI_Gather( &a[n-1], 1, MPI_DOUBLE, last, 1, MPI_DOUBLE, root, comm );
	if(result != MPI_SUCCESS) { return result; }
	if( rank == 0 ) {
	  *answer = 1;
	  for (i = 1; i < n; i++) {
	      if (first[i] > last[i-1]) {
		  answer = 0;
		  break;
	      }
	  }
	}
	result = MPI_Bcast(&answer, 1, MPI_INT, root, comm);
	if(result != MPI_SUCCESS) { return result; }
	return MPI_SUCCESS;
}


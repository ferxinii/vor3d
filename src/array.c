
#include "array.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) {
        if (arr[ii] == entry) return ii;
    }
    fprintf(stderr, "id_where_equal_int: Could not find id.\n"); 
    exit(1);
}


int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) {
        if (arr1[ii] == a) return 1;
    }
    return 0;
}




#include <assert.h>
#include <stdlib.h>

int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) {
        if (arr[ii] == entry) return ii;
    }
    assert(1 == 0 && "Could not find id."); 
    exit(1);
}


int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) {
        if (arr1[ii] == a) return 1;
    }
    return 0;
}


#ifndef ARRAY_OPERATIONS_C
#define ARRAY_OPERATIONS_C

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


void find_unique_pair(const int *arr1, const int *arr2, int N, int *unique1, int *unique2)
{
    int xorAll = 0; // Compute XOR for all elements in both lists.
    for (int ii=0; ii<N; ii++) {
        xorAll ^= arr1[ii];
        xorAll ^= arr2[ii];
    }
    int diffBit = xorAll & -xorAll; // Find rightmost set bit

    *unique1 = 0; *unique2 = 0;
    for (int ii=0; ii<N; ii++) {  // Partition elements based on the set bit.
        if (arr1[ii] & diffBit) {
            *unique1 ^= arr1[ii];
        } else {
            *unique2 ^= arr1[ii];
        }
        if (arr2[ii] & diffBit) {
            *unique1 ^= arr2[ii];
        } else {
            *unique2 ^= arr2[ii];
        }
    }

    if (!inarray(arr1, N, *unique1)) {  // Correct output, flip 1 and 2 if necessary
        int aux = *unique1;
        *unique1 = *unique2;
        *unique2 = aux;
    }
}

#endif

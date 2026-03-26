/* 
 * Header-only dynamic array (arr) implementation with automatic resizing.
 * Stores items by value with configurable item size. Each arr keeps a 
 * block of memory, and if more space is needed, it doubles its size. If
 * items are removed, the space is kept as is, in case the arr grows again.
 *
 * Copyright (c) 2026 Fernando Muñoz
 * MIT License. See bottom of file.
 */

#ifndef HLIBS_DYNARRAY_H
#define HLIBS_DYNARRAY_H

#include <string.h>
#include <stdlib.h>

#define DYNARRAY_INIT_NMAX 10  /* User can change this */

typedef struct dynarray {
	void  *items;        /* pointer to raw bytes */
	size_t item_size;    /* sizeof(element type), DO NOT INPUT MANUALLY TO AVOID ALIGNMENT ISSUES */
	unsigned int N;
	unsigned int Nmax;
} s_dynarray;


/* INTERFACE */
/* All functions returning int, return 0 on ERROR, 1 on SUCCESS. */
static inline s_dynarray dynarray_initialize(size_t item_size, size_t Nmax);  /* (s_dynarray){0} if ERROR */
static inline void dynarray_free(s_dynarray *arr);
static inline void dynarray_memset0(s_dynarray *arr);  /* Sets allocated memory to 0. */
static inline int dynarray_ensure_capacity(s_dynarray *arr, size_t need);  
static inline int dynarray_push(s_dynarray *arr, const void *elem);  
static inline int dynarray_change_entry(s_dynarray *arr, unsigned id, const void *elem);  
static inline void *dynarray_get_ptr(s_dynarray *arr, size_t id);  /* NULL if ERROR (i.e. id>=N), void* if OK */
static inline int dynarray_get_value(const s_dynarray *arr, size_t id, void *out);  
static inline int dynarray_pop(s_dynarray *arr, void *out_elem);


/* IMPLEMENTATION */
static inline s_dynarray dynarray_initialize(size_t item_size, size_t Nmax)
{
	if (item_size == 0) item_size = 1;
	if (Nmax <= 0) Nmax = DYNARRAY_INIT_NMAX;
	s_dynarray out = { 
        .items     = malloc(Nmax * item_size),
		.item_size = item_size,
		.N         = 0,
		.Nmax      = Nmax,
	};
    if (!out.items) return (s_dynarray){0};
	return out;
}

static inline void dynarray_free(s_dynarray *arr) 
{
	if (!arr) return;
	if (arr->items) free(arr->items);
	memset(arr, 0, sizeof(*arr));
}

static inline void dynarray_memset0(s_dynarray *arr)
{
    memset(arr->items, 0, arr->item_size * arr->Nmax);
}

static inline int dynarray_ensure_capacity(s_dynarray *arr, size_t need) 
{   
	if (!arr) return 0;
	if (need < arr->Nmax) return 1;

	size_t Nmax_new = arr->Nmax ? arr->Nmax : 1;
	while (need >= Nmax_new) Nmax_new *= 2;

	void *tmp = realloc(arr->items, Nmax_new * arr->item_size);
	if (!tmp) return 0;

	arr->items = tmp;
	arr->Nmax = Nmax_new;
	return 1;
}

static inline int dynarray_push(s_dynarray *arr, const void *elem) 
{
	if (!arr) return 0;
	if (!dynarray_ensure_capacity(arr, arr->N + 1)) return 0;
	void *dst = (char *)arr->items + arr->N * arr->item_size;
	memcpy(dst, elem, arr->item_size);
	arr->N += 1;
	return 1;
}

static inline int dynarray_change_entry(s_dynarray *arr, unsigned id, const void *elem)
{
    if (!arr || id >= arr->N) return 0;

    void *dst = (uint8_t*)arr->items + id * arr->item_size;
    memmove(dst, elem, arr->item_size);
    return 1;
}

static inline void *dynarray_get_ptr(s_dynarray *arr, size_t id)
{
	if (!arr || id >= arr->N) return NULL;
	return (uint8_t*)arr->items + id * arr->item_size;
}

static inline int dynarray_get_value(const s_dynarray *arr, size_t id, void *out)
{
	if (!arr || !out || id >= arr->N) return 0;

	const void *src = (uint8_t*)arr->items + id * arr->item_size;
	memmove(out, src, arr->item_size);
	return 1;
}

static inline int dynarray_pop(s_dynarray *arr, void *out_elem)
{
	if (!arr || arr->N == 0) return 0;
	arr->N -= 1;
	void *src = (uint8_t*)arr->items + arr->N * arr->item_size;
	if (out_elem) memcpy(out_elem, src, arr->item_size);
	return 1;
}


#endif

/* MIT License.
 *
 * Copyright (c) 2026 Fernando Muñoz.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


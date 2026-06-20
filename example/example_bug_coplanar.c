/* Reproducer for the coplanar-input DT crash documented in ERROR.md.
 *
 * Trigger: 8 cube-corner mesh vertices (all on z=0 or z=1 planes) followed
 * by a 5×5×5 interior grid of 125 points, fed together to construct_dt_3d.
 * Expected: crash / assertion failure at point_id=12 (first interior point
 * inserted after all 8 coplanar corners occupy indices 4–11). */

#include <stdio.h>
#include <stdlib.h>
#include "delaunay.h"

int main(void)
{
    /* 8 cube corners — 4 on z=0, 4 on z=1 (two perfectly coplanar layers) */
    s_point mesh_verts[8] = {
        {{{0.0, 0.0, 0.0}}},
        {{{1.0, 0.0, 0.0}}},
        {{{0.0, 1.0, 0.0}}},
        {{{1.0, 1.0, 0.0}}},
        {{{0.0, 0.0, 1.0}}},
        {{{1.0, 0.0, 1.0}}},
        {{{0.0, 1.0, 1.0}}},
        {{{1.0, 1.0, 1.0}}},
    };

    /* 5×5×5 = 125 interior grid points, away from the faces */
    const int G = 5;
    s_point grid[125];
    int idx = 0;
    for (int ix = 1; ix <= G; ix++)
        for (int iy = 1; iy <= G; iy++)
            for (int iz = 1; iz <= G; iz++)
                grid[idx++] = (s_point){{{ix / (double)(G + 1),
                                          iy / (double)(G + 1),
                                          iz / (double)(G + 1)}}};

    /* Combine: mesh vertices first (will become DT point-ids 4–11),
     * interior grid second (first at point-id 12). */
    const int N = 8 + 125;
    s_point *pts = malloc((size_t)N * sizeof(s_point));
    if (!pts) { fputs("OOM\n", stderr); return 1; }
    for (int i = 0; i < 8;   i++) pts[i]     = mesh_verts[i];
    for (int i = 0; i < 125; i++) pts[8 + i] = grid[i];

    s_points points = {.N = N, .p = pts};

    printf("Calling construct_dt_3d with %d points (%d coplanar mesh verts + %d interior)...\n",
           N, 8, 125);
    fflush(stdout);

    s_scplx dt = construct_dt_3d(&points, NULL, false, 1e-12, NULL);

    if (!dt.head) {
        fputs("construct_dt_3d returned an empty complex (error or crash).\n", stderr);
        free(pts);
        return 1;
    }

    printf("construct_dt_3d succeeded — %d tetrahedra.\n", dt.N_ncells);
    free_complex(&dt);
    free(pts);
    return 0;
}

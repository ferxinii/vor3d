// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "points.h"
#include "convh.h"
#include "scplx.h"
#include "vor3d.h"
#include "volsph.h"
#include "convh.h"
#include "./hlibs/dynarray.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>



void free_vcell(s_vcell *c)
{
    for (int i = 0; i < c->N_pieces; i++) free_convhull(&c->pieces[i]);
    free(c->pieces);
    for (int i = 0; i < c->N_surface; i++) free_trimesh(&c->surface[i]);
    free(c->surface);
    memset(c, 0, sizeof(s_vcell));
}


/* Memory of vdiagram, printing, and other utilities. */
bool vdiagram_is_valid(const s_vdiagram *vd)
{
    return vd->vcells != NULL;
}

void free_vdiagram(s_vdiagram *vd)
{
    if (vd->vcells) {
        for (int i = 0; i < vd->seeds.N; i++) free_vcell(&vd->vcells[i]);
        free(vd->vcells);
    }
    free_points(&vd->seeds);
    free(vd->adj_offset);
    free(vd->adj_neighbor);
    memset(vd, 0, sizeof(s_vdiagram));
}

/* Debugger helpers (e.g. lldb: call print_vdiagram(&vd)); intentionally in no
 * header, and not static so they are kept when unused. */
void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d, N_pieces = %d, vol = %f\n", vcell->seed_id, vcell->N_pieces, vcell->volume);
    if (vcell->N_pieces > 0) {
        printf("Piece 0: Nv = %d\n", vcell->pieces[0].points.N);
        printf("%f, %f, %f\n", vcell->pieces[0].points.p[0].x,
               vcell->pieces[0].points.p[0].y, vcell->pieces[0].points.p[0].z);
        for (int ii=1; ii<vcell->pieces[0].points.N; ii++) {
            printf("%f, %f, %f\n", vcell->pieces[0].points.p[ii].x,
                   vcell->pieces[0].points.p[ii].y, vcell->pieces[0].points.p[ii].z);
        }
    }
    printf("\n");
}

void print_vdiagram(const s_vdiagram *vdiagram)
{
    puts("------- VORONOI DIAGRAM ------");
    for (int ii=0; ii<vdiagram->seeds.N; ii++)
        print_vcell(&vdiagram->vcells[ii]);
    puts("------------------------------");
}

int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double tol_bdry)
{
    for (int ii=0; ii<vd->seeds.N; ii++) {
        for (int p = 0; p < vd->vcells[ii].N_pieces; p++) {
            e_geom_test test = test_point_in_convhull(&vd->vcells[ii].pieces[p], x, tol_bdry);
            if (test == TEST_DEGENERATE || test == TEST_ERROR) fprintf(stderr, "find_inside_which_vcell: TEST_DEGENERATE or TEST_ERROR. Continuing.\n");
            if (test == TEST_IN || test == TEST_BOUNDARY) return ii;
        }
    }
    return -1;
}

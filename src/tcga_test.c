#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "poly.h"


int main(int argc, char ** argv)
{

    //char file[] = "/home/erikw/projects/NucleAI2/src/features.raw";
    char file[] = "/tmp/blca_polygon/TCGA-2F-A9KO-01Z-00-DX1.195576CF-B739-4BD9-B15B-4A70AE287D3E.svs/features.raw";

    FILE * fid = fopen(file, "r");

    double * P = malloc(10000*sizeof(double));

    double * row = malloc(3*sizeof(double));
    size_t id = 0;
    int curr_id = 0;

    size_t ppos = 0;

    size_t nread = fread(row, 8, 3, fid);
    while(curr_id < 1000005)
    {
        ppos = 0;
        if(nread != 3)
        {
            printf("Failed to read file\n");
        }

        while(id == curr_id)
        {
            // printf("%f, %f, %f\n", row[0], row[1], row[2]);

            P[2*ppos] = row[1];
            P[2*ppos+1] = row[2];

            ppos++;
            nread = fread(row, 8, 3, fid);
            if(nread != 3)
            {
                printf("Failed to read file\n");
            }
            id = row[0];
        }
        if(curr_id % 1000 == 0)
        {
            printf(" -> Polygon %d\n", curr_id);
        }
        poly_reverse(P, ppos);
        poly_props * props = poly_measure(P, ppos);
            poly_props_print(stdout, props);

        poly_props_free(&props);
        curr_id++;
    }
    fclose(fid);
}

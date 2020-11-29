#include "fdtd3d.h"

double*** memory_allocate3d( int d1, int d2, int d3, double ini )
{
    double ***array;
    array = new double**[d1];

    for( int i = 0; i < d1; i++ ){
        array[i] = new double*[d2];
        for( int j = 0; j < d2; j++ ){
            array[i][j] = new double [d3];
            for( int k = 0; k < d3; k++ ){
                array[i][j][k] = ini;
            }
        }
    }

    return array;
}
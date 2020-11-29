#include "pml.h"

void pml::set_point_1(int y, int z)
{
    p_j1 = y;
    p_k1 = z;
}

void pml::set_point_2(int y, int z)
{
    p_j2 = y;
    p_k2 = z;
}

void pml::set_point( int v_y1, int v_y2, int v_z1, int v_z2 )
{
    set_point_1(v_y1, v_z1);
    set_point_2(v_y2, v_z2);
}
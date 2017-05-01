#include "Points_3D.h"


Points_3D::Points_3D(void)
{

}

Points_3D::Points_3D(double x1, double y1,double z1)
{
	x=x1;
	y=y1;
    z=z1;
}

Points_3D::~Points_3D(void)
{
}

double Points_3D::get_pos_x ()
{
	return x;
}
double Points_3D::get_pos_y ()
{
	return y;
}
double Points_3D::get_pos_z ()
{
	return z;
}
void Points_3D::set_pos_x (double x1)
{
	x=x1;
}
 void Points_3D::set_pos_y (double y1)
{
	y=y1;
}
 void Points_3D::set_pos_z (double z1)
{
	z=z1;
}

#include "Points_1D.h"


Points_1D::Points_1D(void)
{

}

Points_1D::Points_1D(double x1)
{
	x=x1;
}

Points_1D::~Points_1D(void)
{
}

double Points_1D::get_pos_x ()
{
	return x;
}

void Points_1D::set_pos_x (double x1)
{
	x=x1;
}


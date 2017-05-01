#pragma once
class Points_2D
{

private:
	double x;
	double y;

public:
	Points_2D();
	Points_2D(double,double);
	~Points_2D(void);

	double get_pos_x ();
	double get_pos_y ();

	void set_pos_x (double);
	void set_pos_y (double);
	
};


#pragma once
class Points_3D
{

private:
	double x;
	double y;
	double z;

public:
	Points_3D();
	Points_3D(double,double,double);
	~Points_3D(void);

	double get_pos_x ();
	double get_pos_y ();
	double get_pos_z ();

	void set_pos_x (double);
	void set_pos_y (double);
	void set_pos_z (double);
	
};


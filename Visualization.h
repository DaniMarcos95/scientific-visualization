#ifndef VISUALIZATION_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define VISUALIZATION_H

void rainbow(float,float*,float*,float*);
double RGB_to_lin(double x);
double lin_to_RGB( double y);
void grayscale(float value, float* R,float* G,float* B);

void direction_to_color(float, float, int);
void visualize(void);
void set_colormap(float , int , int , int);
#endif
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GL/glut.h>   
#include "Simulation.h"
#include <rfftw.h>              //the numerical simulation FFTW library

using namespace std;

const int DIM = 50;
extern fftw_real *vx, *vy; 
extern fftw_real *rho;
void visualize(void);
int   winWidth, winHeight;	//size of the graphics window, in pixels 
int   color_dir = 0;           //use direction color-coding or not 
float vec_scale = 1000;			//scaling of hedgehogs 
int   draw_smoke = 0;           //draw the smoke or not 
int   draw_vecs = 1;            //draw the vector field or not 
const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, banded
const int COLOR_RAINBOW=1;
const int COLOR_GRAYSCALE=2;
int   scalar_col = 0;           //method for scalar coloring 
int   frozen = 0;               //toggles on/off the animation 
void rainbow(float,float*,float*,float*);
void set_colormap(float);
void direction_to_color(float, float, int);
double RGB_to_lin(double x);
double lin_to_RGB( double y);
void grayscale(float value, float* R,float* G,float* B);
void draw_color_legend();
void set_colors_in_bar(int color, int x1, int x2);
void rainbow_bar();
void blue_to_yellow_bar();

//rainbow: Implements a color palette, mapping the scalar 'value' to a rainbow color RGB
void rainbow(float value,float* R,float* G,float* B)
{
   const float dx=0.8;
   if (value<0) value=0; 
   if (value>1) value=1;
   value = (6-2*dx)*value+dx;
   *R = max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
   *G = max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
   *B = max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}

//function that converts RGB to linear
double RGB_to_lin(double x) {
	if (x < 0.04045) return x/12.92;
	return pow((x+0.055)/1.055,2.4);
}

//Inverse of RGB_to_lin
double lin_to_RGB( double y) {
	if(y <= 0.003108) return 12.92 * y;
	return 1.055 * pow(y, 1/2.4) - 0.055;
}
//greyscale: Implements a color palette, mallping the scalar 'value' to a greyscale
void grayscale(float value, float* R,float* G,float* B)
{
	const float dx=0.4;
	
	
	if (value<0) value=0; 
	if (value<1) value=1;
	value = (6-2*dx)*value+dx;
	//Get RGB Values
	*R = max(0.0,(3-fabs(value-4) - fabs(value-5))/2);
	*G = max(0.0,(4-fabs(value-2) - fabs(value-4))/2);
	*B = max(0.0,(3-fabs(value-1) - fabs(value-2))/2);
	
	//Convert RGB values to linear grayscale and get the gray value
	double Rlin = 0.2126 * RGB_to_lin(*R/255.0);
	double Glin = 0.7152 * RGB_to_lin(*G/255.0);
	double Blin = 0.0722 * RGB_to_lin(*B/255.0);
	double gray_linear = 0.2126 * Rlin + 0.7152 * Glin + 0.0722 * Blin;
	
	double gray_color = round(lin_to_RGB(gray_linear) * 255);
	
	//Gray value means equal values for R,G,B.
	*R = *G = *B = gray_color;
	
}

//set_colormap: Sets three different types of colormaps
//set_colormap: Sets three different types of colormaps
void set_colormap(float vy)
{
   float R,G,B;

   if (scalar_col==COLOR_BLACKWHITE)
       R = G = B = vy;
   else if (scalar_col==COLOR_RAINBOW)
       rainbow(vy,&R,&G,&B);
   else if (scalar_col==COLOR_GRAYSCALE)
	   
	   grayscale(vy,&R,&G,&B);
	   
       {
          const int NLEVELS = 7;
          vy *= NLEVELS; vy = (int)(vy); vy/= NLEVELS;
	      rainbow(vy,&R,&G,&B);
	   }

   glColor3f(R,G,B);
}


//direction_to_color: Set the current color by mapping a direction vector (x,y), using
//                    the color mapping method 'method'. If method==1, map the vector direction
//                    using a rainbow colormap. If method==0, simply use the white color
void direction_to_color(float x, float y, int method)
{
	float r,g,b,f;
	if (method == 1)
	{
	  f = atan2(y,x) / 3.1415927 + 1;
	  r = f;
	  if(r > 1) r = 2 - r;
	  g = f + .66667;
      if(g > 2) g -= 2;
	  if(g > 1) g = 2 - g;
	  b = f + 2 * .66667;
	  if(b > 2) b -= 2;
	  if(b > 1) b = 2 - b;
	}
	else if (method == 0)
	{ r = g = b = 1; }
	else {f = atan2(y,x) / 3.1415927 + 1;
		  r = f;
		  if(r > 1) r = 2 - r;
		  g = f + .66667;
		  if(g > 2) g -= 2;
		  if(g > 1) g = 2 - g;
		  b = f + 2 * .66667;
		  if(b > 2) b -= 2;
		  if(b > 1) b = 2 - b;
		  float gr = (r+g+b)/3.0;
		  r = g = b = gr;
		
	}
	glColor3f(r,g,b);
	
}

//visualize: This is the main visualization function
void visualize(void)
{
	int        i, j, idx; double px,py;
	fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
	fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

	if (draw_smoke)
	{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_TRIANGLE_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;
		glColor3f(rho[idx],rho[idx],rho[idx]);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;
			set_colormap(rho[idx]);
			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			set_colormap(rho[idx]);
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		set_colormap(rho[idx]);
		glVertex2f(px, py);
		glEnd();
	}
	}

	if (draw_vecs)
	{
	  glBegin(GL_LINES);				//draw velocities
	  for (i = 0; i < DIM; i++)
	    for (j = 0; j < DIM; j++)
	    {
		  idx = (j * DIM) + i;
		  direction_to_color(vx[idx],vy[idx],color_dir);
		  glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
		  glVertex2f((wn + (fftw_real)i * wn) + vec_scale * vx[idx], (hn + (fftw_real)j * hn) + vec_scale * vy[idx]);
	    }
	  glEnd();
	}

	draw_color_legend();
}

void draw_color_legend(){
	switch (color_dir)
	{
	  case 0:  break;

	  case 1:  
	  	rainbow_bar();
	  	break;
	  case 2:  
	  	blue_to_yellow_bar();
	  	break;
	}
	
	
	
}

void rainbow_bar(){
	
		// glBegin(GL_QUAD_STRIP);
	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0.4,0.2,1); //Purple
	glVertex2f(0,0);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(100,0);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(100,30);
	glColor3f(0.4,0.2,1); //Purple
	glVertex2f(0,30);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(100,0);
	glColor3f(0.2,1,0);//Green
	glVertex2f(200,0);
	glColor3f(0.2,1,0);//Green
	glVertex2f(200,30);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(100,30);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0.2,1,0);//Green
	glVertex2f(200,0);
	glColor3f(1,1,0); //Yellow
	glVertex2f(300,0);
	glColor3f(1,1,0); //Yellow
	glVertex2f(300,30);
	glColor3f(0.2,1,0);//Green
	glVertex2f(200,30);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(1,1,0); //Yellow
	glVertex2f(300,0);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(400,0);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(400,30);
	glColor3f(1,1,0); //Yellow
	glVertex2f(300,30);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(400,0);
	glColor3f(1,0,0); //Red
	glVertex2f(500,0);
	glColor3f(1,0,0); //Red
	glVertex2f(500,30);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(400,30);
	glEnd();
}

void blue_to_yellow_bar(){
	
		// glBegin(GL_QUAD_STRIP);
	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0.6,1); 
	glVertex2f(0,0);
	glColor3f(1,1,0); 
	glVertex2f(500,0);
	glVertex2f(500,30);
	glColor3f(0,0.6,1); 
	glVertex2f(0,30);
	glEnd();

}
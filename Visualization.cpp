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
// const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, banded
const int COLOR_RAINBOW=1;
const int COLOR_GRAYSCALE=0;
const int COLOR_BLUEYEL =2;
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
void grayscale_bar();
float min(float x, float y);
int NCOLORS = 255;

//convert RGB values to HSV
void rgb2hsv(float r, float g, float b,
			float& h, float& s, float& v)
{
	float M = max(r,max(g,b));
	float m = min(r,min(g,b));
	float d = M-m;
	v = M;
	s = (M>0.00001)? d/M:0; //saturation
	if (s==0) h = 0;
	else
	{
		if (r==M) h = (g-b)/d;
		else if (g==M) h = 2 + (b-r)/d;
		else h = 4 + (r-g)/d;
		h /= 6;
		if (h <0) h += 1;
	}
	}

// Convert from HSV to RGB
void hsv2rgb(float h, float s, float v, float& r, float& g, float& b)
{
	int hueCase = (int)(h * 6);
	float frac = 6*h-hueCase;
	float lx = v*(1-s);
	float ly = v*(1 - s*frac);
	float lz = v*(1-s*(1-frac));
	switch (hueCase)
	{
		case 0:
		case 6: r=v; g=lz; b=lx; break; //0<hue<1/6
		case 1: r=ly; g=v; b=lx; break; //1/6<hue<2/6
		case 2: r=lx; g=v; b=lz; break; 
		case 3: r=lx; g=ly; b=v; break;
		case 4: r=lz; g=lx; b=v; break;
		case 5: r=v; g=lx; b=ly; break;
		
	}
}
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
	//const float dx=0.8;
	//const float gamma = 2.2;
	//float Y,L;
	float h,s,v;
	if (value<0) value=0; 
	if (value<1) value=1;
	
	//Get RGB Values
	*R = value/3;
	*G = value/3;
	*B = value/3;
	
	rgb2hsv(*R,*G,*B, h,s,v);
	s = 0;
	
	hsv2rgb(h,s,v,*R,*G,*B);
}

//Yellow-Blue Color Spectrum
void blue_yel(float value, float* R, float* G, float* B)
{
   const float dx=0.8;
   //value = (6-2*dx)*value+dx;
   if (value<0) value=0; 
   if (value>1) value=1;
   value = 6*value; //set value to [0,6] range
   
   
   
   *R = max(0.0,(3-fabs(value-6)));
   *G = max(0.0,(6-fabs(value-3)-fabs(value-6))/2);
   *B = max(0.0,(3-fabs(value)));
   /*if (*B>*R && *B>*G) 
   { if (*R>*B) *B=0;
   else *R=0;}
   
   if (*R > *G) (*R = *G);*/
}

//set_colormap: Sets three different types of colormaps
void set_colormap( float value, int scalar_col, int NCOLORS)
{
   float R,G,B;
	
	value *= NCOLORS;
	value = (int)(value); 
	value/= NCOLORS;
   if (scalar_col==0) //WHITE
       {R = value;
	   G = value;
	   B = value;}
   else if (scalar_col==1) //RAINBOW
       {
		   
		   rainbow(value,&R,&G,&B);}      
	   
	else	//grayscale
		blue_yel(value,&R,&G,&B);
	   
       

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
	else if( method == 2)
		{f = atan2(y,x) / 3.1415927 + 1;
		r = f;
		if (r >1 ) {r = 2 -r;}
		g = f + .66667;
		
	    b = f + 2 * .66667;
		if(b>2) {b-= 2;}
		if(b>1) {b = 2 -b;}
		
		if(r>0 && b>0) {
			if(r>b) b=0;
		else r=0;}
		
		if (r > g) r=0;
		}
		
	
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
		  r = gr;
		  g = 0.924 * gr;
		  b = 0.964 * gr;
		
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
			set_colormap(rho[idx], scalar_col,NCOLORS);
			//direction_to_color(vx[idx],vy[idx],color_dir);
			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			//direction_to_color(vx[idx],vy[idx],color_dir);
			set_colormap(rho[idx], scalar_col,NCOLORS);
			//printf("%f\n",rho[idx]);
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		
		set_colormap(rho[idx],scalar_col,NCOLORS);
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
	  case 0:  
	  	grayscale_bar();
	  	break;
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
	glColor3f(0,0.0,1); 
	glVertex2f(0,0);
	glColor3f(0.2,1,0); 
	glVertex2f(250,0);
	glVertex2f(250,30);
	glColor3f(0,0.0,1); 
	glVertex2f(0,30);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0.2,1,0); 
	glVertex2f(250,0);
	glColor3f(1,1,0); 
	glVertex2f(500,0);
	glVertex2f(500,30);
	glColor3f(0.2,1,0); 
	glVertex2f(250,30);
	glEnd();

}

void grayscale_bar(){
	
		// glBegin(GL_QUAD_STRIP);
	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0,0); 
	glVertex2f(0,0);
	glColor3f(1,1,1); 
	glVertex2f(500,0);
	glVertex2f(500,30);
	glColor3f(0,0,0); 
	glVertex2f(0,30);
	glEnd();

}
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GL/glut.h>   
#include "Simulation.h"
#include <rfftw.h>              //the numerical simulation FFTW library
#include <string>

using namespace std;
int DIM = 50;
extern fftw_real *vx, *vy; 
extern fftw_real *fx, *fy;	
extern fftw_real *rho;
void visualize(void);
int   winWidth, winHeight;	//size of the graphics window, in pixels 
int   color_dir = 1;           //use direction color-coding or not 
float vec_scale = 10;			//scaling of hedgehogs 
int   draw_rho = 0;           //draw the smoke or not 
int draw_for = 0;
int   draw_vecs = 0;            //draw the vector field or not 
int draw_vec_mod = 0;
int draw_for_mod = 0;
int diver = 0;
int temp_idx;
// const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, banded
const int COLOR_RAINBOW=1;
const int COLOR_GRAYSCALE=0;
const int COLOR_BLUEYEL =2;
int   scalar_col = 0;           //method for scalar coloring 
int   frozen = 0;               //toggles on/off the animation 
void rainbow(float,float*,float*,float*,int color_bar);
void set_colormap(float , int , int, int color_bar );
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
void compute_codes();
void draw_isolines();

float diverg; 
float array[50][50];
float array2[50][50];
char isocodes[50][50];

void render();
const float dx=0.8;
int NCOLORS = 255;
int dataset_index;
extern float max_rho;
extern float min_rho;
extern float max_f;
extern float min_f;
extern float max_v;
extern float min_v;
extern float max_clamped;
extern float min_clamped;
fftw_real  hn;
extern int scalarIndex;
extern int glyphIndex;
extern int vectorIndex;
extern int numberOfSamples;
extern float isovalue;
extern int isolineselected;






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
void rainbow(float value,float* R,float* G,float* B, int color_bar)
{
    const float dx=0.8;
	
	if(color_bar != 1){
   		if (value<min_clamped) value=min_clamped; 
   		if (value>max_clamped) value=max_clamped;
   		value = (value - min_clamped) / (max_clamped - min_clamped);
    }
	value = (6-2*dx)*value+dx;
	
	*R = max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
	*G = max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
	*B = max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}


//Yellow-Blue Color Spectrum
void blue_yel(float value, float* R, float* G, float* B, int color_bar)
{
   
   if(color_bar != 1){
   		if (value<min_clamped) value=min_clamped; 
   		if (value>max_clamped) value=max_clamped;
   		value = (value - min_clamped) / (max_clamped - min_clamped);
   }
   
   value = 6*value; //set value to [0,6] range
    
   *R = max(0.0,(3-fabs(value-6)));
   *G = max(0.0,(6-fabs(value-3)-fabs(value-6))/2);
   *B = max(0.0,(3-fabs(value)));
   
}

//set_colormap: Sets three different types of colormaps
void set_colormap( float value, int scalar_col, int NCOLORS, int color_bar)
{
	float R,G,B;

	value *= NCOLORS;
	value = (int)(value); 
	value/= NCOLORS;
	if (scalar_col==0){ //WHITE
		if(color_bar != 1){
	   		if (value<min_clamped) value=min_clamped; 
	   		if (value>max_clamped) value=max_clamped;
	   		value = (value - min_clamped) / (max_clamped - min_clamped);
  		}	
	   	R = value;
	   	G = value;
	   	B = value;}
	else if (scalar_col==1) //RAINBOW
	   {
		   
		   rainbow(value,&R,&G,&B, color_bar);}      
	   
	else	//blue-yellow colorbar
		blue_yel(value,&R,&G,&B, color_bar);
	   
	   

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
	  set_colormap(f,scalar_col,NCOLORS,1);
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

 float scaling_clamping(float value,float maxi, float mini){
	//value = NCOLORS * (value - mini) / (maxi - mini);
 	if (value<mini) value=mini; 
	if (value>maxi) value=maxi;
 	//value = (value - mini)/(maxi - mini);
	//printf("%d\n",NCOLORS);
 	return value;
 }

//visualize: This is the main visualization function
void visualize()
{	
	int        i, j, idx; double px,py; float vec_mod, for_mod;
	fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
	hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

	if (draw_rho)
	{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_QUAD_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;

		set_colormap(rho[idx], scalar_col,NCOLORS,0);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;

			set_colormap(rho[idx], scalar_col,NCOLORS,0);

			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);

			set_colormap(rho[idx], scalar_col,NCOLORS,0);
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		set_colormap(rho[idx],scalar_col,NCOLORS,0);
		glVertex2f(px, py);
		glEnd();
	}
	}

	//Drawing the smoke module ||v||- not working at the moment.
	if (draw_vec_mod)
	
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_QUAD_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;
		
		 
		
		vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
		set_colormap(vec_mod, scalar_col,NCOLORS,1);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;
			
			 
			vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
			set_colormap(vec_mod, scalar_col,NCOLORS,1);
			

			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			 
			vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
			set_colormap(vec_mod, scalar_col,NCOLORS,1);
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		
		 
		vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
		set_colormap(vec_mod, scalar_col,NCOLORS,1);
		glVertex2f(px, py);
		glEnd();
	}
	

	if (draw_for_mod) //Draws smoke for the force dataset
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


	for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_QUAD_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;

		for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
		set_colormap(for_mod, scalar_col,NCOLORS,0);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;

			for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
			set_colormap(for_mod, scalar_col,NCOLORS,0);

			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);

			for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
			set_colormap(for_mod, scalar_col,NCOLORS,0);
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
		set_colormap(for_mod, scalar_col,NCOLORS,0);
		glVertex2f(px, py);
		glEnd();
	}
	

	if (draw_vecs){

		glEnable(GL_LIGHTING);
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_LIGHT0);
		
		if(glyphIndex == 0){
			glBegin(GL_LINES);
		}

		for (i = 0; i < DIM-1; i += DIM/numberOfSamples){
	    	for (j = 0; j < DIM-1; j += DIM/numberOfSamples){	

	    		px = wn + (fftw_real)i * wn;
				py = hn + (fftw_real)(j + 1) * hn;

				idx = ((j + 1) * DIM) + i;
				float length = 0;
				switch(scalarIndex){
					case 0: length = rho[idx];
							break;
					case 1: vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
							// length = (vec_mod - min_v)/(max_v - min_v);
							length = vec_mod;
							break;
					case 2: for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
							// length = (for_mod - min_f)/(max_f - min_f);
							length = for_mod;
							break;
				}
				set_colormap(length, scalar_col,NCOLORS,0);
				
				
				float radius = length*vec_scale;
				float angle = atan2(fy[idx],fx[idx])*180/3.1415	;
				GLUquadricObj *quadric = gluNewQuadric();

				switch(glyphIndex){
					case 0:	direction_to_color(vx[idx],vy[idx],scalar_col);
							
						    glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
						    glVertex2f((wn + (fftw_real)i * wn) + vec_scale * vx[idx], (hn + (fftw_real)j * hn) + vec_scale * vy[idx]);
						    
							break;

					case 1: glTranslatef(px, py,1);						
							glRotatef(90-angle,0,0,-1);
							glRotatef(-90,1,0,0);

							glScalef(1.0,1.0,1.0);

							gluCylinder(quadric, radius/3, 0, radius, 4, 4);

							gluDeleteQuadric(quadric);
			                glLoadIdentity();

			                break;

			        case 2:	glTranslatef(px, py,1);						
							glRotatef(90-angle,0,0,-1);
							glRotatef(-90,1,0,0);

							// ARROW
							glTranslated(0,0,radius);
							gluCylinder(quadric, radius/3, 0, radius, 4,4);
							glTranslated(0,0,-radius);
							gluCylinder(quadric, radius/20, radius/20, radius, 4, 4);

							//Line
							glTranslatef(px, py,1);						
							glRotatef(90-angle,0,0,-1);
							glRotatef(-90,1,0,0);
							gluDeleteQuadric(quadric);
							glLoadIdentity();
				}
			}
		}	

		if(glyphIndex == 0){
			glEnd();
		}
	}
	
	//Draw glyphs for the f vector field
	if (draw_for){ 

		glEnable(GL_LIGHTING);
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_LIGHT0);
		
		for (i = 0; i < DIM-1; i++){
	    	for (j = 0; j < DIM-1; j++){	

	    		px = wn + (fftw_real)i * wn;
				py = hn + (fftw_real)(j + 1) * hn;

				idx = ((j + 1) * DIM) + i;

				float length = 0;
				switch(scalarIndex){
					case 0: length = rho[idx];
							break;
					case 1: vec_mod = sqrt(pow(vx[idx],2) + pow(vy[idx],2));
							// length = (vec_mod - min_v)/(max_v - min_v);
							length = vec_mod;
							break;
					case 2: for_mod = sqrt(pow(fx[idx],2) + pow(fy[idx],2));
							// length = (for_mod - min_f)/(max_f - min_f);
							break;
				}

				set_colormap(length, scalar_col,NCOLORS,0);

				float radius = length*vec_scale;
				float angle = atan2(fy[idx],fx[idx])*180/3.1415	;

				GLUquadricObj *quadric = gluNewQuadric();

				// glTranslatef(px, py,1);						
				// glRotatef(90-angle,0,0,-1);
				// glRotatef(-90,1,0,0);

				// glScalef(1.0,1.0,1.0);

				// gluCylinder(quadric, radius/3, 0, radius, 4, 4);

				// gluDeleteQuadric(quadric);
    //             glLoadIdentity();






             	glTranslatef(px, py,1);						
				glRotatef(90-angle,0,0,-1);
				glRotatef(-90,1,0,0);

				// ARROW
				glTranslated(0,0,radius);
				gluCylinder(quadric, radius/3, 0, radius, 4,4);
				glTranslated(0,0,-radius);
				gluCylinder(quadric, radius/20, radius/20, radius, 4, 4);

				//Line
				glTranslatef(px, py,1);						
				glRotatef(90-angle,0,0,-1);
				glRotatef(-90,1,0,0);
				gluDeleteQuadric(quadric);
				glLoadIdentity();
			}
		}	

	}		
	
	
	if (diver){
		switch(vectorIndex){
			case 0: //V vector field
		for(int rowNumber = 0; rowNumber <= DIM-1; ++rowNumber)
		{
			
		
		
		for(int colNumber = 0; colNumber <= DIM-1; ++colNumber)
		{	if(rowNumber == DIM-1){
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber)*DIM + colNumber];
		array2[rowNumber][colNumber] = vy[colNumber] -vy[(rowNumber)*DIM + colNumber];}
		else{
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber*DIM) +colNumber];
		array2[rowNumber][colNumber] = vy[DIM*(rowNumber+1) + colNumber] - vy[(rowNumber*DIM) +colNumber];}
		}
		}
		
		for(int rowNumber = 0; rowNumber <= DIM-1; ++rowNumber)
		{
			
		
		
		for(int colNumber = 0; colNumber <= DIM-1; ++colNumber)
		{	if(rowNumber == DIM-1){
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber)*DIM + colNumber];
		array[rowNumber][colNumber] = vx[DIM*rowNumber] -vx[(rowNumber)*DIM + colNumber];}
		else{
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber*DIM) +colNumber];
		array[rowNumber][colNumber] = vx[DIM*rowNumber + colNumber + 1] - vx[(rowNumber*DIM) +colNumber];}
		}
		}
		
		for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_QUAD_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;
		
		diverg = array[j][i] + array2[j][i];
		set_colormap(diverg, scalar_col,NCOLORS,0);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;

			
		
		diverg = array[j+1][i] + array2[j+1][i];
		set_colormap(diverg, scalar_col,NCOLORS,0);
			

			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			
			diverg = array[j][i+1] + array2[j][i+1];
			set_colormap(diverg, scalar_col,NCOLORS,0);
			
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		
		diverg = array[j+1][DIM-1] + array2[j+1][DIM-1];
		set_colormap(diverg, scalar_col,NCOLORS,0);
		glVertex2f(px, py);
		glEnd();
	}
		
	break;
	case 1: //F-vector field
	for(int rowNumber = 0; rowNumber <= DIM-1; ++rowNumber)
		{
			
		
		
		for(int colNumber = 0; colNumber <= DIM-1; ++colNumber)
		{	if(rowNumber == DIM-1){
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber)*DIM + colNumber];
		array2[rowNumber][colNumber] = vy[colNumber] -vy[(rowNumber)*DIM + colNumber];}
		else{
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber*DIM) +colNumber];
		array2[rowNumber][colNumber] = vy[DIM*(rowNumber+1) + colNumber] - vy[(rowNumber*DIM) +colNumber];}
		}
		}
		
		for(int rowNumber = 0; rowNumber <= DIM-1; ++rowNumber)
		{
			
		
		
		for(int colNumber = 0; colNumber <= DIM-1; ++colNumber)
		{	if(rowNumber == DIM-1){
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber)*DIM + colNumber];
		array[rowNumber][colNumber] = fx[DIM*rowNumber] -fx[(rowNumber)*DIM + colNumber];}
		else{
		//array[rowNumber][colNumber] = vx[DIM*(rowNumber+1) + colNumber] - vx[(rowNumber*DIM) +colNumber];
		array[rowNumber][colNumber] = fx[DIM*rowNumber + colNumber + 1] - fx[(rowNumber*DIM) +colNumber];}
		}
		}
		
		for (j = 0; j < DIM - 1; j++)			//draw smoke
	{
		glBegin(GL_QUAD_STRIP);

		i = 0;
		px = wn + (fftw_real)i * wn;
		py = hn + (fftw_real)j * hn;
		idx = (j * DIM) + i;
		
		diverg = array[j][i] + array2[j][i];
		set_colormap(diverg, scalar_col,NCOLORS,0);
		glVertex2f(px,py);

		for (i = 0; i < DIM - 1; i++)
		{
			px = wn + (fftw_real)i * wn;
			py = hn + (fftw_real)(j + 1) * hn;
			idx = ((j + 1) * DIM) + i;

			
		
		diverg = array[j+1][i] + array2[j+1][i];
		set_colormap(diverg, scalar_col,NCOLORS,0);
			

			glVertex2f(px, py);
			px = wn + (fftw_real)(i + 1) * wn;
			py = hn + (fftw_real)j * hn;
			idx = (j * DIM) + (i + 1);
			
			diverg = array[j][i+1] + array2[j][i+1];
			set_colormap(diverg, scalar_col,NCOLORS,0);
			
			glVertex2f(px, py);
		}

		px = wn + (fftw_real)(DIM - 1) * wn;
		py = hn + (fftw_real)(j + 1) * hn;
		idx = ((j + 1) * DIM) + (DIM - 1);
		
		diverg = array[j+1][DIM-1] + array2[j+1][DIM-1];
		set_colormap(diverg, scalar_col,NCOLORS,0);
		glVertex2f(px, py);
		glEnd();
	}
	
		break;}
}	
	
	if(isolineselected == 0){
		compute_codes();
		draw_isolines();
	}
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
	glVertex2f(470,0);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(470,100);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(500,100);
	glColor3f(0.4,0.2,1); //Purple
	glVertex2f(500,0);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(470,100);
	glColor3f(0.2,1,0);//Green
	glVertex2f(470,200);
	glColor3f(0.2,1,0);//Green
	glVertex2f(500,200);
	glColor3f(0,0.6,1); //Blue
	glVertex2f(500,100);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0.2,1,0);//Green
	glVertex2f(470,200);
	glColor3f(1,1,0); //Yellow
	glVertex2f(470,300);
	glColor3f(1,1,0); //Yellow
	glVertex2f(500,300);
	glColor3f(0.2,1,0);//Green
	glVertex2f(500,200);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(1,1,0); //Yellow
	glVertex2f(470,300);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(470,400);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(500,400);
	glColor3f(1,1,0); //Yellow
	glVertex2f(500,300);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(470,400);
	glColor3f(1,0,0); //Red
	glVertex2f(470,500);
	glColor3f(1,0,0); //Red
	glVertex2f(500,500);
	glColor3f(1,0.6,0); //Orange
	glVertex2f(500,400);
	glEnd();
}

void blue_to_yellow_bar(){
	
		// glBegin(GL_QUAD_STRIP);
	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0.0,1); 
	glVertex2f(470,0);
	glColor3f(0.2,1,0); 
	glVertex2f(470,250);
	glVertex2f(500,250);
	glColor3f(0,0.0,1); 
	glVertex2f(500,0);
	glEnd();

	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0.2,1,0); 
	glVertex2f(470,250);
	glColor3f(1,1,0); 
	glVertex2f(470,500);
	glVertex2f(500,500);
	glColor3f(0.2,1,0); 
	glVertex2f(500,250);
	glEnd();

}

void grayscale_bar(){
	
		// glBegin(GL_QUAD_STRIP);
	glBegin(GL_QUADS);
	glShadeModel(GL_SMOOTH);
	glColor3f(0,0,0); 
	glVertex2f(470,0);
	glColor3f(1,1,1); 
	glVertex2f(470,500);
	glVertex2f(500,500);
	glColor3f(0,0,0); 
	glVertex2f(500,0);
	glEnd();
}

void compute_codes(){
	for(int j = 0; j <= DIM-1; ++j)
		{for(int i = 0; i <= DIM-1; ++i)
		{//Get the coordinates of the cell vertices
		char string1,string2,string3,string4, tempstr;
		int idx1 = (j * DIM) + i;
		int idx2 = (j * DIM) + (i + 1);
		int idx3 = ((j + 1) * DIM) + i;
		int idx4 = ((j + 1) * DIM) + (DIM - 1);
		
		
		string1 = (rho[idx1] > isovalue) ? '1' : '0';
		string2 = (rho[idx2] > isovalue) ?  '1' : '0';
		//empstr = strcat(string1,string2);
		string3 = (rho[idx3] > isovalue) ?   '1' : '0';
		//tempstr = strcat(tempstr,string3);
		string4 = (rho[idx2] > isovalue) ?  '1' : '0';
		
		isocodes[j][i] = string1 + string2 + string3 + string4;
		
		}
		}
		
}



void draw_isolines( ){
	for(int j=0; j < DIM; ++j){
		for(int i=0; i<DIM; ++i){
			//Get cell vectors
		int idx1 = (j * DIM) + i;
		int idx2 = (j * DIM) + (i + 1);
		int idx3 = ((j + 1) * DIM) + (i+1); //DIM -1
		int idx4 = ((j + 1) * DIM) + i;
		 
		 //Get Cell edges
		
		float rat, rat2, rat3, rat4;	
		
		float ed1_x, ed2_x;
		float ed1_y, ed2_y;
		float ed3_x, ed3_y;
		float ed4_x, ed4_y;

			switch(isocodes[i][j]){
				case '0001':
				case '1110':
				rat = ((rho[idx4] - rho[idx1]) != 0) ? ((idx4*(rho[idx4] -isovalue) + idx1*(isovalue-rho[idx1])) /(rho[idx4] - rho[idx1])) : 0; 
				rat2 = ((rho[idx4] - rho[idx3]) != 0) ? ((idx4*(rho[idx4] -isovalue) + idx3*(isovalue-rho[idx3])) /(rho[idx4] - rho[idx3])) : 0; 
				glBegin(GL_LINES);
				ed1_x = ((j+1)*DIM);
				ed1_y = (i+1) + ((i+1) - i)*rat;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = (j+1)*DIM + (i - (i+1))*rat2;
				ed2_y = i;
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '0010':
				case '1101':
				rat = ((rho[idx3] - rho[idx4]) != 0) ? (idx3*(rho[idx3] -isovalue) + idx4*(isovalue-rho[idx4])) /(rho[idx3] - rho[idx4]) : 0; 
				rat2 = ((rho[idx3] - rho[idx2]) != 0) ? (idx3*(rho[idx3] -isovalue) + idx2*(isovalue-rho[idx2]) /(rho[idx3] - rho[idx4])) : 0; 
				glBegin(GL_LINES);
				ed1_x = ((j+1)*DIM);
				ed1_y = (i+1) + rat;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = (j+1)*DIM + rat2;
				ed2_y = i;
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '0011':
				case '1100':
				rat = ((rho[idx4] - rho[idx1]) != 0) ? (idx4*(isovalue-rho[idx4]) + idx1*(rho[idx1] - isovalue))/(rho[idx4] - rho[idx1]) : 0; 
				rat2 = ((rho[idx3] - rho[idx2])  != 0) ? (idx3*(isovalue-rho[idx3]) + idx2*(rho[idx2] -isovalue))/(rho[idx3] - rho[idx2]) : 0; 
				glBegin(GL_LINES);
				
				ed1_x = (j+1)*DIM - rat;
				ed1_y = i;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = (j+1)*DIM - rat2;
				ed2_y = (i+1);
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '0100':
				case '1011':
				rat = ((rho[idx2] - rho[idx1]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx2] - rho[idx1])) : 0; 
				rat2 = ((rho[idx2] - rho[idx3]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx3*(isovalue-rho[idx3]) /(rho[idx2] - rho[idx3])) : 0; 
				glBegin(GL_LINES);
				ed1_x = j*DIM;
				ed1_y = (i+1) - rat;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = (j+1)*DIM - rat2;
				ed2_y = i;
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '0101':
				
				rat = ((rho[idx4] - rho[idx1]) != 0) ? (idx4*(rho[idx4] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx4] - rho[idx1])) : 0; 
				rat2 = ((rho[idx4] - rho[idx3]) != 0) ? (idx4*(rho[idx4] - isovalue) + idx3*(isovalue - rho[idx3]) / (rho[idx4]- rho[idx3])) : 0;
				rat3 = ((rho[idx2] - rho[idx3]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx3*(isovalue-rho[idx3]) /(rho[idx2] - rho[idx3])) : 0; 
				rat4 = ((rho[idx2] - rho[idx1]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx2] - rho[idx1])) : 0;
				glBegin(GL_LINES);
				ed1_x = (j+1)*DIM - rat;
				ed1_y = i;
				glVertex2f(ed1_x, ed1_y);
				
				ed2_x = (j+1)*DIM;
				ed2_y = i + rat2;
				glVertex2f(ed2_x, ed2_y);
				
				ed3_x = j*DIM + rat3;
				ed3_y = i+1;
				glVertex2f(ed3_x, ed3_y);
				
				ed4_x = j*DIM;
				ed4_y = (i+1) - rat4;
				glVertex2f(ed4_x,ed4_y);
				glEnd();
				break;
				
				case '0110':
				case '1001':
				rat = ((rho[idx2] - rho[idx1]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx2] - rho[idx1])) : 0; 
				rat2 = ((rho[idx3] - rho[idx4]) != 0) ? (idx3*(rho[idx3] -isovalue) + idx4*(isovalue-rho[idx4]) /(rho[idx3] - rho[idx4])) : 0; 
				glBegin(GL_LINES);
				ed1_x = j*DIM;
				ed1_y = (i+1) - rat;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = (j+1)*DIM;
				ed2_y = i + rat2;
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '0111':
				case '1000':
				rat = ((rho[idx2] - rho[idx1]) != 0) ? (idx2*(rho[idx2] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx2] - rho[idx1])) : 0; 
				rat2 = ((rho[idx4] - rho[idx1]) != 0) ? (idx4*(rho[idx4] -isovalue) + idx1*(isovalue-rho[idx1]) /(rho[idx4] - rho[idx1])) : 0; 
				glBegin(GL_LINES);
				ed1_x = j*DIM;
				ed1_y = (i+1) - rat;
				glVertex2f(ed1_x,ed1_y);
				
				ed2_x = j*DIM + rat2;
				ed2_y = i;
				glVertex2f(ed2_x,ed2_y);
				glEnd();
				break;
				
				case '1010':
				rat = ((rho[idx1] - rho[idx2]) != 0) ? (idx1*(rho[idx1] -isovalue) + idx2*(isovalue-rho[idx2]) /(rho[idx1] - rho[idx2])) : 0; 
				rat2 = ((rho[idx1] - rho[idx4]) != 0) ? (idx1*(rho[idx1] - isovalue) + idx4*(isovalue - rho[idx4]) / (rho[idx1]- rho[idx4])) : 0;
				rat3 = ((rho[idx3] - rho[idx2]) != 0) ? (idx3*(rho[idx3] -isovalue) + idx2*(isovalue-rho[idx2]) /(rho[idx3] - rho[idx2])) : 0; 
				rat4 = ((rho[idx3] - rho[idx4]) != 0) ? (idx3*(rho[idx3] -isovalue) + idx4*(isovalue-rho[idx4]) /(rho[idx3] - rho[idx4])) : 0;
				glBegin(GL_LINES);
				ed1_x = (j+1)*DIM;
				ed1_y = i + rat;
				glVertex2f(ed1_x, ed1_y);
				
				ed2_x = (j+1)*DIM - rat2;
				ed2_y = i;
				glVertex2f(ed2_x, ed2_y);
				
				ed3_x = j*DIM + rat3;
				ed3_y = i+1;
				glVertex2f(ed3_x, ed3_y);
				
				ed4_x = (j+1)*DIM;
				ed4_y = (i+1) - rat4;
				glVertex2f(ed4_x,ed4_y);
				glEnd();
				break;
				
				
				}
			
		}
	}
}
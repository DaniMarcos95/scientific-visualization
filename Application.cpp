#include <rfftw.h>              //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <math.h>               //for various math functions
#include <GL/glut.h>            //the GLUT graphics library
#include "Visualization.h"
#include "Simulation.h"
#include <GL/glui.h>


using namespace std;

extern int   winWidth, winHeight;      
extern int   color_dir;           
extern float vec_scale;			
extern int   draw_smoke;           
extern int   draw_vecs;   
extern int NCOLORS;         
// const int COLOR_BLACKWHITE = 0;   
const int COLOR_RAINBOW = 1;
const int COLOR_GRAYSCALE = 0;
const int COLOR_BLUEYEL = 2;
extern int   scalar_col;           
extern int   frozen;
extern double dt;	
extern float visc;               
const int DIM = 50;
extern fftw_real *fx, *fy; 
extern fftw_real *rho;
GLUI *glui;
GLUI_Spinner *spinnerNumberColors;
int main_window;
// int  densityCheck  = 1;
// int velocityCheck = 0;
int colorMapIndex;
int datasetIndex;
int numberOfColors = 255;

//display: Handle window redrawing events. Simply delegates to visualize().
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	visualize();
	glFlush();
	glutSwapBuffers();

	// if ( densityCheck ){
 //    	draw_vecs = 0;
	// 	draw_smoke = 1;
	// 	velocityCheck = 0;
 //    }
  	
 //  	if(velocityCheck){
 //  		draw_vecs = 1;
	// 	draw_smoke = 0;
	// 	densityCheck = 0;
 //  	}

	switch(datasetIndex){
		case 0: draw_vecs = 0;
				draw_smoke = 1;
				break;
		case 1: draw_vecs = 1;
				draw_smoke = 0;
				break;
	}

	switch(colorMapIndex){
		case 0: scalar_col = COLOR_GRAYSCALE;
				color_dir = COLOR_GRAYSCALE;
				break;
		case 1: scalar_col = COLOR_RAINBOW;
				color_dir = COLOR_RAINBOW;
				break;
		case 2: scalar_col = COLOR_BLUEYEL;
				color_dir = COLOR_BLUEYEL;
				break;
	}

	NCOLORS = numberOfColors;
}

//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h)
{
 	glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
	winWidth = w; winHeight = h;
}

//keyboard: Handle key presses
void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	  case 't': dt -= 0.001; break;
	  case 'T': dt += 0.001; break;
	  case 'c': color_dir++; if (color_dir>3) color_dir = 0; break;
	  case 'S': vec_scale *= 1.2; break;
	  case 's': vec_scale *= 0.8; break;
	  case 'V': visc *= 5; break;
	  case 'v': visc *= 0.2; break;
	  case 'x': draw_smoke = 1 - draw_smoke;
		    if (draw_smoke==0) draw_vecs = 1; 
		    break;
	  case 'y': draw_vecs = 1 - draw_vecs;
		    if (draw_vecs==0) draw_smoke = 1; 
		    break;
	  case 'm': scalar_col++; color_dir++; if (scalar_col>COLOR_BLUEYEL && color_dir>3) {scalar_col=COLOR_GRAYSCALE; color_dir =0;} break;
	  case 'a': frozen = 1-frozen; break;
	  case 'q': exit(0);
	}
}

// drag: When the user drags with the mouse, add a force that corresponds to the direction of the mouse
//       cursor movement. Also inject some new matter into the field at the mouse location.
void drag(int mx, int my)
{
	int xi,yi,X,Y; double  dx, dy, len;
	static int lmx=0,lmy=0;				//remembers last mouse location

	// Compute the array index that corresponds to the cursor location
	xi = (int)clamp((double)(DIM + 1) * ((double)mx / (double)winWidth));
	yi = (int)clamp((double)(DIM + 1) * ((double)(winHeight - my) / (double)winHeight));

	X = xi; Y = yi;

	if (X > (DIM - 1))  X = DIM - 1; 
	if (Y > (DIM - 1))  Y = DIM - 1;
	if (X < 0) X = 0; 
	if (Y < 0) Y = 0;

	// Add force at the cursor location
	my = winHeight - my;
	dx = mx - lmx; dy = my - lmy;
	len = sqrt(dx * dx + dy * dy);
	if (len != 0.0) {  dx *= 0.1 / len; dy *= 0.1 / len; }
	fx[Y * DIM + X] += dx;
	fy[Y * DIM + X] += dy;
	rho[Y * DIM + X] = 10.0f;
	lmx = mx; lmy = my;
}

void myGlutIdle( void ) {
	glutSetWindow(main_window);
	glutPostRedisplay();         
}

void shutDown(){
	exit(0);
}

int main(int argc, char **argv, int NLEVELS)
{
	printf("Fluid Flow Simulation and Visualization\n");
	printf("=======================================\n");
	printf("Click and drag the mouse to steer the flow!\n");
	printf("T/t:   increase/decrease simulation timestep\n");
	printf("S/s:   increase/decrease hedgehog scaling\n");
	printf("c:     toggle direction coloring on/off\n");
	printf("V/v:   increase decrease fluid viscosity\n");
	printf("x:     toggle drawing matter on/off\n");
	printf("y:     toggle drawing hedgehogs on/off\n");
	printf("m:     toggle thru scalar coloring\n");
	printf("a:     toggle the animation on/off\n");
	printf("q:     quit\n\n");
	//printf("Input number of colors (1-255)\n";
	//cin >> NLEVELS;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(500,500);
	main_window = glutCreateWindow("Real-time smoke simulation and visualization");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(do_one_simulation_step);
	glutKeyboardFunc(keyboard);
	glutMotionFunc(drag);

	
	GLUI *glui = GLUI_Master.create_glui("");
	// GLUI_Master.set_glutIdleFunc(myGlutIdle);
	// glui->set_main_gfx_window( main_window );

	GLUI_Panel *mainPanel = new GLUI_Panel(glui, "Configuration");
	mainPanel->set_alignment(GLUI_ALIGN_LEFT);

	// GLUI_Listbox *listboxDataset = new GLUI_Listbox(objectEdit, "Dataset: ", &datasetIndex,  12);
	// listboxDataset->set_alignment(GLUI_ALIGN_LEFT);
	// listboxDataset->add_item(0, "Density");
	// listboxDataset->add_item(1, "Velocity");
	GLUI_Panel *datasetPanel = new GLUI_Panel(mainPanel, "Dataset");

	GLUI_RadioGroup *radioDataset = new GLUI_RadioGroup(datasetPanel, &datasetIndex);
	GLUI_RadioButton *buttonDensity = new GLUI_RadioButton( radioDataset, "Density" );
	GLUI_RadioButton *buttonVelocity = new GLUI_RadioButton( radioDataset, "Velocity" );


	glui->add_separator_to_panel( mainPanel );

	// GLUI_Listbox *listboxColorMap = new GLUI_Listbox(objectEdit, "Colormap: ", &colorMapIndex,  12);
	// listboxColorMap->set_alignment(GLUI_ALIGN_LEFT);
	// listboxColorMap->add_item(0, "Grayscale");
	// listboxColorMap->add_item(1, "Rainbow");
	// listboxColorMap->add_item(2, "Blue to yellow");
	GLUI_Panel *colorMapPanel = new GLUI_Panel(mainPanel, "Colormap");

	GLUI_RadioGroup *radioColorMap = new GLUI_RadioGroup(colorMapPanel, &colorMapIndex);
	GLUI_RadioButton *buttonGrayScale = new GLUI_RadioButton( radioColorMap, "GrayScale" );
	GLUI_RadioButton *buttonRainbow = new GLUI_RadioButton( radioColorMap, "Rainbow" );
	GLUI_RadioButton *buttonBlueToYellow = new GLUI_RadioButton( radioColorMap, "Blue To Yellow" );

	glui->add_separator_to_panel( mainPanel );

	GLUI_EditText *setNumberOfColors = new GLUI_EditText(mainPanel, "Number of colors(2-255):", &numberOfColors);
	setNumberOfColors->set_int_val(255);
	setNumberOfColors->set_int_limits( 2, 255, GLUI_LIMIT_CLAMP );

	 GLUI_Spinner *spinnerNumberColors = new GLUI_Spinner( mainPanel, "Number of colors:", GLUI_EDITTEXT_INT, &numberOfColors );
	 spinnerNumberColors->set_int_limits(2, 255);
	 spinnerNumberColors->set_int_val(255);
	 NCOLORS = spinnerNumberColors->get_int_val();
	
	init_simulation(DIM);	//initialize the simulation data structures
	glutMainLoop();			//calls do_one_simulation_step, keyboard, display, drag, reshape
	return 0;
}
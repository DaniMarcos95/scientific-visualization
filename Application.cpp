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
extern int draw_vec_mod;   
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
int colorMapIndex;
int datasetIndex;
int numberOfColors = 255;
float max_clamped = 10;
float min_clamped = 0;
extern float max_rho;
extern float min_rho;

//display: Handle window redrawing events. Simply delegates to visualize().
void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	visualize();
	glFlush();
	glutSwapBuffers();

	switch(datasetIndex){
		case 0: draw_vecs = 0;
				draw_smoke = 1;
				draw_vec_mod = 0;
				break;
		case 1: draw_vecs = 0;
				draw_smoke = 0;
				draw_vec_mod = 1;
				break;
		case 2: draw_vecs = 1;
				draw_smoke = 1;
				draw_vec_mod = 0;
				break;
		/*case 3: draw_vecs = 1;
				draw_smoke = 1;
				draw_vec_mod = 0;
				break;*/
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
}

void shutDown(){
	exit(0);
}

int main(int argc, char **argv, int NLEVELS)
{
	/*printf("Fluid Flow Simulation and Visualization\n");
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
	printf("q:     quit\n\n");*/
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

	GLUI_Panel *mainPanel = new GLUI_Panel(glui, "Configuration");
	mainPanel->set_alignment(GLUI_ALIGN_LEFT);

	GLUI_Panel *panel1 = new GLUI_Panel(mainPanel, "");

	GLUI_Panel *datasetPanel = new GLUI_Panel(panel1, "Dataset");

	GLUI_RadioGroup *radioDataset = new GLUI_RadioGroup(datasetPanel, &datasetIndex);
	GLUI_RadioButton *buttonDensity = new GLUI_RadioButton( radioDataset, "rho" );
	GLUI_RadioButton *buttonVecMod = new GLUI_RadioButton( radioDataset, "||v||" );
	GLUI_RadioButton *buttonVelocity = new GLUI_RadioButton( radioDataset, "v" );

	glui->add_column_to_panel(panel1, true);

	GLUI_Panel *colorMapPanel = new GLUI_Panel(panel1, "Colormap");

	GLUI_RadioGroup *radioColorMap = new GLUI_RadioGroup(colorMapPanel, &colorMapIndex);
	GLUI_RadioButton *buttonGrayScale = new GLUI_RadioButton( radioColorMap, "GrayScale" );
	GLUI_RadioButton *buttonRainbow = new GLUI_RadioButton( radioColorMap, "Rainbow" );
	GLUI_RadioButton *buttonBlueToYellow = new GLUI_RadioButton( radioColorMap, "Blue To Yellow" );

	GLUI_Panel *colorPanel = new GLUI_Panel(mainPanel, "Color configuration");

	GLUI_EditText *setNumberOfColors = new GLUI_EditText(colorPanel, "Number of colors(2-255):", &NCOLORS);
	setNumberOfColors->set_int_val(255);
	setNumberOfColors->set_int_limits( 2, 255, GLUI_LIMIT_CLAMP );

	GLUI_Panel *scalingPanel = new GLUI_Panel(mainPanel, "Scaling/Clamping");

	GLUI_EditText *maxClamped = new GLUI_EditText(scalingPanel, "Maximum value for clamping:", &max_clamped);

	GLUI_EditText *minClamped = new GLUI_EditText(scalingPanel, "Minimum value for clamping:", &min_clamped);

	glui->set_main_gfx_window(main_window);
	init_simulation(DIM);	//initialize the simulation data structures
	glutMainLoop();			//calls do_one_simulation_step, keyboard, display, drag, reshape
	return 0;
}
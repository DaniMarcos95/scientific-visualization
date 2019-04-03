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
extern int   draw_rho;           
extern int   draw_vecs;
extern int draw_for;
extern int draw_vec_mod; 
extern int draw_for_mod;   
extern int NCOLORS; 
extern int diver;        
const int COLOR_BLACKWHITE = 0;   
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
int colorMapIndex = 0;
int datasetIndex = 0;
int scalarIndex = 0;
int scalarVectorIndex = 0;
int optionIndex = 0;
int vectorIndex = 0;
int glyphIndex = 0;
int scalingClampingIndex = 0;
int divergenceIndex = 0;
int isolineselected = 0;
int numberOfColors = 255;
float max_clamped = 0;
float min_clamped = 10;
extern float max_rho;
extern float min_rho;
extern float max_v;
extern float min_v;
extern float max_f;
extern float min_f;
extern fftw_real  hn;
GLUI_EditText *maxClamped;
GLUI_EditText *minClamped;
int aux_repetition_scalar = -1;
int aux_repetition_vector = -1;
int numberOfSamples = 50;
float isovalue = 0;

//display: Handle window redrawing events. Simply delegates to visualize().

void drawBitmapText(char *string,float x,float y) 
{  
	char *c;
	glRasterPos2f(x, y);
	int count = 0;

	for (c=string; *c != ' '; c++) 
	{	
		if(count == 5) break;
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
		count++;
	}
}

void render_text(void)
{ 	
	glColor3f(1,1,1);
	

	char array[10];
	snprintf(array, sizeof(array), "%f ", max_rho);
	drawBitmapText(array,winWidth - 70,winHeight - 20);

	snprintf(array, sizeof(array), "%f ", min_rho);
	drawBitmapText(array,winWidth - 70,15);

	snprintf(array, sizeof(array), "%f ", (min_rho + max_rho)/2);
	drawBitmapText(array,winWidth - 70,winHeight/2);

	snprintf(array, sizeof(array), "%f ",min_rho + 3*(fabs(min_rho - max_rho))/4);
	drawBitmapText(array,winWidth - 70,3*winHeight/4);

	snprintf(array, sizeof(array), "%f ", min_rho + (fabs(min_rho - max_rho))/4);
	drawBitmapText(array,winWidth - 70,winHeight/4);
}

void display(void)
{

	float stepSize = (float)(winHeight - 2 *hn) / ((float) NCOLORS+1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	// glLoadIdentity();
	visualize();

	int bar_width = 40;
	int gridWidth = winWidth - bar_width;
	for (int i = 0; i < NCOLORS+1; ++i)
	{
		float value = (float) i/((float) (NCOLORS+1));
		float y1 = hn+stepSize*i;
        float x1 = gridWidth; // 
        float y2 = hn+stepSize*(i+1);
        float x2 = gridWidth + bar_width;
        set_colormap(value, scalar_col, NCOLORS+1,1);
        glRecti(x1,y1,x2,y2);
	}

	render_text();

	glFlush();
	glutSwapBuffers();

	int repetition_scalar = -1;
	int repetition_vector = -1;

	switch(optionIndex){
		case 0: //scalars
			switch(scalarIndex){
				case 0: repetition_scalar = 0;
						draw_rho = 1;
						draw_vec_mod = 0;
						draw_for_mod = 0;
						if(scalingClampingIndex == 0){
								maxClamped->set_float_val(0);
								minClamped->set_float_val(0);
						}else{
							if(repetition_scalar != aux_repetition_scalar){
								if(max_rho < min_rho){
									max_rho = 10;
									min_rho = 0;
								}
								maxClamped->set_float_val(10);
								minClamped->set_float_val(0);
								aux_repetition_scalar = repetition_scalar;
							}
						}
						
						break;
				case 1: repetition_scalar = 1;
						draw_rho = 0;
						draw_vec_mod = 1;
						draw_for_mod = 0;
						if(scalingClampingIndex == 0){
								maxClamped->set_float_val(0);
								minClamped->set_float_val(0);
						}else{
							if(repetition_scalar != aux_repetition_scalar){
								maxClamped->set_float_val(0.02);
								minClamped->set_float_val(0);
								aux_repetition_scalar = repetition_scalar;
							}
						}
						break;
				case 2:	repetition_scalar = 2;
						draw_rho = 0;
						draw_vec_mod = 0;
						draw_for_mod = 1;
						if(scalingClampingIndex == 0){
								maxClamped->set_float_val(0);
								minClamped->set_float_val(0);
						}else{
							if(repetition_scalar != aux_repetition_scalar){
								maxClamped->set_float_val(0.2);
								minClamped->set_float_val(0);
								aux_repetition_scalar = repetition_scalar;
							}
						}
						break;
			}
			break;
		case 1: 
			switch(scalarVectorIndex){
				case 0: repetition_scalar = 0;
						draw_rho = 1;
						draw_vec_mod = 0;
						draw_for_mod = 0;
						if(scalingClampingIndex == 0){
									maxClamped->set_float_val(0);
									minClamped->set_float_val(0);
							}else{
								if(repetition_scalar != aux_repetition_scalar){
									if(max_rho < min_rho){
										max_rho = 10;
										min_rho = 0;
									}
									maxClamped->set_float_val(10);
									minClamped->set_float_val(0);
									aux_repetition_scalar = repetition_scalar;
								}
							}
						break;
				case 1: repetition_scalar = 1;
						draw_rho = 0;
						draw_vec_mod = 1;
						draw_for_mod = 0;
						if(scalingClampingIndex == 0){
									maxClamped->set_float_val(0);
									minClamped->set_float_val(0);
							}else{
								if(repetition_scalar != aux_repetition_scalar){
									maxClamped->set_float_val(0.02);
									minClamped->set_float_val(0);
									aux_repetition_scalar = repetition_scalar;
								}
							}
						break;
				case 2:	repetition_scalar = 2;
						draw_rho = 0;
						draw_vec_mod = 0;
						draw_for_mod = 1;
						if(scalingClampingIndex == 0){
									maxClamped->set_float_val(0);
									minClamped->set_float_val(0);
							}else{
								if(repetition_scalar != aux_repetition_scalar){
									maxClamped->set_float_val(0.2);
									minClamped->set_float_val(0);
									aux_repetition_scalar = repetition_scalar;
								}
							}
						break;
			}

			switch(vectorIndex){
				case 0: repetition_vector = 0;
						draw_vecs = 1;
						draw_for = 0;
						break;
				case 1: repetition_vector = 1;
						draw_vecs = 0;
						draw_for = 1;
						break;
			}

			break;

		case 3: 
			switch(divergenceIndex){
				case 0:
					repetition_vector = 0;
					diver = 1;
					draw_vecs = 1;
					draw_for = 0;
					// if(repetition_vector != aux_repetition_vector){
					// 		 	maxClamped->set_float_val(max_v);
					// 		 	minClamped->set_float_val((-1)*max_v);
					// 		 	aux_repetition_vector = repetition_vector;
					// }
					break;
				case 1:
					repetition_vector = 1;
					diver = 1;
					draw_vecs = 0;
					draw_for = 1;
					// if(repetition_vector != aux_repetition_vector){
					// 		 	maxClamped->set_float_val(max_f);
					// 		 	minClamped->set_float_val((-1)*max_f);
					// 		 	aux_repetition_vector = repetition_vector;
					// }
					break;
			}
		case 4: 
				break;//Height
	}

	// switch(scalarIndex){
	// 	case 0: repetition_scalar = 0;
		
			
	// 			draw_rho = 1;
	// 			draw_vec_mod = 0;
	// 		draw_for_mod = 0;
	// 		diver = 0;
			
	// 			if(repetition_scalar != aux_repetition_scalar){
	// 				if(max_rho < min_rho){
	// 					max_rho = 10;
	// 					min_rho = 0;
	// 				}
	// 				maxClamped->set_float_val(max_rho);
	// 				minClamped->set_float_val(min_rho);
	// 				aux_repetition_scalar = repetition_scalar;
	// 			}
	// 			break;
	// 	case 1: repetition_scalar = 1;
	// 			draw_rho = 0;
	// 			draw_vec_mod = 1;
	// 			draw_for_mod = 0;
	// 			diver = 0;
	// 			if(repetition_scalar != aux_repetition_scalar){
	// 				maxClamped->set_float_val(max_v);
	// 				minClamped->set_float_val(min_v);
	// 				aux_repetition_scalar = repetition_scalar;
	// 			}
	// 			break;
	// 	case 2:	repetition_scalar = 2;
	// 			draw_rho = 0;
	// 			draw_vec_mod = 0;
	// 			draw_for_mod = 1;
	// 			diver = 0;
	// 			if(repetition_scalar != aux_repetition_scalar){
	// 				maxClamped->set_float_val(max_f);
	// 				minClamped->set_float_val(min_f);
	// 				aux_repetition_scalar = repetition_scalar;
	// 			}
	// 			break;
	// 	case 3:	repetition_scalar = 3;
	// 			draw_rho = 0;
	// 			draw_vec_mod = 0;
	// 			draw_for_mod = 0;
	// 			diver = 0;
	// 			if(repetition_scalar != aux_repetition_scalar){
	// 				maxClamped->set_float_val(max_f);
	// 				minClamped->set_float_val(min_f);
	// 				aux_repetition_scalar = repetition_scalar;
	// 			}
	// 			break;
	// 	case 4: repetition_scalar = 4;
	// 			draw_rho = 0;
	// 			draw_vec_mod = 0;
	// 			draw_for_mod = 0;
	// 			diver = 1;
	// 			if(repetition_scalar != aux_repetition_scalar){
	// 						maxClamped->set_float_val(max_f);
	// 						minClamped->set_float_val(min_f);
	// 						aux_repetition_scalar = repetition_scalar;
	// 			}
	// 			break;
	// }

	// switch(vectorIndex){
	// 	case 0: repetition_vector = 0;
	// 			draw_vecs = 1;
	// 			draw_for = 0;
	// 			 if(repetition_vector != aux_repetition_vector){
	// 			 	maxClamped->set_float_val(max_rho);
	// 			 	minClamped->set_float_val(min_rho);
	// 			 	aux_repetition_vector = repetition_vector;
	// 			 }	
	// 			break;
	// 	case 1: repetition_vector = 1;
	// 			draw_vecs = 0;
	// 			draw_for = 1;
	// 			 if(repetition_vector != aux_repetition_vector){
	// 			 	maxClamped->set_float_val(max_rho);
	// 			 	minClamped->set_float_val(min_rho);
	// 			 	aux_repetition_vector = repetition_vector;
	// 			 }
	// 			break;
	// 	case 2: repetition_vector = 2;
	// 			draw_vecs = 0;
	// 			draw_for = 0;
	// 			 if(repetition_vector != aux_repetition_vector){
	// 			 	maxClamped->set_float_val(max_rho);
	// 			 	minClamped->set_float_val(min_rho);
	// 			 	aux_repetition_vector = repetition_vector;
	// 			 }
	// 			break;
	// }

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

void scaleGlyphUp(){
	vec_scale *= 1.2;
}

void scaleGlyphDown(){
	vec_scale *= 0.8;
}
//reshape: Handle window resizing (reshaping) events
void reshape(int w, int h)
{
 	glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
	glMatrixMode(GL_PROJECTION);
	glClearDepth(1.0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glLoadIdentity();
	glOrtho(0.0, (GLdouble)w, 0.0, (GLdouble)h, -1000, 1000);
	glMatrixMode(GL_MODELVIEW);
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
	  case 'x': draw_rho = 1 - draw_rho;
		    if (draw_rho==0) draw_vecs = 1; 
		    break;
	  case 'y': draw_vecs = 1 - draw_vecs;
		    if (draw_vecs==0) draw_rho = 1; 
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

	GLUI_Panel *panel1 = new GLUI_Panel(mainPanel, "General Options");

	GLUI_Panel *optionPanel = new GLUI_Panel(panel1, "");

	GLUI_RadioGroup *radioOptions = new GLUI_RadioGroup(optionPanel, &optionIndex);
	GLUI_RadioButton *optionDensity = new GLUI_RadioButton( radioOptions, "Scalars" );
	GLUI_RadioButton *optionVec = new GLUI_RadioButton( radioOptions, "Vectors" );
	GLUI_RadioButton *optionDIv = new GLUI_RadioButton( radioOptions, "Divergence" );
	GLUI_RadioButton *optionHeight = new GLUI_RadioButton( radioOptions, "Height Plot" );
	

	GLUI_Panel *colorPanel = new GLUI_Panel(panel1, "Color configuration");

	GLUI_EditText *setNumberOfColors = new GLUI_EditText(colorPanel, "Number of colors(2-255):", &NCOLORS);
	setNumberOfColors->set_int_val(255);
	setNumberOfColors->set_int_limits( 2, 255, GLUI_LIMIT_CLAMP );

	GLUI_Panel *scalarPanel = new GLUI_Panel(mainPanel, "Scalars");

	GLUI_RadioGroup *radioScalar = new GLUI_RadioGroup(scalarPanel, &scalarIndex);
	GLUI_RadioButton *buttonDensity = new GLUI_RadioButton( radioScalar, "rho" );
	GLUI_RadioButton *buttonVecMod = new GLUI_RadioButton( radioScalar, "||v||" );
	GLUI_RadioButton *buttonForMod = new GLUI_RadioButton( radioScalar, "||f||" );

	GLUI_Panel *panelVectors = new GLUI_Panel(mainPanel, "Vectors");

	GLUI_Panel *vectorsScalar = new GLUI_Panel(panelVectors, "Scalar");

	GLUI_RadioGroup *vectorScalar = new GLUI_RadioGroup(vectorsScalar, &scalarVectorIndex);
	GLUI_RadioButton *vectorDensity = new GLUI_RadioButton( vectorScalar, "rho" );
	GLUI_RadioButton *vectorVecMod = new GLUI_RadioButton( vectorScalar, "||v||" );
	GLUI_RadioButton *vectorForMod = new GLUI_RadioButton( vectorScalar, "||f||" );
	

	glui->add_column_to_panel(panelVectors, true);

	GLUI_Panel *vectorsFields = new GLUI_Panel(panelVectors, "Vector Fields");

	GLUI_RadioGroup *radioVectors = new GLUI_RadioGroup(vectorsFields, &vectorIndex);
	GLUI_RadioButton *buttonVelocity = new GLUI_RadioButton( radioVectors, "v" );
	GLUI_RadioButton *buttonForces = new GLUI_RadioButton( radioVectors, "f" );
	
	glui->add_column_to_panel(panelVectors, true);

	GLUI_Panel *glyphPanel = new GLUI_Panel(panelVectors, "Glyphs");

	GLUI_RadioGroup *glyphMode = new GLUI_RadioGroup(glyphPanel, &glyphIndex);
	GLUI_RadioButton *lines = new GLUI_RadioButton( glyphMode, "Lines" );
	GLUI_RadioButton *cones = new GLUI_RadioButton( glyphMode, "Cones" );
	GLUI_RadioButton *arrows = new GLUI_RadioButton( glyphMode, "Arrows" );

	glui->add_column_to_panel(panel1, true);

	GLUI_Panel *colorMapPanel = new GLUI_Panel(panel1, "Colormap");

	GLUI_RadioGroup *radioColorMap = new GLUI_RadioGroup(colorMapPanel, &colorMapIndex);
	GLUI_RadioButton *buttonGrayScale = new GLUI_RadioButton( radioColorMap, "GrayScale" );
	GLUI_RadioButton *buttonRainbow = new GLUI_RadioButton( radioColorMap, "Rainbow" );
	GLUI_RadioButton *buttonBlueToYellow = new GLUI_RadioButton( radioColorMap, "Blue To Yellow" );

	GLUI_Panel *scalingPanel = new GLUI_Panel(panel1, "Scaling/Clamping");

	GLUI_RadioGroup *radioScaling = new GLUI_RadioGroup(scalingPanel, &scalingClampingIndex);
	GLUI_RadioButton *buttonScaling = new GLUI_RadioButton( radioScaling, "Scaling" );
	GLUI_RadioButton *buttonClamping = new GLUI_RadioButton( radioScaling, "Clamping" );


	maxClamped = new GLUI_EditText(scalingPanel, "Maximum value for clamping:", &max_clamped);

	minClamped = new GLUI_EditText(scalingPanel, "Minimum value for clamping:", &min_clamped);

	glui->add_column_to_panel(glyphPanel, true);

	GLUI_Button *glyphScaleUp = new GLUI_Button(glyphPanel, "Scale Up", -1, (GLUI_Update_CB)scaleGlyphUp);
	GLUI_Button *glyphScaleDown = new GLUI_Button(glyphPanel, "Scale Down", -1, (GLUI_Update_CB)scaleGlyphDown);


	GLUI_EditText *setNumberOfSamples = new GLUI_EditText(panel1, "Number of Samples (Dimension of Square) ", &numberOfSamples);
	setNumberOfSamples->set_int_val(50);
	setNumberOfSamples->set_int_limits( 1, 50, GLUI_LIMIT_CLAMP );

	GLUI_Panel *diverPanel = new GLUI_Panel(mainPanel, "Divergence");

	GLUI_RadioGroup *radioDivergence = new GLUI_RadioGroup(diverPanel, &divergenceIndex);
	GLUI_RadioButton *vDivergence = new GLUI_RadioButton( radioDivergence, "v" );
	GLUI_RadioButton *fDivergence = new GLUI_RadioButton( radioDivergence, "f" );

	// GLUI_EditText *isovalueText = new GLUI_EditText(mainPanel, "Number of Samples (Dimension of Square) ", &isovalue);
	// setNumberOfSamples->set_int_val(0.0);

	// GLUI_RadioGroup *radioIsolines = new GLUI_RadioGroup(mainPanel, &isolineselected);
	// GLUI_RadioButton *buttonIsolines = new GLUI_RadioButton( radioIsolines, "Isolines" );

	GLUI_Button *Exit = new GLUI_Button(mainPanel, "Exit", -1, (GLUI_Update_CB)exit);

	glui->set_main_gfx_window(main_window);
	init_simulation(DIM);	//initialize the simulation data structures
	glutMainLoop();			//calls do_one_simulation_step, keyboard, display, drag, reshape
	return 0;
}



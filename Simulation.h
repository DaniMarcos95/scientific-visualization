#ifndef SIMULATION_H
#define SIMULATION_H

class Simulation{

	const int DIM = 50;				//size of simulation grid
	double dt;				//simulation time step
	float visc;				//fluid viscosity
	fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
	fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
	fftw_real *fx, *fy;	            //(fx,fy)   = user-controlled simulation forces, steered with the mouse
	fftw_real *rho, *rho0;			//smoke density at the current (rho) and previous (rho0) moment
	rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization

	public:
		Simulation(int);
		void FFT(int, void);
		int clamp(float);
		float max(float, float);
		void solve(int, fftw_real*, fftw_real*, fftw_real*, fftw_real*, fftw_real, fftw_real);
		void diffuse_matter(int, fftw_real, fftw_real, fftw_real, fftw_real, fftw_real dt);
		void set_forces(void);
		void do_one_simulation_step(void);

};
#endif
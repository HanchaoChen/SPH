#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>

const double g = 9.8;
const double rho_0 = 1000.0;
const double gamma = 7;
const double C0 = 20; // maybe should change this later
const double Pi = 3.14;
const double miu = 0.001;

using namespace std;

class SPH_main;


class SPH_particle
{
public:

	double x[2], v[2];				//position and velocity
	double rho, P;					//density and pressure
	double dedv[2]; // dVdT
	double drho;  // drhodt

	static SPH_main *main_data;		//link to SPH_main class so that it can be used in calc_index

	int list_num[2];	//index in neighbour finding array
	
	void init_particle();
	void calc_index(void);
	void cal_P(void);       //calculate the pressure of particle


};


class SPH_main
{
public:
	SPH_main();

	void set_values(void);
	void initialise_grid(void);

	void place_points(double *min, double *max);

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void cal_derivative(SPH_particle *part);      // calculate the dera

	void smooth_density(SPH_particle* part, int t_step);      // smooth the denstity every 10 to 20 times

	void neighbour_iterate(SPH_particle *part);

	void iteration(SPH_particle* part, int t_step,double dt);


	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell
     

};

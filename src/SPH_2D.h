#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>
#include <algorithm>

const double g = 9.8;
const double rho_0 = 1000.0;
const double gamma_0 = 7;
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

	void cal_derivative(SPH_particle *part);      // calculate the derivative

	void smooth_density(SPH_particle* part, int t_step);      // smooth the denstity every 10 to 20 times

	void neighbour_iterate(SPH_particle *part);

	void iteration(SPH_particle* part, int t_step,double dt);
	void cal_interaction(void);
	void update_parameters_fe(double dt,int step);

	void get_cfl_time_step(double h,  double *part_v, double *other_part_v);
	void get_tf_ta_time_step(double h,  double rho, double *dedv);

	void update_min_time_step(void);
	//void get_min_time_step(double h, double P, double rho, double* part_v, double* other_part_v, double* dedv);
	//void update_min_time_step();

	double min_t_cfl = 10, min_tf = 10, min_ta = 10; // Used in adaptive timestepping, default to long value
	double max_dt  = 10;
	double vij; //used to temporarily store the difference in velocity between particle i and j

	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing

	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell

	};

#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>
#include <algorithm>

const double g = 9.8;
const double rho_0 = 1000.0;
const double ga = 7;
const double C0 = 20; // maybe should change this later
const double Pi = 3.14;
const double miu = 0.001;
const double B = rho_0 * C0 * C0 / ga;

using namespace std;

class SPH_main;


class SPH_particle
{
public:

	double x[2], v[2];				//position and velocity
	double rho, P;					//density and pressure
	double dedv[2]; // dVdT
	double drho;  // drhodt
	int if_topped;

	static SPH_main *main_data;		//link to SPH_main class so that it can be used in calc_index

	int list_num[2];	//index in neighbour finding array
	
	void set_particle_deri(void);
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
	void fill_domain();

	void allocate_to_grid(void);			//allocates all the points to the search grid (assumes that index has been appropriately updated)

	void check_if_topped(SPH_particle* part);

	void cal_derivative(SPH_particle *part);      // calculate the derivative

	void smooth_density(SPH_particle* part, int t_step);      // smooth the denstity every 10 to 20 times

	void neighbour_iterate(SPH_particle *part);
	
	//void update_parameters_pc(double dt, int step);
	void update_parameters_fe(int step, double dt);

	void get_cfl_time_step(double* part_v, double* other_part_v);
	void get_tf_ta_time_step(double rho, double *dedv);

	void update_min_time_step();

	int min_t_cfl, min_tf, min_ta; // Used in adaptive timestepping, default to long value
	int max_dt;
	int vij; //used to temporarily store the difference in velocity between particle i and j

	double h;								//smoothing length
	double h_fac;
	double dx;								//particle initial spacing
	double frange;
	double m;


	double min_x[2], max_x[2];				//dimensions of simulation region

	int max_list[2];

	vector<SPH_particle> particle_list;						//list of all the particles

	vector<vector<vector<SPH_particle*> > > search_grid;		//Outer 2 are the grid, inner vector is the list of pointers in each cell

	};

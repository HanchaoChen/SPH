#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <math.h>
#include <algorithm>

// some constant used in simulation

const double g = 9.81;  //gravity
const double rho_0 = 1000.0;  // the initial reference density
const double ga = 7;   // constant used to calculate pressure
const double C0 = 20; // speed of sound in the system (maybe should change this later)
const double miu = 0.001;  // viscosity factor (:Pa*S)
const double B = rho_0 * C0 * C0 / ga;  // constant used to calculate pressure


using namespace std;

class SPH_main;


class SPH_particle
{
public:

    double x[2], v[2];                //position and velocity               
    double rho, P;                    //density and pressure
    double dedv[2];                   //acceleration
    double drho;                      //derivation of density
    int if_topped;                    // the index to indicate if there is any particle overlapped

    static SPH_main *main_data;        //link to SPH_main class so that it can be used in calc_index

    int list_num[2];    //index in neighbour finding array
    
    void set_particle_deri(void);     // initialize the derivative parameters of a particle
    void init_particle();             // initialize the parameters of a particle
    void calc_index(void);            // calculate the grid index of the particle
    void init_particle_topped(void);  // initialize the if_topped value of a particle

};


class SPH_main
{
public:
    SPH_main();

    void set_values(void); // set up the domain sizeand some physical values
    void initialise_grid(void); //Add boundaries and divide the domain to several grids
    
    void place_points(double *min, double *max); // the initial allocation of grids (use it twice)
    void fill_domain();            // Place points for the boundaries and the fluid domain

    void allocate_to_grid(void);  //allocates all the points to the search grid
                                  //(assumes that index has been appropriately updated)

    void check_if_topped(SPH_particle* part);   // function to check if particles topped

    void cal_derivative(SPH_particle *part, int t);  // calculate the derivative

    void smooth_density(SPH_particle* part, int t_step);   // smooth the density every 10 to 20 times
    
    void update_parameters_fe(int step, double dt);      // update the parameters of particles using forward Euler method

    void get_cfl_time_step(double* part_v, double* other_part_v, int t);   // calculate the CFL_time for the purpose of adapative time stepping
    void get_tf_ta_time_step(double rho, double *dedv, int t);           //   calculate the tf_time and ta_time for the purpose of adapative time stepping
    void update_min_time_step(int t);                    //    find the minimum between tf_time, ta_time, and cfl_time for the purpose of adapative time stepping.

    double min_t_cfl, min_tf, min_ta; // Used in adaptive time stepping, default to long value
    double max_dt;                    // adapetive time stepping value
    double vij; //used to temporarily store the difference in velocity between particle i and j

    double h;                                 // smoothing length
    double h_fac;                             // factor used in calculate smoothing length
    double dx;                                // particle initial spacing
    double frange;                            // the range within which addtional force from boundary added
    double m;                                 // the mass of particles (constant)


    double min_x[2], max_x[2];                //dimensions of simulation region

    int max_list[2];                          // maximum of rows and columns in the simulation region

    vector<SPH_particle> particle_list;                        //list of all the particles

    vector<vector<vector<SPH_particle*> > > search_grid;        //Outer 2 are the grid, inner vector is the list of pointers in each cell

    };

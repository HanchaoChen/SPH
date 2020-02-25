#include "SPH_2D.h"
#include "SPH_2D.cpp"
//#include "file_writer.h"

SPH_main domain;

// the repulsive force's range
double frange = domain.dx;


// push back a particle
void pushback(SPH_particle* part)
{
    if (part->x[0] <= 0.0)
    {
        part->x[0] = -part->x[0] + frange;
    }

    if (part->x[0] >= 20.0)
    {
        part->x[0] = 2 * 20.0 - frange;
    }

    if (part->x[1] <= 0.0)
    {
        part->x[1] = -part->x[1] + frange;
    }

    if (part->x[1] >= 10.0)
    {
        part->x[1] = 2 * 10.0 - part->x[1] - frange;
    }
}


int main(void)
{
    double time = 30.0;
    double dt = 0.1;
    int total_step = int(time / dt)+1;

    domain.set_values(); //Set simulation parameters
    domain.initialise_grid(); //initialise simulation grid

    // set up boundary
    int thick = int(2 * domain.h / domain.dx);

    double min_bottom[2] = { -thick * domain.dx, -thick * domain.dx };
    double max_bottom[2] = { 20.0 + thick * domain.dx, -domain.dx };

    double min_left[2] = { -thick * domain.dx, 0.0 };
    double max_left[2] = { -domain.dx, 10.0 };

    double min_right[2] = { 20.0 + domain.dx, 0.0 };
    double max_right[2] = { 20.0 + thick * domain.h, 10.0 };

    double min_top[2] = { -thick * domain.dx, 10.0 + domain.dx };
    double max_top[2] = { 20.0 + thick * domain.dx, 10.0 + thick * domain.dx };

    // set up fluid
    double min_x1[2] = { 0.0, 0.0 };
    double max_x1[2] = { 20.0, 2.0 };

    double min_x2[2] = { 0.0, 2.0 + domain.dx };
    double max_x2[2] = { 3.0, 5.0 };


    domain.place_points(min_bottom, max_bottom);                //places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
    domain.place_points(min_left, max_left);
    domain.place_points(min_right, max_right);
    domain.place_points(min_top, max_top);

    domain.place_points(min_x1, max_x1);
    domain.place_points(min_x2, max_x2);

    domain.allocate_to_grid();                                    //needs to be called for each time step

    //domain.neighbour_iterate(&domain.particle_list[125]);        //finding all the neighbours of the 100th particle in the list - in reality the simulation loop will need to do the calculations for the neighbours of every particle

    // write_file("example.vtp", &domain.particle_list);
    //domain.get_min_time_step(domain.h,domain.rho, part->rho, part->v, other_part->v, part->dedv);
    // for (int t = 0; t < total_step; t++)
    // {
    //     domain.update_parameters_fe(dt,t);
    //     if(t%10==0)
    //         write_file("example.vtp", &domain.particle_list);
    //     cout << t << endl;
    // }

    return 0;
}

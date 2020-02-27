#include "SPH_2D.h"
#include "file_writer.h"
#include "data_output.h"
SPH_main domain;

// the repulsive force's range


// push back a particle

int main(void)
{
    double time = 1.0;
    double dt = 0.1*1.3*0.2/20;
    int total_step = int(time / dt)+1;
    int interval = 100;

    domain.set_values(); //Set simulation parameters
    domain.initialise_grid(); //initialise simulation grid

                                   //needs to be called for each time step

    //domain.get_min_time_step(domain.h,domain.rho, part->rho, part->v, other_part->v, part->dedv);
    
    domain.fill_domain();
    domain.allocate_to_grid();
    for (int p = 0; p < domain.particle_list.size(); p++) // for all particles in the domain
    {
        domain.particle_list[p].init_particle();
    }
    
    cout << "Please input the simulation time output interval you want.(Default is 100)" << endl;
    cin >> interval;
    while (!cin)
    {
        cout << "It is not an integer!" << endl;
        cin.clear();
        cin.ignore(10000, '\n');
        cout << "Please input the simulation time output interval you want.(Default is 100)" << endl;
        cin >> interval;
    }
    cout << "interval = " << interval << endl;


    for (int t = 0; t < total_step; t++)
    {
        if (t % interval == 0)
        {
            cout << t << endl;
            write_file(t, &domain.particle_list);
            //data_output(t, &domain.particle_list);
        }
            
        domain.update_parameters_fe(t, dt);
        //domain.update_parameters_pc(t,dt);
        domain.allocate_to_grid();
    }

    return 0;
}

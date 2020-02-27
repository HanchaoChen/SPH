#define _USE_MATH_DEFINES

#include <cmath>
#include "SPH_2D.hpp"

SPH_main* SPH_particle::main_data;



void SPH_particle::init_particle()
{
    /*Initialize the parameters of a particle. */

    rho = rho_0;
    P = B * (pow((rho / rho_0), ga) - 1);
    v[0] = 0;
    v[1] = 0;
    if_topped = 0;
}

void SPH_particle::set_particle_deri(void)
{
    /*Initialize the derivative parameters of a particle. */

    dedv[0] = 0;
    dedv[1] = -g;
    drho = 0;
}

void SPH_particle::cal_P(void)
{
    /*Calculate pressure of the particle */

    P = B * (pow((rho / rho_0), ga) - 1);
}

void SPH_particle::calc_index(void)
{
    /*Calculate the grid index of the particle */

    for (int i = 0; i < 2; i++)
        list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

SPH_main::SPH_main()
{
    SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
    /*Set up the domain size and some physical values*/
    min_x[0] = 0.0;
    min_x[1] = 0.0;

    max_x[0] = 20.0;
    max_x[1] = 10.0;

    dx = 0.2;

    h_fac = 1.3;
    h = dx * h_fac;
    frange = dx; //Used????
    
    max_dt = 0.1 * h / C0;
    min_ta = 0.1 * h / C0;
    min_tf = 0.1 * h / C0;
    min_t_cfl = 0.1 * h / C0;
    
    m = rho_0 * dx * dx; //global constant varibles
}


void SPH_main::initialise_grid(void)
{
    /*Add boundaries and divide the domain to several grids*/

    for (int i = 0; i < 2; i++)
    {
        min_x[i] -= 2.0 * h;
        max_x[i] += 2.0 * h;                                                //add buffer for virtual wall particles

        max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
    }

    search_grid.resize(max_list[0]);
    for (int i = 0; i < max_list[0]; i++)
        search_grid[i].resize(max_list[1]);
}


void SPH_main::place_points(double* min, double* max) //the initial allocation (use it twice)
{
    /*
    Assign the grid index to particles and put particles into particle_list.

    @param[in] min The minimum value of the domain on the x, y axis
    @param[in] max The maximum value of the domain on the x, y axis

   */
    double x[2] = { min[0], min[1] };
    SPH_particle particle;

    while (x[0] <= max[0])
    {
        x[1] = min[1];
        while (x[1] <= max[1])
        {
            for (int i = 0; i < 2; i++)
                particle.x[i] = x[i];

            particle.calc_index();

            particle_list.push_back(particle);

            x[1] += dx;
        }
        x[0] += dx;
    }
}

void SPH_main::fill_domain()
{
    /*Place points for the boundaries and the fluid domain. */

    int thick = int(2 * h / dx);

    double min_bottom[2] = { -thick * dx, -thick * dx };
    double max_bottom[2] = { 20.0 + thick * dx, -dx };

    double min_left[2] = { -thick * dx, 0.0 };
    double max_left[2] = { -dx, 10.0 };

    double min_right[2] = { 20.0 + dx, 0.0 };
    double max_right[2] = { 20.0 + thick * h, 10.0 };

    double min_top[2] = { -thick * dx, 10.0 + dx };
    double max_top[2] = { 20.0 + thick * dx, 10.0 + thick * dx };

    // set up fluid
    double min_x1[2] = { 0.0, 0.0 };
    double max_x1[2] = { 20.0, 2.0 };

    double min_x2[2] = { 0.0, 2.0 + dx };
    double max_x2[2] = { 3.0, 5.0 };


    place_points(min_bottom, max_bottom);   //places initial points - will need to be modified to include boundary points and the specifics of where the fluid is in the domain
    place_points(min_left, max_left);
    place_points(min_right, max_right);
    place_points(min_top, max_top);

    place_points(min_x1, max_x1);
    place_points(min_x2, max_x2);
}

void SPH_main::allocate_to_grid(void)
{
    /* Update the points in the grids when the particles have their positions updated */

    for (int i = 0; i < max_list[0]; i++)
        for (int j = 0; j < max_list[1]; j++)
            search_grid[i][j].clear();

    for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
    {
        // Set the derivatives back to initial state
        particle_list[cnt].set_particle_deri();
        particle_list[cnt].calc_index();
        search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
    }
}


void SPH_main::check_if_topped(SPH_particle* part)
{
    /* Check if two particles toppled  */

    SPH_particle* other_part;
    size_t size = search_grid[part->list_num[0]][part->list_num[1]].size();

    for (size_t cnt = 0; cnt < size; cnt++)
    {
        other_part = search_grid[part->list_num[0]][part->list_num[1]][cnt];
        if (other_part != part)
        {
            if (other_part->x[0] == part->x[0])
            {
                if (other_part->x[1] == part->x[1])
                    part->if_topped = 1;
            }
        }

    }
}


void SPH_main::cal_derivative(SPH_particle* part, int t)
{
    /*
    Calculate the derivatives on each particle

    @param[in] part The particle calculated

   */
    SPH_particle* other_part;
    double dist;            //distance between particles
    double dn[2];            //vector from 1st to 2nd particle
    double fac = 10. / (7. * M_PI * h * h);
    double dW;              // the derivative of weight
    double dv[2];           // Vij
    double mdv;             // [Vij, eij] dot product

    // calculate the pressure(P) of particle i
    for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
        if (i >= 0 && i < max_list[0])
            for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
                if (j >= 0 && j < max_list[1])
                {
                    for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                    {
                        other_part = search_grid[i][j][cnt];
                        if (part != other_part)    //stops particle interacting with itself
                        {
                            //Calculates the distance between potential neighbors
                            for (int n = 0; n < 2; n++)
                            {
                                dn[n] = part->x[n] - other_part->x[n];
                                dv[n] = part->v[n] - other_part->v[n];
                            }

                            mdv = sqrt(dv[0] * dv[0] + dv[1] * dv[1]); // scalar value of velocity

                            dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                            /*if (h / mdv < min_t_cfl)
                                min_t_cfl = h / mdv;*/
                            double q = dist / h;

                            if (dist < 2. * h)    //only particle within 2h
                            {
                                //TODO: all the interactions between the particles
                                if (dist < h)
                                {
                                    dW = fac/ h * (-3. * q + 9. / 4 * q * q);
                                }
                                else
                                {
                                    dW = -fac / h * 3. / 4 * pow((2 - q), 2);
                                }
                                part->dedv[0] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[0] / dist + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[0] / dist;
                                part->dedv[1] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[1] / dist + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[1] / dist;
                                part->drho += m * dW * (dn[0] * dv[0] + dv[1] * dn[1]) / dist;


                            }
                            get_cfl_time_step(part->v, other_part->v, t);

                        }
                    }
                }

}


void SPH_main::smooth_density(SPH_particle* part, int t_step)
{
    /*
    Smooth density of particles every fixed time steps

    @param[in] part The particle calculated
    @param[in] t_step The time step.

   */
    if (t_step > 0 && t_step % 15 == 0) // every 15 time steps execute this
    {
        double sumW = 0;
        double sumW_rho = 0;
        SPH_particle* other_part;
        double dist;            //distance between particles
        double dn[2];            //vector from 1st to 2nd particle
        double W;               // weight
        double fac = 10. / (7. * Pi * h * h);

        for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
            if (i >= 0 && i < max_list[0]) // Ensure it is in the domain on x axis
                for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
                    if (j >= 0 && j < max_list[1]) // Ensure it is in the domain on y axis
                    {
                        for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
                        {
                            other_part = search_grid[i][j][cnt];
                            //including particle interacting with itself
                            for (int n = 0; n < 2; n++)
                                dn[n] = part->x[n] - other_part->x[n];
                            //Calculates the distance between potential neighbors
                            dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
                            double q = dist / h;

                            if (dist < 2. * h)            //only particle within 2h
                            {
                                //TODO: all the interactions between the particles
                                if (dist < h)
                                {
                                    W = fac * (1 - 3 / 2 * pow(q, 2) + 3 / 4 * pow(q, 3));
                                }
                                else
                                {
                                    W = fac * 1 / 4 * pow((2 - q), 3);
                                }

                                sumW += W;
                                sumW_rho += W / other_part->rho;
                            }

                        }
                    }
        // Get the average density
        part->rho = sumW / sumW_rho;
    }
}


void SPH_main::update_parameters_fe(int step,double dt) // uses forward euler explicit method to update parameters for all particles
{
    /*
    Calculate the new position and velocity for all the particles
    using forward Euler method.

    @param[in] step The time step
    @param[in] dt The interval of time steps.

   */


    //dt = 0.1 * h / C0;
    /*this->min_t_cfl = 0.1 * h / C0;
    this->min_tf = 0.1 * h / C0;
    this->min_ta = 0.1 * h / C0;*/
    for (int p = 0; p < this->particle_list.size(); p++)// for all particles in the domain
    {
        check_if_topped(&(this->particle_list[p]));
        smooth_density(&(this->particle_list[p]), step);
    }

    for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
    {
        cal_derivative(&(this->particle_list[p]), step);
    }

    for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
    {

        // update velocity and position of the particles within the boundary
        if (this->particle_list[p].x[0] > (this->min_x[0] + 2 * this->h - dx) \
            && particle_list[p].x[0] < (this->max_x[0] - 2 * this->h+dx) \
            && particle_list[p].x[1] < (this->max_x[1] - 2 * this->h + dx) \
            && particle_list[p].x[1] > (this->min_x[1] + 2 * this->h - dx))
        {
            //update this->min_tf
            //update this->min_ta

            get_tf_ta_time_step(this->particle_list[p].rho, this->particle_list[p].dedv, step);
            
            for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions)
                                        //for a single particle, using derivatives computed above
            {
                if (this->particle_list[p].x[i] + dt * this->particle_list[p].v[i] < min_x[i] + 2.0 * h - dx)
                {
                    this->particle_list[p].v[i] = abs(this->particle_list[p].v[i] + this->particle_list[p].dedv[i] * dt);
                }
                else if (this->particle_list[p].x[i] + dt * this->particle_list[p].v[i] > max_x[i] - 2.0 * h + dx)
                {
                    this->particle_list[p].v[i] = -abs(this->particle_list[p].v[i] + this->particle_list[p].dedv[i] * dt);
                }
                else
                {
                    double new_x = this->particle_list[p].x[i] + dt * this->particle_list[p].v[i];
                    double new_v = this->particle_list[p].v[i] + dt * this->particle_list[p].dedv[i];
                    this->particle_list[p].x[i] = new_x;
                    this->particle_list[p].v[i] = new_v;
                }
            }
        }

        //update the density and pressures of all the particles
        double new_rho = this->particle_list[p].rho + dt * this->particle_list[p].drho; // updates density
        this->particle_list[p].rho = new_rho;
        this->particle_list[p].P = B * (pow(particle_list[p].rho / rho_0, ga) - 1);

    }
    
    //if (this->max_dt <= 0.0014) cout << this->max_dt << endl;
    
    /*cout << "this->min_t_cfl"<<this->min_t_cfl << "\t";
    cout << "this->min_t_f" << this->min_tf<<"\t";
    cout << "this->min_t_a" << this->min_ta << "\t";*/
    /*dt = min(this->min_t_cfl, this->min_tf);
    dt = min(this->max_dt, this->min_ta);
    cout << "t_step: "<<step <<" dt: "<<dt <<endl;*/

}





void SPH_main::get_cfl_time_step(double* part_v, double* other_part_v, int t)
{
    // this->min_t_cfl = 0.1*h/C0; //Initial time step
    if (t != 0)
    {
        this->vij = sqrt(pow(part_v[1] - other_part_v[1], 2) + pow(part_v[0] - other_part_v[0], 2));
        if (this->h / this->vij < this->min_t_cfl) {
            this->min_t_cfl = h / this->vij;
        }
        
//        if (this->min_t_cfl != 0) cout << "this->min_t_cfl" << endl;
    }
}

void SPH_main::get_tf_ta_time_step(double rho, double* dedv, int t)
{
    //this->min_tf = 0.1*h/C0;
    //this->min_ta = 0.1*h/C0;
    if (t != 0)
    {
        if (sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2))) < this->min_tf) {
            this->min_tf = sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2)));
        }
        if (h / (C0 * sqrt(std::pow(rho / rho_0, ga - 1))) < this->min_ta) {
            this->min_ta = C0 * sqrt(std::pow(rho / rho_0, ga - 1));
        }
        
//        if (this->min_ta != 0) cout << "this->min_ta" << endl;
//        if (this->min_tf != 0) cout << "this->min_ta" << endl;
    }
}

void SPH_main::update_min_time_step(int t)
{
    if (t != 0)
    {
        this->max_dt = min(this->min_t_cfl, this->min_tf);
        this->max_dt = min(this->max_dt, this->min_ta);
    }
}


// ************ GOES IN SPH_Snippet.cpp
bool triangular = true;
    if (triangular == true)
        domain.slope = (double) (max_x1[1] - min_x1[1]) / (double) (max_x1[0] - min_x1[0]); // the line delineating the fluid from the upward-sloping boundary region








// ************ GOES IN SPH_2D.h
// goes in class SPH_main
bool sloped_boundaries; // tells if domain is sloped
double slope; // slope of the boundary








// *************GOES IN SPH_2D.cpp

void SPH_main::update_parameters_fe(double dt,int step) // uses forward euler explicit method to update parameters for all particles
{
	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{
        // update the particles within the boundary
		double min_y = this->min_x1[1] + 2 * this->h;
        
        if (triangular == true)
            min_y += x[0] * this->slope + this->min_x1[1];

        
		if (this->particle_list[p].x[0] >= this->min_x[0] + 2 * this->h && particle_list[p].x[0] <= this->max_x[0] - 2 * this->h && particle_list[p].x[1] <= this->max_x[1] - 2 * this->h && particle_list[p].x[1] >= min_y)
		{

			smooth_density(&(this->particle_list[p]), step);
			cal_derivative(&(this->particle_list[p]));

			for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions) for a single particle, using derivatives computed above
			{
				double new_x = this->particle_list[p].x[i] + dt * this->particle_list[p].v[i];
				this->particle_list[p].x[i] = new_x;
				double new_v = this->particle_list[p].v[i] + dt * this->particle_list[p].dedv[i]; // need to get this from other person's part of assignment
				this->particle_list[p].v[i] = new_v;
			}
			double new_rho = this->particle_list[p].rho + dt * this->particle_list[p].drho; // updates density
			this->particle_list[p].rho = new_rho;
			pushback(&(this->particle_list[p]));
		}
	}
}

void SPH_main::pushback(SPH_particle* part)
{
	if (part->x[0] <= 0.0)
	{
		part->x[0] = -part->x[0] + frange;
	}

	if (part->x[0] >= 20.0)
	{
		part->x[0] = 2 * 20.0 - frange;
	}

	if (part->x[1] >= 10.0)
	{
		part->x[1] = 2 * 10.0 - part->x[1] - frange;
	}

    if (triangular == true)
    {
        if (part->x[1] < this->min_x[1] + 2 * this->h + (x[0] * this->slope + this->min_x1[1])) // for the sloped boundary domain
        {
            part->x[1] = x[0] * this->slope + this->min_x1[1] // push in back *above* the boundary
        }
    }

    else
    {
        if (part->x[1] <= 0.0)
        {
            part->x[1] = -part->x[1] + frange;
        }
    }
    

	
}
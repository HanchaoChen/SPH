#include "SPH_2D.h"

SPH_main* SPH_particle::main_data;

void SPH_particle::init_particle()
{
	rho = rho_0;
	cal_P();
	v[0] = 0;
	v[1] = 0;
}

void SPH_particle::cal_P(void)
{
	const double B = rho_0 * C0 * C0/gam;
	P = B * (pow((rho / rho_0), gam) - 1);
}

void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;
	 
	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = 0.02;

	h_fac = 1.3;
	h = dx * h_fac;
}

void SPH_main::initialise_grid(void)
{
	for (int i = 0; i < 2; i++)
	{
		min_x[i] -= 2.0 * h;
		max_x[i] += 2.0 * h;												//add buffer for virtual wall particles

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0 * h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}


void SPH_main::place_points(double* min, double* max) //the initial allocation (use it twice)
{
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


void SPH_main::allocate_to_grid(void)				//needs to be called each time that all the particles have their positions updated
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}

// calcluate the derivative on each particle
void SPH_main::cal_derivative(SPH_particle* part)
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	//double W;               // weight 
	double fac = 10. / (7. * Pi * h * h);
	double dW;               // the derivative of weight
	double dv[2]; //Vij
	double mdv;   // [Vij, eij] dot product
	
	part->dedv[0] = 0;
	part->dedv[1] = g;
	part->drho = 0;
	double m = rho_0 * dx * dx;

	part->cal_P(); // calculate the pressure(P) of particle i
	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];
						other_part->cal_P(); // calculate the P of particle j
						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								dn[n] = part->x[n] - other_part->x[n];
								dv[n] = part->v[n] - other_part->v[n];
							}

							mdv = sqrt(dv[0] * dv[0] + dv[1] * dv[1]);

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
							if (h / mdv < min_t_cfl)
								min_t_cfl = h / mdv;
							double q = dist / h;

							if (dist < 2. * h)					//only particle within 2h
							{
								//TODO: all the interactions between the particles
								if (dist < h)
								{
									//W = fac * (1 - 3 / 2 * pow(q, 2) + 3 / 4 * pow(q, 3));
									dW = fac/ h * (-3 * q + 9 / 4 * q * q);
								}
								else
								{
									//W = fac * 1 / 4 * pow((2 - q), 3);
									dW = -fac / h * 3 / 4 * pow((2 - q), 2);
								}

								part->dedv[0] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[0] / dist;
								part->dedv[1] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[1] / dist;
								part->drho += m * dW * mdv;
								
							}

						}
					}
				}
	//cout << "drho" << "\t" << part->drho;
	//cout <<"min_t_cfl: "<< min_t_cfl;
}


void SPH_main::smooth_density(SPH_particle* part, int t_step)
{
	if (t_step > 0 && t_step % 15 == 0) // every 15 time steps execute this  
	{
		double sumW = 0;
		double sumW_rho = 0;
		SPH_particle* other_part;
		double dist;			//distance between particles
		double dn[2];			//vector from 1st to 2nd particle
		double W;               // weight 
		double fac = 10. / (7. * Pi * h * h);

		for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
			if (i >= 0 && i < max_list[0])
				for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
					if (j >= 0 && j < max_list[1])
					{
						for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
						{
							other_part = search_grid[i][j][cnt];
							other_part->cal_P(); // calculate the P of particle j
							//including particle interacting with itself
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];
							//Calculates the distance between potential neighbour
							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
							double q = dist / h;

							if (dist < 2. * h)					//only particle within 2h
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
		part->rho = sumW / sumW_rho;
	}
}


// manipulate on single particle
void SPH_main::neighbour_iterate(SPH_particle* part)					//iterates over all particles within 2h of part - can be made more efficient using a stencil and realising that all interactions are symmetric
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle

	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];

						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
								dn[n] = part->x[n] - other_part->x[n];

							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

							if (dist < 2. * h)					//only particle within 2h
							{
								//TODO: all the interactions between the particles

								cout << "dn: " << dn[0] << " " << dn[1] << endl;		//Should be removed from the code - simply here for you to see that it is working
							}

						}
					}
				}
}

void SPH_main::update_parameters_pc(double dt,int step) // uses forward euler explicit method to update parameters for all particles
{
	//dt = 0.1 * h / C0;
	this->min_t_cfl = 0.1 * h / C0;
	this->min_tf = 0.1 * h / C0;
	this->min_ta = 0.1 * h / C0;

	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{
		// update the particles within the boundary
		
		if (this->particle_list[p].x[0] >= this->min_x[0] + 2 * this->h && particle_list[p].x[0] <= this->max_x[0] - 2 * this->h && particle_list[p].x[1] <= this->max_x[1] - 2 * this->h && particle_list[p].x[1] >= this->min_x[1] + 2 * this->h)
		{

			smooth_density(&(this->particle_list[p]), step);
			cal_derivative(&(this->particle_list[p]));
			
			//update this->min_tf
			if (sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2))) < this->min_tf) {
				this->min_tf = sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2)));
			}
			//update this->min_ta
			if (h / (C0 * sqrt(pow(this->particle_list[p].rho / rho_0, gam - 1))) < this->min_ta) {
				this->min_ta = C0 * sqrt(pow(this->particle_list[p].rho / rho_0, gam - 1));
			}

			double x_mid[2], v_mid[2], rho_mid; // stores average values for predictor-corrector method

			for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions) for a single particle, using derivatives computed previously
			{
				x_mid[i] = this->particle_list[p].x[i] + 0.5 * dt * this->particle_list[p].v[i]; 
				v_mid[i] = this->particle_list[p].v[i] + 0.5 * dt * this->particle_list[p].dedv[i];
				
				this->particle_list[p].x[i] = 2 * x_mid[i] - this->particle_list[p].x[i];
				this->particle_list[p].v[i] = 2 * v_mid[i] - this->particle_list[p].v[i];
			}

			rho_mid = this->particle_list[p].rho + 0.5 * dt * this->particle_list[p].drho; // updates densities for a single particle, using derivatives computed previously
        	this->particle_list[p].rho = 2 * rho_mid - this->particle_list[p].rho;

			pushback(&(this->particle_list[p]));
		}
	}
	// cout << "this->min_t_cfl"<<this->min_t_cfl << "\t";
	// cout << "this->min_t_f" << this->min_tf<<"\t";
	// cout << "this->min_t_a" << this->min_ta << "\t";
	/*dt = min(this->min_t_cfl, this->min_tf);
	dt = min(this->max_dt, this->min_ta);*/
	// cout << "t_step: "<<step <<" dt: "<<dt <<endl;
	
}

void SPH_main::update_parameters_fe(double dt,int step) // uses forward euler explicit method to update parameters for all particles
{
	//dt = 0.1 * h / C0;
	this->min_t_cfl = 0.1 * h / C0;
	this->min_tf = 0.1 * h / C0;
	this->min_ta = 0.1 * h / C0;
	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{
		// update the particles within the boundary
		
		if (this->particle_list[p].x[0] >= this->min_x[0] + 2 * this->h && particle_list[p].x[0] <= this->max_x[0] - 2 * this->h && particle_list[p].x[1] <= this->max_x[1] - 2 * this->h && particle_list[p].x[1] >= this->min_x[1] + 2 * this->h)
		{

			smooth_density(&(this->particle_list[p]), step);
			cal_derivative(&(this->particle_list[p]));
			
			//update this->min_tf
			if (sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2))) < this->min_tf) {
				this->min_tf = sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2)));
			}
			//update this->min_ta
			if (h / (C0 * sqrt(pow(this->particle_list[p].rho / rho_0, gam - 1))) < this->min_ta) {
				this->min_ta = C0 * sqrt(pow(this->particle_list[p].rho / rho_0, gam - 1));
			}
			for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions) for a single particle, using derivatives computed above
			{
				double new_x = this->particle_list[p].x[i] + dt * this->particle_list[p].v[i];
				this->particle_list[p].x[i] = new_x;
				double new_v = this->particle_list[p].v[i] + dt * this->particle_list[p].dedv[i]; // need to get this from other person's part of assignment
				this->particle_list[p].v[i] = new_v;
			}
			double new_rho = this->particle_list[p].rho + dt * this->particle_list[p].drho; // updates density
			this->particle_list[p].rho = new_rho;
			//cout << "density:" << this->particle_list[p].rho;
			pushback(&(this->particle_list[p]));
		}
	}
	cout << "this->min_t_cfl"<<this->min_t_cfl << "\t";
	cout << "this->min_t_f" << this->min_tf<<"\t";
	cout << "this->min_t_a" << this->min_ta << "\t";
	/*dt = min(this->min_t_cfl, this->min_tf);
	dt = min(this->max_dt, this->min_ta);*/
	cout << "t_step: "<<step <<" dt: "<<dt <<endl;
	
}

//void SPH_main::get_min_time_step(double h, double* part_v, double* other_part_v, double* dedv) 
//{ // rho_0, C0, gam are globals
//	this->min_t_cfl = 0.1 * h / C0;
//	this->min_tf = 0.1 * h / C0;
//	this->min_ta = 0.1 * h / C0;
//
//	this->vij = sqrt(pow(part_v[1] - other_part_v[1], 2) + pow(part_v[0] - other_part_v[0], 2));
//	if (h / this->vij < this->min_t_cfl) {
//		this->min_t_cfl = h / this->vij;
//	}
//	if (sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2))) < this->min_tf) {
//		this->min_tf = sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2)));
//	}
//	if (h / (C0 * sqrt(pow(rho / rho_0, gam - 1))) < this->min_ta) {
//		this->min_ta = C0 * sqrt(pow(rho / rho_0, gam - 1));
//	}
//}
//void SPH_main::update_min_time_step() {
//	this->max_dt = min(this->min_t_cfl, this->min_tf);
//	this->max_dt = min(this->max_dt, this->min_ta);
//}

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

	if (part->x[1] <= 0.0)
	{
		part->x[1] = -part->x[1] + frange;
	}

	if (part->x[1] >= 10.0)
	{
		part->x[1] = 2 * 10.0 - part->x[1] - frange;
	}
}
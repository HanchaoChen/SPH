#define _USE_MATH_DEFINES

#include <cmath>
#include "SPH_2D.h"

SPH_main* SPH_particle::main_data;

void SPH_particle::init_particle()
{
	rho = rho_0;
	P = B * (pow((rho / rho_0), ga) - 1);
	v[0] = 0;
	v[1] = 0;

}

void SPH_particle::set_particle_deri(void)
{
	dedv[0] = 0;
	dedv[1] = -g;
	drho = 0;
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

	dx = 0.2;

	h_fac = 1.3;
	h = dx * h_fac;
	frange = 0.5*dx;
	max_dt = double(0.1 * h / C0);
	min_t_cfl= double( 0.1 * h / C0);
	

	m = rho_0 * dx * dx; //global constant varibles
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
		particle_list[cnt].set_particle_deri();
		particle_list[cnt].calc_index();
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
}

// to check if two particles toppled 
void SPH_main::check_if_topped(void)
{
	SPH_particle* part;
	SPH_particle* other_part;
	/*for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		part = &particle_list[cnt];
		for (unsigned int i=0;i<)
		particle_list_num[0]][particle_list[cnt].list_num[1]]
		

	}*/
}

// calcluate the derivative on each particle
void SPH_main::cal_derivative(SPH_particle* part)
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	//double W;               // weight 
	double fac = 10. / (7. * M_PI * h * h);
	double dW;               // the derivative of weight
	double dv[2]; //Vij
	double mdv;   // [Vij, eij] dot product

	// calculate the pressure(P) of particle i
	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1])
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];
						//other_part->cal_P(); // calculate the P of particle j
						if (part != other_part)					//stops particle interacting with itself
						{
							//Calculates the distance between potential neighbours
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

							if (dist < 2. * h)					//only particle within 2h
							{
								//TODO: all the interactions between the particles
								if (dist < h)
								{
									//W = fac * (1 - 3 / 2 * pow(q, 2) + 3 / 4 * pow(q, 3));
									dW = fac / h * (-3. * q + 9. / 4 * q * q);
								}
								else
								{
									//W = fac * 1 / 4 * pow((2 - q), 3);
									dW = -fac / h * 3. / 4 * pow((2 - q), 2);
								}
								part->dedv[0] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[0] / dist + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[0] / dist;
								part->dedv[1] += -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[1] / dist + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[1] / dist;
								part->drho += m * dW * (dn[0] * dv[0] + dv[1] * dn[1]) / dist;


							}

						}
					}
				}
	//cout << "drho" << "\t" << part->drho;

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
							//other_part->cal_P(); // calculate the P of particle j
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




void SPH_main::update_parameters_fe(int step, double dt) // uses forward euler explicit method to update parameters for all particles
{
	//dt = 0.1 * h / C0;
	/*this->min_t_cfl = 0.1 * h / C0;
	this->min_tf = 0.1 * h / C0;
	this->min_ta = 0.1 * h / C0;*/
	for (int p = 0; p < this->particle_list.size(); p++)
	{
		smooth_density(&(this->particle_list[p]), step);
	}

	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{


		cal_derivative(&(this->particle_list[p]));
	}

	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{

		// update velocity and position of the particles within the boundary
		if (this->particle_list[p].x[0] > (this->min_x[0] + 2 * this->h - dx) && particle_list[p].x[0] < (this->max_x[0] - 2 * this->h + dx) && particle_list[p].x[1] < (this->max_x[1] - 2 * this->h + dx) && particle_list[p].x[1] > (this->min_x[1] + 2 * this->h - dx))
		{
			//update this->min_tf
			if (sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2))) < this->min_tf) {
				this->min_tf = sqrt(h / sqrt(pow(this->particle_list[p].dedv[1], 2) + pow(this->particle_list[p].dedv[0], 2)));
			}
			//update this->min_ta
			if (h / (C0 * sqrt(pow(this->particle_list[p].rho / rho_0, ga - 1))) < this->min_ta) {
				this->min_ta = C0 * sqrt(pow(this->particle_list[p].rho / rho_0, ga - 1));
			}
			for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions) for a single particle, using derivatives computed above
			{
				//add additional force to the particle which are very close to boundary
				if (this->particle_list[p].x[i]< this->min_x[i] + 2.0 * h - dx+ frange)
				{
					this->particle_list[p].dedv[i] += abs(B * (pow(0.01, ga) - 1) * dx * dx * (frange / abs(particle_list[p].x[i] - min_x[i]) - 1));
				}
				if (this->particle_list[p].x[i] > this->max_x[i] - 2.0 * h + dx - frange)
				{
					this->particle_list[p].dedv[i] -= abs(B * (pow(0.01, ga) - 1) * dx * dx * (frange / abs(particle_list[p].x[i] - min_x[i]) - 1));
				}
				// if this particle is very close to the wall points
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
					double new_v = this->particle_list[p].v[i] + dt * this->particle_list[p].dedv[i]; // need to get this from other person's part of assignment
					this->particle_list[p].x[i] = new_x;
					this->particle_list[p].v[i] = new_v;
				}
			}
			//pushback_position(&(this->particle_list[p]));
		}

		//update the density and pressures of all the particles
		double new_rho = this->particle_list[p].rho + dt * this->particle_list[p].drho; // updates density
		this->particle_list[p].rho = new_rho;
		this->particle_list[p].P = B * (pow(particle_list[p].rho / rho_0, ga) - 1);

	}
	/*cout << "this->min_t_cfl"<<this->min_t_cfl << "\t";
	cout << "this->min_t_f" << this->min_tf<<"\t";
	cout << "this->min_t_a" << this->min_ta << "\t";*/
	/*dt = min(this->min_t_cfl, this->min_tf);
	dt = min(this->max_dt, this->min_ta);
	cout << "t_step: "<<step <<" dt: "<<dt <<endl;*/

}





void SPH_main::get_cfl_time_step(double* part_v, double* other_part_v)
{
	// this->min_t_cfl = 0.1*h/C0; //Initial time step 
	this->vij = sqrt(pow(part_v[1] - other_part_v[1], 2) + pow(part_v[0] - other_part_v[0], 2));
	if (this->h / this->vij < this->min_t_cfl) {
		this->min_t_cfl = h / this->vij;
	}

}

void SPH_main::get_tf_ta_time_step(double rho, double* dedv) {
	//this->min_tf = 0.1*h/C0;
	//this->min_ta = 0.1*h/C0;
	if (sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2))) < this->min_tf) {
		this->min_tf = sqrt(h / sqrt(pow(dedv[1], 2) + pow(dedv[0], 2)));
	}
	if (h / (C0 * sqrt(std::pow(rho / rho_0, ga - 1))) < this->min_ta) {
		this->min_ta = C0 * sqrt(std::pow(rho / rho_0, ga - 1));
	}
}

void SPH_main::update_min_time_step() {
	this->max_dt = min(this->min_t_cfl, this->min_tf);
	this->max_dt = min(this->max_dt, this->min_ta);
}






//////////////////////////////////////////////////////////////////////////

void SPH_main::cal_derivative_second(SPH_particle* part)
{
	SPH_particle* other_part;
	double dist;			//distance between particles
	double dn[2];			//vector from 1st to 2nd particle
	//double W;               // weight 
	double fac = 10. / (7. * Pi * h * h);
	double dW;               // the derivative of weight
	double dv[2]; //Vij
	double mdv;   // [Vij, eij] dot product
	 // calculate the pressure(P) of particle i
	

	// temporary increment while looping the other particles
	double incre_v[2]; 
	
	
	for (int i = part->list_num[0] - 1; i <= part->list_num[0] + 1; i++)
		if (i >= 0 && i < max_list[0])
			for (int j = part->list_num[1]-1; j < part->list_num[1] + 1; j++) // changed for 5 point stencil, used to be: for (int j = part->list_num[1] - 1; j <= part->list_num[1] + 1; j++)
				if (j >= 0 && j < max_list[1] && !((i == part->list_num[0] - 1) && (j == part->list_num[1]))) // added condition to implement 5 point stencil
				{
					for (unsigned int cnt = 0; cnt < search_grid[i][j].size(); cnt++)
					{
						other_part = search_grid[i][j][cnt];
						// calculate the P of particle j
						if (part != other_part)					//stops particle interacting with itself
						{
							get_cfl_time_step(part->v, other_part->v); // compute min cfl dt
							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								dn[n] = part->x[n] - other_part->x[n];
								dv[n] = part->v[n] - other_part->v[n];
							}

							mdv = sqrt(dv[0] * dv[0] + dv[1] * dv[1]);
							dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);
							double q = dist / h;

							if ((dist < 2. * h) && (dist > 0.01*dx) && (dn[1] > 0.01*dx))					//only particle within 2h && (for 5 point stencil, below current particle)
							{
								//TODO: all the interactions between the particles
								if (dist < h)
								{
									//W = fac * (1 - 3 / 2 * pow(q, 2) + 3 / 4 * pow(q, 3));
									dW = fac / h * (-3. * q + 9. / 4. * q * q);
								}
								else
								{
									//W = fac * 1 / 4 * pow((2 - q), 3);
									dW = -fac / h * 3. / 4. * pow((2 - q), 2);
								}
								/*if (abs(dW)> 10000) 
									cout << "dw:"<<dW<<"\t";*/
								if (mdv == 0)
								{
									incre_v[0] = 0;
									incre_v[0] = 0;
								}
								else
								{
									incre_v[0] = -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[0] / dist + miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[0] / dist;
									incre_v[1] = -m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dW * dn[1] /  dist+miu * m * (1 / (part->rho * part->rho) + 1 / (other_part->rho * other_part->rho)) * dW * dv[1] / dist;
								}
								part->dedv[0] += incre_v[0];
								part->dedv[1] += incre_v[1];
								other_part->dedv[0] -= incre_v[0];  // 5 point stencil
								other_part->dedv[1] -= incre_v[1];  // 5 point stencil
								part->drho += m * dW * (dn[0] * dv[0] + dv[1] * dn[1]) / dist;
								other_part->drho += m * dW * (dn[0] * dv[0] + dv[1] * dn[1]) / dist; // 5 point stencil 
							}

						}
					}
					get_tf_ta_time_step(part->rho, part->dedv);
				}
	// When called to iterate through all particles, call update_min_time_step 
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void SPH_main::update_parameters_fe_second(int step,double dt) // uses forward euler explicit method to update parameters for all particles
{
	
	for (int p = 0; p < this->particle_list.size(); p++) // for all particles in the domain
	{
		
		smooth_density(&(this->particle_list[p]), step);
		cal_derivative_second(&(this->particle_list[p]));

		// update the particles within the boundary
		if (this->particle_list[p].x[0] >= this->min_x[0] + 2 * this->h && particle_list[p].x[0] <= this->max_x[0] - 2 * this->h && particle_list[p].x[1] <= this->max_x[1] - 2 * this->h && particle_list[p].x[1] >= this->min_x[1] + 2 * this->h)
		{			
			for (int i = 0; i < 2; i++) // updates position and velocity parameters (in x and y directions) for a single particle, using derivatives computed above
			{
				if (this->particle_list[p].x[i] + dt * this->particle_list[p].v[i] < min_x[i] + 2.0 * h)
				{
					this->particle_list[p].v[i] = abs(this->particle_list[p].v[i] + this->particle_list[p].dedv[i] * dt);
				}
				else if (this->particle_list[p].x[i] + dt * this->particle_list[p].v[i] > max_x[i] - 2.0 * h)
				{
					this->particle_list[p].v[i] = -abs(this->particle_list[p].v[i] + this->particle_list[p].dedv[i] * dt);
				}
				else
				{
					double new_x = this->particle_list[p].x[i] + dt * this->particle_list[p].v[i];
					double new_v = this->particle_list[p].v[i] + dt * this->particle_list[p].dedv[i]; // need to get this from other person's part of assignment
					this->particle_list[p].x[i] = new_x;
					this->particle_list[p].v[i] = new_v;
				}
			}

		}
		//update the density and pressures of all the particles

		double new_rho = this->particle_list[p].rho + dt * this->particle_list[p].drho; // updates density
		this->particle_list[p].rho = new_rho;
		this->particle_list[p].P = B * (pow(new_rho / rho_0, ga) - 1);
		/*if (this->particle_list[p].rho < 0)
		{
			cout << "time step" << step;
			cout << "denstity:" << this->particle_list[p].rho << "\t";
			cout << "De denstity:" << this->particle_list[p].drho << "\t";
		}*/		
		
	}
	update_min_time_step();
	//cout << "dt: max_dt" << this->max_dt;
}


// Improvement Method

//////////////////////////////////////////////////////////////////////////////////////
double SPH_main::rforce(double r)
{
	double force;
	int m = 2;
	double presure = 0.01 * rho_0 * pow(C0, 2) / ga;
	if (r <= frange && r > 0.1 * frange)
	{
		force = presure * (pow((frange / r), m) - 1);
		return force;
	}
	if (r <= 0.1 * frange)
	{
		force = presure * (pow(0.1, m) - 1);
		return force;
	}
	return 0;
}

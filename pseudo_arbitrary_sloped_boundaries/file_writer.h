#pragma once

#include <vector>
#include "SPH_2D.h"

void filename_padLeft(string& file_num, const size_t size = 8,
	       const char padzero = '0');

int write_file(int& iter_num,
	       std::vector<SPH_particle> *particle_list, bool test = false);


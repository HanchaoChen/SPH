#include "SPH_2D.h"
#include "file_writer.h"

int main() {

    int test_num = 111111;
  SPH_particle particle;
  std::vector<SPH_particle> particle_list;

  // Populate particle list
  for (int i=0; i<10; i++) {
    particle.x[0]=i;
    particle.x[1]=i;
    particle.v[0] = i;
    particle.v[1] = -i;
    particle.P = i;
    particle_list.push_back(particle);
  }

  if (particle_list.size() != 10) return 1;

  return write_file(test_num, &particle_list);

}


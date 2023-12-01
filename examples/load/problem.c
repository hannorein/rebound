/**
 * Close encounters with simple MERCURIUS (Hernandez 2023)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "integrator_trace.h"
#include "simulationarchive.h"
#include "particle.h"

void heartbeat(struct reb_simulation* r);

int main(int argc, char* argv[]){
  system("rm -rf restart.bin");
  struct reb_simulation* r = reb_create_simulation();

  r->dt = (8./365.) * 2. *M_PI;
  r->integrator = REB_INTEGRATOR_TRACE;
  r->ri_trace.hillfac = 4;            // By default the switching radius is 4 times the hill radius, from Hernandez 2023
  r->N_active = 2;

  // Initialize masses
  struct reb_particle star = {0};
  star.m = 1;
  struct reb_particle jup = {0};
  jup.m = 0.01 / (star.m - 0.01);

  // velocities
  double a = 5.2;
  double e = 0.0;
  star.x = -(jup.m / (star.m + jup.m)) * (a * (1 + e));
  star.vy = -(jup.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
  reb_add(r, star);

  jup.x = (star.m / (star.m + jup.m)) * (a * (1 + e));
  jup.vy = (star.m / (star.m + jup.m)) * sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
  reb_add(r, jup);

  // Test particle
  struct reb_particle test = {0};
  double xhel = 4.42;
  double vhel = 0.0072 * (365.25) * (1 / (2 * M_PI)); // days to REBOUND years

  test.x = xhel + star.x;
  test.vy = vhel + star.vy;
  reb_add(r, test);
  reb_move_to_com(r);

  // Copy simulation
  struct reb_simulation* r2 = reb_copy_simulation(r);

  reb_integrate(r, 1e3*2.*M_PI);
  reb_integrate(r, 2e3*2.*M_PI);

  reb_integrate(r2, 1e3*2.*M_PI);
  reb_output_binary(r2, "restart.bin");
  reb_free_simulation(r2);
  r2 = NULL;

  struct reb_simulation* r3 = reb_create_simulation_from_binary("restart.bin");
  reb_integrate(r3, 2e3*2.*M_PI);

  printf("dx0: %e\n", r->particles[0].x - r3->particles[0].x);
  printf("dy0: %e\n", r->particles[0].y - r3->particles[0].y);
  printf("dx1: %e\n", r->particles[1].x - r3->particles[1].x);
  printf("dy1: %e\n", r->particles[1].y - r3->particles[1].y);
  printf("dx2: %e\n", r->particles[2].x - r3->particles[2].x);
  printf("dy2: %e\n", r->particles[2].y - r3->particles[2].y);
  printf("dvx0: %e\n", r->particles[0].vx - r3->particles[0].vx);
  printf("dvy0: %e\n", r->particles[0].vy - r3->particles[0].vy);
  printf("dvx1: %e\n", r->particles[1].vx - r3->particles[1].vx);
  printf("dvy1: %e\n", r->particles[1].vy - r3->particles[1].vy);
  printf("dvx2: %e\n", r->particles[2].vx - r3->particles[2].vx);
  printf("dvy2: %e\n", r->particles[2].vy - r3->particles[2].vy);
  printf("%f,%f\n", r->t, r3->t);
  reb_free_simulation(r);
  reb_free_simulation(r3);
}


void heartbeat(struct reb_simulation* r){
}

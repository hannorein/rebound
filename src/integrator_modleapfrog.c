/**
 * integrator_modleapfrog.c: A time regularized Leap Frog integator
 *
 */

#include "rebound.h"
#include "integrator_modleapfrog.h"
#include "binarydata.h"
#include "math.h"
#include "integrator_leapfrog.h"


void reb_integrator_modleapfrog_step(struct reb_simulation* r, void* state);
void* reb_integrator_modleapfrog_create();
void reb_integrator_modleapfrog_free(void* p);
const struct reb_binarydata_field_descriptor reb_integrator_modleapfrog_field_descriptor_list[];

const struct reb_integrator reb_integrator_modleapfrog = {
    .documentation = 
    "Time regularized leapfrog integrator."
    ,
    .step = reb_integrator_modleapfrog_step,
    .create = reb_integrator_modleapfrog_create,
    .free = reb_integrator_modleapfrog_free,
    .field_descriptor_list = reb_integrator_modleapfrog_field_descriptor_list,
};

const struct reb_binarydata_field_descriptor reb_integrator_modleapfrog_field_descriptor_list[] = {
    { "Order of the integrator. Default is 2. Other allowed values are 4, 6 and 8.",
        REB_UINT,        "order",          offsetof(struct reb_integrator_modleapfrog_state, order), 0, 0, 0},
    { 0 }, // Null terminated list
};

void* reb_integrator_modleapfrog_create(){
    struct reb_integrator_modleapfrog_state* modleapfrog = calloc(sizeof(struct reb_integrator_modleapfrog_state),1);
    modleapfrog->order = 2;
    modleapfrog->E_0 = nan("");
    return modleapfrog;
}

void reb_integrator_modleapfrog_free(void* p){
    struct reb_integrator_modleapfrog_state* modleapfrog = p;
    free(modleapfrog);
}


static double potential(struct reb_simulation* r){
    // taken from tools.c's reb_simulation_energy
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX)?N:r->N_active;
    const struct reb_particle* restrict const particles = r->particles;
    double e_pot = 0.;
    size_t N_interact = (r->testparticle_type==0)?N_active:N;
    for (size_t i=0;i<N_active;i++){
        struct reb_particle pi = particles[i];
        for (size_t j=i+1;j<N_interact;j++){
            struct reb_particle pj = particles[j];
            double dx = pi.x - pj.x;
            double dy = pi.y - pj.y;
            double dz = pi.z - pj.z;
            e_pot -= r->G*pj.m*pi.m/sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
    return e_pot;
}

static double kinetic(struct reb_simulation* r){
    // taken from tools.c's reb_simulation_energy
    const size_t N = r->N;
    const size_t N_active = (r->N_active==SIZE_MAX)?N:r->N_active;
    const struct reb_particle* restrict const particles = r->particles;
    double e_kin = 0.;
    size_t N_interact = (r->testparticle_type==0)?N_active:N;
    for (size_t i=0;i<N_interact;i++){
        struct reb_particle pi = particles[i];
        e_kin += 0.5 * pi.m * (pi.vx*pi.vx + pi.vy*pi.vy + pi.vz*pi.vz);
    }
    return e_kin;
}

static double drift(struct reb_simulation* r, struct reb_integrator_modleapfrog_state* modleapfrog, double dt){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;
    double E_0 = modleapfrog->E_0;
    dt /= (kinetic(r) - E_0); // use ficticious time to determine actual time step

#pragma omp parallel for schedule(guided)
    for (size_t i=0;i<N;i++){
        particles[i].x  += dt * particles[i].vx;
        particles[i].y  += dt * particles[i].vy;
        particles[i].z  += dt * particles[i].vz;
    }
    r->t += dt; // kick step advanced time so that force evaluations are correct.
    return dt;
}

static void kick(struct reb_simulation* r, struct reb_integrator_modleapfrog_state* modleapfrog, double dtau){
    const size_t N = r->N;
    struct reb_particle* restrict const particles = r->particles;

    double dt = -1 * dtau/potential(r); // use ficticious time to determine actual time step

#pragma omp parallel for schedule(guided)
    for (size_t i=0;i<N;i++){
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
    }
}

// Leapfrog integrator (Drift-Kick-Drift)
// for non-rotating frame.
void reb_integrator_modleapfrog_step(struct reb_simulation* r, void* state){

    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;
    struct reb_integrator_modleapfrog_state* modleapfrog = state;

    const double dtau = r->dt;
    
    // If we haven't set the initial energy yet, do it.
    if (isnan(modleapfrog->E_0)){
        modleapfrog->E_0 = reb_simulation_energy(r);
    }

    double dt = 0;
    switch (modleapfrog->order){
        case 2: // 1*2 FE
            dt += drift(r, modleapfrog, dtau*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau);
            dt += drift(r, modleapfrog, dtau*0.5);
            break;
        case 4: // 3*2 FE
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf4_a);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*2.*reb_integrator_leapfrog_lf4_a);
            dt += drift(r, modleapfrog, dtau*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*(1.-4.*reb_integrator_leapfrog_lf4_a));
            dt += drift(r, modleapfrog, dtau*(0.5-reb_integrator_leapfrog_lf4_a));
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*2.*reb_integrator_leapfrog_lf4_a);
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf4_a);
            break;
        case 6: // 9*2 FE
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[1]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[2]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[3]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[4]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[3]+reb_integrator_leapfrog_lf6_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[3]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[2]+reb_integrator_leapfrog_lf6_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[2]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[1]+reb_integrator_leapfrog_lf6_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[1]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf6_a[0]+reb_integrator_leapfrog_lf6_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]);
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf6_a[0]*0.5);
            break; 
        case 8: // 17*2 FE
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[1]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[2]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[3]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[4]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[5]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[6]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[7]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[8]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[7]+reb_integrator_leapfrog_lf8_a[8])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[7]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[6]+reb_integrator_leapfrog_lf8_a[7])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[6]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[5]+reb_integrator_leapfrog_lf8_a[6])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[5]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[4]+reb_integrator_leapfrog_lf8_a[5])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[4]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[3]+reb_integrator_leapfrog_lf8_a[4])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[3]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[2]+reb_integrator_leapfrog_lf8_a[3])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[2]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[1]+reb_integrator_leapfrog_lf8_a[2])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[1]);
            dt += drift(r, modleapfrog, dtau*(reb_integrator_leapfrog_lf8_a[0]+reb_integrator_leapfrog_lf8_a[1])*0.5);
            reb_simulation_update_acceleration(r);
            kick(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]);
            dt += drift(r, modleapfrog, dtau*reb_integrator_leapfrog_lf8_a[0]*0.5);
            break;
        default:
            reb_simulation_error(r, "Leapfrog order not supported.");
            return;
    }

    // Assert that we cannot give an exact finish time
    r->exact_finish_time = 0; 
    
    r->dt_last_done = dt;
}

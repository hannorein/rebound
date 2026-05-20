#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <stdbool.h>
#include <string.h>

void setup_sim(char* integrator, int corrector){
    struct reb_simulation* r = reb_simulation_create();
    // Setup constants
    r->dt = 6.0/365.25*2*M_PI;
    r->G = 1.;
    r->exact_finish_time = 0;
    reb_simulation_add_fmt(r, "solarsystem");
    reb_simulation_set_integrator(r, integrator);
    if (strcmp(integrator,"asm512")==0){
        struct reb_integrator_asm512_state* asm512 = r->integrator.state;
        asm512->gr_potential = 0;
        asm512->concatenate_steps = 1e4;
        asm512->corrector = corrector;
    }
    if (strcmp(integrator,"whfast512")==0){
        struct reb_integrator_whfast512_state* whfast512 = r->integrator.state;
        whfast512->gr_potential = 0;
    }
    if (strcmp(integrator,"whfast")==0){
        struct reb_integrator_whfast_state* whfast = r->integrator.state;
        whfast->safe_mode = 0;
        whfast->corrector = corrector;
    }
    double E0 = reb_simulation_energy(r);
    reb_simulation_integrate(r, 1e2*M_PI*2);
    double E1 = reb_simulation_energy(r);
    printf("integrator=%s, corrector=%d: \t%e\n", integrator, corrector, fabs((E0-E1)/E0));
    reb_simulation_free(r);
}



extern void reb_integrator_asm512_kepler_step(struct reb_simulation* const r, int N_steps);
extern uint64_t reb_asm512_counter(struct reb_simulation* r, int test_p);

double next_k = 0.0;
void printhex(char* s, double v){
    uint64_t raw;
    memcpy(&raw, &v, sizeof(double));
    printf("    .quad 0x%016llX",(unsigned long long)raw);
    if (*s=='k'){
        printf("    # Kepler\n");
    }else{
        printf("    # Interaction\n");
    }
}

void corrector_Z(double a, double b){
            printhex("k",(a+next_k));
            printhex("i",(-b));
            printhex("k",(-2.*a));
            printhex("i",(b));
            printhex("i    2x",(2*b));
            next_k = a;
            //printf("k %.16f\n",a);
}
static const double reb_whfast_corrector_a_1 = 0.41833001326703777398908601289259374469640768464934;
static const double reb_whfast_corrector_a_2 = 0.83666002653407554797817202578518748939281536929867;
static const double reb_whfast_corrector_a_3 = 1.2549900398011133219672580386777812340892230539480;
static const double reb_whfast_corrector_a_4 = 1.6733200530681510959563440515703749787856307385973;
static const double reb_whfast_corrector_a_5 = 2.0916500663351888699454300644629687234820384232467;
static const double reb_whfast_corrector_a_6 = 2.5099800796022266439345160773555624681784461078960; 
static const double reb_whfast_corrector_a_7 = 2.9283100928692644179236020902481562128748537925454;
static const double reb_whfast_corrector_a_8 = 3.3466401061363021919126881031407499575712614771947;
static const double reb_whfast_corrector_b_178 = 0.093056103771425958591541059067553547100903397724386; 
static const double reb_whfast_corrector_b_177 = -0.065192863576377893658290760803725762027864651086787; 
static const double reb_whfast_corrector_b_176 = 0.032422198864713580293681523029577130832258806467604; 
static const double reb_whfast_corrector_b_175 = -0.012071760822342291062449751726959664253913904872527; 
static const double reb_whfast_corrector_b_174 = 0.0033132577069380655655490196833451994080066801611459; 
static const double reb_whfast_corrector_b_173 = -0.00063599983075817658983166881625078545864140848560259; 
static const double reb_whfast_corrector_b_172 = 0.000076436355227935738363241846979413475106795392377415; 
static const double reb_whfast_corrector_b_171 = -0.0000043347415473373580190650223498124944896789841432241; 

int main(int argc, char* argv[]) {
    setup_sim("asm512", 0);
    setup_sim("asm512", 17);
    setup_sim("whfast", 0);
    setup_sim("whfast", 17);
    setup_sim("whfast512", 0);

    double inv = 1.;
    double dt = 1;
        corrector_Z(-reb_whfast_corrector_a_8*dt,-inv*reb_whfast_corrector_b_171*dt);
        corrector_Z(-reb_whfast_corrector_a_7*dt,-inv*reb_whfast_corrector_b_172*dt);
        corrector_Z(-reb_whfast_corrector_a_6*dt,-inv*reb_whfast_corrector_b_173*dt);
        corrector_Z(-reb_whfast_corrector_a_5*dt,-inv*reb_whfast_corrector_b_174*dt);
        corrector_Z(-reb_whfast_corrector_a_4*dt,-inv*reb_whfast_corrector_b_175*dt);
        corrector_Z(-reb_whfast_corrector_a_3*dt,-inv*reb_whfast_corrector_b_176*dt);
        corrector_Z(-reb_whfast_corrector_a_2*dt,-inv*reb_whfast_corrector_b_177*dt);
        corrector_Z(-reb_whfast_corrector_a_1*dt,-inv*reb_whfast_corrector_b_178*dt);
        corrector_Z(reb_whfast_corrector_a_1*dt,inv*reb_whfast_corrector_b_178*dt);
        corrector_Z(reb_whfast_corrector_a_2*dt,inv*reb_whfast_corrector_b_177*dt);
        corrector_Z(reb_whfast_corrector_a_3*dt,inv*reb_whfast_corrector_b_176*dt);
        corrector_Z(reb_whfast_corrector_a_4*dt,inv*reb_whfast_corrector_b_175*dt);
        corrector_Z(reb_whfast_corrector_a_5*dt,inv*reb_whfast_corrector_b_174*dt);
        corrector_Z(reb_whfast_corrector_a_6*dt,inv*reb_whfast_corrector_b_173*dt);
        corrector_Z(reb_whfast_corrector_a_7*dt,inv*reb_whfast_corrector_b_172*dt);
        corrector_Z(reb_whfast_corrector_a_8*dt,inv*reb_whfast_corrector_b_171*dt);
            printhex("k",(next_k));
}

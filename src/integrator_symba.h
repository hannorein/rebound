#ifndef INTEGRATOR_SYMBA_H
#define INTEGRATOR_SYMBA_H

extern const struct reb_integrator reb_integrator_symba;

struct reb_integrator_symba_state {
    unsigned int levels;
    double r_outer;
    double ratio;
    unsigned int optimized;
};

#endif
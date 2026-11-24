/**
 * @file 	integrator_modleapfrog.h
 * @brief 	Time regularized leapfrog integrator
 * @author 	Stuart Williamson <stuart.williamson@mail.utoronto.ca>
 */

#ifndef _INTEGRATOR_MODLEAPFROG_H
#define _INTEGRATOR_MODLEAPFROG_H

extern const struct reb_integrator reb_integrator_modleapfrog;

struct reb_integrator_modleapfrog_state {
    unsigned int order;
    double E_0;
    unsigned int regularized;
};


#endif

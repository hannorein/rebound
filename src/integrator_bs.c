/**
 * @file 	integrator.c
 * @brief 	BS integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Gragg-Bulirsch-Stoer integration scheme.  
 *          It is a reimplementation of the fortran code by E. Hairer and G. Wanner.
 *          The starting point was the JAVA implementation in hipparchus:
 *          https://github.com/Hipparchus-Math/hipparchus/blob/master/hipparchus-ode/src/main/java/org/hipparchus/ode/nonstiff/GraggBulirschStoerIntegrator.java
 *
 * @section 	LICENSE
 * Copyright (c) 2021 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright (c) 2004, Ernst Hairer
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 * 
 *  - Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 * BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> // memset
#include <float.h> // for DBL_MAX
#include "rebound.h"
#include "gravity.h"
#include "integrator.h"
#include "integrator_bs.h"
#include "integrator_trace.h"
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

//#define DEBUG 0 // set to 1 to print out debug information (reason for step rejection)
    
// Default configuration parameter. 
// They are hard coded here because it
// is unlikely that these need to be changed by the user.
// static const int maxOrder = 18;// was 18 
static const int sequence_length = 9; // = maxOrder / 2; 
static const double stepControl1 = 0.65;
static const double stepControl2 = 0.94;
static const double stepControl3 = 0.02;
static const double stepControl4 = 4.0;
static const double orderControl1 = 0.8;
static const double orderControl2 = 0.9;
static const double stabilityReduction = 0.5;
static const int maxIter = 2; // maximal number of iterations for which checks are performed
static const int maxChecks = 1; // maximal number of checks for each iteration

void reb_integrator_bs_update_particles(struct reb_simulation* r, const double* y){
    if (r==NULL){
        reb_simulation_error(r, "Update particles called without valid simulation pointer.");
        return;
    }
    if (y==NULL){
        reb_simulation_error(r, "Update particles called without valid y pointer.");
        return;
    }
    for (unsigned int i=0; i<r->N; i++){
        struct reb_particle* const p = &(r->particles[i]);
        p->x  = y[i*6+0];
        p->y  = y[i*6+1];
        p->z  = y[i*6+2];
        p->vx = y[i*6+3];
        p->vy = y[i*6+4];
        p->vz = y[i*6+5];
    }
}


static int tryStep(struct reb_simulation* r, const int Ns, const int k, const int n, const double t0, const double step) {
    struct reb_integrator_bs* ri_bs = &r->ri_bs;
    struct reb_ode** odes = r->odes;
    const double subStep  = step / n;
    double t = t0;
    int needs_nbody = ri_bs->user_ode_needs_nbody;
    if (r->integrator == REB_INTEGRATOR_TRACE){
        needs_nbody = 0; // TRACE does not allow for coupling of N-body and other ODEs
    }

    // LeapFrog Method did not seem to be of any advantage 
    //    switch (method) {
    //        case 0: // LeapFrog
    //            {
    //                // first substep
    //                for (int s=0; s < Ns; s++){
    //                    double* y0 = odes[s].y;
    //                    double* y1 = odes[s].y1;
    //                    const int length = odes[s].length;
    //                    for (int i = 0; i < length; ++i) {
    //                        if (i%6<3){ // Drift
    //                            y1[i] = y0[i] + 0.5*subStep * y0[i+3];
    //                        }
    //                    }
    //                }
    //                t += 0.5*subStep;
    //                for (int s=0; s < Ns; s++){
    //                    odes[s].derivatives(&odes[s], odes[s].yDot, odes[s].y1, t);
    //                }
    //                for (int s=0; s < Ns; s++){
    //                    double* y0 = odes[s].y;
    //                    double* y1 = odes[s].y1;
    //                    double* yDot = odes[s].yDot;
    //                    const int length = odes[s].length;
    //                    for (int i = 0; i < length; ++i) {
    //                        if (i%6>2){ // Kick
    //                            y1[i] = y0[i] + subStep * yDot[i];
    //                        }
    //                    }
    //                }
    //
    //
    //                // other substeps
    //                for (int j = 1; j < n; ++j) {
    //                    t += subStep;
    //                    for (int s=0; s < Ns; s++){
    //                        double* y1 = odes[s].y1;
    //                        const int length = odes[s].length;
    //                        for (int i = 0; i < length; ++i) {
    //                            if (i%6<3){ // Drift
    //                                y1[i] = y1[i] + subStep * y1[i+3];
    //                            }
    //                        }
    //                    }
    //                    for (int s=0; s < Ns; s++){
    //                        odes[s].derivatives(&odes[s], odes[s].yDot, odes[s].y1, t);
    //                    }
    //                    for (int s=0; s < Ns; s++){
    //                        double* y1 = odes[s].y1;
    //                        double* yDot = odes[s].yDot;
    //                        const int length = odes[s].length;
    //                        for (int i = 0; i < length; ++i) {
    //                            if (i%6>2){ // Kick
    //                                y1[i] = y1[i] + subStep * yDot[i];
    //                            }
    //                        }
    //                    }
    //
    //                    // stability check
    //                    //if (performStabilityCheck && (j <= maxChecks) && (k < maxIter)) {
    //                    //    double initialNorm = 0.0;
    //                    //    for (int l = 0; l < length; ++l) {
    //                    //        const double ratio = y0Dot[l] / scale[l];
    //                    //        initialNorm += ratio * ratio;
    //                    //    }
    //                    //    double deltaNorm = 0.0;
    //                    //    for (int l = 0; l < length; ++l) {
    //                    //        const double ratio = (yDot[l] - y0Dot[l]) / scale[l];
    //                    //        deltaNorm += ratio * ratio;
    //                    //    }
    //                    //    //printf("iii   %e %e\n",initialNorm, deltaNorm);
    //                    //    if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
    //                    //        return 0;
    //                    //    }
    //                    //}
    //                }
    //
    //                // correction of the last substep (at t0 + step)
    //                for (int s=0; s < Ns; s++){
    //                    double* y1 = odes[s].y1;
    //                    const int length = odes[s].length;
    //                    for (int i = 0; i < length; ++i) {
    //                        if (i%6<3){ // Drift
    //                            y1[i] = y1[i] + 0.5 * subStep * y1[i+3];
    //                        }
    //                    }
    //                }
    //
    //                return 1;
    //            }
    
    
    // Modified Midpoint method
    // first substep
    t += subStep;
    for (int s=0; s < Ns; s++){
        double* y0 = odes[s]->y;
        double* y1 = odes[s]->y1;
        double* y0Dot = odes[s]->y0Dot;
        const int length = odes[s]->length;
        for (int i = 0; i < length; ++i) {
            y1[i] = y0[i] + subStep * y0Dot[i];
        }
    }

    // other substeps
    if (needs_nbody){
        reb_integrator_bs_update_particles(r, r->ri_bs.nbody_ode->y1);
    }
    for (int s=0; s < Ns; s++){
        odes[s]->derivatives(odes[s], odes[s]->yDot, odes[s]->y1, t);
    }
    for (int s=0; s < Ns; s++){
        double* y0 = odes[s]->y;
        double* yTmp = odes[s]->yTmp;
        const int length = odes[s]->length;
        for (int i = 0; i < length; ++i) {
            yTmp[i] = y0[i];
        }
    }

    for (int j = 1; j < n; ++j) {  // Note: iterating n substeps, not 2n substeps as in Eq. (9.13)
        t += subStep;
        for (int s=0; s < Ns; s++){
            double* y1 = odes[s]->y1;
            double* yDot = odes[s]->yDot;
            double* yTmp = odes[s]->yTmp;
            const int length = odes[s]->length;
            for (int i = 0; i < length; ++i) {
                const double middle = y1[i];
                y1[i]       = yTmp[i] + 2.* subStep * yDot[i];
                yTmp[i]       = middle;
            }
        }

        if (needs_nbody){
            reb_integrator_bs_update_particles(r, r->ri_bs.nbody_ode->y1);
        }
        for (int s=0; s < Ns; s++){
            odes[s]->derivatives(odes[s], odes[s]->yDot, odes[s]->y1, t);
        }

        // stability check
        if (j <= maxChecks && k < maxIter) {
            double initialNorm = 0.0;
            double deltaNorm = 0.0;
            for (int s=0; s < Ns; s++){
                double* yDot = odes[s]->yDot;
                double* y0Dot = odes[s]->y0Dot;
                double* scale = odes[s]->scale;
                const int length = odes[s]->length;
                for (int l = 0; l < length; ++l) {
                    const double ratio1 = y0Dot[l] / scale[l];
                    initialNorm += ratio1 * ratio1;
                    const double ratio2 = (yDot[l] - y0Dot[l]) / scale[l];
                    deltaNorm += ratio2 * ratio2;
                }
            }
            if (deltaNorm > 4 * MAX(1.0e-15, initialNorm)) {
                return 0;
            }
        }

    }

    // correction of the last substep (at t0 + step)
    for (int s=0; s < Ns; s++){
        double* y1 = odes[s]->y1;
        double* yTmp = odes[s]->yTmp;
        double* yDot = odes[s]->yDot;
        const int length = odes[s]->length;
        for (int i = 0; i < length; ++i) {
            y1[i] = 0.5 * (yTmp[i] + y1[i] + subStep * yDot[i]); 
        }
    }

    return 1;
}

static void extrapolate(const struct reb_ode* ode, double * const coeff, const int k) {
    double* const y1 = ode->y1;
    double* const C = ode->C;
    double** const D =  ode->D;
    double const length = ode->length;
    for (int j = 0; j < k; ++j) {
        double xi = coeff[k-j-1];
        double xim1 = coeff[k];
        double facC = xi/(xi-xim1);
        double facD = xim1/(xi-xim1);
        for (int i = 0; i < length; ++i) {
            double CD = C[i] - D[k - j -1][i];
            C[i] = facC * CD;
            D[k - j - 1][i] = facD * CD;
        }
    }
    for (int i = 0; i < length; ++i) {
        y1[i] = D[0][i];
    }
    for (int j = 1; j <= k; ++j) {
        for (int i = 0; i < length; ++i) {
        y1[i] += D[j][i];
        }
    }
}


static void nbody_derivatives(struct reb_ode* ode, double* const yDot, const double* const y, double const t){
    struct reb_simulation* const r = ode->r;
    if (r->t != t) { 
        // Not needed for first step. Accelerations already calculated. Just need to copy them
        reb_integrator_bs_update_particles(r, y);
        reb_simulation_update_acceleration(r);
    }

    for (unsigned int i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        yDot[i*6+0] = p.vx;
        yDot[i*6+1] = p.vy;
        yDot[i*6+2] = p.vz;
        yDot[i*6+3] = p.ax;
        yDot[i*6+4] = p.ay;
        yDot[i*6+5] = p.az;
    }
}


void reb_integrator_bs_part1(struct reb_simulation* r){
    if (r->calculate_megno){
        reb_simulation_error(r, "The BS integrator does currently not support MEGNO.");
    }
            
    struct reb_ode** odes = r->odes;
    int Ns = r->N_odes;
    for (int s=0; s < Ns; s++){
        const int length = odes[s]->length;
        double* y0 = odes[s]->y;
        double* y1 = odes[s]->y1;
        for (int i = 0; i < length; ++i) {
            y1[i] = y0[i];
        }
    }
}

static void allocate_sequence_arrays(struct reb_integrator_bs* ri_bs){
    ri_bs->sequence        = malloc(sizeof(int)*sequence_length);
    ri_bs->cost_per_step     = malloc(sizeof(int)*sequence_length);
    ri_bs->coeff           = malloc(sizeof(double)*sequence_length);
    ri_bs->cost_per_time_unit = malloc(sizeof(double)*sequence_length);
    ri_bs->optimal_step     = malloc(sizeof(double)*sequence_length);

    // step size sequence: 2, 6, 10, 14, ...  // only needed for dense output
     for (int k = 0; k < sequence_length; ++k) {
        ri_bs->sequence[k] = 4 * k + 2;
    }
    
    // step size sequence: 1,2,3,4,5 ...
    //for (int k = 0; k < sequence_length; ++k) {
    //    ri_bs->sequence[k] = 2*( k+1);
    //}

    // initialize the order selection cost array
    // (number of function calls for each column of the extrapolation table)
    ri_bs->cost_per_step[0] = ri_bs->sequence[0] + 1;
    for (int k = 1; k < sequence_length; ++k) {
        ri_bs->cost_per_step[k] = ri_bs->cost_per_step[k - 1] + ri_bs->sequence[k];
    }
    ri_bs->cost_per_time_unit[0]       = 0;

    // initialize the extrapolation tables
    for (int j = 0; j < sequence_length; ++j) {
        double r = 1./((double) ri_bs->sequence[j]);
        ri_bs->coeff[j] = r*r;
    }
}

static void reb_integrator_bs_default_scale(struct reb_ode* ode, double* y1, double* y2, double relTol, double absTol){
    double* scale = ode->scale;
    int length = ode->length;
    for (int i = 0; i < length; i++) {
        scale[i] = absTol + relTol * MAX(fabs(y1[i]), fabs(y2[i]));
    }
}


int reb_integrator_bs_step(struct reb_simulation* r, double dt){
    // return 1 if step was successful
    //        0 if rejected 
    //
    struct reb_integrator_bs* ri_bs = &r->ri_bs;
    
    if (ri_bs->sequence==NULL){
        allocate_sequence_arrays(ri_bs);
    }

    double t = r->t;
    ri_bs->dt_proposed = dt; // In case of early fail

    // initial order selection
    if (ri_bs->target_iter == 0){
        const double tol    = ri_bs->eps_rel;
        const double log10R = log10(MAX(1.0e-10, tol));
        ri_bs->target_iter = MAX(1, MIN(sequence_length - 2, (int) floor(0.5 - 0.6 * log10R)));
    }

    // maxError not used at the moment.
    // double  maxError = DBL_MAX;

    int Ns = r->N_odes; // Number of ode sets
    struct reb_ode** odes = r->odes;
    double error;
    int reject = 0;

    // Check if ODEs have been set up correctly 
    for (int s=0; s < Ns; s++){
        if (!odes[s]->derivatives){
            reb_simulation_error(r,"A user-specified set of ODEs has not been provided with a derivatives function.");
            r->status = REB_STATUS_GENERIC_ERROR;
            return 0;
        }
    }

    for (int s=0; s < Ns; s++){
        // Check if ODEs need pre timestep setup
        if (odes[s]->pre_timestep){
            odes[s]->pre_timestep(odes[s], odes[s]->y);
        }
        // Scaling
        if (odes[s]->getscale){
            odes[s]->getscale(odes[s], odes[s]->y, odes[s]->y); // initial scaling
        }else{
            reb_integrator_bs_default_scale(odes[s], odes[s]->y, odes[s]->y, ri_bs->eps_rel, ri_bs->eps_abs);
        }
    }

    // first evaluation, at the beginning of the step
    for (int s=0; s < Ns; s++){
        odes[s]->derivatives(odes[s], odes[s]->y0Dot, odes[s]->y, t);
    }

    const int forward = (dt >= 0.);

    // iterate over several substep sizes
    int k = -1;
    for (int loop = 1; loop; ) {

        ++k;
        
        // modified midpoint integration with the current substep
        if ( ! tryStep(r, Ns, k, ri_bs->sequence[k], t, dt)) {

            // the stability check failed, we reduce the global step
#if DEBUG
            printf("S");
#endif
            dt  = fabs(dt * stabilityReduction);
            reject = 1;
            loop   = 0;

        } else {
            for (int s=0; s < Ns; s++){
                const int length = odes[s]->length;
                for (int i = 0; i < length; ++i) {
                    double CD = odes[s]->y1[i];
                    odes[s]->C[i] = CD;
                    odes[s]->D[k][i] = CD;
                }
            }

            // the substep was computed successfully
            if (k > 0) {

                // extrapolate the state at the end of the step
                // using last iteration data
                for (int s=0; s < Ns; s++){
                    extrapolate(odes[s], ri_bs->coeff, k);
                    if (odes[s]->getscale){
                        odes[s]->getscale(odes[s], odes[s]->y, odes[s]->y1);
                    }else{
                        reb_integrator_bs_default_scale(odes[s], odes[s]->y, odes[s]->y, ri_bs->eps_rel, ri_bs->eps_abs);
                    }
                }

                // estimate the error at the end of the step.
                error = 0;
                //long int combined_length = 0;
                for (int s=0; s < Ns; s++){
                    const int length = odes[s]->length;
                    //combined_length += length;
                    double * C = odes[s]->C;
                    double * scale = odes[s]->scale;
                    for (int j = 0; j < length; ++j) {
                        const double e = C[j] / scale[j];
                        error = MAX(error, e * e);
                    }
                }
                // Note: Used to be: error = sqrt(error / combined_length). But for N-body applications it might be more consistent to use:
                error = sqrt(error);
                if (isnan(error)) {
                    reb_simulation_error(r, "NaN appearing during ODE integration.");
                    r->status = REB_STATUS_GENERIC_ERROR;
                    return 0;
                }

                if ((error > 1.0e25)){ // TODO: Think about what to do when error increases: || ((k > 1) && (error > maxError))) {
                    // error is too big, we reduce the global step
#if DEBUG
                    printf("R (error= %.5e)",error);
#endif
                    dt  = fabs(dt * stabilityReduction);
                    reject = 1;
                    loop   = 0;
                } else {

                    // Not used at the moment
                    // maxError = MAX(4 * error, 1.0);

                    // compute optimal stepsize for this order
                    const double exp = 1.0 / (2 * k + 1);
                    double fac = stepControl2 / pow(error / stepControl1, exp);
                    const double power = pow(stepControl3, exp);
                    fac = MAX(power / stepControl4, MIN(1. / power, fac));
                    ri_bs->optimal_step[k]     = fabs(dt * fac);
                    ri_bs->cost_per_time_unit[k] = ri_bs->cost_per_step[k] / ri_bs->optimal_step[k];

                    // check convergence
                    switch (k - ri_bs->target_iter) {

                        case -1 : // one before target
                            if ((ri_bs->target_iter > 1) && ! ri_bs->previous_rejected) {

                                // check if we can stop iterations now
                                if (error <= 1.0) {
                                    // convergence have been reached just before target_iter
                                    loop = 0;
                                } else {
                                    // estimate if there is a chance convergence will
                                    // be reached on next iteration, using the
                                    // asymptotic evolution of error
                                    const double ratio = ((double) ri_bs->sequence[ri_bs->target_iter] * ri_bs->sequence[ri_bs->target_iter + 1]) / (ri_bs->sequence[0] * ri_bs->sequence[0]);
                                    if (error > ratio * ratio) {
                                        // we don't expect to converge on next iteration
                                        // we reject the step immediately and reduce order
                                        reject = 1;
                                        loop   = 0;
                                        ri_bs->target_iter = k;
                                        if ((ri_bs->target_iter > 1) &&
                                                (ri_bs->cost_per_time_unit[ri_bs->target_iter - 1] <
                                                 orderControl1 * ri_bs->cost_per_time_unit[ri_bs->target_iter])) {
                                            ri_bs->target_iter -= 1;
                                        }
                                        dt = ri_bs->optimal_step[ri_bs->target_iter];
#if DEBUG
                                        printf("O");
#endif
                                    }
                                }
                            }
                            break;

                        case 0: // exactly on target
                            if (error <= 1.0) {
                                // convergence has been reached exactly at target_iter
                                loop = 0;
                            } else {
                                // estimate if there is a chance convergence will
                                // be reached on next iteration, using the
                                // asymptotic evolution of error
                                const double ratio = ((double) ri_bs->sequence[k + 1]) / ri_bs->sequence[0];
                                if (error > ratio * ratio) {
                                    // we don't expect to converge on next iteration
                                    // we reject the step immediately
#if DEBUG
                                    printf("o");
#endif
                                    reject = 1;
                                    loop = 0;
                                    if ((ri_bs->target_iter > 1) &&
                                            (ri_bs->cost_per_time_unit[ri_bs->target_iter - 1] <
                                             orderControl1 * ri_bs->cost_per_time_unit[ri_bs->target_iter])) {
                                        --ri_bs->target_iter;
                                    }
                                    dt = ri_bs->optimal_step[ri_bs->target_iter];
                                }
                            }
                            break;

                        case 1 : // one past target
                            if (error > 1.0) {
#if DEBUG
                                printf("e");
#endif
                                reject = 1;
                                if ((ri_bs->target_iter > 1) &&
                                        (ri_bs->cost_per_time_unit[ri_bs->target_iter - 1] <
                                         orderControl1 * ri_bs->cost_per_time_unit[ri_bs->target_iter])) {
                                    --ri_bs->target_iter;
                                }
                                dt = ri_bs->optimal_step[ri_bs->target_iter];
                            }
                            loop = 0;
                            break;

                        default :
                            if (ri_bs->first_or_last_step && (error <= 1.0)) {
                                loop = 0;
                            }
                            break;

                    }
                }
            }
        }
    }


    if (! reject) {
#if DEBUG
        printf("."); 
#endif
        // Swap arrays
        for (int s=0; s < Ns; s++){
            double* y_tmp = odes[s]->y;
            odes[s]->y = odes[s]->y1; 
            odes[s]->y1 = y_tmp; 
            // Check if ODEs need post timestep call
            if (odes[s]->post_timestep){
                odes[s]->post_timestep(odes[s], odes[s]->y);
            }
        }

        int optimalIter;
        if (k == 1) {
            optimalIter = 2;
            if (ri_bs->previous_rejected) {
                optimalIter = 1;
            }
        } else if (k <= ri_bs->target_iter) { // Converged before or on target
            optimalIter = k;
            if (ri_bs->cost_per_time_unit[k - 1] < orderControl1 * ri_bs->cost_per_time_unit[k]) {
                optimalIter = k - 1;
            } else if (ri_bs->cost_per_time_unit[k] < orderControl2 * ri_bs->cost_per_time_unit[k - 1]) {
                optimalIter = MIN(k + 1, sequence_length - 2);
            }
        } else {                            // converged after target
            optimalIter = k - 1;
            if ((k > 2) && (ri_bs->cost_per_time_unit[k - 2] < orderControl1 * ri_bs->cost_per_time_unit[k - 1])) {
                optimalIter = k - 2;
            }
            if (ri_bs->cost_per_time_unit[k] < orderControl2 * ri_bs->cost_per_time_unit[optimalIter]) {
                optimalIter = MIN(k, sequence_length - 2);
            }
        }

        if (ri_bs->previous_rejected) {
            // after a rejected step neither order nor stepsize
            // should increase
            ri_bs->target_iter = MIN(optimalIter, k);
            dt = MIN(fabs(dt), ri_bs->optimal_step[ri_bs->target_iter]);
        } else {
            // stepsize control
            if (optimalIter <= k) {
                dt = ri_bs->optimal_step[optimalIter];
            } else {
                if ((k < ri_bs->target_iter) &&
                        (ri_bs->cost_per_time_unit[k] < orderControl2 * ri_bs->cost_per_time_unit[k - 1])) {
                    dt = ri_bs->optimal_step[k] * ri_bs->cost_per_step[optimalIter + 1] / ri_bs->cost_per_step[k];
                } else {
                    dt = ri_bs->optimal_step[k] * ri_bs->cost_per_step[optimalIter] / ri_bs->cost_per_step[k];
                }
            }

            ri_bs->target_iter = optimalIter;

        }
    }

    dt = fabs(dt);

    if (ri_bs->min_dt !=0.0 && dt < ri_bs->min_dt) {
        dt = ri_bs->min_dt;
        reb_simulation_warning(r,"Minimal stepsize reached during ODE integration.");
    }

    if (ri_bs->max_dt !=0.0 && dt > ri_bs->max_dt) {
        dt = ri_bs->max_dt;
        reb_simulation_warning(r,"Maximum stepsize reached during ODE integration.");
    }

    if (! forward) {
        dt = -dt;
    }
    ri_bs->dt_proposed = dt;

    if (reject) {
        ri_bs->previous_rejected = 1;
    } else {
        ri_bs->previous_rejected = 0;
        ri_bs->first_or_last_step = 0;
    }
    return !reject;
}

struct reb_ode* reb_ode_create(struct reb_simulation* r, unsigned int length){
    struct reb_ode* ode = malloc(sizeof(struct reb_ode));
    
    memset(ode, 0, sizeof(struct reb_ode)); // not really necessaery

    if (r->N_allocated_odes <= r->N_odes){
        r->N_allocated_odes += 32;
        r->odes = realloc(r->odes,sizeof(struct reb_ode*)*r->N_allocated_odes);
    }
    
    r->odes[r->N_odes] = ode;
    r->N_odes += 1;


    ode->r = r; // weak reference
    ode->length = length;
    ode->needs_nbody = 1;
    ode->N_allocated = length;
    ode->getscale = NULL;
    ode->derivatives = NULL;
    ode->pre_timestep = NULL;
    ode->post_timestep = NULL;
    ode->D   = malloc(sizeof(double*)*(sequence_length));
    for (int k = 0; k < sequence_length; ++k) {
        ode->D[k]   = malloc(sizeof(double)*length);
    }

    ode->C     = malloc(sizeof(double)*length);
    ode->y     = malloc(sizeof(double)*length);
    ode->y1    = malloc(sizeof(double)*length);
    ode->y0Dot = malloc(sizeof(double)*length);
    ode->yTmp  = malloc(sizeof(double)*length);
    ode->yDot  = malloc(sizeof(double)*length);

    ode->scale = malloc(sizeof(double)*length);

    r->ri_bs.first_or_last_step = 1;

    return ode;
}

void reb_integrator_bs_part2(struct reb_simulation* r){
    struct reb_integrator_bs* ri_bs = &(r->ri_bs);
    
    unsigned int nbody_length = r->N*3*2;
    // Check if particle numbers changed, if so delete and recreate ode.
    if (ri_bs->nbody_ode != NULL){ 
        if (ri_bs->nbody_ode->length != nbody_length){
            reb_ode_free(ri_bs->nbody_ode);
            ri_bs->nbody_ode = NULL;
        }
    }
    if (ri_bs->nbody_ode == NULL){ 
        ri_bs->nbody_ode = reb_ode_create(r, nbody_length);
        ri_bs->nbody_ode->derivatives = nbody_derivatives;
        ri_bs->nbody_ode->needs_nbody = 0; // No need to update unless there's another ode
        ri_bs->first_or_last_step = 1;
    }
    
    for (int s=0; s < r->N_odes; s++){
        if (r->odes[s]->needs_nbody){
            ri_bs->user_ode_needs_nbody = 1;
        }
    }

    double* const y = ri_bs->nbody_ode->y;
    for (unsigned int i=0; i<r->N; i++){
        const struct reb_particle p = r->particles[i];
        y[i*6+0] = p.x;
        y[i*6+1] = p.y;
        y[i*6+2] = p.z;
        y[i*6+3] = p.vx;
        y[i*6+4] = p.vy;
        y[i*6+5] = p.vz;
    }

    int success = reb_integrator_bs_step(r, r->dt);
    if (success){
        r->t += r->dt;
        r->dt_last_done = r->dt;
    }
    r->dt = ri_bs->dt_proposed;

    reb_integrator_bs_update_particles(r, ri_bs->nbody_ode->y);
}

void reb_integrator_bs_synchronize(struct reb_simulation* r){
    // Do nothing.
}

void reb_ode_free(struct reb_ode* ode){
    // Free data array
    free(ode->y);
    ode->y = NULL;
    free(ode->y1);
    ode->y1 = NULL;
    free(ode->C);
    ode->C = NULL;
    free(ode->scale);
    ode->scale = NULL;
    
    if (ode->D){
        for (int k=0; k < sequence_length; k++) {
            free(ode->D[k]);
        }
        free(ode->D);
        ode->D = NULL;
    }
    free(ode->y0Dot);
    ode->y0Dot = NULL;
    free(ode->yTmp);
    ode->yTmp = NULL;
    free(ode->yDot);
    ode->yDot = NULL;
    
    struct reb_simulation* r = ode->r;
    if (r){ // only do this is ode is in a simulation
        struct reb_integrator_bs* ri_bs = &r->ri_bs;
        int shift = 0;
        for (int s=0; s < r->N_odes; s++){
            if (r->odes[s] == ode){
                r->N_odes--;
                shift = 1;
            }
            if (shift && s <= r->N_odes ){
                r->odes[s] = r->odes[s+1];
            }
        }
        if (ri_bs->nbody_ode == ode){
            ri_bs->nbody_ode = NULL;
        }
    }
    free(ode);
}



void reb_integrator_bs_reset(struct reb_simulation* r){
    struct reb_integrator_bs* ri_bs = &(r->ri_bs);
    
    // Delete nbody ode but not others
    if (ri_bs->nbody_ode){
        reb_ode_free(ri_bs->nbody_ode);
        ri_bs->nbody_ode = NULL;
    }

    // Free sequence arrays
    free(ri_bs->sequence);
    ri_bs->sequence = NULL;
    
    free(ri_bs->coeff);
    ri_bs->coeff = NULL;
    free(ri_bs->cost_per_step);
    ri_bs->cost_per_step = NULL;
    free(ri_bs->cost_per_time_unit);
    ri_bs->cost_per_time_unit = NULL;
    free(ri_bs->optimal_step);
    ri_bs->optimal_step = NULL;
    
    
    // Default settings
    ri_bs->eps_abs          = 1e-8;
    ri_bs->eps_rel          = 1e-8;
    ri_bs->max_dt           = 0;
    ri_bs->min_dt           = 0; 
    ri_bs->first_or_last_step  = 1;
    ri_bs->previous_rejected = 0;
    ri_bs->target_iter       = 0;
        
}

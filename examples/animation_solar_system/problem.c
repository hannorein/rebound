/**
 * Animation of the Solar System
 * 
 * This examples show how to use display_settings to
 * programmatically change the visualization of a 
 * REBOUND simulation. This can be used to render movies.
 * To understand what a 4x4 view matrix is, you can read 
 * up on linear algebra for computer graphics, specifically
 * the Model-View-Projection (MVP) paradigm. 
 *
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

struct reb_mat4df view0; // Initial view matrix
double dt0;              // Initial timestep

double inflate_size = 5; // Inflate particle sizes

// Quadratic ease in/out function for smooth animations
double ease_in_out(double x){
    if (x<0.0) x=0.0;
    if (x>1.0) x=1.0;
    return x < 0.5 ? 2.0*x*x : 1.0-(2.0-2.0* x)*(1.0-1.0*x);
}

// Heartbeat function changes the display settings every timestep
// The following animations play in this order:
// - Zoom in
// - Rotation around the x axis
// - Slowing down of simulation
// - Combined zoom in on Earth and rotation 
void heartbeat(struct reb_simulation* const r){
    // Construct a 90 rotation around the x axis
    struct reb_vec3d a = {.x=1.};                               // x axis
    struct reb_rotation rot_x = reb_rotation_init_angle_axis(M_PI_2, a); // pi/2 = 90 degrees
        
    // We start from the original view matrix view0, but we could also 
    // apply consequitive changes to the current view matrix stored in 
    // r->display_settings->view
    struct reb_mat4df view = view0;
                                                                     
    if (r->t < 2.*2.*M_PI){                                         // First 2 years

        // Slerp (interpolate) between no rotation (identity) and a rotation around the x axis
        float t = ease_in_out(r->t/(2.*2.*M_PI));                   // runs from 0 to 1
        struct reb_rotation rot_slerp = reb_rotation_slerp(reb_rotation_identity(), rot_x, t);

        // Convert rotation quaternion to 4d rotation matrix and then
        // operate the rotation matrix on the view matrix.
        struct reb_mat4df rm = reb_rotation_to_mat4df(rot_slerp);
        view = reb_mat4df_multiply(rm, view); 

    }else if (r->t > 2.*2.*M_PI && r->t < 4.*2.*M_PI){              // Year 2 to 4
        // Show orbits as wires
        r->display_settings->wire = 1;
        
        // Increase length of trail to 64
        r->display_settings->breadcrumbs = 64;
                                                       
        // Apply same rotation matrix as before
        struct reb_mat4df rm = reb_rotation_to_mat4df(rot_x);
        view = reb_mat4df_multiply(rm, view); 

        // Also apply a zoom operation
        float t = ease_in_out(((r->t/(2.*M_PI) - 2.0)/2.0));        // runs from 0 to 1
        float s = 1.+30.0*t;                                        // zoom factor 
        view = reb_mat4df_scale(view, s,s,s);                       // zoom in
                                      
    }else if (r->t > 4.*2.*M_PI && r->t < 5.*2.*M_PI){             // Year 4 to 5
    
        
        // Reduce timestep
        r->dt = dt0/10.0; 
        // Hide wires
        r->display_settings->wire = 0;
        if (r->N==9){
            // Add moon (these are not exact parameters, just for illustration)
            struct reb_particle e = r->particles[3];
            reb_simulation_add_fmt(r, "a m r primary", 0.0025695553, 3.6943033e-08, 1.1617812e-05*inflate_size, e);
        }

        // Get Earth-Moon Barycenter
        struct reb_particle em_com = reb_particle_com_of_pair(r->particles[3], r->particles[9]); 
        
        // Apply same rotation and zoom as before (we could cache this)
        struct reb_mat4df rm = reb_rotation_to_mat4df(rot_x);
        view = reb_mat4df_multiply(rm, view); 
        view = reb_mat4df_scale(view, 31.,31.,31.);

        // slowly zoom in
        float ts = ease_in_out((r->t/(2.*M_PI) - 4.0));     // runs from 0 to 1
        float s = 1.+500.0*ts*ts;                                   // zoom factor 
        view = reb_mat4df_scale(view, s,s,s);                   // second zoom operation
        
        // slowly move earth to center
        float tm = ease_in_out((r->t/(2.*M_PI) - 4.0)*3.);    // runs from 0 to 1
        view = reb_mat4df_translate(view, -tm*em_com.x, -tm*em_com.y, -tm*em_com.z);    // translate view
                                                                    
    }else{  // Continue until forever
        
        // Get Earth-Moon Barycenter
        struct reb_particle em_com = reb_particle_com_of_pair(r->particles[3], r->particles[9]); 
        
        // Apply same rotation, zoom, and keep centered on earth
        struct reb_mat4df rm = reb_rotation_to_mat4df(rot_x);
        view = reb_mat4df_multiply(rm, view); 
        view = reb_mat4df_scale(view, 31.,31.,31.);             // We could combine the zoom operations
        view = reb_mat4df_scale(view, 500.,500.,500.);
        view = reb_mat4df_translate(view, -em_com.x, -em_com.y, -em_com.z);    // translate view

        // Show real size of particles
        r->display_settings->spheres = 1;
        // Hide planet trails
        r->display_settings->breadcrumbs = 0;

    }

    // Store the view matrix
    r->display_settings->view = view;
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_simulation_create();

    // Solar System initial conditions from NASA Horizons
    // Units are solar mass, AU, years/2pi.
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 0.9999999999950272, 0.0046524726,
            -0.007784066163300598, -0.003160235321847074, 0.0002085198739177632, 
            0.00029808520873668176, -0.0003933583120820104, -3.0307129383172144e-06);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 1.6601208254808336e-07, 1.6313735e-05,
            -0.0197573020117554, -0.4659490539727674, -0.036512748618975056, 
            1.3072205503890428, 0.04106840511679157, -0.11649038132810724);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 2.447838287784771e-06, 4.0453784e-05,
            -0.31310669579085715, -0.6609240473612075, 0.008792509815427053, 
            1.058788984690766, -0.5006529730168531, -0.06794916369487035);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 3.0404326489511185e-06, 4.2587571e-05,
            -0.7304336184787301, 0.6677833148611707, 0.0001754893039781267, 
            -0.6964579069368259, -0.7370893850360778, 4.193835863036079e-05);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 3.2271560828978514e-07, 2.2702195e-05,
            0.22987379117805698, -1.4195661773490975, -0.035304801938727, 
            0.833220942267478, 0.2040504545094634, -0.01614911011819359);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 0.0009547919099366768, 0.0004778945,
            3.278533805383686, 3.7537759447476775, -0.08892301493418556, 
            -0.3352619068991809, 0.3093717438204439, 0.006217554149365762);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 0.0002858856700231729, 0.0004028667,
            9.053130670591704, -3.529112154855013, -0.2990863565225466, 
            0.09969170490350149, 0.30152611855823136, -0.009212249609912534);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 4.366249613200406e-05, 0.00017085136,
            12.148791396671397, 15.380702941916812, -0.10026566241603418, 
            -0.18110086740935336, 0.1310657378884879, 0.00283290820462697);
    reb_simulation_add_fmt(r, "m r x y z vx vy vz", 5.151383772628957e-05, 0.00016553712,
            29.840954544603573, -1.6778892885804566, -0.6531621601938832, 
            0.00903845442721642, 0.18328001746914682, -0.003982606081774178);

    reb_simulation_move_to_com(r);

    // Inflate sized for illustration
    for(int i=0;i<r->N;i++){
        r->particles[i].r *= inflate_size;
    }
    
    // We use the WHFast integrator and a fixed timestep of 2 days
    r->integrator = REB_INTEGRATOR_WHFAST;
    r->dt = 2./365.25 *2.0*M_PI;

   
    // Starting the REBOUND visualization server. This
    // allows you to visualize the simulation by pointing 
    // your web browser to http://localhost:1234
    reb_simulation_start_server(r, 1234);
    
    // Normally the visualization settings are determined by the 
    // user interface. If we add the display_settings struct to 
    // the simulation itself, it will overwrite any change the 
    // user has made and allows us to programatically change any 
    // settings such as the orientation, zoom, etc. 
    reb_simulation_add_display_settings(r);

    // Initially, rotate to an edge on view
    struct reb_vec3d a = {.x=1.};                               // x axis
    struct reb_rotation rot_x = reb_rotation_init_angle_axis(M_PI_2, a); // pi/2 = 90 degrees
    struct reb_mat4df rm = reb_rotation_to_mat4df(rot_x);
    r->display_settings->view = reb_mat4df_multiply(rm, r->display_settings->view);

    // Store initial view matrix and timestep so we can use it as a reference
    view0 = r->display_settings->view;
    dt0 = r->dt;
    
    // The heartbeat function is called once per timestep and handles the 
    // scripted view changes
    r->heartbeat = heartbeat;

    // Show particles as points (not as spheres with their real size)
    r->display_settings->spheres = 0;
    // Show orbits as planes
    r->display_settings->wire = 2;
    // Show 4 past particle positions
    r->display_settings->breadcrumbs = 4;

    // Slow the simulation down to at most 120 timesteps per second.
    r->usleep = 8333;
    
    // Then keep running forever.
    reb_simulation_integrate(r, INFINITY);

    // Cleanup 
    reb_simulation_free(r);
}


// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"


void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
  
    double K_H, K_P, error, P_out, delta;
    double tot_lander_mass = UNLOADED_LANDER_MASS + FUEL_DENSITY * FUEL_CAPACITY * fuel;


    K_H = 0.018;
    K_P = 5;
    delta = 1.0;

    error = -(K_H * (position.abs()-MARS_RADIUS) + 0.5 + velocity * position.norm());
    P_out = K_P * error;
    
    if (P_out <= -delta) {throttle = 0;
    }
    else if (-delta < P_out < 1 - delta) { throttle = delta + P_out; 
    }
    else { throttle = 1;
    }
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pos. The time step is delta_t (global variable).
{
  
    vector3d thrust, drag, weight, new_pos, acc;
    double density,tot_lander_mass,areaDrag;

    tot_lander_mass = UNLOADED_LANDER_MASS + FUEL_DENSITY*FUEL_CAPACITY*fuel;
    density = atmospheric_density(position);

    double landerDrag = DRAG_COEF_LANDER * M_PI * LANDER_SIZE * LANDER_SIZE;
    double chuteDrag = DRAG_COEF_CHUTE * LANDER_SIZE * LANDER_SIZE * 20; 
    

   
    if (parachute_status == DEPLOYED) {
        areaDrag = chuteDrag + landerDrag;
    }
    else { areaDrag = landerDrag;
    }
    

    thrust = thrust_wrt_world();
    drag = -0.5 * density * areaDrag * velocity.abs2() * velocity.norm();
    weight = -GRAVITY * MARS_MASS * tot_lander_mass * position.norm() / position.abs2();

    acc = (thrust + drag + weight) / tot_lander_mass;


    //select verlet or euler integrator
    if (do_verlet == true) {
        
        //verlet integration
        if (simulation_time == 0) { new_pos = position; }
        else if (simulation_time == delta_t) {
            new_pos = position + velocity * delta_t;
        }
        else {
            new_pos = 2 * position - previous_pos + delta_t * delta_t * acc;
            velocity = 1 / delta_t * (new_pos - position);
        }
        previous_pos = position;
        position = new_pos;
    
    }else { 
        
        //Euler integration
        position += velocity * delta_t;
        velocity += acc * delta_t;
    }




  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();


 
}

void initialize_simulation(void)
{
//reset some global variables for numerical dynamics
 new_pos = previous_pos = vector3d(0, 0, 0);



  // Lander pose initialization - selects one of 10 possible scenarios

  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    break;

  case 6:
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}

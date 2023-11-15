/*
This implements verlet velocity algorithm, to integrate to EOM of 4 planets. 
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <valarray>
using namespace std;
typedef valarray<double> Array;


Array acceleration ( Array x1_vec, Array x2_vec , double m1, double m2 ) {
    Array r = {x2_vec[0] - x1_vec[0],x2_vec[1] - x1_vec[1]}; // displacement vector
    double dist_sqrd = (x2_vec[0] - x1_vec[0]) * (x2_vec[0] - x1_vec[0]) + (x2_vec[1] - x1_vec[1]) * (x2_vec[1] - x1_vec[1]);

    double acc_mag =  (m2) / dist_sqrd;
    Array accel = { r[0] * m2 / (sqrt(dist_sqrd)*dist_sqrd) ,  r[1] * m2 / (sqrt(dist_sqrd) * dist_sqrd) };


    return accel;
}

valarray<Array> Total_acc (valarray<Array> pos, Array m) {

    valarray<Array> acc = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}; // reset

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i != j) {
                 acc[i] += acceleration ( pos[i], pos[j], m[i], m[j] ); //sum gravitys acceleration from each planet
            }
        }
    }
    return acc;
}

//Verlet leap-frog, discussion in notes
void Verlet_lpfg (double dt, valarray<Array>& pos, valarray<Array>& vel, Array m) {

    valarray<Array> old_acc = Total_acc(pos, m);

    for (int i = 0; i < 4; i++) {
        pos[i] = pos[i] + dt * vel[i] + 0.5 * dt*dt * old_acc[i];
    }

    valarray<Array> new_acc = Total_acc(pos, m);

    for (int i = 0; i < 4; i++) {
        vel[i] = vel[i] + 0.5 * (old_acc[i] + new_acc[i]) * dt;
    }

}


int main() {

    valarray<Array> pos = {{-0.5, 0.1}, {-0.6, -0.2}, {0.5, 0.1}, {0.5, 0.4}};
    valarray<Array> vel = {{-0.84, 0.65}, {1.86, 0.7}, {-0.44, -1.5}, {1.15, -1.6}};
    const Array m = {2.2, 0.8, 0.9, 0.4};
    //valarray<Array> acc = Total_acc(pos, m);

    
 
    
    double dt = 0.000001;
    double t = 0;
    int n = 5.0 / dt;
    int i = 0;
    while (i < n) {
        Verlet_lpfg (dt, pos, vel, m);
        t += dt;
        i++;
        //cout << pos[3][0] << setw(10) << pos[3][1] << "\n"; // outputs planet;s trajectory
    }
    for (int i = 0; i < 4; i++) {
        cout << i << " Planet Coordinates and Vel: " << pos[i][0] << setw(10) << pos[i][1] << setw(10) << vel[i][0] << setw(10) << vel[i][1] << "\n"; // outputs final coord/vel

    }
    
    return 0;
}


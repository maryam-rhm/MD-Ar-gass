#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <thread>

const int N = 16;
const double m = 39.948/1000.0/(6.022e23); //kg
const double e = 0.24;
const double s = 3.4;
const double T = 300.0; // K
const double NA = 6.022e23; //Avogadro number
const double kB = 1.38e-23; //J/K
const double dt = 0.1; //s
const double box_size = 20; //A
const double unit = 1.0/(1000.0*NA); // J to KJ/mol


double pos[N][2];
double v[N][2];
double acc[N][2];

void setUp() {
  pos[0][0] = 18.21;
  pos[0][1] = 6.69;

  pos[1][0] = 18.98;
  pos[1][1] = 18.61;

  pos[2][0] = 2.11;
  pos[2][1] = 13.90;
  
  pos[3][0] = 4.85;
  pos[3][1] = 3.55;

  pos[4][0] = 6.14;
  pos[4][1] = 10.61;

  pos[5][0] = 12.09;
  pos[5][1] = 9.97;
  
  pos[6][0] = 11.33;
  pos[6][1] = 3.48;

  pos[7][0] = 12.71;
  pos[7][1] = 17.73;

  pos[8][0] = 10.21;
  pos[8][1] = 14.26;
  
  pos[9][0] = 15.66;
  pos[9][1] = 14.74;

  pos[10][0] = 0.05;
  pos[10][1] = 2.74;

  pos[11][0] = 6.58;
  pos[11][1] = 17.66;
  
  pos[12][0] = 14.70;
  pos[12][1] = 1.21;

  pos[13][0] = 8.20;
  pos[13][1] = 6.25;

  pos[14][0] = 19.03;
  pos[14][1] = 11.29;
  
  pos[15][0] = 2.06;
  pos[15][1] = 7.83;

  v[0][0] = 0.26;
  v[0][1] = 1.28;

  v[1][0] = -1.55;
  v[1][1] = 0.29;

  v[2][0] = -0.72;
  v[2][1] = 4.23;
  
  v[3][0] = 0.65;
  v[3][1] = -1.69;

  v[4][0] = 0.86;
  v[4][1] = -2.02;

  v[5][0] = 1.28;
  v[5][1] = -2.02;

  v[6][0] = -2.05;
  v[6][1] = -2.04;

  v[7][0] = -3.49;
  v[7][1] = -1.92;

  v[8][0] = 2.56;
  v[8][1] = 3.81;

  v[9][0] = 2.27;
  v[9][1] = 2.89;

  v[10][0] = 3.12;
  v[10][1] = 3.02;

  v[11][0] = 3.30;
  v[11][1] = -0.58;

  v[12][0] = -2.69;
  v[12][1] = -3.48;

  v[13][0] = 1.05;
  v[13][1] = -3.91;

  v[14][0] = -3.48;
  v[14][1] = -0.97;

  v[15][0] = -1.34;
  v[15][1] = 3.22;
  
  for (int i=0; i<N; ++i){
     acc[N][0] = 0.0;
     acc[N][1] = 0.0;
  }
}

double potential(double r) {
  return 4*e * (pow((s/r), 12) - pow((s/r), 6));
}

double Fpotential(double r){
  return 4*e * ((12 * pow(s, 12) / pow(r,13)) - (6 * pow(s, 6)/ pow(r, 7)));
}

void Acceleration() {
    for (int i = 0; i<N; ++i) {
        acc[i][0] = 0.0;
        acc[i][1] = 0.0;
    }
    for (int i = 0; i<N-1; ++i) {
        for (int j=i+1; j<N; ++j) {
            double dx = pos[i][0] - pos[j][0];
            double dy = pos[i][1] - pos[j][1];
            const double distance = sqrt(dx*dx + dy*dy);
            double f_magnitude = Fpotential(distance);
            double fx = f_magnitude*dx / distance;
            double fy = f_magnitude*dy / distance;
        
            acc[i][0] += fx;
            acc[i][1] += fy;
            acc[j][0] -= fx;
            acc[j][1] -= fy;
        }
    }
}

double PotentialE(){
    double t_p_e = 0.0; //total potential energy
    for (int i = 0; i<N-1; ++i){
       for (int j=i+1; j<N; ++j){
           double dx = pos[i][0] - pos [j][0];
           double dy = pos[i][1] - pos [j][1];
           const double distance = sqrt(dx*dx + dy*dy);
           t_p_e += potential(distance);
       }
    }
    return t_p_e * unit; 
}

double KineticE(){
    double t_k_e = 0.0;
    for (int i=0; i<N; ++i){
        double speedSq = v[i][0]*v[i][0] + v[i][1]*v[i][1];
        t_k_e += 0.5 * m * speedSq; // J
    }
    return t_k_e * unit;
}

double calculateT () {
    double kinetic = KineticE(); //J
    double temp = (2.0/3.0) * (kinetic / (N * kB)); //K
    return temp;
} 

void thermostat() {
    double currentT = calculateT();
    double l = sqrt(1+ dt * (T / currentT -1));
    for (int i=0; i<N; ++i) {
        v[i][0] *= l;
        v[i][1] *= l;
    }
}

void applyRB(){ //reflective boundries
    for (int i = 0; i < N; ++i){ //x
        if (pos[i][0] < 0){
            pos[i][0] = -pos[i][0];
            v[i][0] = -v[i][0];
        } else if (pos[i][0] > box_size){
            pos[i][0] = 2 * box_size - pos[i][0];
            v[i][0] = -v[i][0]; 
        }
        
        if (pos[i][1] < 0){ //y
            pos[i][1] = -pos[i][1];
            v[i][1] = -v[i][1];
        } else if (pos[i][1] > box_size){
            pos[i][1] = 2 * box_size - pos[i][1];
            v[i][1] = -v[i][1]; 
        }
    }
}


void updateBody() {
  const double dt = 0.0001;
  //update position
  for(int i = 0; i < N; ++i){
      pos[i][0] += v[i][0] * dt + 0.5 * acc[i][0] * dt * dt;
      pos[i][1] += v[i][1] * dt + 0.5 * acc[i][1] * dt * dt;
  }
  //reflective boundries
  applyRB();
  
  //new accelerations    
  double N_acc[N][2];
  Acceleration();
  for (int i=0; i<N; ++i) {
      N_acc[i][0] = acc [i][0];
      N_acc[i][1] = acc [i][1];
  }
  //update velocity
  for(int i=0; i<N; ++i){
      v[i][0] += 0.5*(acc[i][0] + N_acc[i][0])*dt;
      v[i][1] += 0.5*(acc[i][1] + N_acc[i][1])*dt;
  }
  //update acc
  for (int i=0; i<N; ++i){
       acc[i][0] = N_acc[i][0];
       acc[i][1] = N_acc[i][1];
  }
}

void printCSV(int step) {
  std::ofstream out ("result.csv", std::ios::app);
  double t_p_e = PotentialE();
  double t_k_e = KineticE();
  double TotalE = t_p_e + t_k_e;
  double temp = calculateT();
  for (int i=0; i<N; ++i) {
    out << step << ", "
        << pos[i][0]
        << ", "
        << pos[i][1]
        << ", "
        << t_p_e
        << ", "
        << t_k_e
        << ", "
        << TotalE
        << ", "
        << temp
        << "\n";
  }
}

int main() {
    setUp();
    std::ofstream outfile("result.csv");
    outfile << "Step, x, y, Potential Energy(kJ/mol), Kinetic Energy(kJ/mol), Total Energy(kJ/mol), Temperture(K)\n";
    outfile.close();
    
    const int timesteps = 10000;
    const int plotEveryKthStep = 100;
    for (int i=0; i< timesteps; ++i) {
         updateBody();
         thermostat();
         if ((i+1) % plotEveryKthStep ==0) {
             printCSV(i+1); // Please switch off all IO if you do performance tests.
         }      
    }
    
    return 0;
}

#include <iostream>
#include <fstream> // For file handling
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip> // for fixed and setprecision

using namespace std;

// Constants
const double tmax = 10.0; // Maximum simulation time
const double dt = 0.001;   // Time step
const int N = 108;         // Number of particles
const double density = 0.8442; // Desired density
const double temp = 0.728; // Desired temperature (in Kelvin)
const double length = pow((double)N / density, 1.0/3.0); // Length of the simulation box
long idum;
// Particle data
double x[N][3];  // Positions
double v[N][3];  // Velocities
double xm[N][3]; // Previous positions
double f[N][3];  // Force

// Energy data
double kinetic_energy = 0.0;
double potential_energy = 0.0;
double total_energy = 0.0;

// Function prototypes

//Random number generator
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}


// Initialize positions and velocities
void init() {
    double sumv[3] = {0.0, 0.0, 0.0};
    double sumv2 = 0.0;

    // Seed the random number generator
    srand(time(0));

    // Initialize lattice spacing
    int n_per_side = (int)length;  // Number of particles per side of the cubic lattice
    double spacing = 1.0;

    // Loop to initialize positions and velocities
    int index = 0;
    for (int i = 0; i < n_per_side; i++) {
        for (int j = 0; j < n_per_side; j++) {
            for (int k = 0; k < n_per_side; k++) {
                if (index < N) {
                    // Set positions
                    x[index][0] = i * spacing;
                    x[index][1] = j * spacing;
                    x[index][2] = k * spacing;

                    // Set random velocities
                    for (int d = 0; d < 3; d++) {
                        v[index][d] = (ran2(&idum) - 0.5);
                        sumv[d] += v[index][d];
                        sumv2 += v[index][d] * v[index][d];
                    }
                    index++;
                }
            }
        }
    }

    // Calculate mean velocity of center of mass (COM)
    for (int d = 0; d < 3; d++) {
        sumv[d] /= (double)N;
    }
    // Verify line below????????????
    sumv2 /= (double)(3*N-3);

    double fs = sqrt(3.0 * temp / sumv2); // Scale factor of the velocities

    // Adjust velocities to achieve desired temperature
    for (int i = 0; i < N; i++) {
        for (int d = 0; d < 3; d++) {
            v[i][d] = (v[i][d] - sumv[d]) * fs;
            xm[i][d] = x[i][d] - v[i][d] * dt; // Previous position initialization
        }
    }
}

// Calculate forces between particles
double Force() {
    double en = 0; // Energy
    double epsilon = 1.0;//changed from 0.4 to 1.0
    double sigma = 1.0;
    double sigma6 = pow(sigma, 6);
    double rc2 = pow((2.5 * sigma), 2); // Cutoff radius squared
    double ecut = 4 * epsilon * ((sigma6 * sigma6 / pow(rc2, 6)) - (sigma6 / pow(rc2, 3))); // Cut-off energy

    // Initialize forces to zero
    for (int i = 0; i < N; i++) {
        for (int d = 0; d < 3; d++) {
            f[i][d] = 0.0;
        }
    }

    // Loop over all pairs of particles
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double xr[3], r2 = 0.0;

            // Calculate distance between particles i and j in each dimension
            for (int d = 0; d < 3; d++) {
                xr[d] = x[i][d] - x[j][d];
                xr[d] -= length * round(xr[d] / length); // Apply periodic boundary condition
                r2 += xr[d] * xr[d];
            }

            // Check if the distance is within the cutoff radius
            if (r2 < rc2) {
                double r2i = 1.0 / r2;
                double r6i = pow(r2i, 3);
                double ff = 48 * r2i * r6i * (r6i - 0.5); // Lennard-Jones potential

                // Calculate force and update forces for particles i and j
                for (int d = 0; d < 3; d++) {
                    double force = ff * xr[d];
                    f[i][d] += force;
                    f[j][d] -= force;
                }

                // Update potential energy
                en += 4 * epsilon * sigma6 * r6i * (sigma6 * r6i - 1) - ecut;
            }
        }
    }
    return en;
}

// Integrate using the velocity-Verlet algorithm
void integrate(double en) {
    double sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        double xx[3], vi[3];
        for (int d = 0; d < 3; d++) {
            xx[d] = 2.0 * x[i][d] - xm[i][d] + dt * dt * f[i][d];
            vi[d] = (xx[d] - xm[i][d]) / (2.0 * dt);

            // Apply periodic boundary conditions to ensure particles stay within the box
            /*if (xx[d] < 0) {
                xx[d] += length;
            }
            if (xx[d] >= length) {
                xx[d] -= length;
            }*/  //Caused error

            xm[i][d] = x[i][d];
            x[i][d] = xx[d];
            sumv2 += vi[d] * vi[d];
        }
    }

    double ntemp = sumv2 / (3.0 * (double)N);
    kinetic_energy = 0.5 * sumv2 / (double)N; // Kinetic energy
    potential_energy = en / (double)N;
    total_energy = potential_energy + kinetic_energy;
}

// Write function to display kinetic, potential and total energy
/*void write() {
    // Ensure the energies are not NaN or Inf
    if (isnan(potential_energy) || isnan(kinetic_energy) || isnan(total_energy) ||
        isinf(potential_energy) || isinf(kinetic_energy) || isinf(total_energy)) {
        cout << endl;
    } else {
        cout << fixed << setprecision(4) << potential_energy << "\t\t" << kinetic_energy << "\t\t" << total_energy << endl;
    }
}*/
// Print the positions of all particles
void write_positions() {
    for (int i = 0; i < N; ++i) {
        cout << "0" << "\t";
        for (int d = 0; d < 3; ++d) {
            cout << fixed << setprecision(4) << x[i][d] << "\t\t";
        }
        cout << endl;
    }
}
// Print the velocities of all particles
void write_velocities() {
    for (int i = 0; i < N; ++i) {
        cout << "0" << "\t";
        for (int d = 0; d < 3; ++d) {
            cout << fixed << setprecision(4) << v[i][d] << "\t\t";
        }
        cout << endl;
    }
}
void write(){
    cout << fixed << setprecision(3) << potential_energy << "\t\t" << kinetic_energy << "\t\t" << total_energy << endl;
}
/*
int main() {

    idum=-12412424+rand()%100;
    init();
    double t = 0;
    int k = 0;



    while (t < tmax) {
        double e = Force();
        integrate(e);
        cout << k << "\t";
        k++;
        write();
        //write_velocities();
        t += dt;
    }
    return 0;
}
*/


int main() {
    idum = -12412424 + rand() % 100;
    init();
    double t = 0;
    int k = 0;

    // Open file to write energies
    ofstream outfile("energy_data.txt");

    // Check if file is opened successfully
    if (!outfile.is_open()) {
        cerr << "Error opening file for writing." << endl;
        return 1;
    }

    // Write header to the file
    //outfile << "Step\tPotential Energy\tKinetic Energy\tTotal Energy" << endl;

    while (t < tmax) {
        double e = Force();
        integrate(e);

        // Write data to the file
        outfile << k << "\t" << fixed << setprecision(3) << potential_energy << "\t\t" << kinetic_energy << "\t\t" << total_energy << endl;
        k++;
        t += dt;
    }

    // Close the file
    outfile.close();

    return 0;
}

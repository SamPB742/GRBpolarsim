#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

//index into a square 3d matrix
#define VEC_FLD_IDX(x,y,z,s) ((int)((x*s*s) + (y*s) + z))

typedef struct mode {
    double pol;
    double k;
    double phase;
    double theta;
    double phi;
} mode;
typedef struct vector {
    double x;
    double y;
    double z;
} vector;

double *gen1dturb(double sigma, double corr_len, double B0, int iters, double step);
double genRandDouble(double start, double stop);
struct mode* genModes(double min_k, double max_k, int num, int dim);
double G(mode mode1, double corr_len, int dim);
vector *gen3dturb(double sigma, double corr_len, vector B0, int iters, double step, int dim);
vector modeCoords(vector v, mode m);
vector addVectors(vector v, vector w);
vector multVector(double s, vector v);

//
// Created by 4jung on 11/6/2023.
//
int main(int argc, char **argv) {
    //for (int i=0; i<argc; i++)
        //printf("%s\n", argv[i]);
    double sigma = strtod(argv[1], NULL);
    double B0 = strtod(argv[2], NULL);
    vector B0v = {0, 0, B0};
    //printf("vals: %f\n%f\n ", sigma, B0);
    //seed the random generator
    srand(time(NULL));
    double *turb = gen3dturb(sigma, 1, B0v, 1000, .1, 3);
    //TODO figure out how to graph the result
    for (int i = 0; i < 1000; i++) {
        printf("%f\n", turb[i]);
    }
    free(turb);
    return 0;
}


/*
 * Generate 1d turbulence and return an array of sampled values from it
 * of length iters with step size step. The magnetic turbulence variance
 * sigma, the correlation length corr_len, and the background magnetic
 * field strength B0 are all included as parameters.
 *
 * Returns: a pointer to an array of length iters containing the field strength
 * at the sampled points. Each of these is the background field plus the
 * x-component of the calculated dB for that step. Returns NULL if anything goes wrong.
 *TODO add in y component
 *
double *gen1dturb(double sigma, double corr_len, double B0, int iters, double step){
    //TODO allow k generation to be more flexible
    //generate 10 modes with wavenumbers 1-3 (as in Baring turb)
    mode *lower_modes = genModes(1, 3, 10, 1);
    //generate 10 modes with wavenumbers 3-10 (as in Baring turb)
    mode *higher_modes = genModes(3, 10, 10, 1);
    if (lower_modes == NULL || higher_modes == NULL) {
        printf("Mode generation failed");
        return NULL;
    }

    mode *all_modes = malloc(20 * sizeof(*all_modes));
    memcpy(all_modes     , lower_modes, 10 * sizeof(*all_modes));
    memcpy(all_modes + 10, higher_modes, 10 * sizeof(*all_modes));
    free(lower_modes);
    free(higher_modes);

    //calculate sum of G for all modes first
    double G_sum = 0;
    for (int i = 0; i < 10; i++) {
        G_sum += G(all_modes[i], corr_len, 1);
    }

    //go through each of the (x,y,z) points in a cube
    //double *field_vals = malloc(iters * sizeof(*field_vals));
    for (int j = 0; j < iters; j += 1) {
                //(x,y,z) position of current point
                double z = j * step;
                double dB = 0;
                //calculate the element of the total sum from each mode
                //(see G&J eq 3)
                for (int i = 0; i < 20; i++) {
                    double A = sigma * sqrt(G(all_modes[i], corr_len, 1) / G_sum);
                    //TODO currently using cos(pol) for xi hat in G&J eq 3
                    double zeta = cos(all_modes[i].pol);
                    //TODO exponential factor only the real part
                    double exp_factor = cos((all_modes[i].k * z) + all_modes[i].phase);
                    dB += A * exp_factor;
                }
                field_vals[j] = B0 + dB;
            }
        }
    }
    return field_vals;
}*/

/*
 * Generate 3d turbulence and return an array of sampled values from it
 * of size iters^3 with step size step. The magnetic turbulence variance
 * sigma, the correlation length corr_len, and the background magnetic
 * field strength B0 are all included as parameters.
 *
 * Returns: a pointer to an array of size iters^3 containing the field strength
 * at the sampled points. Each of these is the background field plus the
 * x-component of the calculated dB for that step. Returns NULL if anything goes wrong.
 * TODO return a whole 3d matrix of vectors instead of just a list at specific (x,y)
 */
vector *gen3dturb(double sigma, double corr_len, vector B0, int iters, double step, int dim) {
    //TODO allow k generation to be more flexible
    //generate 10 modes with wavenumbers 1-3 (as in Baring turb)
    mode *lower_modes = genModes(1, 3, 10, dim);
    //generate 10 modes with wavenumbers 3-10 (as in Baring turb)
    mode *higher_modes = genModes(3, 10, 10, dim);
    //Error checking
    if (lower_modes == NULL || higher_modes == NULL) {
        printf("Mode generation failed");
        return NULL;
    }

    mode *all_modes = malloc(20 * sizeof(*all_modes));
    memcpy(all_modes     , lower_modes, 10 * sizeof(*all_modes));
    memcpy(all_modes + 10, higher_modes, 10 * sizeof(*all_modes));
    free(lower_modes);
    lower_modes = NULL;
    free(higher_modes);
    higher_modes = NULL;

    //calculate sum of G for all modes first
    double G_sum = 0;
    for (int i = 0; i < 20; i++) {
        G_sum += G(all_modes[i], corr_len, 3);
    }

    //go through each of the (x,y,z) locations and put a vector there
    vector *field_vals = malloc((int) pow(iters, 3) * sizeof(*field_vals));
    //temporarily setting constant values of x and y TODO return a 3d matrix
    for (int j = 0; j < iters; j += 1) {
        for (int k = 0; k < iters; k += 1) {
            for (int l = 0; l < iters; l += 1) {
                vector loc = {j * step, k * step, l * step};
                vector dB = {0, 0, 0};
                //calculate the element of the total sum from each mode
                //(see G&J eq 3)
                for (int i = 0; i < 20; i++) {
                    mode m = all_modes[i];
                    vector vprime = modeCoords(loc, m);
                    double A = sigma * sqrt(G(m, corr_len, dim) / G_sum);
                    vector xp_hat = {cos(m.theta) * cos(m.phi), cos(m.theta) * sin(m.phi), -1 * sin(m.theta)};
                    vector yp_hat = {-1 * sin(m.phi), cos(m.phi), 0};
                    //these two factors times x prime hat and y prime hat respectively are the
                    //real components of what we will add to dB
                    double xph_factor = A * cos(m.pol) * cos(m.k * vprime.z + m.phase);
                    double yph_factor = A * -1 * sin(m.pol) * sin(m.k * vprime.z + m.phase);
                    //total addition to dB from this mode
                    vector res = addVectors(multVector(xph_factor, xp_hat), multVector(yph_factor, yp_hat));
                    dB = addVectors(dB, res);
                }
                field_vals[VEC_FLD_IDX(loc.x,loc.y,loc.z,iters)] = addVectors(B0, dB); //TODO confirm 3d array works
            }
        }
    }
    return field_vals;
}

/*
 * Converts a vector into the coordinates of a particular mode,
 * that is, given an (x,y,z) returns a new vector (x', y', z') for this
 * particular mode according to Eq 5 of G&J
 */
vector modeCoords(vector v, mode m) {
    vector retvec;
    retvec.x = (cos(m.theta) * cos(m.phi) * v.x) + (cos(m.theta) * sin(m.phi) * v.y) + (-1 * sin(m.theta) * v.z);
    retvec.y = (-1 * sin(m.phi) * v.x) + (cos(m.phi) * v.y);
    retvec.z = (sin(m.theta) * cos(m.phi) * v.x) + (sin(m.theta) * sin(m.phi) * v.y) + (cos(m.theta) * v.z);
    return retvec;
}

vector addVectors(vector v, vector w) {
    vector result = {v.x+w.x, v.y+w.y, v.z+w.z};
    return result;
}

vector multVector(double s, vector v) {
    vector result = {s*v.x, s*v.y, s*v.z};
    return result;
}



/*
 * Generates num random mode structs with wavenumbers between min_k and max_k
 * TODO currently sets theta=0 in all of them!
 * Caller is responsible for deallocating returned pointer!
 */
mode *genModes(double min_k, double max_k, int num, int dim){
    mode *modes = malloc(num * sizeof(*modes));
    for (int i = 0; i < num; i++) {
        modes[i].phase = genRandDouble(0, 2*M_PI);
        modes[i].phi = genRandDouble(0, 2*M_PI);
        modes[i].pol = genRandDouble(0, 2 * M_PI);
        if (dim == 1){
            modes[i].theta = 0;
        }
        else if (dim == 3) {
            modes[i].theta = genRandDouble(0, M_PI);
        }
        else {
            printf("Only support dim=1 and dim=3!");
            return NULL;
        }
        modes[i].k = genRandDouble(min_k, max_k);
    }
    return modes;
}

/*
 * Calculates the G of a particular mode (see G&J eq 7)
 * TODO currently implemented for 1d
 */
double G(mode m, double corr_len, int dim) {
   double delk = 1;//TODO I have set delk = 1 for now
   double gamma;
   double delV;
    if (dim == 1){
        gamma = 5.0/3.0;
        delV = delk;
    }
    else if (dim == 3) {
        gamma = 11.0/3.0;
        delV = 4 * M_PI * m.k * m.k * delk;
    }
    else {
        printf("Only support dim=1 and dim=3!");
        return -1.0;
    }
   double denom = 1 + pow(corr_len * m.k, gamma);
   return delV / denom;
}


/*
 * Generates a random double between start and stop
 */
double genRandDouble(double start, double stop){
   double random = rand() / (double)(RAND_MAX);
   double scaled = random * (stop - start);
   return scaled + start;
};
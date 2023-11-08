#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct mode {
    double polarization;
    double wavenumber;
    double phase;
    double theta;
    double phi;
} mode;

double *gen1dturb(double sigma, double corr_len, double B0, int iters, double step);
double genRandDouble(double start, double stop);
struct mode* genModes(double min_k, double max_k, int num);
double G(mode mode1, double corr_len);


//
// Created by 4jung on 11/6/2023.
//
int main(void) {
    double *turb = gen1dturb(.2, 1, 0, 100, 1);
    //TODO figure out how to graph the result
    for (int i = 0; i < 100; i++) {
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
 * at the sampled points. Currently, each of these is the previous plus the
 * x-component of the calculated dB for each step.
 * TODO also return dB y component array?
 */
double *gen1dturb(double sigma, double corr_len, double B0, int iters, double step){
    //TODO allow k generation to be more flexible
    //generate 10 modes with wavenumbers 1-3 (as in Baring turb)
    mode *lower_modes = genModes(1, 3, 10);
    //generate 10 modes with wavenumbers 3-10 (as in Baring turb)
    mode *higher_modes = genModes(3, 10, 10);
    //calculate sum of G for all modes first
    double G_sum = 0;
    for (int i = 0; i < 10; i++) {
        G_sum += G(lower_modes[i], corr_len);
        G_sum += G(higher_modes[i], corr_len);
    }
    mode *all_modes = malloc(20 * sizeof(*all_modes));
    memcpy(all_modes     , lower_modes, 10 * sizeof(*all_modes));
    memcpy(all_modes + 10, higher_modes, 10 * sizeof(*all_modes));
    free(lower_modes);
    free(higher_modes);

    //go through each of the locations in z
    double cumm = B0;
    double *field_vals = malloc(iters * sizeof(*field_vals));
    for (int j = 0; j < iters; j += 1) {
        double z = j * step;
        double dB = 0;
        //calculate the element of the total sum from each mode
        //(see G&J eq 3)
        for (int i = 0; i < 20; i++) {
            double A = sigma * sqrt(G(all_modes[i], corr_len) / G_sum);
            //TODO currently using cos(polarization) for xi hat in G&J eq 3
            double zeta = cos(all_modes[i].polarization);
            //TODO exponential factor only the real part
            double exp_factor = cos((all_modes[i].wavenumber * z) + all_modes[i].phase);
            dB += A * zeta * exp_factor;
        }
       cumm += dB;
       field_vals[j] = cumm;
    }
    return field_vals;
}

/*
 * Generates num random mode structs with wavenumbers between min_k and max_k
 * TODO currently sets theta=0 in all of them!
 * Caller is responsible for deallocating returned pointer!
 */
mode *genModes(double min_k, double max_k, int num){
    mode *modes = malloc(num * sizeof(*modes));
    for (int i = 0; i < num; i++) {
        modes[i].phase = genRandDouble(0, 2*M_PI);
        modes[i].phi = genRandDouble(0, 2*M_PI);
        modes[i].polarization = genRandDouble(0, 2*M_PI);
        modes[i].theta = 0; //TODO change this for 3d
        modes[i].wavenumber = genRandDouble(min_k, max_k);
    }
    return modes;
}

/*
 * Calculates the G of a particular mode (see G&J eq 7)
 * TODO currently implemented for 1d
 */
double G(mode mode1, double corr_len) {
   double delV = 1;//TODO delV = delk for 1d, and I have set delk = 1 for now
   double gamma = 5.0/3.0; //TODO gamma for 1d
   double denom = 1 + pow(corr_len * mode1.wavenumber, gamma);
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
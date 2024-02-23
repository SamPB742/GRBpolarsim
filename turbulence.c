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
typedef struct field_info {
    double sigma;
    double corr_len;
    vector B0;
    mode *modes;
    double G_sum;
} field_info;

double *gen1dturb(double sigma, double corr_len, double B0, int iters, double step);
double genRandDouble(double start, double stop);
struct mode* genModes(double min_k, double max_k, int num, int dim);
double G(mode mode1, double corr_len, int dim);
field_info init_3d_field(double sigma, double corr_len, vector B0);
vector eval_field(vector loc, field_info info);
vector modeCoords(vector v, mode m);
vector addVectors(vector v, vector w);
vector multVector(double s, vector v);
double vecMag(vector v);
double vecDotProd(vector v, vector w);
vector vecCrossProd(vector v, vector w);
void histogram(double *arr, int num_vals, double start, double bin_size, int num_bins);
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
    field_info fieldInfo = init_3d_field(sigma, 1, B0v);

    double helix_rad = 1 * 0.267949; // for 15 deg
    //double helix_rad = 1; // for 45 deg
    //double helix_rad = 1 * 3.73205; //for 75 deg
    double helix_len = 2 * M_PI;
    //for helix centered on z-axis, angle with z axis is
    //arctan(2pi * (rad/len))
    int num_rots = 100;
    double stored_vals[num_rots * 100];
    for (int i = 0; i < num_rots * 100; i++) {
        //helix parameterized by i value
        double j = i / 100.0; //each rotation divided into 100 chunks
        vector loc = {helix_rad * cos(2 * M_PI * j), helix_rad * sin(2 * M_PI * j), helix_len * j};
        vector tan = {-2 * M_PI * helix_rad * sin(2 * M_PI * j), 2 * M_PI * helix_rad * cos(2 * M_PI * j), helix_len};
        vector field = eval_field(loc, fieldInfo);
        //magnitude of cross divided by magnitudes of vectors is sine
       // printf("%f\n", vecMag(vecCrossProd(field, tan)) / (vecMag(field) * vecMag(tan)));
        double sin_theta = (vecMag(vecCrossProd(field, tan)) / (vecMag(field) * vecMag(tan)));
        double B_val = vecMag(field);
        stored_vals[i] = sin_theta;
        //printf("%f\n", sin_theta);

    }

    histogram(stored_vals, num_rots * 100, 0, .001 * B0, 1500);


    //TODO figure out how to graph the result
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
                    // currently using cos(pol) for xi hat in G&J eq 3
                    double zeta = cos(all_modes[i].pol);
                    // exponential factor only the real part
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
 * Initializes a 3d turbulent B-field by generating the modes and precalculating
 * G_sum. Stores relevant information for calculating field vector at an arbitrary
 * point for later use in the returned field_info struct
 *
 * Returns: a field_info struct containing the sigma, correlation length,
 * B0 vector, generated modes, and precalculated G_sum
 */
field_info init_3d_field(double sigma, double corr_len, vector B0) {
    field_info fieldInfo;
    fieldInfo.sigma = sigma;
    fieldInfo.corr_len = corr_len;
    fieldInfo.B0 = B0;
    //TODO allow k generation to be more flexible
    //generate 10 modes with wavenumbers 1-3 (as in Baring turb)
    mode *lower_modes = genModes(1, 3, 10, 3);
    //generate 10 modes with wavenumbers 3-10 (as in Baring turb)
    mode *higher_modes = genModes(3, 10, 10, 3);
    //Error checking
    if (lower_modes == NULL || higher_modes == NULL) {
        printf("Mode generation failed");
        fieldInfo.modes = NULL;
    }

    mode *all_modes = malloc(20 * sizeof(*all_modes));
    memcpy(all_modes     , lower_modes, 10 * sizeof(*all_modes));
    memcpy(all_modes + 10, higher_modes, 10 * sizeof(*all_modes));
    free(lower_modes);
    lower_modes = NULL;
    free(higher_modes);
    higher_modes = NULL;

    fieldInfo.modes = all_modes;

    //calculate sum of G for all modes first
    double G_sum = 0;
    for (int i = 0; i < 20; i++) {
        G_sum += G(all_modes[i], corr_len, 3);
    }
    fieldInfo.G_sum = G_sum;

    return fieldInfo;
}

/*
 * Given vector location, a set of generated modes, and TODO whatever else it needs
 * returns the magnetic field at the provided location.
 * Assumes 3d field
 */
vector eval_field(vector loc, field_info info) {
    vector dB = {0, 0, 0};
    //calculate the element of the total sum from each mode
    //(see G&J eq 3)
    for (int i = 0; i < 20; i++) {
        mode m = info.modes[i];
        vector vprime = modeCoords(loc, m);
        double A = info.sigma * sqrt(G(m, info.corr_len, 3) / info.G_sum);
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
    return addVectors(info.B0, dB);
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

double vecMag(vector v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
};

double vecDotProd(vector v, vector w) {
    return (v.x*w.x) + (v.y*w.y) + (v.z*w.z);
}

vector vecCrossProd(vector v, vector w) {
    vector prod = {(v.y*w.z)-(v.z*w.y), (v.z*w.x)-(v.x*w.z), (v.x*w.y)-(v.y*w.x)};
    return prod;
}
/*
 * Generates num random mode structs with wavenumbers between min_k and max_k
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
}
/*
 * goes through the values in arr and tallies frequency in each of num_bins
 * bins of size bin_size. outputs results in a format for gnuplot histogram plotting
 */
void histogram(double *arr, int num_vals, double start, double bin_size, int num_bins) {
    double max = start + (bin_size * num_bins);
    int freqs[num_bins] = {};

    //calculate frequencies
    for (int i = 0; i < num_vals; i++) {
        if(arr[i] < start) {
            printf("Warning: %f is less than start value of %f, omitting from histogram.\n", arr[i], start);

        }
        else if(arr[i] > max) {
            printf("Warning: %f is greater than max value of %f, omitting from histogram.\n", arr[i], max);
        }
        else {
            freqs[(int) floor((arr[i] - start) / bin_size)] += 1;
        }
    }
    //print bin centers and freqs on successive lines,
    //required form for gnuplot histogram
    for (int i = 0; i < num_bins; i++) {
        double bin_center = start + (bin_size * i) + (bin_size / 2);

        printf("%f %d\n", bin_center, freqs[i]);
    }
}


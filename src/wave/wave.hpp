# ifndef WAVE_FUNCTIONS_HPP
# define WAVE_FUNCTIONS_HPP

/* */
void progressMessage(int timeStep, int nt, float dt, int sId, float sx, float sy, float sz);

/* */
float * rickerGeneration(float delay, float fmax, float dt, int nt);

/* */
void pml_dampers(float * damp1D, float * damp2D, float * damp3D, float factor, int nb);

/* */
void setWaveField(float * U_pas, float * U_pre, float * U_fut, int nPoints);

/* */
void applyWavelet(float * U_pre, float * wavelet, int timeStep, int sId);

/* */
void pml_wavePropagation(float * V, float * U_pas, float * U_pre, float * U_fut, float * damp1D, float * damp2D, float * damp3D, int nb, int nxx, int nyy, int nzz, float dx, float dy, float dz, float dt);

/* */
void wavefieldUpdate(float * U_pas, float * U_pre, float * U_fut, int nPoints);

/* */
void buildSeismogram(float * U_fut, int nxx, int nyy, int nzz, float * seismogram, int timeStep, int nt, int * rx, int * ry, int * rz, int nrec);

# endif
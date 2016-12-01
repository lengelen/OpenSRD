#ifndef MAIN_H
#define MAIN_H

#include "ParticleTracking.h"
#include "ReadAndWrite.h"

//Find intersection of rays with with surfaces
void findIntersection( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 5 coefficients
void findIntersection_2ndorder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 8 coefficients
void findIntersection_3rdorder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 12 coefficients

inline int nancheck(double x); //check if nan - double version
inline int nancheck2(float x); //check if nan - float version

CameraParams Initialize(Settings s, int camera); // Initialization of camera parameters in a CameraParams-var for input of settings and camera number
void calcCorners(Size boardSize, float squareSize, vector<Point3f>& corners, int position); // Calculate theoretical location of 3D points of camera pose pattern

//Optimization of coefficients
vector<double> compute_error(float Lx, float Ly, int ErrorMetric, CameraParams Camera, vector<Point3f> Pixels, vector<Point3f> f, real_1d_array Coeff); // Compute errors Ef for all feature points of one camera
void combined_error3(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of all cameras combined for set of 3 cameras
void combined_error2(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of all cameras combined for set of 2 cameras
void combined_error1(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of single camera
real_1d_array optimizeCoef(real_1d_array Coeff, frameOptimizer Frame); // Finds optimal coefficicients for frame, requires initial guess and FrameOptimizer containing all necessary input
 

#endif

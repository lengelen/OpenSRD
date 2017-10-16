/*
	OpenSRD 1.0.0
    Copyright (C) 2017  Lukas Engelen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef MAIN_H
#define MAIN_H

#include "ParticleTracking.h"
#include "ReadAndWrite.h"
#include "minpack.hpp"
#include "objectTypes.h"

//Find intersection of rays with with surfaces
void findIntersection_constant( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 1 coefficients
void findIntersection_linear( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 3 coefficients
void findIntersection_1storder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 5 coefficients
void findIntersection_2ndorder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 8 coefficients
void findIntersection_3rdorder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface characterized by 12 coefficients
void findIntersection_Paper( int n, double x[], double fvec[], int *iflag, double params[] );

inline int nancheck(double x); //check if nan - double version
inline int nancheck2(float x); //check if nan - float version

CameraParams Initialize(Settings s, int camera); // Initialization of camera parameters in a CameraParams-var for input of settings and camera number
void createFeatureFile(); // Creates text file with theoretical location of reference pattern
/// Main function
vector<uchar> errorfunction(bool flag, float Lx, float Ly, int gridSize, int ParameterAmount, CameraParams Camera, Mat Image, real_1d_array Coeff, Point2f Middle,  uchar i_minThreshold, uchar i_maxThreshold);
Point2f computeMiddle(bool flag, float Lx, float Ly, int gridSize, int ParameterAmount, CameraParams Camera, Mat Image, real_1d_array Coeff, uchar i_minThreshold,  uchar i_maxThreshold);

void combinedErrorFunction(const real_1d_array &x, real_1d_array &fi, void *ptr);
#endif

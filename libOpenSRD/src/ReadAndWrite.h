/*
 * ReadAndWrite.h
 *
 *  Created on: Nov 27, 2016
 *      Author: lengelen
 */

#ifndef READANDWRITE_H_
#define READANDWRITE_H_

#include "objectTypes.h"

//Reading functions
bool readStringList(vector<string>& l, string dir ); // Finds all images in directory and stores their names in vector of strings

// Write functions
void writeMatToFile(cv::Mat& m, const char* file, string  outputDirectory); // Write a Float-matrix to outputDirectory, using filename file
void writeVecToFile(vector<Corner> features, string  path); //write a Float-Vector of corner points to file with path name path
void writePointsToFile(vector<Point3f> features, string  file, string outputDirectory); //write a Float-Vector of 3D feature points to outputDirectory, using filename file
void writeArray(vector<vector<real_1d_array> > plane, string  file, string outputDirectory); // Write an Array of surface coefficients to outputDirectory, using filename file
void writeArrayErrors(vector<vector <double> > errors, string  file, string outputDirectory); // Write an Array of surface coefficients to outputDirectory, using filename file
void saveCameraParams( const string& filename, string outputDirectory, Mat Rotationmatrix, Mat TranslationVector); // Saves the calibration results to outputDirectory, using filename file

#endif /* READANDWRITE_H_ */

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
void writeMatToFile(cv::Mat& m, const char* filename, string  outputDirectory); // Write a Float-matrix to outputdirectory, using filename
void writeVecToFile(vector<Corner> features, string  file); //write a Float-Vector of 3D feature points to file with path name out
void writeArray(vector<vector<real_1d_array> > plane, string  out); // Write an Array of surface coefficients to file with path name out
void writeArrayErrors(vector<vector <double> > errors, string  out); // Write an Array of surface coefficients to file with path name out
void saveCameraParams( const string& filename, Mat Rotationmatrix, Mat TranslationVector); // Saves the results into file with filename

#endif /* READANDWRITE_H_ */

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

#ifndef READANDWRITE_H_
#define READANDWRITE_H_

#include "objectTypes.h"

//Reading functions
bool readStringList(vector<string>& l, string dir ); // Finds all images in directory and stores their names in vector of strings

// Write functions
void writeVecToFile(vector<Corner> features, string  path); //write a Float-Vector of corner points to file with path name path
void writeCentersToFile(vector<vector<Point2f>> features, string  path);
void writePointsToFile(vector<Point3f> features, string  file, string outputDirectory); //write a Float-Vector of 3D feature points to outputDirectory, using filename file
void writeArray(vector<vector<real_1d_array> > plane, string  file, string outputDirectory); // Write an Array of surface coefficients to outputDirectory, using filename file
void saveCameraParams( const string& filename, string outputDirectory, Mat Rotationmatrix, Mat TranslationVector); // Saves the calibration results to outputDirectory, using filename file

#endif /* READANDWRITE_H_ */

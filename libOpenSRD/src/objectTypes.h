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

#ifndef OBJECTTYPES_H_
#define OBJECTTYPES_H_

/// Include header files and external source files

#include <omp.h>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <complex>
#include <iostream>
#include <ctype.h>
#include <sys/param.h>
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "opencv2/core/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/opencv.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <map>
#include <list>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <functional>
#include <unordered_set>

using namespace std;
using namespace cv;


#include "stdafx.h"
#include "optimization.h"

using namespace alglib;

/// Global parameters

const double PI  =3.141592653589793238463;
const float  PI_F=3.14159265358979f;
const double rw=1.336;

/// Template classes

template <class Type> class worldCoords;
template <class Type> class pixelCoordinates;
template <class Type> class Ray;
template <class Type> class Normals_constant;
template <class Type> class Normals_1storder;
template <class Type> class Normals_2ndorder;
template <class Type> class Normals_3rdorder;
template <class Type> class Normals_Paper;
template <class Type> class distance_to_Line;
template <class Type> class distance_to_Point;
template <class Type> class Angle;
template <class Type> class Angle2;
template <class Type> class waterSurface_constant;
template <class Type> class waterSurface_1storder;
template <class Type> class waterSurface_2ndorder;
template <class Type> class waterSurface_3rdorder;
template <class Type> class waterSurface_Paper;
struct arrange1;
struct arrange2;
struct arrange3;
template <class T> struct prod;
template <class T> struct prod2;
template <class T> struct Waterray;
template <class T> struct axis;
template <class T> struct Homogenous;
template <class T> struct takeNorm;
template <class T> struct Rotationmatrix;
template <class T> struct Rotationmin;
template <class T> struct Rotation;
template <class T> struct Distance;
template <class T> struct SubtractNormalized;
template <class T> struct Subtract;
template <class Type> class angleImage1;
template <class Type> class angleImage2;
template <class T> struct reduce_Z;
template <class T> struct createCorneriD;
template <class Type> class angleImageWithPoint;
template <class Type> class findMatches;

/// Object classes

class Corner{

private:
	int iD;
	Point3f coords; // Homog coordinates (u,v,1)
	Point2f vel; // Velocity of feature point in image plane between frames
	bool Found; // Boolean to indicate corner found or not

public:
	// Constructors
	Corner(){};
	Corner(Point3f iCoords, int i_iD){
		coords=iCoords;
		iD=i_iD;
		Found=true;
	};
	// Get private object variables
	int getiD(){
		return iD;
	}
	Point3f getCoords(){
		return coords;
	}
	bool getFound(){
		return Found;
	}
	// Print coordinates to string
	string printCoordinates(){
		return to_string(coords.x)+"\t"+to_string(coords.y);
	}
	// Change object variables
	void setCoords(Point3f iCoords){
		coords=iCoords;
		vel=Point2f(iCoords.x-coords.x, iCoords.y-coords.y);
	};
	void setFound(bool newValue){
		Found=newValue;
	}
	// Predict location in next image
	Point3f predictPoint(){
		return Point3f(coords.x+vel.x,coords.y+vel.y,1);
	}

};
class CameraParams
{ /// Stores all info about camera (extrinsic, intrinsic and initial image points q)
public:
	Mat R;
	Mat T;
	Mat K;
	Point3f c;
	Mat distCoeffs;

public:
    // Constructors
	CameraParams(){
	}
	CameraParams(const CameraParams& cpy){
		R=cpy.R.clone();
		T=cpy.T.clone();
		K=cpy.K.clone();
		c=Point3f(cpy.c.x,cpy.c.y,cpy.c.z);
		distCoeffs=cpy.distCoeffs.clone();
	}
	CameraParams (Mat i_R, Mat i_T, Mat i_K, Point3f i_c, Mat i_distCoeffs) {

        R=i_R;
        T=i_T;
        K=i_K;
        c=i_c;
        distCoeffs=i_distCoeffs;
    }

};
class imageOptimizer
{
public:

	// Camera set-up and surface model parameters
	size_t NumberOfCameras;
	vector<CameraParams> Cameras;
	real_1d_array Scaling;
	int Params;
	float Lx;
	float Ly;
	uchar minThreshold;
	uchar maxThreshold;
	int gridSize;

	// Optimization parameters
	double Epsg;
	double Epsf;
	double Epsx;
	ae_int_t MaxIts;
	double DiffStep;

	// Input per image
	vector<Mat> Images;
	Point2f centerPoint;
	int min_u;
	int max_u;
	int min_v;
	int max_v;

public:
	// Constructors
	imageOptimizer (){
		NumberOfCameras=0;
		vector<CameraParams> Cameras;
		Scaling="[0]";
		Params=0;
		Lx=0;
		Ly=0;
		minThreshold=0;
		maxThreshold=0;
		gridSize=0;

		Epsg=0;
		Epsf=0;
		Epsx=0;
		MaxIts=0;
		DiffStep=0;
		centerPoint=Point2f();
		min_u=0;
		max_u=0;
		min_v=0;
		max_v=0;

	}
	imageOptimizer (const imageOptimizer& cpy){
		NumberOfCameras=cpy.NumberOfCameras;
		Cameras=cpy.Cameras;
		minThreshold=cpy.minThreshold;
		maxThreshold=cpy.maxThreshold;
		gridSize=cpy.gridSize;
		Scaling=cpy.Scaling;
		Params=cpy.Params;
		Lx=cpy.Lx;
		Ly=cpy.Ly;
		Epsg=cpy.Epsg;
		Epsf=cpy.Epsf;
		Epsx=cpy.Epsx;
		MaxIts=cpy.MaxIts;
		DiffStep=cpy.DiffStep;
		Images=cpy.Images;
		centerPoint=cpy.centerPoint;
		min_u=cpy.min_u;
		max_u=cpy.max_u;
		min_v=cpy.min_v;
		max_v=cpy.max_v;
		}
	imageOptimizer (int i_NumberOfCameras, vector<CameraParams> i_Cameras, int i_Params, float i_Lx, float i_Ly, uchar i_minThreshold, uchar i_maxThreshold,
		int i_gridSize, string i_Scaling, double i_epsg, double i_epsf, double i_epsx, int i_maxits, double i_diffStep, int i_min_u, int i_max_u, int i_min_v, int i_max_v) {

		NumberOfCameras=i_NumberOfCameras;
		Cameras=i_Cameras;
		Scaling = i_Scaling.c_str();
		Params=i_Params;
		Lx=i_Lx;
		Ly=i_Ly;
		minThreshold=i_minThreshold;
		maxThreshold=i_maxThreshold;
		gridSize=i_gridSize;

		Epsg=i_epsg;
		Epsf=i_epsf;
		Epsx=i_epsx;
		MaxIts=i_maxits;
		DiffStep=i_diffStep;
		min_u=i_min_u;
		max_u=i_max_u;
		min_v=i_min_v;
		max_v=i_max_v;
	}

	void changeImages(vector<Mat> newImages)
	{
		Images=newImages;
		}
 };
class Settings
{ /// Settings class for reconstruction
public:

	// Change length scale for global optimization
	void changeLengthscales(float i_Lx, float i_Ly){
		Lx=i_Lx;
		Ly=i_Ly;
	}

    //Read serialization for this class
    void read(const FileNode& node)
    { ///Read node in settings file
		node["NumberOfCameras"]  >> NumberOfCameras;
		node["TypeCameraPose"]  >> TypeCameraPose;
        node["ThreadAmount"]  >> ThreadAmount;
        node["SaveCameraPose"]  >> SaveCameraPose;
        node["ShowCorners"]  >> ShowCorners;

        node["InputDirectory1"]  >> InputDirectory1;
        node["InputDirectory2"]  >> InputDirectory2;
        node["InputDirectory3"]  >> InputDirectory3;
        node["InputInitial1"] >> InputInitial1;
        node["InputInitial2"] >> InputInitial2;
        node["InputInitial3"] >> InputInitial3;      
        node["InputReference1"]  >> InputReference1;
        node["InputReference2"]  >> InputReference2;
        node["InputReference3"]  >> InputReference3;

        node["CalibrationFile1"]  >> CalibrationFile1;
        node["CalibrationFile2"]  >> CalibrationFile2;       
     	node["CalibrationFile3"]  >> CalibrationFile3;
     	
     	node["OutputCameraPose1"]  >> OutputCameraPose1;
        node["OutputCameraPose2"]  >> OutputCameraPose2;
        node["OutputCameraPose3"]  >> OutputCameraPose3;
        node["OutputDirectory1"]  >> OutputDirectory1;
        node["OutputDirectory2"]  >> OutputDirectory2;
        node["OutputDirectory3"]  >> OutputDirectory3;
        node["OutputFileName"]  >> OutputFileName;
        node["OutputFileCenters"]  >> OutputFileCenters;
        
        node["Lx"]  >> Lx;
        node["Ly"]  >> Ly;
        node["RefPatternSize_Width" ] >> RefPatternSize.width;
        node["RefPatternSize_Height"] >> RefPatternSize.height;
        node["gridSize" ] >> gridSize;


        node["DiffStep"]  >> DiffStep;
        node["Epsf"]  >> Epsf;
        node["Epsg"]  >> Epsg;        
        node["Epsx"]  >> Epsx;      
        node["InitialGuess"]  >> InitialGuess;
        node["MatchesThreshold"]  >> MatchesThreshold;
        node["MaxIts"]  >> MaxIts;
        node["MinDistance"]  >> MinDistance;
        node["Scaling"]  >> Scaling;
        node["SurfaceModelParameters"]  >> SurfaceModelParameters;
        node["ResponseThreshold"]  >> ResponseThreshold;
        node["ResponseRadius"]  >> ResponseRadius;
        node["minThreshold"]  >> minThreshold;
        node["maxThreshold"]  >> maxThreshold;
        node["min_u"]  >> min_u;
        node["max_u"]  >> max_u;
        node["min_v"]  >> min_v;
        node["max_v"]  >> max_v;

    }

   static bool readStringList(vector<string>& l, string dir )
    { ///Read all files in directory and store in vector string
        l.clear();
        string filepath;
        DIR *dp;
        struct dirent *dirp;
        struct stat filestat;

        dp = opendir( dir.c_str() );
        if (dp == NULL)
        {
            cout << "Error opening " << dir << endl;
            return false;
        }

        while ((dirp = readdir( dp )))
        {
            filepath = dir + "/" + dirp->d_name;

            // If the file is a directory (or is in some way invalid) we'll skip it
            if (stat( filepath.c_str(), &filestat )) continue;
            if (S_ISDIR( filestat.st_mode ))         continue;
            l.push_back((filepath));

        }

        closedir( dp );
        std::sort( l.begin(), l.end() );
        cout <<" Number of images loaded: "<<l.size()<<endl;
        return true;
    }
public:
   

    //General settings
    int NumberOfCameras;		 // Number of cameras's used
   	bool TypeCameraPose;		 // Type of input for camera pose estimation: true(1) for image, false (0) for file
   	int ThreadAmount;			 // Number of threads used to parallelize
   	
   	//Output settings
    bool SaveCameraPose;		 // Save camera pose estimation in file
  	bool ShowCorners;			 // Show Detected corners in each image
   	
   	//Input
 	string InputDirectory1;      // The name of the directory of images - camera 1
    string InputDirectory2;      // The name of the directory of images - camera 2
  	string InputDirectory3;      // The name of the directory of images - camera 3
 	string InputInitial1;		 // The name of the file of vertex locations during camera pose estimation - camera 1
 	string InputInitial2;		 // The name of the file of vertex locations during camera pose estimation - camera 2
 	string InputInitial3;		 // The name of the file of vertex locations during camera pose estimation - camera 3
 	string InputReference1;		 // The name of the reference image (no water) - camera 1
 	string InputReference2;		 // The name of the reference image (no water) - camera 2
 	string InputReference3;		 // The name of the reference image (no water) - camera 3

	//Calibration
    string CalibrationFile1;     // The name of the calibration file used - camera 1
    string CalibrationFile2;     // The name of the calibration file used - camera 2
    string CalibrationFile3;     // The name of the calibration file used - camera 3

	//Output
 	string OutputCameraPose1;	 // The name of the output-file for camera pose estimation - camera 1
 	string OutputCameraPose2;	 // The name of the output-file for camera pose estimation - camera 2
 	string OutputCameraPose3;	 // The name of the output-file for camera pose estimation - camera 3
 	string OutputDirectory1;	 // The name of the output-directory for feature coordinates (pixels) - camera 1
 	string OutputDirectory2;	 // The name of the output-directory for feature coordinates (pixels) - camera 2
 	string OutputDirectory3;	 // The name of the output-directory for feature coordinates (pixels) - camera 3
    string OutputFileName;		 // Name of output-file coefficients
    string OutputFileCenters; 	 // Name of output-file of centers of reconstructed area

    //Detection parameters
    float Lx;					 // Length scale in (lateral) x-direction
    float Ly;					 // Length scale in (streamwise) y-direction 
    Size RefPatternSize;         // The size of the reference board -> Number of items by width and height
    int gridSize;				 // Size of checquerboard squares

    //Optimization parameters
    double DiffStep;			 // Numerical differentiation step (calculation gradient)
    double Epsf;				 // Min function change as stopping condition
    double Epsg;				 // Min Gradient norm as stopping condition
    double Epsx;				 // Min step size as stopping condition
    string InitialGuess;	 	 // Initial guess for coefficients according to model (length=SurfaceModelParameters)
    float MatchesThreshold;		 // Threshold of maximum distance (pixels) between predicted and detected feature point
    int MaxIts;					 // Max amount of iterations for optimization procedure
    float MinDistance;			 // Minimum distance between corner points
    string Scaling;	 			 // Set scale of coefficients (length=SurfaceModelParameters)
    int SurfaceModelParameters;	 // Number of coefficients/parameters in surface model
    int ResponseRadius;			 // Radius for calculation error measure
    float ResponseThreshold;	 // Threshold for corner response to keep interesting points
    uchar minThreshold;			 // Minimum Image intensity Imin
    uchar maxThreshold;			 // Maximum Image intensity Imax
	int min_u;					 // Minimum u-coordinate of reconstructed image area
	int max_u;					 // Maximum u-coordinate of reconstructed image area
	int min_v;					 // Minimum v-coordinate of reconstructed image area
    int max_v;					 // Maximum v-coordinate of reconstructed image area

};

#endif /* OBJECTTYPES_H_ */

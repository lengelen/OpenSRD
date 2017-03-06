/*
 * objectTypes.h
 *
 *  Created on: Nov 26, 2016
 *      Author: lengelen
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
template <class Type> class distance_to_Line;
template <class Type> class distance_to_Point;
template <class Type> class Angle;
template <class Type> class Angle2;
template <class Type> class error_I;
template <class Type> class waterSurface_constant;
template <class Type> class waterSurface_1storder;
template <class Type> class waterSurface_2ndorder;
template <class Type> class waterSurface_3rdorder;
template <class Type> class imagePoints;
template <class Type> class imagePoints2;
struct arrange1;
struct arrange2;
struct arrange3;
template <class T> struct prod;
template <class T> struct prod2;
template <class T> struct shiftPixels;
template <class T> struct error_col;
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
	vector<Point3f> f;
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
		f=cpy.f;
		distCoeffs=cpy.distCoeffs.clone();
	}
	CameraParams (Mat i_R, Mat i_T, Mat i_K, Point3f i_c, Mat i_distCoeffs, vector<Point3f> i_f) {

        R=i_R;
        T=i_T;
        K=i_K;
        c=i_c;
        f=i_f;
        distCoeffs=i_distCoeffs;
    }

};

class frameOptimizer
{ ///Class to store optimization settings, including cameras (with parameters) and vector-list of point coordinates q'
public:

	// Camera set-up and surface model parameters
	size_t NumberOfCameras;
	vector<CameraParams> Cameras;
	int ErrorMetric;
	int Params;
	float Lx;
	float Ly;
	real_1d_array Scaling;

	// Optimization parameters
	double Epsg;
	double Epsf;
	double Epsx;
	ae_int_t MaxIts;
	double DiffStep;

	// Input per image
	vector<vector<Point3f> > Pixels;
	vector<vector<Point3f> > fs;

public:
	// Constructors
	frameOptimizer (){
		NumberOfCameras=0;
		vector<CameraParams> Cameras;
		ErrorMetric=0;
		Params=0;
		Lx=0;
		Ly=0;
		Scaling="[0]";
		Epsg=0;
		Epsf=0;
		Epsx=0;
		MaxIts=0;
		DiffStep=0;
	}
	frameOptimizer (const frameOptimizer& cpy){
		NumberOfCameras=cpy.NumberOfCameras;
		Cameras=cpy.Cameras;
		ErrorMetric=cpy.ErrorMetric;
		Params=cpy.Params;
		Lx=cpy.Lx;
		Ly=cpy.Ly;
		Scaling=cpy.Scaling;
		Epsg=cpy.Epsg;
		Epsf=cpy.Epsf;
		Epsx=cpy.Epsx;
		MaxIts=cpy.MaxIts;
		DiffStep=cpy.DiffStep;
		Pixels=cpy.Pixels;
		fs=cpy.fs;
		}
	frameOptimizer (int i_NumberOfCameras, vector<CameraParams> i_Cameras, int i_ErrorMetric, int i_Params, float i_Lx, float i_Ly, string i_Scaling, double i_epsg, double i_epsf, double i_epsx, int i_maxits, double i_diffStep) {

		NumberOfCameras=i_NumberOfCameras;
		Cameras=i_Cameras;
		ErrorMetric=i_ErrorMetric;
		Params=i_Params;
		Lx=i_Lx;
		Ly=i_Ly;
		Scaling = i_Scaling.c_str();

		Epsg=i_epsg;
		Epsf=i_epsf;
		Epsx=i_epsx;
		MaxIts=i_maxits;
		DiffStep=i_diffStep;
    }
	// Change set of image points for which optimization is run
	void changePixels(vector<vector<Corner> > i_Corners)
	{
		// Reset vector of features f and pixel coordinates
		fs=vector<vector<Point3f> >(NumberOfCameras);
		Pixels=vector<vector<Point3f> >(NumberOfCameras);

		// Fill in updated values
		for(size_t i =0; i<NumberOfCameras; i++){
			for(size_t j=0; j<i_Corners[i].size(); j++){
				if(i_Corners[i][j].getFound()){
					Pixels[i].push_back(i_Corners[i][j].getCoords());
					fs[i].push_back(Cameras[i].f[j]);
				}
			}
		}
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
        node["TypeFeatureInput"]  >> TypeFeatureInput;
        node["ThreadAmount"]  >> ThreadAmount;
        
        node["SaveCameraPose"]  >> SaveCameraPose;
        node["SaveFeatureCoordinates"]  >> SaveFeatureCoordinates;
        node["SaveResiduals"]  >> SaveResiduals;
        node["ShowCorners"]  >> ShowCorners;

        node["InputDirectory1"]  >> InputDirectory1;
        node["InputDirectory2"]  >> InputDirectory2;
        node["InputDirectory3"]  >> InputDirectory3;
        node["InputFeatures1"] >> InputFeatures1;
        node["InputFeatures2"] >> InputFeatures2;
        node["InputFeatures3"] >> InputFeatures3;
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
        node["OutputFileNameResiudals"]  >> OutputFileNameResiudals;
        
        node["FeaturePatternSize_Width" ] >> FeaturePatternSize.width;
        node["FeaturePatternSize_Height"] >> FeaturePatternSize.height;
        node["Lx"]  >> Lx;
        node["Ly"]  >> Ly;
        node["RefPatternSize_Width" ] >> RefPatternSize.width;
        node["RefPatternSize_Height"] >> RefPatternSize.height;

        node["DiffStep"]  >> DiffStep;
        node["Epsf"]  >> Epsf;
        node["Epsg"]  >> Epsg;        
        node["Epsx"]  >> Epsx;      
        node["ErrorMetric"]  >> ErrorMetric;
        node["InitialGuess"]  >> InitialGuess;
        node["MatchesThreshold"]  >> MatchesThreshold;
        node["MaxIts"]  >> MaxIts;
        node["MinDistance"]  >> MinDistance;
        node["Scaling"]  >> Scaling;
        node["SurfaceModelParameters"]  >> SurfaceModelParameters;
        node["ResponseThreshold"]  >> ResponseThreshold;
        node["ResponseRadius"]  >> ResponseRadius;
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
   	bool TypeFeatureInput;		 // Type of input for surface reconstruction: true(1) for images, false (0) for file
   	bool TypeCameraPose;		 // Type of input for camera pose estimation: true(1) for image, false (0) for file
   	int ThreadAmount;			 // Number of threads used to parallelize
   	
   	//Output settings
    bool SaveCameraPose;		 // Save camera pose estimation in file
    bool SaveFeatureCoordinates; // Save feature coordinates (pixels) to text file
 	bool SaveResiduals;			 // Write mean resiudal error over all image points out to file
  	bool ShowCorners;			 // Show Detected corners in each image
   	
   	//Input
 	string InputDirectory1;      // The name of the directory of images - camera 1
    string InputDirectory2;      // The name of the directory of images - camera 2
  	string InputDirectory3;      // The name of the directory of images - camera 3
 	string InputFeatures1;		 // The name of the file of fixed locations of features viewed by camera 1
 	string InputFeatures2;		 // The name of the file of fixed locations of features viewed by camera 2
 	string InputFeatures3;		 // The name of the file of fixed locations of features viewed by camera 3 
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
    string OutputFileNameResiudals; // Name of output-file of residual errors

    //Detection parameters
	Size FeaturePatternSize;     // The size of the feature pattern -> Number of items by width and height
    float Lx;					 // Length scale in (lateral) x-direction
    float Ly;					 // Length scale in (streamwise) y-direction 
    Size RefPatternSize;         // The size of the reference board -> Number of items by width and height

    //Optimization parameters
    double DiffStep;			 // Numerical differentiation step (calculation gradient)
    double Epsf;				 // Min function change as stopping condition
    double Epsg;				 // Min Gradient norm as stopping condition
    double Epsx;				 // Min step size as stopping condition
    int ErrorMetric;		 	 // Type of error metric used
    string InitialGuess;	 	 // Initial guess for coefficients according to model (length=SurfaceModelParameters)
    float MatchesThreshold;		 // Threshold of maximum distance (pixels) between predicted and detected feature point
    int MaxIts;					 // Max amount of iterations for optimization procedure
    float MinDistance;			 // Minimum distance between corner points
    string Scaling;	 			 // Set scale of coefficients (length=SurfaceModelParameters)
    int SurfaceModelParameters;	 // Number of coefficients/parameters in surface model
    int ResponseRadius;			 // Radius for calculation error measure
    float ResponseThreshold;	 // Threshold for corner response to keep interesting points

};

#endif /* OBJECTTYPES_H_ */

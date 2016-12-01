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

#include "minpack.hpp"
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
template <class Type> class Normals;
template <class Type> class Normals_2ndorder;
template <class Type> class distance_to_Line;
template <class Type> class distance_to_Point;
template <class Type> class Angle;
template <class Type> class Angle2;
template <class Type> class error_I;
template <class Type> class waterSurface;
template <class Type> class waterSurface_2ndorder;
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
	size_t CameraAmount;
	vector<CameraParams> Cameras;
	int ErrorMetric;
	int Params;
	float Lx;
	float Ly;

	// Optimization parameters
	double epsg = 0.0000000001;
	double epsf = 0;
	double epsx = 0;
	ae_int_t maxits = 0;
	double diffStep= 0.0001;

	// Input per image
	vector<vector<Point3f> > Pixels;
	vector<vector<Point3f> > fs;

public:
	// Constructors
	frameOptimizer (){
		CameraAmount=0;
		vector<CameraParams> Cameras;
		ErrorMetric=0;
		Params=0;
		Lx=0;
		Ly=0;
	}
	frameOptimizer (int i_CameraAmount, vector<CameraParams> i_Cameras, int i_ErrorMetric, int i_Params, float i_Lx, float i_Ly) {

		CameraAmount=i_CameraAmount;
		Cameras=i_Cameras;
		ErrorMetric=i_ErrorMetric;
		Params=i_Params;
		Lx=i_Lx;
		Ly=i_Ly;
    }

	// Change Stoppingconditions optimization
	void setStoppingConditions(double i_epsg , double i_epsf, double i_epsx, ae_int_t i_maxits)
	{
		epsg = i_epsg;
		epsf = i_epsf;
		epsx = i_epsx;
		maxits = i_maxits;
	}
	// Change differentation-step
	void setDiffStep(double i_diffStep )
	{
		diffStep=i_diffStep;
	}
	// Change set of image points for which optimization is run
	void changePixels(vector<vector<Corner> > i_Corners)
	{
		// Reset vector of features f and pixel coordinates
		fs=vector<vector<Point3f> >(CameraAmount);
		Pixels=vector<vector<Point3f> >(CameraAmount);

		// Fill in updated values
		for(size_t i =0; i<CameraAmount; i++){
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
	//Write serialization for this class
    void write(FileStorage& fs) const
    {
        fs << "{"
        << "MatchesThreshold" << MatchesThreshold
        << "BoardSize_Width"  << boardSize.width
        << "BoardSize_Height" << boardSize.height
        << "PatternSize_Height" << PatternSize.height
        << "PatternSize_Width"  << PatternSize.width
        << "Square_Size"      << squareSize
        << "InputType"      << InputType
        << "InputTypeRef"      << InputTypeRef
        << "InputDirectory1"  << InputDirectory1
        << "InputReference1"  <<InputReference1
        << "InputDirectory2"  << InputDirectory2
        << "InputReference2"  <<InputReference2
        << "InputDirectory3"  << InputDirectory3
        << "InputReference3"  <<InputReference3
        << "OutputCameraPose1"  <<OutputCameraPose1
        << "OutputCameraPose2"  <<OutputCameraPose2
        << "OutputCameraPose3"  <<OutputCameraPose3
        << "OutputDirectory1"  <<OutputDirectory1
        << "OutputDirectory2"  <<OutputDirectory2
        << "OutputDirectory3"  <<OutputDirectory3
        << "InputInitial"  <<InputInitial
        << "ErrorMetric" << ErrorMetric
        << "CameraAmount" << CameraAmount
        << "ResponseThreshold" << ResponseThreshold
        << "MinDistance" << MinDistance
        << "setInitialGuess" << setInitialGuess
        << "setOptim" << setOptim
        << "SurfaceModel" << SurfaceModel
        << "Lx" << Lx
        << "Ly" << Ly
        << "OutputFileName" << OutputFileName
        << "showCorners" << showCorners
        << "saveErrors" << saveErrors
        << "saveCameraPose" << saveCameraPose
        << "OutputFileNameErrors" << OutputFileNameErrors
        << "responseRadius" << responseRadius
        << "saveFeatureCoordinates" << saveFeatureCoordinates
        << "ThreadAmount" << ThreadAmount
        << "}";
    }
    //Read serialization for this class
    void read(const FileNode& node)
    { ///Read node in settings file
        node["InputType"]  >> InputType;
        node["InputTypeRef"]  >> InputTypeRef;
        node["InputDirectory1"]  >> InputDirectory1;
        node["CalibrationFile1"]  >> CalibrationFile1;
        node["InputReference1"]  >> InputReference1;
        node["CalibrationFile2"]  >> CalibrationFile2;
        node["InputDirectory2"]  >> InputDirectory2;
        node["InputReference2"]  >> InputReference2;
        node["InputDirectory3"]  >> InputDirectory3;
     	node["CalibrationFile3"]  >> CalibrationFile3;
        node["InputReference3"]  >> InputReference3;
        node["BoardSize_Width" ] >> boardSize.width;
        node["BoardSize_Height"] >> boardSize.height;
        node["PatternSize_Width" ] >> PatternSize.width;
        node["PatternSize_Height"] >> PatternSize.height;
        node["Square_Size"]  >> squareSize;
        node["Position1"] >> Position1;
        node["Position2"]  >> Position2;
        node["Position3"]  >> Position3;
        node["InputInitial"] >> InputInitial;
        node["ErrorMetric"]  >> ErrorMetric;
        node["CameraAmount"]  >> CameraAmount;
        node["ResponseThreshold"]  >> ResponseThreshold;
        node["MinDistance"]  >> MinDistance;
        node["setInitialGuess"]  >> setInitialGuess;
        node["setOptim"]  >> setOptim;
        node["SurfaceModel"]  >> SurfaceModel;
        node["Lx"]  >> Lx;
        node["Ly"]  >> Ly;
        node["OutputFileName"]  >> OutputFileName;
        node["showCorners"]  >> showCorners;
        node["saveErrors"]  >> saveErrors;
        node["saveCameraPose"]  >> saveCameraPose;
        node["OutputFileNameErrors"]  >> OutputFileNameErrors;
        node["responseRadius"]  >> responseRadius;
        node["OutputCameraPose1"]  >> OutputCameraPose1;
        node["OutputCameraPose2"]  >> OutputCameraPose2;
        node["OutputCameraPose3"]  >> OutputCameraPose3;
        node["OutputDirectory1"]  >> OutputDirectory1;
        node["OutputDirectory2"]  >> OutputDirectory2;
        node["OutputDirectory3"]  >> OutputDirectory3;
        node["MatchesThreshold"]  >> MatchesThreshold;
        node["saveFeatureCoordinates"]  >> saveFeatureCoordinates;
        node["ThreadAmount"]  >> ThreadAmount;

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
   	bool InputType;				 // Type of input for reconstruction: true(1) for images, false (0) for file
   	bool InputTypeRef;			 // Type of input for camera pose estimation: true(1) for image, false (0) for file
    Size boardSize;              // The size of the calibration board -> Number of items by width and height
    Size PatternSize;            // The size of the pattern -> Number of items by width and height
    float squareSize;            // The size of a square in your defined unit (point, millimeter,etc).
    string InputDirectory1;      // The name of the directory of images - camera 1
    string CalibrationFile1;     // The name of the calibration file used - camera 1
    string CalibrationFile2;     // The name of the calibration file used - camera 2
    string CalibrationFile3;     // The name of the calibration file used - camera 3
    string InputDirectory2;      // The name of the directory of images - camera 2
 	string InputDirectory3;      // The name of the directory of images - camera 3
 	string InputReference1;		 // The name of the reference image (no water) - camera 1
 	string InputReference2;		 // The name of the reference image (no water) - camera 2
 	string InputReference3;		 // The name of the reference image (no water) - camera 3
 	string InputInitial;		 // The name of the file of location points f
 	string OutputCameraPose1;	 // The name of output-file for camera pose estimation - camera 1
 	string OutputCameraPose2;	 // The name of output-file for camera pose estimation - camera 2
 	string OutputCameraPose3;	 // The name of output-file for camera pose estimation - camera 3
 	string OutputDirectory1;	 // The name of output-directory for feature coordinates (pixels) - camera 1
 	string OutputDirectory2;	 // The name of output-directory for feature coordinates (pixels) - camera 2
 	string OutputDirectory3;	 // The name of output-directory for feature coordinates (pixels) - camera 3
 	vector<string> imageList1;	 // List of images - camera 1
    vector<string> imageList2;	 // List of images - camera 2
    vector<string> imageList3;	 // List of images - camera 3
    string OutputFileName;		 // Name of output-file coefficients
    int Position1;				 // Position camera 1
    int Position2;				 // Position camera 2
    int Position3;				 // Position camera 3
    int ErrorMetric;		 	 // Type of error metric used
    int CameraAmount;			 // Amount of cameras's used
    float ResponseThreshold;	 // Threshold for corner response to keep interesting points
    float MinDistance;			 // Minimum distance between corner points
    bool setInitialGuess;		 // Allow change to initial guess or not
    bool setOptim;				 // Allow change to optimization parameters or not
    int SurfaceModel;			 // Number of coefficients/parameters in surface model
    float Lx;					 // Length scale in x-direction (crosswise)
    float Ly;					 // Length scale in y-direction (streamwise)
    bool showCorners;			 // Show Detected corners in each image
    bool saveErrors;			 // Write mean error over all image points out to file
    bool saveFeatureCoordinates; // Save feature coordinates (pixels) to text file
    bool saveCameraPose;		 // Save camera pose estimation in file
    string OutputFileNameErrors; // Name of output-file errors
    int responseRadius;			 // Radius for calculation error measure
    float MatchesThreshold;		 // Threshold of maximum distance (pixels) between predicted and detected feature point
    int ThreadAmount;			 // Amount of threads used to parallelize
};

#endif /* OBJECTTYPES_H_ */

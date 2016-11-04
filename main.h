#ifndef MAIN_H
#define MAIN_H

// { /// Inlcude header files and external source files

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
double rw=1.336;


/// Self-defined classes and structures (some are adapted versions of free source code online)
/// Some classes are both defined in double and float variant

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
{ ///Class to store optimization settings, including cameras (with parameters) and vectorlist of point coordinates q'
public:

	// Initialisation of optimisation parameters
	int CameraAmount;
	vector<CameraParams> Cameras;
	int ErrorMetric;
	int Params;
	double epsg = 0.000000001;
	double epsf = 0;
	double epsx = 0;
	double diffStep=0.0001;
	ae_int_t maxits = 0;
	float Lx;
	float Ly;
	vector<vector<Point3f>> Pixels;

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
	// Change Stoppingconditions optimisation
	void setStoppingConditions(double i_epsg , double i_epsf, double i_epsx, ae_int_t i_maxits)
	{
		epsg = i_epsg;
		epsf = i_epsf;
		epsx = i_epsx;
		maxits = i_maxits;
	}
	// Change differention-step
	void setDiffStep(double i_diffStep )
	{
		diffStep=i_diffStep;
	}
	// Change set of image points for which optimisation is runned
	void changePixels(vector<vector<Point3f>> i_Pixels)
	{
		Pixels=i_Pixels;
	}
	// Swap two cameras and remove first
	void swapCameras(int cam1, int cam2)
	{
		Pixels[cam1-1]=Pixels[cam2-1];
		Pixels.erase (Pixels.begin()+cam2-1);
		Cameras[cam1-1]=Cameras[cam2-1];
		Cameras.erase (Cameras.begin()+cam2-1);
		CameraAmount--;
	}
	// Remove 1 camera from camera list
	void removeCamera(int cam)
	{
		Pixels.erase (Pixels.begin()+cam);
		Cameras.erase (Cameras.begin()+cam);
		CameraAmount--;
	}

};
class Settings
{ /// Settings class for reconstruction
public:

	enum {Colinear=1, Disparity =2};
	// Change Lengthscales
	void changeLengthscales(float i_Lx, float i_Ly){
		Lx=i_Lx;
		Ly=i_Ly;
	}
	//Write serialization for this class
    void write(FileStorage& fs) const
    {
        fs << "{"
        << "BoardSize_Width"  << boardSize.width
        << "BoardSize_Height" << boardSize.height
        << "Square_Size"      << squareSize
        << "InputType"      << InputType
        << "InputTypeRef"      << InputTypeRef
        << "InputFilename1"   << InputFilename1
        << "InputDirectory1"  << InputDirectory1
        << "InputReference1"  <<InputReference1
        << "InputFilename2"   << InputFilename2
        << "InputDirectory2"  << InputDirectory2
        << "InputReference2"  <<InputReference2
        << "InputFilename3"   << InputFilename3
        << "InputDirectory3"  << InputDirectory3
        << "InputReference3"  <<InputReference3
        << "OutputCameraPose1"  <<OutputCameraPose1
        << "OutputCameraPose2"  <<OutputCameraPose2
        << "OutputCameraPose3"  <<OutputCameraPose3
        << "InputInitial"  <<InputInitial
        << "ErrorMetric" << ErrorMetric
        << "Opticalflowguess" << Opticalflowguess
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
        << "}";
    }
    //Read serialization for this class
    void read(const FileNode& node)
    { ///Read node in settings file
        node["InputType"]  >> InputType;
        node["InputTypeRef"]  >> InputTypeRef;
        node["InputDirectory1"]  >> InputDirectory1;
        node["InputFilename1"]  >> InputFilename1;
        node["CalibrationFile1"]  >> CalibrationFile1;
        node["InputReference1"]  >> InputReference1;
        node["CalibrationFile2"]  >> CalibrationFile2;
        node["InputFilename2"]  >> InputFilename2;
        node["InputDirectory2"]  >> InputDirectory2;
        node["InputReference2"]  >> InputReference2;
        node["InputDirectory3"]  >> InputDirectory3;
      	node["InputFilename3"]  >> InputFilename3;
     	node["CalibrationFile3"]  >> CalibrationFile3;
        node["InputReference3"]  >> InputReference3;
        node["BoardSize_Width" ] >> boardSize.width;
        node["BoardSize_Height"] >> boardSize.height;
        node["Square_Size"]  >> squareSize;
        node["Position1"] >> Position1;
        node["Position2"]  >> Position2;
        node["Position3"]  >> Position3;
        node["InputInitial"] >> InputInitial;
        node["ErrorMetric"]  >> ErrorMetric;
        node["Opticalflowguess"]  >> Opticalflowguess;
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
        for(size_t i=0;i<l.size();i++)
        	cout <<l[i]<<endl;
        return true;
    }
public:
   	bool InputType;				 // Type of input for reconstruction: true(1) for images, false (0) for file
   	bool InputTypeRef;			 // Type of input for camera pose etsimation: true(1) for image, false (0) for file
    Size boardSize;              // The size of the board -> Number of items by width and height
    float squareSize;            // The size of a square in your defined unit (point, millimeter,etc).
    string InputDirectory1;      // The name of the directory of images - camera 1
    string InputFilename1;       // The name of the input file with image points corners and frame counter - camera 1
    string CalibrationFile1;     // The name of the calibration file used - camera 1
    string CalibrationFile2;     // The name of the calibration file used - camera 2
    string CalibrationFile3;     // The name of the calibration file used - camera 3
    string InputFilename2;		 // The name of the input file with image points corners and frame counter - camera 2
    string InputDirectory2;      // The name of the directory of images - camera 2
    string InputFilename3;       // The name of the input file with image points corners and frame counter - camera 3
 	string InputDirectory3;      // The name of the directory of images - camera 3
 	string InputReference1;		 // The name of the reference image (no water) - camera 1
 	string InputReference2;		 // The name of the reference image (no water) - camera 2
 	string InputReference3;		 // The name of the reference image (no water) - camera 3
 	string InputInitial;		 // The name of the file of location points f
 	string OutputCameraPose1;	 // The name of output-file for camera pose estimation - camera 1
 	string OutputCameraPose2;	 // The name of output-file for camera pose estimation - camera 2
 	string OutputCameraPose3;	 // The name of output-file for camera pose estimation - camera 3
 	vector<string> imageList1;	 // List of images - camera 1
    vector<string> imageList2;	 // List of images - camera 2
    vector<string> imageList3;	 // List of images - camera 3
    string OutputFileName;		 // Name of output-file coefficients
    int Position1;				 // Position camera 1
    int Position2;				 // Position camera 2
    int Position3;				 // Position camera 3
    int ErrorMetric;		 	 // Type of error metric used
    bool Opticalflowguess;		 // Usage of estimate of point with Optical flow if not found
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
    bool saveCameraPose;		 // Save camera pose estimation in file
    string OutputFileNameErrors; // Name of output-file errors
    int responseRadius;			 // Radius for calculation error measure
};
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
template <class T> struct reduce_dist;

/// Self-defined functions
/// Some functions are both defined in double and float variant

//Random functions used
bool point_comparator( cv::Point3f a, cv::Point3f b); // Comparison function of Point3f
int nancheck(double x); // Check if nan - double version
int nancheck2(float x); // Check if nan - float version
bool compare(float a, float b); // Compair-function of two float variables to check for equality (lower than threshold)
bool operator==(const Point2f& pt1, const Point2f& pt2); // Definition of operator to compare two points (2D), necessary for removedupes
void removedupes(std::vector<Point2f> & vec); // Removes doubles in list of 2D points

//Find intersection of rays with with surfaces
void findIntersection( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface charct by 5 coeffs
void findIntersection_2ndorder( int n, double x[], double fvec[], int *iflag, double params[] ); //Find intersection of ray cq' with surface charct by 8 coeffs

//Reading functions
static bool readStringList(vector<string>& l, string dir ); // Finds all images in directory and stores their names in vector of strings
static inline void read(const FileNode& node, Settings& x, const Settings& default_value); // Reads node in settings file and stores it in variable

// Write functions
void writeMatToFile(cv::Mat& m, const char* filename, string  outputDirectory); // Write a Float-matrix to outputdirectory, using filename
void writeVecToFile(vector<Point3f> features, string  file, size_t timestep); //write a Float-Vector of 3D feature points to file with path name out
void writeArray(vector <real_1d_array> plane, string  out); // Write an Array of surface coefficients to file with path name out
void writeArrayErrors(vector <double> errors, string  out); // Write an Array of surface coefficients to file with path name out
static void saveCameraParams( const string& filename, Mat Rotationmatrix, Mat TranslationVector); // Saves the results into file with filename

//Functions to detect corners
float computeResponse(Mat points); // Compute reponse value for center pixel based on points located on circle
Mat getPoints5(Mat image, size_t x, size_t y); // Returns Mat 4x4 with inetnsity value at 16 points in circle with radius 5
Mat getPoints10(Mat image, size_t x, size_t y); // Returns Mat 4x4 with inetnsity value at 16 points in circle with radius 10
Mat corner_detect5(const size_t h, const size_t w,  Mat image); // Compute response function for image with radius 5
Mat corner_detect10(const size_t h, const size_t w,  Mat image); // Compute response function for image with radius 10
vector<Point2f> removeDoubles(vector<Point2f> points, Mat response, float mindist);
// Remove points located closes from each other than mindist
// Results in single individual point for each group of points around corner position
// Individual corner points are weighted average with weight equal to corner response


//Processing images & detecting feature points
vector<Point2f> sortGrid(vector<Point2f> corners, Size boardSize);
// Sort grid of corners points autamatically in ordered pattern
// Columns in streamwise direction
// Rows in crosswise direction

vector<Point2f> create_points(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, bool distortion, Mat cameraMatrix, Mat distCoeffs, bool showCorners); // Do all operations on image to obtain a sorted list of corner points
vector<Point3f> readFeaturesImage(vector<string> imageList, size_t frame, size_t nimages, CameraParams cam, Settings s, string fileName); // Finds all feature points in image, returns vector of sorted imagepoints
CameraParams Initialize(Settings s, int camera); // Initialization of camera parameters in a CameraParams-var for input of settings and camera number
void calcCorners(Size boardSize, float squareSize, vector<Point3f>& corners, int position); // Calculate theoretical location of 3D points of camera pose pattern

//Optimization of coefficients
vector<double> compute_error(float Lx, float Ly, int ErrorMetric, CameraParams Camera, vector<Point3f> Pixels, real_1d_array Coeff); // Compute errors Ef for all feature points of one camera
void combined_error3(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of all cameras combined for set of 3 cameras
void combined_error2(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of all cameras combined for set of 2 cameras
void combined_error1(const real_1d_array &x, real_1d_array &fi, void *ptr); // Compute errors of single camera
real_1d_array optimizeCoef(real_1d_array Coeff, frameOptimizer Frame); // Finds optimal coefficicients for frame, requires initial guess and FrameOptimizer containing all necessary input
 

#endif

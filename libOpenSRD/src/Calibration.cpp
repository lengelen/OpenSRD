#include <iostream>
#include <ctype.h>
#include <sys/param.h>
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <fstream>
#include "opencv2/core/core.hpp"
#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/videoio/videoio.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include <opencv2/imgcodecs.hpp>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <map>
#include <list>
#include <sstream>



#ifndef _CRT_SECURE_NO_WARNINGS
# define _CRT_SECURE_NO_WARNINGS
#endif

using namespace cv;
using namespace std;
RNG rng(12345);
clock_t prevTimestamp = 0;



/*static void help()
{
    cout <<  "This is a camera calibration sample." << endl
    <<  "Usage: calibration configurationFile"  << endl
    <<  "Near the sample file you'll find the configuration file, which has detailed help of "
    "how to edit it.  It may be any OpenCV supported file format XML/YAML." << endl;
}
*/

class Settings
{/// Class to load the predefined settings from settingsfile
public:

     
    void write(FileStorage& fs) const                        //Write serialization for this class
    {
        fs << "{"
        << "BoardSize_Width"  << boardSize.width
        << "BoardSize_Height" << boardSize.height
        << "Square_Size"	<< squareSize
        << "adaptImage"	<< adaptImage
        << "writePerview"   << writePerview
        << "Write_outputFileName"  << outputFileName
        << "InputDirectory"  << InputDirectory
        << "Show_UndistortedImage"	<< showUndistorted
		<< "Showcalib"	<< Showcalib
        << "}";
    }
    void read(const FileNode& node)                          //Read serialization for this class
    {
        node["BoardSize_Width" ] >> boardSize.width;
        node["BoardSize_Height"] >> boardSize.height;
        node["Square_Size"]  >> squareSize;
        node["nframes"] >> nframes;
        node["Write_outputFileName"] >> outputFileName;
        node["Show_UndistortedImage"] >> showUndistorted;
        node["InputDirectory"]  >> InputDirectory;
        node["Showcalib"] >> Showcalib;
        node["writePerview"] >> writePerview;
        node["adaptImage"] >> adaptImage;
    }
      
    static bool readStringList(vector<string>& l )
    {
        l.clear();
        string dir, filepath;
        DIR *dp;
        struct dirent *dirp;
        struct stat filestat;
        
        cout << "dir to get files of: " << flush;
        getline( cin, dir );  // gets everything the user ENTERs
        
        dp = opendir( dir.c_str() );
        if (dp == NULL)
        {
            cout << "Error opening " << dir << endl;
            return false;
        }
        
        while ((dirp = readdir( dp )))
        {
            filepath = dir + "/" + dirp->d_name;
            cout<<filepath;
            
            // If the file is a directory (or is in some way invalid) we'll skip it
            if (stat( filepath.c_str(), &filestat )) continue;
            if (S_ISDIR( filestat.st_mode ))         continue;
            l.push_back((filepath));
            }
        
        closedir( dp );
        return true;
    }
public:
    Size boardSize;              // The size of the board -> Number of items by width and height
    float squareSize;            // The size of a square in your defined unit (point, millimeter,etc).
    bool writePerview;        	 // Write extrinsic parameters
    string outputFileName;       // The name of the file where to write
    string InputDirectory;       // The name of the directory of calibration images 
    bool showUndistorted;        // Show undistorted images after calibration
    int nframes;				 // Number of frames used for calibration
    bool Showcalib;				 // Show calib results or not
    bool adaptImage;			 // Adapt calibration images before calibration procedure
    
    vector<string> imageList;
    size_t atImageList;
};

 static inline void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
{ ///Function to read nodes of settings file
	if(node.empty())
		x = default_value;
	else
		x.read(node);
 }

static inline void write(FileStorage& fs, const String&, const Settings& s )
{ ///Function to write nodes in output file
	s.write(fs);
}
static double computeReprojectionErrors(const vector<vector<Point3f> >& objectPoints,
                                        const vector<vector<Point2f> >& imagePoints,
                                        const vector<Mat>& rvecs, const vector<Mat>& tvecs,
                                        const Mat& cameraMatrix, const Mat& distCoeffs,
                                        vector<float>& perViewErrors )
{ ///Compute reprojection error, based on input objectpoints (theoret), image points and estimate of camera parameters
    vector<Point2f> imagePoints2;
    size_t totalPoints = 0;
	double totalErr = 0, err;
    perViewErrors.resize(objectPoints.size());
    
     for(size_t i = 0; i < objectPoints.size(); ++i )
    {
        projectPoints(Mat(objectPoints[i]), rvecs[i], tvecs[i],
                      cameraMatrix, distCoeffs, imagePoints2); //Compute estimation of image points based on the estimated camera parameters

        err = norm(Mat(imagePoints[i]), Mat(imagePoints2), NORM_L2); //Compute difference in image plane between original and backtraced image points
        size_t n = objectPoints[i].size();
        //Compute errors related to single and total reprojection
        perViewErrors[i] = (float) std::sqrt(err*err/n);
		totalErr        += err*err;
        totalPoints     += n;
    }
    
    return std::sqrt(totalErr/totalPoints);
}

static void calcChessboardCorners(Size boardSize, float squareSize, vector<Point3f>& corners)
{ ///Calculates the theoretical locations of calibration grid
    corners.resize(0);
	for( int j = 0; j < boardSize.height; j++ ) //iterate over boardheight
		for( int i = 0; i < boardSize.width; i++ ){ //iterate ober boardwidth
            corners.push_back(Point3f(float(j*squareSize),
                                      float(i*squareSize), 0));
        
        
    }
}

static bool runCalibration( vector<vector<Point2f> > imagePoints,
                           Size imageSize, Size boardSize,
                           float squareSize,
                           Mat& cameraMatrix, Mat& distCoeffs,
                           vector<Mat>& rvecs, vector<Mat>& tvecs,
                           vector<float>& reprojErrs,
                           double& totalAvgErr )
{ ///Runs the calibration given the input paramters, stores result (camera paremeters AND reprojection errors)
    cameraMatrix = Mat::eye(3, 3, CV_64F);  //initialize 
    distCoeffs = Mat::zeros(5, 1, CV_64F);
    
    vector<vector<Point3f> > objectPoints(1);
    calcChessboardCorners(boardSize, squareSize, objectPoints[0]); //compute theroretical calibration grid
    for( size_t i = 1; i < imagePoints.size(); i++ )
        objectPoints.push_back(objectPoints[0]); //store objects points
    
    TermCriteria criteria = TermCriteria(TermCriteria::COUNT +TermCriteria::EPS, 500, 0.0001 ); //Termination criteria for pose estimation
    calibrateCamera(objectPoints, imagePoints, imageSize, cameraMatrix,
                    distCoeffs, rvecs, tvecs,  0, criteria); //peform the actual callibration
                     
    bool ok = checkRange(cameraMatrix) && checkRange(distCoeffs); //check if values or not NAN or inf
    
    totalAvgErr = computeReprojectionErrors(objectPoints, imagePoints,
                                            rvecs, tvecs, cameraMatrix, distCoeffs, reprojErrs); //compute the reprojection errors
    
    return ok;
}


static void saveCameraParams( const string& filename,
                      Size imageSize, Size boardSize,
                      float squareSize,
                      const Mat& cameraMatrix, const Mat& distCoeffs,
                      const vector<float>& reprojErrs,
                      double totalAvgErr, int nframes  )
{ ///Saves the results into file with filename
    FileStorage fs( filename, FileStorage::WRITE ); //Open filestorage file
    //Time clock
    time_t t;
    time( &t );
    struct tm *t2 = localtime( &t );
    char buf[1024];
    strftime( buf, sizeof(buf)-1, "%c", t2 );
    
    //Writing out the results
    fs << "calibration_time" << buf;
    fs << "image_width" << imageSize.width;
    fs << "image_height" << imageSize.height;
    fs << "board_width" << boardSize.width;
    fs << "board_height" << boardSize.height;
    fs << "squareSize" << squareSize;
    fs << "nframes - max" << nframes;
    fs << "nframes - camera" <<  (int)reprojErrs.size();
      
    fs << "camera_matrix" << cameraMatrix;
    fs << "distortion_coefficients" << distCoeffs;
  
    
    fs << "avg_reprojection_error - camera 1" << totalAvgErr;
    if( !reprojErrs.empty())
    {
    fs << "per_view_reprojection_errors -camera 1" << Mat(reprojErrs);
    }
    fs.release();        
	}

static bool readStringList(vector<string>& l, string dir )
{ //Read vector-string of filenames in directory
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
     
    size_t length=l.size();
	char str[256] = "";

	snprintf(str, sizeof(str), "%zu", length);

	fprintf(stdout, "Images loaded: '");
	fflush(stdout);
	write(fileno(stdout), str, strlen(str));
	fprintf(stdout, "'.\n");
    
    return true;
}

static vector<Point2f> findPoints(Mat new_image, Size board_sz, bool Showcalib)
{ /// Locate the internal corners of the calibration (chessboard) pattern in the image (+possible to show found corners)
    Mat viewGray;
    vector<Point2f> pointbuf;
    if( new_image.rows==0){
        return pointbuf;
    }
    int chessBoardFlags = CALIB_CB_ADAPTIVE_THRESH | CALIB_CB_NORMALIZE_IMAGE; //Flags to change if we want to change calibration settings
    //Find chessboard corners
    bool found = findChessboardCorners( new_image, board_sz, pointbuf, chessBoardFlags );
    if (!found)
		{   cout<<"Corners not found in image "<<endl;
		}
		    
    // improve the found corners' coordinate accuracy
    if(found)
    {   cout<<".";
        cvtColor(new_image, viewGray, cv::COLOR_BGR2GRAY); //change to grayscale because cornerSubPix only works on grayscale images
        
        cornerSubPix( viewGray, pointbuf, Size(11,11),
                     Size(-1,-1), TermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 200, 0.001 )); //refine found corner positions

     if(Showcalib){ //boolean to indicate if detected corners should be shown
            
            drawChessboardCorners( new_image, board_sz, Mat(pointbuf), found ); // Draw corners detected
            
            //Write out amount of corners found (in memory efficient way)
            size_t length=pointbuf.size();
			char str[256] = "";
			snprintf(str, sizeof(str), "%zu", length);
			fprintf(stdout, "Corners found: '");
			fflush(stdout);
			write(fileno(stdout), str, strlen(str));
			fprintf(stdout, "'.\n");

			//Plot the detected corners in the image
            int r = 2;
            for( size_t i = 0; i < pointbuf.size(); i++ )
            { circle(new_image, pointbuf[i], r, Scalar(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255)), -1, 8, 0 ); //Draw circles on the detected corner positions
                }
        
            circle(new_image, pointbuf[0], r*2, Scalar(rng.uniform(0,255), rng.uniform(0,255), rng.uniform(0,255)), -1, 8, 0); //Draw larger circle at origin
            
            //Show image with corners
            namedWindow("Corners found", WINDOW_NORMAL);
			resizeWindow("Corners found", 1500, 800);
            imshow("Corners found", new_image);
            waitKey(0);//found ? 500 : 1000);
            
        }
    destroyWindow("Corners found");
    }
    return pointbuf;//return detected corners
}

static bool runAndSave(const string& outputFilename,
                       const vector<vector<Point2f> >& imagePoints,
                       Size imageSize, Size boardSize, float squareSize,
                       Mat& cameraMatrix,
                       Mat& distCoeffs,
                       bool writePerview, int nframes)
{ ///Gloabl function to run and save the claibration
    vector<Mat> rvecs, tvecs; //Initialize
    vector<float> reprojErrs;
    double totalAvgErr = 0;
    
    bool calib_ok = runCalibration(imagePoints, imageSize, boardSize, squareSize, //Run the calibration, true if no problems occured
                             cameraMatrix, distCoeffs,
                             rvecs, tvecs, reprojErrs, totalAvgErr);
    printf("%s. avg reprojection error = %.2f\n",
           calib_ok ? "Calibration succeeded" : "Calibration failed",
           totalAvgErr); //Print succes or not 
    
   
    if( calib_ok) //If OK save calibration output
        saveCameraParams( outputFilename, imageSize,
                         boardSize, squareSize,
                         cameraMatrix, distCoeffs,
                         writePerview ? reprojErrs : vector<float>(),
                         totalAvgErr, nframes );
    
    return (calib_ok);
    
}
Mat processImage (Mat view, double alpha, int beta)
{ ///Adapt contrast and pixel intensity in image: g(x,y)=alpha*I(x,y)+beta
	
	 Mat new_image = Mat::zeros( view.size(), view.type() ); //initailize result image
			if(beta!=0 && alpha !=1.0){ //If both do not change image intensity nothing to be done
            // Do the operation new_image(i,j) = alpha*image(i,j) + beta
            for( int y = 0; y < view.rows; y++ )
                { for( int x = 0; x < view.cols; x++ )
                    { for( int c = 0; c < 3; c++ )
                        {   new_image.at<Vec3b>(y,x)[c] =
                            saturate_cast<uchar>( alpha*( view.at<Vec3b>(y,x)[c] ) + beta );
                        }
                    }
                }
			}
			else 
			{new_image=view;
				}
			return new_image;
}

int main()
{  ///Main function to run calibration
	
	//Read the settings of the calibration
	Settings s;
    string inputSettings;
    cout << "Give InputsettingsFile: " << flush; //Read which settingsfile has to be used
    getline( cin, inputSettings );  // gets everything the user ENTERs
    string inputSettingsFile="PRE_PROCESS/"+inputSettings;
    FileStorage fs(inputSettingsFile, FileStorage::READ); // Read the settings
    if (!fs.isOpened())
        {
        cout << "Could not open the configuration file: \"" << inputSettingsFile << "\"" << endl;
        return -1;
        }
    fs["Settings"] >> s;
    fs.release();   // close Settings file
   
	// Assign all settings to variables
	Size board_sz = s.boardSize;
    Size imageSize;
    float squareSize = s.squareSize;
    Mat cameraMatrix, distCoeffs;
    int nframes=s.nframes;  
    string InputDirectory=s.InputDirectory;
    bool showUndistorted=s.showUndistorted;
	bool Showcalib=s.Showcalib;
    bool writePerview = s.writePerview;
    bool adaptImage=s.adaptImage;
    string outputFilename = s.outputFileName;
    
    // Declaration of variables
    int nimages;    
    vector<vector<Point2f> > imagePoints;
    vector<string> imageList;
	int beta;
	double alpha;
	
    bool check=readStringList(imageList,InputDirectory); //Read all image names in inputdirectory
    nimages = (int)imageList.size(); //amount of calibration images
    if( !check) return 0; 
    
    for(int i = 0;i<nimages && i<nframes;i=i+1) //Iterate over all images
        {   
			Mat view =imread(imageList[i], 1);  //Read in image   
			
			if(i==0) imageSize =view.size(); //initialize image size
			cout <<"Imagesize "<<imageSize<<endl;
            Mat new_image;
            if(adaptImage ==true) { //If image inensity has to changed
			
			if(i==0)
			{ //Check in first image how the lighting is and choose based on that image appropriate alpha and beta values
				string line;
				namedWindow( "Image", WINDOW_AUTOSIZE );
			 	imshow("Image", view);
				waitKey(0);
				destroyWindow("Image");
				cout << "Give alpha for new image (double so use point e.g. 1.5): " << flush;
				cin >> alpha;
				cout << "Give beta for new image (integer): " << flush;
				cin >> beta;
				}
			new_image= processImage(view, alpha, beta); }//adapt image intensity based on chosen alpha and beta
			else new_image=view; //otherwise keep same image
			                                     
            vector<Point2f> points= findPoints(new_image, board_sz, Showcalib); //Locate corners in processed image
            cout <<"."<<endl;
            if( points.size()== board_sz.width*board_sz.height) //Check if detected amount of corners matches with pattern size
                    {imagePoints.push_back(points);}
       
        }
	
	bool done= runAndSave(outputFilename, imagePoints, imageSize,
                              board_sz, squareSize,
                              cameraMatrix, distCoeffs,
                              writePerview,  nframes); //Run entire calibration
	if(!done) //then error occured during calib
        cout<<"Process failed";
	if(showUndistorted) //show undistorted images
	{  
		Mat view, rview, map1, map2;  
		//Determines rthe undistortion and rectifcation map, for more details-->see openCV website 
		initUndistortRectifyMap(
		cameraMatrix, distCoeffs, Mat(),
		getOptimalNewCameraMatrix(cameraMatrix, distCoeffs, imageSize, 1, imageSize, 0), imageSize,CV_16SC2, map1, map2);

		for(size_t i = 0; i < imageList.size(); i++ ) //iterate over all images
			{
				view = imread(imageList[i], 1); //read in image
				if(view.empty()) //then no image
					continue;
				remap(view, rview, map1, map2, INTER_LINEAR); //unidstort and rectify image
				namedWindow("Image View undistorted", WINDOW_AUTOSIZE ); 
				imshow("Image View undistorted", rview);
				waitKey(0); //wait for keybord press of user
		  }
	destroyWindow("Image View undistorted");
	
	}
	
    return 0;
  }


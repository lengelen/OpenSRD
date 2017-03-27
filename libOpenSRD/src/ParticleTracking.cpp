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
#include "ParticleTracking.h"
#include "ReadAndWrite.h"
/// Definition template classes

template <class Type> class distance_to_Line
{ // Computes distance between 2D Float point and line (charct by start- and end-point),
  // returns (point.x, point.y, dist)

    private:
	double a; double b; double c;

    public:
            // Constructor
        	distance_to_Line (Point2f line_start, Point2f line_end) {

        		// (x- p1X) / (p2X - p1X) = (y - p1Y) / (p2Y - p1Y)
        		 a = (double) line_start.y - line_end.y;
        		 b = (double) line_end.x - line_start.x;
        		 c = (double) line_start.x * line_end.y - line_end.x * line_start.y;
            }

            // The function call
            Type operator ( ) ( Point2f point )
            {
            	double distance =abs(a * point.x + b * point.y + c) / sqrt(a * a + b * b);
            	return Point3f(point.x,point.y,distance);
            }
        };
template <class Type> class distance_to_Point
{ // Computes distance between 2D Float point and one predefined 2D Float point that is initialized in constructor,
  // returns (point.x, point.y, dist)

    private:
	Point2f point;

    public:
            // Constructor
		distance_to_Point (Point2f i_point) {
        		point= i_point;
            }

            // The function call
            Type operator ( ) ( Point2f point2 )
            {
            	double distance = (double)((point2.x - point.x) * (point2.x - point.x)+ (point2.y - point.y) * (point2.y - point.y));
            	return Point3f(point2.x,point2.y,distance);

            }
        };
template <class Type> class angleImage1
{ // Compute  angle between two 2D image points in way 1: -pi ...pi
private:
    Point2f p1;   // inverse matrices
public:
    // Constructor
    angleImage1 (Point2f i_p1) {
        p1=i_p1;
    }

    // The function call
    Type operator ( ) ( Point2f& p2 ) const
    {
    	if(!compare(p1.x,p2.x) ||!compare(p1.y,p2.y) ){
    		return atan2(p2.y-p1.y,p2.x-p1.x);

    	}
    	else return 100; // 2 points are actually the same, return impossible to result to indicate
    	}

};
template <class Type> class angleImage2
{ // Compute  angle between two 2D image points in way 2: 0--2pi
private:
    Point2f p1;   // inverse matrices
public:
    // Constructor
    angleImage2 (Point2f i_p1) {
        p1=i_p1;
    }

    // The function call
    Type operator ( ) ( Point2f& p2 ) const
    {
    	if(!compare(p1.x,p2.x) ||!compare(p1.y,p2.y) ){
    		float ang=atan2(p2.y-p1.y,p2.x-p1.x);
    		if(ang<0)
    			return ang+2*PI_F;
    		else return ang;
    	}
    	else return 100; // 2 points are actually the same, return impossible to result to indicate
    	}

};
template <class T> struct reduce_Z
{ // Throws away z-coordinate of 3D float point and only retains pixel coordinates
    T operator() (Point3f& p) const {
        return Point2f(p.x,p.y);

    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct createCorneriD
{ // Creates corner object from iD and image coordinates
   T operator() (Point2f& p, int& iD) const {

		return Corner(Point3f(p.x,p.y,1),iD);

		}
		typedef T first_argument_type;
		typedef T second_argument_type;
		typedef T result_type;
	};
template <class Type> class angleImageWithPoint
{ // Compute  angle between two 2D image points in way 1: -pi..pi with point coordinate as well
  // result is 3D point with 3rd dim angle
private:
    Point2f p1;   // inverse matrices
public:
    // Constructor
    angleImageWithPoint (Point2f i_p1) {
        p1=i_p1;
    }

    // The function call
    Type operator ( ) ( Point2f& p2 ) const
    {
    	return Point3f(p2.x,p2.y,atan2(p2.y-p1.y,p2.x-p1.x));
    	     }
};
template <class Type> class findMatches
{ // Class object <vector> of point coordinates (locations detected from  image) and threshold value: max distance (in pixels) between dete
  // in case of none or multiple candidates that lie within search radius of predicted location-> set Found on false
  // in case only one falls within search circle-> update coordinates and set Found on true
        private:
            vector<Point2f>candidates;
            float Threshold;
        public:
            // Constructor
            findMatches (vector<Point2f>& i_candidates, float& iThreshold) {
                candidates=i_candidates;
                Threshold=iThreshold;
            }
            // The function call
            Type operator ( ) ( Corner& q) const
            {
            	Point3f temp;
            	size_t count=0;
            	Point3f predCorner=q.predictPoint(); // Predict point location based on stored location from previous time step and velocity
                  for(size_t i=0;i<candidates.size();i++){ // Iterate over candidates
                      	  if((predCorner.x-candidates[i].x)*(predCorner.x-candidates[i].x)+(predCorner.y-candidates[i].y)*(predCorner.y-candidates[i].y)<Threshold){
                		temp=Point3f(candidates[i].x, candidates[i].y, 1);
                		count++;
                	  }}

                  if(count!=1){ // None or more than 1 possible candidates within search circle around predicted point
                	  q.setFound(false);
                	  return q;}
                  else{ // One detected point falls within search circle
                	  q.setCoords(temp);
                	  q.setFound(true);
                	  return q;
                  }

               }
};
struct arrange1
{ ///sort image points according to increasing u
 bool operator() ( Point2f a, Point2f b ){
               return a.x <= b.x;
    }
} sortingu;
struct arrange2
{ ///sort image points according to increasing distance param (3rd elem)
 bool operator() ( Point3f a, Point3f b ){
               return a.z <= b.z;
    }
} sortingZ;
struct arrange3
{ ///sort image points according to increasing distance from origin
 bool operator() ( Point2f a, Point2f b ){
               return pow(a.x,2)+pow(a.y,2) <= pow(b.x,2)+pow(b.y,2);
    }
} sortingOrigin;

/// Simple functions to compare variables etc

inline bool compare(float a, float b)
{ /// Compare-function of two float variables to check for equality (lower than threshold 1/100000)
    if(fabs(a - b) < (1.0 / 100000))
        return true;
    else
        return false;
}
inline bool operator==(const Point2f& pt1, const Point2f& pt2)
{ /// Definition of operator to compare two points (2D)
    return ((pt1.x == pt2.x) && (pt1.y == pt2.y));
}
namespace std
{ // Necessary for removedupes
    template<>
    struct hash<Point2f>
    {
        size_t operator()(Point2f const& pt) const
        {
            return (size_t)(pt.x*100 + pt.y);
        }
    };
}
inline void removedupes(std::vector<Point2f> & vec)
{ /// Removes doubles in list of 2D points
    std::unordered_set<Point2f> pointset;

    auto itor = vec.begin();
    while (itor != vec.end())
    {
        if (pointset.find(*itor) != pointset.end())
        {
            itor = vec.erase(itor);
        }
        else
        {
            pointset.insert(*itor);
            itor++;
        }
    }
}
inline bool point_comparator( cv::Point3f a, cv::Point3f b) {
	//Comparison function of 2 3D float points, true if second point b has a smaller 3rd dim than first point a
    return (b.z < a.z);
}

/// Functions to detect or process feature points

vector<Point2f> sortPattern(vector<Point2f> corners, Size boardSize)
{ /// Sort grid of corners points automatically in ordered pattern
  /// Columns in streamwise direction
  /// Rows in crosswise direction
	vector<Point3f> delta1(corners.size());
	vector<Point3f> delta2(corners.size());
	vector<Point2f> extremes(4);
	vector<Point3f> anglesOrigin(3);
	vector<Point2f> PointsNoOrigin(3);
	vector<Point3f> distances(corners.size());
	vector<Point2f> sortedcorners;


	float stepX1, stepX2, stepY1, stepY2;
	#pragma omp simd
	for(size_t t=0;t<corners.size(); t++){ // Compute angle-measure for all detected image points

		vector<float> anglesImage1(corners.size());
		vector<float> anglesImage2(corners.size());

		transform(corners.begin(), corners.end(), anglesImage1.begin(), angleImage1<float>(corners[t])); // Compute angles between -pi..pi
		transform(corners.begin(), corners.end(), anglesImage2.begin(), angleImage2<float>(corners[t])); // Compute angles between 0..2pi

		anglesImage1.erase(std::remove_if(anglesImage1.begin(), anglesImage1.end(),[](const float& x) {return x==100;}), anglesImage1.end()); // Remove angle that corresponds with same point: indicated with artificial angle 100rad
		anglesImage2.erase(std::remove_if(anglesImage2.begin(), anglesImage2.end(),[](const float& x) {return x==100;}), anglesImage2.end()); // Remove angle that corresponds with same point: indicated with artificial angle 100rad

		// Find min and max in list of angles between considered point and rest of points
		// Store result together with 2D image coordinates (u,v)
		auto temp1= minmax_element(anglesImage1.begin(),anglesImage1.end());
		delta1[t]=Point3f(corners[t].x, corners[t].y, *temp1.second-*temp1.first);
		auto temp2= minmax_element(anglesImage2.begin(),anglesImage2.end());
		delta2[t]=Point3f(corners[t].x, corners[t].y, *temp2.second-*temp2.first);

	}

	// Sort them according to angle-measure: z-coordinate, with smallest first
	sort(delta1.begin(),delta1.end(), sortingZ);
	sort(delta2.begin(),delta2.end(), sortingZ);
	// Only retain best 4
	delta1.resize(4);
	delta2.resize(4);

	// Ideally angle measure is Pi/2, so remove ones that deviate more than 3/4 Pi
	delta1.erase(std::remove_if(delta1.begin(), delta1.end(),[](const Point3f& x) {return x.z > PI_F*3/4;}), delta1.end());
	delta2.erase(std::remove_if(delta2.begin(), delta2.end(),[](const Point3f& x) {return x.z > PI_F*3/4;}), delta2.end());

	if((int)delta1.size()==4) // First angle measure is sufficient to detect all 4 extremes
		transform(delta1.begin(), delta1.end(), extremes.begin(), reduce_Z<Point2f>());
	else if((int)delta2.size()==4)// Second angle measure is sufficient to detect all 4 extremes
			transform(delta2.begin(), delta2.end(), extremes.begin(), reduce_Z<Point2f>());
	else if((int)(delta1.size()+delta2.size())<4){ // Combination of both does not give 4 good ones: something wrong
		cout<<"No 4 outers found"<<endl;
		cout<<"Program will be ended"<<endl;
		exit(1);
	}
	else{

		//Remove z-coordinate (angle measure) from list of interesting points
		vector<Point2f> extremes1(delta1.size());
		vector<Point2f> extremes2(delta2.size());

		transform(delta1.begin(), delta1.end(), extremes1.begin(), reduce_Z<Point2f>());
		transform(delta2.begin(), delta2.end(), extremes2.begin(), reduce_Z<Point2f>());


		extremes1.insert( extremes1.end(), extremes2.begin(), extremes2.end() );
		removedupes(extremes1); //remove all points that already were in extremes 1
		extremes=extremes1;
	}
	if(extremes.size()!=4)
	{
		cout<<"No 4 outers found"<<endl;
		cout<<"Program will be ended"<<endl;
		exit(1);
	}
	//Sort extremes according to distance to upper-left corner and take closest one as origin, sort the rest accordingly
	sort(extremes.begin(), extremes.end(),sortingOrigin);
	transform(extremes.begin()+1, extremes.end(), anglesOrigin.begin(), angleImageWithPoint<Point3f>(extremes[0]));
	sort(anglesOrigin.begin(),anglesOrigin.end(), sortingZ);
	transform(anglesOrigin.begin(), anglesOrigin.end(), PointsNoOrigin.begin(), reduce_Z<Point2f>());

	vector<Point2f> extremesSorted={extremes[0], PointsNoOrigin[2], PointsNoOrigin[0],PointsNoOrigin[1]};

	//Compute steps between two extreme points of lines // y-axis
	stepX1=(extremesSorted[1].x-extremesSorted[0].x)/(boardSize.width-1);
	stepY1=(extremesSorted[1].y-extremesSorted[0].y)/(boardSize.width-1);
	stepX2=(extremesSorted[3].x-extremesSorted[2].x)/(boardSize.width-1);
	stepY2=(extremesSorted[3].y-extremesSorted[2].y)/(boardSize.width-1);


	for(int t=0; t<boardSize.width; t++)
	{
		Point2f point1=Point2f(extremesSorted[0].x+t*stepX1,extremesSorted[0].y+t*stepY1); //edge point left
		Point2f point2=Point2f(extremesSorted[2].x+t*stepX2,extremesSorted[2].y+t*stepY2); //edge point right

		transform(corners.begin(), corners.end(), distances.begin(), distance_to_Line<Point3f> (point1, point2)); //compute distance to line for all points
		sort(distances.begin(), distances.end(), sortingZ); //sort points according to distance
		int k=0;
		vector <Point2f> temp;
		for(size_t j=0; j<corners.size(); j++){
			if(distances[j].z<50 && k<boardSize.height) //threshold of 100 pixels from line and also maximum number of points on one row
				{
				temp.push_back(Point2f(distances[j].x,distances[j].y));
				k++;
				}
			else if(k==boardSize.height){
				sort(temp.begin(), temp.end(), sortingu); //sort all points on one row by increasing u
				sortedcorners.insert(sortedcorners.end(), temp.begin(), temp.end());
				k=0;
				break;
				}
			}
		}
		return sortedcorners;
}
vector<Point2f> undistortCorners(vector<Point2f> distortedCorners, Mat cameraMatrix, Mat distCoeffs){
	//Undistort detected image points of sorted corners
		Mat r;
		undistortPoints(Mat(distortedCorners), r, cameraMatrix, distCoeffs); // undistort sorted image points
		vector<Point2f> points_undistorted(distortedCorners.size());

		float fx = cameraMatrix.at<double>(0, 0);
		float fy = cameraMatrix.at<double>(1, 1);
		float cx = cameraMatrix.at<double>(0, 2);
		float cy = cameraMatrix.at<double>(1, 2);

		#pragma omp simd
		for(int k=0; k<r.rows; k++){
			  // transformation to get correct pixel coordinates (notice ambiguous definition of origin by OpenCV)
			  points_undistorted[k]=Point2f(r.at<Vec2f>(k,0)[0] * fx + cx,r.at<Vec2f>(k,0)[1] * fy + cy);
					  }
	return points_undistorted;
}
vector<Point2f> detectPotentialCorners(Mat img, float ResponseThreshold, float minDistance, int detectionRadius){
	/// Function detects potential corners in image based on threshold values
	/// returns vector of 2D points (u,v)
	Mat response;

	// Calculate corner strength for every pixel, depending on radius of corner measure
	if(detectionRadius==5)
		response = corner_detect5(img.rows, img.cols, img);
	else if(detectionRadius==10)
		response = corner_detect10(img.rows, img.cols, img);
	else{
		cout<<"What is radius of corner response function?"<<endl;
		exit(1);
	}

	// Determine pixels with corner-response larger than ResponeThreshold
	vector<Point2f> corners;
	//#pragma omp parallel for simd collapse(2)
	for(int i=detectionRadius;i<response.cols-detectionRadius;i++)
		for(int j=detectionRadius; j<response.rows-detectionRadius; j++){
			if(response.at<float>(j,i)>ResponseThreshold)
		    		corners.push_back(Point2f(i,j));
					}

	//Remove pixels that correspond to same corner
	return removeDoubles(corners, response, minDistance);

}
vector<Point2f> removeDoubles(vector<Point2f> points, Mat response, float mindist)
{ // Remove points located closes from each other than mindist
  // Results in single individual point for each group of points around corner position
  // Individual corner points are weighted average with weight equal to corner response

	if(points.size()==0){
		cout<<"No corners found";
		exit(1);
	}

	for(size_t i=0; i<points.size()-1;i++){ // For all points until the final-1 (because all previous already checked)
		float N=0;

		Point2f avg=points[i]*response.at<float>(points[i]);

		for(size_t t=i+1; t<points.size();t++){ // For all points after considered point
			//Compute distance
			float dist=sqrt((points[i].x-points[t].x)*(points[i].x-points[t].x)+(points[i].y-points[t].y)*(points[i].y-points[t].y));
			//Remove points located closer than mindist and attach to weighted average
			if(dist<mindist){
				avg+=points[t]*response.at<float>(points[t]);
				N=N+response.at<float>(points[t]);
				points.erase(points.begin() + t);
				t--;
				}
		}
		if(N!=0) points[i]=avg/(response.at<float>(points[i])+N); //Take weighted average

	}
	return points;
}
vector<Point2f> create_points(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, Mat cameraMatrix, Mat distCoeffs, bool ShowCorners)
{ /// Do all operations on image to obtain a sorted list of corner points

	int CornerPoints=boardSize.width*boardSize.height; //number of points to be detected
	Mat response;

	// Calculate corner strength for every pixel, depending on radius of corner measure
	if(detectionRadius==5)
		response = corner_detect5(img.rows, img.cols, img);
	else if(detectionRadius==10)
		response = corner_detect10(img.rows, img.cols, img);
	else{
		cout<<"What is radius of corner response function?"<<endl;
		exit(1);
	}

	// Determine pixels with corner-response larger than ResponeThreshold
	vector<Point2f> corners;
	for(int i=detectionRadius;i<response.cols-detectionRadius;i++)
		for(int j=detectionRadius; j<response.rows-detectionRadius; j++){
		float temp=response.at<float>(j,i);
			if(temp>ResponseThreshold)
		    		corners.push_back(Point2f(i,j));
					}

	//Remove pixels that correspond to same corner
	vector<Point2f> corners2=removeDoubles(corners, response, minDistance);

	if((int)corners2.size()!=CornerPoints){ // Number of individual corner points does not match with grid size
			namedWindow( "Display window", WINDOW_NORMAL );// Create a window for display.
			resizeWindow("Display window", 800, 500);
			cout<<"error corners-> not all found: "<<corners2.size()<<endl;

			for(size_t t=0; t<corners2.size(); t++){
				cout<<t<<" corner "<<corners2[t]<<endl;
				circle( img, corners2[t], 5, Scalar(150, 150, 150), -1, 8, 0 ); // show point located closest
				}
			char k;
			imshow("Display window", img);
			k=waitKey(0);
			destroyWindow("Display window");
	return vector<Point2f>(0); // Return empty vector to alert problem
	}
	vector<Point2f> sortedcorners=sortPattern(corners2, boardSize);
	if(ShowCorners){ // Show the detected corners if users requires it
		namedWindow( "Display window", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display window", 800, 500);
		for(size_t t=0; t<sortedcorners.size(); t++){
					circle( img, sortedcorners[t], 5, Scalar(150, 150, 150), -1, 8, 0 ); // show point located closest
					circle( img, sortedcorners[t], 2
							, Scalar(0, 0, 0), -1, 8, 0 ); // show point located closest

					}
			char k;
			imshow("Display window", img);
			k=waitKey(0);
		destroyWindow("Display window");
	}

	return sortedcorners;
}

/// Functions to compute and process corner strength

inline float computeResponse(Mat points)
{ // Compute response value for center pixel based on points located on circle

	float SR=0;
	float DR=0;

	for(size_t i=0; i<4;i++){
		SR+=abs((points.at<uchar>(0,i)+points.at<uchar>(2,i))-(points.at<uchar>(1,i)+points.at<uchar>(3,i)));
		DR+=abs(points.at<uchar>(0,i)-points.at<uchar>(2,i))+abs(points.at<uchar>(1,i)-points.at<uchar>(3,i));
	}

	return (SR-DR);
}
inline Mat getPoints5(Mat image, size_t x, size_t y)
{ // Returns Mat 4x4 with intensity value at 16 points in circle with radius 5

	 Mat result(4, 4, DataType<int>::type);

	 result.at<uchar>(0,0)=image.at<uchar>(y,x-5);
	 result.at<uchar>(0,1)=image.at<uchar>(y-2,x-5);
	 result.at<uchar>(0,2)=image.at<uchar>(y-4,x-4);
	 result.at<uchar>(0,3)=image.at<uchar>(y-5,x-2);


	 result.at<uchar>(1,0)=image.at<uchar>(y-5,x);
	 result.at<uchar>(1,1)=image.at<uchar>(y-5,x+2);
	 result.at<uchar>(1,2)=image.at<uchar>(y-4,x+4);
	 result.at<uchar>(1,3)=image.at<uchar>(y-2,x+5);

	 result.at<uchar>(2,0)=image.at<uchar>(y,x+5);
	 result.at<uchar>(2,1)=image.at<uchar>(y+2,x+5);
	 result.at<uchar>(2,2)=image.at<uchar>(y+4,x+4);
	 result.at<uchar>(2,3)=image.at<uchar>(y+5,x+2);


	 result.at<uchar>(3,0)=image.at<uchar>(y+5,x);
	 result.at<uchar>(3,1)=image.at<uchar>(y+5,x-2);
	 result.at<uchar>(3,2)=image.at<uchar>(y+4,x-4);
	 result.at<uchar>(3,3)=image.at<uchar>(y+2,x-5);

	 return result;
}
inline Mat getPoints10(Mat image, size_t x, size_t y)
{// Returns Mat 4x4 with intensity value at 16 points in circle with radius 10

	 Mat result(4, 4, DataType<int>::type);

	 result.at<uchar>(0,0)=image.at<uchar>(y,x-10);
	 result.at<uchar>(0,1)=image.at<uchar>(y-4,x-10);
	 result.at<uchar>(0,2)=image.at<uchar>(y-7,x-7);
	 result.at<uchar>(0,3)=image.at<uchar>(y-10,x-4);


	 result.at<uchar>(1,0)=image.at<uchar>(y-10,x);
	 result.at<uchar>(1,1)=image.at<uchar>(y-10,x+4);
	 result.at<uchar>(1,2)=image.at<uchar>(y-7,x+7);
	 result.at<uchar>(1,3)=image.at<uchar>(y-4,x+10);

	 result.at<uchar>(2,0)=image.at<uchar>(y,x+10);
	 result.at<uchar>(2,1)=image.at<uchar>(y+4,x+10);
	 result.at<uchar>(2,2)=image.at<uchar>(y+7,x+7);
	 result.at<uchar>(2,3)=image.at<uchar>(y+10,x+4);


	 result.at<uchar>(3,0)=image.at<uchar>(y+10,x);
	 result.at<uchar>(3,1)=image.at<uchar>(y+10,x-4);
	 result.at<uchar>(3,2)=image.at<uchar>(y+7,x-7);
	 result.at<uchar>(3,3)=image.at<uchar>(y+4,x-10);

	 return result;
}
inline Mat corner_detect5(const size_t h, const size_t w,  Mat image)
{ // Compute response function for image with radius 5
	size_t border=6;
	Mat response = Mat(h,w, CV_32FC1, cvScalar(0));
	for (size_t y = border; y < h - border; y++){
		for (size_t x = border; x < w - border; x++) {
				// Get points on circle radius 5
				// Compute response and save in response Matrix
				response.at<float>(y,x)=computeResponse(getPoints5(image, x, y));

		}}
	return response;
}
inline Mat corner_detect10(const size_t h, const size_t w,  Mat image)
{ // Compute response function for image with radius 10
	size_t border=11;
	Mat response = Mat(h,w, CV_32FC1, cvScalar(0));
		for (size_t y = border; y < h - border; y++)
			for (size_t x = border; x < w - border; x++) {
				// Get points on circle radius 5
				// Compute response and save in response Matrix
				response.at<float>(y,x)=computeResponse(getPoints10(image, x, y));

			}

return response;
}


/// Master functions to detect and update corners in images

vector<Corner> createCornerList(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size PatternSize, Mat cameraMatrix, Mat distCoeffs, bool ShowCorners)
{ /// Do all operations on image to create a sorted list of corner points
  /// Sort them in systematic way

	int CornerPoints=PatternSize.width*PatternSize.height; //number of points to be detected
	vector<Point2f> PotentialCorners=detectPotentialCorners(img, ResponseThreshold, minDistance, detectionRadius); // Detect candidates based on corner measure
	if((int)PotentialCorners.size()!=CornerPoints){ // Number of individual corner points does not match with grid size
			namedWindow( "Display window", WINDOW_NORMAL );// Create a window for display.
			resizeWindow("Display window", 800, 500);
			cout<<"error corners-> not all found: "<<PotentialCorners.size()<<endl;

			for(size_t t=0; t<PotentialCorners.size(); t++){
				cout<<t<<" corner "<<PotentialCorners[t]<<endl;
				circle( img, PotentialCorners[t], 5, Scalar(150, 150, 150), -1, 8, 0 ); // show point located closest
				}
			char k;
			imshow("Display window", img);
			k=waitKey(0);
			destroyWindow("Display window");
	return vector<Corner>(0); // Return empty vector to alert problem
	}
	vector<Point2f> sortedcorners=sortPattern(PotentialCorners, PatternSize); // Sort pattern according to u and v based on patternsize
	if(ShowCorners){ // Show the detected corners if users requires it
		namedWindow( "Display window", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display window", 800, 500);
		for(size_t t=0; t<sortedcorners.size(); t++){
					circle( img, sortedcorners[t], 5, Scalar(150, 150, 150), -1, 8, 0 ); // show point located closest
					circle( img, sortedcorners[t], 2
							, Scalar(0, 0, 0), -1, 8, 0 ); // show point located closest

					}
			char k;
			imshow("Display window", img);
			k=waitKey(0);
		destroyWindow("Display window");
	}

	vector<Point2f> points_undistorted= undistortCorners(sortedcorners, cameraMatrix, distCoeffs ); // Undistort retrieved corner list

	// Create Corner vector from point2f vector
	vector<Corner> cornerList(points_undistorted.size());
	vector<int> iD(points_undistorted.size());
	std::iota (std::begin(iD), std::end(iD), 0);
	transform(points_undistorted.begin(), points_undistorted.end(),iD.begin(), cornerList.begin(), createCorneriD<Corner>());
	return cornerList;
}
vector<Corner> updateCornerlist(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, float MatchesThreshold, Mat cameraMatrix, Mat distCoeffs, bool ShowCorners, vector<Corner> prevCorners){
	/// Update already created corner list prevConers
	/// Use newly detected corners in new image, based on threshold values
	vector<Point2f> PotentialCorners=detectPotentialCorners(img, ResponseThreshold, minDistance, detectionRadius);
	if(ShowCorners){ // Show the detected corners if users requires it
			namedWindow( "Display window", WINDOW_NORMAL );// Create a window for display.
			resizeWindow("Display window", 800, 500);
			for(size_t t=0; t<PotentialCorners.size(); t++){
						circle( img, PotentialCorners[t], 5, Scalar(150, 150, 150), -1, 8, 0 ); // show point located closest
						circle( img, PotentialCorners[t], 2
								, Scalar(0, 0, 0), -1, 8, 0 ); // show point located closest

						}
				char k;
				imshow("Display window", img);
				k=waitKey(0);
			destroyWindow("Display window");
		}
	vector<Point2f> points_undistorted= undistortCorners(PotentialCorners, cameraMatrix, distCoeffs ); // Undistort retrieved corner points

	// Match vector of point coordinates with previous locations
	vector<Corner> nextCorners(prevCorners.size());
	transform(prevCorners.begin(), prevCorners.end(),nextCorners.begin(), findMatches<Corner>(points_undistorted, MatchesThreshold));
	return nextCorners;
}
vector<Corner> readFeaturesImage(string name, CameraParams cam, string OutputName, Settings s, vector<Corner> prevCorners)
{ ///Finds all feature points in image, returns vector of matched corners
  /// If requested write Feature coordinates to file
	if(s.TypeFeatureInput){ // Detect in images

	Mat view = imread(name, IMREAD_GRAYSCALE); //Reading in starts from last image to first to avoid feature loss due to heavy turbulence
	if( !view.data )
		   cout <<"could not load image:"<<name <<endl;
	//else cout<<".."<<endl;

	// Detect features in image and undistort them
	vector<Corner> Corners=updateCornerlist(view, s.ResponseThreshold, s.MinDistance, s.ResponseRadius, s.MatchesThreshold, cam.K, cam.distCoeffs, s.ShowCorners, prevCorners);

	if(s.SaveFeatureCoordinates){ // Save to file if requested

		writeVecToFile(Corners, OutputName);}
	return Corners;
	}
	else{ // Read corners from text-file
	ifstream input;
		vector<Corner> Corners;
		input.open(name, std::ifstream::in);

		//Initialization of temporary variables
		float CoordX, CoordY;
		int iD;
		string line;

		while(std::getline(input, line)){ //Keep reading text-files until end of file

			// Input of camera 1
			istringstream   st(line);
			st>> iD >>CoordX >> CoordY;
			Corners.push_back(Corner(Point3f(CoordX, CoordY, 1),iD)); // Append corner location to list of corner points of that time instance
		}
		return Corners;}
}
vector<Corner> readFeaturesFirstImage(string name, CameraParams cam, string OutputName, Settings s)
{ ///Finds all feature points in image, returns vector of sorted image points
  /// If requested write Feature coordinates to file

	if(s.TypeFeatureInput){// Detect corners in images
	Mat view = imread(name, IMREAD_GRAYSCALE); //Reading in starts from last image to first to avoid feature loss due to heavy turbulence
	if( !view.data )
		   cout <<"could not load image:"<<name <<endl;
	else cout<<".."<<endl;

	// Detect features in image and undistort them
	vector<Corner> Corners= createCornerList(view, s.ResponseThreshold, s.MinDistance, s.ResponseRadius, s.FeaturePatternSize, cam.K, cam.distCoeffs, s.ShowCorners);
	if(s.SaveFeatureCoordinates){ // Write coordinates to text file if requested
		writeVecToFile(Corners, OutputName);}
	return Corners;
	}
	else{ // Read corners from text file
	ifstream input;
	vector<Corner> Corners;
	input.open(name, std::ifstream::in);

	//Initialization of temporary variables
	float CoordX, CoordY;
	int iD;
	string  line;

	while(std::getline(input, line)){ //Keep reading text-files until end of file

		// Input of line
		istringstream   st(line);
		st>> iD >> CoordX >> CoordY;
		Corners.push_back(Corner(Point3f(CoordX, CoordY, 1), iD)); // Append corner location to list of corner points together with iD
	}
	return Corners;}
}

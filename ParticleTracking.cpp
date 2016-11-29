/*
 * PartickeTracking.cpp
 *
 *  Created on: Nov 25, 2016
 *      Author: lengelen
 */
#include "ParticleTracking.h"
#include "ReadAndWrite.h"
/// Definition template classes

template <class Type> class distance_to_Line
{ ///Computes distance between 2D Float point and line (charct by start- and end-point), returns (point.x, point.y, dist)
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
{ ///Computes distance between 2D Float points and one predefined 2D Float point, returns (point.x, point.y, dist)
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
            	Point3f result=Point3f(point2.x,point2.y,distance);
            	return result;
            }
        };
template <class Type> class getResponse
{ // Creates 3D float point with (u,v, Response)
private:
    Mat Resp;
public:
    // Constructor
    getResponse (Mat i_R) {
        Resp=i_R;
    }

    // The function call
    Type operator ( ) ( Point2f& elem ) const
    {	Point3f result =Point3f(elem.x,elem.y, Resp.at<float>(elem));
        return result;
     }
};
template <class Type> class angleImage1
{ ///Compute  angle between two 2D image points in way 1: -pi ...pi
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
    		float ang=atan2(p2.y-p1.y,p2.x-p1.x);
    		return ang;
    	}
    	else return 100; ///2 points are actually the same
    	}

};
template <class Type> class angleImage2
{ ///Compute  angle between two 2D image points in way 2: 0--2pi
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
    	else return 100; ///2 points are actually the same
    	}

};
template <class T> struct reduce_Z
{ ///Throws away distance and only retains pixel coordinates
    T operator() (Point3f& p) const {
        Point2f res = Point2f(p.x,p.y);
        return res;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct createCorneriD
{ ///Creates corner object from iD and image coordinates
   T operator() (Point2f& p, int& iD) const {

		Corner res=Corner(Point3f(p.x,p.y,1),iD);

		return res;}
		typedef T first_argument_type;
		typedef T second_argument_type;
		typedef T result_type;
	};
template <class Type> class angleImageWithPoint
{ ///Compute  angle between two 2D image points in way 1: -pi..pi with point coordinate as well
  ///result is 3D point with 3rd dim angle
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
{
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
            	Point3f predCorner=q.predictPoint();
                  for(size_t i=0;i<candidates.size();i++){
                      	  if((predCorner.x-candidates[i].x)*(predCorner.x-candidates[i].x)+(predCorner.y-candidates[i].y)*(predCorner.y-candidates[i].y)<Threshold){
                		temp=Point3f(candidates[i].x, candidates[i].y, 1);
                		count++;
                	  }}

                  if(count!=1){
                	  q.setFound(false);
                	  return q;}
                  else{
                	  q.setCoords(temp);
                	  q.setFound(true);
                	  return q;
                  }

               }
};
template <class T> struct reduce_Resp
{ ///Throws away Reponse and only retains pixel coordinates
    T operator() (Point3f p) const {
        Point2f res = Point2f(p.x,p.y);
        return res;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
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
{ ///sort image points according to increasing distance param (3rd elem)
 bool operator() ( Point2f a, Point2f b ){
               return pow(a.x,2)+pow(a.y,2) <= pow(b.x,2)+pow(b.y,2);
    }
} sortingOrigin;

/// Simple functions to compare variables etc

inline bool compare(float a, float b)
{ /// Compare-function of two float variables to check for equality (lower than threshold)
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
inline bool point_comparator( cv::Point3f a, cv::Point3f b) { //Comparison function o
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
		Point2f point1=Point2f(extremesSorted[0].x+t*stepX1,extremesSorted[0].y+t*stepY1);//edge point left
		Point2f point2=Point2f(extremesSorted[2].x+t*stepX2,extremesSorted[2].y+t*stepY2);//edge point right

		transform(corners.begin(), corners.end(), distances.begin(), distance_to_Line<Point3f> (point1, point2)); //compute distance to line
		sort(distances.begin(), distances.end(), sortingZ); //sort points accoring to distance
		int k=0;
		vector <Point2f> temp;
		for(size_t j=0; j<corners.size(); j++){
			if(distances[j].z<50 && k<boardSize.height) //threshold of 100 pixels from line and also maximum number of points on one row
				{Point2f P=Point2f(distances[j].x,distances[j].y);
				temp.push_back(P);

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
			  // trasnformation to get correct pixel coordinates-wrong origin by opencv
			  points_undistorted[k]=Point2f(r.at<Vec2f>(k,0)[0] * fx + cx,r.at<Vec2f>(k,0)[1] * fy + cy);
					  }
	return points_undistorted;
}
vector<Point2f> detectPotentialCorners(Mat img, float ResponseThreshold, float minDistance, int detectionRadius){
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
	int number=points.size();
	if(number==0){
		cout<<"No corners found";
		exit(1);
	}

	vector<Point3f> PointsResponse(number);
	vector<Point2f> PointsOrdered(number);

	//Attach response to point coordinates and sort them for decreasing response
	transform(points.begin(), points.end(), PointsResponse.begin(), getResponse<Point3f> (response) );
	sort(PointsResponse.begin(), PointsResponse.end(), point_comparator);

	//Remove Reponse value of points
	transform(PointsResponse.begin(), PointsResponse.end(), PointsOrdered.begin(), reduce_Resp<Point2f>());

	for(size_t i=0; i<PointsOrdered.size()-1;i++){ // For all points until the final-1 (because previous already checked)
		float N=0;

		Point2f avg=PointsOrdered[i]*response.at<float>(PointsOrdered[i]);

		for(size_t t=i+1; t<PointsOrdered.size();t++){ // For all points after considered point
			//Compute distance
			float dist=sqrt((PointsOrdered[i].x-PointsOrdered[t].x)*(PointsOrdered[i].x-PointsOrdered[t].x)+(PointsOrdered[i].y-PointsOrdered[t].y)*(PointsOrdered[i].y-PointsOrdered[t].y));
			//Remove points located closer than mindist and attach to weighted average
			if(dist<mindist){
				avg+=PointsOrdered[t]*response.at<float>(PointsOrdered[t]);
				N=N+response.at<float>(PointsOrdered[t]);
				PointsOrdered.erase(PointsOrdered.begin() + t);
				t--;
				}
		}
		if(N!=0) PointsOrdered[i]=avg/(response.at<float>(PointsOrdered[i])+N); //Take weighetd average

	}
	return PointsOrdered;
}
vector<Point2f> create_points(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, Mat cameraMatrix, Mat distCoeffs, bool showCorners)
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
	if(showCorners){ // Show the detected corners if users requires it
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
{ // Compute reponse value for center pixel based on points located on circle

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

vector<Corner> createCornerList(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size PatternSize, Mat cameraMatrix, Mat distCoeffs, bool showCorners)
{ /// Do all operations on image to obtain a sorted list of corner points

	int CornerPoints=PatternSize.width*PatternSize.height; //number of points to be detected
	vector<Point2f> PotentialCorners=detectPotentialCorners(img, ResponseThreshold, minDistance, detectionRadius);
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
	vector<Point2f> sortedcorners=sortPattern(PotentialCorners, PatternSize);
	if(showCorners){ // Show the detected corners if users requires it
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

	vector<Point2f> points_undistorted= undistortCorners(sortedcorners, cameraMatrix, distCoeffs );
	vector<Corner> cornerList(points_undistorted.size());
	vector<int> iD(points_undistorted.size());
	std::iota (std::begin(iD), std::end(iD), 0);
	transform(points_undistorted.begin(), points_undistorted.end(),iD.begin(), cornerList.begin(), createCorneriD<Corner>());
	return cornerList;
}
vector<Corner> updateCornerlist(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, float MatchesThreshold, Mat cameraMatrix, Mat distCoeffs, bool showCorners, vector<Corner> prevCorners){
	vector<Point2f> PotentialCorners=detectPotentialCorners(img, ResponseThreshold, minDistance, detectionRadius);
	if(showCorners){ // Show the detected corners if users requires it
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
	vector<Point2f> points_undistorted= undistortCorners(PotentialCorners, cameraMatrix, distCoeffs );
	vector<Corner> nextCorners(prevCorners.size());
	transform(prevCorners.begin(), prevCorners.end(),nextCorners.begin(), findMatches<Corner>(points_undistorted, MatchesThreshold));
	return nextCorners;
}
vector<Corner> readFeaturesImage(string name, CameraParams cam, string OutputName, Settings s, vector<Corner> prevCorners)
{ ///Finds all feature points in image, returns vector of sorted imagepoints
	if(s.InputType){
	// Read in image

	Mat view = imread(name, IMREAD_GRAYSCALE); //Reading in starts from last image to first to avoid feature loss due to heavy turbulence
	if( !view.data )
		   cout <<"could not load image:"<<name <<endl;
	//else cout<<".."<<endl;

	// Detect features in image and undistort them
	vector<Corner> Corners=updateCornerlist(view, s.ResponseThreshold, s.MinDistance, s.responseRadius, s.MatchesThreshold, cam.K, cam.distCoeffs, s.showCorners, prevCorners);
	if(s.saveFeatureCoordinates){

		writeVecToFile(Corners, OutputName);}
	return Corners;
	}
	else{
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

	if(s.InputType){
	// Read in image
	Mat view = imread(name, IMREAD_GRAYSCALE); //Reading in starts from last image to first to avoid feature loss due to heavy turbulence
	if( !view.data )
		   cout <<"could not load image:"<<name <<endl;
	else cout<<".."<<endl;

	// Detect features in image and undistort them
	vector<Corner> Corners= createCornerList(view, s.ResponseThreshold, s.MinDistance, s.responseRadius, s.PatternSize, cam.K, cam.distCoeffs, s.showCorners);
	if(s.saveFeatureCoordinates){
		writeVecToFile(Corners, OutputName);}
	return Corners;
	}
	else{
	ifstream input;
	vector<Corner> Corners;
	input.open(name, std::ifstream::in);

	//Initialization of temporary variables
	float CoordX, CoordY;
	int iD;
	string  line;

	while(std::getline(input, line)){ //Keep reading text-files until end of file

		// Input of camera 1
		istringstream   st(line);
		st>> iD >> CoordX >> CoordY;
		Corners.push_back(Corner(Point3f(CoordX, CoordY, 1), iD)); // Append corner location to list of corner points of that time instance
	}
	return Corners;}
}

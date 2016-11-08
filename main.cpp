/// Inlcude header files and external source files
#include "main.h"

/// Self-defined classes and structures (some are adapted versions of free source code online)
/// Some classes are both defined in double and float variant

template <class Type> class worldCoords 
{ ///Determines world coordinates (3D Float) of image point in vector 3-Float format (homog)
private:
    Mat RR, TT, KK;   // inverse matrices
public:
    // Constructor
    worldCoords (Mat i_RR, Mat i_TT, Mat i_KK) {

        RR=i_RR;
        TT=i_TT;
        KK=i_KK;
    }

    // The function call
    Type operator ( ) ( Point3f& elem ) const
    {	Mat q=Mat(elem).reshape(1,3);
        q.convertTo(q, CV_64F);

        Mat result= (RR * (KK*q-TT));
        Mat resultf;
        result.convertTo(resultf, CV_32F);
        Point3f resultt =Point3f(resultf.at<float>(0,0),resultf.at<float>(1,0), resultf.at<float>(2,0));
        return resultt;
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

template <class Type> class pixelCoordinates
{ ///Determines pixel coordinates of 3D point (float), output 3Float (homog)
private:
    Mat R, T, K;   // inverse matrices
public:
    // Constructor
    pixelCoordinates (Mat i_R, Mat i_T, Mat i_K) {

        R=i_R;
        T=i_T;
        K=i_K;
    }

    // The function call
    Type operator ( ) ( Point3f& elem ) const
    {

    	Mat q=Mat(elem).reshape(1,3);
    	q.convertTo(q, CV_64F);
        Mat result= (K * (R*q+T));
        Mat resultf;
        result.convertTo(resultf, CV_32F);
        Point3f resultt= Point3f(resultf.at<float>(0,0)/ resultf.at<float>(2,0),resultf.at<float>(1,0)/ resultf.at<float>(2,0),1.0);
        return resultt;
    }
};
template <class Type> class Ray 
{ ///Computes ray direction (3D Float) between c and image point (not normalized)
        private:
            Point3f c;   // Camera position
        public:
            // Constructor
            Ray (Point3f i_c) {

                c=i_c;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {
                  Point3f result= Point3f(elem.x-c.x,elem.y-c.y, elem.z-c.z);
                  return result;
            }
};
template <class Type> class Normals_2ndorder
{  ///Computes normalized n2 fo set of 8 coeffcieints at position of point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Camera position
        public:
            // Constructor
            Normals_2ndorder (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {

                Mat n= (Mat_<float>(3,1)<< 2*PI_F/Lx*coef[5]*sin(2*PI_F*elem.x/Lx)+PI_F/Lx*coef[7]*sin(PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[1]*sin(PI_F*elem.x/Lx)-coef[3]/Lx, 2*PI_F/Ly*coef[6]*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[7]*cos(PI_F*elem.x/Lx)*sin(PI_F*elem.y/Ly)+PI_F/Ly*coef[2]*sin(PI_F*elem.y/Ly)-coef[4]/Ly, 1);
                float nn =norm(n, NORM_L2);
                Point3f res= Point3f(n.at<float>(0,0)/nn, n.at<float>(1,0)/nn, n.at<float>(2,0)/nn);

                return res;
            }
        };
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
template <class Type> class Angle 
{ ///Computes thetai for given thetadelta and rw
    private:
        double rw;

    public:
            // Constructor
            Angle (double i_rw) {
                rw=i_rw;
            }

            // The function call
            Type operator ( ) ( double& theta ) const
            {
                double res= atan((rw*sin(theta))/(rw*cos(theta)-1));
                return res;
            }
        };
template <class Type> class Angle2 
{ ///Computes thetadelta for given thetai and rw
    private:
        double rw;

    public:
            // Constructor
            Angle2 (double i_rw) {
                rw=i_rw;
            }

            // The function call
            Type operator ( ) ( double& theta ) const
            {

                double res=	theta - asin(sin(theta)/rw);
                return res;
            }
        };
template <class Type> class waterSurface_2ndorder 
{ ///Finds surface point for set of 8 coefficients and ray u=cq' (all 3D Float)

        private:
		real_1d_array coef; Point3f c;  double Lx; double Ly;

        public:
            // Constructor
		waterSurface_2ndorder (Point3f i_c, double i_Lx, double i_Ly, real_1d_array i_coef) {
            	coef=i_coef;
                c=i_c;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& ray ) const
            {

                double *fvec;
                int iflag=1;
                int lwa;
                int n = 4;
                double tol = 0.00001;
                double *wa;
                double *x;
                double *parameters;

                lwa = ( n * ( 3 * n + 13 ) ) / 2;

                fvec = new double[n];

                wa = new double[lwa];
                x = new double[n];
                parameters = new double [16];

                x[0] = 1500;
                x[1] = 0;
                x[2] = 0;
                x[3] = 0;

                parameters[0] = c.x;
                parameters[1] = c.y;
                parameters[2] = c.z;
                parameters[3] = ray.x;
                parameters[4] = ray.y;
                parameters[5] = ray.z;
                parameters[6] = Lx;
                parameters[7] = Ly;
                parameters[8] = coef[0];
                parameters[9] = coef[1];
                parameters[10] = coef[2];
                parameters[11] = coef[3];
                parameters[12] = coef[4];
				parameters[13] = coef[5];
				parameters[14] = coef[6];
				parameters[15] = coef[7];

                findIntersection_2ndorder( n, x, fvec, &iflag, parameters );

                int info = hybrd1 ( findIntersection_2ndorder, n, x, fvec, parameters, tol, wa, lwa );
                Point3f result =Point3f(x[1],x[2],x[3]);
                delete[] fvec;
                delete[] wa;
                delete[] x;
                delete[] parameters;
                return result;
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
{ ///sort image points according to increasing distance param (3rd elem)
 bool operator() ( Point2f a, Point2f b ){
               return pow(a.x,2)+pow(a.y,2) <= pow(b.x,2)+pow(b.y,2);
    }
} sortingOrigin;
template <class T> struct prod 
{ /// computes angle between two vectors by dot product ->for thetad
    T operator() (Point3f& uu, Point3f& vv) const {

    	double res= acos(uu.ddot(vv));

    	return res;}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct prod2 
{/// computes angle between two vectors by dot product but reverses one vector in direction -> for u and n (thetai)
    T operator() (Point3f& u, Point3f& n) const {

    	Point3f uu=Point3f(-u.x,-u.y,-u.z);
    	double res= acos(uu.ddot(n));

    	return res;}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct error_col 
{ ///Computes normal collinearity metric
    T operator() (Point3f& n1, Point3f& n2) const
    {
    	double e= acos(n1.ddot(n2));

    	return e;
    }
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Waterray 
{ ///computes intersection of ray v2 with plane F
    T operator() (float& a, Point3f& v) const
    {
    	Point3f v2= Point3f(a*v.x,a*v.y,a*v.z);

    	return v2;
    }
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct axis 
{ ///computes rotation axis by cross product of û and ^v (order of vectors is important!)
	        T operator() (Point3f& uu, Point3f& vv) const {
        	Point3f res =vv.cross(uu);

        return res;}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Homogenous 
{ ///transforms 2D image point in homogoenous coordinates
        T operator() (Point2f& q) const {return Point3f(q.x,q.y,1);}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct takeNorm 
{ ///computes norm 3D Float image point  (OpenCV gave troubles)
        T operator() (Point3f& p) const {
        	double l =sqrt( p.x*p.x + p.y*p.y + p.z*p.z );
        	return l;}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Rotationmatrix 
{ ///computes rotationmatrix around axis -ax and over theta 
    T operator() (double& theta, Point3f& ax) const {

    	Mat_<double> r=Mat(-ax).reshape(1,3);
        double naxis =norm(r, NORM_L2);
        r=r/naxis;
        double c = cos(theta);
        double s = sin(theta);

        Mat R= (Mat_<double>(3,3)
        << c + (1-c)*r.at<double>(0,0) * r.at<double>(0,0), (1-c) *r.at<double>(0,0) * r.at<double>(1,0) - s*r.at<double>(2,0), (1-c)*r.at<double>(0,0)*r.at<double>(2,0) + s*r.at<double>(1,0),
        (1-c)*r.at<double>(0,0)*r.at<double>(1,0)+s*r.at<double>(2,0), c+(1-c) * r.at<double>(1,0) * r.at<double>(1,0), (1-c)*r.at<double>(2)*r.at<double>(1,0)-s*r.at<double>(0,0),
        (1-c)*r.at<double>(0,0)*r.at<double>(2,0)-s*r.at<double>(1,0), (1-c)*r.at<double>(2,0)*r.at<double>(1,0)+s*r.at<double>(0,0), c+(1-c)*r.at<double>(2,0)* r.at<double>(2,0));

        Mat res;
        R.convertTo(res, CV_32F);
        return res;}
    typedef T first_argument_type;
    typedef T second_argument_type;
    typedef T result_type;
};
template <class T> struct Rotationmin 
{ ///rotates -u with rotationmatrix
      T operator() (Mat& R, Point3f& u) const {

          Mat uu=Mat(-u).reshape(1,3);
          Mat n2=R*uu;
          float nn =norm(n2, NORM_L2);

          Point3f res =Point3f(n2.at<float>(0,0)/nn, n2.at<float>(1,0)/nn, n2.at<float>(2,0)/nn);
          return res;}

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Rotation 
{///rotates u with rotationmatrix
      T operator() (Mat& R, Point3f& u) const {

          Mat uu=Mat(u).reshape(1,3);
          Mat n2=R*uu;
          float nn =norm(n2, NORM_L2);

          Point3f res =Point3f(n2.at<float>(0,0)/nn, n2.at<float>(1,0)/nn, n2.at<float>(2,0)/nn);
          return res;}

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Distance 
{ ///computes distance along v from surface point w towards plane z=0
	    T operator() (Point3f& w, Point3f& v) const {

          float res =-w.z/v.z;
          return res;}

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct SubtractNormalized 
{ ///subtract two 3D Float vectors and gives normalized result
    T operator() (Point3f& v1,Point3f& v2) const {

        Mat res= (Mat_<float>(3,1)<< v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);

        float n =norm(res, NORM_L2);

        Point3f ress = Point3f(res.at<float>(0,0)/n, res.at<float>(1,0)/n, res.at<float>(2,0)/n);
        return ress;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Subtract 
{  ///subtract two 3D Float vectors 
    T operator() (Point3f& v1,Point3f& v2) const {
        Point3f res = Point3f(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
        return res;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct reduce_Z
{ ///Throws away distance and only retains pixel coordinates
    T operator() (Point3f p) const {
        Point2f res = Point2f(p.x,p.y);
        return res;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct AbsoluteAngle
{ ///Takes 3D point with 3rd dim angle and gives abs and returns same 3D point with abs value of angle
    T operator() (Point3f p) const {
        Point3f res = Point3f(p.x,p.y, abs(p.z));
        return res;
    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
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

//Random functions used
bool point_comparator( cv::Point3f a, cv::Point3f b) { //Comparison function o
    return (b.z < a.z);
}
int nancheck(double x) { return x != x; } //check if nan - double version
int nancheck2(float x) { return x != x; } //check if nan - float version
bool compare(float a, float b)
{ /// Compair-function of two float variables to check for equality (lower than threshold)
    if(fabs(a - b) < (1.0 / 100000))
        return true;
    else
        return false;
}
bool operator==(const Point2f& pt1, const Point2f& pt2)
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
void removedupes(std::vector<Point2f> & vec)
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

//Find intersection of rays with with surfaces
void findIntersection( int n, double x[], double fvec[], int *iflag, double params[] )
{///Find intersection of ray cq' with surface charct by 5 coeffs
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double A10 = params[9];
    double A01 = params[10];
    double B = params[11];
    double C = params[12];


    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*x[1]/Lx+C*x[2]/Ly;
    return;
}
void findIntersection_2ndorder( int n, double x[], double fvec[], int *iflag, double params[] )
{///Find intersection of ray cq' with surface charct by 8 coeffs
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double A10 = params[9];
    double A01 = params[10];
    double B = params[11];
    double C = params[12];
	double A20 = params[13];
	double A02 = params[14];
	double A11 = params[15];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*x[1]/Lx+C*x[2]/Ly + A20 * cos(2*PI*x[1]/Lx)+ A02 * cos(2*PI*x[2]/Ly) + A11 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly);
    return;
}

//Reading functions
static bool readStringList(vector<string>& l, string dir )
{///finds all images in directory and stores their names in vector of strings
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
static inline void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
{///reads node in settings file and stores it in variable
    if(node.empty())
        x = default_value;
    else
        x.read(node);
}

// Write functions
void writeMatToFile(cv::Mat& m, const char* filename, string  outputDirectory)
{ ///write a Float-matrix to outputdirectory, using filename
    std::ofstream fout(outputDirectory+filename);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }

    for(int i=0; i<m.rows; i++)
    {
        for(int j=0; j<m.cols; j++)
        {
            fout<<m.at<float>(i,j)<<"\t";
        }
        fout<<endl;
    }

    fout.close();
}
void writeVecToFile(vector<Point3f> features, string  file, size_t timestep)
{///write a Float-Vector of 3D feature points to file with path name out
	std::ofstream out;

	 // std::ios::app is the open mode "append" meaning
	 // new data will be written to the end of the file.
	out.open(file, std::ios_base::app);

    if(!out)
    {
        cout<<"File Not Opened"<<endl;  return;
    }

    for(size_t i=0; i<features.size(); i++)
    {
              out<< timestep<<"\t"<<i<<"\t"<<features[i].x<<"\t"<<features[i].y<<"\t"<<features[i].z<<"\t";
              out <<endl;
    }
    out.close();
}
void writeArray(vector <real_1d_array> plane, string  out)
{///write an Array of surface coefficients to file with path name out
    ofstream fout(out);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }


    for(size_t i=0; i<plane.size(); i++)
    {
    		fout <<plane[i].tostring(2).c_str()<<endl;

    }

    fout.close();
}
void writeArrayErrors(vector <double> errors, string  out)
{///write an Array of surface coefficients to file with path name out
    ofstream fout(out);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }


    for(size_t i=0; i<errors.size(); i++)
    {
    		fout << errors[i]<<endl;

    }

    fout.close();
}
static void saveCameraParams( const string& filename, Mat Rotationmatrix, Mat TranslationVector)
{ ///Saves the results into file with filename
    FileStorage fs( filename, FileStorage::WRITE ); //Open filestorage file

    //Writing out the results
    fs << "Rotationmatrix" << Rotationmatrix;
    fs << "TranslationVector" << TranslationVector;
    fs.release();
	}

//Functions to detect corners
float computeResponse(Mat points)
{ // Compute reponse value for center pixel based on points located on circle

	float SR=0;
	float DR=0;
	//float NM=0; // Old, original definitioon of reponse function
	//float MR=0;
	for(size_t i=0; i<4;i++){
		float temp=abs((points.at<uchar>(0,i)+points.at<uchar>(2,i))-(points.at<uchar>(1,i)+points.at<uchar>(3,i)));
		float temp2=abs(points.at<uchar>(0,i)-points.at<uchar>(2,i))+abs(points.at<uchar>(1,i)-points.at<uchar>(3,i));
		//float temp3=(a[i]+c[i]+b[i]+d[i]);
		SR=SR+temp;
		DR=DR+temp2;
		//NM=NM+temp3;
	}
	//NM=NM/16;
	//MR=abs(NM-LM);
	return (SR-DR);//-16*MR);
}
Mat getPoints5(Mat image, size_t x, size_t y)
{ // Returns Mat 4x4 with inetnsity value at 16 points in circle with radius 5

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
Mat getPoints10(Mat image, size_t x, size_t y)
{// Returns Mat 4x4 with inetnsity value at 16 points in circle with radius 10

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
Mat corner_detect5(const size_t h, const size_t w,  Mat image)
{ // Compute response function for image with radius 5
	size_t border=6;
	Mat response = Mat(h,w, CV_32FC1, cvScalar(0));

	for (size_t y = border; y < h - border; y++)
			for (size_t x = border; x < w - border; x++) {

				Mat points=getPoints5(image, x, y); // Get points on circle radius 5

				 //float LM=(image.at<uchar>(y,x)+image.at<uchar>(y+1,x)+image.at<uchar>(y,x+1)+image.at<uchar>(y-1,x)+image.at<uchar>(y,x-1))/5; // Old, original definitioon of reponse function

				// Compute response and save in response Matrix
				float R=computeResponse(points);
				response.at<float>(y,x) = R;

		}
return response;
}
Mat corner_detect10(const size_t h, const size_t w,  Mat image)
{ // Compute response function for image with radius 10
	size_t border=11;
	Mat response = Mat(h,w, CV_32FC1, cvScalar(0));

	for (size_t y = border; y < h - border; y++)
			for (size_t x = border; x < w - border; x++) {

				Mat points=getPoints10(image, x, y); // Get points on circle radius 5

				 //float LM=(image.at<uchar>(y,x)+image.at<uchar>(y+1,x)+image.at<uchar>(y,x+1)+image.at<uchar>(y-1,x)+image.at<uchar>(y,x-1))/5; // Old, original definitioon of reponse function

				// Compute response and save in response Matrix
				float R=computeResponse(points);
				response.at<float>(y,x) = R;

			}
return response;
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

	//Attach response to point coordinates and sort them for decreasing responnse
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

//Processing images & detecting feature points
vector<Point2f> sortGrid(vector<Point2f> corners, Size boardSize)
{ /// Sort grid of corners points autamatically in ordered pattern
  /// Columns in streamwise direction
  /// Rows in crosswise direction
	vector<Point3f> delta1;
	vector<Point3f> delta2;
	vector<Point2f> extremes(4);
	vector<Point3f> anglesOrigin(3);
	vector<Point2f> PointsNoOrigin(3);
	vector<Point3f> distances(corners.size());
	vector<Point2f> sortedcorners;


	float stepX1, stepX2, stepY1, stepY2;

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
		delta1.push_back(Point3f(corners[t].x, corners[t].y, *temp1.second-*temp1.first));
		auto temp2= minmax_element(anglesImage2.begin(),anglesImage2.end()); 
		delta2.push_back(Point3f(corners[t].x, corners[t].y, *temp2.second-*temp2.first));

	}
	// Sort them according to angle-measure: z-coordinate, with smallest first
	sort(delta1.begin(),delta1.end(), sortingZ);
	sort(delta2.begin(),delta2.end(), sortingZ);
	// Only retain best 4
	delta1.resize(4); 
	delta2.resize(4);

	// Ideallu angle measure is Pi/2, so remove ones that deviate more than 3/4 Pi
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
vector<Point2f> create_points(Mat img, float ResponseThreshold, float minDistance, int detectionRadius, Size boardSize, bool distortion, Mat cameraMatrix, Mat distCoeffs, bool showCorners)
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
	vector<Point2f> sortedcorners=sortGrid(corners2, boardSize);
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

	if(distortion){ //Undistort detected image points of sorted corners if required
	Mat r;
	undistortPoints(Mat(sortedcorners), r, cameraMatrix, distCoeffs); // undistort sorted image points
	vector<Point2f> points_undistorted;

	float fx = cameraMatrix.at<double>(0, 0);
	float fy = cameraMatrix.at<double>(1, 1);
	float cx = cameraMatrix.at<double>(0, 2);
	float cy = cameraMatrix.at<double>(1, 2);

	for(int k=0; k<r.rows; k++){
		  // trasnformation to get correct pixel coordinates-wrong origin by opencv
		  points_undistorted.push_back(Point2f(r.at<Vec2f>(k,0)[0] * fx + cx,r.at<Vec2f>(k,0)[1] * fy + cy));
		          }

	return points_undistorted;}
	else{ return sortedcorners;}
}
vector<Point3f> readFeaturesImage(vector<string> imageList, size_t frame, size_t nimages, CameraParams cam, Settings s, string fileName="/home/lengelen/Desktop/test/coordinates.txt")
{ ///Finds all feature points in image, returns vector of sorted imagepoints

	// Read in image
	Mat view = imread(imageList[nimages-frame-1], IMREAD_GRAYSCALE); //Reading in starts from last image to first to avoid feature loss due to large turbulences
	if( !view.data )
		   cout <<"could not load image:"<<imageList[nimages-frame-1] <<endl;
	else cout<<".."<<endl;

	// Detect features in image and undistort them
	vector <Point2f> points2D=create_points(view, s.ResponseThreshold, s.MinDistance, s.responseRadius, s.boardSize, true, cam.K, cam.distCoeffs, s.showCorners);
	
	// If number of detected and sorted points don't correspond with pattern size give empty vector to indicate problem
	if(static_cast<int>(points2D.size())!=static_cast<int>(s.boardSize.width*s.boardSize.height))
		return vector<Point3f>(0);

	vector<Point3f> pointsHomog(points2D.size());
	transform (points2D.begin(), points2D.end(), pointsHomog.begin(), Homogenous<Point3f>()); //convert to homogeneous coordinates
	//	writeVecToFile(pointsHomog, fileName, frame ); // Write image points to file, used for verification/debugging
	return pointsHomog;
}
CameraParams Initialize(Settings s, int camera)
{ //Initialization of camera parameters in a CameraParams-var for input of settings and camera number
		//Read all settings
	    Mat cameraMatrix, distCoeffs;
	    string CalibrationFile,  InputReference, OutputFile;
	    Mat R, origin, RR, rvec, tvec;
	    int Position;
	    Point3f c;
	    switch(camera){ // Open calibration file and reference image depending on camera
	    	case 1:
	    		 CalibrationFile= s.CalibrationFile1;
	    	     InputReference= s.InputReference1;
	    	     Position=s.Position1;
	    	     OutputFile=s.OutputCameraPose1;
	    	     break;
	    	case 2:
				 CalibrationFile= s.CalibrationFile2;
				 InputReference= s.InputReference2;
	    	     Position=s.Position2;
	    	     OutputFile=s.OutputCameraPose2;
				 break;
	    	case 3:
				 CalibrationFile= s.CalibrationFile3;
				 InputReference= s.InputReference3;
	    	     Position=s.Position3;
	    	     OutputFile=s.OutputCameraPose3;
				 break;
	    }

	    Size board_sz = s.boardSize;
	    float squareSize = s.squareSize;
	    float ResponseThreshold=s.ResponseThreshold;
	    float MinDistance=s.MinDistance;

	    //Read in calibration information:distortion coeffcients and cameramatrix
	    FileStorage fs2(CalibrationFile,FileStorage::READ);
	        if (!fs2.isOpened())
	            {
	                cout << "Could not open the calibration file: \"" << CalibrationFile << "\"" << endl;

	            }
	     cout <<"Calibration file loaded: "<<CalibrationFile <<endl;
	     fs2["distortion_coefficients"] >> distCoeffs;
	     fs2["camera_matrix"] >> cameraMatrix;
	     fs2.release();

	   	// Computation extrinsic parameters
		// Read in reference image for camera position and orientation

		 if(!s.InputTypeRef){
				FileStorage fs3(InputReference,FileStorage::READ);
				 if (!fs3.isOpened())
								{
									cout << "Could not open the Reference file: \"" << InputReference << "\"" << endl;

								}
				 cout <<"Calibration file loaded: "<<InputReference <<endl;
				 fs3["Rotationmatrix"] >> R;
				 fs3["TranslationVector"] >> tvec;
				 fs3.release();
				 RR = R.inv();  // rotation of inverse
				 origin = -RR * tvec; // translation of inverse
				 c=Point3f(origin.at<double>(0,0), origin.at<double>(1,0), origin.at<double>(2,0)); // Camera center
		 }
		 else{
				Mat image=imread(InputReference, IMREAD_GRAYSCALE);
		        if( image.data )
		        cout <<"extra reference image loaded:"<<InputReference <<endl;
		        else exit(1);

		        //Detect corner points in image but don't undistort them
		        vector<Point2f> sortedcorners =create_points(image, ResponseThreshold, MinDistance,s.responseRadius, board_sz, false, cameraMatrix, distCoeffs, s.showCorners);

		        Mat r;
				Mat undistortedimage;
				undistort(image, undistortedimage, cameraMatrix, distCoeffs);//undistort image

		        vector<Point3f> objectPoints(sortedcorners.size());
		        calcCorners(board_sz, squareSize, objectPoints, Position);//theoretical 3D world coordinates
		        bool succes= solvePnP(objectPoints,  sortedcorners,  cameraMatrix, distCoeffs,  rvec, tvec, false); //camera pose estimation
		               if (!succes){

		                   cout<< "Initialization not succeeded" << endl;
		                   exit(1);}
			   Rodrigues(rvec, R);
			   RR = R.inv();  // rotation of inverse
			   origin = -RR * tvec; // translation of inverse
			   c=Point3f(origin.at<double>(0,0), origin.at<double>(1,0), origin.at<double>(2,0)); // Camera center

		 }
		ifstream input;
		input.open(s.InputInitial);

		while(input.fail())
		{
			cout<< "File for initial location of f incorrect" << endl;
			exit(1); ;
		}
		vector<Point3f> f;
		Point3f tmp;
		while (input >> tmp.x && input >> tmp.y && input >> tmp.z)
		{
			f.push_back(tmp);
		};
		if(s.saveCameraPose){

		saveCameraParams(OutputFile, R, tvec);

		}
		CameraParams cam=CameraParams(R, tvec, cameraMatrix, c, distCoeffs, f);
        return cam;
}
void calcCorners(Size boardSize, float squareSize, vector<Point3f>& corners, int position)
{ ///Calculate theoretical location of 3D points of camera pose pattern
    corners.resize(0); //set to empty vector
    if(position==1 || position==4){ //cameras at one side

            for( int i = 0; i < boardSize.width; i++ ){
            	for( int j = 0; j < boardSize.height; j++ ){
                corners.push_back(Point3f(float(i*squareSize),
                                          float(j*squareSize), 0));
                }
            }
    }
     else if(position==2 || position==3){ // cameras at other side-->they see pattern reversed
         for( int i = boardSize.width-1; i>-1; i-- ){
        	for( int j = boardSize.height-1; j > -1; j-- ){
            corners.push_back(Point3f(float(i*squareSize),
                                      float(j*squareSize), 0));
        	}
         }
    }
    else cout<<"Where is the camera positioned???"<<endl;
}

//Optimization of coefficients
vector<double> compute_error(float Lx, float Ly, int ErrorMetric, CameraParams Camera, vector<Point3f> Pixels, real_1d_array Coeff)
{ /// Compute errors Ef for all feature points of one camera


	//Initiliaze all temp vectors
	vector<Point3f> qworld(Pixels.size());
	vector<Point3f> q(Pixels.size());
	vector<Point3f> qImagePlane(Pixels.size());
    vector<Point3f> rays(Pixels.size());
    vector<Point3f> water(Pixels.size());
    vector<Point3f> normals(Pixels.size());
    vector<Point3f> u(Pixels.size());
    vector<Point3f> u2(Pixels.size());
    vector<Point3f> v(Pixels.size());
    vector<Point3f> v2(Pixels.size());
    vector<Point3f> X(Pixels.size());
    vector<Point3f> X2(Pixels.size());
    vector<Point3f> X3(Pixels.size());
    vector<Point3f> normals2(Pixels.size());

    vector<double> thetad(Pixels.size());
    vector<double> E(Pixels.size());
    vector<double> thetai(Pixels.size());
    vector<double> thetai2(Pixels.size());
    vector<double> thetai3(Pixels.size());
    vector<double> thetad2(Pixels.size());
    vector<double> thetar(Pixels.size());
    vector<float> a(Pixels.size());

    vector<Mat> Rotations(Pixels.size());
    vector<Mat> Rotations2(Pixels.size());
    vector<Mat> Rotations3(Pixels.size());
    vector<Point3f> features2(Pixels.size());
    vector<Point3f> uu(Pixels.size());
    vector<Point3f> vv(Pixels.size());
    vector<Point3f> shift(Pixels.size());
    vector<Point3f> shiftImage(Pixels.size());


    Mat RR=Camera.R.inv();
    Mat KK=Camera.K.inv();

    transform(Pixels.begin(), Pixels.end(), qworld.begin(), worldCoords<Point3f> (RR, Camera.T, KK)); //compute world coorinates q'
    transform(qworld.begin(), qworld.end(), rays.begin(), Ray<Point3f> (Camera.c)); //compute rays u=cq'

	//Determine surface points according to surface model and hypthesized Coefff
     transform(rays.begin(), rays.end(), water.begin(),  waterSurface_2ndorder<Point3f> (Camera.c, Lx, Ly, Coeff));

	//Compute normals n2 according to surface model and hypthesized Coefff
    transform(water.begin(), water.end(), normals.begin(), Normals_2ndorder<Point3f> (Coeff, Lx,Ly));

    transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>()); // u in alterative way-->gives the same, I checked

    transform (Camera.f.begin(), Camera.f.end(), water.begin(), v.begin(), SubtractNormalized<Point3f>()); //compute raus v=pf

    transform (u.begin(), u.end(), v.begin(), thetad.begin(),prod<double>()); //compute thetad

    transform (u.begin(), u.end(), v.begin(), X.begin(),axis<Point3f> ()); //compute rotation axis

    transform(thetad.begin(), thetad.end(), thetai.begin(), Angle<double>(rw)); //compute thetai

    transform (thetai.begin(), thetai.end(), X.begin(), Rotations.begin(), Rotationmatrix<Mat>()); //detemine rotationmatrix

    transform (Rotations.begin(), Rotations.end(), u.begin(), normals2.begin(), Rotationmin<Point3f>()); //rotate u with rotationmatrix
    if(ErrorMetric==1) //compute normal colinearity metrix
    transform (normals.begin(), normals.end(),normals2.begin(), E.begin(), error_col<double>());
    else if(ErrorMetric==2){ //compute disparity difference metric
    transform (u.begin(), u.end(), normals.begin(), thetai2.begin(),prod2<double>());
    transform(thetai2.begin(), thetai2.end(), thetad2.begin(), Angle2<double>(rw));
    transform (normals2.begin(), normals2.end(), u.begin(), X2.begin(),axis<Point3f> ());
    transform (thetad2.begin(), thetad2.end(), X2.begin(), Rotations2.begin(), Rotationmatrix<Mat>());
    transform (Rotations2.begin(), Rotations2.end(), u.begin(), v2.begin(), Rotationmin<Point3f>());
    transform (water.begin(), water.end(), v2.begin(), a.begin(), Distance<float>());
    transform (a.begin(), a.end(), v2.begin(), vv.begin(), Waterray<Point3f>());
    transform (water.begin(), water.end(), vv.begin(), features2.begin(), std::plus<Point3f>());
    transform (Camera.f.begin(), Camera.f.end(), features2.begin(), shift.begin(), Subtract<Point3f>());
    transform (shift.begin(), shift.end(),E.begin(), takeNorm<double>());}

     return E;
}
void combined_error3(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of all cameras combined for set of 3 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[1], Frame.Pixels[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }
	//Errors Ef third camera
	vector<double> E3 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[2], Frame.Pixels[2], x);
	for(size_t k=0; k<E3.size(); k++){
		    fi[k+E1.size()+E2.size()] = (nancheck(E3[k]) ? 0 : E3[k]);
		 }

 }
void combined_error2(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of all cameras combined for set of 2 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[1], Frame.Pixels[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }

 }
void combined_error1(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of single camera 

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));
	
	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization
 }
real_1d_array optimizeCoef(real_1d_array Coeff, frameOptimizer Frame)
{ /// Finds optimal coefficicients for frame
  /// Requires initial guess and FrameOptimizer containing all necessary input

		//Optimization configuration
		minlmstate state;
	    minlmreport rep;
	    void *ptr;
	    ptr=&Frame;
		//assign parameters for optimization

		minlmcreatev(Frame.CameraAmount*Frame.Pixels[0].size(), Coeff, Frame.diffStep, state);
		minlmsetcond(state, Frame.epsg, Frame.epsf, Frame.epsx, Frame.maxits);

		if(Frame.Params==1){ //only A00
			real_1d_array bndl = "[-0.05,0.0,0.0,0.0,0.0,0,0,0]"; //lower boundary
			real_1d_array bndu = "[150.05,0,0,0.0,0.0,0,0,0]"; //upper boundary
			alglib::minlmsetbc(state,  bndl, bndu);
			switch(Frame.CameraAmount){ //Determine camera-amount and use appropriate function to compute errors
			case 1:
				alglib::minlmoptimize(state, combined_error1, NULL, ptr); //optimize
				break;
			case 2:
				alglib::minlmoptimize(state, combined_error2, NULL, ptr); //optimize
				break;
			case 3:
				alglib::minlmoptimize(state, combined_error3, NULL, ptr); //optimize
				break;
			default:
				cout<<"Error: Impossible amount of camera's"<<endl;}
		}
		else if(Frame.Params==3){ //only A00, B, C
			real_1d_array bndl = "[-0.05,0.0,0.0,-20.0,-20.0,0,0,0]"; //lower boudnary
			real_1d_array bndu = "[150.05,0,0,20.0,20.0,0,0,0]"; //upper boundary
			alglib::minlmsetbc(state,  bndl, bndu);

			switch(Frame.CameraAmount){ //Determine camera-amount and use appropriate function to compute errors
			case 1:
				alglib::minlmoptimize(state, combined_error1, NULL, ptr); //optimize
				break;
			case 2:
				alglib::minlmoptimize(state, combined_error2, NULL, ptr); //optimize
				break;
			case 3:
				alglib::minlmoptimize(state, combined_error3, NULL, ptr); //optimize
				break;
			default:
				cout<<"Error: Impossible amount of camera's"<<endl;}
		}
		else{ //all 8 coefss
			switch(Frame.CameraAmount){ //Determine camera-amount and use appropriate function to compute errors
			case 1:
				alglib::minlmoptimize(state, combined_error1, NULL, ptr); //optimize
				break;
			case 2:
				alglib::minlmoptimize(state, combined_error2, NULL, ptr); //optimize
				break;
			case 3:
				alglib::minlmoptimize(state, combined_error3, NULL, ptr); //optimize
				break;
			default:
				cout<<"Error: Impossible amount of camera's"<<endl;}}

	    //Print termination code
	    minlmresults(state, Coeff, rep);
		cout <<"TerminationCode: "<< int(rep.terminationtype) <<endl;

	    return Coeff;

	}

//Temporal reconstruction functions
int main ( )
{ /// Main function to reconstruct time-dependent water surface.
	
	//Read all settings
	Settings s;
	FileStorage fs("settings.xml", FileStorage::READ); // Read the settings
	if (!fs.isOpened())
		{
			cout << "Could not open the configuration file: \""<< endl;
			exit(1);
		 }
	fs["Settings"] >> s;
	fs.release();   // close Settings file
	CameraParams cam1, cam2, cam3;

	switch(s.CameraAmount){ // Choice of cameras and intialize each camera seperately

	case 1:
		cam1=Initialize(s, 1);
		break;
	case 2:
		cam1=Initialize(s, 1);
		cam2=Initialize(s, 2);
		break;
	case 3:
		cam1=Initialize(s, 1);
		cam2=Initialize(s, 2);
		cam3=Initialize(s, 3);
		break;
	default:
		cout<<"Invalid amount of camera's"<<endl;
		return 0;
	}

	vector <real_1d_array> CoefficientsList; //Output list of coeffcients
	vector <double> ErrorList; //Output list of corresponding errors

	real_1d_array InitialCoeffs;
	if(s.setInitialGuess){ //Change initial guess of coefficients-> influences convergence rate
		double i[8];

		cout << "Please enter a value for A00: ";
		cin >> i[0];
		cout << "Please enter a value for A10: ";
		cin >> i[1];
		cout << "Please enter a value for A20: ";
		cin >> i[5];
		cout << "Please enter a value for A01: ";
		cin >> i[2];
		cout << "Please enter a value for A02: ";
		cin >> i[6];
		cout << "Please enter a value for A11: ";
		cin >> i[7];
		cout << "Please enter a value for B: ";
		cin >> i[3];
		cout << "Please enter a value for C: ";
		cin >> i[4];
		InitialCoeffs.setcontent(8, i);

	}
	else{
	InitialCoeffs= "[5.05,0.0013,0.001,0.00014,-0.0022,0.0007,-0.0007,0.0007]";
	}
	vector<CameraParams> cams;
	frameOptimizer optim;
	switch(s.CameraAmount){ //Choice of cameras and intialize optimalisation with appropriate camera amount

	case 1:
		cams.push_back(cam1);
		optim=frameOptimizer(1, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
		break;
	case 2:
		cams.push_back(cam1);
		cams.push_back(cam2);
		optim=frameOptimizer(2, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
		break;
	case 3:
		cams.push_back(cam1);
		cams.push_back(cam2);
		cams.push_back(cam3);
		optim=frameOptimizer(3, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
		break;
	default:
		cout<<"Invalid amount of camera's"<<endl;
		return 0;
	}

	if(s.Lx!=((s.boardSize.width-1)*s.squareSize) || s.Ly!=((s.boardSize.height-1)*s.squareSize)){ //Verification of lengths scales
		cout<<"Sure about length scales?? [y/n]"<<endl;
		char in;
		cin >> in;
		if(in=='n'){
			float Lx, Ly;
			cout<<"New value Lx:"<<endl;
			cin >> Lx;
			cout<<"New value Ly:"<<endl;
			cin >> Ly;
			s.changeLengthscales(Lx,Ly);
		}
		else if(in!='y')
			return 0;

	}

	if(s.setOptim){ //Change optimization settings
		double i_epsg ; double i_epsf; double i_epsx; ae_int_t i_maxits;double i_diffStep ;
		cout << "Please enter a value for gradient-based stopping condition (epsg): ";
		cin >> i_epsg;
		cout << "Please enter a value for function-based stopping condition (epsf): ";
		cin >> i_epsf;
		cout << "Please enter a value for step-size based stopping condition (epsx): ";
		cin >> i_epsx;
		cout << "Please enter a value for maximum number of repititions (maxits): ";
		cin >> i_maxits;
		optim.setStoppingConditions( i_epsg ,  i_epsf,  i_epsx,  i_maxits);

		cout << "Please enter a value for differentiation step (diffStep): ";
		cin >> i_diffStep;
		optim.setDiffStep( i_diffStep );
	}

	if(s.InputType){ //Determine corner locations from images
		size_t nimages;
		string InputDirectory1, InputDirectory2, InputDirectory3;
		vector<string> imageList1, imageList2, imageList3;
		// Checks to determine if imagenames were read
		bool check=false;
		bool check1=false;
		bool check2=false;
		bool check3=false;

		switch(s.CameraAmount){ // Choice of cameras and read imagenames in every directory (for each camera)

		case 1:
			InputDirectory1=s.InputDirectory1;
			check1=readStringList(imageList1, InputDirectory1); //read all imagenames
			check=check1;
			break;
		case 2:
			InputDirectory1=s.InputDirectory1;
			InputDirectory2=s.InputDirectory2;
			check1=readStringList(imageList1, InputDirectory1); //read all imagenames
			check2=readStringList(imageList2, InputDirectory2); //read all imagenames
			check=check1 && check2;
			break;
		case 3:
			InputDirectory1=s.InputDirectory1;
			InputDirectory2=s.InputDirectory2;
			InputDirectory3=s.InputDirectory3;
			check1=readStringList(imageList1, InputDirectory1); //read all imagenames
			check2=readStringList(imageList2, InputDirectory2); //read all imagenames
			check3=readStringList(imageList3, InputDirectory3); //read all imagenames
			check=check1 && check2 && check3;
			break;
		}
		
		// Check if imagenaes are succesfully loaded
		if(check){
			cout<<"imageNames loaded"<<endl;
			nimages=imageList1.size(); //number of images
		}
		else{
			cout<<"Images could not be loaded"<<endl;
			exit(1);
		}
		
		// Initialize frames that are skipped on zero
		size_t prev=0;
		//Vector indicating if corner detection was succesfull
		vector<bool> Detection (3, true);
		for(size_t i=0;i<nimages;i++) // Loop over all images
		{

			vector<vector<Point3f>> pixels;
			// Detect corners in image i
			vector<Point3f> ff=readFeaturesImage( imageList1, i, nimages,  cam1,  s); 
			
			// If number of sorted corner points doesn't correspond with pattern size-> skip image
			if(static_cast<int>(ff.size())!=static_cast<int>(s.boardSize.width*s.boardSize.height)){
				Detection[0]=false;
			}
			pixels.push_back(ff); // Add found, sorted corners to list of cornersets
			vector<double> E1, E2, E3, E4, E5, E6; // Errors corresponding to the previous coefficient-sets
			
			if(s.CameraAmount<2){ // Procedure for single camera
				
				optim.changePixels(pixels); // Change FrameOptimizer to cornerlist which requires optimization for this frame
			
				// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
				if(i>4 && Detection[0]){
					E1=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], InitialCoeffs);
					E2=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[i-1-prev]);
					E3=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[i-2-prev]);
					E4=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[i-3-prev]);
					E5=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[i-4-prev]);
					E6=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[i-5-prev]);

				}
			}
			else{ // Procedure for more than 1 camera
				
				// Detect corners in image i
				ff=readFeaturesImage( imageList2, i, nimages,  cam2,  s);
				// If number of sorted corner points doesn't correspond with pattern size-> skip image
				if(static_cast<int>(ff.size())!=static_cast<int>(s.boardSize.width*s.boardSize.height)){
					Detection[1]=false;
				}
				pixels.push_back(ff); // Add found, sorted corners to list of cornersets

				if(s.CameraAmount<3){ //Procedure for 2 cameras

					optim.changePixels(pixels); // Change FrameOptimizer to cornerlist which requires optimization for this frame

					// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
					if(i>4 && Detection[1]){
						vector<double>E11=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], InitialCoeffs);
						E1.insert( E1.end(), E11.begin(), E11.end() );
						vector<double>E22=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[i-1-prev]);
						E2.insert( E2.end(), E22.begin(), E22.end() );
						vector<double>E33=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[i-2-prev]);
						E3.insert( E3.end(), E33.begin(), E33.end() );
						vector<double>E44=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[i-3-prev]);
						E4.insert( E4.end(), E44.begin(), E44.end() );
						vector<double>E55=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[i-4-prev]);
						E5.insert( E5.end(), E55.begin(), E55.end() );
						vector<double>E66=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[i-5-prev]);
						E6.insert( E6.end(), E66.begin(), E66.end() );
					}
				}
				else{ // Procedure for 3 cameras
					
				// Detect corners in image i
				ff=readFeaturesImage( imageList3, i, nimages,  cam3,  s);
				// If number of sorted corner points doesn't correspond with pattern size-> skip image
				if(static_cast<int>(ff.size())!=static_cast<int>(s.boardSize.width*s.boardSize.height)){
					Detection[2]=false;
				}
				pixels.push_back(ff); // Add found, sorted corners to list of cornersets
				
				optim.changePixels(pixels); // Change FrameOptimizer to cornerlist which requires optimization for this frame
			
				// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
				if(i>4 && Detection[2]){
						vector<double>E111=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], InitialCoeffs);
						E1.insert( E1.end(), E111.begin(), E111.end() );
						vector<double>E222=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[i-1-prev]);
						E2.insert( E2.end(), E222.begin(), E222.end() );
						vector<double>E333=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[i-2-prev]);
						E3.insert( E3.end(), E333.begin(), E333.end() );
						vector<double>E444=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[i-3-prev]);
						E4.insert( E4.end(), E444.begin(), E444.end() );
						vector<double>E555=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[i-4-prev]);
						E5.insert( E5.end(), E555.begin(), E555.end() );
						vector<double>E666=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[i-5-prev]);
						E6.insert( E6.end(), E666.begin(), E666.end() );
						}
				}

			}

			//Change frame-optimizer depending on corners found or not
			if(Detection[0]){
				if(Detection[1] &&check2 ){
					if(!Detection[2] &&check3){
						optim.removeCamera(3);
					}
				}
				else{
					if(Detection[2] &&check3){
						optim.swapCameras(2,3);
					}
					else{
						if(check3) optim.removeCamera(3);
						if(check2)optim.removeCamera(2);
					}
				}

			}
			else{
				if(Detection[1] &&check2){
					if(Detection[2] &&check3){
						optim.swapCameras(1,3);
					}
					else{
						if(check2) optim.swapCameras(1,2);
						optim.removeCamera(2);
					}

				}
				else{
					if(Detection[2] &&check3){
						optim.swapCameras(1,3);
						optim.removeCamera(2);
					}
					else{
						cout<<"Non of the cameras found the corner points"<<endl;
						prev++;
						continue;
					}

				}
			}

			// Do optimization with list of sorted cornersets
			real_1d_array coefficients;
			
			if(i<5){ //For images 0-4: use predetermined initial guess
			coefficients = optimizeCoef(InitialCoeffs, optim);
			}
			else{ //For images i>4 -> use initial guess that gives lowest total error
				
				vector<double> S(6); //Vector containing sum of errors
				std::for_each(E1.begin(), E1.end(), [&] (double n) { S[0] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E2.begin(), E2.end(), [&] (double n) { S[1] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E3.begin(), E3.end(), [&] (double n) { S[2] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E4.begin(), E4.end(), [&] (double n) { S[3] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E5.begin(), E5.end(), [&] (double n) { S[4] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E6.begin(), E6.end(), [&] (double n) { S[5] += (nancheck(n) ? 0.0 : n);});
				auto result=min_element(S.begin(),S.end());
				int dist=std::distance(S.begin(), result);
				switch(dist){
					case 0:
						cout<<"Error S[0]: "<<S[0]<<endl;
						coefficients = optimizeCoef(InitialCoeffs, optim);
						break;
					case 1:
						cout<<"Error S[1]: "<<S[1]<<endl;
						coefficients = optimizeCoef(CoefficientsList[i-1-prev], optim);
						break;
					case 2:
						cout<<"Error S[2]: "<<S[2]<<endl;
						coefficients = optimizeCoef(CoefficientsList[i-2-prev], optim);
						break;
					case 3:
						cout<<"Error S[3]: "<<S[2]<<endl;
						coefficients = optimizeCoef(CoefficientsList[i-3-prev], optim);
						break;
					case 4:
						cout<<"Error S[4]: "<<S[4]<<endl;
						coefficients = optimizeCoef(CoefficientsList[i-4-prev], optim);
						break;
					case 5:
						cout<<"Error S[5]: "<<S[5]<<endl;
						coefficients = optimizeCoef(CoefficientsList[i-5-prev], optim);
						break;
					default:
						cout<<"Error sums"<<endl;
						exit(1);
					}
				}
			// Print and store optimized coefficients
			cout<<coefficients.tostring(2).c_str()<<endl; 
			CoefficientsList.push_back(coefficients);
			if(s.saveErrors)
			{ // Save errors corresponding to optimized coefficients
				
				// Compute errors corresponding to optimized coefficients for camera 1
				vector<double> Error;
				if(Detection[0])
					Error=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], coefficients);
				if(check2 && Detection[1])
				{// Compute errors corresponding to optimized coefficients for camera 2
					vector<double> Temp=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
				}
				if(check3 && Detection[2])
				{// Compute errors corresponding to optimized coefficients for camera 3
					vector<double> Temp=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 3 to total error-vector
				}
				double S;
				// Compute mean over all features and save it into errorlist
				std::for_each(Error.begin(), Error.end(), [&] (double n) { S += (nancheck(n) ? 0.0 : n);}); 
				ErrorList.push_back(S/Error.size());
			}

			if(!Detection[0]||!Detection[1]||!Detection[2]){ // Than 1 camera did not detect all corners and was removed
				// Add cameras againfor next image, retaining order of cameras
				switch(s.CameraAmount){ //Choice of cameras and intialize optimalisation with appropriate camera amount

					case 1:
						optim=frameOptimizer(1, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
						break;
					case 2:
						optim=frameOptimizer(2, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
						break;
					case 3:
						optim=frameOptimizer(3, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
						break;
					default:
						cout<<"Invalid amount of camera's"<<endl;
						return 0;
					}
				// Reset boolean detection vector
				vector<bool> Temp (3, true);
				Detection=Temp;
			}

		}

	}
	else{ //Use corner locations from text-files
		// This method assumed that for all text-files of the 3 cameras contain
		// the same amount of detected image points and are taken at corresponding times
		ifstream input1, input2, input3;
		vector<Point3f> ff1, ff2, ff3;

		switch(s.CameraAmount){ // Choice of cameras and open the required input-text files (containing detected corner points)

				case 1:
					input1.open(s.InputFilename1, std::ifstream::in);
					break;
				case 2:
					input1.open(s.InputFilename1, std::ifstream::in);
					input2.open(s.InputFilename2, std::ifstream::in);
					break;

				case 3:
					input1.open(s.InputFilename1, std::ifstream::in);
					input2.open(s.InputFilename2, std::ifstream::in);
					input3.open(s.InputFilename3, std::ifstream::in);
					break;
		}
		
		//Initialisation of temp variables
		int prevTimestep=0;
		int Timestep;
		int Corner;
		float CoordX1, CoordX2, CoordX3;
		float CoordY1, CoordY2, CoordY3;
		string   line1, line2, line3;
		int counter=0; //counter to check number of frames processed
		
		while(std::getline(input1, line1)){ //Keep reading text-files until end of file

			// Input of camera 1
			istringstream   st(line1);
			st>> Timestep >> Corner >>CoordX1 >>CoordY1;
			
			if(Timestep==prevTimestep){ // Than line corresponds to next corner point at same time instance
				ff1.push_back(Point3f(CoordX1, CoordY1, 1)); // Append corner location to list of corner points of that time instance

			if(s.CameraAmount>1){ // Input of more than 1 camera
				
				// Input for camera 2
				std::getline(input2, line2);
				st.str(line2);
				st>> Timestep >> Corner >>CoordX2 >>CoordY2;
				ff2.push_back(Point3f(CoordX2, CoordY2, 1)); // Append corner location to list of corner points of that time instance

			 if(s.CameraAmount>2){ // Input of camera 3

				 std::getline(input3, line3);
				 st.str(line3);
				 st>> Timestep >> Corner >>CoordX3 >>CoordY3;
				 ff3.push_back(Point3f(CoordX3, CoordY3, 1)); // Append corner location to list of corner points of that time instance
			 }
			}
			}
			else{ // Line corresponds to next time instance, for corner point 0

			//append with ff
			vector<vector<Point3f>> pixels;
			pixels.push_back(ff1); // Add corners of camera 1 to vector

			if(s.CameraAmount>1){ 
				pixels.push_back(ff2); // Add corners of camera 2 to vector
			 if(s.CameraAmount>2){
				pixels.push_back(ff3); // Add corners of camera 3 to vector
			 }
			}
			optim.changePixels(pixels);// Change FrameOptimizer to cornerlist which requires optimization for this frame

			vector<double> E1, E2, E3, E4, E5, E6;
			if(s.CameraAmount<2){ // Procedure for single camera
				if(counter>4){ //For images i>4 -> use initial guess that gives lowest total error
					E1=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], InitialCoeffs);
					E2=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[counter-1]);
					E3=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[counter-2]);
					E4=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[counter-3]);
					E5=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[counter-4]);
					E6=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], CoefficientsList[counter-5]);

				}
			}
			else{ // Procedure for more than 1 camera
				if(s.CameraAmount<3){ // Procedure for 2 cameras
					if(counter>4){ //For images i>4 -> use initial guess that gives lowest total error
						vector<double>E11=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], InitialCoeffs);
						E1.insert( E1.end(), E11.begin(), E11.end() );
						vector<double>E22=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[counter-1]);
						E2.insert( E2.end(), E22.begin(), E22.end() );
						vector<double>E33=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[counter-2]);
						E3.insert( E3.end(), E33.begin(), E33.end() );
						vector<double>E44=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[counter-3]);
						E4.insert( E4.end(), E44.begin(), E44.end() );
						vector<double>E55=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[counter-4]);
						E5.insert( E5.end(), E55.begin(), E55.end() );
						vector<double>E66=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], CoefficientsList[counter-5]);
						E6.insert( E6.end(), E66.begin(), E66.end() );
					}
				}
			else{ // Procedure for 3 cameras
				if(counter>4){ //For images i>4 -> use initial guess that gives lowest total error
						vector<double>E111=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], InitialCoeffs);
						E1.insert( E1.end(), E111.begin(), E111.end() );
						vector<double>E222=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[counter-1]);
						E2.insert( E2.end(), E222.begin(), E222.end() );
						vector<double>E333=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[counter-2]);
						E3.insert( E3.end(), E333.begin(), E333.end() );
						vector<double>E444=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[counter-3]);
						E4.insert( E4.end(), E444.begin(), E444.end() );
						vector<double>E555=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[counter-4]);
						E5.insert( E5.end(), E555.begin(), E555.end() );
						vector<double>E666=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], CoefficientsList[counter-5]);
						E6.insert( E6.end(), E666.begin(), E666.end() );
						}

			}

			}

			// Do optimization with list of sorted cornersets
			real_1d_array coefficients;
			
			if(counter<5){  //For images 0-4: use predetermined initial guess
			coefficients = optimizeCoef(InitialCoeffs, optim);
			}
			else{ //For images i>4 -> use initial guess that gives lowest total error
				vector<double> S(6);
				std::for_each(E1.begin(), E1.end(), [&] (double n) { S[0] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E2.begin(), E2.end(), [&] (double n) { S[1] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E3.begin(), E3.end(), [&] (double n) { S[2] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E4.begin(), E4.end(), [&] (double n) { S[3] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E5.begin(), E5.end(), [&] (double n) { S[4] += (nancheck(n) ? 0.0 : n);});
				std::for_each(E6.begin(), E6.end(), [&] (double n) { S[5] += (nancheck(n) ? 0.0 : n);});
				auto result=min_element(S.begin(),S.end());

				int dist=std::distance(S.begin(), result);
				switch(dist){
					case 0:
						cout<<"Error S[0]: "<<S[0]<<endl;
						coefficients = optimizeCoef(InitialCoeffs, optim);
						break;
					case 1:
						cout<<"Error S[1]: "<<S[1]<<endl;
						coefficients = optimizeCoef(CoefficientsList[counter-1], optim);
						break;
					case 2:
						cout<<"Error S[2]: "<<S[2]<<endl;
						coefficients = optimizeCoef(CoefficientsList[counter-2], optim);
						break;
					case 3:
						cout<<"Error S[3]: "<<S[3]<<endl;
						coefficients = optimizeCoef(CoefficientsList[counter-3], optim);
						break;
					case 4:
						cout<<"Error S[4]: "<<S[4]<<endl;
						coefficients = optimizeCoef(CoefficientsList[counter-4], optim);
						break;
					case 5:
						cout<<"Error S[5]: "<<S[5]<<endl;
						coefficients = optimizeCoef(CoefficientsList[counter-5], optim);
						break;
					default:
						cout<<"Error sums"<<endl;
						exit(1);
					}
				}
			
			// Print and store optimized coefficients
			cout<<coefficients.tostring(2).c_str()<<endl;
			CoefficientsList.push_back(coefficients);
			
			if(s.saveErrors)
			{
				// Save errors corresponding to optimized coefficients
				
				// Compute errors corresponding to optimized coefficients for camera 1
				vector<double> Error=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[0], optim.Pixels[0], coefficients);
				if(s.CameraAmount>1){ // Compute errors corresponding to optimized coefficients for camera 2
					vector<double> Temp=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[1], optim.Pixels[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
				}
				if(s.CameraAmount>2){ // Compute errors corresponding to optimized coefficients for camera 3
					vector<double> Temp=compute_error(optim.Lx, optim.Ly, optim.ErrorMetric, optim.Cameras[2], optim.Pixels[2], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 3 to total error-vector
				}
				// Compute mean over all features and save it into errorlist
				double S;
				std::for_each(Error.begin(), Error.end(), [&] (double n) { S += (nancheck(n) ? 0.0 : n);});
				ErrorList.push_back(S/Error.size());
			}
			//Clear vectors and fill first element with line corresponding to next time frame
			vector<Point3f>().swap(ff1);
			ff1.push_back(Point3f(CoordX1, CoordY1, 1));
			vector<Point3f>().swap(ff2);
			ff2.push_back(Point3f(CoordX2, CoordY2, 1));
			vector<Point3f>().swap(ff3);
			ff3.push_back(Point3f(CoordX3, CoordY3, 1));

			counter++; //We have finished a frame-optimization
			prevTimestep=Timestep;
		}
		}
		
	}
	
	//Write results to files
	writeArray(CoefficientsList, s.OutputFileName);
	if(s.saveErrors)
	{writeArrayErrors(ErrorList, s.OutputFileNameErrors);}

	timestamp ( );
    cout << "\n";
    cout << "Surface reconstruction \n";
    cout << "C++ version:\n";
    cout << "By Lukas Engelen:\n";
    return 0;
}

/// Inlcude header files and external source files
#include "main.h"


static inline void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
{///reads node in settings file and stores it in variable
    if(node.empty())
        x = default_value;
    else
        x.read(node);
}
inline string ToString(size_t sz) {

  stringstream ss;
  ss << sz;
  return ss.str();
}
/// Definition template classes

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

inline int nancheck(double x) { return x != x; } //check if nan - double version
inline int nancheck2(float x) { return x != x; } //check if nan - float version

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
		        vector<Point2f> sortedcorners =create_points(image, ResponseThreshold, MinDistance,s.responseRadius, board_sz, cameraMatrix, distCoeffs, s.showCorners);

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
		ifstream input(s.InputInitial);

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

//Optimization of coefficients
vector<double> compute_error(float Lx, float Ly, int ErrorMetric, CameraParams Camera, vector<Point3f> Pixels, vector<Point3f> f, real_1d_array Coeff)
{ /// Compute errors Ef for all feature points of one camera


	//Initialize all temp vectors
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

    transform(Pixels.begin(), Pixels.end(), qworld.begin(), worldCoords<Point3f> (RR, Camera.T, KK)); //compute world coordinates q'
    transform(qworld.begin(), qworld.end(), rays.begin(), Ray<Point3f> (Camera.c)); //compute rays u=cq'

	//Determine surface points according to surface model and hypothesized Coefff
     transform(rays.begin(), rays.end(), water.begin(),  waterSurface_2ndorder<Point3f> (Camera.c, Lx, Ly, Coeff));

	//Compute normals n2 according to surface model and hypothesized Coefff
    transform(water.begin(), water.end(), normals.begin(), Normals_2ndorder<Point3f> (Coeff, Lx,Ly));

    transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>()); // u in alternative way-->gives the same, I checked

    transform (f.begin(), f.end(), water.begin(), v.begin(), SubtractNormalized<Point3f>()); //compute rays v=pf

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
    transform (f.begin(), f.end(), features2.begin(), shift.begin(), Subtract<Point3f>());
    transform (shift.begin(), shift.end(),E.begin(), takeNorm<double>());}

     return E;
}
void combined_error3(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of all cameras combined for set of 3 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[1], Frame.Pixels[1], Frame.fs[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }
	//Errors Ef third camera
	vector<double> E3 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[2], Frame.Pixels[2], Frame.fs[2], x);
	for(size_t k=0; k<E3.size(); k++){
		    fi[k+E1.size()+E2.size()] = (nancheck(E3[k]) ? 0 : E3[k]);
		 }

 }
void combined_error2(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of all cameras combined for set of 2 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[1], Frame.Pixels[1], Frame.fs[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }

 }
void combined_error1(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of single camera 

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr));
	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);


	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization
 }
real_1d_array optimizeCoef(real_1d_array Coeff, frameOptimizer Frame)
{ /// Finds optimal coefficicients for frame
  /// Requires initial guess and frameOptimizer containing all necessary input

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

	   minlmresults(state, Coeff, rep);

	    return Coeff;

	}

//Temporal reconstruction functions
int main ( )
{ /// Main function to reconstruct time-dependent water surface.

	std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

	//Read all settings
	Settings s;
	FileStorage fs("settings.xml", FileStorage::READ); // Read the settings
	if (!fs.isOpened())
		{
			cout << "Could not open the configuration file: "<< endl;
			exit(1);
		 }
	fs["Settings"] >> s;
	fs.release();   // close Settings file
	CameraParams cam1, cam2, cam3;
	int max_threads = omp_get_max_threads();

	if(s.ThreadAmount>max_threads)
	{
		cout<<"Amount of threads unavailable"<<endl;
		cout<<"Changed to maximum amount of threads: "<<max_threads<<endl;
		s.ThreadAmount=max_threads;
	}
	omp_set_dynamic(0);
	omp_set_num_threads(s.ThreadAmount);
	vector<vector<vector<Corner> > > prevPixels(s.ThreadAmount);

	switch(s.CameraAmount){ // Choice of cameras and initialize each camera seperately

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

	vector<vector<real_1d_array> >CoefficientsList(s.ThreadAmount); //Output list of coeffcients
	vector<vector <double>> ErrorList(s.ThreadAmount); //Output list of corresponding errors

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
	InitialCoeffs= "[5.05000,0.000,0.000,0.000,0.000,0.000,0.000,0.000]";
	}

	vector<frameOptimizer> optimizers(s.ThreadAmount);
	for(size_t i=0;i<s.ThreadAmount;i++){
	vector<CameraParams> cams;
	switch(s.CameraAmount){ //Choice of cameras and intialize optimalisation with appropriate camera amount

	case 1:
		cams.push_back(cam1);
		optimizers[i]=frameOptimizer(1, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
		break;
	case 2:
		cams.push_back(cam1);
		cams.push_back(cam2);
		optimizers[i]=frameOptimizer(2, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
		break;
	case 3:
		cams.push_back(cam1);
		cams.push_back(cam2);
		cams.push_back(cam3);
		optimizers[i]=frameOptimizer(3, cams, s.ErrorMetric, s.SurfaceModel, s.Lx, s.Ly);
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
		optimizers[i].setStoppingConditions( i_epsg ,  i_epsf,  i_epsx,  i_maxits);

		cout << "Please enter a value for differentiation step (diffStep): ";
		cin >> i_diffStep;
		optimizers[i].setDiffStep( i_diffStep );
	}
	}

	size_t nimages;
	string InputDirectory1, InputDirectory2, InputDirectory3;
	vector<string> imageList1, imageList2, imageList3;
	// Checks to determine if image-names were read
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

	// Check if image-names are successfully loaded
	if(check){
		cout<<"imageNames/filenames loaded"<<endl;
		nimages=imageList1.size(); //number of images
	}
	else{
		cout<<"Images/files could not be loaded"<<endl;
		exit(1);
	}
	int images_thread=nimages/s.ThreadAmount;
	#pragma omp parallel
	{
		#pragma omp for nowait
		for(size_t j=0;j<s.ThreadAmount;j++){
		for(size_t i=j*images_thread;i<(j+1)*images_thread;i++) // Loop over all images
		//for(size_t i=0;i<nimages;i++) // Loop over all images
		{
			int tid = omp_get_thread_num();

			vector<vector<Corner> > pixels(s.CameraAmount);
			// Detect corners in image i
			switch(i-j*images_thread>0.1){
			
			case true:
				switch(s.CameraAmount){

							case 1:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1, s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s, prevPixels[tid][0]);
								break;
							case 2:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1,  s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s, prevPixels[tid][0]);
								pixels[1]=readFeaturesImage( imageList2[nimages-i-1], cam2,  s.OutputDirectory2+"/CAM2/"+ToString(i)+".txt", s, prevPixels[tid][1]);
								break;
							case 3:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1,  s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s, prevPixels[tid][0]);
								pixels[1]=readFeaturesImage( imageList1[nimages-i-1], cam2,  s.OutputDirectory2+"/CAM2/"+ToString(i)+".txt", s, prevPixels[tid][1]);
								pixels[2]=readFeaturesImage( imageList3[nimages-i-1], cam3,  s.OutputDirectory3+"/CAM3/"+ToString(i)+".txt", s, prevPixels[tid][2]);
								break;
							default:
								cout<<"Error camera amount"<<endl;
								exit(1);
						}
				break;
			case false:
				switch(s.CameraAmount){

							case 1:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1, s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s);
								break;
							case 2:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1, s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s);
								pixels[1]=readFeaturesFirstImage( imageList2[nimages-i-1], cam2, s.OutputDirectory2+"/CAM2/"+ToString(i)+".txt", s);
								break;
							case 3:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1, s.OutputDirectory1+"/CAM1/"+ToString(i)+".txt", s);
								pixels[1]=readFeaturesFirstImage( imageList2[nimages-i-1], cam2, s.OutputDirectory2+"/CAM2/"+ToString(i)+".txt", s);
								pixels[2]=readFeaturesFirstImage( imageList3[nimages-i-1], cam3, s.OutputDirectory3+"/CAM3/"+ToString(i)+".txt", s);
								break;
							default:
								cout<<"Error camera amount"<<endl;
								exit(1);
						}
				break;
			}
			prevPixels[tid]=pixels;
			vector<double> E1, E2, E3, E4, E5, E6; // Errors corresponding to the previous coefficient-sets
			optimizers[tid].changePixels(pixels); // Change frameOptimizer to cornerlist which requires optimization for this frame
			
			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				E1=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], InitialCoeffs);
				E2=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-1]);
				E3=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-2]);
				E4=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-3]);
				E5=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-4]);
				E6=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-5]);

			}

			if(s.CameraAmount>1){ // Procedure for more than 1 camera
				
			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				vector<double>E11=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], InitialCoeffs);
				E1.insert( E1.end(), E11.begin(), E11.end() );
				vector<double>E22=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-1]);
				E2.insert( E2.end(), E22.begin(), E22.end() );
				vector<double>E33=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-2]);
				E3.insert( E3.end(), E33.begin(), E33.end() );
				vector<double>E44=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-3]);
				E4.insert( E4.end(), E44.begin(), E44.end() );
				vector<double>E55=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-4]);
				E5.insert( E5.end(), E55.begin(), E55.end() );
				vector<double>E66=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-5]);
				E6.insert( E6.end(), E66.begin(), E66.end() );
			}
			}
			if(s.CameraAmount>2){ // Procedure for 3 cameras

			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				vector<double>E111=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], InitialCoeffs);
				E1.insert( E1.end(), E111.begin(), E111.end() );
				vector<double>E222=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-1]);
				E2.insert( E2.end(), E222.begin(), E222.end() );
				vector<double>E333=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-2]);
				E3.insert( E3.end(), E333.begin(), E333.end() );
				vector<double>E444=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-3]);
				E4.insert( E4.end(), E444.begin(), E444.end() );
				vector<double>E555=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-4]);
				E5.insert( E5.end(), E555.begin(), E555.end() );
				vector<double>E666=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-5]);
				E6.insert( E6.end(), E666.begin(), E666.end() );
				}
			}

			// Do optimization with list of sorted corner-sets
			real_1d_array coefficients;
			
			if(i-j*images_thread<5){ //For images 0-4: use predetermined initial guess
			coefficients = optimizeCoef(InitialCoeffs, optimizers[tid]);
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
						//cout<<"Error S[0]: "<<S[0]<<endl;
						coefficients = optimizeCoef(InitialCoeffs, optimizers[tid]);
						break;
					case 1:
						//cout<<"Error S[1]: "<<S[1]<<endl;
						coefficients = optimizeCoef(CoefficientsList[tid][i-j*images_thread-1], optimizers[tid]);
						break;
					case 2:
						//cout<<"Error S[2]: "<<S[2]<<endl;
						coefficients = optimizeCoef(CoefficientsList[tid][i-j*images_thread-2], optimizers[tid]);
						break;
					case 3:
						//cout<<"Error S[3]: "<<S[2]<<endl;
						coefficients = optimizeCoef(CoefficientsList[tid][i-j*images_thread-3], optimizers[tid]);
						break;
					case 4:
						//cout<<"Error S[4]: "<<S[4]<<endl;
						coefficients = optimizeCoef(CoefficientsList[tid][i-j*images_thread-4], optimizers[tid]);
						break;
					case 5:
						//cout<<"Error S[5]: "<<S[5]<<endl;
						coefficients = optimizeCoef(CoefficientsList[tid][i-j*images_thread-5], optimizers[tid]);
						break;
					default:
						cout<<"Error sums"<<endl;
						exit(1);
					}
				}
			// Store optimized coefficients
			CoefficientsList[tid].push_back(coefficients);

			if(s.saveErrors)
			{ // Save errors corresponding to optimized coefficients
				
				// Compute errors corresponding to optimized coefficients for camera 1
				vector<double> Error;
				switch (s.CameraAmount){
				case 1:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					break;}
				case 2:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					vector<double> Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					break;}
				case 3:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					vector<double> Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					break;}
				default:
					cout<<"Error camera amount"<<endl;
					exit(1);
				}
			double S;
				// Compute mean over all features and save it into errorlist
				std::for_each(Error.begin(), Error.end(), [&] (double n) { S += (nancheck(n) ? 0.0 : n);}); 
				ErrorList[tid].push_back(S/Error.size());
			}
		}
			}}


	//Write results to files
	writeArray(CoefficientsList, s.OutputFileName);
	if(s.saveErrors)
	{writeArrayErrors(ErrorList, s.OutputFileNameErrors);}
	std::chrono::steady_clock::time_point end1= std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() <<std::endl;

	timestamp ( );
    cout << "\n";
    cout << "Surface reconstruction \n";
    cout << "C++ version:\n";
    cout << "By Lukas Engelen:\n";
    return 0;
}

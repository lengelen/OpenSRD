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
/// Inlcude header files and external source files
#include "main.h"
#include <iterator>     // std::istream_iterator

static inline void read(const FileNode& node, Settings& x, const Settings& default_value = Settings())
{// Reads node in settings file and stores it in variable
    if(node.empty())
        x = default_value;
    else
        x.read(node);
}
inline string ToString(size_t sz) {
// Converts size_t object to string object
  stringstream ss;
  ss << sz;
  return ss.str();
}

/// Definition template classes

template <class Type> class worldCoords
{ // Determines world coordinates (3D Float) of image point in vector 3-Float format (homog)
private:
    Mat RR, T, KK;   // inverse camera matrices
public:
    // Constructor
    worldCoords (Mat i_RR, Mat i_T, Mat i_KK) {
        RR=i_RR;
        T=i_T;
        KK=i_KK;
    }

    // The function call
    Type operator ( ) ( Point3f& elem ) const
    {
    	// Change elem to matrix object to perform calculations
    	Mat q=Mat(elem).reshape(1,3);
        q.convertTo(q, CV_64F);

    	// Compute result
        Mat result= (RR * (KK*q-T));

        // Convert result to Point3f
        Mat resultf;
        result.convertTo(resultf, CV_32F);
        return Point3f(resultf.at<float>(0,0),resultf.at<float>(1,0), resultf.at<float>(2,0));

     }
};
template <class Type> class pixelCoordinates
{ // Determines pixel coordinates of 3D point (float), output 3Float (homog)
private:
    Mat R, T, K;   // camera matrices
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
    	// Change elem to matrix object to perform calculations
    	Mat q=Mat(elem).reshape(1,3);
    	q.convertTo(q, CV_64F);

    	// Compute result
        Mat result= (K * (R*q+T));

        // Convert result to Point3f
        Mat resultf;
        result.convertTo(resultf, CV_32F);
        return Point3f(resultf.at<float>(0,0)/ resultf.at<float>(2,0),resultf.at<float>(1,0)/ resultf.at<float>(2,0),1.0);
        }
};
template <class Type> class Ray
{ // Computes ray direction (3D Float) between c and image point
  // Resulting 3f vector is normalized
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
            	Point3f res=Point3f(elem.x-c.x,elem.y-c.y, elem.z-c.z);
            	return res/norm(res);
            }
};
template <class Type> class Normals_constant
{  // Computes normalized n2 of set of 1 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
        public:
            // Constructor
            Normals_constant (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {

            	return Point3f(0, 0, 1);
            }
        };
template <class Type> class Normals_linear
{  // Computes normalized n2 of set of 3 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
        public:
            // Constructor
            Normals_linear (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {
               Point3f n= Point3f(-coef[1]/Lx,-coef[2]/Ly, 1);
               return n/norm(n);

            }
        };
template <class Type> class Normals_1storder
{  // Computes normalized n2 of set of 8 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
        public:
            // Constructor
            Normals_1storder (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {

               Point3f n= Point3f(PI_F/Lx*coef[3]*sin(PI_F*elem.x/Lx)-coef[1]/Lx, PI_F/Ly*coef[4]*sin(PI_F*elem.y/Ly)-coef[2]/Ly, 1);
               return n/norm(n);

            }
        };
template <class Type> class Normals_2ndorder
{  // Computes normalized n2 of set of 8 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
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

               Point3f n= Point3f(2*PI_F/Lx*coef[5]*sin(2*PI_F*elem.x/Lx)+PI_F/Lx*coef[7]*sin(PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[3]*sin(PI_F*elem.x/Lx)-coef[1]/Lx, 2*PI_F/Ly*coef[6]*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[7]*cos(PI_F*elem.x/Lx)*sin(PI_F*elem.y/Ly)+PI_F/Ly*coef[4]*sin(PI_F*elem.y/Ly)-coef[2]/Ly, 1);
               return n/norm(n);

            }
        };
template <class Type> class Normals_3rdorder
{  // Computes normalized n2 of set of 12 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
        public:
            // Constructor
            Normals_3rdorder (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Lx;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {

               Point3f n= Point3f(3*PI_F/Lx*coef[8]*sin(3*PI_F*elem.x/Lx)+2*PI_F/Lx*coef[5]*sin(2*PI_F*elem.x/Lx)+2*PI_F/Lx*coef[10]*sin(2*PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[7]*sin(PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[3]*sin(PI_F*elem.x/Lx)-coef[1]/Lx, 3*PI_F/Ly*coef[9]*sin(3*PI_F*elem.y/Ly)+2*PI_F/Ly*coef[6]*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[7]*cos(PI_F*elem.x/Lx)*sin(PI_F*elem.y/Ly)+2*coef[11]*cos(PI_F*elem.x/Lx)*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[4]*sin(PI_F*elem.y/Ly)-coef[2]/Ly, 1);
               return n/norm(n);

            }
        };
template <class Type> class Angle
{ // Computes theta_i for given theta_delta and rw
    private:
        double rw; // Refractive index of water

    public:
            // Constructor
            Angle (double i_rw) {
                rw=i_rw;
            }

            // The function call
            Type operator ( ) ( double& theta ) const
            {
                return atan((rw*sin(theta))/(rw*cos(theta)-1));
            }
        };
template <class Type> class Angle2
{ // Computes theta_delta for given theta_i and rw
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

               return (theta - asin(sin(theta)/rw));
            }
        };
template <class Type> class waterSurface_constant
{ // Finds surface point for set of 1 coefficient and ray u=cq' (all 3D Float)

	private:
	real_1d_array coef; Point3f c;  double Lx; double Ly; // Characterization of surface model

	public:
		// Constructor
	waterSurface_constant (Point3f i_c, double i_Lx, double i_Ly, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
		}

		// The function call
		Type operator ( ) ( Point3f& ray ) const
		{
			// Parameters/variables used by Hybrid to solve non-linear problem
			double *fvec;
			int iflag=1;
			int lwa;
			int n = 4;
			double tol = 0.000001;
			double *wa;
			double *x;
			double *parameters;

			lwa = ( n * ( 3 * n + 13 ) ) / 2;

			fvec = new double[n];
			wa = new double[lwa];
			x = new double[n];
			parameters = new double [9];

			// Initial guess
			x[0] = 1500;
			x[1] = 0;
			x[2] = 0;
			x[3] = 0;

			// Additional input information passed by pointer
			parameters[0] = c.x;
			parameters[1] = c.y;
			parameters[2] = c.z;
			parameters[3] = ray.x;
			parameters[4] = ray.y;
			parameters[5] = ray.z;
			parameters[6] = Lx;
			parameters[7] = Ly;
			parameters[8] = coef[0];

			// Initialization
			findIntersection_constant( n, x, fvec, &iflag, parameters );

			// Optimization
			int info = hybrd1 ( findIntersection_constant, n, x, fvec, parameters, tol, wa, lwa );
			Point3f result =Point3f(x[1],x[2],x[3]);
			delete[] fvec;
			delete[] wa;
			delete[] x;
			delete[] parameters;
			return result;
			}
	};
template <class Type> class waterSurface_linear
{ // Finds surface point for set of 3 coefficients and ray u=cq' (all 3D Float)

	private:
	real_1d_array coef; Point3f c;  double Lx; double Ly; // Characterization of surface model

	public:
		// Constructor
	waterSurface_linear (Point3f i_c, double i_Lx, double i_Ly, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
		}

		// The function call
		Type operator ( ) ( Point3f& ray ) const
		{
			// Parameters/variables used by Hybrid to solve non-linear problem
			double *fvec;
			int iflag=1;
			int lwa;
			int n = 4;
			double tol = 0.000001;
			double *wa;
			double *x;
			double *parameters;

			lwa = ( n * ( 3 * n + 13 ) ) / 2;

			fvec = new double[n];
			wa = new double[lwa];
			x = new double[n];
			parameters = new double [11];

			// Initial guess
			x[0] = 1500;
			x[1] = 0;
			x[2] = 0;
			x[3] = 0;

			// Additional input information passed by pointer
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

			// Initialization
			findIntersection_linear( n, x, fvec, &iflag, parameters );

			// Optimization
			int info = hybrd1 ( findIntersection_linear, n, x, fvec, parameters, tol, wa, lwa );
			Point3f result =Point3f(x[1],x[2],x[3]);
			delete[] fvec;
			delete[] wa;
			delete[] x;
			delete[] parameters;
			return result;
			}
	};

template <class Type> class waterSurface_1storder
{ // Finds surface point for set of 8 coefficients and ray u=cq' (all 3D Float)

	private:
	real_1d_array coef; Point3f c;  double Lx; double Ly; // Characterization of surface model

	public:
		// Constructor
	waterSurface_1storder (Point3f i_c, double i_Lx, double i_Ly, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
		}

		// The function call
		Type operator ( ) ( Point3f& ray ) const
		{
			// Parameters/variables used by Hybrid to solve non-linear problem
			double *fvec;
			int iflag=1;
			int lwa;
			int n = 4;
			double tol = 0.000001;
			double *wa;
			double *x;
			double *parameters;

			lwa = ( n * ( 3 * n + 13 ) ) / 2;

			fvec = new double[n];
			wa = new double[lwa];
			x = new double[n];
			parameters = new double [13];

			// Initial guess
			x[0] = 1500;
			x[1] = 0;
			x[2] = 0;
			x[3] = 0;

			// Additional input information passed by pointer
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

			// Initialization
			findIntersection_1storder( n, x, fvec, &iflag, parameters );

			// Optimization
			int info = hybrd1 ( findIntersection_1storder, n, x, fvec, parameters, tol, wa, lwa );
			Point3f result =Point3f(x[1],x[2],x[3]);
			delete[] fvec;
			delete[] wa;
			delete[] x;
			delete[] parameters;
			return result;
			}
        };
template <class Type> class waterSurface_2ndorder
{ // Finds surface point for set of 8 coefficients and ray u=cq' (all 3D Float)

        private:
		real_1d_array coef; Point3f c;  double Lx; double Ly; // Characterization of surface model

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
            	// Parameters/variables used by Hybrid to solve non-linear problem
                double *fvec;
                int iflag=1;
                int lwa;
                int n = 4;
                double tol = 0.000001;
                double *wa;
                double *x;
                double *parameters;

                lwa = ( n * ( 3 * n + 13 ) ) / 2;

                fvec = new double[n];
                wa = new double[lwa];
                x = new double[n];
                parameters = new double [16];

                // Initial guess
                x[0] = 1500;
                x[1] = 0;
                x[2] = 0;
                x[3] = 0;

                // Additional input information passed by pointer
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

				// Initialization
                findIntersection_2ndorder( n, x, fvec, &iflag, parameters );

                // Optimization
                int info = hybrd1 ( findIntersection_2ndorder, n, x, fvec, parameters, tol, wa, lwa );
                Point3f result =Point3f(x[1],x[2],x[3]);
                delete[] fvec;
                delete[] wa;
                delete[] x;
                delete[] parameters;
                return result;
                }
        };
template <class Type> class waterSurface_3rdorder
{ // Finds surface point for set of 12 coefficients and ray u=cq' (all 3D Float)

	private:
	real_1d_array coef; Point3f c;  double Lx; double Ly; // Characterization of surface model

	public:
		// Constructor
	waterSurface_3rdorder (Point3f i_c, double i_Lx, double i_Ly, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
		}

		// The function call
		Type operator ( ) ( Point3f& ray ) const
		{
			// Parameters/variables used by Hybrid to solve non-linear problem
			double *fvec;
			int iflag=1;
			int lwa;
			int n = 4;
			double tol = 0.000001;
			double *wa;
			double *x;
			double *parameters;

			lwa = ( n * ( 3 * n + 13 ) ) / 2;

			fvec = new double[n];
			wa = new double[lwa];
			x = new double[n];
			parameters = new double [20];

			// Initial guess
			x[0] = 1500;
			x[1] = 0;
			x[2] = 0;
			x[3] = 0;

			// Additional input information passed by pointer
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
			parameters[16] = coef[8];
			parameters[17] = coef[9];
			parameters[18] = coef[10];
			parameters[19] = coef[11];

			// Initialization
			findIntersection_3rdorder( n, x, fvec, &iflag, parameters );

			// Optimization
			int info = hybrd1 ( findIntersection_3rdorder, n, x, fvec, parameters, tol, wa, lwa );
			Point3f result =Point3f(x[1],x[2],x[3]);
			delete[] fvec;
			delete[] wa;
			delete[] x;
			delete[] parameters;
			return result;
			}
	};
template <class T> struct prod
{ // Computes angle between two vectors by dot product (used to compute theta_d)
    T operator() (Point3f& uu, Point3f& vv) const {

    	return acos(uu.ddot(vv));
    	}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct prod2
{// Computes angle between two vectors by dot product but reverses one vector in direction (used for u and n (thetai))
    T operator() (Point3f& u, Point3f& n) const {

    	Point3f uu=Point3f(-u.x,-u.y,-u.z);
    	return acos(uu.ddot(n));
    	}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct error_col
{ // Computes normal collinearity metric based on two (normalized) normal directions
    T operator() (Point3f& n1, Point3f& n2) const
    {
    	return 1000*acos(n1.ddot(n2));
    }
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Waterray
{ // Computes intersection of ray v2 with plane F
  // Based on distance a along ray-vector v, return Point3f
    T operator() (float& a, Point3f& v) const
    {
    	return Point3f(a*v.x,a*v.y,a*v.z);
 	    }
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct axis
{ // Computes rotation axis by cross product of û and ^v (order of vectors is important!)
	T operator() (Point3f& uu, Point3f& vv) const {

		return vv.cross(uu);

		}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Homogenous
{ // Transforms 2D image point in homogeneous coordinates
	T operator() (Point2f& q) const {
		return Point3f(q.x,q.y,1);}
	typedef T first_argument_type;
	typedef T second_argument_type;
	typedef T result_type;
};
template <class T> struct takeNorm
{ // Computes norm 3D Float image point  (OpenCV gave troubles)
        T operator() (Point3f& p) const {
        	double l =sqrt( p.x*p.x + p.y*p.y + p.z*p.z );
        	return l;}
        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
    };
template <class T> struct Rotationmatrix
{ // Computes rotationmatrix around axis -axis and over angle theta
    T operator() (double& theta, Point3f& axis) const {

    	Mat_<double> r=Mat(-axis).reshape(1,3);
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
{ // Rotates -u with rotationmatrix
      T operator() (Mat& R, Point3f& u) const {

          Mat uu=Mat(-u).reshape(1,3);
          Mat n2=R*uu;
          float nn =norm(n2, NORM_L2);

          return Point3f(n2.at<float>(0,0)/nn, n2.at<float>(1,0)/nn, n2.at<float>(2,0)/nn);
          }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Rotation
{// Rotates u with rotationmatrix
      T operator() (Mat& R, Point3f& u) const {

          Mat uu=Mat(u).reshape(1,3);
          Mat n2=R*uu;
          float nn =norm(n2, NORM_L2);

          return Point3f(n2.at<float>(0,0)/nn, n2.at<float>(1,0)/nn, n2.at<float>(2,0)/nn);
          }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Distance
{ // Computes distance along v from surface point p towards plane z=0
	    T operator() (Point3f& p, Point3f& v) const {
          return (-p.z/v.z);
          }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct SubtractNormalized
{ // Subtract two 3D Float vectors and gives normalized result
    T operator() (Point3f& v1,Point3f& v2) const {
        Point3f res= Point3f(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
        return res/norm(res);

    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct Subtract
{  // Subtract two 3D Float vectors, gives unnormalized result
    T operator() (Point3f& v1,Point3f& v2) const {
        return Point3f(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);

    }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };
template <class T> struct AbsoluteAngle
{ // Takes 3D point with 3rd dim angle and computes abs of 3rd dim
  // Returns same 3D point with abs value of angle
    T operator() (Point3f p) const {
        return Point3f(p.x,p.y, abs(p.z));
        }

        typedef T first_argument_type;
        typedef T second_argument_type;
        typedef T result_type;
        };

inline int nancheck(double x) { return x != x; } //check if nan - double version
inline int nancheck2(float x) { return x != x; } //check if nan - float version


CameraParams Initialize(Settings s, int camera)
{ // Initialization of camera parameters in a CameraParams-object for input of settings and camera number

	// Initialization and reading of all settings
	Mat cameraMatrix, distCoeffs;
	string in="INPUT/";
	string CalibrationFile,  InputReference, OutputFile, CameraPoseFile, FeaturesFile;
	Mat R, origin, RR, rvec, tvec;
	int Position;
	Point3f c;
	switch(camera){ // Open calibration file and reference image depending on camera
		case 1:
			 CalibrationFile= in+s.CalibrationFile1;
			 InputReference= in+s.InputReference1;
			 CameraPoseFile= in+s.InputInitial1;
			 OutputFile=s.OutputCameraPose1;
			 FeaturesFile=in+s.InputFeatures1;
			 break;
		case 2:
			 CalibrationFile= in+s.CalibrationFile2;
			 InputReference= in+s.InputReference2;
			 CameraPoseFile= in+s.InputInitial2;
			 OutputFile=s.OutputCameraPose2;
			 FeaturesFile=in+s.InputFeatures2;

			 break;
		case 3:
			 CalibrationFile=in+s.CalibrationFile3;
			 InputReference= in+s.InputReference3;
			 CameraPoseFile= in+s.InputInitial3;
			 OutputFile=s.OutputCameraPose3;
			 FeaturesFile=in+s.InputFeatures3;
			 break;
	}

	Size board_sz = s.RefPatternSize;
	float ResponseThreshold=s.ResponseThreshold;
	float MinDistance=s.MinDistance;

	//Read in calibration information:distortion coefficients and cameramatrix
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


	 if(!s.TypeCameraPose){ // Camera pose estimation aleady done, read in external parameters from file
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
	 else{ // Read in reference image for camera position and orientation
			Mat image=imread(InputReference, IMREAD_GRAYSCALE);
			if( image.data )
			cout <<"extra reference image loaded:"<<InputReference <<endl;
			else exit(1);

			// Detect corner points in image but don't undistort them
			vector<Point2f> sortedcorners =create_points(image, ResponseThreshold, MinDistance, s.ResponseRadius, board_sz, cameraMatrix, distCoeffs, s.ShowCorners);

			// Read theoretical 3D world coordinates
			ifstream input(CameraPoseFile);
			while(input.fail())
					{
						cout<< "File for initial location of f incorrect" << endl;
						exit(1); ;
					}
			vector<Point3f> objectPoints;
			Point3f tmp;
			while (input >> tmp.x && input >> tmp.y && input >> tmp.z)
			{

				objectPoints.push_back(tmp);
			};

			if(objectPoints.size()!=board_sz.width*board_sz.height){
				cout<<"File containing camera pose estimation vertices does not match with board_sz"<<endl;
				exit(1);
			}
			// Perform camera pose estimation
			bool succes= solvePnP(objectPoints,  sortedcorners,  cameraMatrix, distCoeffs,  rvec, tvec, false); //camera pose estimation
		    if (!succes){
			   cout<< "Initialization not succeeded" << endl;
			   exit(1);}

		    Rodrigues(rvec, R); // Computation of rotationmatrix from rotationvector with Rodrigues formula
		    RR = R.inv();  // rotation of inverse
		    origin = -RR * tvec; // translation of inverse
		    c=Point3f(origin.at<double>(0,0), origin.at<double>(1,0), origin.at<double>(2,0)); // Camera center

		 }

	 	// Read fixed locations of feature points f on feature plane F
		ifstream input(FeaturesFile);
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
		if(s.SaveCameraPose){ // Save result camera pose estimation if needed

		saveCameraParams(OutputFile,in, R, tvec);

		}
		CameraParams cam=CameraParams(R, tvec, cameraMatrix, c, distCoeffs, f);
        return cam;
}

//Find intersection of rays with with surfaces
void findIntersection_constant( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 1 coefficient
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00;
    return;
}
void findIntersection_linear( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 3 coefficients
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double B = params[9];
    double C = params[10];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 +B*x[1]/Lx+C*x[2]/Ly;
    return;
}

void findIntersection_1storder( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 5 coefficients
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double B = params[9];
    double C = params[10];
    double A10 = params[11];
    double A01 = params[12];


    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*x[1]/Lx+C*x[2]/Ly;
    return;
}
void findIntersection_2ndorder( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 8 coefficients
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double B = params[9];
    double C = params[10];
    double A10 = params[11];
    double A01 = params[12];
	double A20 = params[13];
	double A02 = params[14];
	double A11 = params[15];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*x[1]/Lx+C*x[2]/Ly + A20 * cos(2*PI*x[1]/Lx)+ A02 * cos(2*PI*x[2]/Ly) + A11 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly);
    return;
}
void findIntersection_3rdorder( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 8 coefficients
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double A00 = params[8];
    double B = params[9];
    double C = params[10];
    double A10 = params[11];
    double A01 = params[12];
	double A20 = params[13];
	double A02 = params[14];
	double A11 = params[15];
	double A30 = params[16];
	double A03 = params[17];
	double A21 = params[18];
	double A12 = params[19];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*x[1]/Lx+C*x[2]/Ly + A20 * cos(2*PI*x[1]/Lx)+ A02 * cos(2*PI*x[2]/Ly) + A11 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly)+ A30 * cos(3*PI*x[1]/Lx)+ A03 * cos(3*PI*x[2]/Ly) + A21 * cos(PI*x[1]/Lx*2) * cos(PI*x[2]/Ly)+ A12 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly*2);
    return;
}
//Optimization of coefficients
vector<double> compute_error(float Lx, float Ly, int ErrorMetric, int ParameterAmount, CameraParams Camera, vector<Point3f> Pixels, vector<Point3f> f, real_1d_array Coeff)
{ // Compute errors Ef for all feature points of one camera


	//Initialize all temp vectors
	vector<Point3f> qworld(Pixels.size());
    vector<Point3f> water(Pixels.size());
    vector<Point3f> normals1(Pixels.size());
    vector<Point3f> u(Pixels.size());
    vector<Point3f> v(Pixels.size());
    vector<Point3f> v2(Pixels.size());
    vector<Point3f> X(Pixels.size());
    vector<Point3f> X2(Pixels.size());
    vector<Point3f> normals2(Pixels.size());

    vector<double> thetad(Pixels.size());
    vector<double> E(Pixels.size());
    vector<double> thetai(Pixels.size());
    vector<double> thetai2(Pixels.size());
    vector<double> thetad2(Pixels.size());
    vector<double> thetar(Pixels.size());
    vector<float> a(Pixels.size());

    vector<Mat> Rotations(Pixels.size());
    vector<Mat> Rotations2(Pixels.size());
    vector<Point3f> features2(Pixels.size());
    vector<Point3f> uu(Pixels.size());
    vector<Point3f> vv(Pixels.size());
    vector<Point3f> shift(Pixels.size());

    // Compute inverse of camera matrices
    Mat RR=Camera.R.inv();
    Mat KK=Camera.K.inv();

    // Compute world coordinates q'
    transform(Pixels.begin(), Pixels.end(), qworld.begin(), worldCoords<Point3f> (RR, Camera.T, KK));
    // Compute rays u=cq'
    transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (Camera.c));
	switch(ParameterAmount){
		case 1:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (Camera.c, Lx, Ly, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_constant<Point3f> (Coeff, Lx,Ly));
		break;
		case 3:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (Camera.c, Lx, Ly, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_linear<Point3f> (Coeff, Lx,Ly));
		break;
		case 5:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_1storder<Point3f> (Camera.c, Lx, Ly, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_1storder<Point3f> (Coeff, Lx,Ly));
		break;
		case 8:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_2ndorder<Point3f> (Camera.c, Lx, Ly, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_2ndorder<Point3f> (Coeff, Lx,Ly));
		break;
		case 12:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_3rdorder<Point3f> (Camera.c, Lx, Ly, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_3rdorder<Point3f> (Coeff, Lx,Ly));
		break;
	}
    // Compute rays v=pf
    transform (f.begin(), f.end(), water.begin(), v.begin(), SubtractNormalized<Point3f>());
    // Compute thetad
    transform (u.begin(), u.end(), v.begin(), thetad.begin(),prod<double>());
    // Compute rotation axis
    transform (u.begin(), u.end(), v.begin(), X.begin(),axis<Point3f> ());
    // Compute thetai
    transform(thetad.begin(), thetad.end(), thetai.begin(), Angle<double>(rw));
    // Compute rotation matrix
    transform (thetai.begin(), thetai.end(), X.begin(), Rotations.begin(), Rotationmatrix<Mat>());
    // Rotate u with rotationmatrix
    transform (Rotations.begin(), Rotations.end(), u.begin(), normals1.begin(), Rotationmin<Point3f>());
    if(ErrorMetric==1) // Compute normal collinearity metric
    transform (normals2.begin(), normals2.end(),normals1.begin(), E.begin(), error_col<double>());
    else if(ErrorMetric==2){ // Compute disparity difference metric
    // Compute angle between normals n2 and rays u
    transform (u.begin(), u.end(), normals2.begin(), thetai2.begin(),prod2<double>());
    // Compute second angle theta_delta2
    transform(thetai2.begin(), thetai2.end(), thetad2.begin(), Angle2<double>(rw));
    // Compute rotation axis
    transform (normals2.begin(), normals2.end(), u.begin(), X2.begin(), axis<Point3f> ());
    // Compute rotation matrix
    transform (thetad2.begin(), thetad2.end(), X2.begin(), Rotations2.begin(), Rotationmatrix<Mat>());
    // Rotate vector u to get v2
    transform (Rotations2.begin(), Rotations2.end(), u.begin(), v2.begin(), Rotation<Point3f>());
    // Calculate distance along v2 from surface point to z=0
    transform (water.begin(), water.end(), v2.begin(), a.begin(), Distance<float>());
    // Compute vector from p, along v with length a
    transform (a.begin(), a.end(), v2.begin(), vv.begin(), Waterray<Point3f>());
    // Compute intersections with plane z=0
    transform (water.begin(), water.end(), vv.begin(), features2.begin(), std::plus<Point3f>());
    // Compute disparity shift
    transform (f.begin(), f.end(), features2.begin(), shift.begin(), Subtract<Point3f>());
    // Compute norm of disparity shift
    transform (shift.begin(), shift.end(),E.begin(), takeNorm<double>());
   }
     return E;
}
void combined_error3(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ // Compute errors of all cameras combined for set of 3 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr)); // Get frameoptimizer from pointer

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[1], Frame.Pixels[1], Frame.fs[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }
	//Errors Ef third camera
	vector<double> E3 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params,Frame.Cameras[2], Frame.Pixels[2], Frame.fs[2], x);
	for(size_t k=0; k<E3.size(); k++){
		    fi[k+E1.size()+E2.size()] = (nancheck(E3[k]) ? 0 : E3[k]);
		 }

 }
void combined_error2(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of all cameras combined for set of 2 cameras

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr)); // Get frameoptimizer from pointer

	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization

	//Errors Ef second camera
	vector<double> E2 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[1], Frame.Pixels[1], Frame.fs[1], x);
	for(size_t k=0; k<E2.size(); k++){
		    	fi[k+E1.size()] = (nancheck(E2[k]) ? 0 : E2[k]);

		 }

 }
void combined_error1(const real_1d_array &x, real_1d_array &fi, void *ptr)
{ /// Compute errors of single camera 

	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr)); // Get frameoptimizer from pointer
	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		fi[k] = (nancheck(E1[k]) ? 0.0 : E1[k]) ; //if nan then no contribution to error-> not included in optimization
 }
void combined_errorBF(const real_1d_array &x, double &func, void *ptr)
{ /// Compute errors of single camera
	func=0;
	frameOptimizer Frame=*(static_cast<frameOptimizer*>(ptr)); // Get frameoptimizer from pointer
	//Errors Ef first camera
	vector<double> E1 =compute_error(Frame.Lx, Frame.Ly, Frame.ErrorMetric, Frame.Params, Frame.Cameras[0], Frame.Pixels[0], Frame.fs[0], x);

	for(size_t k=0; k<E1.size(); k++)
		func = (nancheck(E1[k]) ? 0.0 : func+E1[k]) ; //if nan then no contribution to error-> not included in optimization
 }

real_1d_array optimizeCoef(real_1d_array Coeff, frameOptimizer Frame)
{ // Finds optimal coefficients for frame
  // Requires initial guess and frameOptimizer containing all necessary input


		// Optimization configuration parameters
		minlmstate state;
	    minlmreport rep;
	    void *ptr;
	    ptr=&Frame;
	    real_1d_array scale;
		// Initialize optimization
		minlmcreatev(Frame.NumberOfCameras*Frame.Pixels[0].size(), Coeff, Frame.DiffStep, state);
		minlmsetcond(state, Frame.Epsg, Frame.Epsf, Frame.Epsx, Frame.MaxIts);

		switch(Coeff.length()){

			case 1:
				scale="[1]";
				break;
			case 3:
				scale="[1,1,1]";
				break;

			case 5:
				scale="[1,1,1,1,1]";
				break;

			case 8:
				scale="[1,1,1,1,1,1,1,1]";
				break;

			case 12:
				scale="[1,1,1,1,1,1,1,1,1,1,1,1]";
		}

		minlmsetscale(state, scale);



		switch(Frame.NumberOfCameras){ //Determine camera-amount and use appropriate function to compute errors
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
	   minlmresults(state, Coeff, rep);


	   if(Frame.Params != (int) Coeff.length()){
		   real_1d_array Coeff2;
		    string str=Coeff.tostring(5);
			string str2=str.substr(1,str.length()-2);

			Coeff2.setlength(Frame.Params);
			int pos=0;
			int pos2=0;
			int i;
			 std::string::size_type sz;     // alias of size_t
			for(i=0;i<Coeff.length()-1;i++)
			{
				pos2=str2.find(",",pos);
				Coeff2[i]=std::stod (str2.substr(pos,pos2),&sz);
				pos=pos2+1;
			}

			Coeff2[i]=std::stod (str2.substr(pos),&sz);
			for(i=i+1;i<Coeff2.length();i++)
				Coeff2[i]=0;
			cout<<Coeff2.tostring(5)<<endl;


		minlmcreatev(Frame.NumberOfCameras*Frame.Pixels[0].size(), Coeff2, Frame.DiffStep, state);
		minlmsetcond(state, Frame.Epsg, Frame.Epsf, Frame.Epsx, Frame.MaxIts);

		switch(Coeff2.length()){

			case 1:
				scale="[1]";
				break;
			case 3:
				scale="[1,1,1]";
				break;

			case 5:
				scale="[1,1,1,1,1]";
				break;

			case 8:
				scale="[1,1,1,1,1,1,1,1]";
				break;

			case 12:
				scale="[1,1,1,1,1,1,1,1,1,1,1,1]";
				}

		minlmsetscale(state, scale);



		switch(Frame.NumberOfCameras){ //Determine camera-amount and use appropriate function to compute errors
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

		minlmresults(state, Coeff2, rep);
		cout<<Coeff2.tostring(5)<<endl;
	    return Coeff2;}

	   else{
	   return Coeff;}

	}
real_1d_array optimizeCoef2(real_1d_array Coeff, frameOptimizer Frame)
{ // Finds optimal coefficients for frame
  // Requires initial guess and frameOptimizer containing all necessary input

		// Optimization configuration parameters
		minlbfgsstate state;
	    minlbfgsreport rep;
	    void *ptr;
	    ptr=&Frame;

	    real_1d_array scale="[100.0,1,1]";
		// Initialize optimization
	    minlbfgscreatef(3, Coeff, Frame.DiffStep, state);
	    minlbfgssetcond(state, Frame.Epsg, Frame.Epsf, Frame.Epsx, Frame.MaxIts);

		minlbfgssetscale(state, scale);
		minlbfgssetprecscale( state);
		alglib::minlbfgsoptimize(state, combined_errorBF, NULL, ptr);
		minlbfgsresults(state, Coeff, rep);

	   cout<<Coeff.tostring(5).c_str()<<endl;
	    return Coeff;

	}
/// Main function
void Simulate(float std, float std_orig)
{ // Simulate inaccurate feature localization based on std for delta q and std_orig as deviation of initial guess on 'perfect solution'


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
	switch(s.NumberOfCameras){ // Choice of cameras and initialize each camera separately

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

	}

	vector<int> models={1,3,5,8,12}; // Surface models available

	for(int error=1; error<3;error++){ // Iterate over both error metrics

	s.ErrorMetric=error;
	for ( size_t m=1; m<models.size(); m++ ) { // Iterate over chosen surface models

		real_1d_array InitialCoeffs, InitialCoeffs2;
		for(double z=0; z<151;z=z+10){ // Adjust the original coefficient in the chosen surface model
		switch(models[m]){ // Initialize surface model

				case 1:
					InitialCoeffs= "[70.0]";
					break;
				case 3:
					InitialCoeffs= "[70.0,0.000,0.000]";
					break;
				case 5:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000]";
					break;
				case 8:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000]";
					break;
				case 12:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.0,0.0,0.0]";
					break;
				default:
					cout<<"Invalid surface model"<<endl;

			}

		// Adjust the original surface
		InitialCoeffs(0)=z;
		InitialCoeffs2=InitialCoeffs;

		// Choice of length scales
		float Ly=400;
		float Lx=300;

		// Lists (vectors) to store parameters and results
		vector<vector<real_1d_array> >CoefficientsList(s.ThreadAmount); // Output list of coefficients
		//vector<vector <double>> ErrorList(s.ThreadAmount); // Output list of corresponding errors
		//vector<vector <double>>ErrorList2(s.ThreadAmount); // Output list of corresponding errors2
		vector<frameOptimizer> optimizers(s.ThreadAmount); // List of optimizers (one for each thread)

		// Create original image points in image plane
		int rowAmount=40;
		int columnAmount=40;
		vector<Corner> temp;
		vector<Point3f> temp2;
		for(int i=0; i<rowAmount;i++)
			for(int j=0; j<columnAmount;j++){
			temp.push_back(Corner(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1),i*columnAmount+j)); // Create vector with Corner-objects
			temp2.push_back(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1)); // Create vector with Point3f-objects
			}

		// Compute inverse of camera matrices
		Mat RR=cam1.R.inv();
		Mat KK=cam1.K.inv();

		// Initialize required vectors
		vector<Point3f> qworld(temp.size());
		vector<Point3f> rays(temp.size());
		vector<Point3f> water(temp.size());
		vector<Point3f> normals2(temp.size());
		vector<Point3f> u(temp.size());
		vector<Point3f> v(temp.size());
		vector<double> thetai(temp.size());
		vector<Point3f> X(temp.size());
		vector<Mat> Rotations(temp.size());
		vector<double> thetad(temp.size());
		vector<float> a(temp.size());
		vector<Point3f> vv(temp.size());
		vector<Point3f> features(temp.size());

		// Compute world coordinates q'
		transform(temp2.begin(), temp2.end(), qworld.begin(), worldCoords<Point3f> (RR, cam1.T, KK));
		// Compute rays u=cq'
		transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (cam1.c));

		switch(models[m]){
			case 1:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_constant<Point3f> (InitialCoeffs, Lx,Ly));
			break;
			case 3:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_linear<Point3f> (InitialCoeffs, Lx,Ly));
			break;
			case 5:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_1storder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_1storder<Point3f> (InitialCoeffs, Lx,Ly));
			break;
			case 8:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_2ndorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_2ndorder<Point3f> (InitialCoeffs, Lx,Ly));
			break;
			case 12:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_3rdorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_3rdorder<Point3f> (InitialCoeffs, Lx,Ly));
			break;
		}
		// Easy way to compute normalized rays u
		transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>());
		// Compute angle between u and n2
		transform (u.begin(), u.end(), normals2.begin(), thetai.begin(),prod2<double>());
		// Compute theta_d based on Snell's law
		transform(thetai.begin(), thetai.end(), thetad.begin(), Angle2<double>(rw));
		// Compute rotation axis
		transform (normals2.begin(), normals2.end(), u.begin(), X.begin(),axis<Point3f> ());
		// Compute rotation matrix
		transform (thetad.begin(), thetad.end(), X.begin(), Rotations.begin(), Rotationmatrix<Mat>());
		// Rotate u over theta_d to get v
		transform (Rotations.begin(), Rotations.end(), u.begin(), v.begin(), Rotationmin<Point3f>());
		// Find vertical distance between surface point and z=0
		transform (water.begin(), water.end(), v.begin(), a.begin(), Distance<float>());
		// Multiply v with this distance
		transform (a.begin(), a.end(), v.begin(), vv.begin(), Waterray<Point3f>());
		// Surface points p + scaled v's give us intersection on plane F: points f!
		transform (water.begin(), water.end(), vv.begin(), features.begin(), std::plus<Point3f>());

		// Assign feature points to camera
		cam1.f=features;

		// Normal distribution for deviation of initial guess
		std::normal_distribution<float> dist(0, std_orig);

		for(size_t t=0;t<s.ThreadAmount;t++){ // Create per thread a vector with camera-parameters and insert in the corresponding optimizer-object
		vector<CameraParams> cams;
		cams.push_back(cam1);
		//cams.push_back(cam1); // In case of simulations with two cams --> also adjust 1 to 2 in optimizers below
		optimizers[t]=frameOptimizer(1, cams, error, models[m], Lx, Ly, s.Scaling, s.Epsg, s.Epsf, s.Epsx, s.MaxIts, s.DiffStep);
		}

		#pragma omp parallel // Go parallel over multiple threads
			{
			#pragma omp for nowait // Threads don't need each no other: no waiting
			for(size_t j=0;j<s.ThreadAmount;j++){ // Iterate over threads

				// Create vector of vector of image points (Corner objects)
				vector<vector<Corner> > pixels;
				pixels.push_back(temp);
				//pixels.push_back(temp); In case two cameras are used

				int tid = omp_get_thread_num();
				for(size_t i=0;i<20;i++){ // Iterate to get enough total number of simulations
					// Create random generator
					std::random_device rd;
					std::mt19937 gen(rd());

						switch(models[m]){ // Adjust original surface coefficients according to normal distribution as initial guess
								case 1:
									InitialCoeffs2(0)=InitialCoeffs(0)+dist(gen);
									break;
								case 3:
									InitialCoeffs2(0)=InitialCoeffs(0)+dist(gen);
									//InitialCoeffs2(3)=InitialCoeffs(3)+dist(gen);
									//InitialCoeffs(2)=dist(gen);
									break;
								case 5:
									InitialCoeffs2(0)=InitialCoeffs(0)+dist(gen);
									//InitialCoeffs2(3)=InitialCoeffs(3)+dist(gen);
									//InitialCoeffs(2)=dist(gen);
									//InitialCoeffs(3)=dist(gen);
									//InitialCoeffs(4)=dist(gen);
									break;
								case 8:
									InitialCoeffs2(0)=InitialCoeffs(0)+dist(gen);
									//InitialCoeffs2(3)=InitialCoeffs(3)+dist(gen);				//InitialCoeffs(2)=dist(gen);
									//InitialCoeffs(3)=dist(gen);
									//InitialCoeffs(4)=dist(gen);
									//InitialCoeffs(5)=dist(gen);
									//InitialCoeffs(6)=dist(gen);
									//InitialCoeffs(7)=dist(gen);
									break;
								case 12:
									InitialCoeffs2(0)=InitialCoeffs(0)+dist(gen);
									//InitialCoeffs2(3)=InitialCoeffs(3)+dist(gen);				//InitialCoeffs(2)=dist(gen);
									//InitialCoeffs(3)=dist(gen);
									//InitialCoeffs(4)=dist(gen);
									//InitialCoeffs(5)=dist(gen);
									//InitialCoeffs(6)=dist(gen);
									//InitialCoeffs(7)=dist(gen);
									//InitialCoeffs(8)=dist(gen);
									//InitialCoeffs(9)=dist(gen);
									//InitialCoeffs(10)=dist(gen);
									//InitialCoeffs(11)=dist(gen);
									break;
						}
				// Now add random noise to image points
				pixels[0]=ShuffleCorners(temp, std);
				//pixels[1]=ShuffleCorners(temp, std);

				// Change optimizer's image points to shuffled ones
				optimizers[tid].changePixels(pixels);

				// Run optimization
				real_1d_array coefficients= optimizeCoef(InitialCoeffs2, optimizers[tid]);

				// Compute errors corresponding to actual (based on non-shuffled q's) and found solution (based on shufled q's)
				//vector<double> E=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
				//vector<double> E2=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], InitialCoeffs);
			/*	double S,S2;
				// Compute mean over all features and save it into errorlist
				std::for_each(E.begin(), E.end(), [&] (double n) { S += (nancheck(n) ? 0.0 : n);});
				ErrorList[tid].push_back(S/E.size());
				std::for_each(E2.begin(), E2.end(), [&] (double n) { S2 += (nancheck(n) ? 0.0 : n);});
				ErrorList2[tid].push_back(S2/E2.size());*/

				CoefficientsList[tid].push_back(coefficients); // Append to coefficientlist


				}
			}
		}
	// Write out results
	string s="/home/lengelen/Documents/Simulations/COEFA00std1_dev1-FINAL-2cams/";
	string s2="Model-"+to_string(m)+"-ErrorMetric-"+to_string(error)+"-Height-"+to_string(z);
	//string s1="/home/lengelen/Documents/Simulations/C22/errorsBEST/Model-"+to_string(m)+"-ErrorMetric-"+to_string(error)+"-Height-"+to_string(z);
	//string s2="/home/lengelen/Documents/Simulations/C22/errorsINITIAL/Model-"+to_string(m)+"-ErrorMetric-"+to_string(error)+"-Height-"+to_string(z);

	writeArray(CoefficientsList, s2, s);
	//writeArrayErrors(ErrorList, s1);
	//writeArrayErrors(ErrorList2, s2);

	}
	}
}
}
void SimulateFeaturesF(float std)
{ // Simulate calibration errors by shifting feature points f with std

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

	int max_threads = omp_get_max_threads();

	if(s.ThreadAmount>max_threads)
	{
		cout<<"Amount of threads unavailable"<<endl;
		cout<<"Changed to maximum amount of threads: "<<max_threads<<endl;
		s.ThreadAmount=max_threads;
	}
	omp_set_dynamic(0);
	omp_set_num_threads(s.ThreadAmount);

	CameraParams cam1, cam2, cam3,cam4,cam5,cam6,cam7; // For each thread I need seperate camera
	cam1=Initialize(s, 1);
	cam2=Initialize(s, 1);
	cam3=Initialize(s, 1);
	cam4=Initialize(s, 1);
	cam5=Initialize(s, 1);
	cam6=Initialize(s, 1);
	cam7=Initialize(s, 1);

	vector<int> models={1,3,5,8,12}; // Available surface model
	for(int error=1; error<3;error++){ // Iterate over both error metrics

	s.ErrorMetric=error;
	for ( size_t m=0; m<models.size(); m++ ) { // Iterate over chosen surface models

		real_1d_array InitialCoeffs;
		switch(models[m]){ // Initialize surface model

				case 1:
					InitialCoeffs= "[70.0]";
					break;
				case 3:
					InitialCoeffs= "[70.0,0.000,0.000]";
					break;
				case 5:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000]";
					break;
				case 8:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000]";
					break;
				case 12:
					InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.0,0.0,0.0]";
					break;
				default:
					cout<<"Invalid surface model"<<endl;

			}


		for(double z=10; z<101;z=z+10){ // Adjust the original coefficient in the chosen surface model

			// Adjust the original surface
			InitialCoeffs(0)=z;

			// Choice of length scales
			double Ly=400;
			double Lx=300;

			// List (vector) to store parameters and results
			vector<vector<real_1d_array> >CoefficientsList(s.ThreadAmount); // Output list of coefficients
			vector<vector <double>> ErrorList(s.ThreadAmount); // Output list of corresponding errors
			vector<frameOptimizer> optimizers(s.ThreadAmount); // List of optimzers (one for each thread)

			// Create image point in image plane
			int rowAmount=40;
			int columnAmount=40;
			vector<Corner> temp;
			vector<Point3f> temp2;
			for(int i=0; i<rowAmount;i++)
				for(int j=0; j<columnAmount;j++){
				temp.push_back(Corner(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1),i*columnAmount+j)); // Create vector of Corner-objetcs
				temp2.push_back(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1)); // Create vector of Point3f
			}

			// Compute inverse of camera matrices
			Mat RR=cam1.R.inv();
			Mat KK=cam1.K.inv();

			// Initialize required vectors
			vector<Point3f> qworld(temp.size());
			vector<Point3f> rays(temp.size());
			vector<Point3f> water(temp.size());
			vector<Point3f> normals2(temp.size());
			vector<Point3f> u(temp.size());
			vector<Point3f> v(temp.size());
			vector<double> thetai(temp.size());
			vector<Point3f> X(temp.size());
			vector<Mat> Rotations(temp.size());
			vector<double> thetad(temp.size());
			vector<float> a(temp.size());
			vector<Point3f> vv(temp.size());
			vector<Point3f> features(temp.size());

			// Compute world coordinates q'
			transform(temp2.begin(), temp2.end(), qworld.begin(), worldCoords<Point3f> (RR, cam1.T, KK));
			// Compute rays u=cq'
			transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (cam1.c));

			switch(models[m]){
				case 1:
					// Determine surface points according to surface model and hypothesized coefficients
					transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					// Compute normals n2 according to surface model and hypothesized coefficients
					transform(water.begin(), water.end(), normals2.begin(), Normals_constant<Point3f> (InitialCoeffs, Lx,Ly));
				break;
				case 3:
					// Determine surface points according to surface model and hypothesized coefficients
					transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					// Compute normals n2 according to surface model and hypothesized coefficients
					transform(water.begin(), water.end(), normals2.begin(), Normals_linear<Point3f> (InitialCoeffs, Lx,Ly));
				break;
				case 5:
					// Determine surface points according to surface model and hypothesized coefficients
					transform(u.begin(), u.end(), water.begin(),  waterSurface_1storder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					// Compute normals n2 according to surface model and hypothesized coefficients
					transform(water.begin(), water.end(), normals2.begin(), Normals_1storder<Point3f> (InitialCoeffs, Lx,Ly));
				break;
				case 8:
					// Determine surface points according to surface model and hypothesized coefficients
					transform(u.begin(), u.end(), water.begin(),  waterSurface_2ndorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					// Compute normals n2 according to surface model and hypothesized coefficients
					transform(water.begin(), water.end(), normals2.begin(), Normals_2ndorder<Point3f> (InitialCoeffs, Lx,Ly));
				break;
				case 12:
					// Determine surface points according to surface model and hypothesized coefficients
					transform(u.begin(), u.end(), water.begin(),  waterSurface_3rdorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					// Compute normals n2 according to surface model and hypothesized coefficients
					transform(water.begin(), water.end(), normals2.begin(), Normals_3rdorder<Point3f> (InitialCoeffs, Lx,Ly));
				break;
			}

			// Easy way to compute normalized rays u
			transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>());
			// Compute angle between u and n2
			transform (u.begin(), u.end(), normals2.begin(), thetai.begin(),prod2<double>());
			// Compute theta_d based on Snell's law
			transform(thetai.begin(), thetai.end(), thetad.begin(), Angle2<double>(rw));
			// Compute rotation axis
			transform (normals2.begin(), normals2.end(), u.begin(), X.begin(),axis<Point3f> ());
			// Compute rotation matrix
			transform (thetad.begin(), thetad.end(), X.begin(), Rotations.begin(), Rotationmatrix<Mat>());
			// Rotate u over theta_d to get v
			transform (Rotations.begin(), Rotations.end(), u.begin(), v.begin(), Rotationmin<Point3f>());
			// Find vertical distance between surface point and z=0
			transform (water.begin(), water.end(), v.begin(), a.begin(), Distance<float>());
			// Multiply v with this distance
			transform (a.begin(), a.end(), v.begin(), vv.begin(), Waterray<Point3f>());
			// Surface points p + scaled v's give us intersection on plane F: points f!
			transform (water.begin(), water.end(), vv.begin(), features.begin(), std::plus<Point3f>());

			// Assign feature points to camera
			cam1.f=features;
			cam2.f=features;
			cam3.f=features;
			cam4.f=features;
			cam5.f=features;
			cam6.f=features;
			cam7.f=features;

			vector<vector<Corner> > pixels; // create vector<vector of Corner-objects from q's
			pixels.push_back(temp);

			#pragma omp parallel
				{
				#pragma omp for nowait
					for(size_t j=0;j<s.ThreadAmount;j++){ // Iterate over threads

				int tid = omp_get_thread_num();

				for(size_t i=0;i<30;i++){ // To have sufficient of total amount of simulations
					vector<CameraParams> cams;
					switch(tid){ // For each thread shuffle features f with normal distribution based on N(0,std)
					case 1:
					cam1.f=ShufflePoints(features, std);
					cams.push_back(cam1);
					break;
					case 2:
						cam2.f=ShufflePoints(features, std);
						cams.push_back(cam2);
					break;
					case 3:
						cam3.f=ShufflePoints(features, std);
						cams.push_back(cam3);
						break;
					case 4:
						cam4.f=ShufflePoints(features, std);
						cams.push_back(cam4);
					break;
					case 5:
						cam5.f=ShufflePoints(features, std);
						cams.push_back(cam5);
						break;
					case 6:
						cam6.f=ShufflePoints(features, std);
						cams.push_back(cam6);
						break;
					case 0:
						cam7.f=ShufflePoints(features, std);
						cams.push_back(cam7);
						break;}
				// Reset optimizer
				optimizers[tid]=frameOptimizer(1, cams, error, models[m], Lx, Ly, s.Scaling, s.Epsg, s.Epsf, s.Epsx, s.MaxIts, s.DiffStep);
				optimizers[tid].changePixels(pixels);

				// Perform optimization
				real_1d_array coefficients= optimizeCoef(InitialCoeffs, optimizers[tid]);

				// Optional: calculate error of found solution
				//vector<double> E1=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, (int) InitialCoeffs.length(), optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], InitialCoeffs);
				//double sum = std::accumulate(E1.begin(), E1.end(), 0.0);
				//double mean = sum / E1.size();
				//vector<double> E2=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, (int) InitialCoeffs.length(), optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
				//double sum2 = std::accumulate(E2.begin(), E2.end(), 0.0);
				//double mean2 = sum2 / E2.size();

				CoefficientsList[tid].push_back(coefficients); // Append solution to coefficientlist
				}
			}
		}
	// Write out results
	string s="/home/lengelen/Documents/Simulations/ShuffleF_std02/";
	string s2="Model-"+to_string(m)+"-ErrorMetric-"+to_string(error)+"-Height-"+to_string(z);

	writeArray(CoefficientsList, s2, s);
	}
	}
}
}
void DetermineArea()
{ // Determine area over which simulation results need to be evaluated
  // Gives four extreme surface points as return

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

		switch(s.NumberOfCameras){ // Choice of cameras and initialize each camera separately

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

		}

		vector<int> models={1,3,5,8,12}; // Available surface models
		int m=0; // Choice surface model

		real_1d_array InitialCoeffs;
		for(double z=10; z<151;z=z+10){ // Adjust the original coefficient in the chosen surface model
			switch(models[m]){ // Initialize surface model
					case 1:
						InitialCoeffs= "[70.0]";
						break;
					case 3:
						InitialCoeffs= "[70.0,0.000,0.000]";
						break;
					case 5:
						InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000]";
						break;
					case 8:
						InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000]";
						break;
					case 12:
						InitialCoeffs= "[70.0,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.0,0.0,0.0,0.0]";
						break;
					default:
						cout<<"Invalid surface model"<<endl;

				}

			// Adjust original surface
			InitialCoeffs(0)=z;

			// Choice of length scales
			double Ly=400;
			double Lx=300;

			// Create image points in image plane
			int rowAmount=40;
			int columnAmount=40;
			vector<Corner> temp;
			vector<Point3f> temp2;
			for(int i=0; i<rowAmount;i++)
				for(int j=0; j<columnAmount;j++){
					temp.push_back(Corner(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1),i*columnAmount+j));
					temp2.push_back(Point3f(90+z*0.8+j*30,150+z*0.3+i*20,1));
					}

			// Compute inverse of camera matrices
			Mat RR=cam1.R.inv();
			Mat KK=cam1.K.inv();

			// Create required temporary vectors
			vector<Point3f> qworld(temp.size());
			vector<Point3f> rays(temp.size());
			vector<Point3f> water(temp.size());
			vector<Point3f> u(temp.size());

			// Compute world coordinates q'
			transform(temp2.begin(), temp2.end(), qworld.begin(), worldCoords<Point3f> (RR, cam1.T, KK));
			// Compute rays u=cq'
			transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (cam1.c));

			switch(models[m]){
					case 1:
						// Determine surface points according to surface model and hypothesized coefficients
						transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					break;
					case 3:
						// Determine surface points according to surface model and hypothesized coefficients
						transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					break;
					case 5:
						// Determine surface points according to surface model and hypothesized coefficients
						transform(u.begin(), u.end(), water.begin(),  waterSurface_1storder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					break;
					case 8:
						// Determine surface points according to surface model and hypothesized coefficients
						transform(u.begin(), u.end(), water.begin(),  waterSurface_2ndorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));
					break;
					case 12:
						// Determine surface points according to surface model and hypothesized coefficients
						transform(u.begin(), u.end(), water.begin(),  waterSurface_3rdorder<Point3f> (cam1.c, Lx, Ly, InitialCoeffs));

					break;
				}
			// Find 4 extreme corners of grid on water surface
			vector<Point3f> ends;
			ends.push_back(water[0]);
			ends.push_back(water[columnAmount-1]);
			ends.push_back(water[water.size()-columnAmount]);
			ends.push_back(water[water.size()-1]);

			string s="/home/lengelen/Documents/Simulations/SurfaceA00_2/";
			string s2="Height-"+to_string(z);
			writePointsToFile(ends, s2, s);
			}
}
int main ( )
{ // Main function to reconstruct time-dependent water surface.


	// Set chrono-time at begin
	std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

	//Read all settings
	Settings s;
	
	string inputSettingsFile;
	cout << "Give InputsettingsFile: " << flush; //Read which settingsfile has to be used
	getline( cin, inputSettingsFile );  // gets everything the user ENTERs
	FileStorage fs(inputSettingsFile, FileStorage::READ); // Read the settings
	if (!fs.isOpened())
		{
			cout << "Could not open the configuration file: "<< endl;
			exit(1);
		 }
	fs["Settings"] >> s;
	fs.release();   // close Settings file
	CameraParams cam1, cam2, cam3;
	int max_threads = omp_get_max_threads();
	string res="RESULTS/";

	if(s.ThreadAmount>max_threads)
	{
		cout<<"Amount of threads unavailable"<<endl;
		cout<<"Changed to maximum amount of threads: "<<max_threads<<endl;
		s.ThreadAmount=max_threads;
	}
	omp_set_dynamic(0);
	omp_set_num_threads(s.ThreadAmount);
	vector<vector<vector<Corner> > > prevPixels(s.ThreadAmount);

	switch(s.NumberOfCameras){ // Choice of cameras and initialize each camera separately

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

	vector<vector<real_1d_array> >CoefficientsList(s.ThreadAmount); // Output list of coefficients
	vector<vector <double>> ErrorList(s.ThreadAmount); // Output list of corresponding errors

	real_1d_array InitialCoeffs=s.InitialGuess.c_str();
	vector<frameOptimizer> optimizers(s.ThreadAmount);

	for(size_t i=0;i<s.ThreadAmount;i++){
	vector<CameraParams> cams;
	switch(s.NumberOfCameras){ //Choice of cameras and initialize optimization with appropriate camera amount

	case 1:
		cams.push_back(cam1);
		optimizers[i]=frameOptimizer(1, cams, s.ErrorMetric, s.SurfaceModelParameters, s.Lx, s.Ly, s.Scaling, s.Epsg, s.Epsf, s.Epsx, s.MaxIts, s.DiffStep);
		break;
	case 2:
		cams.push_back(cam1);
		cams.push_back(cam2);
		optimizers[i]=frameOptimizer(2, cams, s.ErrorMetric, s.SurfaceModelParameters, s.Lx, s.Ly, s.Scaling, s.Epsg, s.Epsf, s.Epsx, s.MaxIts, s.DiffStep);
		break;
	case 3:
		cams.push_back(cam1);
		cams.push_back(cam2);
		cams.push_back(cam3);
		optimizers[i]=frameOptimizer(3, cams, s.ErrorMetric, s.SurfaceModelParameters, s.Lx, s.Ly, s.Scaling, s.Epsg, s.Epsf, s.Epsx, s.MaxIts, s.DiffStep);
		break;
	default:
		cout<<"Invalid amount of camera's"<<endl;
		return 0;
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

	switch(s.NumberOfCameras){ // Choice of cameras and read image-names in every directory (for each camera)

	case 1:
		InputDirectory1=s.InputDirectory1;
		check1=readStringList(imageList1, InputDirectory1); //read all image-names
		check=check1;
		break;
	case 2:
		InputDirectory1=s.InputDirectory1;
		InputDirectory2=s.InputDirectory2;
		check1=readStringList(imageList1, InputDirectory1); //read all image-names
		check2=readStringList(imageList2, InputDirectory2); //read all image-names
		check=check1 && check2;
		break;
	case 3:
		InputDirectory1=s.InputDirectory1;
		InputDirectory2=s.InputDirectory2;
		InputDirectory3=s.InputDirectory3;
		check1=readStringList(imageList1, InputDirectory1); //read all image-names
		check2=readStringList(imageList2, InputDirectory2); //read all image-names
		check3=readStringList(imageList3, InputDirectory3); //read all image-names
		check=check1 && check2 && check3;
		break;
	}

	// Check if image-names are successfully loaded
	if(check){
		cout<<"image-names/filenames loaded"<<endl;
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
		{
			int tid = omp_get_thread_num();

			vector<vector<Corner> > pixels(s.NumberOfCameras);
			// Detect corners in image i
			if(i-j*images_thread>0.1){
			
			switch(s.NumberOfCameras){

							case 1:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1, res+s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][0]);
								break;
							case 2:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1,  res+ s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][0]);
								pixels[1]=readFeaturesImage( imageList2[nimages-i-1], cam2,   res+s.OutputDirectory2+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][1]);
								break;
							case 3:
								pixels[0]=readFeaturesImage( imageList1[nimages-i-1], cam1,  res+ s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][0]);
								pixels[1]=readFeaturesImage( imageList1[nimages-i-1], cam2,  res+ s.OutputDirectory2+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][1]);
								pixels[2]=readFeaturesImage( imageList3[nimages-i-1], cam3,  res+ s.OutputDirectory3+"/Frame_"+ToString(i)+".txt", s, prevPixels[tid][2]);
								break;
							default:
								cout<<"Error camera amount"<<endl;
								exit(1);
						}
			}else{
				switch(s.NumberOfCameras){

							case 1:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1, res+ s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s);
								break;
							case 2:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1,  res+s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s);
								pixels[1]=readFeaturesFirstImage( imageList2[nimages-i-1], cam2,  res+s.OutputDirectory2+"/Frame_"+ToString(i)+".txt", s);
								break;
							case 3:
								pixels[0]=readFeaturesFirstImage( imageList1[nimages-i-1], cam1, s.OutputDirectory1+"/Frame_"+ToString(i)+".txt", s);
								pixels[1]=readFeaturesFirstImage( imageList2[nimages-i-1], cam2, s.OutputDirectory2+"/Frame_"+ToString(i)+".txt", s);
								pixels[2]=readFeaturesFirstImage( imageList3[nimages-i-1], cam3, s.OutputDirectory3+"/Frame_"+ToString(i)+".txt", s);
								break;
							default:
								cout<<"Error camera amount"<<endl;
								exit(1);
						}
			}
			prevPixels[tid]=pixels;
			vector<double> E1, E2, E3, E4, E5, E6; // Errors corresponding to the previous coefficient-sets
			optimizers[tid].changePixels(pixels); // Change frameOptimizer to corner-list which requires optimization for this frame
			
			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				E1=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, (int) InitialCoeffs.length(), optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], InitialCoeffs);
				E2=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-1]);
				E3=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-2]);
				E4=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-3]);
				E5=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-4]);
				E6=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], CoefficientsList[tid][i-j*images_thread-5]);

			}

			if(s.NumberOfCameras>1){ // Procedure for more than 1 camera
				
			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				vector<double>E11=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, (int) InitialCoeffs.length(),optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], InitialCoeffs);
				E1.insert( E1.end(), E11.begin(), E11.end() );
				vector<double>E22=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-1]);
				E2.insert( E2.end(), E22.begin(), E22.end() );
				vector<double>E33=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-2]);
				E3.insert( E3.end(), E33.begin(), E33.end() );
				vector<double>E44=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-3]);
				E4.insert( E4.end(), E44.begin(), E44.end() );
				vector<double>E55=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-4]);
				E5.insert( E5.end(), E55.begin(), E55.end() );
				vector<double>E66=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], CoefficientsList[tid][i-j*images_thread-5]);
				E6.insert( E6.end(), E66.begin(), E66.end() );
			}
			}
			if(s.NumberOfCameras>2){ // Procedure for 3 cameras

			// Starting from 6th image: find best initial guess out of 5 previous time step and first initial guess (flat surface)
			if(i-j*images_thread>4){
				vector<double>E111=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, (int) InitialCoeffs.length(),optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], InitialCoeffs);
				E1.insert( E1.end(), E111.begin(), E111.end() );
				vector<double>E222=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-1]);
				E2.insert( E2.end(), E222.begin(), E222.end() );
				vector<double>E333=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-2]);
				E3.insert( E3.end(), E333.begin(), E333.end() );
				vector<double>E444=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-3]);
				E4.insert( E4.end(), E444.begin(), E444.end() );
				vector<double>E555=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-4]);
				E5.insert( E5.end(), E555.begin(), E555.end() );
				vector<double>E666=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], CoefficientsList[tid][i-j*images_thread-5]);
				E6.insert( E6.end(), E666.begin(), E666.end() );
				}
			}

			// Do optimization with list of sorted corner-sets
			real_1d_array coefficients;
			
			if(i-j*images_thread<1){ //For images 0-4: use predetermined initial guess
				coefficients = optimizeCoef(InitialCoeffs, optimizers[tid]);
						}

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

			if(s.SaveResiduals)
			{ // Save errors corresponding to optimized coefficients
				
				// Compute errors corresponding to optimized coefficients for camera 1
				vector<double> Error;
				switch (s.NumberOfCameras){
				case 1:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric,optimizers[tid].Params, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					break;}
				case 2:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric,optimizers[tid].Params, optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					vector<double> Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					break;}
				case 3:{
					Error=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric, optimizers[tid].Params,optimizers[tid].Cameras[0], optimizers[tid].Pixels[0], optimizers[tid].fs[0], coefficients);
					vector<double> Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric,optimizers[tid].Params, optimizers[tid].Cameras[1], optimizers[tid].Pixels[1], optimizers[tid].fs[1], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					Temp=compute_error(optimizers[tid].Lx, optimizers[tid].Ly, optimizers[tid].ErrorMetric,optimizers[tid].Params, optimizers[tid].Cameras[2], optimizers[tid].Pixels[2], optimizers[tid].fs[2], coefficients);
					Error.insert( Error.end(), Temp.begin(), Temp.end() ); //Concatenate errors of camera 2 to total error-vector
					break;}
				default:
					cout<<"Error camera amount"<<endl;
					exit(1);
				}
			double S;
				// Compute mean over all features and save it into errorlist
				std::for_each(Error.begin(), Error.end(), [&] (double n) { S += (nancheck(n) ? 0.0 : n);}); 
				ErrorList[tid].push_back(S/Error.size());}

		}
			}}


	//Write results to files

	writeArray(CoefficientsList, s.OutputFileName, res);
	if(s.SaveResiduals)
	{writeArrayErrors(ErrorList, s.OutputFileNameResiudals, res);}
	std::chrono::steady_clock::time_point end1= std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() <<std::endl;

	timestamp ( );
    cout << "\n";
    cout << "Surface reconstruction \n";
    cout << "C++ version:\n";
    cout << "By Lukas Engelen:\n";
    return 0;
}


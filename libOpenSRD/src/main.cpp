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
{  // Computes normalized n* of set of 1 coefficients at position of surface point (3D Float)
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

            	Point3f n =Point3f(-0.16, 0.07, 1);
            	 return n/norm(n);
            }
        };
template <class Type> class Normals_linear
{  // Computes normalized n* of set of 3 coefficients at position of surface point (3D Float)
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
{  // Computes normalized n* of set of 8 coefficients at position of surface point (3D Float)
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
{  // Computes normalized n* of set of 8 coefficients at position of surface point (3D Float)
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
{  // Computes normalized n* of set of 12 coefficients at position of surface point (3D Float)
        private:
            real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
        public:
            // Constructor
            Normals_3rdorder (real_1d_array i_coef, double i_Lx, double i_Ly) {
                coef=i_coef;
                Lx=i_Ly;
                Ly=i_Ly;
            }

            // The function call
            Type operator ( ) ( Point3f& elem ) const
            {


               Point3f n= Point3f(3*PI_F/Lx*coef[8]*sin(3*PI_F*elem.x/Lx)+2*PI_F/Lx*coef[5]*sin(2*PI_F*elem.x/Lx)+2*PI_F/Lx*coef[10]*sin(2*PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[7]*sin(PI_F*elem.x/Lx)*cos(PI_F*elem.y/Ly)+PI_F/Lx*coef[3]*sin(PI_F*elem.x/Lx)-coef[1]/Lx, 3*PI_F/Ly*coef[9]*sin(3*PI_F*elem.y/Ly)+2*PI_F/Ly*coef[6]*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[7]*cos(PI_F*elem.x/Lx)*sin(PI_F*elem.y/Ly)+2*coef[11]*cos(PI_F*elem.x/Lx)*sin(2*PI_F*elem.y/Ly)+PI_F/Ly*coef[4]*sin(PI_F*elem.y/Ly)-coef[2]/Ly, 1);
               return n/norm(n);

            }
        };
template <class Type> class Normals_Paper
{  // Computes normalized n* of set of 2 coefficients at position of surface point (3D Float), as used in Paper
private:
         real_1d_array coef; double Lx; double Ly;   // Characterization of surface model
     public:
         // Constructor
         Normals_Paper (real_1d_array i_coef, double i_Lx, double i_Ly) {
             coef=i_coef;
             Lx=i_Lx;
             Ly=i_Ly;
         }

         // The function call
         Type operator ( ) ( Point3f& elem ) const
         {
        	 // Scaled x-component with 100 because this improved performance of optimization algorithm: changes are same order as changes in water depth
        	 // Fixed y-component, calculated based on average image with linear model
        	 Point3f n= Point3f(-coef[1]/100,0.07, 1);

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
	real_1d_array coef; Point3f c;  double Lx; double Ly; Point2f Middle;// Characterization of surface model

	public:
		// Constructor
	waterSurface_linear (Point3f i_c, double i_Lx, double i_Ly, Point2f i_Middle, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
			Middle=i_Middle;
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
			x[0] = 500;
			x[1] = 0;
			x[2] = 0;
			x[3] = 60;

			// Additional input information passed by pointer
			parameters[0] = c.x;
			parameters[1] = c.y;
			parameters[2] = c.z;
			parameters[3] = ray.x;
			parameters[4] = ray.y;
			parameters[5] = ray.z;
			parameters[6] = Lx;
			parameters[7] = Ly;
			parameters[8] = Middle.x;
			parameters[9] = Middle.y;
			parameters[10] = coef[0];
			parameters[11] = coef[1];
			parameters[12] = coef[2];

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
template <class Type> class waterSurface_Paper
{ // Finds surface point for set of 2 coefficients and ray u=cq' (all 3D Float), as used in Paper

private:
	real_1d_array coef; Point3f c;  double Lx; double Ly; Point2f Middle;// Characterization of surface model

	public:
		// Constructor
	waterSurface_Paper (Point3f i_c, double i_Lx, double i_Ly, Point2f i_Middle, real_1d_array i_coef) {
			coef=i_coef;
			c=i_c;
			Lx=i_Lx;
			Ly=i_Ly;
			Middle=i_Middle;
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
			parameters = new double [12];

			// Initial guess
			x[0] = 1000;
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
			parameters[8] = Middle.x;
			parameters[9] = Middle.y;
			parameters[10] = coef[0];
			parameters[11] = coef[1];

			// Initialization
			findIntersection_Paper( n, x, fvec, &iflag, parameters );

			// Optimization
			int info = hybrd1 ( findIntersection_Paper, n, x, fvec, parameters, tol, wa, lwa );
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
template <class T> struct prod_reverse
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
    	return 980*acos(n1.ddot(n2));
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

struct minX
{ ///sort image points according to increasing u
 bool operator() ( Point3f a, Point3f b ){
               return a.x <= b.x;
    }
} sortingminX;
struct maxX
{ ///sort image points according to increasing u
bool operator() ( Point3f a, Point3f b ){
              return b.x <= a.x;
   }
} sortingmaxX;
struct minY
{ ///sort image points according to increasing u
bool operator() ( Point3f a, Point3f b ){
              return a.y <= b.y;
   }
} sortingminY;
struct maxY
{ ///sort image points according to increasing u
bool operator() ( Point3f a, Point3f b ){
              return b.y <= a.y;
   }
} sortingmaxY;
void createFeatureFile(){ // Creates text file with theoretical location of reference pattern
	ofstream fout("/home/lengelen/Documents/Doctoraat/Measurements/Ms20170508/Ref/VerticesPoseEstimation2-2.txt");
	 if(!fout)
				    {
				        cout<<"File Not Opened"<<endl;  return;
				    }

	 for(int j=0;j<10;j++)
		 for(int i=0;i<39;i++){

    		fout << 390-i*10<<"\t"<<10+j*10<<"\t"<<0<<endl;


		}
	fout.close();
}

CameraParams Initialize(Settings s, int camera)
{ // Initialization of camera parameters in a CameraParams-object for input of settings and camera number

	// Initialization and reading of all settings
	Mat cameraMatrix, distCoeffs;
	string CalibrationFile,  InputReference, OutputFile, CameraPoseFile;
	Mat R, origin, RR, rvec, tvec;
	Point3f c;
	switch(camera){ // Open calibration file and reference image depending on camera
		case 1:
			 CalibrationFile= s.CalibrationFile1;
			 InputReference= s.InputReference1;
			 CameraPoseFile= s.InputInitial1;
			 OutputFile=s.OutputCameraPose1;
			 break;
		case 2:
			 CalibrationFile= s.CalibrationFile2;
			 InputReference= s.InputReference2;
			 CameraPoseFile= s.InputInitial2;
			 OutputFile=s.OutputCameraPose2;
			 break;
		case 3:
			 CalibrationFile=s.CalibrationFile3;
			 InputReference= s.InputReference3;
			 CameraPoseFile= s.InputInitial3;
			 OutputFile=s.OutputCameraPose3;
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
			vector<Point2f> sortedcorners =create_points_OpenCV(image, ResponseThreshold, MinDistance, s.ResponseRadius, board_sz, cameraMatrix, distCoeffs, s.ShowCorners);

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

		if(s.SaveCameraPose){ // Save result camera pose estimation if needed

		saveCameraParams(OutputFile,"", R, tvec);

		}
		CameraParams cam=CameraParams(R, tvec, cameraMatrix, c, distCoeffs);
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
{ // Find intersection of ray cq' with surface model characterized by 3 coefficients and centre of reconstructed area
    double cx = params[0];
    double cy = params[1];
    double cz = params[2];
    double rx = params[3];
    double ry = params[4];
    double rz = params[5];
    double Lx = params[6];
    double Ly = params[7];
    double centerx = params[8];
    double centery = params[9];
    double A00 = params[10];
    double B = params[11];
    double C = params[12];

    fvec[0] = cx+ x[0] * rx -x[1];
    fvec[1] = cy + x[0] * ry -x[2];
    fvec[2] = cz + x[0] * rz -x[3];
    fvec[3] = -x[3] + + A00 +B*(x[1]-centerx)/Lx+C*(x[2]-centery)/Ly;
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
    double Lx = 25;
    double Ly = 25;
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
    fvec[3] = -x[3] + A00 + A10 * cos(PI*x[1]/Lx)+A01*cos(PI*x[2]/Ly)+B*(x[1]-270)/Lx+C*(x[2]-355)/Ly + A20 * cos(2*PI*x[1]/Lx)+ A02 * cos(2*PI*x[2]/Ly) + A11 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly)+ A30 * cos(3*PI*x[1]/Lx)+ A03 * cos(3*PI*x[2]/Ly) + A21 * cos(PI*x[1]/Lx*2) * cos(PI*x[2]/Ly)+ A12 * cos(PI*x[1]/Lx) * cos(PI*x[2]/Ly*2);
    return;
}
void findIntersection_Paper( int n, double x[], double fvec[], int *iflag, double params[] )
{ // Find intersection of ray cq' with surface model characterized by 2 coefficients, as used in Paper
	double cx = params[0];
	double cy = params[1];
	double cz = params[2];
	double rx = params[3];
	double ry = params[4];
	double rz = params[5];
	double Lx = params[6];
	double Ly = params[7];
	double centerx = params[8];
	double centery = params[9];
	double A00 = params[10];
	double C = params[11];


	fvec[0] = cx+ x[0] * rx -x[1];
	fvec[1] = cy + x[0] * ry -x[2];
	fvec[2] = cz + x[0] * rz -x[3];
	fvec[3] = -x[3] + A00;
	return;
}
template <class Type> class crossVerify
{
        private:
         	float gridSize; uchar minThreshold; uchar maxThreshold;
        public:
            // Constructor
         	crossVerify (float i_gridSize, uchar i_minThreshold,  uchar i_maxThreshold) {

         		gridSize=i_gridSize;
                minThreshold=i_minThreshold;
                maxThreshold=i_maxThreshold;
            }

            // The function call
            Type operator ( ) ( Point3f& f, Point3f& pixel ) const
            {
            	if(fmod(floor(f.x/gridSize),2)==fmod(floor((f.y)/gridSize),2)){
            		if(pixel.z<minThreshold){
							return 1;}
            		else return (0);}
				else{
						if(pixel.z>maxThreshold){
							return 1;
						}
						else return 0;
					}
                	}
};
template <class Type> class crossVerify2
{
        private:
         	float gridSize; uchar minThreshold; uchar maxThreshold;
        public:
            // Constructor
         	crossVerify2 (float i_gridSize, uchar i_minThreshold,  uchar i_maxThreshold) {

         		gridSize=i_gridSize;
                minThreshold=i_minThreshold;
                maxThreshold=i_maxThreshold;
            }

            // The function call
            Type operator ( ) ( Point3f& f, Point3f& pixel ) const
            {
            	if(fmod(floor(f.x/gridSize),2)==fmod(floor(f.y/gridSize),2)){
            	            		if(pixel.z<minThreshold){
            								return 0;}
            	            		else return (1);}
            					else{
            							if(pixel.z>maxThreshold){
            								return 0;

            							}
            							else return 1;
            						}

                    	}
        };
vector<uchar> errorfunction(bool flag, float Lx, float Ly, int gridSize, int ParameterAmount, CameraParams Camera, Mat Image, real_1d_array Coeff, Point2f Middle,  uchar i_minThreshold, uchar i_maxThreshold, int min_u, int max_u, int min_v, int max_v)
{
	vector<Point3f> temp;
	vector<Point3f> grayScales;

	for(int i=min_u; i<max_u;i++){
			for(int j=min_v ;j<max_v;j++){
			if(Image.at<uchar>(j,i)< i_minThreshold || Image.at<uchar>(j,i)> i_maxThreshold){
				temp.push_back(Point3f(i,j,1));
				grayScales.push_back(Point3f(i,j,Image.at<uchar>(j,i)));
			}
		}}
	// Compute inverse of camera matrices
	Mat RR=Camera.R.inv();
	Mat KK=Camera.K.inv();

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
	vector<uchar> errors(temp.size());

	// Compute world coordinates q'
	transform(temp.begin(), temp.end(), qworld.begin(), worldCoords<Point3f> (RR, Camera.T, KK));
	// Compute rays u=cq'
	transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (Camera.c));

	switch(Coeff.length()){
			case 1:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (Camera.c, Lx, Ly, Coeff));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_constant<Point3f> (Coeff, Lx,Ly));
			break;
			case 3:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (Camera.c, Lx, Ly, Middle, Coeff));
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
			case 2:
			// Determine surface points according to surface model and hypothesized coefficients
			transform(u.begin(), u.end(), water.begin(),  waterSurface_Paper<Point3f> (Camera.c, Lx, Ly, Middle, Coeff));
			// Compute normals n2 according to surface model and hypothesized coefficients
			transform(water.begin(), water.end(), normals2.begin(), Normals_Paper<Point3f> (Coeff, Lx,Ly));
	}

	// Easy way to compute normalized rays u
	transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>());
	// Compute angle between u and n2
	transform (u.begin(), u.end(), normals2.begin(), thetai.begin(),prod_reverse<double>());
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
	transform (features.begin(), features.end(), grayScales.begin(),errors.begin(), crossVerify<uchar> (gridSize, i_minThreshold,  i_maxThreshold));

	// transform (features.begin(), features.end(), grayScales.begin(),errors.begin(), crossVerify2<uchar> (gridSize, i_minThreshold,  i_maxThreshold));

	//flag=true;
	if(flag){

		Mat error=Mat::zeros(Image.rows, Image.cols, CV_8U);
		Mat projection=Mat::zeros(Image.rows, Image.cols, CV_8U);
		for(int j=0;j<errors.size();j++){
			if(errors[j]!=0) error.at<uchar>(temp[j].y,temp[j].x)=(uchar)errors[j];
			else error.at<uchar>(temp[j].y,temp[j].x)=0;
			if(features[j].x<0 ||features[j].y<0)projection.at<uchar>(temp[j].y,temp[j].x)=0;
					else{
		if(fmod(floor(features[j].x/gridSize),2)==fmod(floor((features[j].y)/gridSize),2)){
					projection.at<uchar>(temp[j].y,temp[j].x)=255;}
							else projection.at<uchar>(temp[j].y,temp[j].x)=100;
		}}
					  // Show our image inside it.
		namedWindow( "Display image 2", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display image 2", 1000, 800);
		imshow( "Display image 2", Image );                   // Show our image inside it.
		namedWindow( "Display error 2", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display error 2", 1000, 800);
		imshow( "Display error 2", error*5 );                   // Show our image inside it.
		namedWindow( "Display projection", WINDOW_NORMAL );// Create a window for display.
			resizeWindow("Display projection", 1000, 800);
			imshow( "Display projection", projection );                   // Show our image inside it.

		waitKey(0);             // Wait for a keystroke in the window

	}

	//double res=accumulate( errors.begin(), errors.end(), 0.0);
	//cout<<res<<endl;
	return errors;

}
Point2f computeMiddle(bool flag, float Lx, float Ly, int gridSize, int ParameterAmount, CameraParams Camera, Mat Image, real_1d_array Coeff, uchar i_minThreshold,  uchar i_maxThreshold, int min_u, int max_u, int min_v, int max_v)
{ // Computes approximate middle of reconstructed area, used in paper during validation on still water

	vector<Point3f> temp;
	vector<Point3f> grayScales;
	float minX=0;
	float minY=0;
	float maxX=0;
	float maxY=0;
	for(int i=min_u; i<max_u;i++){
				for(int j=min_v ;j<max_v;j++){
			if(Image.at<uchar>(j,i)< i_minThreshold || Image.at<uchar>(j,i)> i_maxThreshold){
				temp.push_back(Point3f(i,j,1));
				grayScales.push_back(Point3f(i,j,Image.at<uchar>(j,i)));
			}
		}}
	// Compute inverse of camera matrices
	Mat RR=Camera.R.inv();
	Mat KK=Camera.K.inv();

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
	vector<uchar> errors(temp.size());

	// Compute world coordinates q'
	transform(temp.begin(), temp.end(), qworld.begin(), worldCoords<Point3f> (RR, Camera.T, KK));
	// Compute rays u=cq'
	transform(qworld.begin(), qworld.end(), u.begin(), Ray<Point3f> (Camera.c));

	switch(Coeff.length()){
			case 1:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_constant<Point3f> (Camera.c, Lx, Ly, Coeff));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_constant<Point3f> (Coeff, Lx,Ly));
			break;
			case 3:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_linear<Point3f> (Camera.c, Lx, Ly, Point2f(0,0), Coeff));
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
			case 2:
				// Determine surface points according to surface model and hypothesized coefficients
				transform(u.begin(), u.end(), water.begin(),  waterSurface_Paper<Point3f> (Camera.c, Lx, Ly, Point2f(0,0), Coeff));
				// Compute normals n2 according to surface model and hypothesized coefficients
				transform(water.begin(), water.end(), normals2.begin(), Normals_Paper<Point3f> (Coeff, Lx,Ly));

		}

	// Easy way to compute normalized rays u
	transform (water.begin(), water.end(), qworld.begin(), u.begin(), SubtractNormalized<Point3f>());
	// Compute angle between u and n2
	transform (u.begin(), u.end(), normals2.begin(), thetai.begin(),prod_reverse<double>());
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
	transform (features.begin(), features.end(), grayScales.begin(),errors.begin(), crossVerify<uchar> (gridSize, i_minThreshold,  i_maxThreshold));
	//transform (features.begin(), features.end(), grayScales.begin(),errors.begin(), crossVerify2<uchar> (gridSize, i_minThreshold,  i_maxThreshold));


	if(flag){

		Mat error=Mat::zeros(Image.rows, Image.cols, CV_8U);
		for(size_t j=0;j<errors.size();j++){
			if(errors[j]==1) error.at<uchar>(temp[j].y,temp[j].x)=(uchar)errors[j]*4;
			else error.at<uchar>(temp[j].y,temp[j].x)=1;
		}

					  // Show our image inside it.
		namedWindow( "Display image 2", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display image 2", 2000, 1000);
		imshow( "Display image 2", Image );                   // Show our image inside it.
		namedWindow( "Display error 2", WINDOW_NORMAL );// Create a window for display.
		resizeWindow("Display error 2", 2000, 1000);
		imshow( "Display error 2", error*60 );                   // Show our image inside it.

		waitKey(0);             // Wait for a keystroke in the window

	}

	sort(water.begin(), water.end(), sortingminX); //sort all points on one row by increasing u
	minX=water[0].x;
	sort(water.begin(), water.end(), sortingminY); //sort all points on one row by increasing u
	minY=water[0].y;
	sort(water.begin(), water.end(), sortingmaxX); //sort all points on one row by increasing u
	maxX=water[0].x;
	sort(water.begin(), water.end(), sortingmaxY); //sort all points on one row by increasing u
	maxY=water[0].y;

	float midX=(minX+maxX)/2;
	float midY=(minY+maxY)/2;
	Point2f m=Point2f(midX,midY);
	//cout<<m<<endl;
	return m;

}
void combinedErrorFunction(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
	imageOptimizer Frame=*(static_cast<imageOptimizer*>(ptr)); // Get ImageOptimzer from pointer
	vector<uchar> E1=errorfunction(false, Frame.Lx, Frame.Ly, Frame.gridSize, Frame.Params, Frame.Cameras[0], Frame.Images[0], x, Frame.centerPoint, Frame.minThreshold, Frame.maxThreshold, Frame.min_u, Frame.max_u, Frame.min_v, Frame.max_v);
	vector<uchar> E2;
	vector<uchar> E3;
	for(size_t t=0;t<E1.size();t++)
		fi[t]=E1[t];
	if(Frame.NumberOfCameras>1){
	E2=errorfunction(false, Frame.Lx, Frame.Ly, Frame.gridSize, Frame.Params, Frame.Cameras[1], Frame.Images[1], x, Frame.centerPoint, Frame.minThreshold, Frame.maxThreshold, Frame.min_u, Frame.max_u, Frame.min_v, Frame.max_v);

	for(size_t t=E1.size();t<E1.size()+E2.size();t++)
		fi[t]=E2[t-E1.size()];}

	if(Frame.NumberOfCameras>2){
		E3=errorfunction(false, Frame.Lx, Frame.Ly, Frame.gridSize, Frame.Params, Frame.Cameras[2], Frame.Images[2], x, Frame.centerPoint, Frame.minThreshold, Frame.maxThreshold, Frame.min_u, Frame.max_u, Frame.min_v, Frame.max_v);

		for(size_t t=E1.size()+E2.size();t<E1.size()+E2.size()+E3.size();t++)
			fi[t]=E3[t-E1.size()-E2.size()];}

}
int main ( )
{ // Main function to reconstruct time-dependent water surface.
	//createFeatureFile();
	//return 0;

	//Read all settings
	Settings set;
	FileStorage fss("settings.xml", FileStorage::READ); // Read the settings
	if (!fss.isOpened())
		{
			cout << "Could not open the configuration file: "<< endl;
			exit(1);
		 }
	fss["Settings"] >> set;
	fss.release();   // close Settings file
	CameraParams cam1, cam2, cam3;

	// Reconstructed area currently same for every camera used
	int min_u=set.min_u;
	int max_u=set.max_u;
	int min_v=set.min_v;
	int max_v=set.max_v;


	switch(set.NumberOfCameras){ // Choice of cameras and initialize each camera separately

		case 1:
			cam1=Initialize(set, 1);
			break;
		case 2:
			cam1=Initialize(set, 1);
			cam2=Initialize(set, 2);
			break;
		case 3:
			cam1=Initialize(set, 1);
			cam2=Initialize(set, 2);
			cam3=Initialize(set, 3);
			break;
		default:
			cout<<"Invalid amount of camera's"<<endl;
			return 0;
		}
	vector<imageOptimizer> optimizers(set.ThreadAmount);

	for(size_t i=0;i<set.ThreadAmount;i++){
	vector<CameraParams> cams;
	switch(set.NumberOfCameras){ //Choice of cameras and initialize optimization with appropriate camera amount

		case 1:
			cams.push_back(cam1);
			optimizers[i]=imageOptimizer(1, cams, set.SurfaceModelParameters, set.Lx, set.Ly, set.minThreshold, set.maxThreshold, set.gridSize, set.Scaling, set.Epsg, set.Epsf, set.Epsx, set.MaxIts, set.DiffStep, min_u, max_u, min_v, max_v);
			break;
		case 2:
			cams.push_back(cam1);
			cams.push_back(cam2);
			optimizers[i]=imageOptimizer(2, cams, set.SurfaceModelParameters, set.Lx, set.Ly, set.minThreshold, set.maxThreshold, set.gridSize, set.Scaling, set.Epsg, set.Epsf, set.Epsx, set.MaxIts, set.DiffStep, min_u, max_u, min_v, max_v);
			break;
		case 3:
			cams.push_back(cam1);
			cams.push_back(cam2);
			cams.push_back(cam3);
			optimizers[i]=imageOptimizer(3, cams, set.SurfaceModelParameters, set.Lx, set.Ly, set.minThreshold, set.maxThreshold, set.gridSize, set.Scaling, set.Epsg, set.Epsf, set.Epsx, set.MaxIts, set.DiffStep, min_u, max_u, min_v, max_v);
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

	switch(set.NumberOfCameras){ // Choice of cameras and read image-names in every directory (for each camera)

		case 1:
			InputDirectory1=set.InputDirectory1;
			check1=readStringList(imageList1, InputDirectory1); //read all image-names
			check=check1;
			break;
		case 2:
			InputDirectory1=set.InputDirectory1;
			InputDirectory2=set.InputDirectory2;
			check1=readStringList(imageList1, InputDirectory1); //read all image-names
			check2=readStringList(imageList2, InputDirectory2); //read all image-names
			check=check1 && check2;
			break;
		case 3:
			InputDirectory1=set.InputDirectory1;
			InputDirectory2=set.InputDirectory2;
			InputDirectory3=set.InputDirectory3;
			check1=readStringList(imageList1, InputDirectory1); //read all image-names
			check2=readStringList(imageList2, InputDirectory2); //read all image-names
			check3=readStringList(imageList3, InputDirectory3); //read all image-names
			check=check1 && check2 && check3;
			break;
		}

		// Check if image-names are successfully loaded
		if(check){
			cout<<"image-names loaded"<<endl;
			nimages=imageList1.size(); //number of images
		}
		else{
			cout<<"Images could not be loaded"<<endl;
			exit(1);
		}
	omp_set_dynamic(0);
	omp_set_num_threads(set.ThreadAmount);
	vector<vector<real_1d_array> >CoefficientsList(set.ThreadAmount); // Output list of coefficients
	vector<real_1d_array> prevCoefficientsList(set.ThreadAmount);


	size_t images_thread=nimages/set.ThreadAmount;
	double prevprev;

	#pragma omp parallel
	{
		#pragma omp for nowait
		for(size_t j=0;j<set.ThreadAmount;j++){


			for(size_t i=j*images_thread;i<(j+1)*images_thread;i++){ // Loop over all images
				vector<Mat> Images;
				Mat image, image2, image3;
				switch(set.NumberOfCameras){

						case 1:
							image=imread(imageList1[i],0);
							undistort(image, image2, cam1.K, cam1.distCoeffs, noArray() );
							Images.push_back(image2);
							break;
						case 2:
							image=imread(imageList1[i],0);
							undistort(image, image2, cam1.K, cam1.distCoeffs, noArray() );
							Images.push_back(image2);
							image=imread(imageList2[i],0);
							image2=Mat();
							undistort(image, image2, cam2.K, cam2.distCoeffs, noArray() );
							Images.push_back(image2);
							break;
						case 3:
							image=imread(imageList1[i],0);
							undistort(image, image2, cam1.K, cam1.distCoeffs, noArray() );
							Images.push_back(image2);
							image=imread(imageList2[i],0);
							image2=Mat();
							undistort(image, image2, cam2.K, cam2.distCoeffs, noArray() );
							Images.push_back(image2);
							image=imread(imageList3[i],0);
							image3=Mat();
							undistort(image, image3, cam3.K, cam3.distCoeffs, noArray() );
							Images.push_back(image3);
							break;
						default:
							cout<<"Error camera amount"<<endl;
							exit(1);
				}
				optimizers[j].changeImages(Images);
				minlmstate state;
				minlmreport rep;
				void *ptr;
				ptr=&optimizers[j];

			if(i-j*images_thread>0.1){

				vector<uchar> E1=errorfunction(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[0], optimizers[j].Images[0], prevCoefficientsList[j], optimizers[j].centerPoint, set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);
				vector<uchar> E2;
				if(set.NumberOfCameras>1)
					E2=errorfunction(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[1], optimizers[j].Images[1], prevCoefficientsList[j], Point2f(0,0), set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);

				double temp=prevCoefficientsList[j][0];
				if(i>1){
					cout<<prevprev<<endl;
					temp=prevCoefficientsList[j][0]-prevprev;
								}
				prevprev=prevCoefficientsList[j][0];
				if(temp>0)
					prevCoefficientsList[j][0]=prevCoefficientsList[j][0]+1;
				else prevCoefficientsList[j][0]=prevCoefficientsList[j][0]-1;
				cout<<prevCoefficientsList[j][0]<<endl;


				//bndl[0]=prevCoefficientsList[j][0]-5;
				//bndu[0]=prevCoefficientsList[j][0]+5;

				minlmcreatev(E1.size() +E2.size(), prevCoefficientsList[j], optimizers[j].DiffStep, state);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				minlmcreatev(E1.size() +E2.size(), prevCoefficientsList[j], optimizers[j].DiffStep/2, state);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				minlmcreatev(E1.size() +E2.size(), prevCoefficientsList[j], optimizers[j].DiffStep/10, state);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				// Repeat refining differentiation step until difference in water depth during refinement is less than hard-coded limit (0.05 mm)
				bool repeat=true;
				double num=20;
				while(repeat){
					double temporary=prevCoefficientsList[j][0];
					minlmcreatev(E1.size() +E2.size(), prevCoefficientsList[j], optimizers[j].DiffStep/num, state);
					minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
					alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
					minlmresults(state, prevCoefficientsList[j], rep);
					num=num+20;
					cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;
					if(abs(prevCoefficientsList[j][0]-temporary)<0.05)repeat=false;
				}
			CoefficientsList[j].push_back(prevCoefficientsList[j]);
			}
			else{

				/* Example to compute middle of approximate reconstructed surface area
				   Used in paper during validation on still water

				prevCoefficientsList[j]="[118.1555]";
				real_1d_array bndl="[85]";
				real_1d_array bndu = "[90]";
				vector<uchar> E=errorfunction(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[0], optimizers[j].Images[0], prevCoefficientsList[j], Point2f(0,0), set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);

				minlmcreatev(E.size() , prevCoefficientsList[j], optimizers[j].DiffStep, state);
				//alglib::minlmsetbc(state, bndl, bndu);
				//minlmsetscale(state, scale);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				//cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				minlmcreatev(E.size(), prevCoefficientsList[j], optimizers[j].DiffStep/2, state);
				//alglib::minlmsetbc(state, bndl, bndu);
				//minlmsetscale(state, scale);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				//cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				minlmcreatev(E.size(), prevCoefficientsList[j], optimizers[j].DiffStep/4, state);
				//alglib::minlmsetbc(state, bndl, bndu);
				//minlmsetscale(state, scale);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				//cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				minlmcreatev(E.size(), prevCoefficientsList[j], optimizers[j].DiffStep/8, state);
				//alglib::minlmsetbc(state, bndl, bndu);
				//minlmsetscale(state, scale);
				minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
				alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
				minlmresults(state, prevCoefficientsList[j], rep);
				cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

				bndl="[85,-3,-3]";

				bndu = "[75,3,3]";
				Point2f M=computeMiddle(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[0], optimizers[j].Images[0], prevCoefficientsList[j], set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);

				optimizers[j].centerPoint=M;
				ptr=&optimizers[j];
				cout<<"Approximate centre reconstructed area: "<<M<<endl;
			*/

			prevCoefficientsList[j]=set.InitialGuess.c_str();
			vector<uchar> E1=errorfunction(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[0], optimizers[j].Images[0], prevCoefficientsList[j], Point2f(0,0), set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);
			vector<uchar> E2;
			if(set.NumberOfCameras>1)
				E2=errorfunction(false, set.Lx, set.Ly, optimizers[j].gridSize, optimizers[j].Params, optimizers[j].Cameras[1], optimizers[j].Images[1], prevCoefficientsList[j], Point2f(0,0), set.minThreshold, set.maxThreshold, min_u, max_u, min_v, max_v);

			minlmcreatev(E1.size() +E2.size() , prevCoefficientsList[j], optimizers[j].DiffStep, state);
			minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
			alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
			minlmresults(state, prevCoefficientsList[j], rep);
			cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

			minlmcreatev(E1.size() +E2.size() , prevCoefficientsList[j], optimizers[j].DiffStep/5, state);
			minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
			alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
			minlmresults(state, prevCoefficientsList[j], rep);
			cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

			minlmcreatev(E1.size() +E2.size() , prevCoefficientsList[j], optimizers[j].DiffStep/15, state);
			minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
			alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
			minlmresults(state, prevCoefficientsList[j], rep);
			cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;

			minlmcreatev(E1.size() +E2.size() , prevCoefficientsList[j], optimizers[j].DiffStep/50, state);
			minlmsetcond(state, optimizers[j].Epsg, optimizers[j].Epsf, optimizers[j].Epsx, optimizers[j].MaxIts);
			alglib::minlmoptimize(state, combinedErrorFunction, NULL, ptr); //optimize
			minlmresults(state, prevCoefficientsList[j], rep);
			cout<<imageList1[i]<<"\t"<<prevCoefficientsList[j].tostring(5).c_str()<<endl;
			CoefficientsList[j].push_back(prevCoefficientsList[j]);

			}
		}
		}

	}

	writeArray(CoefficientsList, set.OutputFileName, "");
	timestamp ( );
    cout << "\n";
    cout << "Surface reconstruction \n";
    cout << "C++ version:\n";
    cout << "By Lukas Engelen:\n";
    return 0;
}


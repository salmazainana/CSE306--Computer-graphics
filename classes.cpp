#define _CRT_SECURE_NO_WARNINGS 1

#include <algorithm>
#include <cmath>
#include <limits>

#include <random>
static std::default_random_engine engine(10) ;
static std::uniform_real_distribution<double> uniform (0 , 1) ;


/*------------------------------------------ VECTOR CLASS ------------------------------------------*/

// Implement path tracer TD1 

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
    Vector& operator+=(const Vector& V) {
        data[0] += V[0]; 
        data[1] += V[1]; 
        data[2] += V[2]; 
        return *this ;
    }
    
    double norm2() {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	double norm() {
		return sqrt(norm2());
	}
    void normalize() {
		double n = norm();
		data[0] /= n;
		data[1] /= n;
		data[2] /= n;
	}
    const double& operator []( int i) const { return data[i] ;}
    double& operator []( int i) { return data[i]; }
    double data[3] ;
};

Vector operator+(const Vector& V, const Vector& U) {
    return Vector(V[0] + U[0] , V[1] + U[1] , V[2] + U[2]) ;
}
Vector div(Vector b, Vector a)
{
    return Vector(b[0] / a[0], b[1] / a[1], b[2] / a[2]);
}
Vector operator-(const Vector& V, const Vector& U){
    return Vector(V[0] - U[0] , V[1] - U[1] , V[2] - U[2]) ;
}
Vector operator*(const double alpha, const Vector& V) {
	return Vector(alpha*V[0], alpha*V[1], alpha*V[2]);
}
Vector operator*(const Vector& V, const double alpha) {
	return Vector(alpha*V[0], alpha*V[1], alpha*V[2]);
}
Vector operator/(const Vector& V, const double alpha) {
	return V*(1/alpha);
}

Vector cross(const Vector& V, const Vector &U) {
    return Vector(V[1] * U[2] - V[2] * U[1], V[2] * U[0] - V[0] * U[2], V[0] * U[1] - V[1] * U[0]);
}

double dot( const Vector& V, const Vector& U) { 
    return V[0] * U[0] + V[1] * U[1] + V[2] * U[2];
}

Vector operator*(const Vector& V, const Vector& U) {
	return Vector(V[0]*U[0], V[1]*U[1], V[2]*U[2]);       
}

Vector pow( Vector& V,  double n) {
    V[0]= std::pow(V[0], n);
    V[1]= std::pow(V[1], n);
    V[2]= std::pow(V[2], n);
	return Vector(V[0], V[1], V[2]);       
}

// Vector normalize(const Vector& V){
//     Vector U = V / V.norm();
//     return U ;
// }
/*------------------------------------------ AUX Functions------------------------------------------*/

Vector random_cos(const Vector& N) {

    //Indirect lightning Auxillary function

    double r1 = uniform(engine) ;
    double r2 = uniform(engine) ;
    double x = sqrt(1- r1 )*cos( 2 * M_PI*r2 );
    double y = sqrt(1- r1 )*sin( 2 * M_PI*r2);
    double z = sqrt(r1) ;

    Vector T1; 
    if ( std::abs(N[0])<= std::abs(N[1]) && std::abs(N[0])<= std::abs(N[2])){
        T1 = Vector(0,N[2],(-1)*N[1]);
    }
    else {
        if (std::abs(N[1])<= std::abs(N[0]) && std::abs(N[1])<= std::abs(N[2])){
            T1= Vector(N[2],0,(-1)*N[0]);
		}
        else{
            T1= Vector((-1)*N[1],N[0],0 );
        } 
	}
    T1.normalize();
	Vector T2;
    T2 = cross(T1, N); 
    return x*T1 + y*T2 + z*N;
}

void boxMuller(double stdev, double &x, double &y)
{
    double r1 = uniform(engine);
    double r2 = uniform(engine);

    x = sqrt(-2 * log(r1)) * cos(2 * M_PI * r2) * stdev;
    y = sqrt(-2 * log(r1)) * sin(2 * M_PI * r2) * stdev;
}

/*------------------------------------------ RAY CLASS ---------------------------------*/

class Ray{
    public:
        Vector origin;
        Vector direction;

        explicit Ray(Vector O, Vector u) {
            origin = O;
            direction = u;
        }
};

/* ------------------------------------------ CAMERA CLASS ------------------------------------------ */
class Camera{
public : 
    Vector position;
    double fov;

    explicit Camera(Vector Q, double alpha){
        position = Q;
        fov = alpha;
    }

    Ray normalized_ray_direction(int i,int j, bool transformed){
        Vector ray_u(0,0,0);

        int W = 512;
        int H = 512;
        int x = j; 
        int y = H - i -1; 

        ray_u[0] = x - W/2. + 0.5;
        ray_u[1] = y -  H/2. + 0.5;
        ray_u[2] = - W/(2*tan(fov/2));
        

        if (transformed){ // Monte carlo transformation 
            double a;
            double b;
            boxMuller(0.5, a, b);
            ray_u[0] += a;
            ray_u[1] += b;
        }

        ray_u.normalize();
        
        return Ray(position,ray_u );

    }   
};

class Geometry;
/*------------------------------------------ Intersection CLASS ---------------------------------*/

class Intersection{
public: 
    bool is_intersection;
    Vector point;
    Vector normal; 
    double distance;
    Vector rho;
    Geometry* geometry;
    Intersection(){};
    Intersection(bool is_inter, double t , Vector P, Vector NP, Vector rho, Geometry* geometry){
        is_intersection = is_inter; 
        point = P; 
        normal = NP; 
        distance = t;
        rho = rho;
        this->geometry = geometry;
    }
};

/* ------------------------------------------ Geometry CLASS ------------------------------------------ */

class Geometry{
public:
    virtual Intersection intersect(Ray &r) = 0;
    Vector albedo;
    Vector details;
};

/* ------------------------------------------ SPHERE CLASS ------------------------------------------ */
class Sphere : public Geometry{
public:
    Vector center;
    double radius;

    explicit Sphere(Vector C, double R, Vector color, Vector info){
        center = C; // center of the sphere
        radius = R; // Radius
        albedo = color; //   [0,1]^3 !!
        details = info; // [mirror_effect, is_transparent, refrac_ind] 1 == True && 0 == False
    }

    Intersection intersect(Ray& r){
    
        Vector u = r.direction; 
        Vector O = r.origin;
        
        Vector C = center;
        double R = radius;
        double is_inter;

        Vector OC = O - C;
        Vector CO = C - O;

        double delta = dot(u, OC) * dot(u, OC) - OC.norm2() + R*R ;

        if (delta < 0){  // No solution case
            is_inter = 0.; 
            double t = 1e15;
            return Intersection(is_inter, t, Vector(0,0,0), Vector(0,0,0),albedo, nullptr);
        }
        else {  // One or two solutions case 
            double t1 = -dot(u, OC) - sqrt(delta);
            double t2 = -dot(u, OC) + sqrt(delta);

            if (t2<0){ // THE RAY DOESN'T INTERSECT THE SPHERE
                is_inter = 0.; 
                double t = 1e15;
                return Intersection(is_inter, t, Vector(0,0,0), Vector(0,0,0),albedo, nullptr);
            }
            else{
                double t;
                if (t1 < 0){
                    t = t2; 
                }
                else{
                    t = t1;
                }    
                is_inter = 1.;
                Vector P = O + u * t; 
                Vector PC= P-C;  
                // Vector NaP =  normalize(PC);  // Normal @ P 
                PC.normalize();
                Vector NaP = PC;
                return Intersection(is_inter, t, P, NaP, albedo, this);
            }
        }
    }
};




/*-----------------------------------------BoundingBox Class-----------------------------------------------*/

class Boundingbox{
public:
    Vector B_min;
    Vector B_max;
    Boundingbox(){};
    Boundingbox(Vector &b_min,Vector &b_max){
        B_min = b_min;
        B_max = b_max;
    }

    bool is_intersect(Ray &r, double &d_inter)
    {

        Vector tmin = div(B_min - r.origin ,r.direction);
        Vector tmax = div(B_max - r.origin ,r.direction);

        for (int i =0; i<3; i++){
            if (tmin[i]>tmax[i]){
                std::swap(tmin[i], tmax[i]);
            }
        }

        double max = std::max(tmin[0], std::max(tmin[1], tmin[2])); // le max des mins
        double min = std::min(tmax[0], std::min(tmax[1], tmax[2])); // le miin des maxs

        if (min < max ){ return false;}

        d_inter = max;
        return true;
    }
};

/*----------------------------------------- Node Class-----------------------------------------------*/

class node{
public:
    node(){};
    int starting_triangle;
    int ending_triangle;
    Boundingbox Bbox;
    node* child_left;
    node* child_right;
    bool is_leaf;
};
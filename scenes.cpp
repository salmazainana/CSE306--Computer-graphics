#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

/* ------------------------------------------ SCENE CLASS ------------------------------------------ */
class Scene{
public:
    Vector light; 
    //std::vector<Sphere> objects;
    std::vector<Geometry*> objects;

    explicit Scene( Vector star){
        light = star;
    }

 	void add_geometry(Geometry* S) {
		this->objects.push_back(S);
	}
    
    Intersection intersection(Ray R){
        Vector rho; 
        Intersection closest = objects[0]->intersect(R) ;
        for (int idx=1; idx<objects.size(); idx++){
            Intersection sphere_inter = objects[idx]->intersect(R);
            if (sphere_inter.is_intersection && sphere_inter.distance <= closest.distance){
                closest = sphere_inter;
            }
        }
        return closest;
    } 

    Vector get_color(Ray r, int ray_depth, bool indirect_light){
        // Inspired in majority from the professor code and formulas from lecture notes  
        double I = 2e10; 
        Intersection inter = intersection(r);

        if (ray_depth < 0) { return Vector(0.,0.,0.);}
            
        if (!inter.is_intersection){ return Vector(0.,0.,0.);}

        Vector info = inter.geometry->details; // Vector (mirror, transp,n_idx)
        //std::cout << "idx is" << inter.idx <<std::endl;
        // The sphere is a mirror
        if (info[0]==1.)
        {
            //std::cout << "mirror" << std::endl;
            Vector reflexion_origine = inter.point + 1e-3*inter.normal;
            Vector mirror_direction =  r.direction - 2 * inter.normal * dot(inter.normal, r.direction);
            Ray reflected_ray = Ray(reflexion_origine,mirror_direction);
            return get_color(reflected_ray, ray_depth - 1, indirect_light);
        }

        // The sphere is transparent 
        else if(info[1] == 1.)
        {          
            //std::cout << "transparent" << std::endl;
            double n1 = 1; // refraction_index of the scene 
            double n2 = info[2];  // refraction_index of the sphere_idx
            double n = n1 / n2;
            Vector Ntransp = inter.normal;
            if (dot(r.direction, Ntransp) > 0) {
                std::swap(n1, n2);
                n = 1/n;
                Ntransp = (-1)*Ntransp;
            }
            
            double k_0 = pow(n1-n2,2)/pow(n1+n2,2);
            double R = k_0 + (1-k_0) * pow(1. - abs(dot(inter.normal,r.direction)),5);
            double u = uniform(engine);
            // Under the assumption that scene's light is fresnel light. 
            if ( u < R ){
                Ray reflected_ray = Ray(inter.point + 0.001* inter.normal, r.direction - 2*dot(r.direction,inter.normal)*inter.normal) ;
                return get_color(reflected_ray , ray_depth - 1, indirect_light) ;
            }
            Vector tTangent, tNormal;
            double cos = dot(r.direction, Ntransp);
            tTangent = n * (r.direction - cos * Ntransp);
            double rad = 1 - n*n*(1 - cos*cos);
            if (rad <= 0) {  
                //reflexion
                Vector direction = r.direction - 2 * dot(r.direction, Ntransp)*Ntransp ;
                Ray reflected_ray = Ray(inter.point - 0.001 * inter.normal,direction);
                return get_color(reflected_ray, ray_depth - 1,indirect_light);
            }
            else {
                //Refraction 
                tNormal = Ntransp * (-1) * sqrt(rad) ;
                Ray refractiveRay( inter.point - 0.001 * Ntransp, tTangent + tNormal); 
                return get_color(refractiveRay, ray_depth - 1,indirect_light); 
            }
        }

        else 
        {   //Lambertian model 
            //DIRECT COMPONENT 
            Vector color= Vector(0,0,0); 
            double d2 = (light - inter.point).norm2();
            Vector lightdir = (light - inter.point); 
            lightdir.normalize();
            Ray shadowRay(inter.point + 0.001 * inter.normal, lightdir);
            bool in_shadow = false;
            Intersection inter_shadow = intersection(shadowRay);
            if (inter_shadow.is_intersection){
                if ((inter_shadow.distance * inter_shadow.distance) < d2) {
                    in_shadow = true;
                }
            }
            if (!in_shadow) {
                color += (I / (4 * M_PI * d2)) * (inter.geometry->albedo / M_PI) *
                         std::max(0.0, dot(inter.normal, lightdir));
            }

            //INDIRECT COMPONENT 
            if (indirect_light){
                Ray indirect_ray(inter.point+0.001*inter.normal, random_cos(inter.normal));
                color += inter.geometry->albedo * get_color(indirect_ray, ray_depth-1,indirect_light);
            }
            return color;
        }
    }
}; 
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <chrono>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "classes.cpp"
#include "scenes.cpp"
#include "TriangleMesh.cpp"

// Time taken for chat 1: 158367 milliseconds with bvh 
// Time taken: 172639 milliseconds ~ 172,639 seconds ~ 2,87731667minutes - chat3 without optional features but bvh on 
// Time taken: 3248245 milliseconds for chat 4
/* ------------------------------------------ MAIN FUNCTION ------------------------------------------ */

int main() {

    // Timing 
    auto start_time = std::chrono::high_resolution_clock::now();

    // Scene version
    const int version = 1;

    int W = 512;
    int H = 512;

    //Light :
    Vector L(-10,20,40);
    double I = 2e10;

    //Camera
    Vector Q(0,0,55);
    double alpha = 60.*M_PI/180.;
    Camera C(Q, alpha);
    
    // Scene 
    Scene scene = Scene(L);
    
    // Scene parameters
    int n_rays = 100;
    bool indirect_light = true; // image false
    int ray_depth = 5;
    bool antialiasing = true;
        // The walls 
    Sphere* blue_wall = new Sphere(Vector(0,-1000,0), 990., Vector(0,0,1), Vector(0,0,1)); //albedos are random
    Sphere* pink_wall = new Sphere( Vector(0,0,1000), 940., Vector(0,1,0), Vector(0,0,1));
    Sphere* red_wall = new Sphere(Vector(0,1000,0), 940., Vector(1,0,1), Vector(0,0,1));
    Sphere* green_wall = new Sphere( Vector(0,0,-1000), 940., Vector(1,1,1), Vector(0,0,1));

    scene.add_geometry(blue_wall);
    scene.add_geometry(pink_wall);
    scene.add_geometry(red_wall);
    scene.add_geometry(green_wall);

    if (version==1){

        // The spheres
        Sphere* S1=new Sphere(Vector(-22,0,0), 10., Vector(0.5,0.5,0.5), Vector(1.,0.,1.4));
        Sphere* S2=new Sphere(Vector(0,0,0), 10., Vector(1,0,0), Vector(0.,1.,1.5));
        Sphere* S3=new Sphere(Vector(22,0,0), 10., Vector(1,0,1), Vector(0.,0.,1.4));

        // add objects
        scene.add_geometry(S1);
        scene.add_geometry(S2);
        scene.add_geometry(S3);

    }
    if (version==2)
    {   
        // Parameters 
        bool BVH = true;
        Vector color = Vector(0.3, 0.2, 0.25);
        double scaling = 0.6;   
        Vector translation = Vector(0,-10,0);
        Vector details = Vector(1., 1.,0.3);

        TriangleMesh *Meshes = new TriangleMesh(color, scaling, translation, details); 
        Meshes->readOBJ("cat.obj");

        for (int i = 0; i < Meshes->vertices.size(); i++)
        {
            Meshes->vertices[i] = Meshes->vertices[i] * scaling;
            Meshes->vertices[i] = Meshes->vertices[i] + translation;
        }
        if (BVH)
            Meshes->BVH(Meshes->root, 0, Meshes->indices.size());
        else 
        {
            Meshes->root->Bbox = Meshes->compute_bbox(0, Meshes->indices.size());
            Meshes->root->starting_triangle = 0;
            Meshes->root->ending_triangle = Meshes->indices.size();
            Meshes->root->is_leaf = true;
        }
        scene.add_geometry(Meshes);
    }

    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i= 0; i < H; i++){
        for (int j = 0; j < W; j++) {
            //rays from camera
            Vector colors = Vector(0.,0.,0.);
            for (int k = 0; k < n_rays;k++){
                Ray r = C.normalized_ray_direction(i,j,antialiasing);
                colors += scene.get_color(r,ray_depth,indirect_light);
            }
            colors = colors/n_rays;
            colors = pow(colors, 1/(2.2));

            image[(i * W + j) * 3] = std::min(255., colors[0]);
            image[(i * W + j) * 3 + 1] = std::min(255., colors[1]);
            image[(i * W + j) * 3 + 2] = std::min(255., colors[2]);
        }
    } 

    stbi_write_png("lab3.png", W, H, 3, &image[0], 0);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}
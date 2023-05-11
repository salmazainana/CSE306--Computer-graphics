#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>

//Copied from the link given in lecture notes
class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
	};
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

class TriangleMesh: public Geometry {
public:
	double scaling;
	Vector translation;
	node *root;
	Vector details; // [mirror_effect, is_transparent, refrac_ind] 1 == True && 0 == False

    ~TriangleMesh() {};

	TriangleMesh(Vector color, double factor, Vector trans, Vector details){
		albedo = color;
		scaling = factor;
		translation = trans;
		root = new node();
		details = details;
	}
		
	Boundingbox compute_bbox(int starting_triangle, int ending_triangle)
	{
		Boundingbox Bbox;
		double max_double = std::numeric_limits<double>::max();
		Vector B_max = (-1)*Vector(max_double,max_double,max_double) ;
		Vector B_min = Vector(max_double,max_double,max_double);

		for (int i = starting_triangle; i < ending_triangle; i++) // the ending triangle excluded 
		{
			TriangleIndices triangle_idx = indices[i];
			//std::list<Vector> dots;
			for (int j : {triangle_idx.vtxi,triangle_idx.vtxj,triangle_idx.vtxk}) {
				Vector point = vertices[j];

				B_max[0] = std::max(B_max[0], point[0]);
				B_min[0]= std::min(B_min[0], point[0]);
				B_max[1] = std::max(B_max[1], point[1]);
				B_min[1] = std::min(B_min[1], point[1]);
				B_max[2] = std::max(B_max[2], point[2]);
				B_min[2] = std::min(B_min[2], point[2]);
			}
		}
		Bbox.B_min = B_min;
		Bbox.B_max = B_max;
		return Bbox;
	}

	// Lecture page 49 BVH algorithm
	void BVH(node *node, int starting_triangle, int ending_triangle)
    {
		node->Bbox = compute_bbox(starting_triangle, ending_triangle ); 
		node->starting_triangle = starting_triangle ;
		node->ending_triangle = ending_triangle ;
		Boundingbox bbox = node->Bbox;
		Vector diag = bbox.B_max - bbox.B_min;
		Vector middle_diag = bbox.B_min + diag * 0.5 ;
		int longest_axis = diag[0] > diag[1] && diag[0] > diag[2] ? 0 : (diag[1] > diag[2] ? 1 : 2);
		int pivot_index = starting_triangle ;
		for ( int i=starting_triangle ; i<ending_triangle ; i++) 
		{
			TriangleIndices triangle = indices[i];
    		Vector barycenter = (1./3) * Vector(vertices[triangle.vtxi] + vertices[triangle.vtxj] + vertices[triangle.vtxk]);
			if ( barycenter[longest_axis] < middle_diag[longest_axis ] ) {
				std::swap( indices[i] , indices[pivot_index] ) ;
				pivot_index++;
			}
		}
		if ( pivot_index<=starting_triangle || pivot_index>=(ending_triangle-1) || (ending_triangle-starting_triangle)<5 ) {
            node->child_left = nullptr;
            node->child_right = nullptr;
			node->is_leaf = true;
			return ;
        }
		node->is_leaf = false;
        node->child_left = new class node();
        node->child_right = new class node();
		BVH( node->child_left,starting_triangle,pivot_index);
		BVH( node->child_right,pivot_index,ending_triangle);
	}

	// Ray-BVH intersection from page 49 lecture notes 		
	Intersection intersect(Ray &r)
	{
		Intersection intersection;
		intersection.is_intersection = false;
		double d_inter;
		Boundingbox bbox = root->Bbox;
		if (!bbox.is_intersect(r, d_inter)){
			return intersection;
		}

		std::list<node*> nodes_to_visit;
		nodes_to_visit.push_front(root);
		double d_max = std::numeric_limits<double>::max();
		Vector v_i, v_j, v_k, cp;
		Vector diff1; 
		Vector diff2;
		while (!nodes_to_visit.empty())
		{
			node *curNode = nodes_to_visit.back();
			nodes_to_visit.pop_back();
			if (!curNode->is_leaf)
			{	
				if (curNode->child_left->Bbox.is_intersect(r, d_inter)){
					if (d_inter < d_max){
						nodes_to_visit.push_back(curNode->child_left);
					}
				}
				if (curNode->child_right->Bbox.is_intersect(r, d_inter)){
					if (d_inter < d_max){
						nodes_to_visit.push_back(curNode->child_right);
					}
				}
			}
			else
			{
				for (int i = curNode->starting_triangle; i < curNode->ending_triangle; i++)
				{
					TriangleIndices triangle = indices[i];
					v_i = vertices[triangle.vtxi];
					v_j= vertices[triangle.vtxj];
					v_k = vertices[triangle.vtxk];
					diff1 = v_j - v_i;
					diff2 = v_k - v_i;
					cp = cross(diff1, diff2);
					d_inter = dot(cp, v_i - r.origin) / dot(cp, r.direction);
					if (0 < d_inter && d_inter < d_max){
						double temp =  dot(cp, r.direction);
						double beta = (1/ temp) * dot(diff2, cross(v_i - r.origin, r.direction));
						double gamma = (-1/ temp) * dot(diff1, cross(v_i - r.origin, r.direction));
						double alpha = 1 - beta - gamma;
						if (beta > 0 && gamma > 0 && alpha > 0)
						{
							d_max = d_inter;
							intersection.distance = d_max;
							intersection.normal = cp;
							intersection.normal.normalize();
							intersection.point = v_i + diff1*beta + diff2*gamma;
                            intersection.geometry = this;
							intersection.is_intersection = true;
						}
					}
				}
			}
		}
		return intersection;
	}

	void readOBJ(const char* obj) {

		char matfile[255];
		char grp[255];

		FILE* f;
		f = fopen(obj, "r");
		int curGroup = -1;
		while (!feof(f)) {
			char line[255];
			if (!fgets(line, 255, f)) break;

			std::string linetrim(line);
			linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
			strcpy(line, linetrim.c_str());

			if (line[0] == 'u' && line[1] == 's') {
				sscanf(line, "usemtl %[^\n]\n", grp);
				curGroup++;
			}

			if (line[0] == 'v' && line[1] == ' ') {
				Vector vec;

				Vector col;
				if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
					col[0] = std::min(1., std::max(0., col[0]));
					col[1] = std::min(1., std::max(0., col[1]));
					col[2] = std::min(1., std::max(0., col[2]));

					vertices.push_back(vec);
					vertexcolors.push_back(col);

				} else {
					sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
					vertices.push_back(vec);
				}
			}
			if (line[0] == 'v' && line[1] == 'n') {
				Vector vec;
				sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				normals.push_back(vec);
			}
			if (line[0] == 'v' && line[1] == 't') {
				Vector vec;
				sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
				uvs.push_back(vec);
			}
			if (line[0] == 'f') {
				TriangleIndices t;
				int i0, i1, i2, i3;
				int j0, j1, j2, j3;
				int k0, k1, k2, k3;
				int nn;
				t.group = curGroup;

				char* consumedline = line + 1;
				int offset;

				nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
				if (nn == 9) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
					if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
					if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
					if (nn == 6) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
						if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
						if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
						if (nn == 3) {
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							indices.push_back(t);
						} else {
							nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
							if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
							if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
							if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
							if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
							if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
							if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
							indices.push_back(t);
						}
					}
				}

				consumedline = consumedline + offset;

				while (true) {
					if (consumedline[0] == '\n') break;
					if (consumedline[0] == '\0') break;
					nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
					TriangleIndices t2;
					t2.group = curGroup;
					if (nn == 3) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
						if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
						if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
						indices.push_back(t2);
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						k2 = k3;
					} else {
						nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
							if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
							if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
							consumedline = consumedline + offset;
							i2 = i3;
							j2 = j3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
							if (nn == 2) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
								if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
								if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
								consumedline = consumedline + offset;
								i2 = i3;
								k2 = k3;
								indices.push_back(t2);
							} else {
								nn = sscanf(consumedline, "%u%n", &i3, &offset);
								if (nn == 1) {
									if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
									if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
									if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
									consumedline = consumedline + offset;
									i2 = i3;
									indices.push_back(t2);
								} else {
									consumedline = consumedline + 1;
								}
							}
						}
					}
				}
			}
		}
		fclose(f);
	}

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};

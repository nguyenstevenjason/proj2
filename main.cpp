#include <stack>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include "Transform.h"
#include <vector>
//#include <GL/glew.h>
//#include <GL/glut.h>

using namespace std ; 

float probMirror = 34.0f/100.0f;
float probRand = 34.0f/100.0f;
float probDirect = 31.0f/100.0f;
float probEmit = 1.0f/100.0f;

vec3 center; 
vec3 eye; // The (regularly updated) vector coordinates of the eye location 
vec3 up;  // The (regularly updated) vector coordinates of the up location 

const float ip_corners[4][3] = { {-1, 1, -1}, {1, 1, -1}, {1, -1, -1}, {-1, -1, -1}};

int width, height;  // assigned from size later
const int maxwidth = 500;
const int maxheight = 500;
int depth = 5; // maximum recursive depth
float depth_of_field = 1;
//enum shape {sphere, triangle} ; // geometry
int size[2] ;   // size of the image
//vec3 lights [];
int numlights = 0;
float attenuation[3] = {1, 0, 0};

string filename;     // if a filename is specified, this will be it
float camera[10];  // 0-2: lookfrom, 3-5: lookat, 6-8: up 9: fov

vec3 vertarray[1000]; // array storing every single vertex, numbered 0 through 999
int facearray[2000][3]; // array storing faces (triangles, up to 2000), which consist of 3 verticies each
int maxverts; // just a convenience variable
int maxvertnorms; // max num of verticies w/ normals
float vertexnormalarray[1000][6]; // x, y, z, nx, ny, nz
int vertexcounter = 0;
int trianglecounter = 0;
const int maxobjects = 50;

vec3 g_ambient;
vec3 g_diffuse;
vec3 g_specular;
vec3 g_emission; 
float g_shininess ;
float refraction;

vec3 current_translate;
vec4 current_rotate;
vec3 current_scale;

const int directional = 0;
const int point = 1;

float shadowoffset = 0.2;


struct Color {
	int r, g, b;
};

class Light {
public:
	int type; // 0 directional, 1 point
	vec3 pos; // Either the position for a point light, or direction for a directional
	//Color color;
	vec3 intensity;
	// attenuation too!
        float radius;
}lights[50];
int lightcounter = 0;

class Ray {
	public:
		vec4 pos;
		vec4 dir;
		bool is_shadowray;
		float t_min, t_max;   // these might be necessary for intersection calculations
		void RayThruPixel(int,int); // and a camera too
		float current_index; // initially 1
                void RandRayThruPixel(int,int);
};

void Ray::RayThruPixel(int i, int j) {
	// i,j are position of pixel, i -> height, w -> width
	// for vec4's, w term: 0 indicates vector (direction), 1 indicates point
	is_shadowray = 0;
	current_index = 1;

	vec3 a = eye - center;
	vec3 b = up;
	vec3 w = glm::normalize(a);
	vec3 u = glm::normalize(glm::cross(b, w)); 
	vec3 v = glm::cross(w, u);

	// Handling proper aspect ratio
	float aspect = ((float) width / (float) height);
	float fovy = camera[9];
	float fovx = fovy;
	// Lecture implementation, radians conversion
	fovy = fovy *(pi/180);
	fovx = fovx *(pi/180) ;

	//For depth of field, we randomly nudge i & j by a random amount, and redefine the focal plane

	float alpha = tan(aspect*fovx/2)*((float) (j - (width/2)) / (width/2));
	float beta = tan(fovy/2)*( (float) ( (height/2) - i) / (height/2)); 
	pos = vec4(eye,1);  // now they are set as vec4's
	//vec3 temp = glm::normalize(alpha*u + beta*v - w);

	if (depth_of_field >= 0) {
		vec3 rnd(( ((float) (rand() % 10)) / 50), ( ((float) (rand() % 10)) / 50), 0);
		vec3 neye = eye + rnd; // new eye in 3D
		vec4 neweye = pos + vec4(rnd, 0);  // new eye in 4D

		vec3 temp = alpha*u + beta*v - w;
		dir = glm::normalize(vec4(temp,0));
		float D = 1; // eye to image plane (always 1)
		float d = sqrt(pow(temp[0], 2) + pow(temp[1], 2) + pow(temp[2], 2));  // dist of eye to pixel
		// image plane is 1 away from eye 		
		
		vec4 p = pos + (d/ (D / (D + 4.5)) )*dir; // the intersection point of the ray with the focal point
		vec4 newdir = glm::normalize(p - neweye);
		pos = neweye;
		dir = newdir;
	}
	else {
		vec3 temp = alpha*u + beta*v - w;
		dir = glm::normalize(vec4(temp,0));
	}
}

void Ray::RandRayThruPixel(int i, int j) {
    	is_shadowray = 0;
	current_index = 1;
        
        float randPhi = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
        float randTheta = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
        float randPsi = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
        
        mat3 rotation(cos(randTheta)*cos(randPsi),
                              -cos(randTheta)*sin(randPsi),
                              sin(randTheta),
                              cos(randPhi)*sin(randPsi) + sin(randPhi)*sin(randTheta)*cos(randPsi),
                              cos(randPhi)*cos(randPsi) - sin(randPhi)*sin(randTheta)*sin(randPsi),
                              -sin(randPhi)*cos(randTheta),
                              sin(randPhi)*sin(randPsi) - cos(randPhi)*sin(randTheta)*cos(randPsi),
                              sin(randPhi)*cos(randPsi) + cos(randPhi)*sin(randTheta)*sin(randPsi),
                              cos(randPhi)*cos(randTheta)
                        );
        
	vec3 a = eye - center;
	vec3 b = up;
	vec3 w = glm::normalize(a);
	vec3 u = glm::normalize(glm::cross(b, w)); 
	vec3 v = glm::cross(w, u);

	// Handling proper aspect ratio
	float aspect = ((float) width / (float) height);
	float fovy = camera[9];
	float fovx = fovy;
	// Lecture implementation, radians conversion
	fovy = fovy *(pi/180);
	fovx = fovx *(pi/180) ;

	//For depth of field, we randomly nudge i & j by a random amount, and redefine the focal plane

	float alpha = tan(aspect*fovx/2)*((float) (j - (width/2)) / (width/2));
	float beta = tan(fovy/2)*( (float) ( (height/2) - i) / (height/2)); 
	pos = vec4(eye,1);  // now they are set as vec4's
	//vec3 temp = glm::normalize(alpha*u + beta*v - w);

	if (depth_of_field >= 0) {
		vec3 rnd(( ((float) (rand() % 10)) / 50), ( ((float) (rand() % 10)) / 50), 0);
		vec3 neye = eye + rnd; // new eye in 3D
		vec4 neweye = pos + vec4(rnd, 0);  // new eye in 4D

		vec3 temp = alpha*u + beta*v - w;
		dir = glm::normalize(vec4(temp,0));
		float D = 1; // eye to image plane (always 1)
		float d = sqrt(pow(temp[0], 2) + pow(temp[1], 2) + pow(temp[2], 2));  // dist of eye to pixel
		// image plane is 1 away from eye 		
		
		vec4 p = pos + (d/ (D / (D + 4.5)) )*dir; // the intersection point of the ray with the focal point
		vec4 newdir = glm::normalize(p - neweye);
		pos = neweye;
		dir = newdir;
	}
	else {
		vec3 temp = alpha*u + beta*v - w;
		dir = glm::normalize(vec4(temp,0));
	}
        dir = vec4(glm::normalize(rotation * vec3(dir.x, dir.y, dir.z)), 0);
}
// Shape is a basic abstract class, holds intersect method
class Shape {
public:
	Color color;    // RGB color, from scale of 0 to 1
	float intersect(Ray); 
	vec3 ambient ;  // Lighting properties, as vec3's
	vec3 diffuse ; 
	vec3 specular ;
	vec3 emission ; 
	float shininess ; // simple float
	mat4 matrix ;  // Transformation matrix, mat4
	vec4 normal; // to get around messy transformation issues
	float refractive_index; // n, index of refraction
};

class Sphere: public Shape {
	
public:
	float args[4];  //float x, y, z, radius;
	float intersect(Ray);
}spheres[100];
int spherecounter = 0;    // Tracks how many spheres in scene

float Sphere::intersect(Ray r) {
	// actual implementation goes here
	float ret = 999999;
	vec4 pos = r.pos;
	vec4 dir = r.dir;
	//if (r.is_shadowray == true) {
	//	r.t_min = 0.001; // arbitrary?
	//}
	r.t_min = 0;
	// In the case of a miss, r.t_min is 0, and returned t is 999999 (6 9's)
	vec4 cen = vec4(args[0], args[1], args[2],1); //center of the sphere
	float A = glm::dot(dir,dir);
	float B = 2*glm::dot(dir, pos - cen);
	float C = glm::dot(pos-cen, pos-cen) - args[3]*args[3];
	float det = B*B - 4*A*C;
	if (det >= 0) {
		float root1 = (-B + sqrt(det))/2*A;
		float root2 = (-B - sqrt(det))/2*A;
		float min;
		if (root1 < root2) {
			min = root1;
		} else {
			min = root2;
		}
		if (det = 0) {
			ret = root1;
		} else if (root1 >= 0 && root2 >= 0) {
			ret = min;
		} else if (root1 >=0 && root2 < 0) {
			ret = root1;
		} else if (root1 < 0 && root2 >= 0) {
			ret = root2;
		}
	}
	//if (ret >= r.t_min) {
		return ret;
	//}
	//else {
	//	return 999999;
	//}
}

class Triangle: public Shape {
public:
	vec3 vertexarray[3];  // 3 verticies of a triangle
	float intersect(Ray);
}triangles[500];

float Triangle::intersect(Ray r) {
	// actual implementation goes here
	float ret = 999999;
	vec4 pos = r.pos;
	vec4 dir = r.dir;
	r.t_min = 0;
	vec3 vert_a = vertexarray[0];
	vec3 vert_b = vertexarray[1];
	vec3 vert_c = vertexarray[2];
	normal = normal = glm::normalize(vec4(glm::normalize(glm::cross(vert_c - vert_a,vert_b - vert_a)),0));
	float a = vert_a[0] - vert_b[0]; // Xa - Xb
	float b = vert_a[1] - vert_b[1]; // Ya - Yb
	float c = vert_a[2] - vert_b[2]; // Za - Zb
	float d = vert_a[0] - vert_c[0];
	float e = vert_a[1] - vert_c[1];
	float f = vert_a[2] - vert_c[2];
	float g = dir[0];
	float h = dir[1];
	float i = dir[2];
	float j =  vert_a[0] - pos[0];
	float k =  vert_a[1] - pos[1];
	float l =  vert_a[2] - pos[2];
	float beta;
	float gamma;
	float t;
	float M = a*(e*i-h*f)+b*(g*f-d*i)+c*(d*h-e*g);
	beta = (j*(e*i-h*f)+k*(g*f-d*i)+l*(d*h-e*g))/M;
	gamma = (i*(a*k-j*b)+h*(j*c-a*l)+g*(b*l-k*c))/M;
	t = -(f*(a*k-j*b)+e*(j*c-a*l)+d*(b*l-k*c))/M;
	vec4 n = vec4(glm::normalize(glm::cross(vertexarray[2]-vertexarray[0],vertexarray[1]-vertexarray[0])),0);
	if (glm::dot(dir,n) != 0 && (t >= r.t_min) && (gamma >= 0) && (gamma <= 1) && (beta >= 0) && (beta <= 1 - gamma)) {
		ret = t;
	}
	return ret;
}

class Image {      // stores the color values and has method writeIamge?
public:
	Color ** colors;
	void writeImage();
	void initialize();
};

void Image::initialize() {  // magic memory allocation bullshit
	int size_x = height;
	int size_y = width;
	colors = (Color**) malloc (size_x * sizeof(Color *));
	for (int i = 0; i < size_x; i++) {
		colors[i] = (Color*) malloc(size_y * sizeof(Color));
	}
	//colors = (Color*) malloc (10000);
}

void Image::writeImage() {
	ofstream myfile; //("output.ppm"); //say we had an output file
	myfile.open ("output.ppm");
	if (myfile.is_open()) {
		myfile << "P3 \n";
		myfile << width << " " << height << " \n";
		myfile << "255 \n";
		//for (int i = height - 1; i >= 0; i--) {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {   // order is changed since ppm starts from top left corner
				myfile << colors[i][j].r << " " << colors[i][j].g << " " << colors[i][j].b << " ";
				//myfile << 0 << " " << 0 << " " << 0 << " ";
			}
			myfile << "\n";
		}
		myfile.close();
	} 
}

// **************************** Raytracer stuff

class Intersection {
	public:
		Shape s;
		Ray r;
		vec4 normal;
		Intersection (Ray);
		bool flag;
		vec4 point; // The intersection point, for use in transformations
};

Intersection::Intersection(Ray ray) {
	flag = false;
	float min = 999999;
	//r.t_min = min;
	for (int i = 0; i < spherecounter; i++) {
		Sphere sphere = spheres[i];
		
		Ray sray;
		sray.pos = ray.pos*glm::transpose(glm::inverse(sphere.matrix));
		vec4 tempdir = ray.dir*glm::transpose(glm::inverse(sphere.matrix));
		sray.dir = glm::normalize(tempdir);
		float t = sphere.Sphere::intersect(sray); //we now have the transformed intersection "t"
		vec4 newpoint = sray.pos + t*sray.dir;
		float temp = t;   // temp holds the object-space t value
		t = t/sqrt(glm::dot(tempdir, tempdir)); // returns a proper t in world coordinates. It seems to make everything really slow...
		
		//float t = sphere.Sphere::intersect(ray);
		if (t < min && temp < 999999) { // the min has to be world coordinate min, so here the problem comes in
			min = t; // now are we returning the correct t?
			point = newpoint*glm::transpose(sphere.matrix); // intersection point p turns into Mp
			s = spheres[i];
			flag = true;
			r = sray;
			r.current_index = ray.current_index;
			r.dir = sray.dir*glm::transpose(sphere.matrix);
			r.t_min = min;
			normal = glm::normalize(glm::normalize( (sray.pos + sray.dir*temp) - vec4(sphere.args[0], sphere.args[1], sphere.args[2],1))*glm::inverse(sphere.matrix));   // Normals are hit with inverse transpose
			//normal = glm::normalize(((ray.pos + ray.dir*t) - vec4(sphere.args[0], sphere.args[1], sphere.args[2],1)));
		}
	}
	for (int i = 0; i < trianglecounter; i++) {	
		Triangle triangle = triangles[i];
		Ray tray;
		
		tray.pos = ray.pos*glm::transpose(glm::inverse(triangle.matrix));
		vec4 tempdir = ray.dir*glm::transpose(glm::inverse(triangle.matrix));
		tray.dir = glm::normalize(tempdir);
		float t = triangle.Triangle::intersect(tray);
		float temp = t;
		vec4 newpoint = tray.pos + t*tray.dir;
		t = t/sqrt(glm::dot(tempdir, tempdir));
		
		//float t = triangle.Triangle::intersect(ray);
		if (t < min && temp < 999999) {
			min = t;
			point = newpoint*glm::transpose(triangle.matrix); // intersection point p turns into Mp
			s = triangles[i];
			flag = true;
			//r = ray;
			r = tray;
			r.current_index = ray.current_index;
			r.dir = tray.dir*glm::transpose(triangle.matrix);
			r.t_min = min;
			normal = -glm::normalize(triangle.normal*glm::inverse(triangle.matrix));
		}
	}
	
}

struct Colorf {
	float r, g, b;
};

class Shading {
public: 
	Colorf shadecolor;
	void DirectShade(Intersection, int); // also must take in all lights later
        void InDirectRandShade(Intersection, int);
        void MirrorShade(Intersection, int);
        void EmitShade(Intersection, int);
        
};
void Shading::DirectShade(Intersection hit, int d) {
	if (d < depth && (hit.flag == true)) { //recursive depth 
                Shape shape = hit.s;
		vec3 veccolor = vec3(shape.ambient[0] + shape.emission[0], shape.ambient[1] + shape.emission[1], 
			shape.ambient[2] + shape.emission[2]);
                vector<Ray> shadowRays;
		for (int i = 0; i < lightcounter; i++) {
			Ray shadow; // generate a shadow ray
                        shadow.is_shadowray = true;
                        shadow.current_index = 1;
			//shadow.pos = (hit.r.pos + hit.r.dir*hit.r.t_min + hit.normal*0.001);
			shadow.pos = (hit.point + hit.normal*0.001);  // now, replace calculating point with hit.point

			vec4 lightposition;
			vec4 lightdirection;
			vec4 lightside;
			vec4 lightup;
                        float area;
			if (lights[i].type == 1) {  // point light
				vec4 lightpositiontest = vec4(lights[i].pos,1);
				vec4 lightdirectiontest = glm::normalize(lightpositiontest - shadow.pos);
				
				lightside = glm::normalize(vec4(glm::cross(vec3(lightdirectiontest), up),0)); //might not work for lights directly above or below a point
				lightup = glm::normalize(vec4(glm::cross(vec3(lightside),vec3(lightdirectiontest)),0)); //fixing lightup will fix this as well
				//lightposition = lightpositiontest + unifRand()*lightup*shadowoffset + unifRand()*lightside*shadowoffset;  // Randomized shadows
				lightposition = lightpositiontest;
				lightdirection = glm::normalize(lightposition - shadow.pos); 
                                area = 1.0f;

			} else if (lights[i].type == 0) { //directional
				lightposition = vec4(lights[i].pos,1); // the 3d position in space
				lightdirection = glm::normalize(vec4(vec3(lightposition),0));
                                area = 1.0f;
			} else if (lights[i].type == 2) { //area
                                vec4 lightpositiontest = vec4(lights[i].pos,1);
				vec4 lightdirectiontest = glm::normalize(lightpositiontest - shadow.pos);
				
				lightside = glm::normalize(vec4(glm::cross(vec3(lightdirectiontest), up),0)); //might not work for lights directly above or below a point
				lightup = glm::normalize(vec4(glm::cross(vec3(lightside),vec3(lightdirectiontest)),0)); //fixing lightup will fix this as well
				//lightposition = lightpositiontest + unifRand()*lightup*shadowoffset + unifRand()*lightside*shadowoffset;  // Randomized shadows
				lightposition = lightpositiontest;
				lightdirection = glm::normalize(lightposition - shadow.pos);
                                area = 4.0f*pi*pow(lights[i].radius,2.0f);
                        }
                        
			vec4 eyedirn = glm::normalize(vec4(eye,1) - shadow.pos);
			vec4 half = glm::normalize(lightdirection + eyedirn);
                        
                        shadow.dir = lightdirection;
                        Intersection intersect(shadow);
                        float visibility;
                        if (lights[i].type == 1 && intersect.flag == true) {
                            visibility = (float)(sqrt(glm::dot(lightposition-shadow.pos, lightposition-shadow.pos))-sqrt(glm::dot(lightposition-intersect.point, lightposition-intersect.point)))/sqrt(glm::dot(lightposition-shadow.pos, lightposition-shadow.pos));
                        } else if (lights[i].type == 0 && intersect.flag == true) {
                            visibility = (float)0.0f;
                        } else if (lights[i].type == 2) {
                            float numshadows = 10;
                            visibility = (float)((float)1.0f/(float)(numshadows+(float)1.0f))*(sqrt(glm::dot(lightposition-shadow.pos, lightposition-shadow.pos))-sqrt(glm::dot(lightposition-intersect.point, lightposition-intersect.point)))/sqrt(glm::dot(lightposition-shadow.pos, lightposition-shadow.pos));
                            for (int s = 0; s < numshadows; s++) {
                                float randPhi = (float) ((float) rand() / RAND_MAX) * 180.0f-90.0f;
                                float randTheta = (float) ((float) rand() / RAND_MAX) * 180.0f-90.0f;
                                float randPsi = (float) ((float) rand() / RAND_MAX) * 180.0f-90.0f;
        
                                mat3 rotation(cos(randTheta)*cos(randPsi),
                                        -cos(randTheta)*sin(randPsi),
                                        sin(randTheta),
                                        cos(randPhi)*sin(randPsi) + sin(randPhi)*sin(randTheta)*cos(randPsi),
                                        cos(randPhi)*cos(randPsi) - sin(randPhi)*sin(randTheta)*sin(randPsi),
                                        -sin(randPhi)*cos(randTheta),
                                        sin(randPhi)*sin(randPsi) - cos(randPhi)*sin(randTheta)*cos(randPsi),
                                        sin(randPhi)*cos(randPsi) + cos(randPhi)*sin(randTheta)*sin(randPsi),
                                        cos(randPhi)*cos(randTheta)
                                );
                
                                vec4 tempdir = -shadow.dir;
                                vec4 newdir = glm::normalize(vec4(rotation*vec3(tempdir[0],tempdir[1],tempdir[2]),0));
                                vec4 templightpos = lightposition + newdir*lights[i].radius;
                                Ray newshadow;
                                newshadow.is_shadowray = true;
                                newshadow.current_index = 1;
                                newshadow.pos = (hit.point + hit.normal*0.001);
                                newshadow.dir = glm::normalize(templightpos - newshadow.pos);
                                Intersection newintersect(newshadow);
                                if (newintersect.flag == true) {
                                        visibility += (float)((float)1.0f/(float)(numshadows+(float)1.0f))*(float)(sqrt(glm::dot(templightpos-newshadow.pos, templightpos-newshadow.pos))-sqrt(glm::dot(templightpos-newintersect.point, templightpos-newintersect.point)))/sqrt(glm::dot(templightpos-newshadow.pos, templightpos-newshadow.pos));
                                } else {
                                        visibility += (float)((float)1.0f/(float)(numshadows+(float)1.0f))*(float)1.0f;
                                }
                            }
                        } else {
                            visibility = (float)1.0f;
                        }
                                
                        //visibility = 0.5f;
                        veccolor[0] = (veccolor[0] + visibility*area*lights[i].intensity[0]*(shape.diffuse[0]*max((float) glm::dot(hit.normal,lightdirection), (float) 0) + shape.specular[0]*pow((float) max((float) glm::dot(hit.normal,half), (float) 0), shape.shininess) ));
                        veccolor[1] = (veccolor[1] + visibility*area*lights[i].intensity[1]*(shape.diffuse[1]*max((float) glm::dot(hit.normal,lightdirection), (float) 0) + shape.specular[1]*pow((float) max((float) glm::dot(hit.normal,half), (float) 0), shape.shininess) ));
                        veccolor[2] = (veccolor[2] + visibility*area*lights[i].intensity[2]*(shape.diffuse[2]*max((float) glm::dot(hit.normal,lightdirection), (float) 0) + shape.specular[2]*pow((float) max((float) glm::dot(hit.normal,half), (float) 0), shape.shininess) ));
			
		}
		// Reflection rays, sent after all other rays are cast. Does this make sense?
		Ray rray;
                
		rray.dir = glm::normalize(hit.r.dir - (2*glm::dot(hit.r.dir,hit.normal))*hit.normal);
		rray.pos = hit.point + rray.dir*0.001;

		if (veccolor[0] > 1) {
				veccolor[0] = 1;
		}
		if (veccolor[0] < 0) {
				veccolor[0] = 0;
		}
		if (veccolor[1] > 1) {
				veccolor[1] = 1;
		}
		if (veccolor[1] < 0) {
			veccolor[1] = 0;
		}
		if (veccolor[2] > 1) {
				veccolor[2] = 1;
		}
		if (veccolor[2] < 0) {
				veccolor[2] = 0;
		}  
		shadecolor.r =  (1/(1-probDirect)) * veccolor[0];
		shadecolor.g =  (1/(1-probDirect)) * veccolor[1];
		shadecolor.b =  (1/(1-probDirect)) * veccolor[2];
	}


	else {
		shadecolor.r =  0;
		shadecolor.g =  0;
		shadecolor.b =  0;
	}
}

void Shading::InDirectRandShade(Intersection hit, int d) {
	if (d < depth && (hit.flag == true)) { //recursive depth 
                Shape shape = hit.s;
		vec3 veccolor = vec3(shape.ambient[0] + shape.emission[0], shape.ambient[1] + shape.emission[1], 
			shape.ambient[2] + shape.emission[2]);
		// Reflection rays, sent after all other rays are cast. Does this make sense?
                Ray randRay;
                float randPhi = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
                float randTheta = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
                float randPsi = (float) ((float) rand() / RAND_MAX) * 180.0f - 90.0f;
        
                mat3 rotation(cos(randTheta)*cos(randPsi),
                              -cos(randTheta)*sin(randPsi),
                              sin(randTheta),
                              cos(randPhi)*sin(randPsi) + sin(randPhi)*sin(randTheta)*cos(randPsi),
                              cos(randPhi)*cos(randPsi) - sin(randPhi)*sin(randTheta)*sin(randPsi),
                              -sin(randPhi)*cos(randTheta),
                              sin(randPhi)*sin(randPsi) - cos(randPhi)*sin(randTheta)*cos(randPsi),
                              sin(randPhi)*cos(randPsi) + cos(randPhi)*sin(randTheta)*sin(randPsi),
                              cos(randPhi)*cos(randTheta)
                        );
                
                vec3 newRraydir = vec3(-hit.normal.x, -hit.normal.y, hit.normal.z);
                vec3 rotRraydir = glm::normalize(rotation * newRraydir);
                
                randRay.dir = glm::normalize(vec4(rotRraydir, 1));
                randRay.pos = hit.point + randRay.dir*0.001;

                Intersection randreflect(randRay);
                Shading randReflectShade;
                
                Shading thisShade;
                thisShade.DirectShade(hit, d+1);
                
                float aveg = (thisShade.shadecolor.r + thisShade.shadecolor.g + thisShade.shadecolor.b)/3.0f;
                
                float randChoice = ((float) rand() / (float) RAND_MAX);
                
                if (aveg < randChoice) {
                    float randN = 1.0f;
                    while (randN > (probDirect + probEmit)) {
                        randN = ((float) rand() / (float) RAND_MAX);
                    }
                    if (randN < probDirect) {
                        randReflectShade.DirectShade(randreflect, d+1);
                    } else {
                        randReflectShade.EmitShade(randreflect, d+1);
                    }
                    /* WE EMIT */
                } else {
                    float randN = 1.0f;
                    while (randN > (probRand + probMirror)) {
                        randN = ((float) rand() / (float) RAND_MAX);
                    }
                    if (randN < probRand) {
                        randReflectShade.InDirectRandShade(randreflect, d+1);
                    } else {
                        randReflectShade.MirrorShade(randreflect, d+1);
                    }
                }
                /*float randN = ((float) rand() / (float) RAND_MAX);
                if (randN < probDirect) {
                    randReflectShade.DirectShade(randreflect,d+1);
                }
                if ((randN >= probDirect) && (randN < (probDirect+probRand))) {
                    randReflectShade.InDirectRandShade(randreflect, d+1);
                }
                if ((randN >= (probDirect+probRand)) && (randN < (probDirect+probRand+probMirror))) {
                    randReflectShade.MirrorShade(randreflect, d+1);
                }
                if ((randN >= (probDirect+probRand+probMirror)) && (randN < (probDirect+probRand+probMirror+probEmit))) {
                    randReflectShade.EmitShade(randreflect, d+1);
                    //randReflectShade.DirectShade(randreflect, d+1);
                }*/
                //vec3 cosineWeight = 1.0f;
                float dotVec = glm::dot(randRay.dir, hit.normal);
                float randRayMag = sqrt(glm::dot(randRay.dir, randRay.dir));
                float hitMag = sqrt(glm::dot(hit.normal, hit.normal));
                float cosTerm = dotVec / ((randRayMag) * (hitMag));
                cosTerm = 1.0f;
                //vec4 temp = glm::normalize(randRay.dir * hit.normal);
                veccolor[0] = shape.specular[0]* cosTerm *((float) randReflectShade.shadecolor.r);
                veccolor[1] = shape.specular[1]* cosTerm *((float) randReflectShade.shadecolor.g);
                veccolor[2] = shape.specular[2]* cosTerm *((float) randReflectShade.shadecolor.b);
                if (veccolor[0] > 1) {
				veccolor[0] = 1;
		}
		if (veccolor[0] < 0) {
				veccolor[0] = 0;
		}
		if (veccolor[1] > 1) {
				veccolor[1] = 1;
		}
		if (veccolor[1] < 0) {
			veccolor[1] = 0;
		}
		if (veccolor[2] > 1) {
				veccolor[2] = 1;
		}
		if (veccolor[2] < 0) {
				veccolor[2] = 0;
		}
		shadecolor.r =  (1/(1-probRand)) * veccolor[0];
		shadecolor.g =  (1/(1-probRand)) * veccolor[1];
		shadecolor.b =  (1/(1-probRand)) * veccolor[2];
	}


	else {
		shadecolor.r =  0;
		shadecolor.g =  0;
		shadecolor.b =  0;
	}
}
void Shading::MirrorShade(Intersection hit, int d) {
	if (d < depth && (hit.flag == true)) { //recursive depth 
                Shape shape = hit.s;
		vec3 veccolor = vec3(shape.ambient[0] + shape.emission[0], shape.ambient[1] + shape.emission[1], 
			shape.ambient[2] + shape.emission[2]);
		// Reflection rays, sent after all other rays are cast. Does this make sense?
		Ray rray;
               
		rray.dir = glm::normalize(hit.r.dir - (2*glm::dot(hit.r.dir,hit.normal))*hit.normal);
		rray.pos = hit.point + rray.dir*0.001;
                
                Intersection reflect(rray);
                Shading reflectshade;
                float randN = ((float) rand() / (float) RAND_MAX);
                if (randN < probDirect) {
                    reflectshade.DirectShade(reflect,d+1);
                }
                if ((randN >= probDirect) && (randN < (probDirect+probRand))) {
                    reflectshade.InDirectRandShade(reflect, d+1);
                }
                if ((randN >= (probDirect+probRand)) && (randN < (probDirect+probRand+probMirror))) {
                    reflectshade.MirrorShade(reflect, d+1);
                }
                if ((randN >= (probDirect+probRand+probMirror)) && (randN < (probDirect+probRand+probMirror+probEmit))) {
                    reflectshade.EmitShade(reflect, d+1);
                    //reflectshade.DirectShade(reflect, d+1);
                }
                //reflectshade.DirectShade(reflect, d+1);
                veccolor[0] = shape.specular[0]*((float) reflectshade.shadecolor.r);
                veccolor[1] = shape.specular[1]*((float) reflectshade.shadecolor.g);
                veccolor[2] = shape.specular[2]*((float) reflectshade.shadecolor.b);
                if (veccolor[0] > 1) {
        		veccolor[0] = 1;
		}
                
		if (veccolor[0] < 0) {
				veccolor[0] = 0;
		}
		if (veccolor[1] > 1) {
				veccolor[1] = 1;
		}
		if (veccolor[1] < 0) {
			veccolor[1] = 0;
		}
		if (veccolor[2] > 1) {
				veccolor[2] = 1;
		}
		if (veccolor[2] < 0) {
				veccolor[2] = 0;
		}
                
		shadecolor.r =  (1/(1-probMirror)) * veccolor[0];
		shadecolor.g =  (1/(1-probMirror)) * veccolor[1];
		shadecolor.b =  (1/(1-probMirror)) * veccolor[2];
	}


	else {
		shadecolor.r =  0;
		shadecolor.g =  0;
		shadecolor.b =  0;
	}
}
void Shading::EmitShade(Intersection hit, int d) {
	if (d < depth && (hit.flag == true)) { //recursive depth 
                Shape shape = hit.s;
		vec3 veccolor = vec3(shape.ambient[0] + shape.emission[0], shape.ambient[1] + shape.emission[1], 
			shape.ambient[2] + shape.emission[2]);
                
		shadecolor.r =  (1/(1-probEmit)) * veccolor[0];
		shadecolor.g =  (1/(1-probEmit)) * veccolor[1];
		shadecolor.b =  (1/(1-probEmit)) * veccolor[2];
	}


	else {
		shadecolor.r =  0;
		shadecolor.g =  0;
		shadecolor.b =  0;
	}
}

void Parser (const char * filename) {
	stack <mat4> transfstack ; 
	transfstack.push(mat4(1.0)) ; //sets initial value to identity
  string str, ret = "" ; 
  ifstream in ; 
  in.open(filename) ; 
  if (in.is_open()) {
    getline (in, str) ; 
	int n = 0;
    while (in) { 
		if ((str.find_first_not_of("\t\r\n") != string::npos) && (str[0] != '#')) {
			string cmd;
			stringstream s(str);
			s >> cmd;
				if (cmd == "directional" && lightcounter < 50) {   // change later
					lights[lightcounter].type = 0; // signifies directional 
					for (int i = 0; i < 6; i++) {  // xyz rgb
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							float x;
							s >> x;
							if (i < 3) {
								lights[lightcounter].pos[i] = x;
							} else {
								lights[lightcounter].intensity[i-3] = x;
							}
                                                        lights[lightcounter].radius = 0.001;
						}
					}
					lightcounter++;
				}
				else if (cmd == "point" && lightcounter < 50) {   // change later
					//Light l = lights[lightcounter];
					lights[lightcounter].type = 1; // signifies directional 
					for (int i = 0; i < 6; i++) {  // xyz rgb
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							float x;
							s >> x;
							if (i < 3) {
								lights[lightcounter].pos[i] = x;
							}
							else {
								lights[lightcounter].intensity[i-3] = x;
							}
                                                        lights[lightcounter].radius = 0.001;
						}
					}
					lightcounter++;
				}
                                else if (cmd == "area" && lightcounter < 50) {   // change later
					//Light l = lights[lightcounter];
					lights[lightcounter].type = 2; // signifies area 
					for (int i = 0; i < 7; i++) {  // xyz rgb
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							float x;
							s >> x;
							if (i < 3) {
								lights[lightcounter].pos[i] = x;
							}
							else if (i < 6) {
								lights[lightcounter].intensity[i-3] = x;
							}
                                                        else {
                                                            lights[lightcounter].radius = x;
                                                        }
						}
					}
					lightcounter++;
				}
				else if (cmd == "size") {   // refers to size of image, width and height
					for (int i = 0; i < 2; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> size[i];
						}
					}
				}
				else if (cmd == "attenuation") {   // refers to size of image, width and height
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> attenuation[i];
						}
					}
				}
				else if (cmd == "maxdepth") {   
					if (s.str().empty()) {
						cerr << "Not enough arguments to 'size'\n";
						throw 2;
					} else {
							s >> depth;
					}
				}
				else if (cmd == "camera") {
					for (int i = 0; i < 10; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> camera[i] ;
						}
					}
				}
				else if (cmd == "maxverts") {   
					if (s.str().empty()) {
						cerr << "Not enough arguments to 'size'\n";
						throw 2;
					} else {
							s >> maxverts;
					}
				}
				else if (cmd == "maxvertnorms") {   
					if (s.str().empty()) {
						cerr << "Not enough arguments to 'size'\n";
						throw 2;
					} else {
							s >> maxvertnorms;
					}
				}
				else if (cmd == "vertex") {
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> vertarray[vertexcounter][i] ;
						}
					}
					vertexcounter++;
				}
				else if (cmd == "tri") {
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							int x;
							s >> x;
							triangles[trianglecounter].vertexarray[i] = vertarray[x];  //conversion from 1 based indexing
						}
					}
					triangles[trianglecounter].ambient = g_ambient;
					triangles[trianglecounter].diffuse = g_diffuse;
					triangles[trianglecounter].emission = g_emission;
					triangles[trianglecounter].specular = g_specular;
					triangles[trianglecounter].shininess = g_shininess;
					triangles[trianglecounter].matrix = transfstack.top();
					triangles[trianglecounter].refractive_index = refraction;
					trianglecounter++;
					
				}
				else if (cmd == "sphere") {          // temporarily make sphere red
					for (int i = 0; i < 4; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> spheres[spherecounter].args[i] ;
						}
					}
					spheres[spherecounter].ambient = g_ambient;
					spheres[spherecounter].diffuse = g_diffuse;
					spheres[spherecounter].emission = g_emission;
					spheres[spherecounter].specular = g_specular;
					spheres[spherecounter].shininess = g_shininess;
					spheres[spherecounter].matrix = transfstack.top();
					spheres[spherecounter].refractive_index = refraction;
					spherecounter++;
				}
				//          Lighting
				else if (cmd == "ambient") {            // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> g_ambient[i];
						}
					}
				}
				else if (cmd == "diffuse") {          // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> g_diffuse[i] ;
						}
					}
				}
				else if (cmd == "specular") {            // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> g_specular[i] ;
						}
					}
				}
				else if (cmd == "emission") {               // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> g_emission[i] ;
						}
					}
				}
				else if (cmd == "shininess") {            // change later?
					if (s.str().empty()) {
						cerr << "Not enough arguments to 'size'\n";
						throw 2;
					} else {
						s >> g_shininess;
					}
				}
				else if (cmd == "refraction") {            // change later?
					if (s.str().empty()) {
						cerr << "Not enough arguments to 'size'\n";
						throw 2;
					} else {
						s >> refraction;
					}
				}
				else if (cmd == "pushTransform") {
					//flag = true ;
					mat4 top = transfstack.top() ;
					transfstack.push(top) ;
				}
				else if (cmd == "popTransform") {
					//flag = false ;
					
					transfstack.pop() ;
				}
				
				else if (cmd == "translate") {               // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> current_translate[i] ;
						}
					}
					mat4 M = Transform::translate(current_translate[0], current_translate[1], current_translate[2]) ;
					mat4 & T = transfstack.top() ;
					M = glm::transpose(M) ;
					T = T * M ;
				}
				else if (cmd == "rotate") {               // change later?
					for (int i = 0; i < 4; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> current_rotate[i] ;
						}
					}
					mat4 M (0) ;
					vec3 axis(current_rotate[0], current_rotate[1], current_rotate[2]) ;
					mat3 R = Transform::rotate(current_rotate[3], axis) ;

					//Setting up M
					
					M[0][0] = R[0][0] ;
					M[0][1] = R[0][1] ;
					M[0][2] = R[0][2] ;
					M[1][0] = R[1][0] ;
					M[1][1] = R[1][1] ;
					M[1][2] = R[1][2] ;
					M[2][0] = R[2][0] ;
					M[2][1] = R[2][1] ;
					M[2][2] = R[2][2] ;	
					M[0][3] = 0 ;
					M[1][3] = 0 ;
					M[2][3] = 0 ;
					M[3][0] = 0 ;
					M[3][1] = 0 ;
					M[3][2] = 0 ;
					M[3][3] = 1 ;

					mat4 & T = transfstack.top() ;
					M = glm::transpose(M) ;
					T = T * M ;
				}
				else if (cmd == "scale") {               // change later?
					for (int i = 0; i < 3; i++) {
						if (s.str().empty()) {
							cerr << "Not enough arguments to 'size'\n";
							throw 2;
						} else {
							s >> current_scale[i] ;
						}
					}
					mat4 M = Transform::scale(current_scale[0], current_scale[1], current_scale[2]) ;
					
					mat4 & T = transfstack.top() ;
					M = glm::transpose(M) ;
					T = T * M ;
				}
		}
		getline (in, str) ; 
    }
  }
  else {
    cerr << "Unable to Open File " << filename << "\n" ; 
    throw 2 ; 
  }
}

int main () {	
	
	Parser("scene1.test");
	
	//eye = vec3((camera[0] - camera[3]), (camera[1] - camera[4]), (camera[2] - camera[5])) ; 
	eye = vec3(camera[0], camera[1], camera[2]) ;  // lookfrom
	up = vec3(camera[6], camera[7], camera[8]) ; 
	center = vec3(camera[3], camera[4], camera[5]);  // lookat
	width = size[0];
	height = size[1];
	
	
	Image image;
	image.initialize();
	
    for (int i = 0; i < height; i++) {            
        for (int j = 0; j < width; j++) {
            float red = 0;
            float green = 0;
            float blue = 0;
            
            float numz2 = 2.0f;
            float numz = pow(numz2, 3) + 4;
/*            numz = 10;
            float y = i;
            float x = j;*/

            //float flat = 4.0f;
            for (float y = i; y < i + 1.0; y += 1.0f/numz2) {
                for (float x = j; x < j + 1.0; x += 1.0f/numz2) {
                    
                    float nDRays = 1.0f;
                    float nIDRays = 1.0f;
                    float nMRays = 1.0f;
                    float nERays = 1.0f;
                    for (int kzoo = 0; kzoo < numz2; kzoo ++) {
                                float randN = ((float) rand() / (float) RAND_MAX);
                                if (randN < probDirect) {
                                    nDRays += 1.0f;
                                }
                                else if ((randN >= probDirect) && (randN < (probDirect+probRand))) {
                                    nIDRays += 1.0f;
                                }
                                else if ((randN >= (probDirect+probRand)) && (randN < (probDirect+probRand+probMirror))) {
                                    nMRays += 1.0f;
                                }
                                else if ((randN >= (probDirect+probRand+probMirror)) && (randN < (probDirect+probRand+probMirror+probEmit))) {
                                    nERays += 1.0f;
                                }
                    }

                    float FnDRays = nDRays;
                    float rnDRays = 0;
                    float bnDRays = 0;
                    float gnDRays = 0;
                    for (int DRAYS = 0; DRAYS < nDRays; DRAYS ++) {
                        Ray ray;
                        ray.RayThruPixel(y, x); // initializes the ray
                        Intersection hit (ray);  // stores the hit object
                        if (hit.flag == true) {
                            Shading sh;
                            sh.DirectShade(hit,0);
                            if ((sh.shadecolor.r == 0) && (sh.shadecolor.g == 0) && (sh.shadecolor.b == 0) && (FnDRays > 1)) {
                                FnDRays -= 1;
                            } else {
                                rnDRays += sh.shadecolor.r*255;
                                gnDRays += sh.shadecolor.g*255;
                                bnDRays += sh.shadecolor.b*255;
                            }
                            }
                        }
                        red += (nDRays/numz) * (1/FnDRays) * rnDRays;
                        blue += (nDRays/numz) * (1/FnDRays) * bnDRays;
                        green += (nDRays/numz) * (1/FnDRays) * gnDRays;

                        float FnIDRays = nIDRays;
                        float rnIDRays = 0;
                        float bnIDRays = 0;
                        float gnIDRays = 0;
                        for (int IDRAYS = 0; IDRAYS < nIDRays; IDRAYS ++) {
                            Ray ray;
                            ray.RayThruPixel(y, x); // initializes the ray
                            Intersection hit (ray);  // stores the hit object
                            if (hit.flag == true) {
                                float randN = ((float) rand() / (float) RAND_MAX);
                                Shading sh;
                                sh.InDirectRandShade(hit,0);
                                if ((sh.shadecolor.r == 0) && (sh.shadecolor.g == 0) && (sh.shadecolor.b == 0) && (FnIDRays > 1)) {
                                    FnIDRays -= 1;
                                } else {
                                    rnIDRays += sh.shadecolor.r*255;
                                    gnIDRays += sh.shadecolor.g*255;
                                    bnIDRays += sh.shadecolor.b*255;
                                }
                            }
                        }
                        red += (nIDRays/numz) * (1/FnIDRays) * rnIDRays;
                        blue += (nIDRays/numz) * (1/FnIDRays) * bnIDRays;
                        green += (nIDRays/numz) * (1/FnIDRays) * gnIDRays;

                        float FnMRays = nMRays;
                        float rnMRays = 0;
                        float bnMRays = 0;
                        float gnMRays = 0;
                        for (int MRAYS = 0; MRAYS < nMRays; MRAYS ++) {
                            Ray ray;
                            ray.RayThruPixel(y, x); // initializes the ray
                            Intersection hit (ray);  // stores the hit object
                            if (hit.flag == true) {
                                Shading sh;
                                sh.MirrorShade(hit,0);
                                if ((sh.shadecolor.r == 0) && (sh.shadecolor.g == 0) && (sh.shadecolor.b == 0) && (FnMRays > 1)) {
                                    FnMRays -= 1;
                                } else {
                                    rnMRays += sh.shadecolor.r*255;
                                    gnMRays += sh.shadecolor.g*255;
                                    bnMRays += sh.shadecolor.b*255;
                                }
                            }
                        }
                        red += (nMRays/numz) * (1/FnMRays) * rnMRays;
                        blue += (nMRays/numz) * (1/FnMRays) * bnMRays;
                        green += (nMRays/numz) * (1/FnMRays) * gnMRays;

                        float FnERays = nERays;
                        float rnERays = 0;
                        float bnERays = 0;
                        float gnERays = 0;
                        for (int ERAYS = 0; ERAYS < nERays; ERAYS ++) {
                            Ray ray;
                            ray.RayThruPixel(y, x); // initializes the ray
                            Intersection hit (ray);  // stores the hit object
                            if (hit.flag == true) {
                                Shading sh;
                                sh.EmitShade(hit,0);
                                if ((sh.shadecolor.r == 0) && (sh.shadecolor.g == 0) && (sh.shadecolor.b == 0) && (FnERays > 1)) {
                                    FnERays -= 1;
                                } else {
                                    rnERays += sh.shadecolor.r*255;
                                    gnERays += sh.shadecolor.g*255;
                                    bnERays += sh.shadecolor.b*255;
                                }
                            }
                        }
                        red += (nERays/numz) * (1/FnERays) * rnERays;
                        blue += (nERays/numz) * (1/FnERays) * bnERays;
                        green += (nERays/numz) * (1/FnERays) * gnERays;
                }
            }
            
            int newRed = (int) red;
            int newGreen = (int) green;
            int newBlue = (int) blue;
            if (newRed > 255){
                newRed = 255;
            }
            if (newRed < 0) {
                newRed = 0;
            }
            if (newGreen > 255) {
                newGreen = 255;
            }
            if (newGreen < 0) {
                newGreen = 0;
            }
            if (newBlue > 255) {
                newBlue = 255;
            }
            if (newBlue < 0) {
                newBlue = 0;
            }
            
            image.colors[i][j].r = newRed;
            image.colors[i][j].g = newGreen;
            image.colors[i][j].b = newBlue;
        }
    }
	image.writeImage();
	return 0;
	
}
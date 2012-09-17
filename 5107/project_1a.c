#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159265

typedef struct {
	float x, y, z;
} Point;

typedef struct {
	float x, y, z;
} Vector;

typedef struct {
	float r, g, b;
} ColorType;

typedef struct {
	float x, y, z;
	float dir_x, dir_y, dir_z;
} RayType;

typedef struct {
	float x, y, z;
	float r;
	ColorType m;
} SphereType;



ColorType Shade_Ray(Point point, ColorType color)
{
    ColorType computed_color;
    
    computed_color = color;
    
    return computed_color;
}

ColorType Trace_Ray(RayType ray, SphereType sphere, ColorType background)
{
	ColorType return_color;
	double a, b, c;
    double dis; //b^2 - 4ac
    double t1, t2, t;
    int hit = 0;
    
    a = 1.0;
    b = 2*(ray.dir_x*(ray.x-sphere.x) + ray.dir_y*(ray.y-sphere.y) + ray.dir_z*(ray.z-sphere.z))
	c = (ray.x-sphere.x)*(ray.x-sphere.x) + (ray.y-sphere.y)*(ray.y-sphere.y) + (ray.z-sphere.z)*(ray.z-sphere.z) - sphere.r*sphere.r;
    dis = b*b - 4*a*c;
    
    if (dis > 0) {//ray pierces the sphere
        t1 = (-1*b + sqrt(dis))/(2*a);
        t2 = (-1*b - sqrt(dis))/(2*a);
        if (t1 <= 0) {
            if (t2 <= 0) {
                hit = 0;
            } else {
                t = t2;
                hit = 1;
            }
        } else {
            if (t2 <= 0) {
                t = t1;
                hit = 1;
            } else {
                if (t1 > t2) {
                    t = t2;
                    hit = 1;
                } else {
                    t = t1;
                    hit = 1;
                }
            }
        }
    } else if (dis = 0) {//ray grazes the sphere
        t = (-1*b)/(2*a);
        if (t > 0) {
            hit = 1;
        }
    } else {//ray misses the sphere
        hit = 0;
    }
    
    if (hit == 1) {
        Point intersection;
        intersection.x = ray.x + ray.dir_x*t;
        intersection.y = ray.y + ray.dir_y*t;
        intersection.z = ray.z + ray.dir_z*t;
        
        return_color = Shade_Ray(Point intersection, ColorType sphere.m);
    } else {
        return_color = background;
    }
    
    return (return_color);
}

int main(int argc, char* const argv[])
{
	char junkchar[20]; //stores the useless strings in input file such as "eye", "viewdir", etc.
	Point eye, ul, ur, ll, lr;
	SphereType sphere;
	Vector vdir, unit_vdir, up, unit_up;
	Vector u, unit_u, v, unit_v;
	ColorType bkg, material;
	float fov_h, aspect_ratio, view_dist;
	int pixwidth, pixheight; //the image width in pixel units
	double width, height; //width and height of the viewing window
	int i, j;
	ColorType **image;
    
	FILE *inFile = NULL;
	FILE *outFile = NULL;
    
	inFile = fopen("input.txt", "r");
	/* read scene description from input file */
	fscanf(inFile, "%s %f %f %f", junkchar, &eye.x, &eye.y, &eye.z);
	fscanf(inFile, "%s %f %f %f", junkchar, &vdir.x, &vdir.y, &vdir.z);    
	fscanf(inFile, "%s %f %f %f", junkchar, &up.x, &up.y, &up.z);
	fscanf(inFile, "%s %f", junkchar, &fov_h);
	fscanf(inFile, "%s %f", junkchar, &aspect_ratio);
	fscanf(inFile, "%s %f", junkchar, &view_dist);
	fscanf(inFile, "%s %d", junkchar, &pixwidth);
	fscanf(inFile, "%s %f %f %f", junkchar, &bkg.r, &bkg.g, &bkg.b);
	fscanf(inFile, "%s %f %f %f", junkchar, &material.r, &material.g, &material.b);
	fscanf(inFile, "%s %f %f %f %f", junkchar, &sphere.x, &sphere.y, &sphere.z, &sphere.r);
    fclose(inFile);
    
	/* determine the width and height of the viewing window */
	width = 2*view_dist*tan((fov_h/2)*PI/180);
	height = width/aspect_ratio;
	
	//printf("%lf %lf\n", width, height);
	/* determine the vector and unit vector in the horizontal direction of the viewing window */
	u.x = vdir.y*up.z - vdir.z*up.y;
	u.y = vdir.z*up.x - vdir.x*up.z;
	u.z = vdir.x*up.y - vdir.y*up.x;
    
	//printf("%f %f %f\n", u.x, u.y, u.z);
    
	unit_u.x = u.x/(sqrt(u.x*u.x+u.y*u.y+u.z*u.z));
	unit_u.y = u.y/(sqrt(u.x*u.x+u.y*u.y+u.z*u.z));
	unit_u.z = u.z/(sqrt(u.x*u.x+u.y*u.y+u.z*u.z));
    
	//printf("%f %f %f\n", unit_u.x, unit_u.y, unit_u.z);
	
	/* determine the vector and unit vector in the vertical direction of the viewing window */
	v.x = unit_u.y*vdir.z - unit_u.z*vdir.y;
	v.y = unit_u.z*vdir.x - unit_u.x*vdir.z;
	v.z = unit_u.x*vdir.y - unit_u.y*vdir.x;
	//printf("%f %f %f\n", v.x, v.y, v.z);
    
	unit_v.x = v.x/(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
	unit_v.y = v.y/(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
	unit_v.z = v.z/(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
	//printf("%f %f %f\n", unit_vdir.x, unit_vdir.y, unit_vdir.z);		
    
	/* determine the unit vector in the view direction and the four corners of the viewing window */
	unit_vdir.x = vdir.x/(sqrt(vdir.x*vdir.x+vdir.y*vdir.y+vdir.z*vdir.z));
	unit_vdir.y = vdir.y/(sqrt(vdir.x*vdir.x+vdir.y*vdir.y+vdir.z*vdir.z));
	unit_vdir.z = vdir.z/(sqrt(vdir.x*vdir.x+vdir.y*vdir.y+vdir.z*vdir.z));
    
	ul.x = eye.x + view_dist*unit_vdir.x + (height/2)*unit_v.x - (width/2)*unit_u.x;
	ul.y = eye.y + view_dist*unit_vdir.y + (height/2)*unit_v.y - (width/2)*unit_u.y;
	ul.z = eye.z + view_dist*unit_vdir.z + (height/2)*unit_v.z - (width/2)*unit_u.z;
    
	ur.x = eye.x + view_dist*unit_vdir.x + (height/2)*unit_v.x + (width/2)*unit_u.x;
	ur.y = eye.y + view_dist*unit_vdir.y + (height/2)*unit_v.y + (width/2)*unit_u.y;
	ur.z = eye.z + view_dist*unit_vdir.z + (height/2)*unit_v.z + (width/2)*unit_u.z;
    
 	ll.x = eye.x + view_dist*unit_vdir.x - (height/2)*unit_v.x - (width/2)*unit_u.x;
	ll.y = eye.y + view_dist*unit_vdir.y - (height/2)*unit_v.y - (width/2)*unit_u.y;
	ll.z = eye.z + view_dist*unit_vdir.z - (height/2)*unit_v.z - (width/2)*unit_u.z;
    
	lr.x = eye.x + view_dist*unit_vdir.x - (height/2)*unit_v.x + (width/2)*unit_u.x;
	lr.y = eye.y + view_dist*unit_vdir.y - (height/2)*unit_v.y + (width/2)*unit_u.y;
	lr.z = eye.z + view_dist*unit_vdir.z - (height/2)*unit_v.z + (width/2)*unit_u.z;
	
	/* determine the height of the out put image */
	pixheight = pixwidth/aspect_ratio;
	//printf("%d %d\n", pixheight, pixwidth);
	
	/* initialize pixel array to store output image */
	image = malloc(sizeof *image * pixheight);
	for (i = 0; i < pixheight; i++)
	{	
		image[i] = malloc(sizeof *image[i] * pixwidth);
	}
	
	for (i = 0; i < pixheight; i++)
		for (j = 0; j < pixwidth; j++)
		{	
			image[i][j] = bkg;
			//image[i][j].r = bkg.r;
			//image[i][j].g = bkg.g;
			//image[i][j].b = bkg.b;
			//printf("%f %f %f\n", image[i][j].r, image[i][j].g, image[i][j].b);
		}
    
	/*	determine the horizontal and vertical offset vectors delta_h, delta_v to use 
     in stepping from point to point across the viewing window from left to right
     and from up to bottom*/
	Vector delta_h, delta_v;
    
	delta_h.x = (ur.x-ul.x)/(pixwidth-1);
	delta_h.y = (ur.y-ul.y)/(pixwidth-1);
	delta_h.z = (ur.z-ul.z)/(pixwidth-1);	
    
	delta_v.x = (ll.x-ul.x)/(pixheight-1);
	delta_v.y = (ll.y-ul.y)/(pixheight-1);
	delta_v.z = (ll.z-ul.z)/(pixheight-1);	
    
	/* initialize rays */
	RayType ray;
	ray.x = eye.x;
	ray.y = eye.y;
	ray.z = eye.z;
    
	/* for each pixel in the image array: */
	for (i = 0; i < pixheight; i++)
		for (j = 0; j < pixwidth; j++)
		{		
            Point window;
			window.x = (ul.x+delta_h.x*j+delta_v.x*i);
			window.y = (ul.y+delta_h.y*j+delta_v.y*i);
			window.z = (ul.z+delta_h.z*j+delta_v.z*i);
            
            Vector ray_dir;
            
            ray_dir.x = window.x - ray.x;
            ray_dir.y = window.y - ray.y;
            ray_dir.z = window.z - ray.z;
            
            ray.dir_x = ray_dir.x/(sqrt(ray_dir.x*ray_dir.x + ray_dir.y*ray_dir.y + ray_dir.z*ray_dir.z));
            ray.dir_z = ray_dir.z/(sqrt(ray_dir.x*ray_dir.x + ray_dir.y*ray_dir.y + ray_dir.z*ray_dir.z));
            ray.dir_z = ray_dir.z/(sqrt(ray_dir.x*ray_dir.x + ray_dir.y*ray_dir.y + ray_dir.z*ray_dir.z));
            
            sphere.m = material;
			image[i][j] = Trace_Ray(ray, sphere, bkg);
		}
    
    outFile = fopen("output_1.ppm", "w");
    fprintf(outFile, "P3\n");
    fprintf(outFile, "%d %d\n", pixwidth, pixheight);
    fprintf(outFile, "# created by Jian Liu's program\n");
    fprintf(outFile, "255\n");
    
    for (i = 0; i < pixheight; i++) {
        for (j = 0; j < pixwidth; j++) {
            fprintf(outFile, "%f %f %f\n", image[i][j].r, image[i][j].g, image[i][j].b);
        }
    }
    
    fclose(outFile);
    /*
     printf("%f %f %f\n", ul.x, ul.y, ul.z);
     printf("%f %f %f\n", ur.x, ur.y, ur.z);	
     printf("%f %f %f\n", ll.x, ll.y, ll.z);
     printf("%f %f %f\n", lr.x, lr.y, lr.z);
     */
    
    /*
     printf("%f %f %f\n", eye.x, eye.y, eye.z);
     printf("%f %f %f\n", vdir.x, vdir.y, vdir.z);
     printf("%f %f %f\n", up.x, up.y, up.z);
     printf("%f\n", fov_h);
     printf("%f\n", aspect_ratio);
     printf("%f\n", view_dist);
     printf("%d\n", pixwidth);
     printf("%d %d %d\n", bkg_r, bkg_g, bkg_b);
     printf("%d %d %d\n", material_r, material_g, material_b);
     printf("%f %f %f %f\n", sphere.x, sphere.y, sphere.z, sphere_r);
     */
    
	return 0;
}

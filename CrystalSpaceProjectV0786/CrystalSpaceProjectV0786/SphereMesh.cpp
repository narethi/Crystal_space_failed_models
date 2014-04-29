#ifndef PI4
#define PI4 atan(1) //1/4 of pi (45)
#endif

#ifndef PI2
#define PI2 atan(1)*2 //half of pi (90)
#endif



#ifndef __CS_PLUGINS_MESH_PLANET_MESH_CPP__
#define __CS_PLUGINS_MESH_PLANET_MESH_CPP__

#include <iostream>
#include <gl/glut.h>
#include <gl/GL.h>
#include <noise/noise.h>
//#include <noise/noiseutils/noiseutils.h>

//using namespace std;
using namespace noise;

/*I am now introducing a sphere structure to keep track of all the elements within the sphere
This will help keep track of how the whole structure is kept together */
struct sphere
{
	double radius, /* Size is used for the draw arrays*/size; 
	int resolution;

	//colors is for testing purposes
	double* colors;

	//This is where all the points will be stored
	double* plane0;
	double* plane1;
	double* plane2;
	double* plane3;
	double* plane4;
	double* plane5;

	double* normal0;
	double* normal1;
	double* normal2;
	double* normal3;
	double* normal4;
	double* normal5;

	//These are the locations of the intersection of planes
	double* intersection0;
	double* intersection1;
	double* intersection2;
	double* intersection3;
	double* intersection4;
	double* intersection5;
};

sphere solveIntersection(double radius, int resolution, sphere world);

//All of these variables are used exclusively for the stand alone version
GLint WIDTH = 600, HEIGHT = 600, /*cRadius is the radius of the camera viewing angle*/
	cRadius = 15, offsetX = 0, offsetY = 0, offsetZ = 0, up = 1;
GLdouble axis[] = {-60,0,0, 60,0,0, 0,-60,0, 0,60,0, 0,0,-60, 0,0,60};
GLdouble colors[] = {1,1,1, 1,0,0, 1,1,1, 0,1,0, 1,1,1, 0,0,1, 1,1,1, 1,0,1, 1,1,1, 1,1,0, 1,1,1, 0,1,1, 1,1,1, 
					.5,.5,.5, 1,1,1, 1,.5,.5, 1,1,1, .5,1,.5, 1,1,1, .5,.5,1, 1,1,1};
GLdouble eye[] = {-300,-300,300};
GLdouble center[] = {0,0,0};
GLdouble cameraAngle = 0, cameraVerticle = 0;
//greg is the test sphere
sphere greg;
//test is the resolution of the test sphere
int test = 0;
//meowsers is the radius of the test sphere
double meowsers = 10;

//moves along the x axis
void camera(void)
{
	
	eye[0] = cRadius*sin(cameraAngle)*cos(cameraVerticle),
	eye[1] = cRadius*sin(cameraVerticle),
	eye[2] = cRadius*cos(cameraAngle)*cos(cameraVerticle);
	if(cameraVerticle < 0) 
		cameraVerticle += 6.28, up = -1;
	else if(cameraVerticle > 6.28) 
		cameraVerticle -= 6.28, up = 1;
	else if(cameraVerticle < 3.14/2 || cameraVerticle > 3*3.14/2) up = 1;
	else if(cameraVerticle > 3.14/2) up = -1;

}

void drawAxis(void)
{

	glVertexPointer(3, GL_DOUBLE, 0, axis);
	glColorPointer(3,GL_DOUBLE,0, colors);
	glDrawArrays(GL_LINES,0,6);

}

sphere initializeSphere(int resolution, double radius)
{
	sphere world;
	//each plane requires and array that can store 5 pieces
	//of information for 2^(resolution + 1) + 2 points of
	//intersection, the 5 pieces of information are the x,y,z
	//values of the point of intersection and where the intersection
	//occurs in the loop
	//The order will be pitch, roll, x,y,z
	int count = (pow(2,resolution+1)+2)*5;
	world.resolution = resolution;
	world.radius = radius;
	world.intersection0 = new double[count];
	world.intersection1 = new double[count];
	world.intersection2 = new double[count];
	world.intersection3 = new double[count];
	world.intersection4 = new double[count];
	world.intersection5 = new double[count];

	count = pow(pow(2, resolution)+1, 2)*6;

	world.colors = new double[count];

	world.plane0 = new double[count];
	world.plane1 = new double[count];
	world.plane2 = new double[count];
	world.plane3 = new double[count];
	world.plane4 = new double[count];
	world.plane5 = new double[count];

	
	world.normal0 = new double[count];
	world.normal1 = new double[count];
	world.normal2 = new double[count];
	world.normal3 = new double[count];
	world.normal4 = new double[count];
	world.normal5 = new double[count];

	world = solveIntersection(radius, resolution, world);

	return world;
}

/*

Generate grid currently builds a full grid of vertices using the
radius and resolution as defined by the user. The grid is made of
6 sides that are expanded to appear to form a spherical surface

Known bugs:

	Currently the biggest issue is that after the half way point the
	values for the edge vertices begin to overlap

*/
sphere generateGrid(double radius, int resolution /*,heightMap*/)
{
	sphere world = initializeSphere(resolution, radius);
	//first fix from previous version

	/*earlier I had made the mistake of converting from degrees to radians
	every step so I decided to just start in radians, so step is started in radians*/
	double step = PI2/pow(2, resolution), check = 0;
	int count = 0, stageIntersection = (pow(2, resolution) + 1)*5;

	/*
	
	The plan

	2 nested for loops that are span 90 degrees, making a grid that is exactly the size
	of one side of a sphere

	Another fix to be done later for the heightmapping is during the generation of the intersection values
	an average of the heightmap values will taken into account to make the move from 1 plane to the other
	more fluid and less noticable.

	*/

	for(double i = 0; i < pow(2, resolution); i++)
	{
		for(double j = 0; j < pow(2, resolution) + 1; j++)
		{
			if(j == world.intersection0[int(i*5) + 1] || j == world.intersection0[int(i*5) + stageIntersection + 1])
			{
				world.colors[count]		= j/(pow(2, resolution) + 1);
				world.colors[count + 1] = 0;
				world.colors[count + 2] = 0;

				world.colors[count + 3] = (j + 1)/(pow(2, resolution) + 1);
				world.colors[count + 4] = 0;
				world.colors[count + 5] = 0;

				if((int(i)%2 == 0 && j == world.intersection0[int(i*5) + 1])
					|| (int(i)%2 != 0 && j == world.intersection0[int(i*5) + stageIntersection + 1])) check = 0;
				else check = stageIntersection;
				
				/*Plane0*/
				world.plane0[count]		= world.intersection0[int(i*5 + check) + 2];
				world.plane0[count + 1] = world.intersection0[int(i*5 + check) + 3];
				world.plane0[count + 2] = world.intersection0[int(i*5 + check) + 4];

				world.plane0[count + 3] = world.intersection0[int(i*5  + check) + 7];
				world.plane0[count + 4] = world.intersection0[int(i*5 + check) + 8];
				world.plane0[count + 5] = world.intersection0[int(i*5 + check) + 9];
		
				/*Plane1*/
				world.plane1[count]		= -world.intersection1[int(i*5 + check) + 2];
				world.plane1[count + 1] = world.intersection1[int(i*5 + check) + 3];
				world.plane1[count + 2] = world.intersection1[int(i*5 + check) + 4];

				world.plane1[count + 3] = -world.intersection1[int(i*5 + check) + 7];
				world.plane1[count + 4] = world.intersection1[int(i*5 + check) + 8];
				world.plane1[count + 5] = world.intersection1[int(i*5 + check) + 9];

				/*Plane2*/
				world.plane2[count]		= world.intersection2[int(i*5 + check) + 2];
				world.plane2[count + 1] = world.intersection2[int(i*5 + check) + 3];
				world.plane2[count + 2] = -world.intersection2[int(i*5 + check) + 4];

				world.plane2[count + 3] = world.intersection2[int(i*5 + check) + 7];
				world.plane2[count + 4] = world.intersection2[int(i*5 + check) + 8];
				world.plane2[count + 5] = -world.intersection2[int(i*5 + check) + 9];

				/*Plane3*/
				world.plane3[count]		= world.intersection3[int(i*5 + check) + 2];
				world.plane3[count + 1] = world.intersection3[int(i*5 + check) + 3];
				world.plane3[count + 2] = world.intersection3[int(i*5 + check) + 4];

				world.plane3[count + 3] = world.intersection3[int(i*5 + check) + 7];
				world.plane3[count + 4] = world.intersection3[int(i*5 + check) + 8];
				world.plane3[count + 5] = world.intersection3[int(i*5 + check) + 9];

				/*Plane5*/
				world.plane5[count]		= world.intersection5[int(i*5 + check) + 2];
				world.plane5[count + 1] = -world.intersection5[int(i*5 + check) + 3];
				world.plane5[count + 2] = world.intersection5[int(i*5 + check) + 4];

				world.plane5[count + 3] = world.intersection5[int(i*5 + check) + 7];
				world.plane5[count + 4] = -world.intersection5[int(i*5 + check) + 8];
				world.plane5[count + 5] = world.intersection5[int(i*5 + check) + 9];

				if((int(i)%2 == 0 && j == world.intersection0[int(i*5) + 1])
					|| (int(i)%2 != 0 && j == world.intersection0[int(i*5) + stageIntersection + 1])) check = stageIntersection;
				else check = 0;

				/*Plane4*/
				world.plane4[count]		= world.intersection4[int(i*5 + check)+ 2];
				world.plane4[count + 1] = world.intersection4[int(i*5 + check) + 3];
				world.plane4[count + 2] = world.intersection4[int(i*5 + check) + 4];

				world.plane4[count + 3] = world.intersection4[int(i*5 + check) + 7];
				world.plane4[count + 4] = world.intersection4[int(i*5 + check) + 8];
				world.plane4[count + 5] = world.intersection4[int(i*5 + check) + 9];

				count += 6;
				//This is neccesary for the asymmetric nature of this problem
				//Since there is spaces missing if you use j expressly (the removal of the
				//overlap)
			}
			else if(j > world.intersection0[int(i*5) + 1] && j < world.intersection0[int(i*5) + stageIntersection + 1])
			{

				world.colors[count]		= j/(pow(2, resolution) + 1);
				world.colors[count + 1] = 0;
				world.colors[count + 2] = 0;

				world.colors[count + 3] = (j + 1)/(pow(2, resolution) + 1);
				world.colors[count + 4] = 0;
				world.colors[count + 5] = 0;

				//Plane 0

				if(int(i)%2 == 0) check = j*step;
				else check = (pow(2, resolution) - j)*step;

				world.plane0[count]		= radius*sin(check + PI4)*sin(i*step + 7*PI4);
				world.plane0[count + 1] = radius*sin(check + PI4)*cos(i*step + 7*PI4); 
				world.plane0[count + 2] = radius*cos(check + PI4);

				world.plane0[count + 3] = radius*sin(check + PI4)*sin((i + 1)*step + 7*PI4);
				world.plane0[count + 4] = radius*sin(check + PI4)*cos((i + 1)*step + 7*PI4); 
				world.plane0[count + 5] = radius*cos(check + PI4);

				//Plane 1
				world.plane1[count]		= radius*sin(check + PI4)*sin(i*step + 3*PI4);
				world.plane1[count + 1] = radius*sin(check + PI4)*cos(i*step + 3*PI4); 
				world.plane1[count + 2] = radius*cos(check + PI4);

				world.plane1[count + 3] = radius*sin(check + PI4)*sin((i + 1)*step + 3*PI4);
				world.plane1[count + 4] = radius*sin(check + PI4)*cos((i + 1)*step + 3*PI4); 
				world.plane1[count + 5] = radius*cos(check + PI4);

				//Plane 2
				world.plane2[count]		= radius*sin(check + PI4)*cos(i*step + 7*PI4);
				world.plane2[count + 1] = radius*cos(check + PI4); 
				world.plane2[count + 2] = radius*sin(check + PI4)*sin(i*step + 7*PI4);

				world.plane2[count + 3] = radius*sin(check + PI4)*cos((i + 1)*step + 7*PI4);
				world.plane2[count + 4] = radius*cos(check + PI4); 
				world.plane2[count + 5] = radius*sin(check + PI4)*sin((i + 1)*step + 7*PI4);

				//Plane 3
				world.plane3[count]		= radius*sin(check + PI4)*cos(i*step + 3*PI4);
				world.plane3[count + 1] = radius*cos(check + PI4); 
				world.plane3[count + 2] = radius*sin(check + PI4)*sin(i*step + 3*PI4);

				world.plane3[count + 3] = radius*sin(check + PI4)*cos((i + 1)*step + 3*PI4);
				world.plane3[count + 4] = radius*cos(check + PI4); 
				world.plane3[count + 5] = radius*sin(check + PI4)*sin((i + 1)*step + 3*PI4);

				//Plane 4
				world.plane4[count] = radius*cos(check + PI4);
				world.plane4[count + 1] = radius*sin(check + PI4)*cos(i*step + PI4); 
				world.plane4[count + 2] = radius*sin(check + PI4)*sin(i*step + PI4);

				world.plane4[count + 3] = radius*cos(check + PI4);
				world.plane4[count + 4] = radius*sin(check + PI4)*cos((i + 1)*step + PI4); 
				world.plane4[count + 5] = radius*sin(check + PI4)*sin((i + 1)*step + PI4);

				//Plane 5
				world.plane5[count] = radius*cos(check + 5*PI4);
				world.plane5[count + 1] = radius*sin(check + 5*PI4)*cos(i*step + PI4); 
				world.plane5[count + 2] = radius*sin(check + 5*PI4)*sin(i*step + PI4);

				world.plane5[count + 3] = radius*cos(check + 5*PI4);
				world.plane5[count + 4] = radius*sin(check + 5*PI4)*cos((i + 1)*step + PI4); 
				world.plane5[count + 5] = radius*sin(check + 5*PI4)*sin((i + 1)*step + PI4);

				count += 6;
			}
		}
	}
	world.size = count;
	return world;
}

/*

Pitch and Roll, which tell the map where each intersection occurs
(currently these values are k and pitch in the current version of the
solver). So the pitch tells which row the intersection occurs and the
roll tells which segment the intersection occurs

x,y,z which tells where the intersection occurs this will replace the end
points that are currently being used at the intersection point. A predicted
bug being that when the terrain is mapped and the radius at that point is no
longer default. The rememdy of this is to double the radius of the closest adjacent
point and then average the values of the intersecting point and its adjacent point

known bugs:
	
	Resolution 2 is still broken

	The intersection is not tracked as effectively as possible and therefore
	higher resolutions take a long time to finish

*/
sphere solveIntersection(double radius, int resolution, sphere world)
{
	double x,y,z, x1, y1, z1, incX, incY, incZ, incX1, incY1, incZ1;
	double nextIntersection = 0, endOfCheck = pow(2, resolution - 1), end = 1;
	//init is the variable that will be used to track where the 
	//first intersection occurs. This is important for the initial
	//resolutions because the intersections occur in the
	//first line segments

	bool init = false, intersect = true;
	int count = (pow(2,resolution+1)+2)*5, position;

	//This is the intersection machine currently it only
	//Solves the equations of all of the lines that make up
	//The sphere at each resolution
	for(int pitch = 0; pitch < pow(2, resolution - 1) + 1; pitch++)
	{
		for(int res = 0; res <  pow(2, resolution - 1); res++)
		{
			incX = (radius*cos(3*PI4 - (nextIntersection + 1)*(PI2/pow(2, resolution))) - radius*cos(3*PI4 - nextIntersection*(PI2/pow(2, resolution))))/(500);
			incY = (radius*sin(3*PI4 - (nextIntersection + 1)*(PI2/pow(2, resolution)))*cos(PI4) - radius*sin(3*PI4 - nextIntersection*(PI2/pow(2, resolution)))*cos(PI4))/(500);
			incZ = (radius*sin(3*PI4 - (nextIntersection + 1)*(PI2/pow(2, resolution)))*sin(PI4) - radius*sin(3*PI4 - nextIntersection*(PI2/pow(2, resolution)))*sin(PI4))/(500);

			incX1 = (radius*sin(PI4 + res*(PI2/pow(2, resolution)))*sin(7*PI4 + pitch*PI2/pow(2, resolution)) - radius*sin(PI4 + (res + 1)*(PI2/pow(2, resolution)))*sin(7*PI4 + pitch*PI2/pow(2, resolution)))/(500);
			incY1 = (radius*sin(PI4 + res*(PI2/pow(2, resolution)))*cos(7*PI4 + pitch*PI2/pow(2, resolution)) -  radius*sin(PI4 + (res + 1)*(PI2/pow(2, resolution)))*cos(7*PI4 + pitch*PI2/pow(2, resolution)))/(500);
			incZ1 = (radius*cos(PI4 + res*(PI2/pow(2, resolution))) - radius*cos(PI4 + (res + 1)*(PI2/pow(2, resolution))))/(500);

			x = radius*cos(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution)))
				, y = radius*sin(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution)))*cos(PI4), 
				z = radius*sin(PI4)*sin(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution)));
			x1 = radius*sin(PI4 + res*(PI2/pow(2, resolution)))*sin(7*PI4 + pitch*PI2/pow(2, resolution)), 
				y1 = radius*sin(PI4 + res*(PI2/pow(2, resolution)))*cos(7*PI4 + pitch*PI2/pow(2, resolution)), 
				z1 = radius*cos(PI4 + res*(PI2/pow(2, resolution)));

			glColor3f((res+1)%3,(res+1)%2,0);
			//This part of the code looks for the initial intersection of the planes and then maps them
			//to the sphere object
			if(init) end = 500;
			for(int j = 0; j < end; j++)
			{
				for(double i = 0; i < 500; i++)
				{
					if((abs(x - x1) < (radius/(500*(resolution + 1)))) && 
						(abs(y - y1) < (radius/(500*(resolution + 1))) && 
						(abs(z - z1) < (radius/500*(resolution + 1)))))
					{
						init = true;
						//Plane 0
						//first corner
						position = (pow(2, resolution) - pitch)*5;
						world.intersection0[pitch*5] = pitch, world.intersection0[pitch*5 + 1] = res, 
							world.intersection0[pitch*5 + 2] = x, world.intersection0[pitch*5 + 3] = y, world.intersection0[pitch*5 + 4] = z;

						world.intersection0[position] = pow(2, resolution) - pitch, world.intersection0[position + 1] = res, 
							world.intersection0[position + 2] = -x, world.intersection0[position + 3] = y, world.intersection0[position + 4] = z;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection0[position] = pitch, world.intersection0[position + 1] = pow(2, resolution) - res, 
							world.intersection0[position + 2] = x, world.intersection0[position + 3] = y, world.intersection0[position + 4] = -z;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection0[position] = pow(2, resolution) - pitch, world.intersection0[position + 1] = pow(2, resolution) - res, 
							world.intersection0[position + 2] = -x, world.intersection0[position + 3] = y, world.intersection0[position + 4] = -z;

						//Plane 1
						position = (pow(2, resolution) - pitch)*5;
				
						world.intersection1[pitch*5] = pitch, world.intersection1[pitch*5 + 1] = res, 
							world.intersection1[pitch*5 + 2] = x, world.intersection1[pitch*5 + 3] = -y, world.intersection1[pitch*5 + 4] = z;

						world.intersection1[position] = pow(2, resolution) - pitch, world.intersection1[position + 1] = res, 
							world.intersection1[position + 2] = -x, world.intersection1[position + 3] = -y, world.intersection1[position + 4] = z;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection1[position] = pitch, world.intersection1[position + 1] = pow(2, resolution) - res, 
							world.intersection1[position + 2] = x, world.intersection1[position + 3] = -y, world.intersection1[position + 4] = -z;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection1[position] = pow(2, resolution) - pitch, world.intersection1[position + 1] = pow(2, resolution) - res, 
							world.intersection1[position + 2] = -x, world.intersection1[position + 3] = -y, world.intersection1[position + 4] = -z;

								
						//Plane 2
						position = (pow(2, resolution) - pitch)*5;

						world.intersection2[pitch*5] = pitch, world.intersection2[pitch*5 + 1] = res, 
							world.intersection2[pitch*5 + 2] = z, world.intersection2[pitch*5 + 3] = y, world.intersection2[pitch*5 + 4] = -x;

						world.intersection2[position] = pow(2, resolution) - pitch, world.intersection2[position + 1] = res, 
							world.intersection2[position + 2] = z, world.intersection2[position + 3] = y, world.intersection2[position + 4] = x;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection2[position] = (pow(2, resolution) + 1) + pitch, world.intersection2[position + 1] = pow(2, resolution) - res, 
							world.intersection2[position + 2] = z, world.intersection2[position + 3] = -y, world.intersection2[position + 4] = -x;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection2[position] = (pow(2, resolution + 1) + 1) - pitch, world.intersection2[position + 1] = pow(2, resolution) - res, 
							world.intersection2[position + 2] = z, world.intersection2[position + 3] = -y, world.intersection2[position + 4] = x;

						//Plane 3
						position = (pow(2, resolution) - pitch)*5;

						world.intersection3[pitch*5] = pitch, world.intersection3[pitch*5 + 1] = res, 
							world.intersection3[pitch*5 + 2] = -z, world.intersection3[pitch*5 + 3] = y, world.intersection3[pitch*5 + 4] = -x;

						world.intersection3[position] = pow(2, resolution) - pitch, world.intersection3[position + 1] = res, 
							world.intersection3[position + 2] = -z, world.intersection3[position + 3] = y, world.intersection3[position + 4] = x;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection3[position] = pitch, world.intersection3[position + 1] = pow(2, resolution) - res, 
							world.intersection3[position + 2] = -z, world.intersection3[position + 3] = -y, world.intersection3[position + 4] = -x;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection3[position] = pow(2, resolution) - pitch, world.intersection3[position + 1] = pow(2, resolution) - res, 
							world.intersection3[position + 2] = -z, world.intersection3[position + 3] = -y, world.intersection3[position + 4] = x;

						//Plane 4
						position = (pow(2, resolution) - pitch)*5;

						world.intersection4[pitch*5] = pitch, world.intersection4[pitch*5 + 1] = res, 
							world.intersection4[pitch*5 + 2] = -y, world.intersection4[pitch*5 + 3] = -x, world.intersection4[pitch*5 + 4] = z;

						world.intersection4[position] = pow(2, resolution) - pitch, world.intersection4[position + 1] = res, 
							world.intersection4[position + 2] = -y, world.intersection4[position + 3] = x, world.intersection4[position + 4] = z;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection4[position] = pitch, world.intersection4[position + 1] = pow(2, resolution) - res, 
							world.intersection4[position + 2] = y, world.intersection4[position + 3] = -x, world.intersection4[position + 4] = z;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection4[position] = pow(2, resolution) - pitch, world.intersection4[position + 1] = pow(2, resolution) - res, 
							world.intersection4[position + 2] = y, world.intersection4[position + 3] = x, world.intersection4[position + 4] = z;

						//Plane 5
						position = (pow(2, resolution) - pitch)*5;

						world.intersection5[pitch*5] = pitch, world.intersection5[pitch*5 + 1] = res, 
							world.intersection5[pitch*5 + 2] = -y, world.intersection5[pitch*5 + 3] = -x, world.intersection5[pitch*5 + 4] = -z;

						world.intersection5[position] = pow(2, resolution) - pitch, world.intersection5[position + 1] = res, 
							world.intersection5[position + 2] = -y, world.intersection5[position + 3] = x, world.intersection5[position + 4] = -z;

						position = (pow(2, resolution) + pitch + 1)*5;

						world.intersection5[position] = pitch, world.intersection5[position + 1] = pow(2, resolution) - res, 
							world.intersection5[position + 2] = y, world.intersection5[position + 3] = -x, world.intersection5[position + 4] = -z;

						position = (pow(2, resolution + 1) - pitch + 1)*5;

						world.intersection5[position] = pow(2, resolution) - pitch, world.intersection5[position + 1] = pow(2, resolution) - res, 
							world.intersection5[position + 2] = y, world.intersection5[position + 3] = x, world.intersection5[position + 4] = -z;

						res = pow(2, resolution - 1);
						j = 500;
						i = 500;

					}
					if (init) x+=incX, y+=incY, z+=incZ;
					else x+=incX, y+=incY, z+=incZ, x1-=incX1, y1-=incY1, z1-=incZ1;
				}
				x1-=incX1, y1-=incY1, z1-=incZ1;
				x = radius*cos(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution))), 
				y = radius*sin(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution)))*cos(PI4), 
				z = radius*sin(PI4)*sin(3*PI4 - (nextIntersection)*(PI2/pow(2, resolution)));
			}
			/*

			This is where the interesting stuff happens, essentially this will
			be an asymmetric nested for loop, where the computer checks each 
			line segments individually

			*/
			/*This conditional is very important this is where the asymmetry of the problem
			becomes visible. Essentially after this section res and nextintersection diverge
			so that each position on the curve in the intersecting plane can
			be compared to the intersecting arc*/
			if(!init) nextIntersection++;
			//This statement forces 1 step back if no intersection occurs on the current pitch value
			//ALL pitch values have an intersection
			else if(res == pow(2, resolution - 2) - 1)
			{
				intersect = false;
			}

		}
		//This is the part in the main loop that forces the current pitch to be searched in the previous section
		if(intersect || resolution == 2) nextIntersection++;
		else
		{
			intersect = true;
			pitch--;
			nextIntersection--;
		}
	}
	return world;
}


void display(void)
{

	glutInitDisplayMode(GL_DEPTH|GL_FLOAT|GL_RGBA);
	glEnable(GL_DEPTH_TEST);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
	glLoadIdentity();

	glClearColor(.5,.5,.5,0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glViewport(0, 0, WIDTH, HEIGHT);
	glFrustum(-cRadius, cRadius, -cRadius, cRadius, 2, 2);
	glOrtho(offsetX-cRadius, offsetX+cRadius, offsetY-cRadius, offsetY+cRadius, offsetZ-cRadius*2, offsetZ+cRadius*2);
	gluPerspective(0,1,-cRadius,cRadius);
	gluLookAt(eye[0]+offsetX,eye[1]+offsetY,eye[2]+offsetZ,center[0]+offsetX,center[1]+offsetY,center[2]+offsetZ,0,up,0);//Be careful if the scene is on the lense then it will not be transformable

	glVertexPointer(3, GL_DOUBLE, 0, greg.plane0);
	glColorPointer(3, GL_DOUBLE, 0, greg.colors);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));
	glVertexPointer(3, GL_DOUBLE, 0, greg.plane1);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));
	glVertexPointer(3, GL_DOUBLE, 0, greg.plane2);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));
	glVertexPointer(3, GL_DOUBLE, 0, greg.plane3);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));
	glVertexPointer(3, GL_DOUBLE, 0, greg.plane4);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));
	glVertexPointer(3, GL_DOUBLE, 0, greg.plane5);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, (greg.size/3));

	drawAxis();

	glFlush();
	glutPostRedisplay();

}

void reshape(int x, int y)
{
	WIDTH = x;
	HEIGHT = x;
}

void keyboard(unsigned char key, int g,int h)
{
	if(key =='b')
	{
		meowsers++;
		greg = generateGrid(meowsers, test);
		std::cout << "Resolution:  " << greg.intersection0[int(pow(2, test + 1) + 2)*5] << std::endl;
	}
	else if(key == 'n')
	{
		meowsers--;
		greg = generateGrid(meowsers, test);
		std::cout << "Resolution:  " << greg.intersection0[int(pow(2, test + 1) + 2)*5] << std::endl;
	}
	if(key == '3') cRadius += 1;
	else if(key == '4') cRadius -= 1;

	if(key == '+') 
	{
		test++;
		greg = generateGrid(meowsers, test);
		std::cout << "Resolution:  " << greg.resolution << std::endl;
	}

	else if(key == '-' && test > 0)
	{
		test--;
		greg = generateGrid(meowsers, test);
		std::cout << "Resolution:  " << greg.intersection0[int(pow(2, test + 1) + 2)*5] << std::endl;
	}

	if(key == 'u') cameraAngle+=.01;
	else if(key == 'y') cameraAngle -=.01;
	else if(key == 'o') cameraVerticle -= .01;
	else if(key == 'i') cameraVerticle+=.01;

	if(key == 'a') offsetX--;
	else if(key == 'd') offsetX++;

	if(key == 'w') offsetY--;
	else if(key == 's') offsetY++;

	camera();
	
}

void main(int arguV, char* argc[], char** argv)
{
	greg = generateGrid(meowsers, test);
	camera();

	module::Perlin myModule;
	double value = myModule.GetValue(5, 0.75, 0.50);
	std::cout << value << std::endl;
	myModule.SetFrequency(myModule.GetValue(5, 0.75, 0.50));
	value = myModule.GetValue(5, 0.75, 0.50);
	std::cout << value << std::endl;
	value = myModule.GetValue(3, 0.75, 0.50);
	std::cout << value << std::endl;

	glutInit(&arguV,argc);
	glutInitWindowSize(WIDTH,HEIGHT);
	glutCreateWindow("Crystal Space Planet Simulation V0.7");
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);	
	glutReshapeFunc(reshape);
	glutMainLoop();
}

#endif
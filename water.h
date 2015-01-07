#ifndef __WATER__H_
#define __WATER__H_

#include<GL/glew.h>
#include<GLFW/glfw3.h>

/* Structures  */
struct Particle{
	int id;
	GLfloat px,py,pz;
	GLfloat vx,vy,vz;
	GLfloat fx,fy,fz;
	GLfloat rho;
	GLfloat lambda;
};
typedef struct Particle Particle;
typedef struct Particle* Particle_t;

struct Bubble{
	int type;
	float radius;
	GLfloat x,y,z;
	GLfloat vx,vy,vz;
	GLfloat lifetime;
};
typedef struct Bubble Bubble;
typedef struct Bubble* Bubble_t;

struct Triple{
	GLfloat x,y,z;
};
typedef struct Triple Triple;
typedef struct Triple* Triple_t;

struct Indices{
	int ix,iy,iz;
};
typedef struct Indices Indices;
typedef struct Indices* Indices_t;

#endif




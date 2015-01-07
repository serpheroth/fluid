#ifndef __GRIDS__
#define __GRIDS__

#include<vector>
#include<set>
#include"water.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Grid{
public:
Grid();
Grid(Grid& g);
Grid(GLfloat xmin,GLfloat ymin,GLfloat zmin,GLfloat xmax,GLfloat ymax,GLfloat zmax,GLfloat gap,GLfloat h);
~Grid();
void initialize(GLfloat xmin,GLfloat ymin,GLfloat zmin,GLfloat xmax,GLfloat ymax,GLfloat zmax,GLfloat gap,GLfloat h);
void add(Particle_t t);
void remove(Particle_t t);
void getNeighbors(Particle_t t,std::vector<Indices>& neighbors);
inline std::set<Particle_t>* get(int i,int j,int k){return grid[i][j][k];};
void adjust();
glm::vec3 getNormal(Particle_t t);
int getDensity(float x,float y,float z);
glm::vec3 getVelocity(float x,float y,float z);
private:
void normalize(int& a,int& b,int& c);


std::set<Particle_t> ****grid; 
int xnum,ynum,znum;
GLfloat gap,xmin,xmax,ymin,ymax,zmin,zmax,h;
};


#endif

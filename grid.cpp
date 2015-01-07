#include <cmath>
#include "grid.h"
#include "sph.h"

Grid::Grid(){
	xmin=xmax=ymin=ymax=zmin=zmax=gap=h=0.0f;
	grid=NULL;		
	xnum=ynum=znum=0;
}

void Grid::normalize(int& a,int& b,int& c){
	if(a<0)
		a=0;
	if(a>xnum-1)
		a=xnum-1;
	if(b<0)
		b=0;
	if(b>ynum-1)
		b=ynum-1;
	if(c<0)
		c=0;
	if(c>znum-1)
		c=znum-1;
}

Grid::Grid(Grid& g):xmin(g.xmin),xmax(g.xmax),ymin(g.ymin),ymax(g.ymax),zmin(g.zmin),zmax(g.zmax),h(g.h){}

Grid::Grid(GLfloat xmin,GLfloat ymin,GLfloat zmin,GLfloat xmax,GLfloat ymax,GLfloat zmax,GLfloat gap,GLfloat h){
	initialize(xmin,ymin,zmin,xmax,ymax,zmax,gap,h);
}

void Grid::initialize(GLfloat xmin,GLfloat ymin,GLfloat zmin,GLfloat xmax,GLfloat ymax,GLfloat zmax,GLfloat gap,GLfloat h){
	this->xmin=xmin;
	this->xmax=xmax;
	this->ymin=ymin;
	this->ymax=ymax;
	this->zmin=zmin;
	this->zmax=zmax;
	this->xnum=ceil((xmax-xmin)/gap);
	this->ynum=ceil((ymax-ymin)/gap);
	this->znum=ceil((zmax-zmin)/gap);
	this->gap=gap;	
	this->h=h;
	this->grid=new std::set<Particle_t>***[xnum];
	for(int i=0;i<xnum;i++){
		this->grid[i]=new std::set<Particle_t>**[ynum];
		for(int j=0;j<ynum;j++){
			this->grid[i][j]=new std::set<Particle_t>*[znum];
			for(int k=0;k<znum;k++){
				this->grid[i][j][k]=new std::set<Particle_t>();
			}
		}
	}
}

Grid::~Grid(){
	for(int i=0;i<xnum;i++){
		for(int j=0;j<ynum;j++){
			for(int k=0;k<znum;k++){
				grid[i][j][k]->clear();
				delete grid[i][j][k];
			}
			delete[] grid[i][j];
		}
		delete[] grid[i];
	}
	delete[] grid;
}


void Grid::adjust(){
	for(int i=0;i<xnum;i++){
		for(int j=0;j<ynum;j++){
			for(int k=0;k<znum;k++){
				std::set<Particle_t> *m=grid[i][j][k];
				for(std::set<Particle_t>::iterator iter=m->begin();iter!=m->end();){
					int ix=((*iter)->px-xmin)/gap;
					int iy=((*iter)->py-ymin)/gap;
					int iz=((*iter)->pz-zmin)/gap;
					normalize(ix,iy,iz);
					if(ix!=i||iy!=j||iz!=k){
						Particle_t particle=*iter;
						m->erase(iter++);
						grid[ix][iy][iz]->insert(particle);	
					}else
						iter++;
						
				}	
			}
		}			
	}
}


void Grid::add(Particle_t t){
	int ix=(t->px-xmin)/gap;
	int iy=(t->py-ymin)/gap;
	int iz=(t->pz-zmin)/gap;

	normalize(ix,iy,iz);
	
	grid[ix][iy][iz]->insert(t);	

}


void Grid::remove(Particle_t t){
	int ix=(t->px-xmin)/gap;
	int iy=(t->py-ymin)/gap;
	int iz=(t->pz-zmin)/gap;
	
	normalize(ix,iy,iz);
	
	grid[ix][iy][iz]->erase(grid[ix][iy][iz]->find(t));	
}

void Grid::getNeighbors(Particle_t t,std::vector<Indices>& neighbors){
	int ix=(t->px-xmin)/gap;
	int iy=(t->py-ymin)/gap;
	int iz=(t->pz-zmin)/gap;

	normalize(ix,iy,iz);

	int ixmin=(t->px-h-xmin)/gap;
	int ixmax=(t->px+h-xmin)/gap;
	int iymin=(t->py-h-ymin)/gap;
	int iymax=(t->py+h-ymin)/gap;
	int izmin=(t->pz-h-zmin)/gap;
	int izmax=(t->pz+h-zmin)/gap;


	normalize(ixmin,iymin,izmin);
	normalize(ixmax,iymax,izmax);

	
	for(int i=ixmin;i<ixmax+1;i++){
		for(int j=iymin;j<iymax+1;j++){
			for(int k=izmin;k<izmax+1;k++){
				Indices indices;
				indices.ix=i;indices.iy=j;indices.iz=k;
				neighbors.push_back(indices);	
			}
		}
	}
}

glm::vec3 Grid::getNormal(Particle_t t){
	int d=getDensity(t->px,t->py,t->pz);
	int xplus=getDensity(t->px+h,t->py,t->pz);
	int yplus=getDensity(t->px,t->py+h,t->pz);
	int zplus=getDensity(t->px,t->py,t->pz+h);
	return glm::normalize(glm::vec3(xplus-d,yplus-d,zplus-d));
}


int Grid::getDensity(float x,float y,float z){
	if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax)
		return 0;
	int ix=(x-xmin)/gap;
	int iy=(y-ymin)/gap;
	int iz=(z-zmin)/gap;
	
	
	int ixmin=(x-h-xmin)/gap;
	int ixmax=(x+h-xmin)/gap;
	int iymin=(y-h-ymin)/gap;
	int iymax=(y+h-ymin)/gap;
	int izmin=(z-h-zmin)/gap;
	int izmax=(z+h-zmin)/gap;


	normalize(ixmin,iymin,izmin);
	normalize(ixmax,iymax,izmax);

	int result=0;
	for(int i=ixmin;i<ixmax+1;i++){
		for(int j=iymin;j<iymax+1;j++){
			for(int k=izmin;k<izmax+1;k++){
				std::set<Particle_t>* s=grid[i][j][k];
				for(std::set<Particle_t>::iterator iter=s->begin();iter!=s->end();iter++){
					Particle_t p=*iter;
					float dx=p->px-x;
					float dy=p->py-y;
					float dz=p->pz-z;
					float r=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
					if(h>r)
						result++;
				}
			}
		}
	}
	return result;
}

glm::vec3 Grid::getVelocity(float x,float y,float z){
	if(x<xmin||x>xmax||y<ymin||y>ymax||z<zmin||z>zmax)
		return glm::vec3(0.0f,0.0f,0.0f);
	int ix=(x-xmin)/gap;
	int iy=(y-ymin)/gap;
	int iz=(z-zmin)/gap;

	int ixmin=(x-h-xmin)/gap;
	int ixmax=(x+h-xmin)/gap;
	int iymin=(y-h-ymin)/gap;
	int iymax=(y+h-ymin)/gap;
	int izmin=(z-h-zmin)/gap;
	int izmax=(z+h-zmin)/gap;


	normalize(ixmin,iymin,izmin);
	normalize(ixmax,iymax,izmax);

	float d=0.0f;
	glm::vec3 n(0.0f,0.0f,0.0f);
	for(int i=ixmin;i<ixmax+1;i++){
		for(int j=iymin;j<iymax+1;j++){
			for(int k=izmin;k<izmax+1;k++){
				std::set<Particle_t>* s=grid[i][j][k];
				for(std::set<Particle_t>::iterator iter=s->begin();iter!=s->end();iter++){
					Particle_t p=*iter;
					float dx=p->px-x;
					float dy=p->py-y;
					float dz=p->pz-z;
					float rsqr=pow(dx,2)+pow(dy,2)+pow(dz,2);
					if(hsqr>rsqr){
						float tmp=wk* pow(hsqr - rsqr, 3);
						n=n+glm::vec3((*iter)->vx*tmp,(*iter)->vy*tmp,(*iter)->vz*tmp);
						d=d+wk * pow(hsqr - rsqr, 3);
					}
				}
			}
		}
	}
	if(d==0.0f)
		return glm::vec3(0.0f,0.0f,0.0f);
	n/=d;
	return n;
	
}








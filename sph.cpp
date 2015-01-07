#include<list>
#include<iostream>
#include<mutex>
#include"sph.h"
#include"util.h"


// for bubbles
float rmin=0.001;
float rmax=0.008;


/* Variables */
float *pPos=NULL;
Particle_t *particles=NULL;
Grid grid;
std::vector<std::thread> threads;
std::list<Bubble> bubbles;

// for walls
float d=0.2;
float wxmin=-d,wxmax=2*d,wymin=-2*d,wymax=d,wzmin=-d,wzmax=d;

std::mutex print_mutex;

// function declarations
Triple getOmega(int i, int j);
float getC(int i);
float calcW(int i, int j);
Triple calcDW(int i, int j);
void predictPos();
void updateVelocityPos();
void detect_collision();
void computeRho();
void computeLambda();
void computeDeltap();
void computeVorticity();
void updateBubbles();
void generateBubbles();
float Clamp(float v,float min,float max);

void start(){
	grid.initialize(wxmin,wymin,wzmin,wxmax,wymax,wzmax,h,h);
	particles=new Particle_t[ptotalNum];
	pPos=new float[ptotalNum*3];
	for(int i=0;i<pxnum;i++){
		for(int j=0;j<pynum;j++){
			for(int k=0;k<pznum;k++){
				Particle_t t=new Particle();
				t->px=pstartx+i*pdist;t->py=pstarty+j*pdist;t->pz=pstartz+k*pdist;
				t->fx=0;t->fy=-gravity;t->fz=0;
				t->vx=t->vy=t->vz=0;	
				t->rho=t->lambda=0;
				int index=i*pynum*pznum+j*pznum+k;
				t->id=index;
				particles[index]=t;	
				pPos[3*index]=t->px;pPos[3*index+1]=t->py;pPos[3*index+2]=t->pz;
				grid.add(t);
			}
		}
	}	
}

void predictPos(){
    for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {
				particles[j]->vx += particles[j]->fx*tstep;
				particles[j]->vy += particles[j]->fy*tstep;
				particles[j]->vz += particles[j]->fz*tstep;

				particles[j]->fx = 0; particles[j]->fy = -gravity; particles[j]->fz = 0;

				particles[j]->px = pPos[3 * j] + particles[j]->vx*tstep;
				particles[j]->py = pPos[3 * j + 1] + particles[j]->vy*tstep;
				particles[j]->pz = pPos[3 * j + 2] + particles[j]->vz*tstep;
			}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}

void updateVelocityPos(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {


		particles[j]->vx = (particles[j]->px - pPos[3 * j]) / tstep;
		particles[j]->vy = (particles[j]->py - pPos[3 * j + 1]) / tstep;
		particles[j]->vz = (particles[j]->pz - pPos[3 * j + 2]) / tstep;


		std::vector<Indices> neighbors;
		neighbors.reserve(100);
		grid.getNeighbors(particles[j], neighbors);
		for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
			int ix = n->ix; int iy = n->iy; int iz = n->iz;
			std::set<Particle_t>* g = grid.get(ix, iy, iz);
			for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
				int k = (*p)->id;
				float w = calcW(j, k);
				float vijx = particles[k]->vx - particles[j]->vx;
				float vijy = particles[k]->vy - particles[j]->vy;
				float vijz = particles[k]->vz - particles[j]->vz;
				float c = 0.001;
				particles[j]->vx += c*vijx*w;
				particles[j]->vy += c*vijy*w;
				particles[j]->vz += c*vijz*w;
			}
		}



		pPos[3 * j] = particles[j]->px;
		pPos[3 * j + 1] = particles[j]->py;
		pPos[3 * j + 2] = particles[j]->pz;
	}
		}));
	}
		for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}

void update(){
	predictPos();
	computeRho();
	for (int i = 0; i<3; i++){
		computeLambda();
		computeDeltap();
		detect_collision();
	}
	updateVelocityPos();
	
	//computeVorticity();
	grid.adjust();
	generateBubbles();
	updateBubbles();
}

float Clamp(float v,float min,float max){
	if(v>max)
		v=max;
	if(v<min)
		v=min;
	return (v-min)/(max-min);
}

void updateBubbles(){
	std::list<Bubble>::iterator iter=bubbles.begin();
	while(iter!=bubbles.end()){
		
		if(iter->lifetime<0){
			std::list<Bubble>::iterator tmp=iter;
			iter++;
			bubbles.erase(tmp);
			continue;
		}
		
		if(iter->type==SPRAY){
			iter->vy+=-gravity*tstep;
			iter->x+=iter->vx*tstep;
			iter->y+=iter->vy*tstep;
			iter->z+=iter->vz*tstep;
		}else if(iter->type==FOAM){
			glm::vec3 v=grid.getVelocity(iter->x,iter->y,iter->z);
			iter->x+=v.x*tstep;
			iter->y+=v.y*tstep;
			iter->z+=v.z*tstep;
		}else{
			glm::vec3 v=grid.getVelocity(iter->x,iter->y,iter->z);
			iter->vx+=kd*(v.x-iter->vx);
			iter->vy+=kb*gravity*tstep+kd*(v.y-iter->vy);
			iter->vz+=kd*(v.z-iter->vz);
			iter->x+=iter->vx*tstep;
			iter->y+=iter->vy*tstep;
			iter->z+=iter->vz*tstep;
		}
		
		if(iter->x<wxmin){
			iter->x=wxmin;
			iter->vx=-iter->vx;
		}
		if(iter->x>wxmax){
			iter->x=wxmax;
			iter->vx=-iter->vx;
		}
		if(iter->y<wymin){
			iter->y=wymin;
			iter->vy=-iter->vy;
		}
		if(iter->y>wymax){
			iter->y=wymax;
			iter->vy=-iter->vy;
		}
		if(iter->z<wzmin){
			iter->z=wzmin;
			iter->vz=-iter->vz;
		}
		if(iter->z>wzmax){
			iter->z=wzmax;
			iter->vz=-iter->vz;
		}
		
		int density=grid.getDensity(iter->x,iter->y,iter->z);
		if(density<LIMIT1)
			iter->type=SPRAY;
		else if(density<LIMIT2)
			iter->type=FOAM;
		iter->lifetime-=tstep;
		iter++;
	}
}

void generateBubbles(){
	for(int i=0;i<ptotalNum;i++){
		
		float vel=sqrt(pow(particles[i]->vx,2)+pow(particles[i]->vy,2)+pow(particles[i]->vz,2));
		if(vel==0.0f)
			continue;
		float delta=vel*tstep;

		float ita=0.0f;
		glm::vec3 ni=grid.getNormal(particles[i]);
		
		std::vector<Indices> neighbors;
		neighbors.reserve(100);
		grid.getNeighbors(particles[i], neighbors);
		for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
			int ix = n->ix; int iy = n->iy; int iz = n->iz;
			std::set<Particle_t>* g = grid.get(ix, iy, iz);
			for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
				if((*p)->id==i){
					continue;
				}
				float dx=particles[i]->px-(*p)->px;
				float dy=particles[i]->py-(*p)->py;
				float dz=particles[i]->pz-(*p)->pz;
				float r = sqrt(pow(dx, 2)+pow(dy,2)+pow(dz,2));
				if(r==0)
					continue;
				if (h < r)
					continue;
				// for trapped air
				float dvx=particles[i]->vx-(*p)->vx;
				float dvy=particles[i]->vy-(*p)->vy;
				float dvz=particles[i]->vz-(*p)->vz;
				float dv=sqrt(pow(dvx,2)+pow(dvy,2)+pow(dvz,2));
				ita+=dv*(1-(dx/r*dvx/dv+dy/r*dvy/dv+dz/r*dvz/dv))*(1-r/h);
			}
		}
		
		
		ita=Clamp(ita,minta,maxta);
		float ik=0.5*mass*(pow(particles[i]->vx,2)+pow(particles[i]->vy,2)+pow(particles[i]->vz,2));
		ik=Clamp(ik,mink,maxk);
		int nd=ik*kta*ita*tstep;
		
		// generate bubbles
		for(int j=0;j<nd;j++){
			Bubble m;
			
			m.radius=randFloat(rmin,rmax);
			
			m.x=randFloat(particles[i]->px-delta,particles[i]->px+tstep);
			m.y=randFloat(particles[i]->py-delta,particles[i]->py+tstep);
			m.z=randFloat(particles[i]->pz-delta,particles[i]->pz+tstep);
			
			m.vx=randFloat(0,particles[i]->vx);
			m.vy=randFloat(0,particles[i]->vy);
			m.vz=randFloat(0,particles[i]->vz);
			
			m.lifetime=ik*maxlifetime;
			m.type=BUBBLE;
			bubbles.push_back(m);
		}
	}
}


void cleanup(){
	for(int i=0;i<ptotalNum;i++){
		delete[] particles[i];
	}
	delete[] particles;
	delete[] pPos;
}

void detect_collision(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {
				if (particles[j]->px < wxmin) {
					particles[j]->px = wxmin;
					particles[j]->vx = -particles[j]->vx;
				}
				if (particles[j]->px> wxmax) {
					particles[j]->px = wxmax;
					particles[j]->vx = -particles[j]->vx;
				}
				if (particles[j]->py < wymin) {
					particles[j]->py = wymin;
					particles[j]->vy = -particles[j]->vy;
				}
				if (particles[j]->py > wymax) {
					particles[j]->py = wymax ;
					particles[j]->vy = -particles[j]->vy;
				}
				if (particles[j]->pz < wzmin) {
					particles[j]->pz = wzmin;
					particles[j]->vz = -particles[j]->vz;
				}
				if (particles[j]->pz > wzmax) {
					particles[j]->pz = wzmax;
					particles[j]->vz = -particles[j]->vz;
				}
			
		}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}


void computeRho(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {
				particles[j]->rho = 0;
				std::vector<Indices> neighbors;
				neighbors.reserve(100);
				grid.getNeighbors(particles[j], neighbors);
				for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
					int ix = n->ix; int iy = n->iy; int iz = n->iz;
					std::set<Particle_t>* g = grid.get(ix, iy, iz);
					for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
						particles[j]->rho += pmass*calcW(j, (*p)->id);

					}
				}
			}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}

float getC(int i){
	return particles[i]->rho / prho0 - 1;
}

void computeLambda(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {
				float gradCSqr = 0;
				std::vector<Indices> neighbors;
				neighbors.reserve(100);
				grid.getNeighbors(particles[j], neighbors);
				glm::vec3 sum;
				sum.x = sum.y = sum.z = 0;
				for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
					int ix = n->ix; int iy = n->iy; int iz = n->iz;
					std::set<Particle_t>* g = grid.get(ix, iy, iz);
					for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
						int k = (*p)->id;
						Triple t = calcDW(j, k);
						gradCSqr += (pow(t.x, 2) + pow(t.y, 2) + pow(t.z, 2));
						sum.x += t.x; sum.y += t.y; sum.z += t.z;
					}
				}
				gradCSqr += (pow(sum.x, 2) + pow(sum.y, 2) + pow(sum.z, 2));
				gradCSqr*=pow(pmass/prho0,2);
				particles[j]->lambda = -getC(j) / (gradCSqr + sigma);
			}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}


void computeDeltap(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {
				glm::vec3 sum;
				sum.x = sum.y = sum.z = 0;
				std::vector<Indices> neighbors;
				neighbors.reserve(100);
				grid.getNeighbors(particles[j], neighbors);
				for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
					int ix = n->ix; int iy = n->iy; int iz = n->iz;
					std::set<Particle_t>* g = grid.get(ix, iy, iz);
					for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
						int k = (*p)->id;
						float scorr = -0.01 * pow(calcW(j, k) / scorrk, 4);
						Triple t = calcDW(k, j);
						float c = particles[j]->lambda + particles[k]->lambda + scorr;
						t.x *= c; t.y *= c; t.z *= c;
						sum.x += t.x; sum.y += t.y; sum.z += t.z;
					}
				}
				float t=pmass/prho0;
				sum.x *= t; sum.y *=t; sum.z *=t;
				particles[j]->px += sum.x; particles[j]->py += sum.y; particles[j]->pz += sum.z;
			}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();
}

Triple getOmega(int i,int j){
	Triple result;
	float vijx=particles[j]->vx-particles[i]->vx;
	float vijy=particles[j]->vy-particles[i]->vy;
	float vijz=particles[j]->vz-particles[i]->vz;
	float dx=particles[i]->px-particles[j]->px;float dy=particles[i]->py-particles[j]->py;
	float dz=particles[i]->pz-particles[j]->pz;
	float r=sqrt(pow(dx,2)+pow(dy,2)+pow(dz,2));
	r+=1e-5;
	float c=dwk*3*pow(h-r,2)/r;
	result.x=vijy*dz-vijz*dy;
	result.y=vijz*dx-vijx*dz;
	result.z=vijx*dy-vijy*dx;
	result.x*=c;result.y*=c;result.z*=c;
	return result;
}

void computeVorticity(){
	for(int i=0;i<threadNum;i++){
		threads.push_back(std::thread([&,i](){
			for (int j = i*workload; j < (i+1)*workload&&j<ptotalNum; j++) {


		glm::vec3 omega, newOmegax, newOmegay, newOmegaz;
		omega.x = omega.y = omega.z = 0.0f;
		newOmegax.x = newOmegax.y = newOmegax.z = 0.0f;
		newOmegay.x = newOmegay.y = newOmegay.z = 0.0f;
		newOmegaz.x = newOmegaz.y = newOmegaz.z = 0.0f;
		std::vector<Indices> neighbors;
		neighbors.reserve(100);
		grid.getNeighbors(particles[j], neighbors);
		for (std::vector<Indices>::iterator n = neighbors.begin(); n != neighbors.end(); n++){
			int ix = n->ix; int iy = n->iy; int iz = n->iz;
			std::set<Particle_t>* g = grid.get(ix, iy, iz);
			for (std::set<Particle_t>::iterator p = g->begin(); p != g->end(); p++){
				int k = (*p)->id;
				Triple o = getOmega(j, k);
				omega.x += o.x; omega.y += o.y; omega.z += o.z;

				particles[j]->px += 0.1;
				Triple dox = getOmega(j, k);
				newOmegax.x += dox.x; newOmegax.y += dox.y; newOmegax.z += dox.z;
				particles[j]->px -= 0.1;

				particles[j]->py += 0.1;
				Triple doy = getOmega(j, k);
				newOmegay.x += doy.x; newOmegay.y += doy.y; newOmegay.z += doy.z;
				particles[j]->py -= 0.1;

				particles[j]->pz += 0.1;
				Triple doz = getOmega(j, k);
				newOmegaz.x += doz.x; newOmegaz.y += doz.y; newOmegaz.z += doz.z;
				particles[j]->pz -= 0.1;
			}
		}
		float normal = sqrt(pow(omega.x, 2) + pow(omega.y, 2) + pow(omega.z, 2));
		float newNormalx = sqrt(pow(newOmegax.x, 2) + pow(newOmegax.y, 2) + pow(newOmegax.z, 2));
		float newNormaly = sqrt(pow(newOmegay.x, 2) + pow(newOmegay.y, 2) + pow(newOmegay.z, 2));
		float newNormalz = sqrt(pow(newOmegaz.x, 2) + pow(newOmegaz.y, 2) + pow(newOmegaz.z, 2));
		float etax = newNormalx - normal; float etay = newNormaly - normal; float etaz = newNormalz - normal;
		float etaNormal = sqrt(pow(etax, 2) + pow(etay, 2) + pow(etaz, 2));
		etaNormal += 1e-5;
		float nx = etax / etaNormal; float ny = etay / etaNormal; float nz = etaz / etaNormal;

		float fx = ny*omega.z - nz*omega.y;
		float fy = nz*omega.x - nx*omega.z;
		float fz = nx*omega.y - ny*omega.x;

		float c = 0.00001f;
		particles[j]->fx += c*fx; particles[j]->fy += c*fy; particles[j]->fz += c*fz;
			}
		}));
	}
	for(int i=0;i<threadNum;i++)
		threads[i].join();
	threads.clear();			
}

float calcW(int i,int j){
	float dx=particles[j]->px-particles[i]->px;
	float dy=particles[j]->py-particles[i]->py;
	float dz=particles[j]->pz-particles[i]->pz;
	float rsqr = pow(dx, 2)+pow(dy,2)+pow(dz,2);
	if (hsqr < rsqr)
		return 0.0f;
	return wk * pow(hsqr - rsqr, 3);
}

Triple calcDW(int i,int j){
	Triple result;
	result.x=result.y=result.z=0;
	if (i == j)
		return result;

	float dx=particles[j]->px-particles[i]->px;
	float dy=particles[j]->py-particles[i]->py;
	float dz=particles[j]->pz-particles[i]->pz;
	float r = sqrt(pow(dx, 2)+pow(dy,2)+pow(dz,2));
	if (h < r)
		return result;
	float c=3 * dwk* pow(h - r, 2) / (r + 1e-10);
	result.x=dx*c;result.y=dy*c;result.z=dz*c;
	return result;
}

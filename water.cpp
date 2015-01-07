#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <list>


#include"sph.h"
#include"util.h"


// for walls
extern int margin;

extern std::list<Bubble> bubbles;

// for mouse
int mxOrigin=-1;

// for camera
float cangle=0.0f;
float cx0=0.25f;
float cx=cx0,cy=-0.4f,cradius=0.8,cz=cradius;
float cdeltaAngle=0;
float near=0.00001f,far=500.0f;

// for OpenGL
int width=600,height=600;
float ratio=(float)width/(float)height;
int iterNum=0;
GLuint vertexbuffer, colorbuffer, radiibuffer;
GLuint programID, MVPID,lightID;
glm::mat4 MVP, Projection, View;
GLFWwindow* window;
glm::vec3 lightdir=glm::normalize(glm::vec3(1,0,0));

// data for rendering
extern float *pPos;

typedef std::chrono::duration<int, std::ratio<1, 24>> frame_duration;

void bindUniformVariables(){
	glUniformMatrix4fv(MVPID, 1, GL_FALSE, &MVP[0][0]);
	glUniform3fv(lightID,1, &lightdir[0]);
}

void renderScene(void){
	glClear( GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glUseProgram(programID);
	bindUniformVariables();
	
	std::vector<glm::vec3> colors(ptotalNum);
	for(int i=0;i<ptotalNum;i++){
		colors[i]=glm::vec3(1.0,1.0,1.0);
	}
	std::vector<float> radii(ptotalNum);
	for(int i=0;i<ptotalNum;i++){
		radii[i]=0.002;
	}
	
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*ptotalNum*3, pPos, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*ptotalNum, &colors[0], GL_STATIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
	glEnableVertexAttribArray(2);
	glBindBuffer(GL_ARRAY_BUFFER, radiibuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float)*ptotalNum, &radii[0], GL_STATIC_DRAW);
	glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
	
	
	glDrawArrays(GL_POINTS, 0, ptotalNum);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	
	if(bubbles.size()>0){
		colors.resize(bubbles.size());
		std::vector<glm::vec3> positions;
		radii.resize(bubbles.size());
		int i=0;
		for(std::list<Bubble>::iterator iter=bubbles.begin();iter!=bubbles.end();iter++){
			if(iter->type==BUBBLE)
				colors[i]=glm::vec3(0.0f,1.0f,0.0f);
			else if(iter->type==SPRAY)
				colors[i]=glm::vec3(1.0f,1.0f,0.0f);
			else
				colors[i]=glm::vec3(0.0f,0.0f,1.0f);
			positions.push_back(glm::vec3(iter->x,iter->y,iter->z));
			radii[i]=iter->radius;
			i++;
		}
	
		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*positions.size(), &positions[0], GL_STATIC_DRAW);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
		glEnableVertexAttribArray(1);
		glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(glm::vec3)*colors.size(), &colors[0], GL_STATIC_DRAW);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	
		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, radiibuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float)*radii.size(), &radii[0], GL_STATIC_DRAW);
		glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);
		
		glDrawArrays(GL_POINTS, 0, positions.size());
		glDisableVertexAttribArray(0);
		glDisableVertexAttribArray(1);
	}
}


void keyPress(GLFWwindow* window, int key, int scancode, int action, int mods){
	if(action!=	GLFW_PRESS&&action!=GLFW_REPEAT)
		return;
	switch (key){
		case GLFW_KEY_LEFT:
		case GLFW_KEY_A:
			cangle-=0.05f;
			cx = cx0+cradius*sin(cangle);
			cz = cradius*cos(cangle);
			break;
		case GLFW_KEY_RIGHT:
		case GLFW_KEY_D:
			cangle+=0.05f;
			cx = cx0+cradius*sin(cangle);
			cz = cradius*cos(cangle);
			break;
		case GLFW_KEY_UP:
		case GLFW_KEY_W:
			cy += 0.1f;
			break;
		case GLFW_KEY_DOWN:
		case GLFW_KEY_S:
			cy -= 0.1f;
			break;
	}
	Projection=glm::perspective(90.0f,ratio,near,far);
	View=glm::lookAt(glm::vec3(cx,cy,cz), glm::vec3(cx,cy,0), glm::vec3(0,1,0));
	MVP = Projection * View;
}

void computeMatricesFromInputs(){
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	float ratio=1.0*width/height;
	Projection=glm::perspective(90.0f,ratio,near,far);
	View=glm::lookAt(glm::vec3(cx,cy,cz), glm::vec3(cx,cy,0), glm::vec3(0,1,0));
	MVP = Projection * View;
}


int main( void )
{
	if( !glfwInit() )
	{
		fprintf( stderr, "Failed to initialize GLFW\n" );
		return -1;
	}
	
	glfwWindowHint(GLFW_SAMPLES, 1);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	
	window = glfwCreateWindow( width, height, "Water - Pu", NULL, NULL);
	if( window == NULL ){
		fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	
	glewExperimental=GL_TRUE;
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	
	start();
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_PROGRAM_POINT_SIZE);
	glClearColor(0.9f, 0.9f, 0.9f, 0.0f);

	GLuint VertexArrayID;
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);
	
	
	glGenBuffers(1, &vertexbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glPointSize(10.0);
	
	glGenBuffers(1, &colorbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, colorbuffer);
	
	glGenBuffers(1, &radiibuffer);
	glBindBuffer(GL_ARRAY_BUFFER, radiibuffer);
	
	
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	glfwSetKeyCallback(window, keyPress);
	
	programID = LoadShaders( "partShader.vertex","partShader.fragment" );
	glUseProgram(programID);
	MVPID = glGetUniformLocation(programID, "MVP");
	lightID  = glGetUniformLocation(programID, "lightdir");

	do{
		auto start_time = std::chrono::steady_clock::now();
		auto end_time = start_time + frame_duration(1);
		computeMatricesFromInputs();
		update();
		renderScene();
		glfwSwapBuffers(window);
		glfwPollEvents();
		std::this_thread::sleep_until(end_time);
	}while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
		   glfwWindowShouldClose(window) == 0 );

	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1,&colorbuffer);
	glDeleteBuffers(1,&radiibuffer);
	glDeleteVertexArrays(1, &VertexArrayID);
	glDeleteProgram(programID);
	
	glfwTerminate();

	return 0;
}

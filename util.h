#ifndef __UTIL_H_
#define __UTIL_H_

#include"water.h"
#include <fstream>
#include <vector>

GLuint LoadShaders(const char * vertex_file_path,const char * fragment_file_path);
float randFloat(float s,float e);
float randFloat();

#endif

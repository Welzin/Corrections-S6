// --------------------------------------------------------------------------
// gMini,
// a minimal Glut/OpenGL app to extend                              
//
// Copyright(C) 2007-2009                
// Tamy Boubekeur
//                                                                            
// All rights reserved.                                                       
//                                                                            
// This program is free software; you can redistribute it and/or modify       
// it under the terms of the GNU General Public License as published by       
// the Free Software Foundation; either version 2 of the License, or          
// (at your option) any later version.                                        
//                                                                            
// This program is distributed in the hope that it will be useful,            
// but WITHOUT ANY WARRANTY; without even the implied warranty of             
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)           
// for more details.                                                          
//                                                                          
// --------------------------------------------------------------------------

#include "Shader.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <fcntl.h>

#include <sys/types.h>
#include <unistd.h>

#ifdef _WIN32 /*[*/
#include <io.h>
#endif /*]*/

using namespace std;

#define printOpenGLError() printOglError(__FILE__, __LINE__)

/// Returns 1 if an OpenGL error occurred, 0 otherwise.
static int printOglError (const char * file, int line) {
  GLenum glErr;
  int    retCode = 0;
  glErr = glGetError ();
  while (glErr != GL_NO_ERROR) {
    printf ("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
    retCode = 1;
    glErr = glGetError ();
  }
  return retCode;
}


Shader::Shader () : shaderProgram (0), vertexShader (0), fragmentShader (0),
                    vertexShaderSize (0), fragmentShaderSize (0) {}


Shader::~Shader () {
  if (hasVertexShader ())
    glDeleteShader (vertexShader);
  if (hasFragmentShader ())
    glDeleteShader (fragmentShaderSize);
  glDeleteProgram (shaderProgram);
}


void Shader::loadFromFile (const string & vertexShaderFilename, const string & fragmentShaderFilename) {
  const GLchar * vertexShaderSource;
  const GLchar * fragmentShaderSource;
  if (vertexShaderFilename != "")
    vertexShaderSource = readShaderSource (vertexShaderFilename, vertexShaderSize);
  if (fragmentShaderFilename != "")
    fragmentShaderSource = readShaderSource (fragmentShaderFilename, fragmentShaderSize);

  shaderProgram = glCreateProgram ();
  if (hasVertexShader () == true) {
    compileAttach (vertexShader, GL_VERTEX_SHADER, &vertexShaderSource);
    delete [] vertexShaderSource;
  }
  if (hasFragmentShader () == true) {
    compileAttach (fragmentShader, GL_FRAGMENT_SHADER, &fragmentShaderSource);
    delete [] fragmentShaderSource;
  }
  glLinkProgram (shaderProgram);
  printOpenGLError ();
	
  GLint linked;
  glGetProgramiv (shaderProgram, GL_LINK_STATUS, &linked);
  printProgramInfoLog (shaderProgram);
  if (!linked)
    throw ShaderException ("Error: Shaders not linked");
}


void Shader::bind () {
  glUseProgram (shaderProgram);
}


void Shader::unbind () {
  glUseProgram (0);
}


// ------------------
// Protected methods.
// ------------------

unsigned int Shader::getShaderSize (const string & filename) {
  int fd;
  unsigned int count = 0;
#ifdef _WIN32
  fd = _open (filename.c_str (), _O_RDONLY);
  if (fd != -1) {
    count = static_cast<unsigned int>(_lseek (fd, 0, SEEK_END) + 1);
    _close(fd);
  } else 
    throw ShaderException (string ("getShaderSize: bad Shader File Name") + filename); 
#else 
  fd = open (filename.c_str (), O_RDONLY);
  if (fd != -1) {
    count = static_cast<unsigned int>(lseek (fd, 0, SEEK_END) + 1);
    close(fd);
  } else 
    throw ShaderException (string ("getShaderSize: bad Shader File Name") + filename); 
#endif
  return count;
}


GLchar * Shader::readShaderSource(const string & shaderFilename, unsigned int & shaderSize) {
  shaderSize = getShaderSize (shaderFilename);
  FILE * fh = fopen (shaderFilename.c_str (), "r");
  if (!fh)
    throw ShaderException (string ("readShaderSource : Bad Shader Filename ") + shaderFilename);
    
  GLchar * shaderSource = new GLchar[shaderSize];
  fseek (fh, 0, SEEK_SET);
  int count = fread (shaderSource, 1, shaderSize, fh);
  shaderSource[count] = '\0';
  fclose (fh);
  //  if (ferror (fh))
  //  throw ShaderException (string ("readShaderSource : error while reading file") + shaderFilename);

  return shaderSource;
}


GLint Shader::getUniLoc (GLuint program, const GLchar *name){
  GLint loc = glGetUniformLocation (program, name);
  //if (loc == -1)
    //throw ShaderException (string ("No such uniform named ") + string (name));
  printOpenGLError(); 
  return loc;
}


void Shader::printShaderInfoLog (GLuint shader) {
  int infologLength = 0;
  int charsWritten  = 0;
  GLchar *infoLog;
  
  printOpenGLError ();  // Check for OpenGL errors
  glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &infologLength);
  printOpenGLError ();  // Check for OpenGL errors

  if (infologLength > 0) {
    infoLog = new GLchar[infologLength];
    if (infoLog == NULL) {
      printf("ERROR: Could not allocate InfoLog buffer\n");
      exit(1);
    }
    glGetShaderInfoLog (shader, infologLength, &charsWritten, infoLog);
    //    cerr << "InfoLog:" << endl << infoLog << endl << endl;
    delete [] infoLog;
  }
  printOpenGLError();  // Check for OpenGL errors
}

void Shader::printProgramInfoLog (GLuint program) {
  int infologLength = 0;
  int charsWritten  = 0;
  GLchar *infoLog;
  
  printOpenGLError ();  // Check for OpenGL errors
  glGetProgramiv (program, GL_INFO_LOG_LENGTH, &infologLength);
  printOpenGLError ();  // Check for OpenGL errors

  if (infologLength > 0) {
    infoLog = new GLchar[infologLength];
    if (infoLog == NULL) {
      printf("ERROR: Could not allocate InfoLog buffer\n");
      exit(1);
    }
    glGetProgramInfoLog (program, infologLength, &charsWritten, infoLog);
    //cerr << "InfoLog:" << endl << infoLog << endl << endl;
    delete [] infoLog;
  }
  printOpenGLError();  // Check for OpenGL errors
}

void Shader::compileAttach (GLuint & shader, GLenum type, const GLchar ** source) {
  GLint shaderCompiled;
  shader = glCreateShader (type);
  glShaderSource (shader, 1, source, NULL);
  glCompileShader (shader);
  printOpenGLError ();  // Check for OpenGL errors
  glGetShaderiv (shader, GL_COMPILE_STATUS, &shaderCompiled);
  printOpenGLError ();  // Check for OpenGL errors
  printShaderInfoLog (shader);
  if (!shaderCompiled)
    throw ShaderException ("Error: shader not compiled");
  glAttachShader (shaderProgram, shader);
}

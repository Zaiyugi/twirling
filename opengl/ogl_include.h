/* Zachary Shore
 * 2014-01-29
 * Includes for OpenGL and Camera
 */

#ifndef __OGL_INCLUDE_H__
#define __OGL_INCLUDE_H__

#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include "OGL_Matrix.h"
#include "OGL_Camera.h"

#endif

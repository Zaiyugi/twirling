/* Zachary Shore
 * 2014-01-03
 * OpenGL Functions and Callbacks
 */

#ifndef __OGL_H__
#define __OGL_H__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

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

#define PI180 M_PI/180.0 

using namespace std;

int windowWidth, windowHeight;

/* ---------------------------------- */
bool update_display = true;
int update_timestep = 25;
int drawMode = 0;
int mouseInWindow = 0;

OGL_Camera* cam;
/* ---------------------------------- */

void setWindowSize(int w, int h)
{ 
   windowWidth = w;
   windowHeight = h;
}

void setCamera(Vector3d eye, Vector3d view, Vector3d up)
{
   cam = new Camera(eye, view up);
}

void cleanup() { delete cam; }

void IdleCB(){}

void timerCB(int millisec)
{
   glutTimerFunc(millisec, timerCB, millisec);

   glutPostRedisplay();
}

void reshapeCB(int w, int h)
{
   windowWidth = w;
   windowHeight = h;
   glViewport(0, 0, (GLsizei)w, (GLsizei)h);

   //float aspectRatio = (float)w / h;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glMatrixMode(GL_MODELVIEW);
}

void keyboardCB(unsigned char key, int x, int y)
{
    switch(key)
    {
    case 27: // ESCAPE
        exit(0);
        break;

    case 'r': // switch rendering modes (fill -> wire -> point)
    case 'R':
        drawMode = ++drawMode % 3;
        if(drawMode == 0)        // fill mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glEnable(GL_DEPTH_TEST);
        }
        else if(drawMode == 1)  // wireframe mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDisable(GL_DEPTH_TEST);
        }
        else                    // point mode
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
            glDisable(GL_DEPTH_TEST);
        }
        break;

    default:
        ;
    }

    glutPostRedisplay();
}

void mouseEntryCB(int state)
{
   if(state == GLUT_ENTERED)
   {
      /*printf("------------------------------\n");
      printf("!! MOUSE ENTERED WINDOW %d !!\n",state);
      printf("------------------------------\n");*/
      mouseInWindow = state;
   } else if(state == GLUT_LEFT) {
      /*printf("------------------------------\n");
      printf("!! MOUSE LEFT WINDOW %d !!\n", state);
      printf("------------------------------\n");*/
      mouseInWindow = state;
   }

}

void mouseMotionCB(int x, int y)
{
   cam->HandleMouseMotion(x, y);
   glutPostRedisplay();

}

void mouseCB(int button, int state, int x, int y)
{
   cam->HandleMouseEvent(button, state, x, y);
}

int initGLUT(int argc, char **argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
   glutInitWindowSize(windowWidth, windowHeight);
   glutInitWindowPosition(100, 100);
   int handle = glutCreateWindow(argv[0]);

   glutDisplayFunc(draw);
   glutReshapeFunc(reshapeCB);
   glutKeyboardFunc(keyboardCB);
   glutMouseFunc(mouseCB);
   glutMotionFunc(mouseMotionCB);
   glutEntryFunc(mouseEntryCB);
   glutIdleFunc(IdleCB);
   glutTimerFunc(update_timestep, timerCB, update_timestep);

return handle;
}

void initGL()
{
   glShadeModel(GL_SMOOTH);

   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
   glEnable(GL_DEPTH_TEST);

   //glEnable(GL_BLEND);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   glClearColor(0, 0, 0, 0);                   // background color
   glClearDepth(1.0f);                         // 0 is near, 1 is far
   glDepthFunc(GL_LEQUAL);

   initLights();
}

void initLights()
{
   glLightModelf(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

   GLfloat lightKa[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat lightKd[] = {0.7f, 0.7f, 0.7f, 1.0f};  // diffuse light
   GLfloat lightKs[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
   glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);

   float lightPos[4] = {0, 5, 0, 1}; // positional light
   float lightDir[4] = {0, -1, 0};
   glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
   glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, lightDir);

   float spotVars[2] = {90, 180};
   glLightfv(GL_LIGHT0, GL_SPOT_CUTOFF, &spotVars[0]);
   //glLightfv(GL_LIGHT0, GL_SPOT_EXPONENT, &spotVars[1]);

   //glEnable(GL_LIGHTING);
   //glEnable(GL_LIGHT0);
}

#endif

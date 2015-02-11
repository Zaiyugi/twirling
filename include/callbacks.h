/* Zachary Shore
 * 2014-01-29
 * CPSC-881-001
 * Lab 1: Curve Modeling
 * Callbacks for OpenGL
 */
#include "opengl/ogl_include.h"
#include <AntTweakBar.h>

int windowHandle;
int display_update = 10;
int windowWidth, windowHeight;
int windowPositionX, windowPositionY;

extern ogl::OGL_Camera *cam;

void draw();
void ogl_cleanup();

void IdleCB() {}

void timerCB(int ms)
{
   glutTimerFunc(ms, timerCB, ms);

   glutPostRedisplay();
}

void reshapeCB(int w, int h)
{
   windowWidth = w;
   windowHeight = h;
   glViewport(0, 0, (GLsizei)w, (GLsizei)h);

   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glMatrixMode(GL_MODELVIEW);

   TwWindowSize(windowWidth, windowHeight);
}

/*
void keyboardCB(unsigned char key, int x, int y)
{
   switch(key)
   {
      case 27:
         ogl_cleanup();
         exit(0);
         break;

      case 'v':
      case 'V':
         model_state = VIEWING;
         break;

      case 'e':
      case 'E':
         model_state = EDITING;
         break;

      default:
         ;
   }

   glutPostRedisplay();
}
*/

void mouseEntryCB(int state)
{
   if(state == GLUT_ENTERED)
      ;
   else if(state == GLUT_LEFT)
      ;
}

/*
void mouseMotionCB(int x, int y)
{
   if( !TWEventMouseMotionGLUT(x, y) )
   {
      cam->HandleMouseMotion(x, y);
      glutPostRedisplay();
   }
}
*/

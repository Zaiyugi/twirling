/* Name: Zachary Shore
 * Date: July 22, 2013
 * Desc: Torus Knot
 */
#include "global_prototypes.h"
#include "OGL_Matrix.h"
#include "OGL_Camera.h"

#include <cmath>
#define PI180 0.0174532925

using namespace std;

/* ---------------------------------- */
bool pauseP = true;
int drawMode = 0;
int focalPosX, focalPosY;
int mouseInWindow = 0;
int winWidth(850), winHeight(850);

int p = 3;
int q = 2;
int segments = 20;

Vector3d RIGHT(1, 0, 0);
Vector3d    UP(0, 1, 0);
Vector3d FRONT(0, 0, 1);
/* ---------------------------------- */

Camera *cam;

void draw()
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   cam->PerspectiveDisplay(winWidth, winHeight);

   glPushMatrix();

   // First material is always the wall color
   glPointSize(2.0);
/*
   glColor3f(1, 1, 1);
   glBegin(GL_LINES);
      glVertex3d(0, 0, 0);
      glVertex3d(5, 0, 0);

      glVertex3d(0, 0, 0);
      glVertex3d(0, 5, 0);

      glVertex3d(0, 0, 0);
      glVertex3d(0, 0, 5);
   glEnd();
*/

   double i = 0;
   double di = 1.0 / (double)segments;
   double phi = 0.0;
   double r, x, y, z;
   Vector3d f;

   glBegin(GL_LINE_LOOP);
   for(i; i < 1.0; i += di)
   {
      phi = 2.0 * M_PI * i;
      r = cos((double)q * phi) + 2;
      x = r * cos((double)p * phi);
      y = r * sin((double)p * phi);
      z = -sin((double)q * phi);

      f = x * RIGHT.normalize() + y * UP.normalize() + z * FRONT.normalize();

      glColor3f(1.0 - x / r, 1.0 - y / r, 1.0 - z / r);
      glVertex3d(f[0], f[1], f[2]);
   }
   glEnd();

   glPopMatrix();
   glutSwapBuffers();
}

int main(int argc, char **argv)
{
   if(argc == 4)
   {
      p = atoi(argv[1]);
      q = atoi(argv[2]);
      segments = atoi(argv[3]);
   } else {
      printf("Usage: TorusKnot <p> <q>\np and q should be relatively prime\n");
      exit(0);
   }

   // Begin OpenGL
   initGLUT(argc, argv);
   initGL();

   glutMainLoop();

return 0;
}

void setCamera()
{
   cam = new Camera(Vector3d(0, 0, 1), Vector3d(0, 0, 0), Vector3d(0, 1, 0));
}

void IdleCB()
{
}

void timerCB(int millisec)
{
   glutTimerFunc(millisec, timerCB, millisec);

   glutPostRedisplay();
}

void reshapeCB(int w, int h)
{
   winWidth = w;
   winHeight = h;
   glViewport(0, 0, (GLsizei)w, (GLsizei)h);

   float aspectRatio = (float)w / h;
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   //gluPerspective(60.0f, (float)(w)/h, 1.0f, 1000.0f);
   cam->PerspectiveDisplay(winWidth, winHeight);

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
   glutInitWindowSize(winWidth, winHeight);
   glutInitWindowPosition(100, 100);
   int handle = glutCreateWindow(argv[0]);

   glutDisplayFunc(draw);
   glutReshapeFunc(reshapeCB);
   glutKeyboardFunc(keyboardCB);
   glutMouseFunc(mouseCB);
   glutMotionFunc(mouseMotionCB);
   glutEntryFunc(mouseEntryCB);
   glutIdleFunc(IdleCB);
   glutTimerFunc(25, timerCB, 25);

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
   setCamera();
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

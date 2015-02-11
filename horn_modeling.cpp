/* Zachary Shore
 * Created: 2014-03-31
 * Last Edit: 2014-03-31
 * CPSC-881-001
 * Final: Twisting Scupltures
 */
#include <cstdlib>
#include <cstdio>
#include <cmath> 
#include <string>
#include <vector>
#include <algorithm>

#include "opengl/ogl_include.h"
#include "callbacks.h"
#include <AntTweakBar.h>

#include "Vector.h"
#include "Curves.h"
#include "Surfaces.h"

using namespace std;
using namespace ogl;

const static unsigned int VIEW = 0;
const static unsigned int EDIT = 1;
const static unsigned int CREATE = 2;
const static unsigned int SELECT = 3;

void initGLUT(int argc, char **argv);
void initGL();

ogl::OGL_Camera *cam;
TwBar *main_bar;
TwBar *horn_bar;
int model_state = 0;

lux::Vector selected_point;

/* Used in projecting window coords to object coords */
/* Used for creating control polygon */
GLint viewport[4];
GLdouble modelview[16];
GLdouble projection[16];
GLfloat winX, winY, winZ;
GLdouble posX, posY, posZ;
/* ------------------------- */

/* ------------------------- */
/*
 * Curve Variables
 * Types:
 *		Bezier_Sampled
 *		Bezier_Subdiv
 * 	Uniform_3D_bspline
 */
//int samples = 30;
//int subdivs = 5;
int spline_order = 0;
bool showControlPolygons = true;
int curve_type = cagd::curves::BEZIER_SAMPLED;

std::vector<lux::Vector> global_curve;
std::vector<lux::Vector> control_polygon;
std::vector<lux::Vector> point_glyph;
//std::vector<cagd::curves::Frame> local_coord_frame;
/* ------------------------- */

/* ------------------------- */
/* Surface Variables */
/* Types:
 * 	Extrusion
 * 	Revolution
 * 	Sweep
 * 	Loft
 */

cagd::surf::param_t surface_params;
//int surface_type = cagd::surf::EXTRUSION;
double surface_axis[3]{0, 0, 1};
double surface_axis_orig[3]{0, 0, 0};

cagd::surf::Surface* global_surface = nullptr;

/* ------------------------- */

int selected_surface = -1;
int selected_slice = -1;
int selected_ndx = -1;
std::string IO_Surface_Name = "models/tet";

int sweep_select_slice = -1;

lux::Vector OGLToLux(ogl::Vector3d V)
{ return lux::Vector(V[0], V[1], V[2]); }

/* ------------------------- */

void updateSelectedCurveProperties()
{
   if(spline_order != surface_params.order)
   {
      if(surface_params.closed)
      {
         for(int i = 0; i < spline_order; i++)
            control_polygon.pop_back();
         for(int i = 0; i < surface_params.order; i++)
            control_polygon.push_back(control_polygon[i]);
      }
   }
   spline_order = surface_params.order;
}

void regenerateCurve()
{
   global_curve.clear();

   if(curve_type == cagd::curves::BEZIER_SAMPLED)
   {
      if(surface_params.samples > 0)
      {
         cagd::curves::deCasteljauSampling(
            control_polygon,
            global_curve,
            surface_params.samples
         );
      }

   } else if(curve_type == cagd::curves::UNIFORM_3D_BSPLINE) {
      if(surface_params.samples > 0 && surface_params.order > 0)
      {
         cagd::curves::BSplineUniformKnotsArbitraryOrder(
            control_polygon,
            global_curve,
            surface_params.order,
            surface_params.samples
         );
      }

   }
/*
   if(global_curve.size() > 0)
   {
      // Build Rotation Frames for Curve
      cagd::curves::Frame U_0;
      U_0.r = lux::Vector(1, 0, 0);
      U_0.s = lux::Vector(0, 1, 0);
      U_0.t = lux::Vector(0, 0, 1);

      lux::Vector target = (global_curve[1] - global_curve[0]).unitvector();
      lux::Vector rot_axis = target ^ U_0.t;
      double rot_angle = std::acos(U_0.t * target);
      U_0.transform(cagd::surf::rot_mat(rot_axis, rot_angle));

      cagd::curves::generateFrame(global_curve, U_0, local_coord_frame);
   }*/

}

lux::Vector unproject(int x, int y)
{
   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
   glGetDoublev( GL_PROJECTION_MATRIX, projection );
   glGetIntegerv( GL_VIEWPORT, viewport );

   winX = static_cast<float>(x);
   winY = static_cast<float>(viewport[3] - y);
   winZ = 1.0;

/*
   printf("\tUnproject\n");
   printf("\tX: %d | Y: %d\n", x, (int)winY);
   printf("\tWinX: %f | WinY: %f | WinZ: %f\n", winX, winY, winZ);
*/

   gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ );

   return lux::Vector(posX, posY, posZ);
}

void updateCurve(int x, int y)
{
   lux::Vector p_obj = unproject(x, y);

   lux::Vector c_aim(cam->Aim[0], cam->Aim[1], cam->Aim[2]);
   lux::Vector c_pos(cam->Pos[0], cam->Pos[1], cam->Pos[2]);
   lux::Vector c_dir = (c_aim - c_pos).unitvector();

   lux::Vector p_far = c_pos + c_dir * cam->FarPlane;
   double farToAim = (c_aim - p_far).magnitude();
   double farToPos = cam->FarPlane;
   double aimRatio = farToAim / farToPos;

   double objToPos = (c_pos - p_obj).magnitude();
   double objToAim = aimRatio * objToPos;

   lux::Vector pos = p_obj + (c_pos - p_obj).unitvector() * objToAim;

   control_polygon[selected_ndx] = pos;
   if(surface_params.closed)
   {
      int N_cp = control_polygon.size();
      if(selected_ndx < surface_params.order)
         control_polygon[N_cp - (surface_params.order - selected_ndx)] = pos;
      else if(selected_ndx > (N_cp - surface_params.order) - 1)
         control_polygon[selected_ndx - (N_cp - surface_params.order)] = pos;
   }

   regenerateCurve();
 
}

void keyboardCB(unsigned char key, int x, int y)
{
   int ret = TwEventKeyboardGLUT(key, x, y);

   if( ret )
   {
      if(control_polygon.size() > 0 && global_curve.size() > 0)
      {
         updateSelectedCurveProperties();
         regenerateCurve();
      }
   } else if( !ret ) {

      switch(key)
      {
         case 'v':
         case 'V':
            model_state = VIEW;
            break;

         case 'e':
         case 'E':
            model_state = EDIT;
            break;

         case 'c':
         case 'C':
            model_state = CREATE;
            break;

         case 's':
         case 'S':
            model_state = SELECT;
            break;

         default:
            ;
      }
   }

   glutPostRedisplay();
}

void TW_CALL GenerateCurve(void *)
{
   if(control_polygon.size() > 0)
   {
      regenerateCurve();
   }

}

void TW_CALL ClearCurve(void *)
{
   global_curve.clear();
   control_polygon.clear();
}

void TW_CALL ImportSurface(void *)
{
   if(global_surface != nullptr)
   {
      delete global_surface;
      global_surface = nullptr;
   }

   global_surface = new cagd::surf::Surface(surface_params);
   int ret = cagd::surf::ReadSurface(global_surface, IO_Surface_Name);

   global_surface->generate();
}

void TW_CALL ExportSurface(void *)
{
   if(global_surface != nullptr)
      int ret = cagd::surf::WriteSurface(global_surface, IO_Surface_Name);
}

void TW_CALL DeleteSurface(void *)
{
   if(global_surface != nullptr)
   {
      delete global_surface;
      global_surface = nullptr;
   }

}

void TW_CALL CloseCurve(void *)
{
   if(control_polygon.size() > 0)
   {
      if(surface_params.closed)
      {
         for(int i = 0; i < surface_params.order; i++)
            control_polygon.pop_back();
         surface_params.closed = false;
      } else {
         for(int i = 0; i < surface_params.order; i++)
            control_polygon.push_back(control_polygon[i]);
         surface_params.closed = true;
      }

      regenerateCurve();
   }
}

void TW_CALL DisplayControlPolygons(void *)
{
   showControlPolygons = !showControlPolygons;
}

void TW_CALL ClearSelection(void *)
{
   if(global_surface)
      global_surface->clearSelection();
}

void TW_CALL SelectAll(void *)
{
   if(global_surface)
      global_surface->selectAll();
}

void TW_CALL SelectAllTriangles(void *)
{
   if(global_surface)
      global_surface->selectAllTriangles();
}

void TW_CALL Undo(void *)
{
   if(global_surface)
      global_surface->Undo();
}

void TW_CALL AssignPath(void *)
{
   if(global_surface)
      if(global_curve.size() > 0)
         global_surface->assignPath(global_curve);
}

void TW_CALL AddSubdivision(void *)
{
   surface_params.axis = lux::Vector(surface_axis[0], surface_axis[1], surface_axis[2]).unitvector();
   surface_params.orig = lux::Vector(surface_axis_orig[0], surface_axis_orig[1], surface_axis_orig[2]);

   if(global_surface)
      global_surface->AddSubdivision(surface_params);
}

void TW_CALL AddExtrusion(void *)
{
   surface_params.axis = lux::Vector(surface_axis[0], surface_axis[1], surface_axis[2]).unitvector();
   surface_params.orig = lux::Vector(surface_axis_orig[0], surface_axis_orig[1], surface_axis_orig[2]);

   if(global_surface)
      global_surface->AddExtrusion(surface_params);
}

void TW_CALL ExitProgram(void *)
{
   ogl_cleanup();
   exit(0);
}

void TW_CALL SetMyStdStringCB(const void *value, void *)
{
   const std::string *srcPtr = static_cast<const std::string *>(value);
   IO_Surface_Name = *srcPtr;
}

void TW_CALL GetMyStdStringCB(void *value, void *)
{
   std::string *destPtr = static_cast<std::string *>(value);
   TwCopyStdStringToLibrary(*destPtr, IO_Surface_Name);
}

void mouseMotionCB(int x, int y)
{
   int ret = TwEventMouseMotionGLUT(x, y);
   if( ret )
   {
      if(control_polygon.size() > 0 && global_curve.size() > 0)
      {
         updateSelectedCurveProperties();
         regenerateCurve();
      }
   } else if( !TwEventMouseMotionGLUT(x, y) ) {

      if( model_state == VIEW )
      {
         cam->HandleMouseMotion(x, y);
      } else if( model_state == EDIT ) {
         if(showControlPolygons)
         {
            if(control_polygon.size() > 0 && selected_ndx > -1)
               updateCurve(x, y);
         }

      }

   }

   glutPostRedisplay();
}

void mouseButtonCB(int button, int state, int x, int y)
{
   int ret = TwEventMouseButtonGLUT(button, state, x, y);
   if( ret )
   {
      if(control_polygon.size() > 0 && global_curve.size() > 0)
      {
         updateSelectedCurveProperties();
         regenerateCurve();
      }
   } else if( !ret ) {

      if( model_state == VIEW )
      {
         cam->HandleMouseEvent(button, state, x, y);
      } else if( model_state == EDIT && state == GLUT_DOWN && button == GLUT_LEFT_BUTTON ) {
         if(showControlPolygons)
         {
            lux::Vector p_far = unproject(x, y);

            lux::Vector r_pos = OGLToLux(cam->Pos);
            lux::Vector r_dir = (p_far - r_pos).unitvector();

            // First, test each surface
            double best_distance = (p_far - r_pos).magnitude();
            int closest_ndx = -1;

            int ndx;
            if(cagd::curves::testCurve(control_polygon, surface_params.closed, surface_params.order, r_pos, r_dir, ndx))
            {
               double curve_dist = (control_polygon[ndx] - r_pos).magnitude();
               if(curve_dist < best_distance)
               {
                  closest_ndx = ndx;
               }
            }

            selected_ndx = closest_ndx;

         }

      } else if( model_state == CREATE && state == GLUT_DOWN && button == GLUT_LEFT_BUTTON ) {

         lux::Vector p_obj = unproject(x, y);

         lux::Vector c_aim(cam->Aim[0], cam->Aim[1], cam->Aim[2]);
         lux::Vector c_pos(cam->Pos[0], cam->Pos[1], cam->Pos[2]);
         lux::Vector c_dir = (c_aim - c_pos).unitvector();

         lux::Vector p_far = c_pos + c_dir * cam->FarPlane;
         double farToAim = (c_aim - p_far).magnitude();
         double farToPos = cam->FarPlane;
         double aimRatio = farToAim / farToPos;

         double objToPos = (c_pos - p_obj).magnitude();
         double objToAim = aimRatio * objToPos;

         lux::Vector pos = p_obj + (c_pos - p_obj).unitvector() * objToAim;

         control_polygon.push_back(pos);

      } else if( model_state == SELECT && state == GLUT_DOWN && button == GLUT_LEFT_BUTTON ) {

         lux::Vector p_obj = unproject(x, y);

         lux::Vector c_aim(cam->Aim[0], cam->Aim[1], cam->Aim[2]);
         lux::Vector c_pos(cam->Pos[0], cam->Pos[1], cam->Pos[2]);
         lux::Vector c_dir = (c_aim - c_pos).unitvector();

         lux::Vector p_far = c_pos + c_dir * cam->FarPlane;
         double farToAim = (c_aim - p_far).magnitude();
         double farToPos = cam->FarPlane;
         double aimRatio = farToAim / farToPos;

         double objToPos = (c_pos - p_obj).magnitude();
         double objToAim = aimRatio * objToPos;

         lux::Vector pos = p_obj + (c_pos - p_obj).unitvector() * objToAim;

         if(global_surface)
            if(global_surface->select(cam, x, y, windowWidth, windowHeight))
               std::printf("Face Selection SUCCESS:\n");
      }
   }

   glutPostRedisplay();
}

void draw()
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   cam->PerspectiveDisplay(windowWidth, windowHeight);

   glPointSize(3.0);

   lux::Vector c_view = OGLToLux( (cam->Aim - cam->Pos).normalize() );
   lux::Vector c_up = OGLToLux( cam->Up ).unitvector();
   lux::Vector c_right = (c_view ^ c_up).unitvector();

   // Draw partial control polygons

   if(showControlPolygons)
   {
      for(int i = 0; i < static_cast<int>(control_polygon.size()); i++)
      {
         lux::Vector pt = control_polygon[i];

         if(i == selected_ndx)
            glColor4f(0.0, 0.0, 1.0, 1.0);
         else
            glColor4f(1.0, 1.0, 1.0, 1.0);
         glBegin(GL_LINE_LOOP);
            for(size_t j = 0; j < point_glyph.size(); j++)
            {
               lux::Vector X = pt
                             + c_right * point_glyph[j][0]
                             + c_up * point_glyph[j][1];

               glVertex3f(X[0], X[1], X[2]);
            }
         glEnd();

         glBegin(GL_POINTS);
            glNormal3f(0, 1, 0);
            glVertex3f(pt[0], pt[1], pt[2]);
         glEnd();
      
      }

      glColor4f(0.35, 0.35, 0.35, 1.0);
      glBegin(GL_LINE_STRIP);
      for_each(control_polygon.begin(), control_polygon.end(), [] (lux::Vector pt)
      {
         glNormal3f(0, 1, 0);
         glVertex3f(pt[0], pt[1], pt[2]);
      } );
      glEnd();

   }

   // Draw curves

   glColor4f(1.0, 1.0, 0.0, 1.0);
   glBegin(GL_LINE_STRIP);
      for_each(global_curve.begin(), global_curve.end(), [] (lux::Vector pt)
      {
         glVertex3f(pt[0], pt[1], pt[2]);
      } );
   glEnd();

   // Draw Curve Rotation Frame

/*   glBegin(GL_LINES);
      for(size_t i = 0; i < local_coord_frame.size(); i++)
      {
         lux::Vector orig = global_curve[i];

         // R
         lux::Vector offset = 0.25 * local_coord_frame[i].r;
         glColor4f(1.0, 0.0, 0.0, 1.0);
         glVertex3f(orig[0], orig[1], orig[2]);
         glVertex3f(
            orig[0] + offset[0], 
            orig[1] + offset[1], 
            orig[2] + offset[2]
         );

         // S
         offset = 0.25 * local_coord_frame[i].s;
         glColor4f(0.0, 1.0, 0.0, 1.0);
         glVertex3f(orig[0], orig[1], orig[2]);
         glVertex3f(
            orig[0] + offset[0], 
            orig[1] + offset[1], 
            orig[2] + offset[2]
         );

         // T
         offset = 0.25 * local_coord_frame[i].t;
         glColor4f(0.0, 0.0, 1.0, 1.0);
         glVertex3f(orig[0], orig[1], orig[2]);
         glVertex3f(
            orig[0] + offset[0], 
            orig[1] + offset[1], 
            orig[2] + offset[2]
         );
      }
   glEnd();*/

   // Draw surfaces
   if(global_surface != nullptr)
   {
      global_surface->draw();
   }

   TwDraw();

   glutSwapBuffers();
}

void ogl_cleanup()
{
   delete cam;
   TwTerminate();
}

int main(int argc, char **argv)
{
   windowWidth = 960;
   windowHeight = 540;
   surface_params.horn_params.hn = 0.5;
   surface_params.horn_params.sn = 0.5;
   surface_params.surf_type = 3;

   initGLUT(argc, argv);
   initGL();

   // Init AntTweakBar
   TwInit(TW_OPENGL, NULL);
   TwWindowSize(windowWidth, windowHeight);

   /* MAIN CONTROLS */
   main_bar = TwNewBar("General");
   TwDefine(" 'General' size='300 280' valueswidth=fit contained=true ");

   TwEnumVal modelingStatesEV[4] = { {0, "View"}, {1, "Edit"}, {2, "Create"}, {3, "Select"} };
   TwType modelStateType = TwDefineEnum("modelStateType", modelingStatesEV, 4);
   TwAddButton(main_bar, "Exit", ExitProgram, NULL, " group='Main Controls' ");
   TwAddVarRW(main_bar, "Current Mode", modelStateType, &model_state,
              " keyIncr='<' keyDecr='>' help='Select modeling mode' group='Main Controls' ");
   TwAddButton(main_bar, "Show/Hide Control Polytopes", 
               DisplayControlPolygons, NULL, 
               " group='Main Controls' ");
   TwAddButton(main_bar, "Clear Selection", ClearSelection, NULL, " group='Main Controls' ");
   TwAddButton(main_bar, "Select All", SelectAll, NULL, " group='Main Controls' ");
   TwAddButton(main_bar, "Select All Triangles", SelectAllTriangles, NULL, " group='Main Controls' ");

   TwAddButton(main_bar, "Undo last operation", Undo, NULL, " group='Main Controls' ");
   TwAddButton(main_bar, "Apply Extrusion", AddExtrusion, NULL, " group='Main Controls' ");
   TwAddButton(main_bar, "Apply Subdivision", AddSubdivision, NULL, " group='Main Controls' ");
   TwAddButton(main_bar, "Assign path to surface", AssignPath, NULL, " group='Main Controls' ");

   TwAddVarCB(main_bar, "Filename", TW_TYPE_STDSTRING, SetMyStdStringCB, 
              GetMyStdStringCB, NULL, " label='Surface Name' group='I/O Controls' ");
   TwAddButton(main_bar, "Import Surface", ImportSurface, NULL, " group='I/O Controls' ");
   TwAddButton(main_bar, "Export Surface", ExportSurface, NULL, " group='I/O Controls' ");
   TwAddButton(main_bar, "Delete Surface", DeleteSurface, NULL, " group='I/O Controls' ");


   /* CURVE MODELING CONTROLS */
   TwEnumVal curveEV[2] = {
      {cagd::curves::BEZIER_SAMPLED, "Bezier (Sampled)"},
      {cagd::curves::UNIFORM_3D_BSPLINE, "B-Spline (Uniform Knots)"},
   };
   TwType curveEVType = TwDefineEnum("curveEVType", curveEV, 2);
   TwAddVarRW(main_bar, "Curve Type", curveEVType, &curve_type,
              " help='Select curve type' group='Curve Modeling' ");

   TwAddVarRW(main_bar, "Curve Samples", TW_TYPE_INT32, &surface_params.samples,
              " min=1 max=100 help='Number of samples for c(t)' group='Curve Modeling' ");
   TwAddVarRW(main_bar, "B-Spline order", TW_TYPE_INT32, &surface_params.order,
              " min=1 max=10 help='Order for B-Spline' group='Curve Modeling' ");
   TwAddButton(main_bar, "Generate Curve", GenerateCurve, NULL, " group='Curve Modeling' ");
   TwAddButton(main_bar, "Clear Curve", ClearCurve, NULL, " group='Curve Modeling' ");

   /* SURFACE MODELING CONTROLS */
   TwEnumVal surfTypeEV[1] = {
      {cagd::surf::DOOSABIN, "Doo-Sabin"},
   };
   TwType surfEVType = TwDefineEnum("surfEVType", surfTypeEV, 1);
   TwAddVarRW(main_bar, "Subdivision Type", surfEVType, &surface_params.surf_type,
              " help='Select subdivision type' group='Surface Modeling' ");
   TwAddVarRW(main_bar, "Subdivisons", TW_TYPE_INT32, &surface_params.subdivs,
              " min=0 max=10 help='Number of subdivisions' group='Surface Modeling' ");

   horn_bar = TwNewBar("Horn Controls");
   TwDefine(" 'Horn Controls' size='300 140' valueswidth=fit contained=true ");

   TwEnumVal extruTypeEV[2] = {
      {cagd::surf::EXTRUDE_ALONG_NORMAL, "Extrude along normal"},
      {cagd::surf::EXTRUDE_ALONG_PATH, "Extrude along path"},
   };
   TwType extruEVType = TwDefineEnum("extruEVType", extruTypeEV, 2);
   TwAddVarRW(horn_bar, "Extrusion Type", extruEVType, &surface_params.extru_type,
              " help='Select extrusion type' group='Horn Modeling' ");

   TwAddVarRW(horn_bar, "Steps", TW_TYPE_INT32, &surface_params.horn_params.steps,
              " min=0 help='Number of steps to extrude' group='Horn Modeling' ");
   TwAddVarRW(horn_bar, "Step Size", TW_TYPE_DOUBLE, &surface_params.horn_params.hn,
              " help='Amount of extrusion per step' group='Horn Modeling' ");
   TwAddVarRW(horn_bar, "Step Scale", TW_TYPE_DOUBLE, &surface_params.horn_params.sn,
              " help='Amount of scaling per step' group='Horn Modeling' ");
   TwAddVarRW(horn_bar, "Step Rotation", TW_TYPE_DOUBLE, &surface_params.horn_params.theta,
              " help='Amount of rotation per step' group='Horn Modeling' ");


   cam = new OGL_Camera(ogl::Vector3d(0, 0, 1), ogl::Vector3d(0, 0, 0), ogl::Vector3d(0, 1, 0));

   double glyph_radius = 0.1;
   size_t glyph_segs = 16;

   double step = (2.0 * M_PI) / static_cast<double>(glyph_segs);
   point_glyph.resize(glyph_segs);
   for(size_t i = 0; i < glyph_segs; i++)
   {
      double angle = (step * i);
      point_glyph[i] = lux::Vector(glyph_radius * std::cos(angle), glyph_radius * std::sin(angle), 0.0);
   }

   glutMainLoop();

   return 0;
}

void initGLUT(int argc, char **argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
   glutInitWindowSize(windowWidth, windowHeight);
   glutInitWindowPosition(windowPositionX, windowPositionY);
   windowHandle = glutCreateWindow("Surface Modeling");

   glutDisplayFunc(draw);
   glutReshapeFunc(reshapeCB);
   glutKeyboardFunc(keyboardCB);
   glutMouseFunc(mouseButtonCB);
   glutMotionFunc(mouseMotionCB);
   glutIdleFunc(IdleCB);
   glutTimerFunc(display_update, timerCB, display_update);

   glutSpecialFunc((GLUTspecialfun)TwEventSpecialGLUT);
   TwGLUTModifiersFunc(glutGetModifiers);
}

void initLights()
{
   glLightModelf(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);

   GLfloat light0_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light0_Kd[] = {0.7f, 0.0f, 0.0f, 1.0f};  // diffuse light
   GLfloat light0_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT0, GL_AMBIENT, light0_Ka);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_Kd);
   glLightfv(GL_LIGHT0, GL_SPECULAR, light0_Ks);
   float light0_Pos[4] = {15, 0, 0, 1}; // positional light
   glLightfv(GL_LIGHT0, GL_POSITION, light0_Pos);

   GLfloat light1_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light1_Kd[] = {0.0f, 0.7f, 0.0f, 1.0f};  // diffuse light
   GLfloat light1_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT1, GL_AMBIENT, light1_Ka);
   glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_Kd);
   glLightfv(GL_LIGHT1, GL_SPECULAR, light1_Ks);
   float light1_Pos[4] = {-15, 0, 0, 1}; // positional light
   glLightfv(GL_LIGHT1, GL_POSITION, light1_Pos);


   GLfloat light2_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light2_Kd[] = {0.0f, 1.0f, 1.0f, 1.0f};  // diffuse light
   GLfloat light2_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT2, GL_AMBIENT, light2_Ka);
   glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_Kd);
   glLightfv(GL_LIGHT2, GL_SPECULAR, light2_Ks);
   float light2_Pos[4] = {0, 15, 0, 1}; // light position
   float light2_Dir[4] = {0,  1, 0, 1}; // light direction
   glLightfv(GL_LIGHT2, GL_POSITION, light2_Pos);
   glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light2_Dir);

   GLfloat light3_Ka[] = {0.0f, 0.0f, 0.0f, 1.0f};  // ambient light
   GLfloat light3_Kd[] = {1.0f, 1.0f, 1.0f, 1.0f};  // diffuse light
   GLfloat light3_Ks[] = {1, 1, 1, 1};           // specular light
   glLightfv(GL_LIGHT3, GL_AMBIENT, light3_Ka);
   glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_Kd);
   glLightfv(GL_LIGHT3, GL_SPECULAR, light3_Ks);
   float light3_Pos[4] = {0, 15, 0, 1}; // light position
   float light3_Dir[4] = {0, -1, 0, 1}; // light direction
   glLightfv(GL_LIGHT3, GL_POSITION, light3_Pos);
   glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, light3_Dir);

   //glEnable(GL_LIGHT0);
   //glEnable(GL_LIGHT1);
   //glEnable(GL_LIGHT2);
   glEnable(GL_LIGHT3);
   //glEnable(GL_LIGHTING);
}

void initGL()
{
   glShadeModel(GL_SMOOTH);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
   glEnable(GL_DEPTH_TEST);

   glClearColor(0, 0, 0, 0);
   glClearDepth(1.0f);
   glDepthFunc(GL_LEQUAL);

   glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
   glEnable(GL_COLOR_MATERIAL);

   initLights();
}

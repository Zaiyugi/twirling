#ifndef BOX3D
#define BOX3D
#ifndef SPHERE
#define SPHERE
#ifndef COLLPAIR
#define COLLPAIR
#include "global_prototypes.h"
using namespace std;

extern double ***sphereFrame;

class Material
{
   public:
      int ID;
      int illum;

      float *Kd;
      float *Ka;
      float *Ks;
      float *Ke;
      float *Ns;

      Material()
      {
         illum = -1;
         Kd = NULL;
         Ka = NULL;
         Ks = NULL;
         Ke = NULL;
         Ns = NULL;
      }

      void applyMaterial()
      {
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
         glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
         glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, Ke);
         glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, Ns);
      }

      void loadMaterialfromFile(FILE *mat_fd)
      {
         char *junk = new char[50];

         // Diffuse
         Kd = new float[4];
         fscanf(mat_fd, "%s %f %f %f %f", junk, Kd, (Kd + 1), (Kd + 2), (Kd + 3));

         // Ambient
         Ka = new float[4];
         fscanf(mat_fd, "%s %f %f %f %f", junk, Ka, (Ka + 1), (Ka + 2), (Ka + 3));

         // Specular
         Ks = new float[4];
         fscanf(mat_fd, "%s %f %f %f %f", junk, Ks, (Ks + 1), (Ks + 2), (Ks + 3));

         // Emissive
         Ke = new float[4];
         fscanf(mat_fd, "%s %f %f %f %f", junk, Ke, (Ke + 1), (Ke + 2), (Ke + 3));

         // Shininess
         Ns = new float;
         fscanf(mat_fd, "%s %f", junk, Ns);

         delete [] junk;
      }
};

// Wall Class is currently unused.
class Wall
{
   private:
      Vec3<double> v0, v1, v2, v3;
      Vec3<double> normal;

   public:
      double length, width;

      Wall()
      {
         length = width = 1;
         v0.set(0, 0, 0);
         v1.set(v0.x + length, v0.y, v0.z);
         v2.set(v0.x + length, v0.y, v0.z + width);
         v3.set(v0.x, v0.y, v0.z + width);
         normal.set(0, 1, 0);
      }

      Wall(Vec3<double> v0, Vec3<double> v1, 
           Vec3<double> v2, Vec3<double> v3)
      {
         this->v0 = v0;			// V0 is the smallest vertex
         this->v1 = v1;
         this->v2 = v2;			// V2 is the biggest vertex
         this->v3 = v3;

         Vec3<double> temp1 = v0 - v1;
         Vec3<double> temp2 = v2 - v1;
         normal = temp1.crossP(temp2);
      }

      void draw()
      {
         glBegin(GL_TRIANGLES);

            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v0.x, v0.y, v0.z);
            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v3.x, v3.y, v3.z);
            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v2.x, v2.y, v2.z);

            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v2.x, v2.y, v2.z);
            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v1.x, v1.y, v1.z);
            glNormal3d(normal.x, normal.y, normal.z);
            glVertex3d(v0.x, v0.y, v0.z);

         glEnd();
      }

      Vec3<double> getPt(int index)
      {
         switch(index)
         {
            case 0:
               return v0;

            case 1:
               return v1;

            case 2:
               return v2;

            case 3:
               return v3;

         }
      }

};

class CollPairLocal
{
   public:
      int indexA;
      int indexB;
      char typeA[10];
      char typeB[10];

      CollPairLocal()
      {
         indexA = -1;
         indexB = -1;
         strcpy(typeA, "Sphere");
         strcpy(typeB, "Sphere");
      }

      CollPairLocal(int a, int b, char aT[10], char bT[10])
      {
         indexA = a;
         indexB = b;
         strcpy(typeA, aT);
         strcpy(typeB, bT);
      }

};

class Sphere
{
   private:
      double radius;
      Vec3<double> up;
      bool resting;

   public:
      int mat_index;
      double mass;
      Vec3<double> velocity;
      Vec3<double> center;
      Vec3<double> new_velocity;
      Vec3<double> new_center;

      Sphere()
      {
         resting = false;
         mat_index = -1;
         center = Vec3<double>(0, 0, 0);
         up = Vec3<double>(0, 1, 0);
         radius = 1;
      }

      Sphere(double r, Vec3<double> c, Vec3<double> up, int mi)
      {
         resting = false;
         mat_index = mi;
         radius = r;
         center = c;
         this->up = up;
      }

      void draw()
      {
         glMatrixMode(GL_MODELVIEW);
         glPushMatrix();

         glTranslatef(center.x, center.y, center.z);

         glBegin(GL_POINTS);
               glNormal3d(up.x, up.y, up.z);
               glVertex3d(center.x, center.y, center.z);
         glEnd();

         glTranslatef(-center.x, -center.y, -center.z);

         glPopMatrix();
      }

      Vec3<double> getUpVector()
      {
         return up;
      }

      bool isResting()
      {
         return resting;
      }

      void setResting(bool r)
      {
         resting = r;
      }

};

class Box3D
{
   private:
      Vec3<double> pts[8];
      Vec3<double> up;

   public:
      double length, width, height;

      Box3D()
      {
         length = 1;
         height = 1;
         width = 1;
         up = Vec3<double>(0, 1, 0);
      }

      Box3D(double l, double w, double h, Vec3<double> v0, Vec3<double> up)
      {
         length = l;
         width = w;
         height = h;
         this->up = up;
         pts[0] = v0;
         pts[7] = v0 - Vec3<double>(length, height, width);
      }

      void fillInVertices()
      {
         pts[1] = Vec3<double>(pts[0].x - length, pts[0].y,          pts[0].z);
         pts[2] = Vec3<double>(pts[0].x - length, pts[0].y - height, pts[0].z);
         pts[3] = Vec3<double>(pts[0].x,          pts[0].y - height, pts[0].z);

         pts[4] = Vec3<double>(pts[7].x + length, pts[7].y,          pts[7].z);
         pts[5] = Vec3<double>(pts[7].x + length, pts[7].y + height, pts[7].z);
         pts[6] = Vec3<double>(pts[7].x,          pts[7].y + height, pts[7].z);
         
      }

      void draw()
      {
         glBegin(GL_TRIANGLES);
            // Bottom face
            glNormal3d(up.x, up.y, up.z);
            glVertex3f(pts[7].x, pts[7].y, pts[7].z);
            glVertex3f(pts[4].x, pts[4].y, pts[4].z);
            glVertex3f(pts[3].x, pts[3].y, pts[3].z);

            glVertex3f(pts[3].x, pts[3].y, pts[3].z);
            glVertex3f(pts[2].x, pts[2].y, pts[2].z);
            glVertex3f(pts[7].x, pts[7].y, pts[7].z);
            // End bottom

            // Left face
            glNormal3d(1, 0, 0);
            glVertex3f(pts[1].x, pts[1].y, pts[1].z);
            glVertex3f(pts[6].x, pts[6].y, pts[6].z);
            glVertex3f(pts[7].x, pts[7].y, pts[7].z);

            glVertex3f(pts[7].x, pts[7].y, pts[7].z);
            glVertex3f(pts[2].x, pts[2].y, pts[2].z);
            glVertex3f(pts[1].x, pts[1].y, pts[1].z);
            // End Left face

            // Back face
            glNormal3d(0, 0, 1);
            glVertex3f(pts[4].x, pts[4].y, pts[4].z);
            glVertex3f(pts[7].x, pts[7].y, pts[7].z);
            glVertex3f(pts[6].x, pts[6].y, pts[6].z);

            glVertex3f(pts[6].x, pts[6].y, pts[6].z);
            glVertex3f(pts[5].x, pts[5].y, pts[5].z);
            glVertex3f(pts[4].x, pts[4].y, pts[4].z);
            // End back face

            // Top face
            glNormal3d(0, -1, 0);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            glVertex3f(pts[5].x, pts[5].y, pts[5].z);
            glVertex3f(pts[6].x, pts[6].y, pts[6].z);

            glVertex3f(pts[6].x, pts[6].y, pts[6].z);
            glVertex3f(pts[1].x, pts[1].y, pts[1].z);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            // End top face

            // Right face
            glNormal3d(-1, 0, 0);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            glVertex3f(pts[3].x, pts[3].y, pts[3].z);
            glVertex3f(pts[4].x, pts[4].y, pts[4].z);

            glVertex3f(pts[4].x, pts[4].y, pts[4].z);
            glVertex3f(pts[5].x, pts[5].y, pts[5].z);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            // End right face

            // Front face
            glNormal3d(0, 0, -1);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            glVertex3f(pts[1].x, pts[1].y, pts[1].z);
            glVertex3f(pts[2].x, pts[2].y, pts[2].z);

            glVertex3f(pts[2].x, pts[2].y, pts[2].z);
            glVertex3f(pts[3].x, pts[3].y, pts[3].z);
            glVertex3f(pts[0].x, pts[0].y, pts[0].z);
            // End front face

         glEnd();
      }

      Vec3<double> getPt(int index)
      {
         return pts[index];
      }

      bool testResting(Sphere* ball, Vec3<double> accel)
      {
         double d, epsilon;
         Vec3<double> X = ball->new_center;
         Vec3<double> V = ball->new_velocity;

         d = (Vec3<double>(X.x, X.y, X.z) - pts[7]).dotP(up);
         epsilon = 0.0025;
         if(fabs(d) < epsilon)
         {
            if(fabs(V.dotP(up)) < epsilon)
            {
               if(accel.dotP(up) < epsilon)
               {
                  ball->setResting(true);

                  return true;
               }
            }
         }

         return false;
      }

      bool ballBoxCollision(int index, Sphere* ball, double dt, double COR, double COF, Vec3<double> accel)
      {
         /* Indexing scheme for walls:
          * 0 = Floor
          * 1 = Ceiling
          * 2 = Left Wall
          * 3 = Right Wall
          * 4 = Back Wall
          * 5 = Front Wall
          */

         Vec3<double> V_new = ball->new_velocity;
         Vec3<double> X_new = ball->new_center;
         Vec3<double> V_c, X_c, V_p, V_t, V_cp;

         switch(index)
         {
            case 0:
               if(X_new.y < getPt(7).y)
               {
                  //printf("Collide with bottom!\n");
                  double f = (ball->center.y - getPt(7).y);
                  f = f / (ball->center.y - X_new.y);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(0, 1, 0) * V_c.dotP(Vec3<double>(0, 1, 0));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;

            case 1:
               if(X_new.y > getPt(0).y)
               {
                  //printf("Collide with top!\n");
                  double f = (getPt(0).y - ball->center.y);
                  f = f / (X_new.y - ball->center.y);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(0, -1, 0) * V_c.dotP(Vec3<double>(0, -1, 0));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;

            case 2:
               if(X_new.x < getPt(7).x)
               {
                  //printf("Collide with left wall!\n");
                  double f = (ball->center.x - getPt(7).x);
                  f = f / (ball->center.x - X_new.x);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(1, 0, 0) * V_c.dotP(Vec3<double>(1, 0, 0));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;

            case 3:
               if(X_new.x > getPt(0).x)
               {
                  //printf("Collide with right wall!\n");
                  double f = (getPt(0).x - ball->center.x);
                  f = f / (X_new.x - ball->center.x);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(-1, 0, 0) * V_c.dotP(Vec3<double>(-1, 0, 0));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;

            case 4:
               if(X_new.z < getPt(7).z)
               {
                  //printf("Collide with back wall!\n");
                  double f = (ball->center.z - getPt(7).z);
                  f = f / (ball->center.z - X_new.z);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(0, 0, 1) * V_c.dotP(Vec3<double>(0, 0, 1));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;

            case 5:
               if(X_new.z > getPt(0).z)
               {
                  //printf("Collide with front wall!\n");
                  double f = (getPt(0).z - ball->center.z);
                  f = f / (X_new.z - ball->center.z);
                  V_c = ball->velocity + accel * f * dt;
                  X_c = ball->center + ball->velocity * f * dt;

                  // Frictional Reflect
                  V_p = Vec3<double>(0, 0, -1) * V_c.dotP(Vec3<double>(0, 0, -1));
                  V_t = V_c - V_p;
                  V_cp = V_p * -COR + V_t * (1 - COF);
                  V_c = V_cp;

                  V_new = V_c + accel * (1 - f) * dt;
                  X_new = X_c + V_c * (1 - f) * dt;
               }
               break;
         }

         ball->new_velocity = V_new;
         ball->new_center = X_new;
      }

      Vec3<double> getCenter()
      {
         return (pts[0] + pts[7]) / 2;
      }

      Vec3<double> getSmallestPt()
      {
         return pts[7];
      }

      Vec3<double> getBiggestPt()
      {
         return pts[0];
      }

      Vec3<double> getUpVector()
      {
         return up;
      }
};

#endif
#endif
#endif

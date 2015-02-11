#ifndef BLOCK
#define BLOCK
#ifndef VERTICES
#define VERTICES
#ifndef QUATERNION
#define QUATERNION
#ifndef VEC3
#define VEC3

#include "global_prototypes.h"
#define PIOVER180 0.0174532925
using namespace std;

template <class T>
class Vec3
{
   public:
      T x;
      T y;
      T z;

   Vec3() {
      x = y = z = NULL;
   }

   Vec3(T x, T y, T z) {
      this->x = x;
      this->y = y;
      this->z = z;
   }

   void set(T x, T y, T z)
   {
      this->x = x;
      this->y = y;
      this->z = z;
   }

   T mag()
   {
      T temp; 
      temp = x * x + y * y + z * z;

   return sqrt(temp);
   }

   void normalize()
   {
      T mag;
      mag = sqrt(x*x + y*y + z*z);

      x = x / mag;
      y = y / mag;
      z = z / mag;
   }

   T dotP(Vec3<T> v2)
   {
      T ans;
      ans = x * v2.x + y * v2.y + z * v2.z;
      return ans;
   }

   Vec3<T> crossP(Vec3<T> v2)
   {
      Vec3<T> ans;
      ans.x = y * v2.z - z * v2.y;
      ans.y = z * v2.x - x * v2.z;
      ans.z = x * v2.y - y * v2.x;

      return ans;
   }

   Vec3<T> operator*(const T &scalar) const
   {
      Vec3<T> v3;
      v3.set(scalar * this->x, scalar * this->y, scalar * this->z);

      return v3;
   }

   Vec3<T> operator/(const T &scalar) const
   {
      Vec3<T> v3;
      v3.set(this->x / scalar, this->y / scalar, this->z / scalar);

      return v3;
   }

   Vec3<T> operator+(const Vec3<T> &v2) const
   {
      Vec3<T> v3;
      v3.set(this->x + v2.x, this->y + v2.y, this->z + v2.z);

   return v3;
   }

   Vec3<T> operator-(const Vec3<T> &v2) const
   {
      Vec3<T> v3;
      v3.set(this->x - v2.x, this->y - v2.y, this->z - v2.z);

   return v3;
   }

   bool operator>=(const Vec3<T> &v2) const
   {
      if(this->x >= v2.x)
         if(this->y >= v2.y)
            if(this->z >= v2.z)
               return true;

   return false;
   }

   bool operator<=(const Vec3<T> &v2) const
   {
      if(this->x <= v2.x)
         if(this->y <= v2.y)
            if(this->z <= v2.z)
               return true;

   return false;
   }

   bool operator==(const Vec3<T> &v2) const
   {
      if(this->x == v2.x)
         if(this->y == v2.y)
            if(this->z == v2.z)
               return true;

   return false;
   }

   bool operator!=(const Vec3<T> &v2) const
   {
      if(this->x != v2.x)
         return true;
      if(this->y != v2.y)
         return true;
      if(this->z != v2.z)
         return true;

   return false;
   }

};

class Quaternion
{
   public:
      float w;
      float x;
      float y;
      float z;

      // Default Quaternion
      // Used in Block Definition
      // Defines Quaternion with 0 rotation
      Quaternion()
      {
         w = 1; x = 0; y = 0; z = 0;
      }

      // Build quaternion from axis and angle
      // q = (s,v) 
      Quaternion(float angle, Vec3<float> axis)
      {
         float rad = PIOVER180 * angle;
         rad = rad / 2;

         this->w = cos(rad);
         this->x = axis.x * sin(rad);
         this->y = axis.y * sin(rad);
         this->z = axis.z * sin(rad);
      }

      // Set previous declared quaternion 
      // q = (s,v) 
      void set(float angle, Vec3<float> axis)
      {
         float rad = PIOVER180 * angle;
         rad = rad / 2;

         this->w = cos(rad);
         this->x = axis.x * sin(rad);
         this->y = axis.y * sin(rad);
         this->z = axis.z * sin(rad);
      }

      Quaternion conjugate()
      {
         Quaternion q;
         q.w = w;
         q.x = -x;
         q.y = -y;
         q.z = -z;

         return q;
      }

      // Check if quaternion is within tolerance
      // If not, normalize the quaternion
      void normalize()
      {
         float TOL = 0.00005;
         float mag = sqrt(w*w + x*x + y*y + z*z);

         if(fabs(mag - 1) > TOL)
         {
            w = w / mag;
            x = x / mag;
            y = y / mag;
            z = z / mag;
         }
      }

      // Concatenate two quaternions by multiplication
      // q2 * q1 = (s2*s1 - dotp(v2,v1), s2*v1 + s1*v2 + crossp(v2,v1))
      Quaternion operator*(Quaternion &q2)
      {
         Quaternion q3;

         q3.w = (this->w * q2.w) - (this->x * q2.x) - (this->y * q2.y) - (this->z * q2.z);
         q3.x = (this->w * q2.x) + (this->x * q2.w) + (this->y * q2.z) - (this->z * q2.y);
         q3.y = (this->w * q2.y) - (this->x * q2.z) + (this->y * q2.w) + (this->z * q2.x);
         q3.z = (this->w * q2.z) + (this->x * q2.y) - (this->y * q2.x) + (this->z * q2.w);

      return q3;
      }
};

class Normals
{
   private:
      Vec3<float> *normals;

   public:

      Normals()
      {
         normals = new Vec3<float>[6];
         normals[0].set(0, 0, 1);
         normals[1].set(1, 0, 0);
         normals[2].set(0, 1, 0);
         normals[0].set(-1, 0, 0);
         normals[0].set(0, -1, 0);
         normals[0].set(0, 0, -1);
      }

      Vec3<float> getNormal(char axis[3])
      {
         if(axis[0] == '+')
         {
            switch(axis[1])
            {
               case 'x':
                  return normals[1];
                  break;

               case 'y':
                  return normals[2];
                  break;

               case 'z':
                  return normals[0];
                  break;
            }
         } else if(axis[0] == '-') {
            switch(axis[1])
            {
               case 'x':
                  return normals[3];
                  break;

               case 'y':
                  return normals[4];
                  break;

               case 'z':
                  return normals[5];
                  break;
            }
         }
      }

};

class Camera_Mine
{
   public:
      float distance;
      float pan;
      float elev;

      Vec3<float> Center;
      Vec3<float> Eye;
      Vec3<float> Up;
      Vec3<float> Right;

      Vec3<float> localXAxis;
      Vec3<float> localZAxis;

      float prevy;
      float prevx;

      Quaternion camRotX;
      Quaternion camRotY;

      Camera_Mine()
      {
         Eye.set(0, 0, 1);
         Center.set(0, 0, 0);
         Up.set(0, 1, 0);
         Right.set(1, 0, 0);

         distance = pan = elev = 0;
         prevy = prevx = 0.0;

      }

      void moveRight(float amt)
      {
         pan += amt;
         //localXAxis = localXAxis + localXAxis * amt;
      }

      void moveLeft(float amt)
      {
         pan -= amt;
         //localXAxis = localXAxis + localXAxis * amt;
      }

      void Zoom(float amt)
      {
         distance += amt;
         //localZAxis = localZAxis + localZAxis * amt;
      }

      void snapTo(float x, float y)
      {
         prevy = y;
         prevx = x;
      }

      void updateAngleAboutX(float deltaY, float scale, float y)
      {
         Quaternion local;

         // Set Rotation around Camera's Right axis (WORLD)
         local.set(deltaY * scale * 1.0, Up);

         local.normalize();
         camRotX = local * camRotX;

         prevy = y;
      }

      void updateAngleAboutY(float deltaX, float scale, float x)
      {
         Quaternion local;

         // Set Rotation around Camera's Up axis (WORLD)
         local.set(deltaX * scale * 1.0, Right);

         local.normalize();
         camRotY = local * camRotY;

         prevx = x;
      }

      void setModelView()
      {
         double M[16];

         construct(M, camRotY);
         glMultMatrixd(M);
         construct(M, camRotX);
         glMultMatrixd(M);

         glTranslatef(pan, elev, distance);
      }

      void invertModelView()
      {
         double M[16];
         Quaternion inverse_x;
         Quaternion inverse_y;

         inverse_x = camRotX.conjugate();
         inverse_y = camRotY.conjugate();

         glTranslatef(-pan, -elev, -distance);

         construct(M, inverse_x);
         glMultMatrixd(M);
         construct(M, inverse_y);
         glMultMatrixd(M);
      }

      void construct(double mat[16], Quaternion Q1)
      {
         double w,x,y,z;
         w = Q1.w; x = Q1.x; y = Q1.y; z = Q1.z;

         mat[0] = (w*w)+(x*x)-(y*y)-(z*z); mat[4] = (2*x*y)+(2*w*z);         mat[8] = (2*x*z)-(2*w*y);          mat[12] = 0;
         mat[1] = (2*x*y)-(2*w*z); 	   mat[5] = (w*w)-(x*x)+(y*y)-(z*z); mat[9] = (2*y*z)+(2*w*x);          mat[13] = 0;
         mat[2] = (2*x*z)+(2*w*y);         mat[6] = (2*y*z)-(2*w*x);         mat[10] = (w*w)-(x*x)-(y*y)+(z*z); mat[14] = 0;
         mat[3] = 0;                       mat[7] = 0;                       mat[11] = 0;                       mat[15] = 1;
      }

};

class Vertices
{
   private:
      float x;
      float y;
      float z;

   public:
      float v0[3];
      float v1[3];
      float v2[3];
      float v3[3];
      float v4[3];
      float v5[3];
      float v6[3];
      float v7[3];

      Vertices()
      {
         x = y = z = 0;
      }

      // Set vertices to be 0.5 units from center of mass
      // Creates a cube of size 1x1x1 units
      void vertUpdate(float x, float y, float z)
      {
         this->x = x; this->y = y; this->z = z;

         v0[0] = x + 0.5; v0[1] = y + 0.5; v0[2] = z + 0.5;
         v1[0] = x - 0.5; v1[1] = y + 0.5; v1[2] = z + 0.5;
         v2[0] = x - 0.5; v2[1] = y - 0.5; v2[2] = z + 0.5;
         v3[0] = x + 0.5; v3[1] = y - 0.5; v3[2] = z + 0.5;
         v4[0] = x + 0.5; v4[1] = y - 0.5; v4[2] = z - 0.5;
         v5[0] = x + 0.5; v5[1] = y + 0.5; v5[2] = z - 0.5;
         v6[0] = x - 0.5; v6[1] = y + 0.5; v6[2] = z - 0.5;
         v7[0] = x - 0.5; v7[1] = y - 0.5; v7[2] = z - 0.5;
      }

      Vec3<float> getVertex(int n)
      {
         Vec3<float> temp;
         switch(n)
         {
            case 0:
               temp.set(v0[0], v0[1], v0[2]);
               break;

            case 7:
               temp.set(v7[0], v7[1], v7[2]);
         }

         return temp;
      }

};

class Block
{
   private:
      Vertices v;		// Vertices of cube
      Quaternion quat;		// Quaternion for rotation

      // Color: RGB 
      float r;
      float g;
      float b;

   public:
      Vec3<float> vel;		// Velocity
      Vec3<float> c;		// Center of Mass
      bool quatFlag;
      bool picked;
      bool active;
      short int id;
      short int stackLvl;

      Block()
      {
         id = stackLvl = -1;
         c.set(0,0,0);
         v.vertUpdate(0,0,0);
         vel.set(0,0,0);
         quatFlag = picked = active = false;
      }

      Block(int id, float x, float y, float z, float r, float g, float b)
      {
         this->id = id;
         stackLvl = -1;
         c.set(x,y,z);
         v.vertUpdate(c.x, c.y, c.z);
         vel.set(0,0,0);
         quatFlag = picked = active = false;
         this->r = r; this->g = g; this->b = b;
      }

      // Update vertices to reflect new center of mass
      void update()
      {
         v.vertUpdate(c.x, c.y, c.z);
      }

      Vec3<float> getVertex(int n)
      {
         return v.getVertex(n);
      }

      // Render a 1x1x1 cube
      // Note: Should be updated to better method
      // Located in 'render.cpp'
      void renderBlock();

      // Check if ray point is inside of block
      // Used by the 'select(x, y)'
      bool checkVertices(Vec3<float> ray)
      {
//         cout << "checking vertices inside block func...\n";
//         printf("IS PT(%f %f %f) ", ray.x, ray.y, ray.z);
//         printf(" INSIDE PT(%f %f %f) ", v.v0[0], v.v0[1], v.v0[2]);
//         printf(" AND PT(%f %f %f)\n", v.v7[0], v.v7[1], v.v7[2]);
         if(ray.x <= v.v0[0] && ray.x >= v.v7[0]) {
            if(ray.y <= v.v0[1] && ray.y >= v.v7[1])
               if(ray.z <= v.v0[2] && ray.z >= v.v7[2])
                  return true;
         } else
            return false;
      }

      // Update current quaternion to reflect new rotation
      void updateQuat(float angle, Vec3<float> axis)
      {
         // Create quaternion from axis-angle representation
         Quaternion local(angle, axis);
         local.normalize();

         // Multiply local rotation and current rotation
         // Uses operator overload
         quat = local * quat;

         // Set flag to reflect rotated status
         quatFlag = true;
      }

      void rotate()
      {
         float quatMat[16];

         // Construct rotation matrix from quaternion
         // Matrix is in Column-Major form for glMultMatrix
         construct(quatMat);

         glMatrixMode(GL_MODELVIEW);
         glPushMatrix();

         glTranslatef(c.x, c.y, c.z);
         glMultMatrixf(quatMat);
         glTranslatef(-c.x, -c.y, -c.z);

         renderBlock();

         glPopMatrix();
      }

      void construct(float mat[16])
      {
         float w,x,y,z;
         w = quat.w; x = quat.x; y = quat.y; z = quat.z;

         mat[0] = (w*w)+(x*x)-(y*y)-(z*z); mat[4] = (2*x*y)+(2*w*z);         mat[8] = (2*x*z)-(2*w*y);          mat[12] = 0;
         mat[1] = (2*x*y)-(2*w*z); 	   mat[5] = (w*w)-(x*x)+(y*y)-(z*z); mat[9] = (2*y*z)+(2*w*x);          mat[13] = 0;
         mat[2] = (2*x*z)+(2*w*y);         mat[6] = (2*y*z)-(2*w*x);         mat[10] = (w*w)-(x*x)-(y*y)+(z*z); mat[14] = 0;
         mat[3] = 0;                       mat[7] = 0;                       mat[11] = 0;                       mat[15] = 1;
      }
};

class CollPair
{
   public:
      Block *a;
      Block *b;
      double distance;

   CollPair()
   {
      a = new Block();
      b = new Block();
      distance = 0;
   }

   CollPair(Block *b1, Block *b2)
   {
      a = b1;
      b = b2;
      distance = 0;
   }

   ~CollPair()
   {
      delete a;
      delete b;
   }

};

#endif
#endif
#endif
#endif


/* Zachary Shore
 * CPSC-881-001
 * 2014-02-21
 * Geometric Modeling
 * Control Polytope class
 */
#ifndef __CONTROL_POLYTOPE_H__
#define __CONTROL_POLYTOPE_H__
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cmath>

#include <vector>
#include <algorithm>

#include "opengl/ogl_include.h"
#include "Vector.h"

namespace cagd
{

class ControlPolytope
{
   public:

      ControlPolytope()
      {
         glyph_radius = 0.025;
         glyph_segs = 16;

         double step = (2.0 * M_PI) / static_cast<double>(glyph_segs);
         glyph.resize(glyph_segs);
         for(int i = 0; i < glyph_segs; i++)
         {
            double angle = (step * i);
            glyph[i] = lux::Vector(glyph_radius * std::cos(angle), glyph_radius * std::sin(angle), 0.0);
         }

         order = 1;

         closed = false;

         points.resize(1);
      };

      ~ControlPolytope() {};

      void setFirstSlice(std::vector<lux::Vector>& slice)
      { points[0] = slice; }

      void setSlice(int i, std::vector<lux::Vector>& slice)
      { points[i] = slice; }

      void addSlice(std::vector<lux::Vector>& slice)
      { points.push_back(slice); }

      void draw(lux::Vector right, lux::Vector up, lux::Vector color)
      {
         int pointsPerSlice = points[0].size();
         int slices = points.size();
         for(int i = 0; i < slices; i++)
         {
            for(int j = 0; j < pointsPerSlice; j++)
            {
               lux::Vector pt = points[i][j];

               if(i == selected_slice && j == selected_ndx)
                  glColor4f(1.0 - color[0], 1.0 - color[1], 1.0 - color[2], 1.0);
               else
                  glColor4f(color[0], color[1], color[2], 1.0);

               glBegin(GL_POINTS);
                  glNormal3f(0, 1, 0);
                  glVertex3f(pt[0], pt[1], pt[2]);
               glEnd();

               glBegin(GL_LINE_LOOP);
                  for(size_t k = 0; k < glyph.size(); k++)
                  {
                     lux::Vector X = pt
                                   + right * glyph[k][0]
                                   + up * glyph[k][1];

                     glNormal3f(0, 1, 0);
                     glVertex3f(X[0], X[1], X[2]);
                  }
               glEnd();

            }
         }

         glColor4f(0.35, 0.35, 0.35, 1.0);
         for(int i = 0; i < slices; i++)
         {
            glBegin(GL_LINE_STRIP);

            for(int j = 0; j < pointsPerSlice; j++)
            {
               lux::Vector pt = points[i][j];
               glNormal3f(0, 1, 0);
               glVertex3f(pt[0], pt[1], pt[2]);
            }

            glEnd();
         }

         for(int i = 0; i < pointsPerSlice; i++)
         {
            glBegin(GL_LINE_STRIP);

            for(int j = 0; j < slices; j++)
            {
               lux::Vector pt = points[j][i];
               glNormal3f(0, 1, 0);
               glVertex3f(pt[0], pt[1], pt[2]);
            }

            glEnd();
         }

      }

      bool testPoints(lux::Vector r_pos, lux::Vector r_dir, int& slice, int& ndx)
      {
         double radius = glyph_radius;
         double best_t = INT_MAX;
         double best_ndx = -1;
         double best_slice = -1;
         size_t start = (closed) ? order : 0;

         for(size_t i = 0; i < points.size(); i++)
         {
            for(size_t j = start; j < points[i].size(); j++)
            {

               lux::Vector pt = points[i][j];

               double A = r_dir * r_dir;
               double B = 2.0 * (r_pos - pt) * r_dir;
               double C = (r_pos - pt) * (r_pos - pt) - radius * radius;

               double disc = B * B - 4.0 * A * C;
               if(disc < 0.0)
                  continue;

               double discSqrt = std::sqrt(disc);
               double Q;
               if(B < 0)
                  Q = (-B - discSqrt) / 2.0;
               else
                  Q = (-B + discSqrt) / 2.0;

               double t0 = Q / A;
               double t1 = C / Q;

               if(t0 > t1)
               {
                  double temp = t0;
                  t0 = t1;
                  t1 = temp;
               }

               if(t1 < 0)
                  continue;

               double t = t0;
               if(t0 < 0)
                  t = t1;

               if(t == t0 || t == t1)
                  if(t < best_t)
                  {
                     best_t = t0;
                     best_slice = i;
                     best_ndx = j;
                  }
            }
         }

         ndx = best_ndx;
         slice = best_slice;
         if(slice > -1 && ndx > -1)
            return true;

         return false;
      }

      void clearSelected()
      { 
         selected_slice = -1;
         selected_ndx = -1;
      }

      void setSelected(int ss, int sn)
      {
         selected_slice = ss;
         selected_ndx = sn;
      }

      lux::Vector getSelected()
      {
         return points[selected_slice][selected_ndx];
      }

      int close()
      { 
         if(closed) return 0;

         for(size_t i = 0; i < points.size(); i++)
            for(int j = 0; j < order; j++)
               points[i].push_back(points[i][j]);

         closed = true;

         return 1;
      }

      int open()
      {
         if(!closed) return 0;

         for(size_t i = 0; i < points.size(); i++)
            for(int j = 0; j < order; j++)
               points.pop_back();

         closed = false;

         return 1;
      }

      void updateOrder(int o)
      {
         //int ret = close();
         //if(ret)
         //   open();
         order = o;
         //if(ret)
         //   close();
      };

      void setControlPoint(int i, int j, lux::Vector P)
      {
         points[i][j] = P;
         /*
         if(closed)
         {
            size_t N_cp = points.size();
            if(i < order)
               points[N_cp - (order - i)] = P;
            else if(i > (N_cp - order) - 1)
               points[i - (N_cp - order)] = P;
         }
         */
      }

      lux::Vector getControlPoint(int i, int j) { return points[i][j]; }

      int getSlices() { return points.size(); }
      int getPointsPerSlice() { return points[0].size(); }

      std::vector<lux::Vector>& operator[] ( size_t n ) { return points[n]; }

      lux::Vector getCentroid(int s)
      {
         lux::Vector cog(0, 0, 0);
         for(size_t i = 0; i < points[s].size(); i++)
            cog = cog + points[s][i];
         return cog * (1.0 / static_cast<double>(points[s].size()));
      }

      lux::Vector getNormal(int s)
      {
         lux::Vector normal(0, 0, 0);
         lux::Vector cog = getCentroid(s);

         int count = 0;
         for(size_t i = 0; i < points[s].size() - 1; i++)
         {
            lux::Vector v1 = (points[s][i] - cog);
            lux::Vector v2 = (points[s][i+1] - cog);
            normal = normal + (v1^v2);
            count++;
         }

         normal = normal * (1.0 / static_cast<double>(count));
         return normal.unitvector();
      }

      void clear() { points.clear(); }

      std::vector< std::vector<lux::Vector> > points;

   protected:

      int glyph_segs;
      double glyph_radius;
      std::vector<lux::Vector> glyph;

      int order;

      int selected_slice, selected_ndx;
      bool closed;
};

}

#endif

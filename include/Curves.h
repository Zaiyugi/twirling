/* Zachary Shore
 * CPSC-881-001
 * Created: 2014-01-30
 * Edited: 2014-04-15
 * Geometric Modeling
 * Curve classes
 */
#ifndef __CURVES_H__
#define __CURVES_H__
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cmath>

#include <vector>
#include <algorithm>

#include "opengl/ogl_include.h"
#include "LinearAlgebra.h"

namespace cagd
{

namespace curves
{

const unsigned int NULL_CURVE_TYPE = 0;
const unsigned int BEZIER_SAMPLED = 1;
const unsigned int BEZIER_SUBDIV = 2;
const unsigned int UNIFORM_3D_BSPLINE = 3;

class Frame
{
   public:
     lux::Vector r, s, t;

     Frame() {};
     Frame(lux::Vector _r, lux::Vector _s, lux::Vector _t) : r(_r), s(_s), t(_t) {};

     void transform(lux::Matrix mat)
     {
         r = r * mat;
         s = s * mat;
         t = t * mat;
     }
};

/* ~ Generate Rotation Minimizing Frames ~ */
void generateFrame(std::vector<lux::Vector>& X, Frame U_0, std::vector<Frame>& U)
{
   //U.resize(X.size());
   U.clear();
   U.push_back(U_0);

   Frame U_ip1;
   for(size_t i = 0; i < X.size()-1; ++i)
   {
      std::printf("Frame: %lu\n", i);

      lux::Vector v1 = X[i+1] - X[i];
      double c1 = v1 * v1;
      U_ip1.t = (X[i+1] - X[i]).unitvector();

      lux::Vector r_L = (U[i].r - (2.0/c1) * (v1 * U[i].r) * v1).unitvector();
      lux::Vector t_L = (U[i].t - (2.0/c1) * (v1 * U[i].t) * v1).unitvector();
      lux::Vector v2 = U_ip1.t - t_L;
      double c2 = v2 * v2;

      U_ip1.r = (r_L - (2.0/c2) * (v2 * r_L) * v2).unitvector();
      U_ip1.s = (U_ip1.t ^ U_ip1.r);

      U.push_back(U_ip1);
   }

}

/* ~ Control Polygon Intersection Test ~ */
bool testCurve(std::vector<lux::Vector> CP, bool closed, int order, lux::Vector r_pos, lux::Vector r_dir, int& x)
{
   double radius = 0.1;
   double best_t = INT_MAX;
   double best_ndx = -1;
   size_t start = (closed) ? order : 0;

   x = -1;
   for(size_t i = start; i < CP.size(); i++)
   {
      lux::Vector pt = CP[i];

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
            best_ndx = i;
         }
   }

   x = best_ndx;
   if(x > -1)
      return true;

   return false;
}

/* ~ Curve generation functions ~ */

/* Bezier - de Casteljau Sampling */
lux::Vector dC_Sample(std::vector<lux::Vector> P, double t)
{
   int n = P.size();
   while(n > 0)
   {
      n -= 1;
      for(int i = 0; i < n; i++)
         P[i] = (1.0 - t) * P[i] + t * P[i+1];
   }

   return P[0];
}

void deCasteljauSampling(std::vector<lux::Vector>& cp, std::vector<lux::Vector>& curve, int samples)
{
   curve.resize(samples);

   int i = 0;
   double s = static_cast<double>(samples-1);
   for(; i < samples-1; i++)
   {
      double t = static_cast<double>(i) / s;
      curve[i] = dC_Sample(cp, t);
   }
   curve[i] = dC_Sample(cp, 1.0);
}

/* Bezier - de Casteljau Subdivision */
std::vector<lux::Vector> dC_Subdivide(std::vector<lux::Vector> P, int n, double u)
{
   std::vector<lux::Vector> L;
   std::vector<lux::Vector> R;

   if(n == 1)
      return P;

   size_t N = P.size();
   while(N > 1)
   {
      L.push_back(P[0]); // L = L . P_0

      R.push_back(P[N-1]); // R = R . P_N

      N -= 1;
      for(size_t i = 0; i < N; i++)
         P[i] = (1.0 - u) * P[i] + u * P[i+1];
   }
   L.push_back(P[0]);
   R.push_back(P[0]);

   std::vector<lux::Vector> concat;

   L = dC_Subdivide(L, n-1, u);
   for_each(L.begin(), L.end(), [&] (lux::Vector X)
   {
      concat.push_back(X);
   } );

   concat.push_back(P[0]);

   R = dC_Subdivide(R, n-1, u);
   for_each(R.rbegin(), R.rend(), [&] (lux::Vector X)
   {
      concat.push_back(X);
   } );

   return concat;
}

void deCasteljauSubdivide(std::vector<lux::Vector>& cp, std::vector<lux::Vector>& curve, int subdiv)
{
   double u = 0.5;
   curve = dC_Subdivide(cp, subdiv, u);
}

/* B-Spline - Uniform knot vector w/ arbitrary order */
double knots(int i) { return static_cast<double>(i); }

double NormalizedBasis(int i, int d, double t)
{
   if(d == 1)
   {
      if(knots(i) <= t && t < knots(i+1))
         return 1;

      return 0;
   }

   double A = (t - knots(i))
            / (knots(i+d-1) - knots(i));

   double B = (knots(i+d) - t) 
            / (knots(i+d) - knots(i+1));

   double ret = A * NormalizedBasis(i, d-1, t) + B * NormalizedBasis(i+1, d-1, t);

   return ret;
}

void BSplineUniformKnotsArbitraryOrder(std::vector<lux::Vector>& cp, std::vector<lux::Vector>& curve, int order, int samples)
{
   size_t n = cp.size();
   curve.clear();

   size_t ndx = 0;
   for(size_t a = order-1; a < n; a++)
   {
      double t_range = knots(a + 1) - knots(a);
      double dt = t_range / static_cast<double>(samples);

      for(size_t it = 0; it < samples; ++it)
      {
         double t = knots(a) + dt * static_cast<double>(it);
         //printf("\tt: %f\n", t);
         lux::Vector pt;

         for(size_t i = 0; i < n; i++)
         {
            pt = pt + cp[i] * NormalizedBasis(i, order, t);
         }

         //printf("\tp: %f %f %f\n", pt[0], pt[1], pt[2]);
         curve.push_back(pt);
         ndx++;
      }
   }
}

/* ~ Curve classes ~ */
class Curve
{
   public:

      Curve()
      {
         glyph_radius = 0.025;
         glyph_segs = 16;

         double step = (2.0 * M_PI) / static_cast<double>(glyph_segs);
         point_glyph.resize(glyph_segs);
         for(int i = 0; i < glyph_segs; i++)
         {
            double angle = (step * i);
            point_glyph[i] = lux::Vector(glyph_radius * std::cos(angle), glyph_radius * std::sin(angle), 0.0);
         }

         order = 1;

         closed = false;
      };

      virtual ~Curve() {};

      virtual void generateCurve() {};
      virtual void regenerateCurve() {};
      virtual void drawControlPolygon(lux::Vector right, lux::Vector up) {};

      void drawCurve()
      {
         glBegin(GL_LINE_STRIP);
            for_each(curve.begin(), curve.end(), [] (lux::Vector pt)
            {
               glVertex3f(pt[0], pt[1], pt[2]);
            } );
         glEnd();

      }

      bool testCurve(lux::Vector r_pos, lux::Vector r_dir, int& x)
      {
         double radius = glyph_radius;
         double best_t = INT_MAX;
         double best_ndx = -1;
         size_t start = (closed) ? order : 0;

         for(size_t i = start; i < control_polygon.size(); i++)
         {
            lux::Vector pt = control_polygon[i];

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
                  best_ndx = i;
               }
         }

         x = best_ndx;
         if(x > -1)
            return true;

         return false;
      }

      void clearSelected() { selected = -1; };
      void setSelected(int s) { selected = s; };
      lux::Vector getSelected() { return control_polygon[selected]; }

      int closeCurve()
      { 
         if(closed) return 0;

         for(int i = 0; i < order; i++)
            control_polygon.push_back(control_polygon[i]);

         closed = true;

         return 1;
      }

      int openCurve()
      {
         if(!closed) return 0;

         for(int i = 0; i < order; i++)
            control_polygon.pop_back();

         closed = false;

         return 1;
      }

      void updateSamples(int s) { samples = s; };
      void updateSubdivs(int s) { subdivs = s; };
      virtual void updateOrder(int o) { order = 1; };

      int getOrder() { return order; }
      int getSamples() { return samples; }
      int getSubdivs() { return subdivs; }
      virtual int getType() { return NULL_CURVE_TYPE; }

      void setControlPoint(int i, lux::Vector P)
      {
         control_polygon[i] = P;
         if(closed)
         {
            int N_cp = control_polygon.size();
            if(i < order)
               control_polygon[N_cp - (order - i)] = P;
            else if(i > (N_cp - order) - 1)
               control_polygon[i - (N_cp - order)] = P;
         }
      }

      lux::Vector getControlPoint(int i) { return control_polygon[i]; }

      size_t getCurveLength() { return curve.size(); }
      lux::Vector getCurvePoint(int i)
      {
         if(i > -1 && i < curve.size())
            return curve[i];
         else
            return lux::Vector();
      }

      lux::Vector centroid()
      {
         lux::Vector cog(0, 0, 0);
         for(size_t i = 0; i < control_polygon.size(); i++)
            cog = cog + control_polygon[i];
         return cog * (1.0 / static_cast<double>(control_polygon.size()));
      }

      std::vector<lux::Vector> control_polygon;

   protected:

      std::vector<lux::Vector> curve;
      int samples;
      int subdivs;
      int order;

      int glyph_segs;
      double glyph_radius;
      std::vector<lux::Vector> point_glyph;

      int selected;
      bool closed;
};

class Bezier_DeCasteljauSampling : public Curve
{
   public:

      Bezier_DeCasteljauSampling() {};
      Bezier_DeCasteljauSampling(std::vector<lux::Vector> cp, int sa, int su)
      {
         control_polygon = cp;
         samples = sa;
         subdivs = su;

         selected = -1;
      };

      ~Bezier_DeCasteljauSampling()
      {
         control_polygon.clear();
         curve.clear();
      }

      void generateCurve()
      {
         curve.resize(samples);

         int i = 0;
         double s = static_cast<double>(samples-1);
         for(; i < samples-1; i++)
         {
            double t = static_cast<double>(i) / s;
            curve[i] = deCasteljau(control_polygon, t);
         }
         curve[i] = deCasteljau(control_polygon, 1.0);
      }

      void regenerateCurve()
      {
         generateCurve();
      }

      lux::Vector deCasteljau(std::vector<lux::Vector> P, double t)
      {
         //printf("\tSize of P: %d\n", P.size());
         //
         int n = P.size();
         while(n > 0)
         {
            n -= 1;
            for(int i = 0; i < n; i++)
               P[i] = (1.0 - t) * P[i] + t * P[i+1];
            //printf("\t\tPoint: <%.4e, %.4e, %.4e>\n", Q[i][0], Q[i][1], Q[i][2]);

         }

         return P[0];
      }

      void drawControlPolygon(lux::Vector right, lux::Vector up)
      {
         int cp_size = control_polygon.size();
         for(int i = 0; i < cp_size; i++)
         {
            lux::Vector pt = control_polygon[i];

            if(i == selected)
               glColor4f(0.0, 0.0, 1.0, 1.0);
            else
               glColor4f(1.0, 0.0, 0.0, 1.0);

            glBegin(GL_POINTS);
               glVertex3f(pt[0], pt[1], pt[2]);
            glEnd();

            glBegin(GL_LINE_LOOP);
               for(size_t j = 0; j < point_glyph.size(); j++)
               {
                  lux::Vector X = pt
                                + right * point_glyph[j][0]
                                + up * point_glyph[j][1];

                  glVertex3f(X[0], X[1], X[2]);
               }
            glEnd();

         }

         glColor4f(0.35, 0.35, 0.35, 1.0);
         glBegin(GL_LINE_STRIP);
            for_each(control_polygon.begin(), control_polygon.end(), [] (lux::Vector pt)
            {
               glVertex3f(pt[0], pt[1], pt[2]);
            } );
         glEnd();
      }

      int getType() { return BEZIER_SAMPLED; }
};

class Bezier_DeCasteljauSubdiv : public Curve
{
   public:

      Bezier_DeCasteljauSubdiv() {};
      Bezier_DeCasteljauSubdiv(std::vector<lux::Vector> cp, int sa, int su)
      {
         control_polygon = cp;
         samples = sa;
         subdivs = su;

         selected = -1;
      };

      ~Bezier_DeCasteljauSubdiv()
      {
         control_polygon.clear();
         curve.clear();
      }

      void generateCurve()
      {
         double u = 0.5;
         curve = Subdivide(control_polygon, subdivs, u);
      }

      void regenerateCurve()
      {
         generateCurve();
      }

      std::vector<lux::Vector> Subdivide(std::vector<lux::Vector> P, int n, double u)
      {
         std::vector<lux::Vector> L;
         std::vector<lux::Vector> R;

         if(n == 1)
            return P;

         size_t N = P.size();
         while(N > 1)
         {
            L.push_back(P[0]); // L = L . P_0

            R.push_back(P[N-1]); // R = R . P_N

            N -= 1;
            for(size_t i = 0; i < N; i++)
               P[i] = (1.0 - u) * P[i] + u * P[i+1];
         }
         L.push_back(P[0]);
         R.push_back(P[0]);

         std::vector<lux::Vector> concat;

         L = Subdivide(L, n-1, u);
         for_each(L.begin(), L.end(), [&] (lux::Vector X)
         {
            concat.push_back(X);
         } );

         concat.push_back(P[0]);

         R = Subdivide(R, n-1, u);
         for_each(R.rbegin(), R.rend(), [&] (lux::Vector X)
         {
            concat.push_back(X);
         } );

         return concat;
         //printf("\t\tPoint: <%.4e, %.4e, %.4e>\n", Q[i][0], Q[i][1], Q[i][2]);
      }

      void drawControlPolygon(lux::Vector right, lux::Vector up)
      {
         int cp_size = control_polygon.size();
         for(int i = 0; i < cp_size; i++)
         {
            lux::Vector pt = control_polygon[i];

            if(i == selected)
               glColor4f(0.0, 0.0, 1.0, 1.0);
            else
               glColor4f(0.0, 1.0, 0.0, 1.0);

            glBegin(GL_POINTS);
               glVertex3f(pt[0], pt[1], pt[2]);
            glEnd();

            glBegin(GL_LINE_LOOP);
               for(size_t j = 0; j < point_glyph.size(); j++)
               {
                  lux::Vector X = pt
                                + right * point_glyph[j][0]
                                + up * point_glyph[j][1];

                  glVertex3f(X[0], X[1], X[2]);
               }
            glEnd();

         }

         glColor4f(0.35, 0.35, 0.35, 1.0);
         glBegin(GL_LINE_STRIP);
            for_each(control_polygon.begin(), control_polygon.end(), [] (lux::Vector pt)
            {
               glVertex3f(pt[0], pt[1], pt[2]);
            } );
         glEnd();
      }

      int getType() { return BEZIER_SUBDIV; }

};

class UniformBSpline : public Curve
{
   public:

      UniformBSpline() {};
      UniformBSpline(std::vector<lux::Vector> cp, int o, int sa, int su)
      {
         control_polygon = cp;
         samples = sa;
         subdivs = su;
         order = o;

         selected = -1;
      };

      ~UniformBSpline()
      {
         control_polygon.clear();
         curve.clear();
      }

      void updateOrder(int o)
      {
         int ret = openCurve();
         order = o;
         if( ret )
            closeCurve();
      }

      double knots(int i) { return static_cast<double>(i); }

      void generateCurve()
      {
         size_t n = control_polygon.size();
         curve.clear();

         size_t ndx = 0;
         for(size_t a = order-1; a < n; a++)
         {
            double t_range = knots(a + 1) - knots(a);
            double dt = t_range / static_cast<double>(samples);

            for(double t = knots(a); t <= knots(a + 1); t += dt)
            {
               //printf("\tt: %f\n", t);
               lux::Vector pt;

               for(size_t i = 0; i < n; i++)
               {
                  pt = pt + control_polygon[i] * NormalizedBasis(i, order, t);
               }

               //printf("\tp: %f %f %f\n", pt[0], pt[1], pt[2]);
               curve.push_back(pt);
               ndx++;
            }
         }
      }

      void regenerateCurve()
      {
         generateCurve();
      }

      double NormalizedBasis(int i, int d, double t)
      {
         if(d == 1)
         {
            if(knots(i) <= t && t < knots(i+1))
               return 1;

            return 0;
         }

         double A = (t - knots(i))
                  / (knots(i+d-1) - knots(i));

         double B = (knots(i+d) - t) 
                  / (knots(i+d) - knots(i+1));

         double ret = A * NormalizedBasis(i, d-1, t) + B * NormalizedBasis(i+1, d-1, t);

         return ret;
      }

      void drawControlPolygon(lux::Vector right, lux::Vector up)
      {
         int cp_size = control_polygon.size();
         for(int i = 0; i < cp_size; i++)
         {
            lux::Vector pt = control_polygon[i];

            if(i == selected)
               glColor4f(0.0, 0.0, 1.0, 1.0);
            else
               glColor4f(0.0, 1.0, 1.0, 1.0);

            glBegin(GL_POINTS);
               glVertex3f(pt[0], pt[1], pt[2]);
            glEnd();

            glBegin(GL_LINE_LOOP);
               for(size_t j = 0; j < point_glyph.size(); j++)
               {
                  lux::Vector X = pt
                                + right * point_glyph[j][0]
                                + up * point_glyph[j][1];

                  glVertex3f(X[0], X[1], X[2]);
               }
            glEnd();

         }

         glColor4f(0.35, 0.35, 0.35, 1.0);
         glBegin(GL_LINE_STRIP);
            for_each(control_polygon.begin(), control_polygon.end(), [] (lux::Vector pt)
            {
               glVertex3f(pt[0], pt[1], pt[2]);
            } );
         glEnd();
      }

      int getType() { return UNIFORM_3D_BSPLINE; }

};

}

} 
#endif

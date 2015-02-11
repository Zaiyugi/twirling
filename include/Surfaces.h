/* Zachary Shore
 * CPSC-881-001
 * Created: 2014-02-21
 * Last Edit: 2014-03-05
 * Geometric Modeling
 * Surface classes
 */
#ifndef __SURFACES_H__
#define __SURFACES_H__
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <string>

#include <vector>
#include <set>
#include <algorithm>

// OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

// Personal
#include "opengl/ogl_include.h"
#include "Vector.h"
#include "Matrix.h"
#include "LinearAlgebra.h"

#include "Curves.h"
#include "ControlPolytope.h"

namespace cagd
{

struct MyTraits : public OpenMesh::DefaultTraits
{
   VertexAttributes(    OpenMesh::Attributes::Normal  | OpenMesh::Attributes::Status );
   FaceAttributes(      OpenMesh::Attributes::Normal  | OpenMesh::Attributes::Status );
   HalfedgeAttributes(  OpenMesh::Attributes::Normal  | OpenMesh::Attributes::Status );
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits> MyMesh;

typedef OpenMesh::VPropHandleT< MyMesh::VertexHandle > VProp_VH;

typedef OpenMesh::HPropHandleT< MyMesh::VertexHandle > HProp_VH;

typedef OpenMesh::FPropHandleT< std::vector<MyMesh::VertexHandle> > FProp_VH;
typedef OpenMesh::FPropHandleT< MyMesh::Point > FProp_Centroid;

namespace surf
{

const unsigned int NULL_SURFACE = 0;

/* SURFACE CONSTRUCTIONS */
const unsigned int EXTRUSION = 1;
const unsigned int REVOLUTION = 2;
const unsigned int SWEEP = 3;
const unsigned int LOFT = 4;

/* SURFACE TYPES */
const unsigned int DOOSABIN = 3;
const unsigned int CATMULLCLARK = 4;
const unsigned int LOOP = 5;

/* SURFACE QUALIFIER */
const unsigned int DIRECT = 1;
const unsigned int NUMERIC = 2;
const unsigned int SUBDIV = 3;

/* EXTRUSION TYPES */
const unsigned int EXTRUDE_ALONG_NORMAL = 1;
const unsigned int EXTRUDE_ALONG_PATH = 2;
const unsigned int EXTRUDE_REGION = 3;

// -------------------------------------
typedef struct HornParameters
{
   int steps;

   double hn;
   double sn;
   double theta;
} horn_t;

class param_t
{
   public:
      param_t() {} 

      std::string name;

      int construct;
      int curve_type;
      int surf_type;
      int extru_type;

      lux::Vector axis;
      lux::Vector orig;
      double distance;
      double rotation;

      int slices;

      int samples;
      int subdivs;
      int order;
      bool closed;

      lux::Vector color;

      lux::Matrix rot;

      bool fromFile;
      bool triangulated;

      horn_t horn_params;

      void copy(param_t& src)
      {
         this->construct = src.construct;
         this->curve_type = src.curve_type;
         this->surf_type = src.surf_type;

         this->axis = src.axis;
         this->orig = src.orig;

         this->distance = src.distance;
         this->rotation = src.rotation;

         this->slices = src.slices;

         this->samples = src.samples;
         this->subdivs = src.subdivs;
         this->order = src.order;

         // Update Horn Properties
         this->extru_type = src.extru_type;
         this->horn_params.steps = src.horn_params.steps;
         this->horn_params.hn = src.horn_params.hn;
         this->horn_params.sn = src.horn_params.sn;
         this->horn_params.theta = src.horn_params.theta;
      }

};

class SurfaceState
{
   public:
      SurfaceState() {};
      SurfaceState(MyMesh _m, param_t _p) : M(_m), P(_p) {};

      MyMesh M;
      param_t P;
};

void buildParametersFromCurve(param_t& P, cagd::curves::Curve* C)
{
   P.curve_type = C->getType();
   P.samples = C->getSamples();
   P.subdivs = C->getSubdivs();
   P.order = C->getOrder();
}

lux::Matrix rot_mat(lux::Vector axis, double angle)
{
   lux::Matrix op;
   lux::outer_product(axis, axis, op);
   lux::Matrix cpm(0,       -axis[2], axis[1],
                   axis[2],  0,      -axis[0],
                  -axis[1],  axis[0], 0);

   lux::Matrix iden(1, 0, 0, 
                    0, 1, 0, 
                    0, 0, 1);

   lux::Matrix rot = iden * std::cos(angle)
                   + std::sin(angle) * cpm
                   + (1.0 - std::cos(angle)) * op;

   return rot;
}

// -------------------------------------

// -------------------------------------
void DooSabin_Subdivide(MyMesh& mesh, param_t& props)
{

   MyMesh mesh2;

   FProp_Centroid centroid;
   HProp_VH heh_vertices;
   mesh.add_property(centroid);
   mesh.add_property(heh_vertices);

   MyMesh::FaceIter f_it, f_end(mesh.faces_end());
   MyMesh::FaceVertexIter fv_it;

   // Centroids
   for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
   {
      mesh.property(centroid, *f_it).vectorize(0.0f);
      float valence = 0.0f;

      for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
      {
         mesh.property(centroid, *f_it) += mesh.point(*fv_it);
         ++valence;
      }

      mesh.property(centroid, *f_it) /= valence;

      MyMesh::Point C = mesh.property(centroid, *f_it);
   }

   MyMesh::FaceHandle test_fh;
   std::vector<MyMesh::VertexHandle> vhandles;

   // Loop through all faces of mesh
   // Iterate over halfedges, compute required vertices
   // Build face faces
   for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
   {
      MyMesh::HalfedgeHandle heh_init = mesh.halfedge_handle(*f_it);
      MyMesh::HalfedgeHandle heh_cur, heh_next;
      MyMesh::Point V, Ve1, Ve2, Vf;

      heh_cur = heh_init;
      heh_next = mesh.prev_halfedge_handle(heh_cur);

      MyMesh::VertexHandle vh0 = mesh.to_vertex_handle(heh_cur);
      MyMesh::VertexHandle vh1 = mesh.from_vertex_handle(heh_next);

      V = mesh.point(mesh.from_vertex_handle(heh_cur));
      Ve1 = (mesh.point(vh0) + V) / 2.0;
      Ve2 = (mesh.point(vh1) + V) / 2.0;
      Vf = (V + Ve1 + Ve2 + mesh.property(centroid, *f_it)) / 4.0;

      vhandles.push_back(mesh2.add_vertex(Vf));
      mesh.property(heh_vertices, heh_cur) = vhandles.back();

      heh_cur = heh_next;
      while(heh_cur != heh_init)
      {
         heh_next = mesh.prev_halfedge_handle(heh_cur);

         vh0 = mesh.to_vertex_handle(heh_cur);
         vh1 = mesh.from_vertex_handle(heh_next);

         V = mesh.point(mesh.from_vertex_handle(heh_cur));
         Ve1 = (mesh.point(vh0) + V) / 2.0;
         Ve2 = (mesh.point(vh1) + V) / 2.0;
         Vf = (V + Ve1 + Ve2 + mesh.property(centroid, *f_it)) / 4.0;

         vhandles.push_back(mesh2.add_vertex(Vf));
         mesh.property(heh_vertices, heh_cur) = vhandles.back();

         heh_cur = heh_next;
      }

      std::reverse(vhandles.begin(), vhandles.end());
      mesh2.add_face(vhandles);
      vhandles.clear();
   }

   // Build Edge faces
   MyMesh::EdgeIter e_it, e_end(mesh.edges_end());
   for(e_it = mesh.edges_begin(); e_it != e_end; ++e_it)
   {
      MyMesh::HalfedgeHandle heh0 = mesh.halfedge_handle(*e_it, 0);
      MyMesh::HalfedgeHandle heh1 = mesh.halfedge_handle(*e_it, 1);
      MyMesh::FaceHandle fh0 = mesh.face_handle(heh0);
      MyMesh::FaceHandle fh1 = mesh.face_handle(heh1);

      if(fh0 == MyMesh::InvalidFaceHandle || fh1 == MyMesh::InvalidFaceHandle)
         continue;

      vhandles.push_back( mesh.property(heh_vertices, heh0 ) );
      vhandles.push_back( mesh.property(heh_vertices, mesh.next_halfedge_handle(heh0) ) );

      vhandles.push_back( mesh.property(heh_vertices, heh1 ) );
      vhandles.push_back( mesh.property(heh_vertices, mesh.next_halfedge_handle(heh1) ) );

      std::reverse(vhandles.begin(), vhandles.end());
      test_fh = mesh2.add_face(vhandles);

      vhandles.clear();
   }

   // Build Vertex faces
   MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
   for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
   {
      if(mesh.is_boundary(*v_it))
         continue;

      MyMesh::VertexOHalfedgeIter voh_it;
      for(voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); --voh_it)
         vhandles.push_back( mesh.property(heh_vertices, *voh_it) );

      test_fh = mesh2.add_face(vhandles);

      vhandles.clear();
   }

   // Clear Mesh, deleting all including properties
   // Assign new mesh
   mesh.clear();
   mesh = mesh2;

}

void extrudeSelectedFaces(MyMesh& mesh, std::vector<bool>& selected, horn_t& H)
{
   std::printf("Create Horns!\n");

   std::vector<MyMesh::FaceHandle> next_fh;
   VProp_VH mesh2_vh;
   mesh.add_property(mesh2_vh);

   MyMesh mesh2;

   MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
   for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
      mesh.property(mesh2_vh, *v_it) = mesh2.add_vertex(mesh.point(*v_it));

   MyMesh::FaceIter f_it, f_end(mesh.faces_end());
   MyMesh::FaceVertexIter fv_it;

   for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
   {
      std::printf("FACE IDX: %d\n", f_it->idx());
      size_t f_val = mesh.valence(*f_it);
      size_t i = 0;

      MyMesh::Point Cen(0,0,0);
      std::vector<MyMesh::VertexHandle> local, extr_local;
      local.resize(f_val);

      for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
      {
         Cen += mesh.point(*fv_it);
         local[i++] = mesh.property(mesh2_vh, *fv_it);
      }
      Cen /= static_cast<float>(f_val);

      // If this face is not selected, add to mesh2
      if( !selected.at(f_it->idx()) )
      {
         mesh2.add_face(local);
         continue;
      }
      // If this face is selected, extrude it

      MyMesh::Normal fn = mesh.normal(*f_it);

      extr_local.resize(f_val);

      lux::Matrix rotM = rot_mat(lux::Vector(fn[0], fn[1], fn[2]), H.theta * M_PI / 180.0);
      for(i = 0; i < local.size(); ++i)
      {
         MyMesh::Point P = mesh2.point(local[i]);

         // Translate to origin
         MyMesh::Point tP = P - Cen;
         // Scale
         tP = tP * H.sn;
         // Rotate about face normal
         lux::Vector rotP = rotM * lux::Vector(tP[0], tP[1], tP[2]);
         // Translate back, with the addition extrude along normal
         MyMesh::Point Pf = MyMesh::Point(rotP[0], rotP[1], rotP[2]) + Cen + fn * H.hn;

         extr_local[i] = mesh2.add_vertex(Pf);
      }

      for(i = 0; i < extr_local.size(); ++i)
      {
         MyMesh::VertexHandle vh0 = extr_local[i];
         MyMesh::VertexHandle vh1 = extr_local[(i+1)%f_val];
         MyMesh::VertexHandle vh2 = local[i];
         MyMesh::VertexHandle vh3 = local[(i+1)%f_val];

         mesh2.add_face(vh2, vh3, vh1, vh0);

      }

      next_fh.push_back(mesh2.add_face(extr_local));
   }

   mesh = mesh2;
   selected.clear();
   selected.resize(mesh.n_faces());

   for(size_t i = 0; i < next_fh.size(); ++i)
      selected.at(next_fh[i].idx()) = true;

   printf("EXIT HORNS\n");
}

void extrudeSelectedFacesAlongPath(MyMesh& mesh, std::vector<bool>& selected, std::vector<lux::Vector>& path, horn_t& H)
{
   std::printf("Create Horns along Path!\n");

   std::vector<MyMesh::FaceHandle> next_fh;
   VProp_VH mesh2_vh;
   mesh.add_property(mesh2_vh);

   MyMesh mesh2;

   MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
   for(v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
      mesh.property(mesh2_vh, *v_it) = mesh2.add_vertex(mesh.point(*v_it));

   lux::Vector path_cen;
   for(size_t i = 0; i < path.size(); i++)
      path_cen += path[i];
   path_cen /= static_cast<double>(path.size());

   MyMesh::FaceIter f_it, f_end(mesh.faces_end());
   MyMesh::FaceVertexIter fv_it;

   for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
   {
      std::printf("FACE IDX: %d\n", f_it->idx());
      size_t f_val = mesh.valence(*f_it);
      size_t i = 0;

      MyMesh::Point Cen(0,0,0);
      std::vector<MyMesh::VertexHandle> local, extr_local;
      local.resize(f_val);

      for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
      {
         Cen += mesh.point(*fv_it);
         local[i++] = mesh.property(mesh2_vh, *fv_it);
      }
      Cen /= static_cast<float>(f_val);

      // If this face is not selected, add to mesh2
      if( !selected.at(f_it->idx()) )
      {
         mesh2.add_face(local);
         continue;
      }
      // If this face is selected, extrude it

      MyMesh::Normal fn = mesh.normal(*f_it);
      MyMesh::Point ft = mesh.point(*fv_it) - Cen;

      lux::Matrix rotT, rotN, rotM;
      lux::Vector T(ft[0], ft[1], ft[2]);
      lux::Vector N(fn[0], fn[1], fn[2]);
      lux::Vector face_cen(Cen[0], Cen[1], Cen[2]);

      T = T.unitvector();
      N = N.unitvector();

      // Generate local coord frames for path
      // Initial frame: local coord frame for face
      cagd::curves::Frame U_face(T, N ^ T, N);
      lux::Vector target = (path[1] - path[0]).unitvector();
      lux::Vector rot_axis = (U_face.t ^ target).unitvector();
      double rot_angle = std::acos(U_face.t * target);
      rotT = rot_mat(rot_axis, rot_angle);

      std::vector<lux::Vector> local_path; local_path.resize(path.size());
      for(size_t k = 0; k < path.size(); ++k)
      {
         lux::Vector P = path[k] - path[0];
         P = P * rotT;
         local_path[k] = P + face_cen;
      }

      std::vector<cagd::curves::Frame> local_frames;
      cagd::curves::generateFrame(local_path, U_face, local_frames);

      for(size_t k = 1; k < local_frames.size(); ++k)
      {
         // Calculate rotation to tangent
         rot_axis = local_frames[k].t ^ U_face.t;
         if(local_frames[k].t * U_face.t >= 1.0)
            rot_angle = 0.0;
         else
            rot_angle = std::acos(local_frames[k].t * U_face.t);
         rotT = rot_mat(rot_axis.unitvector(), rot_angle);
         if(rot_axis.magnitude() > 0.0)
            U_face.transform(rotT);

         // Calculate rotation to reference
         rot_axis = local_frames[k].r ^ U_face.r;
         if(local_frames[k].r * U_face.r >= 1.0)
            rot_angle = 0.0;
         else
            rot_angle = std::acos(local_frames[k].r * U_face.r);
         rotN = rot_mat(rot_axis.unitvector(), rot_angle);
         if(rot_axis.magnitude() > 0.0)
            U_face.transform(rotN);
         
         extr_local.resize(f_val);
         face_cen = local_path[k-1];

         rotM = rot_mat(local_frames[k-1].t, H.theta * M_PI / 180.0);
         for(size_t i = 0; i < local.size(); ++i)
         {
            MyMesh::Point P = mesh2.point(local[i]);

            // Translate to origin
            MyMesh::Point tP = P - MyMesh::Point(face_cen[0], face_cen[1], face_cen[2]);
            // Scale
            tP = tP * H.sn;
            // Rotate to tangent, reference, then about tangent
            lux::Vector rotP(tP[0], tP[1], tP[2]);
            rotP = rotP * rotT;
            rotP = rotP * rotN;
            rotP = rotP * rotM;
            // Translate to new position
            lux::Vector ntP = rotP + local_path[k];
            MyMesh::Point Pf = MyMesh::Point(ntP[0], ntP[1], ntP[2]);

            extr_local[i] = mesh2.add_vertex(Pf);
         }

         for(size_t i = 0; i < extr_local.size(); ++i)
         {
            MyMesh::VertexHandle vh0 = extr_local[i];
            MyMesh::VertexHandle vh1 = extr_local[(i+1)%f_val];
            MyMesh::VertexHandle vh2 = local[i];
            MyMesh::VertexHandle vh3 = local[(i+1)%f_val];

            //std::printf("Create face from: %d %d %d %d\n", vh2, vh3, vh1, vh0);
            mesh2.add_face(vh2, vh3, vh1, vh0);
            //std::printf("Face %d added\n", i);
         }

         local = extr_local;
      }

      next_fh.push_back(mesh2.add_face(extr_local));
   }

   mesh = mesh2;
   selected.clear();
   selected.resize(mesh.n_faces());

   for(size_t i = 0; i < next_fh.size(); ++i)
      selected.at(next_fh[i].idx()) = true;

   printf("EXIT HORNS\n");
}

class Surface
{
   public:

      Surface() : existence(false) {}

      Surface(param_t& P)
      {
         properties = P;
         properties.fromFile = false;
         existence = false;
      }

      ~Surface() {};

      void initialize()
      {
         mesh.request_vertex_normals();
         mesh.request_halfedge_normals();
         mesh.request_face_normals();

         mesh.update_normals();

         selected.resize(mesh.n_faces());
         selectNone();
      }

      void AddSubdivision(param_t& surface_params)
      {
         properties.copy(surface_params);

         if(properties.subdivs == 0)
         {
            std::printf("No subdivisions specificed\n");
            return;
         }

         mesh.update_normals();
         for(int i = 0; i < properties.subdivs; i++)
         {
            DooSabin_Subdivide(mesh, properties);
         }

         selected.resize(mesh.n_faces());
         selectNone();
         mesh.update_normals();
         history.push_back(new SurfaceState(mesh, properties));
         std::printf("Subdivided; added to history\n");
      }

      void AddExtrusion(param_t& surface_params)
      {
         properties.copy(surface_params);

         size_t i = 0;
         while(i < selected.size() && !selected[i]) ++i;
         if(i == selected.size())
         {
            std::printf("Nothing selected\n");
            return;
         }

         if(properties.extru_type == EXTRUDE_ALONG_PATH)
         {
            extrudeSelectedFacesAlongPath(mesh, selected, path, properties.horn_params);
         } else if(properties.horn_params.steps > 0) {
            if(properties.extru_type == EXTRUDE_ALONG_NORMAL)
            {
               std::printf("EXTRUDE SELECTED FACES!\n");
               int cur = properties.horn_params.steps;
               for(int i = 0; i < std::abs(cur); ++i)
               {
                  mesh.update_normals();
                  extrudeSelectedFaces(mesh, selected, properties.horn_params);
               }

            }

         }

         mesh.update_normals();
         history.push_back(new SurfaceState(mesh, properties));
         std::printf("Extrusion done; added to history\n");
      }

      void Undo()
      {
         if(history.size() > 1)
         {
            SurfaceState* recent = history.back();
            history.pop_back();
            delete recent;

            mesh = history.back()->M;
            properties.copy(history.back()->P);
            mesh.update_normals();

            selected.resize(mesh.n_faces());
            selectNone();
         } else {
            std::printf("History only contains the original mesh\n");
         }

         std::printf("UNDO complete\n");
      }

      void selectNone()
      {
         for(MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
            selected.at(f_it->idx()) = false;   
      }

      void selectAll()
      {
         for(MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
            selected.at(f_it->idx()) = true;
      }

		void selectAllTriangles()
		{
         for(MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
				if(mesh.valence(*f_it) == 3)
            	selected.at(f_it->idx()) = true;
		}

      void assignPath(std::vector<lux::Vector>& pts)
      {
         path = pts;

         std::printf("Path assigned to surface\n");
      }

      void generate()
      {
         std::printf("Exit Generate!\n");
      }

      void regenerate()
      {
         generate();
      }

      void draw()
      {
         glColor4f(1.0, 1.0, 1.0, 1.0);

         MyMesh::FaceIter f_it, f_end(mesh.faces_end());
         MyMesh::FaceVertexIter fv_it;
         MyMesh::Point pt;
         MyMesh::Normal n;

         for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
         {
            bool useNormalColoring = true;
            if( selected.at(f_it->idx()) )
               useNormalColoring = false;

            // Surface
            glEnable(GL_POLYGON_OFFSET_FILL);
            glPolygonOffset(1.0, 1.0);

            MyMesh::Point cen(0,0,0);
            int val = 0;
            glBegin(GL_POLYGON);

            for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
            {
               pt = mesh.point(*fv_it);
               n = mesh.normal(*fv_it);
               cen += pt;
               val++;

               float r = (n[0] + 1.0) * 0.5;
               float g = (n[1] + 1.0) * 0.5;
               float b = (n[2] + 1.0) * 0.5;
               float alpha = 0.5;
               if(useNormalColoring)
                  alpha = 1.0;
               
               glColor4f(r * alpha, g * alpha, b * alpha, alpha);
               glNormal3fv(&n[0]);
               glVertex3fv(&pt[0]);
            }
            glEnd();

            glPolygonOffset(0.0, 0.0);
            glDisable(GL_POLYGON_OFFSET_FILL);

            // Wireframe
            glColor4f(0.0, 0.0, 0.0, 1.0);
            glBegin(GL_LINE_LOOP);
            for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
            {
               pt = mesh.point(*fv_it);
               glVertex3f(pt[0], pt[1], pt[2]);
            }
            glEnd();

         }

         glColor4f(1.0, 1.0, 1.0, 1.0);

      }

      void select_draw()
      {
         MyMesh::FaceIter f_it, f_end(mesh.faces_end());
         MyMesh::FaceVertexIter fv_it;
         MyMesh::Point pt;

         for(f_it = mesh.faces_begin(); f_it != f_end; ++f_it)
         {
            glPushName(f_it->idx());
            glBegin(GL_POLYGON);
            for(fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
            {
               pt = mesh.point(*fv_it);
               glVertex3fv(&pt[0]);
            }
            glEnd();
            glPopName();
         }

      }

      bool select(ogl::OGL_Camera* cam, int x, int y, int w, int h)
      {
         GLint viewport[4];
         // BUFSIZE = 512
         GLuint selectBuf[512];

         // Initialize/setup Selection buffer and render mode
         glSelectBuffer(512, selectBuf);
         glRenderMode(GL_SELECT);

         glMatrixMode(GL_PROJECTION);
         glPushMatrix();
         glLoadIdentity();

         glGetIntegerv(GL_VIEWPORT, viewport);
         gluPickMatrix(x, viewport[3] - y, 2, 2, viewport);
         gluPerspective(cam->Fov, (float)w/(float)h, cam->NearPlane, cam->FarPlane);
         gluLookAt(
            cam->Pos.x, cam->Pos.y, cam->Pos.z, 
            cam->Aim.x, cam->Aim.y, cam->Aim.z,
            cam->Up.x, cam->Up.y, cam->Up.z
         );
         glMatrixMode(GL_MODELVIEW);
         glInitNames();

         // Draw all mesh items with appropriate names
         select_draw();

         // Gather number of hits
         // If 0, nothing selected
         int hits;
         glMatrixMode(GL_PROJECTION);
         glPopMatrix();
         glMatrixMode(GL_MODELVIEW);
         glFlush();

         hits = glRenderMode(GL_RENDER);
         if(hits == 0)
            return false;

         // Process hits
         std::printf("Hits: %d\n", hits);
         uint32_t *ptr = (uint32_t*)selectBuf;
         uint32_t smallest = 0xffffffff;
         int small_ndx = -1;

         for(int i = 0; i < hits; ++i)
         {
            uint32_t n_names = ptr[i * 4];
            uint32_t min = ptr[i * 4 + 1];
            uint32_t max = ptr[i * 4 + 2];
            uint32_t name = ptr[i * 4 + 3];

            if(min < smallest)
            {
               smallest = min;
               small_ndx = name;
            }

            std::printf("Number of names stored for hit record %d: %d\n", i, n_names);
            std::printf("\tMin depth: %u\n", min);
            std::printf("\tMax depth: %u\n", max);
            std::printf("\tName: %d\n", name);

            std::printf("\n");
         }

         if(small_ndx != -1)
            selected.at( mesh.face_handle(small_ndx).idx() ) = true;

         return true;
      }

      std::string getName() { return properties.name; }

      void updateSamples(int s) { properties.samples = s; };
      void updateSubdivs(int s) { properties.subdivs = s; };
      void updateOrder(int o)
      {
         properties.order = o;
      };

      void clearSelection()
      {
         std::replace(selected.begin(), selected.end(), true, false);
      }

      bool do_I_exist() { return existence; }

      friend int WriteSurface();
      friend int ReadSurface();

      param_t properties;

      bool existence;

      MyMesh mesh;
      std::vector<SurfaceState*> history;

      std::vector<bool> selected;

      /* Path and Local coord frame */
      std::vector<lux::Vector> path;

};

// -------------------------------------
int WriteSurface(Surface* S, std::string name)
{
   std::string filename = name;
   filename += ".obj";

   try
   {
      if(!OpenMesh::IO::write_mesh(S->mesh, filename.c_str()))
      {
         fprintf(stderr, "Cannot write mesh to %s\n", filename.c_str());
         return 0;
      } else {
         fprintf(stderr, "%s: Mesh write successful\n", filename.c_str());
      }
   } catch( std::exception& x) {
      fprintf(stderr, "%s\n", x.what());
      return 0;
   }

   return 1;
}

int ReadSurface(Surface* S, std::string name)
{
   std::string filename = name;
   filename += ".obj";

   try
   {
      if(!OpenMesh::IO::read_mesh(S->mesh, filename.c_str()))
      {
         fprintf(stderr, "Cannot read mesh to %s\n", filename.c_str());
         return 0;
      } else {
         fprintf(stderr, "%s: Mesh read successful\n", filename.c_str());
      }
   } catch( std::exception& x) {
      fprintf(stderr, "%s\n", x.what());
      return 0;
   }

   S->properties.fromFile = true;
   S->existence = true;
   S->history.clear();
   S->history.push_back(new SurfaceState(S->mesh, S->properties));
   S->initialize();

   return 1;
}
// -------------------------------------

}

}

#endif

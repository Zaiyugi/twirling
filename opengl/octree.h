#ifndef OCTREE
#define OCTREE

#include "global_prototypes.h"
#include "global_classes.h"

#define MAX_INDICES 4
#define MIN_INDICES 2
using namespace std;

extern vector<Sphere*> balls;
extern set<CollPairLocal*> collSet;

class Octree
{
   private:
      Octree* children[2][2][2];	// Pointers to child nodes
      Vec3<double> CoM; 			// Center of mass point
      Vec3<double> v0;			// Largest vertex
      Vec3<double> v7;			// Smallest vertex
      set<int> node;			// Contains the set of ball indices for this space
      bool hasChild;			// Does this node have children
      int nIndices;			// Number of ball indices starting from here
      int depth;			// Depth within (this) tree
      int max_depth;			// Max Depth for (this) tree
      double side_len;			// side length of node
    
      void fileIndex(int index, bool addIndex) 
      {
         //printf("ENTERED FILING...\n");
         Vec3<double> pos = balls[index]->center;

         double x, y, z;
         int i, j, k;
         x = pos.x;
         y = pos.y;
         z = pos.z;

         if(x <= CoM.x)
            i = 0;
         else
            i = 1;
    
         if(y <= CoM.y)
            j = 0;
         else
            j = 1;
    
         if(z <= CoM.z)
            k = 0;
         else
            k = 1;
  
         if(addIndex) 
            children[i][j][k]->add(index);
         else
            children[i][j][k]->remove(index);
      }

      void buildChild()
      {
         double sl2 = side_len / 2;

         /* -------------------------------------------------------------------------------------------- */
            children[0][0][0] = new Octree(Vec3<double>(v7.x, v7.y, v7.z), 
                                           Vec3<double>(CoM.x, CoM.y, CoM.z), 
                                           depth + 1, max_depth - 1, sl2);

            children[1][0][0] = new Octree(Vec3<double>(v7.x + sl2, v7.y, v7.z), 
                                           Vec3<double>(CoM.x + sl2, CoM.y, CoM.z), 
                                           depth + 1, max_depth - 1, sl2);

            children[0][0][1] = new Octree(Vec3<double>(v7.x, v7.y, v7.z + sl2), 
                                           Vec3<double>(CoM.x, CoM.y, CoM.z + sl2), 
                                           depth + 1, max_depth - 1, sl2);

            children[1][0][1] = new Octree(Vec3<double>(v7.x + sl2, v7.y, v7.z + sl2), 
                                           Vec3<double>(CoM.x + sl2, CoM.y, CoM.z + sl2), 
                                           depth + 1, max_depth - 1, sl2);
         /* -------------------------------------------------------------------------------------------- */

         /* -------------------------------------------------------------------------------------------- */
            children[1][1][1] = new Octree(Vec3<double>(CoM.x, CoM.y, CoM.z), 
                                           Vec3<double>(v0.x, v0.y, v0.z), 
                                           depth + 1, max_depth - 1, sl2);

            children[0][1][1] = new Octree(Vec3<double>(CoM.x - sl2, CoM.y, CoM.z), 
                                           Vec3<double>(v0.x - sl2, v0.y, v0.z), 
                                           depth + 1, max_depth - 1, sl2);

            children[1][1][0] = new Octree(Vec3<double>(CoM.x, CoM.y, CoM.z - sl2), 
                                           Vec3<double>(v0.x, v0.y, v0.z - sl2), 
                                           depth + 1, max_depth - 1, sl2);

            children[0][1][0] = new Octree(Vec3<double>(CoM.x - sl2, CoM.y, CoM.z - sl2), 
                                           Vec3<double>(v0.x - sl2, v0.y, v0.z - sl2), 
                                           depth + 1, max_depth - 1, sl2);
         /* -------------------------------------------------------------------------------------------- */

         set<int>::iterator del_it;
         for(set<int>::iterator it = node.begin(); it != node.end(); )
         {
            del_it = it;
            ++it;
            fileIndex((*del_it), true);
            node.erase((*del_it));
         }

         hasChild = true;
      }

      void collectIndices(set<int> &ns) 
      {
         if(hasChild)
         {
            for(int i = 0; i < 2; i++)
               for(int j = 0; j < 2; j++)
                  for(int k = 0; k < 2; k++)
                     children[i][j][k]->collectIndices(ns);
         } else {
            for(set<int>::iterator it = node.begin(); it != node.end(); it++)
               ns.insert((*it));
         }
      }
 
      void deconChild()
      {         
         collectIndices(node);
         for(int i = 0; i < 2; i++)
            for(int j = 0; j < 2; j++)
               for(int k = 0; k < 2; k++)
                  delete children[i][j][k];

         hasChild = false;
      }
 
   public:

      Octree()
      {
         hasChild = false;
         nIndices = 0;
         depth = 0;
         max_depth = 0;
      }

      Octree(Vec3<double> min, Vec3<double> max, int d, int md, int sl)
      {
         v7 = min;
         v0 = max;
         CoM = (v0 + v7) / 2.0;
         depth = d;
         max_depth = md;
         side_len = sl;

         hasChild = false;
         nIndices = 0;

         for(int i = 0; i < 2; i++)
            for(int j = 0; j < 2; j++)
               for(int k = 0; k < 2; k++)
                  children[i][j][k] = NULL;
      }

      ~Octree()
      {
         if(hasChild)
            deconChild();
      }

      void cntChildren(int &cnt)
      {
         bool isNull;
         if(hasChild)
         {
            for(int i = 0; i < 2; i++)
               for(int j = 0; j < 2; j++)
                  for(int k = 0; k < 2; k++)
                  {
                     children[i][j][k]->cntChildren(cnt);
                  }
         } else {
            cnt += node.size();
         }
      }

      void add(int index)
      {
         nIndices += 1;
         //printf("Node size: %d\n", (int)node.size());
         if(!hasChild && node.size() > MAX_INDICES - 1 && side_len > 1 && max_depth > -1) {
            //printf("NODE LIMITS -> MIN: %.2lf %.2lf %.2lf | MAX: %.2lf %.2lf %.2lf\n", 
            //       v7.x, v7.y, v7.z, v0.x, v0.y, v0.z);
            //printf("CALLING BUILDCHILD() AT DEPTH %d W/ MAXDEPTH %d\n", depth, max_depth);
            buildChild();
            //printf("FINISHED BUILDCHILD() AT DEPTH %d W/ MAXDEPTH %d\n", depth, max_depth);
         }

         if(hasChild)
            fileIndex(index, true);
         else if(node.size() < MAX_INDICES)
            node.insert(index);
      }

      void deleteIndex(int index)
      {
         nIndices--;

         if(hasChild && nIndices < MIN_INDICES)
         {
            deconChild();
         }

         if(hasChild)
            fileIndex(index, false);
         else
            node.erase(index);
      }

      void remove(int index)
      {
         deleteIndex(index);
      }

      void moveIndice(int index)
      {
         remove(index);
         removeCollPairs(index);

         add(index);
      }

      void updateCollPairs()
      {
         if(hasChild)
         {
            for(int i = 0; i < 2; i++)
               for(int j = 0; j < 2; j++)
                  for(int k = 0; k < 2; k++)
                     children[i][j][k]->updateCollPairs();

         } else if(node.size() > 1) {
            CollPairLocal* cp = new CollPairLocal;
            for(set<int>::iterator it = node.begin(); it != node.end(); it++)
            {
               for(set<int>::iterator it2 = node.begin(); it2 != node.end(); it2++)
               {
                  if((*it) > (*it2))
                  {
                     cp = new CollPairLocal;
                     cp->indexA = *it;
                     cp->indexB = *it2;
                     strcpy(cp->typeA, "Sphere");
                     strcpy(cp->typeB, "Sphere");
                     collSet.insert(cp);
                  }
               }
            }
         }
      }

      void removeCollPairs(int index)
      {
         set<CollPairLocal*>::iterator it = collSet.begin();
         set<CollPairLocal*>::iterator del_it;
         while(it != collSet.end())
         {
            del_it = it;
            ++it;
            if((*del_it)->indexA == index || (*del_it)->indexB == index)
               collSet.erase(*del_it);
         }

      }

      bool checkSameSpace(Vec3<double> ballPos, Vec3<double> oldPos)
      {
         if(oldPos <= v0 && oldPos >= v7)
         {
            if(hasChild)
            {
               for(int i = 0; i < 2; i++)
                  for(int j = 0; j < 2; j++)
                     for(int k = 0; k < 2; k++)
                     {
                        if(children[i][j][k]->checkSameSpace(ballPos, oldPos))
                        {
                           return true;
                        }
                     }
            } else {
               if(ballPos <= v0 && ballPos >= v7)
               {
                  return true;
               } else {
                  return false;
               }
            }
         }
         return false;
      }

      void drawStruct()
      {
         if(hasChild)
         {
            for(int i = 0; i < 2; i++)
               for(int j = 0; j < 2; j++)
                  for(int k = 0; k < 2; k++)
                     children[i][j][k]->drawStruct();
         }
         
         glMatrixMode(GL_MODELVIEW);
         glPointSize(2.0);
         glBegin(GL_POINTS);
            glVertex3f(CoM.x, CoM.y, CoM.z);
         glEnd();
         glPushMatrix();
         glTranslatef(CoM.x, CoM.y, CoM.z);
         glutWireCube(side_len);
         glTranslatef(-CoM.x, -CoM.y, -CoM.z);
         glPopMatrix();
      }

      void drawObjects()
      {
         if(hasChild)
         {
            for(int i = 0; i < 2; i++)
               for(int j = 0; j < 2; j++)
                  for(int k = 0; k < 2; k++)
                     children[i][j][k]->drawObjects();
         } else {
            for(set<int>::iterator it = node.begin(); it != node.end(); it++)
            {
               glColor4f(1.0, 0.0, 0.0, 1.0);
               balls[(*it)]->draw();
            }
         }
      }
};

#endif

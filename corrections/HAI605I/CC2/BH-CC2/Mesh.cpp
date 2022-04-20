// --------------------------------------------------------------------------
// gMini,
// a minimal Glut/OpenGL app to extend                              
//
// Copyright(C) 2007-2009                
// Tamy Boubekeur
//                                                                            
// All rights reserved.                                                       
//                                                                            
// This program is free software; you can redistribute it and/or modify       
// it under the terms of the GNU General Public License as published by       
// the Free Software Foundation; either version 2 of the License, or          
// (at your option) any later version.                                        
//                                                                            
// This program is distributed in the hope that it will be useful,            
// but WITHOUT ANY WARRANTY; without even the implied warranty of             
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              
// GNU General Public License (http://www.gnu.org/licenses/gpl.txt)           
// for more details.                                                          
//                                                                          
// --------------------------------------------------------------------------

#include "Mesh.h"
#include "Vec3D.h"
#include <algorithm>

using namespace std;

void Mesh::clear () {
    clearTopology ();
    clearGeometry ();
}

void Mesh::clearGeometry () {
    vertices.clear ();
}

void Mesh::clearTopology () {
    triangles.clear ();
}

void Mesh::computeTriangleNormals (vector<Vec3> & triangleNormals) {
    for (vector<Triangle>::const_iterator it = triangles.begin ();
         it != triangles.end ();
         it++) {
        Vec3 e01 (vertices[it->getVertex (1)].getPosition ()
                - vertices[it->getVertex (0)].getPosition ());
        Vec3 e02 (vertices[it->getVertex (2)].getPosition ()
                - vertices[it->getVertex (0)].getPosition ());
        Vec3 n (Vec3::crossProduct (e01, e02));
        n.normalize ();
        triangleNormals.push_back (n);
    }
}

void Mesh::recomputeSmoothVertexNormals (unsigned int normWeight) {
    vector<Vec3> triangleNormals;
    computeTriangleNormals (triangleNormals);
    for (vector<Vertex>::iterator it = vertices.begin ();
         it != vertices.end ();
         it++)
        it->setNormal (Vec3 (0.0, 0.0, 0.0));
    vector<Vec3>::iterator itNormal = triangleNormals.begin ();
    vector<Triangle>::iterator it = triangles.begin ();
    for ( ; it != triangles.end (); it++, itNormal++)
        for (unsigned int  j = 0; j < 3; j++) {
            Vertex & vj = vertices[it->getVertex (j)];
            float w = 1.0; // uniform weights
            Vec3 e0 = vertices[it->getVertex ((j+1)%3)].getPosition ()
                    - vj.getPosition ();
            Vec3 e1 = vertices[it->getVertex ((j+2)%3)].getPosition ()
                    - vj.getPosition ();
            if (normWeight == 1) { // area weight
                w = Vec3::crossProduct (e0, e1).getLength () / 2.0;
            } else if (normWeight == 2) { // angle weight
                e0.normalize ();
                e1.normalize ();
                w = (2.0 - (Vec3::dotProduct (e0, e1) + 1.0)) / 2.0;
            }
            if (w <= 0.0)
                continue;
            vj.setNormal (vj.getNormal () + (*itNormal) * w);
        }
    Vertex::normalizeNormals (vertices);
}

//one-ring of each vertex, i.e. a list of vertices with which it shares an edge
void Mesh::collectOneRing (vector<vector<unsigned int> > & oneRing) const {
    //Initialisation le vecetur de o_one_ring de la taille du vecteur vertices
    oneRing.resize (vertices.size ());
    //Parcours les triangles et ajout les voisins dans le 1-voisinage
    //Tous les points opposés dans le triangle sont reliés
    for (unsigned int i = 0; i < triangles.size (); i++) {
        const Triangle & ti = triangles[i];
        for (unsigned int j = 0; j < 3; j++) {
            unsigned int vj = ti.getVertex (j);
            for (unsigned int k = 1; k < 3; k++) {
                unsigned int vk = ti.getVertex ((j+k)%3);
                if (find (oneRing[vj].begin (), oneRing[vj].end (), vk) == oneRing[vj].end ())
                    oneRing[vj].push_back (vk);
            }
        }
    }
}

void Mesh::addNoise(){
    for( unsigned int i = 0 ; i < vertices.size() ; i ++ )
    {
        float factor = 0.03;
        const Vec3 & p = vertices[i].position;
        const Vec3 & n = vertices[i].normal;
        vertices[i].position = Vec3( p[0] + factor*((double)(rand()) / (double)(RAND_MAX))*n[0],
                p[1] + factor*((double)(rand()) / (double)(RAND_MAX))*n[1],
                p[2] + factor*((double)(rand()) / (double)(RAND_MAX))*n[2]);

    }

}

//Fonction à completer Exercice 1
void Mesh::smooth( float lambda ){
    //Utiliser la fonction collectOneRing pour récuperer le 1-voisinage

    //(1) (a) Initialisation d'un vecteur pour les positions moyennes avec des coordonnées nulles

    vector<vector<unsigned int>> ring;
    collectOneRing(ring);
    cout << lambda << endl;

    vector<Vec3> new_positions;
    new_positions.resize(vertices.size());

    for (int i = 0; i < vertices.size(); ++i) {
        new_positions[i] = Vec3(0, 0, 0);

        for (int j = 0; j < ring[i].size(); ++j) {
            new_positions[i] += vertices[ring[i][j]].position;
        }

        new_positions[i] /= ring[i].size();
        new_positions[i] -= vertices[i].position;
    }

    for (int i = 0; i < vertices.size(); ++i) {
        vertices[i].position += lambda * new_positions[i];
    }


    //Parcourir les sommets


    //Parcourir le vecteur des positions resultante et diviser par le nombre de voisins (taille du vecteur de OneRing du sommet courant)


    // Parcourir les sommet et mettre à jour les positions en faisant un déplacement vers le barycentre du 1-voisinage

        // (1) (b) Calcul du vecteur de déplacement (position du sommet courant  vers la position moyenne)

        // (2) Déplacer le sommet de factor dans la direction calculée (nouvelle position = position courante + factor * vecteur)


}

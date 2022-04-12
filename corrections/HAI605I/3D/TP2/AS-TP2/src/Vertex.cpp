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
   
#include "Vertex.h"

#include <cmath>
#include <algorithm>

#include "Vec3D.h"

using namespace std;

static const unsigned int SIZE_OF_VERTEX = 10;

ostream & operator<< (ostream & output, const Vertex & v) {
    output << v.getPos () << endl 
           << v.getNormal () << endl;
    return output;
}

void Vertex::interpolate (const Vertex & u, const Vertex & v, float alpha) {
    setPos (Vec3Df::interpolate (u.getPos (), v.getPos (), alpha));
    Vec3Df normal = Vec3Df::interpolate (u.getNormal (), v.getNormal (), alpha);
    normal.normalize ();
    setNormal (normal);
}

// ------------------------------------
// Static Members Methods.
// ------------------------------------

void Vertex::computeAveragePosAndRadius (std::vector<Vertex> & vertices, 
                                         Vec3Df & center, float & radius) {
    center = Vec3Df (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < vertices.size (); i++)
        center += vertices[i].getPos ();
    center /= float (vertices.size ());
    radius = 0.0;
    for (unsigned int i = 0; i < vertices.size (); i++) {
        float vDistance = Vec3Df::distance (center, vertices[i].getPos ());
        if (vDistance > radius)
            radius = vDistance;
    }
}

void Vertex::scaleToUnitBox (vector<Vertex> & vertices, 
                             Vec3Df & center, float & scaleToUnit) {
    computeAveragePosAndRadius (vertices, center, scaleToUnit);
    for (unsigned int i = 0; i < vertices.size (); i++) 
        vertices[i].setPos (Vec3Df::segment (center, vertices[i].getPos ()) / scaleToUnit);
}

void Vertex::normalizeNormals (vector<Vertex> & vertices) {
    for (std::vector<Vertex>::iterator it = vertices.begin (); 
         it != vertices.end (); 
         it++) {
        Vec3Df n = it->getNormal ();
        if (n != Vec3Df (0.0, 0.0, 0.0)) {
            n.normalize ();
            it->setNormal (n);
        }
    }    
}


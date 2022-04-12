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
   
#pragma once

#include <iostream>
#include <vector>

#include "Vec3D.h"

class Vertex {
public:
    inline Vertex () 
        : pos (Vec3Df (0.0,0.0,0.0)), normal (Vec3Df (0.0, 0.0, 1.0)) {}
    inline Vertex (const Vec3Df & pos) 
        : pos (pos), normal (Vec3Df (0.0, 0.0, 1.0)) {}
    inline Vertex (const Vec3Df & pos, const Vec3Df & normal) : pos (pos), normal (normal) {}
    inline Vertex (const Vertex & v) : pos (v.pos), normal (v.normal) {}
    inline virtual ~Vertex () {}
    inline Vertex & operator= (const Vertex & vertex) {
        pos = vertex.pos;
        normal = vertex.normal;
        return (*this);
    }
    inline const Vec3Df & getPos () const { return pos; }
    inline const Vec3Df & getNormal () const { return normal; }  
    inline void setPos (const Vec3Df & newPos) { pos = newPos; }
    inline void setNormal (const Vec3Df & newNormal) { normal = newNormal; }
    inline bool operator== (const Vertex & v) { return (v.pos == pos && v.normal == normal); }
    void interpolate (const Vertex & u, const Vertex & v, float alpha = 0.5);

    static void computeAveragePosAndRadius (std::vector<Vertex> & vertices, 
                                            Vec3Df & center, float & radius);
    static void scaleToUnitBox (std::vector<Vertex> & vertices, 
                                Vec3Df & center, float & scaleToUnitBox);
    static void normalizeNormals (std::vector<Vertex> & vertices);

private:
    Vec3Df pos;
    Vec3Df normal;
};

extern std::ostream & operator<< (std::ostream & output, const Vertex & v);

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:

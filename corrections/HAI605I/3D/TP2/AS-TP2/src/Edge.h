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

#include <map>

// -------------------------------------------------
// Intermediate Edge structure for hashed adjacency
// -------------------------------------------------

struct Edge {
public:
  inline Edge (unsigned int v0, unsigned int v1) { 
    if (v0 < v1) {v[0] = v0; v[1] = v1; } else {v[0] = v1; v[1] = v0; }  
  }
  inline Edge (const Edge & e) { v[0] = e.v[0]; v[1] = e.v[1]; }
  inline virtual ~Edge () {}
  inline Edge & operator= (const Edge & e) { v[0] = e.v[0]; v[1] = e.v[1]; return (*this); }
  inline bool operator== (const Edge & e) { return (v[0] == e.v[0] && v[1] == e.v[1]); }
  inline bool operator< (const Edge & e) { return (v[0] < e.v[0] || (v[0] == e.v[0] && v[1] < e.v[1])); }
  inline bool contains (unsigned int i) const { return (v[0] == i || v[1] == i); }
  unsigned int v[2];
};

struct compareEdge {
  inline bool operator()(const Edge e1, const Edge e2) const { 
    if (e1.v[0] < e2.v[0])
      return true;
    if (e1.v[0] > e2.v[0])
      return false;
    if (e1.v[1] > e2.v[1])
      return true;
    return false;
  }
};

typedef std::map<Edge, unsigned int, compareEdge> EdgeMapIndex;

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:

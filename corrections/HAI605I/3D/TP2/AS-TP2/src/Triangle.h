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

class Triangle {
public:
  inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2) { init (v0, v1, v2); }
  inline Triangle (unsigned int * vp) { init (vp[0], vp[1], vp[2]); }
  inline Triangle (const Triangle & it) { init (it.v[0], it.v[1], it.v[2]); }
  inline virtual ~Triangle () {}
  inline Triangle & operator= (const Triangle & it) { init (it.v[0], it.v[1], it.v[2]); return (*this); }
  inline bool operator== (const Triangle & t) const { return (v[0] == t.v[0] && v[1] == t.v[1] && v[2] == t.v[2]); }
  inline unsigned int getVertex (unsigned int i) const { return v[i]; }
  inline void setVertex (unsigned int i, unsigned int vertex) { v[i] = vertex; }
  inline bool contains (unsigned int vertex) const { return (v[0] == vertex || v[1] == vertex || v[2] == vertex); }
  
protected:
  inline void init (unsigned int v0, unsigned int v1, unsigned int v2) {
	v[0] = v0; v[1] = v1; v[2] = v2;  
  }
  
private:
  unsigned int v[3];
};

extern std::ostream & operator<< (std::ostream & output, const Triangle & t);

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:

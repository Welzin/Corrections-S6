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

#include <cmath>
#include <iostream>

template<typename T> class Vec3D;

template <class T> bool operator!= (const Vec3D<T> & p1, const Vec3D<T> & p2) {
    return (p1[0] != p2[0] || p1[1] != p2[1] || p1[2] != p2[2]);
}

template <class T> const Vec3D<T> operator* (const Vec3D<T> & p, float factor) {
    return Vec3D<T> (p[0] * factor, p[1] * factor, p[2] * factor);
}

template <class T> const Vec3D<T> operator* (float factor, const Vec3D<T> & p) {
    return Vec3D<T> (p[0] * factor, p[1] * factor, p[2] * factor);
}

template <class T> const Vec3D<T> operator* (const Vec3D<T> & p1, const Vec3D<T> & p2) {
    return Vec3D<T> (p1[0] * p2[0], p1[1] * p2[1], p1[2] * p2[2]);
}

template <class T> const Vec3D<T> operator+ (const Vec3D<T> & p1, const Vec3D<T> & p2) {
    return Vec3D<T> (p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]);
}

template <class T> const Vec3D<T> operator- (const Vec3D<T> & p1, const Vec3D<T> & p2) {
    return Vec3D<T> (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
}

template <class T> const Vec3D<T> operator- (const Vec3D<T> & p) {
    return Vec3D<T> (-p[0], -p[1], -p[2]);
}

template <class T> const Vec3D<T> operator/ (const Vec3D<T> & p, float divisor) {
    return Vec3D<T> (p[0]/divisor, p[1]/divisor, p[2]/divisor);
}

template <class T> bool operator== (const Vec3D<T> & p1, const Vec3D<T> & p2) {
    return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]);
}

template <class T> bool operator< (const Vec3D<T> & a, const Vec3D<T> & b) {
    return (a[0] < b[0] && a[1] < b[1] && a[2] < b[2]);
}

template <class T> bool operator>= (const Vec3D<T> & a, const Vec3D<T> & b) {
    return (a[0] >= b[0] || a[1] >= b[1] || a[2] >= b[2]);
}


/**
 * Vector in 3 dimensions, with basics operators overloaded.
 */
template <typename T>
class Vec3D {
public:
    inline Vec3D (void)	{
        p[0] = p[1] = p[2] = T ();
    }
    inline Vec3D (T p0, T p1, T p2) {
        p[0] = p0;
        p[1] = p1;
        p[2] = p2;
    };
    inline Vec3D (const Vec3D & v) {
        init (v[0], v[1], v[2]);
    }
    inline Vec3D (T* pp) {
        p[0] = pp[0];
        p[1] = pp[1];
        p[2] = pp[2];
    };
    // ---------
    // Operators
    // ---------
    inline T& operator[] (int Index) {
        return (p[Index]);
    };
    inline const T& operator[] (int Index) const {
        return (p[Index]);
    };
    inline Vec3D& operator= (const Vec3D & P) {
        p[0] = P[0];
        p[1] = P[1];
        p[2] = P[2];
        return (*this);
    };
    inline Vec3D& operator+= (const Vec3D & P) {
        p[0] += P[0];
        p[1] += P[1];
        p[2] += P[2];
        return (*this);
    };
    inline Vec3D& operator-= (const Vec3D & P) {
        p[0] -= P[0];
        p[1] -= P[1];
        p[2] -= P[2];
        return (*this);
    };
    inline Vec3D& operator*= (const Vec3D & P) {
        p[0] *= P[0];
        p[1] *= P[1];
        p[2] *= P[2];
        return (*this);
    };
    inline Vec3D& operator*= (T s) {
        p[0] *= s;
        p[1] *= s;
        p[2] *= s;
        return (*this);
    };
    inline Vec3D& operator/= (const Vec3D & P) {
        p[0] /= P[0];
        p[1] /= P[1];
        p[2] /= P[2];
        return (*this);
    };
    inline Vec3D& operator/= (T s) {
        p[0] /= s;
        p[1] /= s;
        p[2] /= s;
        return (*this);
    };

    //---------------------------------------------------------------

    inline Vec3D & init (T x, T y, T z) {
        p[0] = x;
        p[1] = y;
        p[2] = z;
        return (*this);
    };
    inline T getSquaredLength() const {
        return (dotProduct (*this, *this));
    };
    inline T getLength() const {
        return (T)sqrt (getSquaredLength());
    };
    /// Return length after normalization
    inline T normalize (void) {
        T length = getLength();
        if (length == 0.0f)
            return 0;
        T rezLength = 1.0f / length;
        p[0] *= rezLength;
        p[1] *= rezLength;
        p[2] *= rezLength;
        return length;
    };
    inline void fromTo (const Vec3D & P1, const Vec3D & P2) {
        p[0] = P2[0] - P1[0];
        p[1] = P2[1] - P1[1];
        p[2] = P2[2] - P1[2];
    };
    inline float transProduct (const Vec3D & v) const {
        return (p[0]*v[0] + p[1]*v[1] + p[2]*v[2]);
    }
    inline void getTwoOrthogonals (Vec3D & u, Vec3D & v) const {
        if (fabs(p[0]) < fabs(p[1])) {
            if (fabs(p[0]) < fabs(p[2]))
                u = Vec3D (0, -p[2], p[1]);
            else
                u = Vec3D (-p[1], p[0], 0);
        } else {
            if (fabs(p[1]) < fabs(p[2]))
                u = Vec3D (p[2], 0, -p[0]);
            else
                u = Vec3D(-p[1], p[0], 0);
        }
        v = crossProduct (*this, u);
    }
    inline Vec3D projectOn (const Vec3D & N, const Vec3D & P) const {
        T w = dotProduct (((*this) - P), N);
        return (*this) - (N * w);
    }
    static inline Vec3D segment (const Vec3D & a, const Vec3D & b) {
        Vec3D r;
        r[0] = b[0] - a[0];
        r[1] = b[1] - a[1];
        r[2] = b[2] - a[2];
        return r;
    };
    static inline Vec3D crossProduct(const Vec3D & a, const Vec3D & b) {
        Vec3D result;
        result[0] = a[1] * b[2] - a[2] * b[1];
        result[1] = a[2] * b[0] - a[0] * b[2];
        result[2] = a[0] * b[1] - a[1] * b[0];
        return(result);
    }
    static inline T dotProduct(const Vec3D & a, const Vec3D & b) {
        return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    }
    static inline T squaredDistance (const Vec3D &v1, const Vec3D &v2) {
        Vec3D tmp = v1 - v2;
        return (tmp.getSquaredLength());
    }
    static inline T distance (const Vec3D &v1, const Vec3D &v2) {
        Vec3D tmp = v1 - v2;
        return (tmp.getLength());
    }
    static inline Vec3D interpolate (const Vec3D & u, const Vec3D & v, T alpha) {
        return (u * (1.0f - alpha) + v * alpha);
    }

    // cartesion to polar coordinates
    // result:
    // [0] = length
    // [1] = angle with z-axis
    // [2] = angle of projection into x,y, plane with x-axis
    static inline Vec3D cartesianToPolar (const Vec3D &v) {
        Vec3D polar;
        polar[0] = v.getLength();
        if (v[2] > 0.0f)
            polar[1] = (T) atan (sqrt (v[0] * v[0] + v[1] * v[1]) / v[2]);
        else if (v[2] < 0.0f)
            polar[1] = (T) atan (sqrt (v[0] * v[0] + v[1] * v[1]) / v[2]) + M_PI;
        else
            polar[1] = M_PI * 0.5f;
        if (v[0] > 0.0f)
            polar[2] = (T) atan (v[1] / v[0]);
        else if (v[0] < 0.0f)
            polar[2] = (T) atan (v[1] / v[0]) + M_PI;
        else if (v[1] > 0)
            polar[2] = M_PI * 0.5f;
        else
            polar[2] = -M_PI * 0.5;
        return polar;
    }

    // polar to cartesion coordinates
    // input:
    // [0] = length
    // [1] = angle with z-axis
    // [2] = angle of projection into x,y, plane with x-axis
    static inline Vec3D polarToCartesian (const Vec3D & v) {
        Vec3D cart;
        cart[0] = v[0] * (T) sin (v[1]) * (T) cos (v[2]);
        cart[1] = v[0] * (T) sin (v[1]) * (T) sin (v[2]);
        cart[2] = v[0] * (T) cos (v[1]);
        return cart;
    }
    static inline Vec3D projectOntoVector (const Vec3D & v1, const Vec3D & v2) {
        return v2 * dotProduct (v1, v2);
    }
    inline Vec3D transformIn (const Vec3D & pos, const Vec3D & n, const Vec3D & u, const Vec3D & v) const {
        Vec3D q = (*this) - pos;
        return Vec3D (u[0]*q[0] + u[1]*q[1] + u[2]*q[2],
                      v[0]*q[0] + v[1]*q[1] + v[2]*q[2],
                      n[0]*q[0] + n[1]*q[1] + n[2]*q[2]);
    }

protected:
    T p[3];
};

template <class T> inline Vec3D<T> swap (Vec3D<T> & P, Vec3D<T> & Q) {
    Vec3D<T> tmp = P;
    P = Q;
    Q = tmp;
}

template <class T> std::ostream & operator<< (std::ostream & output, const Vec3D<T> & v) {
    output << v[0] << " " << v[1] << " " << v[2];
    return output;
}

template <class T> std::istream & operator>> (std::istream & input, Vec3D<T> & v) {
    input >> v[0] >> v[1] >> v[2];
    return input;
}

typedef Vec3D<float> Vec3Df;
typedef Vec3D<double> Vec3Dd;
typedef Vec3D<int> Vec3Di;

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:

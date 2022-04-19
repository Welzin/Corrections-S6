#ifndef MESH_H
#define MESH_H


#include <vector>
#include <string>
#include "Vec3.h"
#include "Skeleton.h"

#include <cmath>

#include <GL/glut.h>


// -------------------------------------------
// Basic Mesh class
// -------------------------------------------

struct MeshVertex {
    inline MeshVertex () {
        weights.clear();
    }
    inline MeshVertex (const Vec3 & _p, const Vec3 & _n) : position (_p), normal (_n) {
        weights.clear();
    }
    inline MeshVertex (const MeshVertex & vertex) : position (vertex.position), normal (vertex.normal) , weights(vertex.weights) {
    }
    inline virtual ~MeshVertex () {}
    inline MeshVertex & operator = (const MeshVertex & vertex) {
        position = vertex.position;
        normal = vertex.normal;
        weights = vertex.weights;
        return (*this);
    }
    // membres :
    Vec3 position; // une position
    Vec3 normal; // une normale
    std::vector< double > weights; // skinning weights
};

struct MeshTriangle {
    inline MeshTriangle () {
        v[0] = v[1] = v[2] = 0;
    }
    inline MeshTriangle (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
    }
    inline MeshTriangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;   v[1] = v1;   v[2] = v2;
    }
    inline virtual ~MeshTriangle () {}
    inline MeshTriangle & operator = (const MeshTriangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
        return (*this);
    }
    // membres :
    unsigned int v[3];
};




class Mesh {
public:
    std::vector<MeshVertex> vertices;
    std::vector<MeshTriangle> triangles;

    void loadOFF (const std::string & filename);

    Vec3 scalarToRGB( float scalar_value ) const;

    void recomputeNormals ();

    void computeSkinningWeights( Skeleton & skeleton, int n );

    void draw( int displayed_bone ) const ;

    void drawTransformedMesh( SkeletonTransformation & transfo ) const ;
};



#endif

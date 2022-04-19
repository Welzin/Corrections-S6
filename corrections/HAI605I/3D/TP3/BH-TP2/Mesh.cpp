#include "Mesh.h"
#include "Skeleton.h"
#include "Vec3.h"
#include <iostream>
#include <fstream>

Vec3 Mesh::scalarToRGB( float scalar_value ) const //Scalar_value ∈ [0, 1]
{
    Vec3 rgb;
    float H = scalar_value*360., S = 1., V = 0.85,
            P, Q, T,
            fract;

    (H == 360.)?(H = 0.):(H /= 60.);
    fract = H - floor(H);

    P = V*(1. - S);
    Q = V*(1. - S*fract);
    T = V*(1. - S*(1. - fract));

    if      (0. <= H && H < 1.)
        rgb = Vec3( V,  T,  P);
    else if (1. <= H && H < 2.)
        rgb = Vec3( Q,  V,  P);
    else if (2. <= H && H < 3.)
        rgb = Vec3( P,  V,  T);
    else if (3. <= H && H < 4.)
        rgb = Vec3( P,  Q,  V);
    else if (4. <= H && H < 5.)
        rgb = Vec3( T,  P,  V);
    else if (5. <= H && H < 6.)
        rgb = Vec3( V,  P,  Q);
    else
        rgb = Vec3( 0.,  0.,  0.);

    return rgb;
}

void Mesh::computeSkinningWeights( Skeleton & skeleton, int n ) {
    //---------------------------------------------------//
    //---------------------------------------------------//
    // code to change :

    // Indications:
    // you should compute weights for each vertex w.r.t. the skeleton bones
    // so each vertex will have B weights (B = number of bones)
    // Don't forget to normalize weights
    // these weights shoud be stored in vertex.weights:

    for( unsigned int i = 0 ; i < vertices.size() ; ++i ) {
        MeshVertex& vertex = vertices[i];
        vertex.weights.clear();

        double total_weights = 0;
        for (size_t j = 0; j < skeleton.bones.size(); ++j) {
            Bone& bone = skeleton.bones[j];

            Vec3& a = skeleton.articulations[bone.joints[0]].position;
            Vec3& b = skeleton.articulations[bone.joints[1]].position;
            Vec3& c = vertex.position;

            Vec3 ab = b - a;
            Vec3 ac = c - a;

            float norm = Vec3::dot(ac, ab) / Vec3::dot(ab, ab);
            Vec3 d = a + norm * ab;

            if (norm > 1) d = b;
            if (norm < 0) d = a;

            Vec3 cd = d - c;

            double w = std::pow(1 / cd.length(), n);
            vertex.weights.push_back(w);
            total_weights += w;
        }

        float sum = 0;
        // On normalise
        for (size_t j = 0; j < skeleton.bones.size(); ++j) {
            vertex.weights[j] /= total_weights;
            sum += vertex.weights[j];
        }

        if (abs(1 - sum) > 1e-6) {
            std::cout << sum << std::endl;
        }
    }

    //---------------------------------------------------//
    //---------------------------------------------------//
    //---------------------------------------------------//
}

void Mesh::draw( int displayed_bone ) const {

    glEnable(GL_LIGHTING);
    glBegin (GL_TRIANGLES);
    for (unsigned int i = 0; i < triangles.size (); i++)
        for (unsigned int j = 0; j < 3; j++) {
            const MeshVertex & v = vertices[triangles[i].v[j]];
            if( displayed_bone >= 0 && v.weights.size() > 0 ){
                // code to change :

                // Indications:
                //Call the function scalarToRGB so that you get the same coloring as slide 51
                //Update the color from the Vec3 resulting color
                
                Vec3 rgb = scalarToRGB(1. - v.weights[displayed_bone]);
                glColor3f (rgb[0], rgb[1], rgb[2]);
            }
            glNormal3f (v.normal[0], v.normal[1], v.normal[2]);
            glVertex3f (v.position[0], v.position[1], v.position[2]);
        }
    glEnd ();
}

void Mesh::drawTransformedMesh( SkeletonTransformation & transfo ) const {
    std::vector< Vec3 > new_positions( vertices.size() );

    //---------------------------------------------------//
    //---------------------------------------------------//
    // code to change :
    for( unsigned int i = 0 ; i < vertices.size() ; ++i ) {
        Vec3 p = vertices[i].position;

        // Indications:
        // you should use the skinning weights to blend the transformations of the vertex position by the bones.
        // to update the position use the weight and the bone transformation
        // for each bone p'=R*p+t
        new_positions[i] = Vec3(0, 0, 0);

        for (size_t j = 0; j < transfo.bone_transformations.size(); ++j) {
            BoneTransformation& bone_transfo = transfo.bone_transformations[j];

            new_positions[i] += vertices[i].weights[j] * (bone_transfo.world_space_rotation * p + bone_transfo.world_space_translation);
        }

    }
    //---------------------------------------------------//
    //---------------------------------------------------//
    //---------------------------------------------------//

    glEnable(GL_LIGHTING);
    glBegin (GL_TRIANGLES);
    for (unsigned int i = 0; i < triangles.size (); i++)
        for (unsigned int j = 0; j < 3; j++) {
            const MeshVertex & v = vertices[triangles[i].v[j]];
            Vec3 p = new_positions[ triangles[i].v[j] ];
            glNormal3f (v.normal[0], v.normal[1], v.normal[2]);
            glVertex3f (p[0], p[1], p[2]);
        }
    glEnd ();
}

void Mesh::loadOFF (const std::string & filename) {
    std::ifstream in (filename.c_str ());
    if (!in)
        exit (EXIT_FAILURE);
    std::string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    vertices.resize (sizeV);
    triangles.resize (sizeT);
    for (unsigned int i = 0; i < sizeV; i++)
        in >> vertices[i].position;
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
            in >> triangles[i].v[j];
    }
    in.close ();
    recomputeNormals ();
}

void Mesh::recomputeNormals () {
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal = Vec3 (0.0, 0.0, 0.0);
    for (unsigned int i = 0; i < triangles.size (); i++) {
        Vec3 e01 = vertices[triangles[i].v[1]].position -  vertices[triangles[i].v[0]].position;
        Vec3 e02 = vertices[triangles[i].v[2]].position -  vertices[triangles[i].v[0]].position;
        Vec3 n = Vec3::cross (e01, e02);
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
            vertices[triangles[i].v[j]].normal += n;
    }
    for (unsigned int i = 0; i < vertices.size (); i++)
        vertices[i].normal.normalize ();
}

// author: Robin Boanca
// file  : tp3.cpp
// date  : 19/04/2022


// ----------------
//           tp.cpp
// ----------------

// pour afficher les os
// dans la fonction key():
void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    // ...

        case 'x':
            displayed_bone++;
            break;
        case 'X':
            displayed_bone = max(-1, displayed_bone-1);
            break;

    // ...

    }
}

// ----------------
//         Mesh.cpp
// ----------------

// ----------------
//       Question 1
// ----------------

// dans la fonction computeSkinningWeights (question 1)


// Correction du prof
void Mesh::computeSkinningWeights( Skeleton & skeleton ) {
    for (unsigned int i =0; i < vertices.size(); ++i) {
        MeshVertex & vertex = vertices[i];
        Vec3 p = vertex.position;
        vertex.weights.resize(skeleton.bones.size());
        double sum_w_ij = 0.0;
        for (unsigned int j = 0; j < skeleton.bones.size(); ++j) {
            Bone & b = skeleton.bones[j];
            Vec3 p0 = skeleton.articulations[b.joints[0]].position;
            Vec3 p1 = skeleton.articulations[b.joints[1]].position;

            Vec3 proj = p0 + std::min<double>(1.0, 
                        std::max<double>(0.0,
                        Vec3::dot(p - p0, p1 - p0) / 
                        Vec3::dot(p1 - p0, p1 - p0))) 
                        * (p1 - p0);
            double alpha_ij = 1.0 / (proj-p).squareLength();
            alpha_ij = pow(alpha_ij, 4);
            vertex.weights[j] = alpha_ij;
            sum_w_ij += alpha_ij;
        }
        for (unsigned int j=0; j < skeleton.bones.size(); ++j) {
                vertex.weights[j] /= sum_w_ij;
        }
    }
}

// la même, mais version "je fais tous les calculs moi même"
void Mesh::computeSkinningWeights( Skeleton & skeleton, int factor_n ) {
    for( unsigned int i = 0 ; i < vertices.size() ; ++i ) {
        MeshVertex & vertex = vertices[ i ];
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            Bone & bone = skeleton.bones[j];
            Vec3 A, B, C;
            A = skeleton.articulations[bone.joints[0]].position;
            B = skeleton.articulations[bone.joints[1]].position;
            C = vertex.position;

            Vec3 AB, AC, ACP, u;
            AB = B - A;
            AC = C - A;
            double nAB = sqrt(AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]);
            double nAC = sqrt(AC[0]*AC[0] + AC[1]*AC[1] + AC[2]*AC[2]);
            double nACP = (AB[0]*AC[0] + AB[1]*AC[1] + AB[2]*AC[2]) / nAC;
            u = AB/nAB;
            Vec3 off;
            off[0] = u[0]*nACP; off[1] = u[1]*nACP; off[2] = u[2]*nACP;
            ACP = A + off;
            double dist = sqrt( (C[0]-ACP[0])*(C[0]-ACP[0]) +
                                (C[1]-ACP[1])*(C[1]-ACP[1]) +
                                (C[2]-ACP[2])*(C[2]-ACP[2]));
            if (dist != 0.0) {
                double poids = pow(1/dist,factor_n);
                vertex.weights.push_back(poids);
            } else {
                vertex.weights.push_back(0);
            }
        }
        // normaliser les poids
        double somme = 0;
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            somme+= vertex.weights[j];
        }
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            vertex.weights[j] /= somme;
        }
    }
}


// la même, mais version "j'utilise des algos tout faits"
void Mesh::computeSkinningWeights( Skeleton & skeleton , int factor_n) {
    for( unsigned int i = 0 ; i < vertices.size() ; ++i ) {
        MeshVertex & vertex = vertices[ i ];
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            Bone & bone = skeleton.bones[j];
            Vec3 A, B, C;
            A = skeleton.articulations[bone.joints[0]].position;
            B = skeleton.articulations[bone.joints[1]].position;
            C = vertex.position;

            //Distance avec la ligne
            Vec3 res;
            double r = Vec3::dot(B-A, B-A);
            if (fabs(r) < 1e-12) res = A;
            else {
                r = Vec3::dot(C-A, B-A)/r;
                if (r < 0) res = A;
                else if (r > 1) res = B;
                else {
                    Vec3 D = B-A;
                    D *= r;
                    res = A + D;
                }                
            }
            double dist = sqrt(Vec3::dot(C-res, C-res));

            if (dist != 0.0) {
                double poids = pow(1/dist,factor_n);
                vertex.weights.push_back(poids);
            } else {
                vertex.weights.push_back(0);
            }
        }
        // normaliser les poids
        double somme = 0;
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            somme+= vertex.weights[j];
        }
        for (unsigned int j=0; j<skeleton.bones.size(); ++j) {
            vertex.weights[j] /= somme;
        }
    }
}

// -------------
//    Question 2
// -------------

// dans draw():
void Mesh::draw( int displayed_bone ) const {
    // ...
            if( displayed_bone >= 0 && v.weights.size() > 0 ){
                Vec3 rgb = scalarToRGB(1 - v.weights[displayed_bone]);
                glColor3f(rgb[0], rgb[1], rgb[2]);
            }
    // ... 
}

// -------------
//    Question 3
// -------------

// dans drawTransformedMesh
void Mesh::drawTransformedMesh( SkeletonTransformation & transfo ) const {
    // ...

    for( unsigned int i = 0 ; i < vertices.size() ; ++i ) {
        Vec3 p = vertices[i].position;
        Vec3 somme;
        for (unsigned int j=0; j<transfo.bone_transformations.size(); ++j) {
            double poids = vertices[i].weights[j];
            Mat3 rotation = transfo.bone_transformations[j].world_space_rotation;
            Vec3 translation = transfo.bone_transformations[j].world_space_translation;
            if (j == 0) somme = poids * (rotation*p + translation);
            else somme += poids * (rotation*p + translation);
        }

        new_positions[ i ] = somme;
    }

    // ...
}

// -------------
//    Question 4
// -------------

// Pour faire varier n, j'ai ajouté un paramètre factor_n à la fonction computeSkinningWeights
// Dans la puissance, on utilise ce n
// comme ici : double poids = pow(1/dist,factor_n);
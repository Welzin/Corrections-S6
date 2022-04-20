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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>

#define GLEW_STATIC 1
#include <GL/glew.h>
#include <GL/glut.h>

#include "src/Shader.h"
#include "src/Vec3D.h"
#include "src/Vertex.h"
#include "src/Triangle.h"
#include "src/Mesh.h"
#include "src/Camera.h"

using namespace std;

class PhongShader : public Shader {
public:
    PhongShader () { init ("shader.vert", "shader.frag"); }
    inline virtual ~PhongShader () {}

    void setAmbientRef (float s) {
        glUniform1fARB (ambientRefLocation, s);
    }

    void setDiffuseRef (float s) {
        glUniform1fARB (diffuseRefLocation, s);
    }

    void setSpecularRef (float s) {
        glUniform1fARB (specularRefLocation, s);
    }

    void setShininess (float s) {
        glUniform1fARB (shininessLocation, s);
    }

    void setLevels (int l) {
        glUniform1iARB (levelsLocation, l);
    }
    
private:
    void init (const std::string & vertexShaderFilename,
               const std::string & fragmentShaderFilename) {
        loadFromFile (vertexShaderFilename, fragmentShaderFilename);
        bind ();
        ambientRefLocation = getUniLoc ("ambientRef");
        diffuseRefLocation = getUniLoc ("diffuseRef");
        specularRefLocation = getUniLoc ("specularRef");
        shininessLocation = getUniLoc ("shininess");
        levelsLocation = getUniLoc ("levels");
    }
    GLint ambientRefLocation;
    GLint diffuseRefLocation;
    GLint specularRefLocation;
    GLint shininessLocation;
    GLint levelsLocation;
};

static GLint window;
static unsigned int SCREENWIDTH = 1024;
static unsigned int SCREENHEIGHT = 768;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static unsigned int FPS = 0;
static bool fullScreen = false;

static PhongShader * phongShader;

static Mesh current_mesh;
static Mesh mesh_pose_0;
static Mesh mesh_pose_1;
static Mesh mesh_pose_2;
static GLuint glID;

static int levels = 4;
static float ambientRef = 0.2f;
static float diffuseRef = 0.8f;
static float specularRef = 0.5f;
static float shininess = 16.0f;
static float factor = 0.5f;
static float interpolant1 = 0.f;
static float interpolant2 = 0.f;
static bool animate = false;

typedef enum {Wire, Phong, Solid} RenderingMode;
static RenderingMode mode = Phong;
void initGLList ();
void openOFF (const std::string filename, Mesh &mesh, unsigned int normWeight) {
    vector<Vertex> V;
    vector<Triangle> T;

    ifstream in (filename.c_str ());
    if (!in)
        exit (EXIT_FAILURE);
    string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    for (unsigned int i = 0; i < sizeV; i++) {
        Vec3 p;
        in >> p;
        V.push_back (Vertex (p));
    }
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        unsigned int v[3];
        for (unsigned int j = 0; j < 3; j++)
            in >> v[j];
        T.push_back (Triangle (v[0], v[1], v[2]));
    }
    in.close ();

    Vec3 center;
    float radius;
    Vertex::scaleToUnitBox (V, center, radius);
    mesh = Mesh (V, T);
    mesh.recomputeSmoothVertexNormals (normWeight);

}

inline void glVertexVec3Df (const Vec3 & v) {
    glVertex3f (v[0], v[1], v[2]);
}

inline void glNormalVec3Df (const Vec3 & n) {
    glNormal3f (n[0], n[1], n[2]);
}

inline void glDrawPoint (const Vec3 & pos, const Vec3 & normal) {
    glNormalVec3Df (normal);
    glVertexVec3Df (pos);
}

inline void glDrawPoint (const Vertex & v) { 
    glDrawPoint (v.getPosition (), v.getNormal ());
}

//Lissage Exercice 1
void smooth(){
    current_mesh.smooth(factor);

    //Recalcule des normales et mise à jour de l'affichage
    current_mesh.recomputeSmoothVertexNormals(0);
    initGLList ();
}


//A completer Exercice 1 (Bonus)
void taubinSmooth(float lambda, float mu){
    current_mesh.smooth(lambda);
    current_mesh.smooth(mu);

    current_mesh.recomputeSmoothVertexNormals(0);
    initGLList ();
}

//A completer Exercice 3
void updateInterpolation (){
    //Récupérer la position des sommets du maillage courant à mettre à jour
    vector<Vertex> & V = current_mesh.getVertices ();

    //(1) Faire l'interpolation linéaire entre les positions de mesh_pose_0 et mesh_pose_1
    //Affecter le résultat aux positions de current_mesh

    const vector<Vertex> & V0 = mesh_pose_0.getVertices ();
    const vector<Vertex> & V1 = mesh_pose_1.getVertices ();

    for (int i = 0; i < V.size(); ++i) {
        V[i].position = (1 - interpolant1) * V0[i].position + interpolant1 * V1[i].position;
    }

    ////////////////////////////////////////////////////
    //(2) Faire la moyenne pondérée entre les positions de mesh_pose_0, mesh_pose_1 et Mesh_pose_2
    //Affecter le résultat aux positions de current_mesh
    const vector<Vertex> & V2 = mesh_pose_2.getVertices ();

    //Les variables interpolant1 et interpolant2 sont comprises entre 0 et 1
    //Calule les poids w0 tel que la somme de w0, w1, w2 soit égale à 1
    float w0 = 1. - interpolant1 - interpolant2;
    float w1 = interpolant1;
    float w2 = interpolant2;
    // Normaliser les poids : i.e. diviser chaque poids par la somme des poids


    for (int i = 0; i < V.size(); ++i) {
        V[i].position = w0 * V0[i].position + w1 * V1[i].position + w2 * V2[i].position;
    }

    //Mettre à jour la position des sommets V
    //v_i = w0*V0 + w1*V1 +w2*V2


    //Recalcule des normales et mise à jour de l'affichage
    current_mesh.recomputeSmoothVertexNormals(0);
    initGLList ();
}

void setShaderValues () {
    phongShader->setAmbientRef(ambientRef);
    phongShader->setDiffuseRef (diffuseRef);
    phongShader->setSpecularRef (specularRef);
    phongShader->setShininess (shininess);
    phongShader->setLevels(levels);
}

void drawMesh (bool flat) {
    const vector<Vertex> & V = current_mesh.getVertices ();
    const vector<Triangle> & T = current_mesh.getTriangles ();
    glBegin (GL_TRIANGLES);
    for (unsigned int i = 0; i < T.size (); i++) {
        const Triangle & t = T[i];
        if (flat) {
            Vec3 normal = Vec3::crossProduct (V[t.getVertex (1)].getPosition ()
                    - V[t.getVertex (0)].getPosition (),
                    V[t.getVertex (2)].getPosition ()
                    - V[t.getVertex (0)].getPosition ());
            normal.normalize ();
            glNormalVec3Df (normal);
        }
        for (unsigned int j = 0; j < 3; j++)
            if (!flat) {
                glNormalVec3Df (V[t.getVertex (j)].getNormal ());
                glVertexVec3Df (V[t.getVertex (j)].getPosition ());
            } else
                glVertexVec3Df (V[t.getVertex (j)].getPosition ());
    }
    glEnd ();
}

void drawSolidModel () {
    glEnable (GL_LIGHTING);
    glEnable (GL_COLOR_MATERIAL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glPolygonOffset (1.0, 1.0);
    glEnable (GL_POLYGON_OFFSET_FILL);
    glShadeModel (GL_FLAT);
    phongShader->bind ();
    drawMesh (true);
    glPolygonMode (GL_FRONT, GL_LINE);
    glPolygonMode (GL_BACK, GL_FILL);
    glColor3f (0.0, 0.0, 0.0);
    drawMesh (true);
    glDisable (GL_POLYGON_OFFSET_FILL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    glDisable (GL_COLOR_MATERIAL);
    glDisable (GL_LIGHTING);
    glShadeModel (GL_SMOOTH);
}

void drawPhongModel () {
    glCallList (glID);
}

void initLights () {
    float cam_pos_x, cam_pos_y, cam_pos_z;
    camera.getPosition(cam_pos_x, cam_pos_y, cam_pos_z);

    GLfloat light_position_0[4] = {42, 374, 161, 0};
    GLfloat light_position_1[4] = {473, -351, -259, 0};
    //GLfloat light_position_2[4] = {-438, 167, -48, 0};
    GLfloat light_position_2[4] = {cam_pos_x, cam_pos_y, cam_pos_z, 0};

    GLfloat direction_0[3] = {-42, -374, -161};
    GLfloat direction_1[3] = {-473, 351, 259};
    GLfloat direction_2[3] = {438, -167, 48};

    GLfloat diffuse_color_0[4] = {1.0, 1.0, 1.0, 1};
    GLfloat diffuse_color_1[4] = {0.28, 0.39, 1.0, 1};
    GLfloat diffuse_color_2[4] = {1.0, 0.69, 0.23, 1};

    GLfloat specular_color_0[4] = {0.8, 0.0, 0.0, 1};
    GLfloat specular_color_1[4] = {0.0, 0.8, 0.0, 1};
    GLfloat specular_color_2[4] = {0.0, 0.0, 0.8, 1};

    GLfloat ambient[4] = {0.4f, 0.4f, 0.4f, 1.f};

    glLightfv (GL_LIGHT0, GL_POSITION, light_position_0);
    glLightfv (GL_LIGHT0, GL_SPOT_DIRECTION, direction_0);
    glLightfv (GL_LIGHT0, GL_DIFFUSE, diffuse_color_0);
    glLightfv (GL_LIGHT0, GL_SPECULAR, specular_color_0);

    glLightfv (GL_LIGHT1, GL_POSITION, light_position_1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction_1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, diffuse_color_1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, specular_color_1);

    glLightfv (GL_LIGHT2, GL_POSITION, light_position_2);
    glLightfv (GL_LIGHT2, GL_SPOT_DIRECTION, direction_2);
    glLightfv (GL_LIGHT2, GL_DIFFUSE, diffuse_color_2);
    glLightfv (GL_LIGHT2, GL_SPECULAR, specular_color_2);

    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);

    glEnable (GL_LIGHTING);
}

void setSunriseLight () {
    glDisable (GL_LIGHT0);
    glDisable (GL_LIGHT1);
    glDisable (GL_LIGHT2);
}

void setSingleSpotLight () {
    glEnable (GL_LIGHT0);
    glDisable (GL_LIGHT1);
    glDisable (GL_LIGHT2);
}

void setDefaultMaterial () {
    GLfloat material_color[4] = {1.0, 1.0, 1., 1.0f };
    GLfloat material_specular[4] = {0.5, 0.5, 0.5, 1.0 };
    GLfloat material_ambient[4] = {1.0, 0.0, 0.0, 1.0};

    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);

    glDisable (GL_COLOR_MATERIAL);
}

void initGLList () {
    glID = glGenLists (1);
    glNewList (glID, GL_COMPILE);
    drawMesh (false);
    glEndList ();
}

void init () {
    glewInit();
    if (glewGetExtension ("GL_ARB_vertex_shader")        != GL_TRUE ||
            glewGetExtension ("GL_ARB_shader_objects")       != GL_TRUE ||
            glewGetExtension ("GL_ARB_shading_language_100") != GL_TRUE) {
        cerr << "Driver does not support OpenGL Shading Language" << endl;
        exit (EXIT_FAILURE);
    }
    if (glewGetExtension ("GL_ARB_vertex_buffer_object") != GL_TRUE) {
        cerr << "Driver does not support Vertex Buffer Objects" << endl;
        exit (EXIT_FAILURE);
    }

    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    glClearColor (0.5, 0.5, 0.5, 1.0);

    initLights ();
    setSunriseLight ();
    setDefaultMaterial ();
    openOFF(std::string("./data/camel.off"), current_mesh, 0);
    //current_mesh.addNoise();
    openOFF(std::string("./data/camel.off"), mesh_pose_0, 0);
    openOFF(std::string("./data/camel_pose_1.off"), mesh_pose_1, 0);
    openOFF(std::string("./data/camel_pose_2.off"), mesh_pose_2, 0);
    initGLList ();
    
    try {
        phongShader = new PhongShader;
        phongShader->bind ();
        setShaderValues ();
    } catch (ShaderException e) {
        cerr << e.getMessage () << endl;
        exit (EXIT_FAILURE);
    }
}

void clear () {
    delete phongShader;
    glDeleteLists (glID, 1);
}

void reshape(int w, int h) {
    camera.resize (w, h);
}

void display () {
    if (animate) {
        static bool animation_revert = false;
        interpolant1 += 0.01f * (animation_revert ? -1 : 1);
        if (abs(1 - interpolant1) < 1e-6 || interpolant1 < 1e-6) {
            animation_revert = ! animation_revert;
        }
        updateInterpolation();
    }

    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    if (mode == Solid)
        drawSolidModel ();
    else if (mode == Phong || mode == Wire )
        drawPhongModel ();
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char FPSstr [128];
        unsigned int numOfTriangles = current_mesh.getTriangles ().size ();
        if (mode == Solid)
            sprintf (FPSstr, "HAI60I - Examen TP: %d tri. - solid shading - %d FPS.",
                     numOfTriangles, FPS);
        else if (mode == Phong)
            sprintf (FPSstr, "HAI60I - Examen TP: %d tri. - Phong shading - %d FPS.",
                     numOfTriangles, FPS);
        glutSetWindowTitle (FPSstr);
        lastTime = currentTime;

    }

    glutPostRedisplay ();
}

void printUsage () {
    cerr << endl
         << "--------------------------------------" << endl
         << "HAI60I - Examen TP" << endl
         << "--------------------------------------" << endl
         << "USAGE: ./Main <file>.off" << endl
         << "--------------------------------------" << endl
         << "Keyboard commands" << endl
         << "--------------------------------------" << endl
         << " ?: Print help" << endl
         << " w: Toggle wireframe Mode" << endl
         << " f: Toggle full screen mode" << endl
         << " A/a: Increase/Decrease ambient reflection" << endl
         << " D/d: Increase/Decrease diffuse reflection" << endl
         << " S/s: Increase/Decrease specular reflection" << endl
         << " +/-: interpolation step" << endl
         << " <drag>+<left button>: rotate model" << endl
         << " <drag>+<right button>: move model" << endl
         << " <drag>+<middle button>: zoom" << endl
         << " q, <esc>: Quit" << endl << endl
         << "--------------------------------------" << endl;
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'P':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;
    case 'q':
    case 27:
        clear ();
        exit (0);
        break;
    case 'w':
        if( mode == Wire ){
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
            phongShader->bind ();
            mode = Phong;
        } else {
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
            phongShader->bind ();
            mode = Wire;
        }
        break;
    case 'D':
        diffuseRef += 0.1f;
        phongShader->setDiffuseRef (diffuseRef);
        break;
    case 'd':
        if (diffuseRef > 0.1f) {
            diffuseRef -= 0.1f;
            phongShader->setDiffuseRef (diffuseRef);
        }
        break;
    case 'S':
        specularRef += 0.1f;
        phongShader->setSpecularRef (specularRef);
        break;
    case 's':
        if (specularRef > 0.1f) {
            specularRef -= 0.1f;
            phongShader->setSpecularRef (specularRef);
        }
        break;
    case 'A':
        ambientRef += 0.1f;
        phongShader->setAmbientRef (ambientRef);
        break;
    case 'a':
        if (ambientRef > 0.1f) {
            ambientRef -= 0.1f;
            phongShader->setAmbientRef (ambientRef);
        }
        break;

    case 'N':
        levels ++;
        phongShader->setLevels(levels);
        break;
    case 'n':
        if (levels > 1) {
            levels --;
            phongShader->setLevels(levels);
        }
        break;
    case 'F':
        factor = std::min( factor + 0.01f, 1.0f);
        break;
    case 'f':
        factor = std::max( factor - 0.01f, 0.0f);
        break;
    case '+':
        interpolant1 = std::min( interpolant1 + 0.01f, 1.0f);
        if( interpolant1 + interpolant2 > 1.0f ) interpolant2 -= 0.01f;
        updateInterpolation();
        break;
    case '-':
        interpolant1 = std::max( interpolant1 - 0.01f, 0.0f);
        updateInterpolation();
        break;
    case 'M':
        interpolant2 = std::min( interpolant2 + 0.01f, 1.0f);
        if( interpolant1 + interpolant2 > 1.0f ) interpolant1 -= 0.01f;
        updateInterpolation();
        break;
    case 'm':
        interpolant2 = std::max( interpolant2 - 0.01f, 0.0f);
        updateInterpolation();
        break;
    case 'l':
        smooth();
        break;
    case 't':
        taubinSmooth(factor, -factor-0.01);
        break;
    case 'g':
        animate = !animate;
        break;
    case '?':
    default:
        printUsage ();
        break;
    }
    setShaderValues ();
    idle ();
}

void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }
    idle ();
}

void motion (int x, int y) {
    if (mouseRotatePressed == true)
        camera.rotate (x, y);
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH),
                     (lastY-y)/static_cast<float>(SCREENHEIGHT),
                     0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}

void usage () {
    printUsage ();
    exit (EXIT_FAILURE);
}



int main (int argc, char ** argv) {
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ( "HAI60I - Examen TP");


    init ();

    glCullFace (GL_BACK);
    glEnable (GL_CULL_FACE);
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);

    key ('?', 0, 0);

    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);

    phongShader->bind ();
    glutMainLoop ();
    return EXIT_SUCCESS;
}


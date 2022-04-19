

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <GL/glut.h>
#include "src/Vec3.h"
#include "src/Camera.h"
#include "src/Mesh.h"
#include "src/Skeleton.h"



using namespace std;


// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 640;
static unsigned int SCREENHEIGHT = 480;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;

Mesh mesh;
Skeleton skeleton;
SkeletonTransformation skeleton_transformation;

int displayed_bone = -1;
int n = 4;

int displayMode = 0;
// 0 : show rest pose with weights
// 1 : show procedural anim

void printUsage () {
    cerr << endl
         << "Usage : ./main" << endl
         << "Keyboard commands" << endl
         << "------------------" << endl
         << " ?: Print help" << endl
         << " w: Toggle Wireframe Mode" << endl
         << " f: Toggle full screen mode" << endl
         << " <drag>+<left button>: rotate model" << endl
         << " <drag>+<right button>: move model" << endl
         << " <drag>+<middle button>: zoom" << endl << endl;
}

void usage () {
    printUsage ();
    exit (EXIT_FAILURE);
}



// ------------------------------------

void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glEnable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);

}

void draw () {

    glColor3f(0.8, 0.8, 0.8);

    if( displayMode == 0 ) {
        mesh.draw( displayed_bone );
        skeleton.draw( displayed_bone );
    }

    if( displayMode == 1 ) {
        mesh.drawTransformedMesh( skeleton_transformation );
        skeleton.drawTransformedSkeleton( displayed_bone , skeleton_transformation );
    }
}


void display () {
    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();
    glFlush ();
    glutSwapBuffers ();
}

void updateProceduralAnim(){
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    skeleton.computeProceduralAnimation( 0.001 * currentTime , skeleton_transformation );
}

void idle () {
    updateProceduralAnim();
    glutPostRedisplay ();
}

void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;

    case 'd':
        displayMode = (displayMode + 1) % 2;
        break;

    case 'w':
        GLint polygonMode[2];
        glGetIntegerv(GL_POLYGON_MODE, polygonMode);
        if(polygonMode[0] != GL_FILL)
            glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        else
            glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        break;

    case '+':
        n++;
        n %= 10;
        std::cout << n << std::endl;
        mesh.computeSkinningWeights( skeleton, n );
        break;

    case '-':
        n--;
        n += 10;
        n %= 10;
        std::cout << n << std::endl;
        mesh.computeSkinningWeights( skeleton, n );
        break;

    default:
        printUsage ();
        break;
    }
    idle ();
}
void specialKey(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP :
        if( skeleton.bones.size() > 0 ) {
            ++displayed_bone;
            displayed_bone = std::min<int>( displayed_bone , skeleton.bones.size()-1 );
        }
        break;
    case GLUT_KEY_DOWN :
        --displayed_bone;
        displayed_bone = std::max( displayed_bone , -1 );
        break;
    }
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
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}




int main (int argc, char ** argv) {
    if (argc > 2) {
        printUsage ();
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("anim");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutSpecialFunc(specialKey);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);


    mesh.loadOFF("models/Draco.off");
    skeleton.load("models/Draco.skel");
    mesh.computeSkinningWeights( skeleton, n );
    updateProceduralAnim();

    glutMainLoop ();
    return EXIT_SUCCESS;
}


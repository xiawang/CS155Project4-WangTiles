//
// Author: Xiaotian Wang
//

#include <stdlib.h>
#include <iostream>
#include <glut/glut.h>
#include <math.h>


using namespace std;

#define PI 3.14159

// function prototypes
void display(void);
void reshape(int width, int height);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void init(void);
void shadowMatrixPointLight(float matrix[16], float groundPlane[4], float lightPosition[4]);
void drawCube(void);
void drawFloor(void);

// viewpoint
double theta=1.0, phi=1.0, distanceToViewpoint=15;
// window
int width=400;
int height=400;

// global yuck
float pointLightPosition[]={-3.0f,2.0f,0.0f,1.0f};

int
main(int argc, char **argv)
{
    
    // set up window
    glutInitWindowSize(500, 500);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL);
    glutCreateWindow("Shadows");
    
    // register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMotionFunc(motion);
    
    // initalize opengl parameters
    init();
    
    // loop until something happens
    glutMainLoop();
    return 0;
}

void init()
{
    // initialize viewing system
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, 1.0, 1.0, 1000.0);
    glMatrixMode(GL_MODELVIEW);
    
    // initialize clear colors
    glClearColor(0,0,0,0);
    
    // initialize shade model to smooth
    glShadeModel(GL_SMOOTH);
    
    // initialize background color to white
    glClearColor(1.0,1.0,1.0,0.0);
    
    // enable depth buffering
    glEnable(GL_DEPTH_TEST);
    
    // create some ambient light
    float ambient[]={0.5f,.5f,.5f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
    
    // create a point light
    float pointLightColor[]={1.0f,1.0f,1.0f,1.0f};
    
    glLightfv(GL_LIGHT0, GL_DIFFUSE, pointLightColor);
    
    
    // enable lighting and the point light
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    
    // set polygon mode
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    
}

void motion(int x, int y)
{
    
    static int currX=-1;
    static int currY=-1;
    
    if (currX>0  && abs(x-currX) < width/6 && abs(y-currY) < height/6) {
        phi += (double) (x-currX)/ (double) width * 3.14159;
        theta += (double) (y-currY)/ (double) height * 3.14159;
        // limit theta to -4pi/9 and 4pi/9
        if (theta < -4*3.14159/9.0)
            theta = -4*3.14159/9.0;
        if (theta > 4*3.14159/9.0)
            theta = 4*3.14159/9.0;
    }
    currX = x;
    currY = y;
    glutPostRedisplay();
}

void reshape(int width, int height)
{
    
    glViewport(0,0,width,height);
}

void display()
{
    
    // check for openGL errors
    GLenum errCode;
    if ((errCode = glGetError()) != GL_NO_ERROR)
        cout << gluErrorString(errCode) << endl;
    
    // set up modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // compute eye
    float eyex = (float) -distanceToViewpoint * cos(theta) * sin(phi);
    float eyey = (float) distanceToViewpoint * sin(theta);
    float eyez = (float) distanceToViewpoint * cos(theta) * cos(phi);
    
    gluLookAt(eyex,eyey,eyez, 0,0,0, 0, 1, 0);
    glPushMatrix();
    
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    
    // set point light position
    glLightfv(GL_LIGHT0, GL_POSITION, pointLightPosition);
    
    
    // draw floor
    float floorDiffuse[]={1.0,0.0,0.0};
    float floorAmbient[]={1.0,0.0,0.0};
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, floorDiffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, floorDiffuse);
    drawFloor();
    
    glutSwapBuffers();
    
}


void drawCube(void)
{
    glPushMatrix();
    glTranslatef(1,.5,0);
    glutSolidCube( 1.0);
    glPopMatrix();
}

void drawFloor(void)
{
    glNormal3f(0,1,0);
    glBegin(GL_QUADS);
    glVertex3f(-5,0,-5);
    glVertex3f(-5,0,5);
    glVertex3f(5,0,5);
    glVertex3f(5,0,-5);
    glEnd();
}


/* Create a matrix that will project the desired shadow. */
void
shadowMatrixPointLight(GLfloat shadowMat[16],
                       GLfloat groundplane[4],
                       GLfloat lightpos[4])
{
    GLfloat dot;
    
    
    /* Find dot product between light position vector and ground plane normal. */
    dot = groundplane[0] * lightpos[0] +
    groundplane[1] * lightpos[1] +
    groundplane[2] * lightpos[2] +
    groundplane[3] * lightpos[3];
    
    shadowMat[0] = dot - lightpos[0] * groundplane[0];
    shadowMat[4] = 0.f - lightpos[0] * groundplane[1];
    shadowMat[8] = 0.f - lightpos[0] * groundplane[2];
    shadowMat[12] = 0.f - lightpos[0] * groundplane[3];
    
    shadowMat[1] = 0.f - lightpos[1] * groundplane[0];
    shadowMat[5] = dot - lightpos[1] * groundplane[1];
    shadowMat[9] = 0.f - lightpos[1] * groundplane[2];
    shadowMat[13] = 0.f - lightpos[1] * groundplane[3];
    
    shadowMat[2] = 0.f - lightpos[2] * groundplane[0];
    shadowMat[6] = 0.f - lightpos[2] * groundplane[1];
    shadowMat[10] = dot - lightpos[2] * groundplane[2];
    shadowMat[14] = 0.f - lightpos[2] * groundplane[3];
    
    shadowMat[3] = 0.f - lightpos[3] * groundplane[0];
    shadowMat[7] = 0.f - lightpos[3] * groundplane[1];
    shadowMat[11] = 0.f - lightpos[3] * groundplane[2];
    shadowMat[15] = dot - lightpos[3] * groundplane[3];
}
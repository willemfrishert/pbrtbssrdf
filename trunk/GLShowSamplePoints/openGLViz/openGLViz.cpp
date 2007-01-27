/* -*- Mode: c; c-indentation-style: stroustrup; c-basic-offset: 4 -*- */

#include <Windows.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

struct Point3d 
{
	Point3d(float x, float y, float z):
	x(x), y(y), z(z) {}
	float x;
	float y;
	float z;
};


GLdouble angle;
GLdouble wobble;
GLfloat xRot = 0;
GLfloat yRot = 0.0f;
GLfloat xRotOld = 0.0f;
GLfloat yRotOld = 0.0f;
int mouseState = 0;
int xCenter = 0;
int yCenter = 0;

bool initpoints = true;
bool repelpoints = true;
bool endpoints = true;

float numberOfIterations;
float numberOfSamples;
float numberOfForcePoints;

vector<Point3d> stratifiedPoints;
vector<Point3d> initialPoints;
vector<Point3d> endPoints;
vector<Point3d> turkPoints;

#define M_ROTATE_XY     1


void ParseFile( string fileName, vector<Point3d>& iPoints, bool storeMetaData )
{
	fstream input( fileName.c_str(), fstream::in);
	float nrOfIterations;
	float nrOfSamples;

	input >> nrOfIterations;
	input >> nrOfSamples;
	
	if (storeMetaData)
	{
		numberOfIterations = nrOfIterations;
		numberOfSamples = nrOfSamples;
		numberOfForcePoints = numberOfIterations*numberOfSamples;
	}

	while(!input.eof())
	{
		Point3d p(0,0,0);
		input >> p.x >> p.y >> p.z;
		iPoints.push_back(p);
	}

	iPoints.pop_back();

	input.close();
}

void displayVec(vector<Point3d>& points)
{
	glPushMatrix();
	glTranslatef(0.0, -1.5, 0.0);
	if (initpoints)
	{
		glColor3f(1, 0, 0);

		glEnableClientState( GL_VERTEX_ARRAY );
		glVertexPointer( 3, GL_FLOAT, 0, &initialPoints[0].x );
		glDrawArrays(GL_POINTS, 0, numberOfSamples);
	
		glDisableClientState( GL_VERTEX_ARRAY );
	}

	if (repelpoints)
	{
		glColor3f(0, 0, 1);
		glPushMatrix();
		glEnableClientState( GL_VERTEX_ARRAY );
		glVertexPointer( 3, GL_FLOAT, 0, &turkPoints[0].x );
		glDrawArrays(GL_POINTS, 0, numberOfForcePoints);
	
		glDisableClientState( GL_VERTEX_ARRAY );
	}

	if (endpoints)
	{
		glColor3f(0, 1, 0);
		glPushMatrix();
		glEnableClientState( GL_VERTEX_ARRAY );
		glVertexPointer( 3, GL_FLOAT, 0, &endPoints[0].x );
		glDrawArrays(GL_POINTS, 0, numberOfSamples);
	
		glDisableClientState( GL_VERTEX_ARRAY );
	}
	glPopMatrix();
}

void init() {
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glPointSize(2.0);
	glColor3f(1, 0, 0);

	glClearColor(0, 0, 0, 0);
	angle = 0;
}

void anim() {
	glutPostRedisplay();
}

void display(void) {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glPushMatrix();
	glTranslatef(0.0, 0, -5.0);

		glPushMatrix();
			/* "World" rotation, controlled by mouse */
			glRotatef(xRot, 1, 0, 0);
			glRotatef(yRot, 0, 1, 0);

			displayVec( initialPoints );

		glPopMatrix();

	glPopMatrix();

	glColor3f(0, 0, 0);
	glBegin(GL_QUADS);
	glVertex3f(-3,-3,0);
	glVertex3f(3,-3,0);
	glVertex3f(3,3,0);
	glVertex3f(-3,3,0);
	glEnd();
	glPopMatrix();

	glutSwapBuffers();
}

void reshape(int w, int h) {
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, (GLdouble)w/(GLdouble)h, 0.1, 100);
	glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {

	switch (key) {
	 case 3:             /* Ctrl-C */
	 case 27:            /* ESC */
		 exit(0);
	 case 'r':
		 repelpoints = !repelpoints;
		 break;
	 case 'i':
		 initpoints = !initpoints;
		 break;
	 case 'e':
		 endpoints = !endpoints;
	}
}

void mouse(int button, int state, int x, int y) {
	xCenter = x;
	yCenter = y;

	if (state == GLUT_DOWN) {
		if (button == GLUT_LEFT_BUTTON) {
			mouseState = M_ROTATE_XY;
			xRotOld = xRot;
			yRotOld = yRot;
		}
	} else {
		mouseState = 0;
	}
}

void motion(int x, int y) {
	if (mouseState == M_ROTATE_XY) {
		xRot = xRotOld + (float)(y - yCenter) / 4.0;
		yRot = yRotOld + (float)(x - xCenter) / 4.0;
	}
}

int main(int argc, char **argv) 
{

	ParseFile("../../scenes/stratified2D.txt", stratifiedPoints, false);
	ParseFile("../../scenes/mapped3D.txt", initialPoints, false );
	ParseFile("../../scenes/turk3D.txt", turkPoints, true );
	ParseFile("../../scenes/final3D.txt", endPoints, false );
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(800, 600);
	glutCreateWindow("Test Stratified Points");
	init();
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutIdleFunc(anim);
	glutReshapeFunc(reshape);
	glutMainLoop();
}
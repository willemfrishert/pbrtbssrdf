/* -*- Mode: c; c-indentation-style: stroustrup; c-basic-offset: 4 -*- */

#include <Windows.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;


GLdouble angle;
GLdouble wobble;
GLfloat xRot = 0;
GLfloat yRot = 0.0f;
GLfloat xRotOld = 0.0f;
GLfloat yRotOld = 0.0f;
int mouseState = 0;
int xCenter = 0;
int yCenter = 0;

#define M_ROTATE_XY     1

struct Point3d 
{
	Point3d(float x, float y, float z):
		x(x), y(y), z(z) {}
	float x;
	float y;
	float z;
};

void ParseFile( string fileName, vector<Point3d>& iPoints )
{
	fstream input( fileName.c_str(), fstream::in);

	while(!input.eof())
	{
		Point3d p(0,0,0);
		input >> p.x >> p.y >> p.z;
		iPoints.push_back(p);
	}

	iPoints.pop_back();

	input.close();
}

vector<Point3d> stratifiedPoints;
vector<Point3d> initialPoints;
vector<Point3d> turkPoints;

void displayVec(vector<Point3d>& points)
{
	vector<Point3d>::iterator it = initialPoints.begin();
	
	glColor3f(1, 0, 0);
	glPushMatrix();
	glTranslatef(-0.75, -0.75, 0.75);
	glBegin(GL_POINTS);
		for (; it != initialPoints.end(); it++)
		{
			Point3d p = *it;
			glVertex3f(p.x, p.y, p.z);
		}
	glEnd();
	glPopMatrix();


	it = turkPoints.begin();

	glColor3f(0, 0, 1);
	glPushMatrix();
	glTranslatef(-0.75, -0.75, 0.75);
	glBegin(GL_POINTS);
	for (; it != turkPoints.end(); it++)
	{
		Point3d p = *it;
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();
	glPopMatrix();



	glColor3f(1,1,1);
//	glutWireCube( 1.5 );

	glColor3f(0,0,0);
	glutSolidCube( 1.499 );

}

void init() {
	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glPointSize(5.0);
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
	glTranslatef(0, 0, -4);

	glPushMatrix();
	/* "World" rotation, controlled by mouse */
	glRotatef(xRot, 1, 0, 0);
	glRotatef(yRot, 0, 1, 0);

	displayVec( initialPoints );

	glPopMatrix();

	glPushMatrix();
	glColor3f(1, 0, 0);
	glTranslatef(0.75, -0.5, 1.0);
	glBegin(GL_POINTS);

	vector<Point3d>::iterator it = stratifiedPoints.begin();
	for (; it != stratifiedPoints.end(); it++)
	{
		Point3d p = *it;
		glVertex3f(p.x, p.y, p.z);
	}
	glEnd();

	glColor3f(1, 1, 1);
	glBegin(GL_LINE_LOOP);
	glVertex3f(0,0,0);
	glVertex3f(1.0,0,0);
	glVertex3f(1.0,1.0,0);
	glVertex3f(0,1.0,0);
	glEnd();
	glPopMatrix();

//	glutSolidCube(1.0);

	glutSwapBuffers();
}

void reshape(int w, int h) {
	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, (GLdouble)w/(GLdouble)h, 0.1, 10);
	glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {

//	switch (key) {
////	 case 3:             /* Ctrl-C */
////	 case 27:            /* ESC */
////		 exit(0);
//	}
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

int main(int argc, char **argv) {

	//points.push_back( Point3d(0, 0, 0) );
	//points.push_back( Point3d(0, 1, 0) );
	//points.push_back( Point3d(1, 0, 0) );

	ParseFile("../../scenes/stratified2D.txt", stratifiedPoints);
	ParseFile("../../scenes/mapped3D.txt", initialPoints );
	ParseFile("../../scenes/turk3D.txt", turkPoints );
	

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
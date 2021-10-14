#include <iostream>
#include <gl/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <Windows.h>
#include "NurbsSurface_qj.h"
#include "NurbsSurface_triangle.h"


#pragma comment(linker, "/subsystem:console") 
#pragma comment( lib, "glew32.lib")
#pragma comment( lib, "glu32.lib")
#pragma comment( lib, "glut32.lib")
#pragma comment( lib, "opengl32.lib")
std::vector< Member> mVec;//ȫ�ֱ������洢member
std::vector< Node> Node_uv;


const int vertexMax =100;//���漫ֵ������ֵ
const int vertexMin = -80;//���漫ֵ�����Сֵ
bool surface_flag=FALSE;
int meshtype=1;
GLUnurbsObj *theNurb1;

GLfloat ctrlpoints[6][6][3] =
{
	{ { 0, 0, 0 },
	{ 20, 0, 0 },
	{ 40, 0, 0 },
	{ 60, 0, 0 },
	{ 80, 0, 0 },
	{ 100,0, 0 } },

	{ { 0, 20, 0 },
	{ 20, 20, 10 },
	{ 40, 20, 15 },
	{ 60, 20, 20 },
	{ 80, 20, 10 },
	{ 100, 20, 0 } },

	{ { 0, 40, 0 },
	{ 20, 40, 5 },
	{ 40, 40, 20 },
	{ 60, 40, vertexMax },
	{ 80, 40, 20 },
	{ 100, 40, 0 } },

	{ { 0, 60, 0 },
	{ 20, 60, 0 },
	{ 40, 60, 10 },
	{ 60, 60, 20 },
	{ 80, 60, 15 },
	{ 100, 60, 0 } },

	{ { 0, 80, 0 },
	{ 20, 80, vertexMin },
	{ 40, 80, 0 },
	{ 60, 80, 5 },
	{ 80, 80, 10 },
	{ 100, 80, 0 } },

	{ { 0, 100, 0 },
	{ 20, 100, 0 },
	{ 40, 100, 0 },
	{ 60, 100, 0 },
	{ 80, 100, 0 },
	{ 100, 100, 0 } }

};//���Ƶ�
//GLfloat ctrlpoints[6][6][3] =
//{
//	{ { 0, 0, 0 },
//	{ 20,-10, 0 },
//	{ 40, -20, 0 },
//	{ 60, -20, 0 },
//	{ 80, -10, 0 },
//	{ 100, 0, 0 } },
//
//	{ { -10, 20, 0 },
//	{ 20, 15, 10 },
//	{ 40, 10, 15 },
//	{ 60, 10, 20 },
//	{ 80, 15, 10 },
//	{ 110, 20, 0 } },
//
//	{ { -15, 40, 0 },
//	{ 20, 30, 5 },
//	{ 40, 25, 20 },
//	{ 60, 25, vertexMax },
//	{ 80, 30, 20 },
//	{ 115, 40, 0 } },
//
//	{ { -20, 60, 0 },
//	{ 20, 60, 0 },
//	{ 40, 60, 10 },
//	{ 60, 60, 20 },
//	{ 80, 60, 15 },
//	{ 120, 60, 0 } },
//
//	{ { -10, 80, 0 },
//	{ 20, 90, vertexMin },
//	{ 40, 95, 0 },
//	{ 60, 95, 5 },
//	{ 80, 90, 10 },
//	{ 110, 80, 0 } },
//
//	{ { 0, 100, 0 },
//	{ 20, 110, 0 },
//	{ 40, 120, 0 },
//	{ 60, 120, 0 },
//	{ 80, 110, 0 },
//	{ 100, 100, 0 } }
//
//};//���Ƶ�

GLfloat mat_diffuse[] = { 1.0, 0.5, 0.1, 1.0 };
GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess[] = { 100.0 };
GLfloat light_position[] = { 0.0, -10.0, 0.0, 1.0 };


//����ͼ����ɫ�����յ�OpenGL
void myInit(void)
{
	glClearColor(1.0, 1.0, 1.0, 0.0);//���ñ���ɫ
	glEnable(GL_DEPTH_TEST);

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);//����ǿ�ȣ������䣩����ʾ�������䵽�ò����ϣ�������������γɵĹ���ǿ�ȣ���ɫ��
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);//����ǿ�ȣ����淴�䣩����ʾ�������䵽�ò����ϣ��������淴����γɵĹ���ǿ�ȣ���ɫ
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);//����ָ������ȡֵ��Χ��0��128����ֵԽС����ʾ����Խ�ֲڣ����Դ����Ĺ������䵽���棬Ҳ���Բ����ϴ������
	glLightfv(GL_FRONT, GL_POSITION, light_position);//���ù�Դ����
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);//���ù���ģ�Ͳ���

	//glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LEQUAL);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glEnable(GL_BLEND);
	glFrontFace(GL_CW);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	theNurb1 = gluNewNurbsRenderer();//����NURBS����theNurb1

	gluNurbsProperty(theNurb1, GLU_SAMPLING_METHOD, GLU_PARAMETRIC_ERROR);
	gluNurbsProperty(theNurb1, GLU_PARAMETRIC_TOLERANCE, 25);
	gluNurbsProperty(theNurb1, GLU_DISPLAY_MODE, GLU_OUTLINE_POLYGON);


}

//ͼ����תƽ�����ſ���
static void myKey(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'x'://��x����ת
		glRotatef(3, 1.0, 0.0, 0.0);
		glutPostRedisplay();
		break;

	case 'y'://��y����ת
		glRotatef(3, 0.0, 1.0, 0.0);
		glutPostRedisplay();
		break;

	case 'z'://��z����ת

		glRotatef(3, 0.0, 0.0, 1.0);
		glutPostRedisplay();
		break;

	case '='://�Ŵ�
		glScalef(1.1, 1.1, 1.1);
		glutPostRedisplay();
		break;

	case '-'://��С
		glScalef(0.9, 0.9, 0.90);
		glutPostRedisplay();
		break;


	case 'a'://����ƽ��
		glTranslatef(-5, 0, 0);
		glutPostRedisplay();
		break;

	case 'd'://����ƽ��

		glTranslatef(5, 0, 0);
		glutPostRedisplay();
		break;

	case 'w':///����ƽ��

		glTranslatef(0, 5, 0);
		glutPostRedisplay();
		break;

	case 's'://����ƽ��

		glTranslatef(0, -5, 0);
		glutPostRedisplay();
		break;
	case 27:
		exit(0);
	default:
		break;
	}
}

void myDisplay()//��ͼ
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GLfloat knots[11] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0 };

	Matrix_HPoint3Dd cont(6, 6);

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cont(i, j) = PLib::HPoint3Dd(ctrlpoints[i][j][0], ctrlpoints[i][j][1], ctrlpoints[i][j][2], 1);
		}
	}
	
	int degree = 4;//����

	Vector_DOUBLE knots1;//Vector_DOUBLE��typedef��double����
	knots1.resize(11);
	knots1[0] = 0;
	knots1[1] = 0;
	knots1[2] = 0;
	knots1[3] = 0;
	knots1[4] = 0;
	knots1[5] = 0.5;
	knots1[6] = 1.0;
	knots1[7] = 1.0;
	knots1[8] = 1.0;
	knots1[9] = 1.0;
	knots1[10] = 1.0;

	if (surface_flag == TRUE) {
		glPushMatrix();
		glTranslatef(-20.0, -55, 0);
		glRotatef(45.0, 1.0, 0.0, 0.0);
		glRotatef(40.0, 1.0, 0.0, 1.0);
		glRotatef(30.0, 0.0, 1.0, 0.0);//xyz��ת
		gluBeginSurface(theNurb1);

		glColor3f(0, 0, 0);
		gluNurbsSurface(theNurb1, 11, knots, 11, knots, 6 * 3, 3, &ctrlpoints[0][0][0], 5, 5, GL_MAP2_VERTEX_3);
		gluEndSurface(theNurb1);
		glPopMatrix();
	}

	NurbsSurface_qj nurbsSurface1(degree, degree, knots1, knots1, cont, 0, 0);//***************************************************��Ҫ�޸�
	//NurbsSurface_triangle nurbsSurface1(degree, degree, knots1, knots1, cont,0,0);//***************************************************��Ҫ�޸�

	glPushMatrix();
	glTranslatef(-20.0, -55, 0);
	glRotatef(45.0, 1.0, 0.0, 0.0);
	glRotatef(40.0, 1.0, 0.0, 1.0);
	glRotatef(30.0, 0.0, 1.0, 0.0);//xyz��ת

	glPointSize(3.0);
	glColor3f(1.0, 0.0, 0.0);
	
	
	Node *tempNode1;
	Node *tempNode2;
	glLineWidth(2);
	glBegin(GL_LINES);

	for (size_t i = 0; i < mVec.size(); ++i)
	{
		tempNode1 = nurbsSurface1.GetSurfacePointatUandV(Node_uv[mVec[i].getInode()].getX(), Node_uv[mVec[i].getInode()].getY());
		tempNode2 = nurbsSurface1.GetSurfacePointatUandV(Node_uv[mVec[i].getJnode()].getX(), Node_uv[mVec[i].getJnode()].getY());		
		glVertex3f(GLfloat((*tempNode1).getX()), GLfloat((*tempNode1).getY()), GLfloat((*tempNode1).getZ()));
		glVertex3f(GLfloat((*tempNode2).getX()), GLfloat((*tempNode2).getY()), GLfloat((*tempNode2).getZ()));
		delete tempNode1;
		delete tempNode2;		
	}
	glEnd();
	glPopMatrix();
	glutSwapBuffers();

}

void myReshape(GLsizei w, GLsizei h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-60, 60, -60 * h / w, 60 * h / w, -10000, 10000);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScalef(0.8, 0.8, 0.8);
	glutSwapBuffers();

}

int main(int argc, char ** argv)
{
	Matrix_HPoint3Dd cont(6, 6);
	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 6; j++)
		{
			cont(i, j) = PLib::HPoint3Dd(ctrlpoints[i][j][0], ctrlpoints[i][j][1], ctrlpoints[i][j][2], 1);
		}
	}

	int degree = 4;

	Vector_DOUBLE knots;
	knots.resize(11);
	knots[0] = 0;
	knots[1] = 0;
	knots[2] = 0;
	knots[3] = 0;
	knots[4] = 0;
	knots[5] = 0.5;
	knots[6] = 1.0;
	knots[7] = 1.0;
	knots[8] = 1.0;
	knots[9] = 1.0;
	knots[10] = 1.0;
	const double  USplitLength = 5;
	const double  VSplitLength = 5;


	if (meshtype == 1)
	{
		NurbsSurface_qj nurbsSurface(degree, degree, knots, knots, cont, 0, 0);//��������
		mVec = nurbsSurface.surfaceAverageTotalSpit(0, 1, 0, 1, 5, 5, Node_uv);//��������*************************************�ĸ˳�
	}
	else if (meshtype == 2)
	{
		NurbsSurface_triangle nurbs_triangle(degree, degree, knots,knots, cont);
		mVec = nurbs_triangle.SurfaceAverageSpit(0, 1, 0, 1,10,10, Node_uv);//*************************************************�ĸ˳�
	}
	
	
	//NurbsSurface_triangle nurbs_triangle(degree, degree, knots,knots, cont);
	//mVec = nurbs_triangle.SurfaceAverageSpit(0, 1, 0, 1,10,10, Node_uv);//*************************************************�ĸ˳�


	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(600, 400);
	glutInitWindowPosition(200, 200);

	glutCreateWindow("NURBS surface");

	myInit();
	glutKeyboardFunc(myKey);
	glutReshapeFunc(myReshape);
	glutDisplayFunc(myDisplay);

	glutMainLoop();
	
	return(0);
}




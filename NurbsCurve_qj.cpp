#include "NurbsCurve_qj.h"
#include "Node.h"


NurbsCurve_qj::NurbsCurve_qj(void)
{
}


NurbsCurve_qj::~NurbsCurve_qj(void)
{
}


//void NurbsCurve_qj::display()
//{
//	//
//
//
//	//Init
//	glClearColor(1.0,1.0,1.0,0.0);//���ñ���ɫ 
//	theNurb = gluNewNurbsRenderer();//����NURBS����
//	gluNurbsProperty(theNurb,GLU_SAMPLING_TOLERANCE,10);
//	//
//}





NurbsCurve_qj::NurbsCurve_qj(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree)//���㣺�������Ƶ㣬�ڵ㣬�ȣ���������
{

	

	this->Degree = degree;

	this->Controlpoint.resize(control_points.size());//�������Ƶ����Ŀ
	this->Knotpoint.resize(knots.size());//�����ڵ����Ŀ

	//****************�����Ƶ㸳ֵ******************
	for(int i = 0;i < control_points.size();i++)
	{
		this->Controlpoint[i] = control_points[i];

	}
	//****************�����Ƶ㸳ֵ******************

	//*****************���ڵ㸳ֵ*******************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//*****************���ڵ㸳ֵ*******************

	this->nurbscurve.reset(control_points,knots,degree);//��������


}


NurbsCurve_qj::NurbsCurve_qj(Vector_HPoint3Dd curve_points,int degree)//���㣺���������ϵĲ����㣬�ȣ���������
{
	
	

	this->Degree = degree;//���ȸ�ֵ

	//******************�������ϵĵ㸳ֵ********************
	for(int i = 0;i < curve_points.size();i++)
	{
		this->Curvepoint[i] = curve_points[i];
	}
	//******************�������ϵĵ㸳ֵ********************

	

	this->nurbscurve.globalInterpH(Curvepoint,Degree) ;//��������

}




NurbsCurve_qj::NurbsCurve_qj(std::vector<Node*> control_points,std::vector<double> knots,int degree)//����:�����Ŀ��Ƶ���	Node��ʽ���ڵ㣬�ȣ���������
{
	this->Degree = degree;//���ȸ�ֵ



	this->Controlpoint.resize(control_points.size());//�������Ƶ�ĸ���
	this->Knotpoint.resize(knots.size());//�����ڵ�ĸ���

	Vector_HPoint3Dd thecontrol_points;//����Vector_HPoint3Dd��ʽ�Ŀ��Ƶ�
	

	

	//*********************�����Ƶ㸳ֵ*********************
	for(int i = 0;i < control_points.size();i++)
	{
		(thecontrol_points[i]).data[0] = (control_points[i]->getX());//�����Ƶ��X���긳ֵ
		(thecontrol_points[i]).data[1] = (control_points[i]->getY());//�����Ƶ��Y���긳ֵ
		(thecontrol_points[i]).data[2] = (control_points[i]->getZ());//�����Ƶ��Z���긳ֵ
		(thecontrol_points[i]).data[3] = 1;//�����Ƶ��W���긳ֵ
		
		this->Controlpoint[i] = thecontrol_points[i];

	}
	//*********************�����Ƶ㸳ֵ*********************

	//******************���ڵ㸳ֵ********************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//******************���ڵ㸳ֵ********************

	
	this->nurbscurve.reset(Controlpoint,Knotpoint,Degree);//��������

}



NurbsCurve_qj::NurbsCurve_qj(std::vector<Node*> curve_points,int degree)//���㣺���������ϵĲ�������Node��ʽ���ȣ���������
{
	Vector_HPoint3Dd thecurve_points;
	int Degree;

	this->Degree = degree;//���ȸ�ֵ

	//********************�������ϵĲ����㸳ֵ***********************
	for(int i = 0;i <= curve_points.size();i++)
	{
		
		(thecurve_points[i]).data[0] = (curve_points[0]->getX());
		(thecurve_points[i]).data[1] = (curve_points[0]->getY());
		(thecurve_points[i]).data[2] = (curve_points[0]->getZ());
		(thecurve_points[i]).data[3] = 1;

		this->Curvepoint[i] = thecurve_points[i];
	}
	//********************�������ϵĲ����㸳ֵ***********************



	this->nurbscurve.globalInterpH(Curvepoint,Degree) ;//��������


}



Node* NurbsCurve_qj::getCruvePointatU(double u)//�������ߺͱ������ӣ��������ϵĵ�
{
	


	PLib::HPoint3Dd hpoint = this->nurbscurve.pointAt(u);//�������ߺͱ������ӣ������ά�ռ��ϵĵ�
	PLib::Point3Dd point = project(hpoint);//����ά�ռ��ϵĵ�ͶӰ����ά�ռ�
	//********************�õ�Node��ʽ�ĵ�********************
	Node* thepoint = new Node();
	thepoint->setX((point.data[0]));
	thepoint->setY((point.data[1]));
	thepoint->setZ((point.data[2]));
	//********************�õ�Node��ʽ�ĵ�********************

	return thepoint;
}




int NurbsCurve_qj::getControlPoint(std::vector<Node*>& controlpoint)//������������Ƶ�
{


	
	
	int num = this->nurbscurve.ctrlPnts().size();
	//controlpoint.resize(num);
	PLib::Vector<PLib::HPoint3Dd> hpoint = nurbscurve.ctrlPnts();
	PLib::Point3Dd point;
	for(int i = 0; i < num; i++)
	{
		Node* thenode = new Node();
		point = project(hpoint[i]);
		thenode->setX(point.data[0]) ;
		thenode->setY(point.data[1]) ;
		thenode->setZ(point.data[2]) ;
		controlpoint.push_back(thenode);

	}
	
	return  num;

}






int NurbsCurve_qj::getKnotPoint(double* knotpoint)
{
	


	int num = this->nurbscurve.knot().size();
	PLib::Vector<double> a = nurbscurve.knot();

	for (int i = 0; i < num; i++)
	{
		knotpoint[i] = a[i] ;
	}

	return num;


}









void NurbsCurve_qj::drawNurbsCurve()
{
	//myInit
	glClearColor(1.0,1.0,1.0,0.0);//���ñ���ɫ 

	
	std::vector<Node*> controlpointtemp;
	int controlpointsize = this->getControlPoint(controlpointtemp);
	GLfloat** ctrlpoints;
	ctrlpoints=(GLfloat**)(new GLfloat[controlpointsize]);
	for(int i=0;i!=controlpointsize;i++)
	{
		ctrlpoints[i]=new GLfloat[3];

	}



	for (int j = 0;j<controlpointsize;++j)
	{
		
		
		ctrlpoints[j][0] = GLfloat(controlpointtemp[j]->getX());
		ctrlpoints[j][1] = GLfloat(controlpointtemp[j]->getY());
		ctrlpoints[j][2] = GLfloat(controlpointtemp[j]->getZ());

	}
	
	
	
	
	
	double *knotstemp;
    knotstemp = new double[this->Knotpoint.size()];

	int knotsSize = this->getKnotPoint(knotstemp);
	GLfloat *knots = new GLfloat[knotsSize];
	for (int i = 0; i< knotsSize; ++i)
	{
		knots[i] = GLfloat(knotstemp[i]);
	}
	
	
	Node *thenode = this->getCruvePointatU(0.5);//!!!!!!!!!!!!!!!!
	GLfloat nodef[3] = { GLfloat(thenode->getX()),  GLfloat(thenode->getY()),GLfloat(thenode->getZ())};

	GLUnurbsObj *theNurb; 
	theNurb = gluNewNurbsRenderer();//����NURBS����
	gluNurbsProperty(theNurb,GLU_SAMPLING_TOLERANCE,10);//����Nurbs��������

	//******************display**********************
	int i; 
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
	glColor3f(0.0,0.0,0.0); 
	glLineWidth(3.0); 
	/*��������*/ 
	gluBeginCurve(theNurb); 
	
	gluNurbsCurve(theNurb,8/*this->Knotpoint.size()*/,knots,3,&ctrlpoints[0][0],3/*this->Degree+1*/,GL_MAP1_VERTEX_3); //!!!!!!!!!!!!!!!!
	//gluNurbsCurve(theNurb,8,knots,4,&color[0][0],3,GL_MAP1_COLOR_4); 
	gluEndCurve(theNurb); 

	/*���Ƶ�*/ 
	glColor3f(1.0,0.0,0.0); 
	glPointSize(5.0); 
	glBegin(GL_POINTS); 
	for(i = 0;i < this->Controlpoint.size();i++)
		glVertex3f(ctrlpoints[i][0],ctrlpoints[i][1],ctrlpoints[i][2]); 

	glColor3f(0.0,1,0.0); 
	glVertex3f(nodef[0],nodef[1],nodef[2]); //!!!!!!!!!!!!!!!!

	glEnd(); 
	glutSwapBuffers();
	//******************display**********************



}



void NurbsCurve_qj::reset(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree)
{
	this->Degree = degree;

	this->Controlpoint.resize(control_points.size());//�������Ƶ����Ŀ
	this->Knotpoint.resize(knots.size());//�����ڵ����Ŀ

	//****************�����Ƶ㸳ֵ******************
	for(int i = 0;i < control_points.size();i++)
	{
		this->Controlpoint[i] = control_points[i];

	}
	//****************�����Ƶ㸳ֵ******************

	//*****************���ڵ㸳ֵ*******************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//*****************���ڵ㸳ֵ*******************

	this->nurbscurve.reset(control_points,knots,degree);//��������
}
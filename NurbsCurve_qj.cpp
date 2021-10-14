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
//	glClearColor(1.0,1.0,1.0,0.0);//设置背景色 
//	theNurb = gluNewNurbsRenderer();//创建NURBS对象
//	gluNurbsProperty(theNurb,GLU_SAMPLING_TOLERANCE,10);
//	//
//}





NurbsCurve_qj::NurbsCurve_qj(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree)//正算：给定控制点，节点，度，构造曲线
{

	

	this->Degree = degree;

	this->Controlpoint.resize(control_points.size());//给定控制点的数目
	this->Knotpoint.resize(knots.size());//给定节点的数目

	//****************给控制点赋值******************
	for(int i = 0;i < control_points.size();i++)
	{
		this->Controlpoint[i] = control_points[i];

	}
	//****************给控制点赋值******************

	//*****************给节点赋值*******************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//*****************给节点赋值*******************

	this->nurbscurve.reset(control_points,knots,degree);//构造曲线


}


NurbsCurve_qj::NurbsCurve_qj(Vector_HPoint3Dd curve_points,int degree)//反算：给定曲线上的采样点，度，构造曲线
{
	
	

	this->Degree = degree;//给度赋值

	//******************给曲线上的点赋值********************
	for(int i = 0;i < curve_points.size();i++)
	{
		this->Curvepoint[i] = curve_points[i];
	}
	//******************给曲线上的点赋值********************

	

	this->nurbscurve.globalInterpH(Curvepoint,Degree) ;//构造曲线

}




NurbsCurve_qj::NurbsCurve_qj(std::vector<Node*> control_points,std::vector<double> knots,int degree)//正算:给定的控制点是	Node形式，节点，度，构造曲线
{
	this->Degree = degree;//给度赋值



	this->Controlpoint.resize(control_points.size());//给定控制点的个数
	this->Knotpoint.resize(knots.size());//给定节点的个数

	Vector_HPoint3Dd thecontrol_points;//定义Vector_HPoint3Dd形式的控制点
	

	

	//*********************给控制点赋值*********************
	for(int i = 0;i < control_points.size();i++)
	{
		(thecontrol_points[i]).data[0] = (control_points[i]->getX());//给控制点的X坐标赋值
		(thecontrol_points[i]).data[1] = (control_points[i]->getY());//给控制点的Y坐标赋值
		(thecontrol_points[i]).data[2] = (control_points[i]->getZ());//给控制点的Z坐标赋值
		(thecontrol_points[i]).data[3] = 1;//给控制点的W坐标赋值
		
		this->Controlpoint[i] = thecontrol_points[i];

	}
	//*********************给控制点赋值*********************

	//******************给节点赋值********************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//******************给节点赋值********************

	
	this->nurbscurve.reset(Controlpoint,Knotpoint,Degree);//构造曲线

}



NurbsCurve_qj::NurbsCurve_qj(std::vector<Node*> curve_points,int degree)//反算：给定曲线上的采样点是Node形式，度，构造曲线
{
	Vector_HPoint3Dd thecurve_points;
	int Degree;

	this->Degree = degree;//给度赋值

	//********************给曲线上的采样点赋值***********************
	for(int i = 0;i <= curve_points.size();i++)
	{
		
		(thecurve_points[i]).data[0] = (curve_points[0]->getX());
		(thecurve_points[i]).data[1] = (curve_points[0]->getY());
		(thecurve_points[i]).data[2] = (curve_points[0]->getZ());
		(thecurve_points[i]).data[3] = 1;

		this->Curvepoint[i] = thecurve_points[i];
	}
	//********************给曲线上的采样点赋值***********************



	this->nurbscurve.globalInterpH(Curvepoint,Degree) ;//构造曲线


}



Node* NurbsCurve_qj::getCruvePointatU(double u)//给定曲线和比例因子，求曲线上的点
{
	


	PLib::HPoint3Dd hpoint = this->nurbscurve.pointAt(u);//根据曲线和比例因子，求出四维空间上的点
	PLib::Point3Dd point = project(hpoint);//将四维空间上的点投影到三维空间
	//********************得到Node形式的点********************
	Node* thepoint = new Node();
	thepoint->setX((point.data[0]));
	thepoint->setY((point.data[1]));
	thepoint->setZ((point.data[2]));
	//********************得到Node形式的点********************

	return thepoint;
}




int NurbsCurve_qj::getControlPoint(std::vector<Node*>& controlpoint)//给定曲线算控制点
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
	glClearColor(1.0,1.0,1.0,0.0);//设置背景色 

	
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
	theNurb = gluNewNurbsRenderer();//创建NURBS对象
	gluNurbsProperty(theNurb,GLU_SAMPLING_TOLERANCE,10);//设置Nurbs对象属性

	//******************display**********************
	int i; 
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
	glColor3f(0.0,0.0,0.0); 
	glLineWidth(3.0); 
	/*绘制曲线*/ 
	gluBeginCurve(theNurb); 
	
	gluNurbsCurve(theNurb,8/*this->Knotpoint.size()*/,knots,3,&ctrlpoints[0][0],3/*this->Degree+1*/,GL_MAP1_VERTEX_3); //!!!!!!!!!!!!!!!!
	//gluNurbsCurve(theNurb,8,knots,4,&color[0][0],3,GL_MAP1_COLOR_4); 
	gluEndCurve(theNurb); 

	/*绘制点*/ 
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

	this->Controlpoint.resize(control_points.size());//给定控制点的数目
	this->Knotpoint.resize(knots.size());//给定节点的数目

	//****************给控制点赋值******************
	for(int i = 0;i < control_points.size();i++)
	{
		this->Controlpoint[i] = control_points[i];

	}
	//****************给控制点赋值******************

	//*****************给节点赋值*******************
	for(int j = 0;j < knots.size();j++ )
	{
		this->Knotpoint[j] = knots[j];

	}
	//*****************给节点赋值*******************

	this->nurbscurve.reset(control_points,knots,degree);//构造曲线
}
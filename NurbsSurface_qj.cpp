#include "NurbsSurface_qj.h"
//#include <nurbsS.h>

bool compare1_x( PointStruct p1,  PointStruct p2)//�Ƚϴ�С,��С����
{
	return p1.point.x() < p2.point.x();
}
bool compare2_x(PointStruct p1, PointStruct p2)//�Ƚϴ�С,�Ӵ�С
{
	return p1.point.x() > p2.point.x();
}
bool compare3_y(PointStruct p1, PointStruct p2)//�Ƚϴ�С,��С����
{
	return p1.point.y() < p2.point.y();
}
bool compare4_y(PointStruct p1, PointStruct p2)//�Ƚϴ�С,�Ӵ�С
{
	return p1.point.y() > p2.point.y();
}

bool edge_compare1(PointStruct *p1, PointStruct *p2)//���ڲ�����߽��Ե��ڲ�������
{
	return p1->point.x() < p2->point.x();
}
bool edge_compare2(PointStruct *p1, PointStruct *p2)//���ڲ�����߽��Ե��ڲ�������
{
	return p1->point.x() > p2->point.x();
}
bool edge_compare3(PointStruct *p1, PointStruct *p2)//���ڲ�����߽��Ե��ڲ�������
{
	return p1->point.y() < p2->point.y();
}
bool edge_compare4(PointStruct *p1, PointStruct *p2)//���ڲ�����߽��Ե��ڲ�������
{
	return p1->point.y() > p2->point.y();
}


bool compare_CPoint2D(CPoint2D p1, CPoint2D p2)//�Ƚϴ�С,�Ӵ�С
{
	return p1.y > p2.y;
}



NurbsSurface_qj::NurbsSurface_qj(void)
{
}


NurbsSurface_qj::~NurbsSurface_qj(void)
{
}



NurbsSurface_qj::NurbsSurface_qj(Vector_HPoint3Dd ccurve_points, int cdegree, Vector_HPoint3Dd tcurve_points, int tdegree)//��������C�Ĳ�����Ͷȣ���������T�Ĳ�����Ͷȣ�������
{

	this->CDegree = cdegree;//���ȸ�ֵ
	CCurvepoint.resize(ccurve_points.size());

	//******************�������ϵĵ㸳ֵ********************
	for (int i = 0; i < ccurve_points.size(); i++)
	{
		this->CCurvepoint[i] = ccurve_points[i];
	}
	//******************�������ϵĵ㸳ֵ********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//��������C

	this->TDegree = tdegree;//���ȸ�ֵ
	TCurvepoint.resize(tcurve_points.size());

	//******************�������ϵĵ㸳ֵ********************
	for (int i = 0; i < tcurve_points.size(); i++)
	{
		this->TCurvepoint[i] = tcurve_points[i];
	}
	//******************�������ϵĵ㸳ֵ********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//��������T

	nurbssurface.sweep(tnurbscurve, cnurbscurve, 5);
}

NurbsSurface_qj::NurbsSurface_qj(std::vector<Node*> Surface_points)//��������Ĳ����㣬������
{
	Matrix_Point3Dd Pts(4, 5);
	int i, j;
	//using namespace PLib;
	PLib::Point3Dd point;
	for (i = 0; i < Pts.rows(); ++i)
	{
		for (j = 0; j < Pts.cols(); ++j)
		{
			for (int k = 0; k < 20; ++k)
			{

				Pts(i, j) = PLib::Point3Dd(Surface_points[k]->getX(), Surface_points[k]->getY(), Surface_points[k]->getZ());
			}
		}
	}
	nurbssurface.globalInterp(Pts, 3, 3);
}

NurbsSurface_qj::NurbsSurface_qj(std::vector<Node*> ccurve_points, int cdegree, std::vector<Node*> tcurve_points, int tdegree)//��������C�Ĳ�������Node��ʽ�Ͷȣ���������T�Ĳ�������Node��ʽ�Ͷȣ�������	
{
	Vector_HPoint3Dd theccurve_points;
	int CDegree;

	this->CDegree = cdegree;//���ȸ�ֵ

	//********************�������ϵĲ����㸳ֵ***********************
	for (int i = 0; i <= ccurve_points.size(); i++)
	{

		(theccurve_points[i]).data[0] = (ccurve_points[0]->getX());
		(theccurve_points[i]).data[1] = (ccurve_points[0]->getY());
		(theccurve_points[i]).data[2] = (ccurve_points[0]->getZ());
		(theccurve_points[i]).data[3] = 1;

		this->CCurvepoint[i] = theccurve_points[i];
	}
	//********************�������ϵĲ����㸳ֵ***********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//��������

	Vector_HPoint3Dd thetcurve_points;
	int TDegree;
	this->TDegree = tdegree;//���ȸ�ֵ

	//********************�������ϵĲ����㸳ֵ***********************
	for (int i = 0; i <= tcurve_points.size(); i++)
	{

		(thetcurve_points[i]).data[0] = (tcurve_points[0]->getX());
		(thetcurve_points[i]).data[1] = (tcurve_points[0]->getY());
		(thetcurve_points[i]).data[2] = (tcurve_points[0]->getZ());
		(thetcurve_points[i]).data[3] = 1;

		this->TCurvepoint[i] = thetcurve_points[i];
	}
	//********************�������ϵĲ����㸳ֵ***********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//��������
	nurbssurface.sweep(tnurbscurve, cnurbscurve, 5);
}


//���������uv�������ڵ�ʸ�������Ƶ�
NurbsSurface_qj::NurbsSurface_qj(int  uDeg, int vDeg, const Vector_DOUBLE& uKnots, const Vector_DOUBLE& vKnots, const Matrix_HPoint3Dd& controlPoints, unsigned long long int n_Id, unsigned long long int m_Id)
{
	nId = n_Id;
	mId = m_Id;
	PlNurbsSurfaced tempNurbsSurface(uDeg, vDeg, uKnots, vKnots, controlPoints);
	nurbssurface = tempNurbsSurface;
}


Node* NurbsSurface_qj::GetSurfacePointatUandV(double u, double v)//��������ͱ���u��v,�������ϵĵ�
{
	PLib::HPoint3Dd hpoint = this->nurbssurface.pointAt(u, v);
	PLib::Point3Dd point = project(hpoint);

	Node* thepoint = new Node();
	thepoint->setX((point.data[0]));
	thepoint->setY((point.data[1]));
	thepoint->setZ((point.data[2]));

	return thepoint;
}
void NurbsSurface_qj::reset(Vector_HPoint3Dd ccurve_points, int cdegree, Vector_HPoint3Dd tcurve_points, int tdegree)
{
	this->CDegree = cdegree;//���ȸ�ֵ
	CCurvepoint.resize(ccurve_points.size());

	//******************�������ϵĵ㸳ֵ********************
	for (int i = 0; i < ccurve_points.size(); i++)
	{
		this->CCurvepoint[i] = ccurve_points[i];
	}
	//******************�������ϵĵ㸳ֵ********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//��������C
	this->TDegree = tdegree;//���ȸ�ֵ
	TCurvepoint.resize(tcurve_points.size());

	//******************�������ϵĵ㸳ֵ********************
	for (int i = 0; i < tcurve_points.size(); i++)
	{
		this->TCurvepoint[i] = tcurve_points[i];
	}
	//******************�������ϵĵ㸳ֵ********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//��������T

	nurbssurface.sweep(tnurbscurve, cnurbscurve, 5);

}



void NurbsSurface_qj::reset(Vector_HPoint3Dd Surface_points)
{
	Matrix_Point3Dd Pts(4, 5);
	int i, j;
	//using namespace PLib;

	int k = 0;
	for (i = 0; i < Pts.rows() - 1; ++i)
	{

		for (j = 0; j < Pts.cols() - 1; ++j)
		{

			Pts(i, j) = PLib::Point3Dd(Surface_points[k].data[0], Surface_points[k].data[1], Surface_points[k].data[2]);
			++k;

		}
	}
	nurbssurface.globalInterp(Pts, 3, 3);
}





int NurbsSurface_qj::getControlPoint(std::vector<Node*>& controlpoint)
{
	PLib::Matrix<PLib::HPoint3Dd> controlpointmatrix;
	controlpointmatrix = nurbssurface.ctrlPnts();

	int rows = 0;
	rows = controlpointmatrix.rows();

	int cols = 0;
	cols = controlpointmatrix.cols();

	int num = rows*cols;
	
	PLib::Point3Dd point;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			Node* thenode = new Node();
			point = project(controlpointmatrix(i, j));
			thenode->setX(point.data[0]);
			thenode->setY(point.data[1]);
			thenode->setZ(point.data[2]);
			controlpoint.push_back(thenode);
		}

	}
	return  num;
}

int NurbsSurface_qj::getKnotUPoint(double* knotpoint)
{
	int num = this->nurbssurface.knotU().size();
	PLib::Vector<double> knotu = nurbssurface.knotU();

	for (int i = 0; i < num; i++)
	{
		knotpoint[i] = knotu[i];
	}
	return num;
}



int NurbsSurface_qj::getKnotVPoint(double* knotpoint)
{
	int num = this->nurbssurface.knotV().size();
	PLib::Vector<double> knotv = nurbssurface.knotU();

	for (int i = 0; i < num; i++)
	{
		knotpoint[i] = knotv[i];
	}
	return num;
}


//��ȡ������һ�����ߵļ�ֵ��,vֵ�̶�,u ֵ�ɱ䣬(u����Сֵ��u�����ֵ��vֵ��
std::vector<double> NurbsSurface_qj::getUExtremums(double min, double max, const double v)//��ȡ������һ�����ߵļ�ֵ��,vֵ�̶�,u ֵ�ɱ�
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//�趨����
	double u = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 1, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	if (mat(1, 0).z() >= 1e-6)
	{
		preFlagZ = 1;
	}
	else if (mat(1, 0).z() < -1e-6)
	{
		preFlagZ = -1;
	}
	while (u <= max)
	{
		u += deltaD;
		nurbssurface.deriveAtH(u, v, 1, mat);
		node2 = GetSurfacePointatUandV(u, v);
		if (mat(1, 0).z() >= 1e-6)
		{
			flagZ = 1;
		}
		else if (mat(1, 0).z() < -1e-6)
		{
			flagZ = -1;
		}

		if (flagZ*preFlagZ < 0)
		{
			if ((*node1 - *node2).GetLength() < 0.001)
			{
				vec.push_back(u);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				u = u - deltaD;
				deltaD = deltaD / 2;
				delete node2;
				continue;
			}

		}
		delete node1;
		node1 = node2;

	}
	return vec;
}


//��ȡ������һ�����ߵļ�ֵ��,uֵ�̶�,v ֵ�ɱ�,(v����Сֵ��v�����ֵ��uֵ��
std::vector<double> NurbsSurface_qj::getVExtremums(double min, double max, const double u)//��ȡ������һ�����ߵļ�ֵ��,uֵ�̶�,v ֵ�ɱ�
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//�趨����
	double v = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 1, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	if (mat(0, 1).z() >= 1e-6)
	{
		preFlagZ = 1;
	}
	else if (mat(0, 1).z() < -1e-6)
	{
		preFlagZ = -1;
	}
	while (v <= max)
	{
		v += deltaD;
		nurbssurface.deriveAtH(u, v, 1, mat);
		node2 = GetSurfacePointatUandV(u, v);
		if (mat(0, 1).z() >= 1e-6)
		{
			flagZ = 1;
		}
		else if (mat(0, 1).z() < -1e-6)
		{
			flagZ = -1;
		}

		if (flagZ*preFlagZ < 0)
		{
			if ((*node1 - *node2).GetLength() < 0.001)
			{
				vec.push_back(v);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				v = v - deltaD;
				deltaD = deltaD / 2;
				delete node2;
				continue;
			}

		}
		delete node1;
		node1 = node2;

	}
	return vec;
}


//��ȡu����Ĺյ�,(u����Сֵ��u�����ֵ��vֵ��
std::vector<double> NurbsSurface_qj::getUInflectionPoints(double min, double max, const double v)//��ȡu����Ĺյ�
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//�趨����
	double u = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 2, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	if (mat(2, 0).z() >= 1e-6)
	{
		preFlagZ = 1;
	}
	else if (mat(2, 0).z() < -1e-6)
	{
		preFlagZ = -1;
	}
	while (u <= max)
	{
		u += deltaD;
		nurbssurface.deriveAtH(u, v, 2, mat);
		node2 = GetSurfacePointatUandV(u, v);
		if (mat(2, 0).z() >= 1e-6)
		{
			flagZ = 1;
		}
		else if (mat(2, 0).z() < -1e-6)
		{
			flagZ = -1;
		}

		if (flagZ*preFlagZ < 0)
		{
			if ((*node1 - *node2).GetLength() < 0.001)
			{
				vec.push_back(u);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				u = u - deltaD;
				deltaD = deltaD / 2;
				delete node2;
				continue;
			}

		}
		delete node1;
		node1 = node2;

	}
	return vec;
}


//��ȡv����Ĺյ�,(v����Сֵ��v�����ֵ��uֵ��
std::vector<double> NurbsSurface_qj::getVInflectionPoints(double min, double max, const double u)
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//�趨����
	double v = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 2, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	if (mat(0, 2).z() >= 1e-6)
	{
		preFlagZ = 1;
	}
	else if (mat(0, 2).z() < -1e-6)
	{
		preFlagZ = -1;
	}
	while (v <= max)
	{
		v += deltaD;
		nurbssurface.deriveAtH(u, v, 2, mat);
		node2 = GetSurfacePointatUandV(u, v);
		if (mat(2, 0).z() >= 1e-6)
		{
			flagZ = 1;
		}
		else if (mat(2, 0).z() < -1e-6)
		{
			flagZ = -1;
		}

		if (flagZ*preFlagZ < 0)
		{
			if ((*node1 - *node2).GetLength() < 0.001)
			{
				vec.push_back(v);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				v = v - deltaD;
				deltaD = deltaD / 2;
				delete node2;
				continue;
			}

		}
		delete node1;
		node1 = node2;

	}
	return vec;
}



std::vector<double> NurbsSurface_qj::getUParallelPoints(double min, double max, const double v)
{
	std::vector<double> vec;
	Node *node1 = GetSurfacePointatUandV(min, v);
	Node *node2 = GetSurfacePointatUandV(max, v);
	Node *node3 = GetSurfacePointatUandV((min + max) / 2, v);
	Node node4 = (*node1 - *node3)*(*node2 - *node3);//���Ϊ(min, v)�㣬(max, v)����((min + max)/2, v)��������ƽ��ķ�����l
	Node node5 = (*node1 - *node2)*node4;//*node5Ϊ��(min, v)�㣬(max, v)��ͷ�����l������ƽ��ķ�����
	double length = max - min;
	double deltaD = length / 20;


	Node tempNode1;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	double u = min;

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	nurbssurface.deriveAtH(u, v, 1, mat);
	tempNode1 = Node(mat(1, 0).x(), mat(1, 0).y(), mat(1, 0).z());

	if ((tempNode1 | node5) >= 1e-6)
	{
		preFlagZ = 1;
	}
	else
	{
		preFlagZ = -1;
	}

	while (u <= max)
	{
		u += deltaD;
		nurbssurface.deriveAtH(u, v, 1, mat);
		tempNode1 = Node(mat(1, 0).x(), mat(1, 0).y(), mat(1, 0).z());


		if ((tempNode1 | node5) >= 1e-6)
		{
			flagZ = 1;
		}
		else
		{
			flagZ = -1;
		}
		if (flagZ*preFlagZ < 0)
		{
			if (abs(tempNode1 | node5) < 0.001)
			{
				vec.push_back(u);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				u = u - deltaD;
				deltaD = deltaD / 2;
				continue;
			}

		}
	}
	
	delete node1;
	delete node2;
	delete node3;

	return vec;

}

std::vector<double> NurbsSurface_qj::getVParallelPoints(double min, double max, const double u)
{
	std::vector<double> vec;
	Node *node1 = GetSurfacePointatUandV(u, min);
	Node *node2 = GetSurfacePointatUandV(u, max);
	Node *node3 = GetSurfacePointatUandV(u, (min + max) / 2);
	Node node4 = (*node1 - *node3)*(*node2 - *node3);//���Ϊ(u, min)�㣬(u, max)����(u, (min + max)/2)��������ƽ��ķ�����l
	Node node5 = (*node1 - *node2)*node4;//*node5Ϊ��(u, min)�㣬(u, max)��ͷ�����l������ƽ��ķ�����
	double length = max - min;
	double deltaD = length / 20;


	Node tempNode1;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	double v = min;

	int preFlagZ = 0;//��¼�ϴα�־ֵ������ϴα�־ֵ�ķ�����ôεĲ�ͬ����������������֮��ľ������㾫��Ҫ������Ϊ������е��Ǽ�ֵ��
	int flagZ = 0;//��¼�˴α�־

	nurbssurface.deriveAtH(u, v, 1, mat);
	tempNode1 = Node(mat(1, 0).x(), mat(1, 0).y(), mat(1, 0).z());

	if ((tempNode1 | node5) >= 1e-6)
	{
		preFlagZ = 1;
	}
	else
	{
		preFlagZ = -1;
	}

	while (v <= max)
	{
		v += deltaD;
		nurbssurface.deriveAtH(u, v, 1, mat);
		tempNode1 = Node(mat(1, 0).x(), mat(1, 0).y(), mat(1, 0).z());


		if ((tempNode1 | node5) >= 1e-6)
		{
			flagZ = 1;
		}
		else
		{
			flagZ = -1;
		}
		
		if (flagZ*preFlagZ < 0)
		{
			if (abs(tempNode1 | node5) < 0.001)
			{
				vec.push_back(u);
				preFlagZ = flagZ;
				deltaD = length / 20;
			}
			else
			{
				v = v - deltaD;
				deltaD = deltaD / 2;
				continue;
			}

		}
	}
	
	delete node1;
	delete node2;
	delete node3;

	return vec;

}


//vֵ�̶�����u�������߰����ҳ��ָ�,
std::vector<double> NurbsSurface_qj::getUAveragePerSection(double min, double max, const double v, int splitNum)//vֵ�̶�����u�������߰����ҳ��ָ�
{
	double *u = new double[splitNum + 1];
	double *d = new double[splitNum];
	double *r = new double[splitNum];
	Node *tempNode1 = GetSurfacePointatUandV(min, v);
	Node *tempNode2 = GetSurfacePointatUandV(max, v);
	double totalStringLength = 0;
	double L = 0;
	double initD = (max - min) / splitNum;
	double D = initD;

	u[0] = min;
	for (int i = 0; i < splitNum; ++i)//����ÿ�ε��ҳ������ҳ������ҳ�β�����ҳ�֮�ͣ�
	{
		u[i + 1] = u[i] + D;
		tempNode1 = GetSurfacePointatUandV(u[i], v);
		tempNode2 = GetSurfacePointatUandV(u[i + 1], v);
		L = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		totalStringLength += L;
		delete tempNode1;
		delete tempNode2;
	}
	double R = totalStringLength / splitNum;//����ƽ���ҳ�,���Դ���Ϊ��ʼ�İ뾶
	double deltaR = R / 100;

	double preR = R;
	int r_flag = 0;
	int d_flag = 0;

	double deltaD = D / (4 * splitNum);
	double l = 0;
	int count = 0;

	do
	{
		for (int i = 0; i < splitNum - 1; ++i)
		{
			u[i + 1] = u[i] + D;
			deltaD = D / (4 * splitNum);
			d_flag = 0;
			do
			{
				tempNode1 = GetSurfacePointatUandV(u[i], v);
				tempNode2 = GetSurfacePointatUandV(u[i + 1], v);

				r[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
					(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
					(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
				delete tempNode1;
				delete tempNode2;
				if (r[i] - R > 0.0001)
				{
					if (d_flag == 1)
					{
						deltaD /= 2;
					}
					u[i + 1] = u[i + 1] - deltaD;
					d_flag = -1;
				}
				else if (R - r[i] > 0.0001)
				{
					if (d_flag == -1)
					{
						deltaD /= 2;
					}
					u[i + 1] = u[i + 1] + deltaD;
					d_flag = 1;
				}

			} while (abs(r[i] - R) > 0.0001);
		}

		if (u[splitNum - 1] > max)
		{
			if (r_flag == 1 && deltaR > 1e-4)
			{
				deltaR = deltaR / 2;
			}
			R = R - deltaR;
			r_flag = -1;
			continue;
		}

		tempNode1 = GetSurfacePointatUandV(u[splitNum - 1], v);
		tempNode2 = GetSurfacePointatUandV(max, v);

		l = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		delete tempNode1;
		delete tempNode2;
		if (R - l > 0.001)
		{
			preR = R;
			if (r_flag == 1 && deltaR > 1e-8)
			{
				deltaR = deltaR / 2;
			}
			R = R - deltaR;
			r_flag = -1;
		}
		else if (l - R > 0.001)
		{
			preR = R;
			if (r_flag == -1 && deltaR > 1e-8)
			{
				deltaR = deltaR / 2;
			}
			R = R + deltaR;
			r_flag = 1;
		}

		count++;
	} while (abs(l - R) > 0.001);
	
	std::vector<double> vec;
	
	for (int i = 0; i < splitNum + 1; ++i)
	{
		vec.push_back(u[i]);
	}

	delete[] u;
	delete[] d;
	delete[] r;
	return vec;
}


//uֵ�̶�����v�������߰����ҳ��ָ�
std::vector<double> NurbsSurface_qj::getVAveragePersection(double min, double max, const double u, int splitNum)//uֵ�̶�����v�������߰����ҳ��ָ�
{

	double *v = new double[splitNum + 1];
	double *d = new double[splitNum];
	double *r = new double[splitNum];
	Node *tempNode1 = GetSurfacePointatUandV(u, min);
	Node *tempNode2 = GetSurfacePointatUandV(u, max);
	double totalStringLength = 0;
	double L = 0;
	double initD = (max - min) / splitNum;
	double D = initD;

	v[0] = min;
	for (int i = 0; i < splitNum; ++i)//����ÿ�ε��ҳ������ҳ������ҳ�β�����ҳ�֮�ͣ�
	{
		v[i + 1] = v[i] + D;
		tempNode1 = GetSurfacePointatUandV(u, v[i]);
		tempNode2 = GetSurfacePointatUandV(u, v[i + 1]);
		L = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		totalStringLength += L;
		delete tempNode1;
		delete tempNode2;
	}
	double R = totalStringLength / splitNum;//����ƽ���ҳ�,���Դ���Ϊ��ʼ�İ뾶
	double deltaR = R / 100;

	double preR = R;
	int r_flag = 0;
	int d_flag = 0;

	double deltaD = D / (4 * splitNum);
	double l = 0;
	int count = 0;

	do
	{
		for (int i = 0; i < splitNum - 1; ++i)
		{
			v[i + 1] = v[i] + D;
			deltaD = D / (4 * splitNum);
			d_flag = 0;
			do
			{
				tempNode1 = GetSurfacePointatUandV(u, v[i]);
				tempNode2 = GetSurfacePointatUandV(u, v[i + 1]);

				r[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
					(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
					(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
				delete tempNode1;
				delete tempNode2;
				if (r[i] - R > 0.0001)
				{
					if (d_flag == 1)
					{
						deltaD /= 2;
					}
					v[i + 1] = v[i + 1] - deltaD;
					d_flag = -1;
				}
				else if (R - r[i] > 0.0001)
				{
					if (d_flag == -1)
					{
						deltaD /= 2;
					}
					v[i + 1] = v[i + 1] + deltaD;
					d_flag = 1;
				}

			} while (abs(r[i] - R) > 0.0001);
		}

		if (v[splitNum - 1] > max)
		{
			if (r_flag == 1 && deltaR > 1e-4)
			{
				deltaR = deltaR / 2;
			}
			R = R - deltaR;
			r_flag = -1;
			continue;
		}

		tempNode1 = GetSurfacePointatUandV(u, v[splitNum - 1]);
		tempNode2 = GetSurfacePointatUandV(u, max);

		l = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		delete tempNode1;
		delete tempNode2;
		if (R - l > 0.001)
		{
			preR = R;
			if (r_flag == 1 && deltaR > 1e-8)
			{
				deltaR = deltaR / 2;
			}
			R = R - deltaR;
			r_flag = -1;
		}
		else if (l - R > 0.001)
		{
			preR = R;
			if (r_flag == -1 && deltaR > 1e-8)
			{
				deltaR = deltaR / 2;
			}
			R = R + deltaR;
			r_flag = 1;
		}

		count++;
	} while (abs(l - R) > 0.001);

	std::vector<double> vec;
	for (int i = 0; i < splitNum + 1; ++i)
	{
		vec.push_back(v[i]);
	}
	
	delete[] v;
	delete[] d;
	delete[] r;
	return vec;
}

std::vector<double> NurbsSurface_qj::getUAverage(double min, double max, const double v, int splitNum)
{
	std::vector<double> result;
	std::vector<double> tmp1;//����ÿһ���ĵ��uֵ
	std::vector<double> tmp2;//����ÿһ���ĵ��uֵ
	tmp1 = getUExtremums(min, max, v);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//����ֵ�����
	tmp1.clear();
	tmp1 = getUInflectionPoints(min, max, v);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//���յ����
	tmp1.clear();
	tmp2.push_back(min);//����ʼ�����
	tmp2.push_back(max);//����ֹ�����
	std::sort(tmp2.begin(), tmp2.end());//����ʼ�㡢��ֹ�㡢�յ㰴��С��������
	result.insert(result.end(), tmp2.begin(), tmp2.end());
	for (size_t i = 0; i < result.size() - 1; ++i)
	{
		tmp1 = getUParallelPoints(result[i], result[i + 1], v);//��ȡ�������յ�ƽ�еĵ�
		tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());
		tmp1.clear();
	}
	std::sort(tmp2.begin(), tmp2.end());//����ʼ�㡢��ֹ�㡢�յ��Լ�ƽ�е㰴��С��������


	int * splitNumPerSection = new int[tmp2.size()];//�洢ÿ�ΰ������ֶ��ٶ�
	Node *tempNode1;
	Node *tempNode2;
	double *l = new double[tmp2.size()];//��¼ÿ�ε��ҳ�
	double totalLength = 0;//��¼�ܵĳ���

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		tempNode1 = GetSurfacePointatUandV(tmp2[i], v);
		tempNode2 = GetSurfacePointatUandV(tmp2[i + 1], v);
		l[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		//��������������tempNode1��tempNode2�ľ���
		totalLength += l[i];
		delete tempNode1;
		delete tempNode2;
	}

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		splitNumPerSection[i] = int(l[i] / totalLength*splitNum + 0.5);//����ÿ���ҳ�ռ�ܳ��ȵı�����ȷ���ָ���ٴ�
	}

	result.clear();
	for (size_t i = 0; i < tmp2.size() - 1; ++i)
	{

		tmp1 = getUAveragePerSection(tmp2[i], tmp2[i + 1], v, splitNumPerSection[i]);
		//��ȡÿ�εȷֵ������

		result.insert(result.end(), tmp1.begin(), tmp1.end());
	}
	delete[] l;
	delete[] splitNumPerSection;
	return result;

}

//��vֵ�̶���u�����ϵ����߰����ҳ����֣�type=0��ʾ���µ��ϣ�type = 1��ʾ���ϵ���
//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�u�������ߵȷ�
std::vector<double> NurbsSurface_qj::getUAverageLength(double min, double max, const double v, double lengthPerSection, int type)//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�u�������ߵȷ�
{
	double length = max - min;
	std::vector<double> vec;
	double nextU = min;
	switch (type)
	{
	case 0://0��3����ʾ��������Ѱ�ҵȷֵ�
	case 3:
		nextU = min;
		while (nextU <= max)
		{
			vec.push_back(nextU);
			nextU = nextUPoint(nextU, v, lengthPerSection, type);
		}
		break;
	case 1:
	case 2:
		nextU = max;
		while (nextU >= min)
		{
			vec.push_back(nextU);
			nextU = nextUPoint(nextU, v, lengthPerSection, type);
		}
		break;

	default:
		break;
	}


	return vec;
}


//��uֵ�̶���v�����ϵ����߰����ҳ����֣�type=0��ʾ���µ��ϣ�type = 1��ʾ���ϵ���
//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�v�������ߵȷ�
std::vector<double> NurbsSurface_qj::getVAverageLength(double min, double max, const double u, double lengthPerSection, int type)//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�v�������ߵȷ�
{
	double length = max - min;
	std::vector<double> vec;
	double nextV;
	switch (type)
	{
	case 0://��������
	case 1:
		nextV = min;
		while (nextV <= max)
		{
			vec.push_back(nextV);
			nextV = nextVPoint(u, nextV, lengthPerSection, type);
		}
		break;
	case 2:
	case 3:
		nextV = max;
		while (nextV >= min)
		{
			vec.push_back(nextV);
			nextV = nextVPoint(u, nextV, lengthPerSection, type);
		}
	default:
		break;
	}

	return vec;
}


//�趨���ֶ��ٷݣ�uֵ�̶�����v�������߰����ҳ��ָ�
std::vector<double> NurbsSurface_qj::getVAverage(double min, double max, const double u, int splitNum)
{
	std::vector<double> result;
	std::vector<double> tmp1;//����ÿһ���ĵ��uֵ
	std::vector<double> tmp2;//����ÿһ���ĵ��uֵ
	tmp1 = getVExtremums(min, max, u);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//����ֵ�����
	tmp1.clear();
	tmp1 = getVInflectionPoints(min, max, u);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//���յ����
	tmp1.clear();
	tmp2.push_back(min);//����ʼ�����
	tmp2.push_back(max);//����ֹ�����
	std::sort(tmp2.begin(), tmp2.end());//����ʼ�㡢��ֹ�㡢�յ㰴��С��������
	result.insert(result.end(), tmp2.begin(), tmp2.end());
	for (size_t i = 0; i < result.size() - 1; ++i)
	{
		tmp1 = getVParallelPoints(result[i], result[i + 1], u);//��ȡ�������յ�ƽ�еĵ�
		tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());
		tmp1.clear();
	}
	std::sort(tmp2.begin(), tmp2.end());//����ʼ�㡢��ֹ�㡢�յ��Լ�ƽ�е㰴��С��������


	int * splitNumPerSection = new int[tmp2.size()];//�洢ÿ�ΰ������ֶ��ٶ�
	Node *tempNode1;
	Node *tempNode2;
	double *l = new double[tmp2.size()];//��¼ÿ�ε��ҳ�
	double totalLength = 0;//��¼�ܵĳ���

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		tempNode1 = GetSurfacePointatUandV(u, tmp2[i]);
		tempNode2 = GetSurfacePointatUandV(u, tmp2[i + 1]);
		l[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		//��������������tempNode1��tempNode2�ľ���
		totalLength += l[i];
		delete tempNode1;
		delete tempNode2;
	}

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		splitNumPerSection[i] = int(l[i] / totalLength*splitNum + 0.5);//����ÿ���ҳ�ռ�ܳ��ȵı�����ȷ���ָ���ٴ�
	}

	result.clear();
	for (size_t i = 0; i < tmp2.size() - 1; ++i)
	{

		tmp1 = getVAveragePersection(tmp2[i], tmp2[i + 1], u, splitNumPerSection[i]);
		//��ȡÿ�εȷֵ������

		result.insert(result.end(), tmp1.begin(), tmp1.end());
	}
	delete[] l;
	delete[] splitNumPerSection;
	return result;
}


//Ѱ������ļ���ֵ��
std::vector<PLib::Point2Dd> NurbsSurface_qj::GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV)//Ѱ������ļ���ֵ��
{
	double lengthU = maxU - minU;
	double lenghtV = maxV - minV;

	double deltaU = lengthU / 50; //u�����ʼ����
	double deltaV = lenghtV / 50; //v�����ʼ����

	Node** record = new Node*[50];//��¼u��v���꣬�Լ���Ӧ��zֵ��
	for (int i = 0; i < 50; ++i)
	{
		record[i] = new Node[50];
	}

	double u = minU;
	double v = minV;
	Node * tempNode;

	for (int i = 0; i < 50; ++i)
	{
		u = minU+i * deltaU;

		for (int j = 0; j < 50; ++j)
		{
			v = minV+j * deltaV;
			tempNode = GetSurfacePointatUandV(u, v);
			record[i][j] = Node(u, v, tempNode->getZ());
			delete tempNode;

		}
	}
	std::vector<Node *> vec;
	for (int i = 1; i <49; ++i)
	{
		for (int j = 1; j < 49; ++j)
		{
			if ((record[i][j].getZ() > record[i - 1][j - 1].getZ()) && (record[i][j].getZ() > record[i][j - 1].getZ()) && (record[i][j].getZ() > record[i + 1][j - 1].getZ()) &&
				(record[i][j].getZ() > record[i - 1][j].getZ()) && (record[i][j].getZ() > record[i + 1][j].getZ()) &&
				(record[i][j].getZ() > record[i - 1][j + 11].getZ()) && (record[i][j].getZ() > record[i][j + 1].getZ()) && (record[i][j].getZ() > record[i + 1][j + 11].getZ())
				)//����ֵ��zֵ������Χ8�����zֵ
			{
				vec.push_back(&record[i][j]);
			}
		}
	}
	std::vector<PLib::Point2Dd> result;
	for (size_t i = 0; i < vec.size(); ++i)
	{
		double deltaZ = 0;
		u = vec[i]->getX();
		v = vec[i]->getY();
		minU = vec[i]->getX() - deltaU;
		minV = vec[i]->getY() - deltaV;
		maxU = u + deltaU;
		maxV = v + deltaV;
		deltaU = deltaU / 4;
		deltaV = deltaV / 4;

		double Z = vec[i]->getZ();
		double maxZ = vec[i]->getZ();

		int changedNum = 0;
		while ((u >= minU && u <= maxU) && (v >= minV && v <= maxV))
		{
			changedNum = 0;
			for (double tempU = minU; tempU < maxU; tempU += deltaU)
			{
				for (double tempV = minV; tempV < maxV; tempV += deltaV)
				{
					tempNode = GetSurfacePointatUandV(tempU, tempV);
					if (tempNode->getZ() > maxZ)
					{
						changedNum++;
						maxZ = tempNode->getZ();
						u = tempU;
						v = tempV;
					}
					delete tempNode;
				}
			}
			deltaZ = maxZ - Z;
			Z = maxZ;

			if (changedNum == 0)
			{
				deltaU = deltaU / 2;
				deltaV = deltaV / 2;
			}
			else
			{
				if (deltaZ <= 0.0001)
				{
					result.push_back(PLib::Point2Dd(u, v));
					break;
				}
			}

		}
	}

	for (int i = 0; i < 50; ++i)//���ٶ�̬����Ķ�ά����
	{
		delete[] record[i];
		record[i] = NULL;
	}
	delete[] record;
	record = NULL;
	return result;
}


//Ѱ������ļ�Сֵ��
std::vector<PLib::Point2Dd> NurbsSurface_qj::GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV)//Ѱ������ļ���ֵ��
{
	double lengthU = maxU - minU;
	double lenghtV = maxV - minV;

	double deltaU = lengthU / 50; //u���򲽳�
	double deltaV = lenghtV / 50; //v���򲽳�

	Node** record = new Node*[50];//��¼u��v���꣬�Լ���Ӧ��zֵ��
	for (int i = 0; i < 50; ++i)
	{
		record[i] = new Node[50];
	}

	double u = minU;
	double v = minV;
	Node * tempNode;
	//�Ƚ�������Ի��֣�ȷ����ֵ��Ĵ��·�Χ��Ȼ�������ҵ��ķ�Χ���ö��ַ��Ҽ�ֵ��
	for (int i = 0; i < 50; ++i)//�Ƚ����滮�֣����ֵ������洢�������Է���ȷ����ֵ��ķ�Χ
	{
		u = minU+i * deltaU;

		for (int j = 0; j < 50; ++j)
		{
			v = minV+j * deltaV;
			tempNode = GetSurfacePointatUandV(u, v);
			record[i][j] = Node(u, v, tempNode->getZ());
			delete tempNode;

		}
	}
	std::vector<Node *> vec;
	for (int i = 1; i < 49; ++i)
	{
		for (int j = 1; j < 49; ++j)
		{
			if ((record[i][j].getZ() < record[i - 1][j - 1].getZ()) && (record[i][j].getZ() < record[i][j - 1].getZ()) && (record[i][j].getZ() < record[i + 1][j - 1].getZ()) &&
				(record[i][j].getZ() < record[i - 1][j].getZ()) && (record[i][j].getZ() < record[i + 1][j].getZ()) &&
				(record[i][j].getZ() < record[i - 1][j + 11].getZ()) && (record[i][j].getZ() < record[i][j + 1].getZ()) && (record[i][j].getZ() < record[i + 1][j + 11].getZ())
				)//��Сֵ���zֵС����Χ8�����zֵ
			{
				vec.push_back(&record[i][j]);
			}
		}
	}
	std::vector<PLib::Point2Dd> result;
	double deltaZ = 0;
	for (size_t i = 0; i != vec.size(); ++i)
	{
		u = vec[i]->getX();
		v = vec[i]->getY();
		minU = vec[i]->getX() - deltaU;
		minV = vec[i]->getY() - deltaV;
		maxU = u + deltaU;
		maxV = v + deltaV;
		deltaU = deltaU / 4;
		deltaV = deltaV / 4;

		double Z = vec[i]->getZ();
		double minZ = vec[i]->getZ();

		int changedNum = 0;
		while ((u >= minU && u <= maxU) && (v >= minV && v <= maxV))
		{
			changedNum = 0;
			for (double tempU = minU; tempU < maxU; tempU += deltaU)
			{
				for (double tempV = minV; tempV < maxV; tempV += deltaV)
				{
					tempNode = GetSurfacePointatUandV(tempU, tempV);
					if (tempNode->getZ() < minZ)
					{
						changedNum++;
						minZ = tempNode->getZ();
						u = tempU;
						v = tempV;
					}
					delete tempNode;
				}
			}
			deltaZ = Z - minZ;
			Z = minZ;

			if (changedNum == 0)
			{
				deltaU = deltaU / 2;
				deltaV = deltaV / 2;
			}
			else
			{
				if (deltaZ <= 0.0001)
				{
					result.push_back(PLib::Point2Dd(u, v));
					break;
				}
			}
		}
	}

	for (int i = 0; i < 50; ++i)//���ٶ�̬����Ķ�ά����
	{
		delete[] record[i];
		record[i] = NULL;
	}
	delete[] record;
	record = NULL;
	return result;
}


//�����򣺻������棬���ظ����ֿ���ڲ���ͱ߽�㣬�ȴ��߽紦��
std::vector< Member> NurbsSurface_qj::surfaceAverageTotalSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> &Node_uv)
{
	
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ�Сֵ��
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ���ֵ��
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//�洢��ֵ��

	int num = 0;//��¼��ֵ����

	for (std::vector<PLib::Point2Dd>::iterator itr = pMinExtrmums.begin(); itr != pMinExtrmums.end(); ++itr)//����Сֵ���뼫ֵ������
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	
	for (std::vector<PLib::Point2Dd>::iterator itr = pMaxExtrmums.begin(); itr != pMaxExtrmums.end(); ++itr)//������ֵ���뼫ֵ������
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	std::sort(extrmums, extrmums + num, compare_CPoint2D);//�Լ�ֵ�����鰴vֵ�Ӵ�С����compare_CPoint2D�Ǹ�������
	//std::sort�Զ���ṹ
	std::vector <PointStruct> uvEdge;//���ڼ�¼�߽���uvֵ
	//std::vector <PointStruct> vVec;//���ڼ�¼��v����ı߽�㣬���ǵı�ʶ��value��¼uֵ��С
	//equalDivPoints divPoint;		//����ֵ�����ڵ��ݺ�ֱ���ϵ����еȷֵ�洢����,equalDivPoints��value����ʶ�����������
	PointStruct divPoint;//������ʱ��¼��ֵ��	
	divPoint.Id = -1;
	//Node AddSpace;//��������UVForOut
	std::vector<double> temp;//������ʱ��¼�ȷֵ�
	//����߽�ȷֵ�
	//�±߽�
	temp = getUAverageLength(minU, maxU,minV, splitLengthU, 0);//��������
	divPoint.point.y() = minV;
	//AddSpace.setY(minV);
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.x() = temp[i];	
		//divPoint.Id = nId;
		uvEdge.push_back(divPoint);
		//AddSpace.setX(temp[i]);
		//AddSpace.setID(nId);
		//UVForOut.push_back(AddSpace);
		//++nId;
	}
	if ((maxU - (uvEdge.end() - 1)->point.x()) < 0.001)//�ж������Ƿ񿿵�̫��
	{
		uvEdge.pop_back();
	}

	//�ұ߽�
	temp = getVAverageLength(minV, maxV, maxU, splitLengthV, 1);//��������	
	divPoint.point.x() = maxU;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.y() = temp[i];	
		uvEdge.push_back(divPoint);
	}
	if ((maxV - (uvEdge.end() - 1)->point.y()) < 0.001)//�ж������Ƿ񿿵�̫��
	{
		uvEdge.pop_back();
	}

	//�ϱ߽�
	temp = getUAverageLength(minU,maxU,maxV, splitLengthU, 1);//��������	
	divPoint.point.y() = maxV;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.x() = temp[i];
		uvEdge.push_back(divPoint);
	}
	if (((uvEdge.end() - 1)->point.x()-minU ) < 0.001)//�ж������Ƿ񿿵�̫��
	{
		uvEdge.pop_back();
	}

	//��߽�
	temp = getVAverageLength(minV, maxV, minU, splitLengthV, 2);//��������
	divPoint.point.x() =minU;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.y() = temp[i];
		uvEdge.push_back(divPoint);
	}
	if (((uvEdge.end() - 1)->point.y()-minV ) < 0.001)//�ж������Ƿ񿿵�̫��
	{
		uvEdge.pop_back();
	}

	for (int i = 0; i < num; ++i)//��ֵ���vֵ�̶�����u����ĵȷֵ�ȫ����¼�������Է��㻮���߽�Ĵ���num��ʾ��ֵ������
	{
		divPoint.point.y() = extrmums[i].y;//��¼��ֵ���vֵ
		temp = getUAverageLength(minU, extrmums[i].x, extrmums[i].y, splitLengthU, 1);//��������ȷ֣���Ӧtype=1��2,��¼�ȷֵ㣬��std::vector
		temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ���滹�����һ��
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.x() = temp[j];//�����temp��u���PointStruct���������uvEdge��
			uvEdge.push_back(divPoint);//���ȷֱ߽�����uvEdge
		}

		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());//��temp���뵽divPoint��

		temp = getUAverageLength(extrmums[i].x, maxU, extrmums[i].y, splitLengthV, 0);//�������ҵȷ֣���Ӧtype=0,3
		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.x() = temp[j];//���Ұ�temp��u���PointStruct���������uvEdge��			
			uvEdge.push_back(divPoint);
		}
		//std::sort(uvEdge[i].points.begin(), uvEdge[i].points.end());//��������uֵ��С����
	}

	
	//�˴�compare��surface.cpp���ж��壬��˸�Ϊcompare1
	for (int i = 0; i < num; ++i)//��ֵ���uֵ�̶�����v����ĵȷֵ�ȫ����¼�������Է��㻮���߽�Ĵ���num��ʾ��ֵ������
	{
		divPoint.point.x() = extrmums[i].x;//��¼��ֵ���vֵ
		temp = getVAverageLength(minV, extrmums[i].y, extrmums[i].x, splitLengthU, 2);//��������ȷ֣���Ӧtype=1��2,��¼�ȷֵ㣬��std::vector
		temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ�����Ѿ����
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.y() = temp[j];//���°�temp��v���PointStruct���������vVec��
			uvEdge.push_back(divPoint);//���ȷֱ߽�����vVec
		}

		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());//��temp���뵽divPoint��

		temp = getVAverageLength(extrmums[i].y, maxV, extrmums[i].x, splitLengthV, 0);//�������ҵȷ֣���Ӧtype=0,3
		temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ�����Ѿ����
		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.y() = temp[j];//���ϰ�temp��v���PointStruct���������vVec��			
			uvEdge.push_back(divPoint);
		}

	}
	PLib::Point2Dd * extrmums2Dd = new PLib::Point2Dd[num];

	for (int i = 0; i < num; ++i)//����ֵ�������ת��ΪPointDd�洢
	{
		extrmums2Dd[i].x() = extrmums[i].x;
		extrmums2Dd[i].y() = extrmums[i].y;
	}
	//std::sort(extrmums2Dd, extrmums2Dd + num, compare_CPoint2D);//����ֵ�㰴��v�Ӵ�С�����Է��㻮���߽�Ĵ���
	PBSTreeNode BSTree = NULL;
	BiTree_Create(&BSTree, extrmums2Dd, num);//���ݼ�ֵ�㽨�����Ҷ�����
	rectangleRegion * rec = new rectangleRegion[3 * num + 1];//num����ֵ��Ҫ����3*num+1������
	PreOrder(BSTree, rec, 3 * num + 1, minU, maxU, minV, maxV);//����ֵ�㽫���飬�������Ϣ������rec������

	std::vector< Member> mVec;//���ڼ�¼�ܵĸ����ֿ��ڵ�member
	
	std::vector<Member> mVecTemp;//��ʱ�����ڼ�¼ѭ���зֿ��ڵ�member���β�
	//std::vector<Node> nVecTemp;//��ʱ���βΣ�Node

	for (int i =0; i < 3*num+1; ++i)//������������еȷ֣������ȷֵĵ�͸���Ϣ�ֱ�洢��nVec��mVec�У�3*num+1Ϊ�ֿ���
		//***************************************************************************************************7���޸�
	{
		mVecTemp = surfaceAverageSpit(rec[i], splitLengthU, splitLengthU, uvEdge);//��i����еȷ�
		(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//��mVecTemp�ӵ�mVec��

		//mVec.push_back(mVecTemp);//����ᴫ�ݵ�main��
		//nVec.push_back(nVecTemp);//����ᴫ�ݵ�main��
		mVecTemp.clear();
		//nVecTemp.clear();
	}
	
	//�߽���и˼�����*********************************************
	std::vector<PLib::Point2Dd> extrmumVec;
	for (size_t i = 0; i < num; i++)
	{
		extrmumVec.push_back(extrmums2Dd[i]);
	}
	mVecTemp = EdgeTotal_Member(uvEdge, extrmumVec, num, minU, maxU, minV, maxV, splitLengthU, splitLengthV);
	(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//��mVecTemp�ӵ�mVec��
	//�߽���и˼�����*********************************************
	
	delete[] extrmums;
	delete[] rec;
	delete[] extrmums2Dd;
	
	Node_uv = UVForOut;
	return mVec;

}


//����һ���u��v���꣬����v���䣬Ѱ��u�������ҳ�����uChordLength��uֵ,type��ʾ��������
double NurbsSurface_qj::nextUPoint(double u, const double v, double uChordLength, int type, double error )
{

	double deltaU = 0.05;
	double nextU = u;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	switch (type)
	{
	case 0://�����µ�����
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r�ϴ���Ҫ��С
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}
			else if (uChordLength - r > error)//r��С����Ҫ����
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}

		} while (abs(r - uChordLength) > error);

		break;
	case 1://�����µ�����
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r�ϴ�u��Ҫ����
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}
			else if (uChordLength - r > error)//r��С��u��Ҫ��С
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}

		} while (abs(r - uChordLength) > error);
		break;
	case 2: //�����ϵ�����
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r�ϴ�u��Ҫ����
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}
			else if (uChordLength - r > error)//r��С��u��Ҫ��С
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}

		} while (abs(r - uChordLength) > error);
		break;
	case 3://�����ϵ�����
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r�ϴ���Ҫ��С
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}
			else if (uChordLength - r > error)//r��С����Ҫ����
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}

		} while (abs(r - uChordLength) > error);
		break;
	default:
		break;
	}


	return nextU;
}


//����һ���u��v���꣬Ѱ��v�������ҳ�����vChordLength��uֵ,type��ʾ��������
double NurbsSurface_qj::nextVPoint(const double u, double v, double vChordLength, int type, double error )
{
	double deltaV = 0.05;
	double nextV = v;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	switch (type)
	{
	case 0://�������������ҵ�
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r�ϴ���Ҫ��С
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}
			else if (vChordLength - r > error)//r��С����Ҫ����
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}

		} while (abs(r - vChordLength) > error);
		break;
	case 1://������������Ѱ��
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r�ϴ���Ҫ��С
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}
			else if (vChordLength - r > error)//r��С����Ҫ����
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}

		} while (abs(r - vChordLength) > error);
		break;
		break;
	case 2://�����ϵ�����
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r�ϴ�v��Ҫ����
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}
			else if (vChordLength - r > error)//r��С��v��Ҫ��С
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}

		} while (abs(r - vChordLength) > error);
		break;

	case 3://�����ϵ�����
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r�ϴ�v��Ҫ����
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}
			else if (vChordLength - r > error)//r��С��v��Ҫ��С
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}

		} while (abs(r - vChordLength) > error);
		break;
		break;
	default:
		break;
	}


	return nextV;
}


//�������㣬����uChordLength��vChordLength������㣨u��v����type��ʾ��������
PLib::Point2Dd NurbsSurface_qj::nextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type,double error)
//
{
	if ((point2.x() == -1 && point2.y() == -1) || (point1.x() == -1 && point1.y() == -1))
	{
		return PLib::Point2Dd(-1, -1);
		
	}
	//cout<<point1.x()<<"	"<<point1.y()<<"	"<<point2.x()<<"	"<<point2.y()<<endl;
	int Uflag = 0;
	int Vflag = 0;
	int preUFlag = 0;
	int preVFlag = 0;
	double u = point2.x();
	double v = point1.y();

	double deltaU = 0.005;
	double deltaV = 0.005;
	Node *tempNode1;
	Node *tempNode2;
	Node *tempNode3;
	double r1 = 0;
	double r2 = 0;
	tempNode1 = GetSurfacePointatUandV(point1.x(), point1.y());
	tempNode2 = GetSurfacePointatUandV(point2.x(), point2.y());
	unsigned int count = 0;
	switch (type)
	{
	case 0://�������������ҵ�
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2�ϴ���Ҫ��С
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= 0)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < 1e-10 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 > 0)//r2��С����Ҫ����
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1�ϴ���Ҫ��С
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= 0)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < 1e-10 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 > 0)//r2��С����Ҫ����
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 1://������������Ѱ��
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2�ϴ���Ҫ��С
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 > 0)//r2��С����Ҫ����
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1��С����Ҫ��Сuʹr1����
				   {
					   if (Uflag == 1)//˵����һ����Ҫ����u�����˵���Ҫ��Сu�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= -0.3)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < -0.3 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 < 0)//r1�ϴ���Ҫ����uʹr1��С
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 2://�����ϵ�����
	{

			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2��С����Ҫ��Сvʹr2����
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 < 0)//r2�ϴ���Ҫ����vʹr2��С
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1��С����Ҫ��Сuʹr1����
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= -0.3)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < -0.3 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 < 0)//r1�ϴ���Ҫ����u��r1��С
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 3://�����ϵ�����
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2��С����Ҫ��Сvʹr2����
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 < 0)//r2�ϴ���Ҫ����vʹr2��С
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1�ϴ���Ҫ��С
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= 0)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < 1e-10 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 > 0)//r2��С����Ҫ����
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;
						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;

	default:
		break;
	}
	delete tempNode1;
	delete tempNode2;
	return PLib::Point2Dd(u, v);
}
//extern std::vector<equalDivPoints> uvEdge;	//��¼u������ֵ�����
//extern std::vector<equalDivPoints> vVec;	//��¼v������ֵ�����



//�ֿ��ڲ����ҳ����֣�pDisplayVec���ڴ洢Ҫ��ʾ�ĵ㣬���ݵ�surfaceAverageTotalSpit�е�nVecTemp��Ȼ���ٴ��ݵ�main�е�nVec
//mDisplayVec���ڴ洢Ҫ��ʾ�ĸˣ����ݵ�surfaceAverageTotalSpit�е�mVecTemp��Ȼ���ٴ��ݵ�main�е�mVec
//����mDisplayVec
std::vector<Member> NurbsSurface_qj::surfaceAverageSpit(rectangleRegion &rec, double splitLengthU, double splitLengthV, std::vector<PointStruct> &uvEdge)
{
	double maxU = rec.maxU + 0.3;
	//***********�˴�0.3�д��Ż�***************************************************
	double maxV = rec.maxV + 0.3;
	double minU = rec.minU - 0.3;
	double minV = rec.minV - 0.3;

	PointStruct ** points = new PointStruct *[1000];
	for (int i = 0; i < 1000; ++i)
	{
		points[i] = new PointStruct[1000];
		for (int j = 0; j < 1000; ++j)
		{
			points[i][j].point = PLib::Point2Dd(-1, -1);
			points[i][j].Id=-1;
		}
	}
	int row0 = 0;
	int col0 = 0;
	int row1 = 1;
	int col1 = 0;

	double nextU = 0;
	double nextV = 0;
	PLib::Point2Dd nextCrossP;
	Node *tempNode1;
	Node *tempNode2;
	Node *tempNode3;
	double cosValue = 0;
	double u = 0;
	double v = 0;
	switch (rec.type)
	{
	case 0://�������������ҵ�
		u = rec.minU;
		v = rec.minV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u <= maxU)//����һ�е�������ҵ�
		{
			nextU = nextUPoint(u, v, splitLengthU);
			points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
			col0 = col0 + 1;
			u = nextU;
		}

		v = rec.minV;

		for (row0 = 0, row1 = 1; (row0 < 1000) && (v <= maxV); row0++, row1++)
		{
			u = rec.minU;
			col0 = 0;
			col1 = 0;
			nextV = nextVPoint(u, v, splitLengthV);//Ѱ��ÿ�еĵ�һ����v����
			if (nextV <= maxV)
			{
				points[row1][col1].point = PLib::Point2Dd(u, nextV);
			}

			while (u < maxU)
			{
				tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());
				tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());
				tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());

				cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
				delete tempNode1;
				delete tempNode2;


				nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV);

				if (nextCrossP.x() == -1 && nextCrossP.y() == -1 && points[row1][col1].point.x() != -1 && points[row1][col1].point.y() != -1
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//���1 2 ���������������ҵ�������Ϊ-1 -1����������
				{
					if (cosValue < -0.866)//����1 2�������ҵ�����û���ҵ��������㹹�ɵ������εĽǶȴ���150�ȣ�����������ε������������һ��
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV);
						if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						continue;
					}
					else
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV);
						if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}

				}
				if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
				{
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));
					points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
				}
				else
				{
					break;
				}
				u = points[row0][col0 + 1].point.x();
				col0 = col0 + 1;
				col1 = col1 + 1;
			}
			v = nextV;
		}

		break;

	case 1://������������Ѱ��
		u = rec.maxU;
		v = rec.minV;

		points[0][0].point = PLib::Point2Dd(u, v);
		while (u >= minU)//����һ�е�������ҵ�
		{
			nextU = nextUPoint(u, v, splitLengthU, 1);
			points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
			col0 = col0 + 1;
			u = nextU;
		}
		v = rec.minV;

		for (row0 = 0, row1 = 1; (row0 < 1000) && (v <= maxV); row0++, row1++)
		{
			u = rec.maxU;
			col0 = 0;
			col1 = 0;
			nextV = nextVPoint(u, v, splitLengthV, 1);//Ѱ��ÿ�еĵ�һ����v����
			if (nextV <= maxV)
			{
				points[row1][col1].point = PLib::Point2Dd(u, nextV);
			}

			while (u >= minU)
			{
				tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());
				tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());
				tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());

				cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
				delete tempNode1;
				delete tempNode2;

				nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 1);

				if (nextCrossP.x() == -1 && nextCrossP.y() == -1 && points[row1][col1].point.x() != -1 && points[row1][col1].point.y() != -1
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//���1 2 ���������������ҵ�������Ϊ-1 -1����������
				{
					if (cosValue < -0.866)//����1 2�������ҵ�����û���ҵ��������㹹�ɵ������εĽǶȴ���150�ȣ�����������ε������������һ��
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 1);
						if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						//col0 = col0 + 1;
						continue;
					}
					else
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 1);
						if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						col0 = col0 + 2;
						u = points[row0][col0 + 2].point.x();
						continue;
					}

				}
				if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (rec.minV <= nextCrossP.y() && nextCrossP.y() <= maxV))
				{
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));
					points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
				}
				else
				{
					break;
				}
				u = points[row0][col0 + 1].point.x();
				col0 = col0 + 1;
				col1 = col1 + 1;
			}
			v = nextV;
		}

		break;

	case 2://��������������
		u = rec.maxU;
		v = rec.maxV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u >= minU)//����һ�е�������ҵ�
		{
			nextU = nextUPoint(u, v, splitLengthU, 2);
			points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
			//points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			col0 = col0 + 1;
			u = nextU;
		}
		v = rec.maxV;

		for (row0 = 0, row1 = 1; (row0 < 1000) && (v >= minV); row0++, row1++)
		{
			u = rec.maxU;
			col0 = 0;
			col1 = 0;
			nextV = nextVPoint(u, v, splitLengthV, 2);//Ѱ��ÿ�еĵ�һ����v����
			if (nextV >= minV)
			{
				points[row1][col1].point = PLib::Point2Dd(u, nextV);
				//points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			}

			while (u >= minU)
			{
				tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());
				tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());
				tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());

				cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
				delete tempNode1;
				delete tempNode2;

				nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 2);

				if (nextCrossP.x() == -1 && nextCrossP.y() == -1 && points[row1][col1].point.x() != -1 && points[row1][col1].point.y() != -1
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//���1 2 ���������������ҵ�������Ϊ-1 -1����������
				{
					if (cosValue < -0.866)//����1 2�������ҵ�����û���ҵ��������㹹�ɵ������εĽǶȴ���150�ȣ�����������ε������������һ��
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 2);
						if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (minV<= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						//col0 = col0 + 1;
						continue;
					}
					else
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 2);
						if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (minV<= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}

					//if ((nextCrossP.x() == -1 && nextCrossP.y() == -1 && points[row1][col1].x() != -1 && points[row1][col1].y() != -1
					//	&& points[row1][col1].x() != -1 && points[row1][col1].y() != -1))//�����������Ҳ��������㣬����1�����������u����������һ��
					//{
					//	nextCrossP.x() = points[row0][col0].x();
					//	nextCrossP.y() = nextVPoint(points[row0][col0].x(), points[row1][col1].y(), splitLengthV);

					//}

				}
				if ((minU <= nextCrossP.x() && nextCrossP.x() <= rec.maxU) && (minV <= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
				{
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));
					points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
				}
				else
				{
					break;
				}
				u = points[row0][col0 + 1].point.x();
				col0 = col0 + 1;
				col1 = col1 + 1;
			}
			v = nextV;
		}

		break;

	case 3://��������������
		u = rec.minU;
		v = rec.maxV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u <= maxU)//����һ�е�������ҵ�
		{
			nextU = nextUPoint(u, v, splitLengthU, 3);
			points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
			col0 = col0 + 1;
			u = nextU;
		}
		v = rec.maxV;

		for (row0 = 0, row1 = 1; (row0 < 1000) && (v >= minV); row0++, row1++)
		{
			u = rec.minU;
			col0 = 0;
			col1 = 0;
			nextV = nextVPoint(u, v, splitLengthV, 3);//Ѱ��ÿ�еĵ�һ����v����
			if (nextV >= minV)
			{
				points[row1][col1].point = PLib::Point2Dd(u, nextV);
				//points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			}

			while (u <= maxU)
			{
				tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());
				tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());
				tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());

				cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
				delete tempNode1;
				delete tempNode2;

				nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 3);

				if (nextCrossP.x() == -1 && nextCrossP.y() == -1 && points[row1][col1].point.x() != -1 && points[row1][col1].point.y() != -1
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//���1 2 ���������������ҵ�������Ϊ-1 -1����������
				{
					if (cosValue < -0.866)//����1 2�������ҵ�����û���ҵ��������㹹�ɵ������εĽǶȴ���150�ȣ�����������ε������������һ��
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 3);
						if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (minV <= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						continue;
					}
					else
					{
						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 3);
						if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (minV <= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
						{
							points[row1][col1 + 1].point = nextCrossP;
							points[row0][col0].next.push_back(pointCoordinate(row1, col1));
							points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
							points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						}
						else
						{
							break;
						}
						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}

				}
				if ((rec.minU <= nextCrossP.x() && nextCrossP.x() <= maxU) && (minV <= nextCrossP.y() && nextCrossP.y() <= rec.maxV))
				{
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));
					points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
				}
				else
				{
					break;
				}
				u = points[row0][col0 + 1].point.x();
				col0 = col0 + 1;
				col1 = col1 + 1;
			}
			v = nextV;
		}
		break;
	default:
		break;
	}


	//Node ** nodeForDisplay = new Node *[1000];//�洢��ʾ������㣬��¼u��vֵ�����չ��ܲ���Ҫ��������pDisplayVec
	////nodeForSave���û��
	//Node ** nodeForSave = new Node *[1000];//�洢����㣬Ϊ�˷���˲��ң����ID˳������������е�˳��һ��
	//for (int i = 0; i < 1000; ++i)
	//{
	//	nodeForDisplay[i] = new Node[1000];
	//	nodeForSave[i] = new Node[1000];
	//	for (int j = 0; j < 1000; ++j)
	//	{
	//		nodeForDisplay[i][j].setX(-1);
	//		nodeForDisplay[i][j].setY(-1);
	//		nodeForDisplay[i][j].setID(-1);
	//		nodeForSave[i][j].setX(-1);
	//		nodeForSave[i][j].setY(-1);
	//		nodeForSave[i][j].setZ(-1);
	//	}
	//}


	//unsigned long long int pID = 0;//��ID
	//unsigned long long int mID = 0;//��ID
	//std::vector<Node> pSaveVec;//�洢������������˳����ID��ȫһ�£�Node��x��y��z�����ʾʵ��ֵ
	std::vector<Member> mDisplayVec; //��Ϣ�����洢�˵��������˵�˳�����ID��ȫһ��
	std::vector<Member> mSaveVec; //�洢�˵��������˵�˳�����ID��ȫһ��
	Node addspace;//��������NodeForOut�����Ĵ�С
	//std::vector<double> uDirPoints;//��¼u�����ϲ�����߽紦�ı߽�ȷֵ�
	//std::vector<double> vDirPoints;//��¼v�����ϲ�����߽紦�ı߽�ȷֵ�
	//std::vector<PointStruct> uv_edge;//��ʱ��¼�÷ֿ�ı߽��,�Ա���С������Χ
	std::vector<PointStruct> uvedge_up, uvedge_down, uvedge_left, uvedge_right;//�洢�ֿ��ıߵı߽��
	std::vector<PointStruct*> points_up, points_down,points_left,points_right;//�洢���������߽���ڲ����ָ��
	std::vector<PointStruct*> point_corner;//��¼�ǵ�


	switch (rec.type)
	{
	case 0://�����µ�����
	{			  
		//�Ը÷ֿ�ı߽�����������**********************************************************************************
		//���߽���uv����Ŷ������������UVForOut��

		//���طֿ���������ұ߽������������Ա߽����������ţ����߽����뵽�����
		Edge_udlr(uvEdge, rec,uvedge_up,uvedge_down,uvedge_left,uvedge_right);
		std::sort(uvedge_down.begin(), uvedge_down.end(), compare1_x);//�Ա߽���С��������
		std::sort(uvedge_up.begin(), uvedge_up.end(), compare1_x);
		std::sort(uvedge_left.begin(), uvedge_left.end(), compare3_y);
		std::sort(uvedge_right.begin(), uvedge_right.end(), compare3_y);
					
		//�Էֿ��ڲ��ĵ����������,���벻����߽�������ڲ���,�������벻����߽磨���ϱ߽磩������ڲ�������
		for (int i = 1; i < 999; ++i)
		{
			for (int j = 1; j < 999; ++j)
			{
				if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				{
					if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() < rec.maxV))//*****************************��䶯
						//������Ҫ��������ڵĵ�˳��洢	
					{
						//�ж��Ƿ�Ϊ������������߽�ĵ�**********************************��Ķ�
						if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() <rec.maxU) && (points[i][j + 1].point.y() <rec.maxV))
						{					
							points[i][j].Id = nId;
							addspace.setX(points[i][j].point.x());
							addspace.setY(points[i][j].point.y());
							addspace.setID(nId);
							UVForOut.push_back(addspace);
							++nId;
							continue;
						}
						//�ж��Ƿ�Ϊ�����������ϱ߽�ĵ�********************************��Ķ�
						//��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��,���Ǽ����ֻ����ͨ�����
						if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() >= rec.maxU))
						{
							point_corner.push_back ( &points[i][j]);//����ǵ�
							point_corner.push_back(&points[i][j-1]);//1
							point_corner.push_back(&points[i-1][j]);//2
							point_corner.push_back(&uvedge_right[uvedge_right.size() - 2]);//3
							point_corner.push_back(&uvedge_right.back());//4
							point_corner.push_back(&uvedge_up[uvedge_up.size() - 2]);//5
							point_corner.push_back(&uvedge_up.back());//6
						}
						else if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() < rec.maxU))
						{
							points_up.push_back(&points[i][j]);//�ϱ߽總�����ڲ���
						}
						else if ((points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() >= rec.maxU))
						{
							points_right.push_back(&points[i][j]);//�ұ߽總�����ڲ���
						}
						//�ж��Ƿ�Ϊ�����������ұ߽�ĵ�********************************��Ķ�
						//��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��
						
					}
				}				
			}
		}
		
		std::sort(points_up.begin(), points_up.end(), edge_compare1);//��point_up�еĵ㰴u��С��������
		std::sort(points_right.begin(), points_right.end(), edge_compare3);//��point_right�еĵ㰴v��С��������
		

		//�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
		for (int j = 0; j < 999; ++j)
		{
			if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
			{
				if (points[0][j].point.x() <=rec.maxU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				{
					for (size_t k = 0; k < uvedge_down.size(); ++k)
					{
						if ((points[0][j].point.x() == uvedge_down[k].point.x()) && (points[0][j].point.y() == uvedge_down[k].point.y()))
						{
							points[0][j].Id = uvedge_down[k].Id;
							break;
						}
					}
				}
			}
		}
		//�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
		for (int i = 1; i < 999; ++i)
		{
			if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
			{
				if (points[i][0].point.y() <= rec.maxV)//ȷ�����ڷֿ��ڣ�Ҫ�ߵ㣬������=//*****************************��䶯
				{
					for (size_t k = 0; k < uvedge_left.size(); ++k)
					{
						if ((points[i][0].point.x() == uvedge_left[k].point.x()) && (points[i][0].point.y() == uvedge_left[k].point.y()))
						{
							points[i][0].Id = uvedge_left[k].Id;
							break;
						}
					}
				}
			}
		}

		//�Ը÷ֿ�ı߽�����������***************************************************************************************
		//����Ĺ�ϵ��nId��¼����********************************************************************************************

		//double up_value = 2, down_value = 2;//��¼�߽����ڵ���߽��ľ���,2�ǳ�ʼֵ
		//int up_index = 0, down_index = 0;//��¼�߽����ڵ�������ı߽��ı��,�˱��Ϊuv_edge���±�
		Member m;//���ڼ�¼����Ϣ
		std::vector<Member> edge_m;//��¼������߽紦�ĸ˼�

		//���ȴ���������Ͻǵ�ı߽磬��ʱ��Ϊ�ϱ߽磬�������Ա�֤�ǵ�����
		//�ϱ߽紦��
		edge_m = Edge_Member(uvEdge,uvedge_up, points_up, 1);//x��������
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
		//�ұ߽紦��
		edge_m = Edge_Member(uvEdge, uvedge_right, points_right,  3);//y��������
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
		//�ǵ㴦��
		edge_m = EdgeCorner_Member(point_corner);
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
		
		//���±߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
		for (int j = 1; j < 999; ++j)
		{
			if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
			{
				if ((points[0][j].point.x() < rec.maxU) && (points[1][j].point.x()<rec.maxU))//ȷ����[0],[1]�ڷֿ��ڣ�//*****************************��䶯
				{
					for (size_t k = 0; k < points[0][j].next.size(); ++k)
					{
						if (points[0][j].next[k].x == 1)
						{
							m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//���������γɸ˼�
							mDisplayVec.push_back(m);//�˼�����mDisplayVec						
						}
					}
				}
			}
		}

		//����߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
		for (int i = 1; i < 999; ++i)
		{
			if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
			{
				if ((points[i][0].point.y() < rec.maxV)&&(points[i][1].point.y()<rec.maxV))//ȷ����[0],[1]�ڷֿ��ڣ�//*****************************��䶯
				{
					for (size_t k = 0; k < points[i][0].next.size(); ++k)
					{
						if (points[i][0].next[k].y == 1)
						{
							m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//���������γɸ˼�
							mDisplayVec.push_back(m);//�˼�����mDisplayVec		
						}
					}
				}
			}
		}
		
		//�����ڲ��㣬���������±߽�ĵ㣬�������������ϱ߽���ڲ��㣨�����������ϵ��Ҫ�������˴�i,j��ʼΪ1
		for (int i = 1; i < 999; ++i)
		{
			for (int j = 1; j < 999; ++j)
			{
				if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				{
					if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() < rec.maxV))//�жϵ��Ƿ�������ֿ���//*****************************��䶯
					{
						//Ѱ�����ұ߽�������ڲ���						
						if (points[i][j + 1].point.x() >= rec.maxU)//�ж�j��u��vС���ұ߽磬j+1��u���ڵ����ұ߽�//*****************************��䶯
						{							
							if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//�жϵ���Ϸ����Ƿ��ڷֿ��ڣ��ж���������
							{	
								m = Add_2Member(points[i][j], points[i+1][j]);//���������γɸ˼�
								mDisplayVec.push_back(m);//�˼�����mDisplayVec																						
							}
						}
						if (points[i + 1][j].point.y() >= rec.maxV)//Ѱ�����ϱ߽�������ڲ���//*****************************��䶯
						{							
							if ((points[i][j + 1].point.x() < rec.maxU) && (points[i][j + 1].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//�жϵ���ҷ����Ƿ��ڷֿ���
							{
								m = Add_2Member(points[i][j], points[i][j + 1]);//���������γɸ˼�
								mDisplayVec.push_back(m);//�˼�����mDisplayVec								
							}
						}
						if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() <rec.maxU) && (points[i][j + 1].point.y() <rec.maxV))
						{
							for (size_t k = 0; k < points[i][j].next.size(); ++k)
							{
								m = Add_2Member(points[i][j], points[points[i][j].next[k].x][points[i][j].next[k].y]);
								mDisplayVec.push_back(m);								
							}
						}
					}
				}
			}
		}
			   
		
			   break;
	}
	//����Ĺ�ϵ��nId��¼����********************************************************************************************



	case 1://�����µ�����
	{
			   //�Ը÷ֿ�ı߽�����������**********************************************************************************
			   //���߽���uv����Ŷ������������UVForOut��
			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare2_x);//�Ա߽��Ӵ�С����
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare2_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare3_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare3_y);

			   //�Էֿ��ڲ��ĵ����������
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() >rec.minU) && (points[i][j].point.y() < rec.maxV))//*****************************��䶯
							   //������Ҫ��������ڵĵ�˳��洢	
						   {							 
							   if ((points[i + 1][j].point.x()> rec.minU) && (points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y() <rec.maxV))
							   {
								   points[i][j].Id = nId;
								   addspace.setX(points[i][j].point.x());
								   addspace.setY(points[i][j].point.y());
								   addspace.setID(nId);
								   UVForOut.push_back(addspace);
								   ++nId;
							   }
							   //�ж��Ƿ�Ϊ�����������ϱ߽�ĵ�********************************��Ķ�
							   //��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��,���Ǽ����ֻ����ͨ�����
							   //if (points[i + 1][j].point.y() >= rec.maxV)
							   //{
								  // points_up.push_back(&points[i][j]);
							   //}
							   ////�ж��Ƿ�Ϊ�����������ұ߽�ĵ�********************************��Ķ�
							   ////��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��
							   //if (points[i][j + 1].point.x() <= rec.minU)
							   //{
								  // points_left.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   //point_corner = &points[i][j];//����ǵ�
								   point_corner.push_back(&points[i][j]);//����ǵ�
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_left[uvedge_left.size() - 2]);//3
								   point_corner.push_back(&uvedge_left.back());//4
								   point_corner.push_back(&uvedge_up[uvedge_up.size() - 2]);//5
								   point_corner.push_back(&uvedge_up.back());//6
							   }
							   else if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() >rec.minU))
							   {
								   points_up.push_back(&points[i][j]);//�ϱ߽總�����ڲ���
							   }
							   else if ((points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() < rec.minU))
							   {
								   points_left.push_back(&points[i][j]);//��߽總�����ڲ���
							   }
						   }
					   }
				   }
			   }
			   //�Բ�����߽��Ե��ڲ��������������
			   std::sort(points_up.begin(), points_up.end(),edge_compare2);//��point_up�еĵ㰴u�Ӵ�С����
			   std::sort(points_left.begin(), points_left.end(), edge_compare3);//��point_right�еĵ㰴v��С��������
			   

			   //�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[0][j].point.x() >= rec.minU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_down.size(); ++k)
						   {
							   if ((points[0][j].point.x() == uvedge_down[k].point.x()) && (points[0][j].point.y() == uvedge_down[k].point.y()))
							   {
								   points[0][j].Id = uvedge_down[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�ж��ұ߽��ϵĵ�,��uvEdge�е������Ÿ����ұ߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[i][0].point.y() <= rec.maxV)//ȷ�����ڷֿ���//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_right.size(); ++k)
						   {
							   if ((points[i][0].point.x() == uvedge_right[k].point.x()) && (points[i][0].point.y() == uvedge_right[k].point.y()))
							   {
								   points[i][0].Id = uvedge_right[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�Ը÷ֿ�ı߽�����������***************************************************************************************



			   //����Ĺ�ϵ��nId��¼����********************************************************************************************

			  // double up_value = 2, down_value = 2;//��¼�߽����ڵ���߽��ľ���,2�ǳ�ʼֵ
			   //int up_index = 0, down_index = 0;//��¼�߽����ڵ�������ı߽��ı��
			   Member m;//���ڼ�¼����Ϣ
				std::vector<Member> edge_m;//��¼������߽紦�ĸ˼�
			 
			   //���ȴ���������Ͻǵ�ı߽磬��ʱ��Ϊ��߽磬�������Ա�֤�ǵ�����
			   //��߽紦��
				edge_m = Edge_Member(uvEdge, uvedge_left, points_left,3);//y��������
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //�ϱ߽紦��
			   edge_m = Edge_Member(uvEdge, uvedge_up, points_up, 2);//x��������
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //�ǵ㴦��
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
			 
			   //���±߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int j = 1; j < 999; ++j)
			   {				   
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if ((points[0][j].point.x() > rec.minU) && (points[1][j].point.x()>rec.minU))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //���ұ߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {

					   if ((points[i][0].point.y() < rec.maxV) && (points[i][1].point.y()<rec.maxV))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //�����ڲ��㣬���������±߽�ĵ㣬�˴�i,j��ʼΪ1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() > rec.minU) && (points[i][j].point.y() < rec.maxV))//�жϵ��Ƿ�������ֿ���//*****************************��䶯
						   {
							   //Ѱ������߽�������ڲ���						
							   if (points[i][j + 1].point.x() <= rec.minU)//�ж�j��u��vС����߽磬j+1��u���ڵ�����߽�//*****************************��䶯
							   {
								   if ((points[i + 1][j].point.x() >rec.minU) && (points[i + 1][j].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//�жϵ���Ϸ����Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() >= rec.maxV)//Ѱ�����ϱ߽�������ڲ���//*****************************��䶯
							   {								   
								   if ((points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//�жϵ���󷽵��Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i][j+1]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else
							   {
								   for (size_t k = 0; k < points[i][j].next.size(); ++k)
								   {
									   m = Add_2Member(points[i][j], points[points[i][j].next[k].x][points[i][j].next[k].y]);
									   mDisplayVec.push_back(m);
								   }
							   }
						   }
					   }
				   }
			   }
			  
			   break;
	}
		//����Ĺ�ϵ��nId��¼����********************************************************************************************			   

	case 2://�����ϵ�����
	{
			   //�Ը÷ֿ�ı߽�����������**********************************************************************************
			   //���߽���uv����Ŷ������������UVForOut��

			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare2_x);//�Ա߽���С��������
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare2_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare4_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare4_y);

			   //�Էֿ��ڲ��ĵ����������
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() >rec.minU) && (points[i][j].point.y() > rec.minV))//*****************************��䶯
							   //������Ҫ��������ڵĵ�˳��洢	
						   {							  
							   if ((points[i + 1][j].point.x()>rec.minU) && (points[i + 1][j].point.y() >rec.minV) && (points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y()>rec.minV))
							   {
								   points[i][j].Id = nId;
								   addspace.setX(points[i][j].point.x());
								   addspace.setY(points[i][j].point.y());
								   addspace.setID(nId);
								   UVForOut.push_back(addspace);
								   ++nId;
							   }
							   //�ж��Ƿ�Ϊ�����������ϱ߽�ĵ�********************************��Ķ�
							   //��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��,���Ǽ����ֻ����ͨ�����
							   //if (points[i + 1][j].point.y() <= rec.minV)
							   //{
								  // points_down.push_back(&points[i][j]);
							   //}
							   ////�ж��Ƿ�Ϊ�����������ұ߽�ĵ�********************************��Ķ�
							   ////��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��
							   //if (points[i][j + 1].point.x() <= rec.minU)
							   //{
								  // points_left.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   //point_corner = &points[i][j];//����ǵ�
								   point_corner.push_back(&points[i][j]);//����ǵ�
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_left[uvedge_left.size()-2]);//3
								   point_corner.push_back(&uvedge_left.back());//4
								   point_corner.push_back(&uvedge_down[uvedge_down.size()-2]);//5
								   point_corner.push_back(&uvedge_down.back());//6
							   }
							   else if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() > rec.minU))
							   {
								   points_down.push_back(&points[i][j]);//�±߽總�����ڲ���
							   }
							   else if ((points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   points_left.push_back(&points[i][j]);//��߽總�����ڲ���
							   }
						   }
					   }
				   }
			   }
			   //�Բ�����߽��Ե��ڲ�������������£�
			   std::sort(points_down.begin(), points_down.end(), edge_compare2);//��point_down�еĵ㰴u�Ӵ�С����
			   std::sort(points_left.begin(), points_left.end(), edge_compare4);//��point_left�еĵ㰴v�Ӵ�������
			  

			   //�ж��ϱ߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[0][j].point.x() >= rec.minU)//ȷ�����ڷֿ��ڣ�Ҫ�ߵ㣬������=//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_up.size(); ++k)
						   {
							   if ((points[0][j].point.x() == uvedge_up[k].point.x()) && (points[0][j].point.y() == uvedge_up[k].point.y()))
							   {
								   points[0][j].Id = uvedge_up[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�ж��ұ߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[i][0].point.y() >= rec.minV)//ȷ�����ڷֿ��ڣ�Ҫ�ߵ㣬������=//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_right.size(); ++k)
						   {
							   if ((points[i][0].point.x() == uvedge_right[k].point.x()) && (points[i][0].point.y() == uvedge_right[k].point.y()))
							   {
								   points[i][0].Id = uvedge_right[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�Ը÷ֿ�ı߽�����������***************************************************************************************



			   //����Ĺ�ϵ��nId��¼����********************************************************************************************

			   //double up_value = 2, down_value = 2;//��¼�߽����ڵ���߽��ľ���,2�ǳ�ʼֵ
			   //int up_index = 0, down_index = 0;//��¼�߽����ڵ�������ı߽��ı��
			   Member m;//���ڼ�¼����Ϣ
			   std::vector<Member> edge_m;//��¼������߽紦�ĸ˼�

			   //���ȴ���������½ǵ�ı߽磬��ʱ��Ϊ�±߽磬�������Ա�֤�ǵ�����
			   //�±߽紦��
			   edge_m = Edge_Member(uvEdge, uvedge_down, points_down, 2);//x��������
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //��߽紦��
			   edge_m = Edge_Member(uvEdge, uvedge_left, points_left, 4);//y��������
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //�ǵ㴦��
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
			   

			   //���±߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int j = 1; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if ((points[0][j].point.x() > rec.minU) && (points[1][j].point.x()>rec.minU))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //���ұ߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if ((points[i][0].point.y() > rec.minV) && (points[i][1].point.y()> rec.minV))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //�����ڲ��㣬���������±߽�ĵ㣬�˴�i,j��ʼΪ1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() > rec.minU) && (points[i][j].point.y() > rec.minV))//�жϵ��Ƿ�������ֿ���//*****************************��䶯
						   {
							   //Ѱ������߽�������ڲ���						
							   if (points[i][j + 1].point.x() <= rec.minU)//�ж�j��u��vС����߽磬j+1��u���ڵ�����߽�//*****************************��䶯
							   {								  
								   if ((points[i + 1][j].point.x() >rec.minU) && (points[i + 1][j].point.y()>rec.minV) && (points[i][j].next.size() != 0))//�жϵ���·����Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() <= rec.minV)//Ѱ�����ϱ߽�������ڲ���//*****************************��䶯
							   {								   
								   if ((points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y()>rec.minV) && (points[i][j].next.size() != 0))//�жϵ���󷽵��Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i][j+1]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else
							   {
								   for (size_t k = 0; k < points[i][j].next.size(); ++k)
								   {
									   m = Add_2Member(points[i][j], points[points[i][j].next[k].x][points[i][j].next[k].y]);
									   mDisplayVec.push_back(m);
								   }
							   }
						   }
					   }
				   }
			   }

			   break;
	}

	case 3://�����ϵ�����
	{
			   //�Ը÷ֿ�ı߽�����������**********************************************************************************
			   //���߽���uv����Ŷ������������UVForOut��
			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare1_x);//�Ա߽���С��������
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare1_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare4_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare4_y);

			   //�Էֿ��ڲ��ĵ����������
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() <rec.maxU) && (points[i][j].point.y() > rec.minV))//*****************************��䶯
							   //������Ҫ��������ڵĵ�˳��洢	
						   {
							   //�ж��Ƿ�Ϊ������������߽�ĵ�**********************************��Ķ�
							   if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() <rec.maxU) && (points[i][j + 1].point.y() >rec.minV))
							   {
								   points[i][j].Id = nId;
								   addspace.setX(points[i][j].point.x());
								   addspace.setY(points[i][j].point.y());
								   addspace.setID(nId);
								   UVForOut.push_back(addspace);
								   ++nId;
							   }
							   //�ж��Ƿ�Ϊ�����������ϱ߽�ĵ�********************************��Ķ�
							   //��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��,���Ǽ����ֻ����ͨ�����
							   //if (points[i + 1][j].point.y() <= rec.minV)
							   //{
								  // points_down.push_back(&points[i][j]);
							   //}
							   ////�ж��Ƿ�Ϊ�����������ұ߽�ĵ�********************************��Ķ�
							   ////��[i][j]���ڵ����㶼�п����ǳ�Խ�߽��
							   //if (points[i][j + 1].point.x() >= rec.maxU)
							   //{
								  // points_right.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() >= rec.maxU))
							   {
								  // point_corner = &points[i][j];//����ǵ�
								   point_corner.push_back(&points[i][j]);//����ǵ�
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_right[uvedge_right.size()-2]);//3
								   point_corner.push_back(&uvedge_right.back());//4
								   point_corner.push_back(&uvedge_down[uvedge_down.size()-2]);//5
								   point_corner.push_back(&uvedge_down.back());//6
							   }
							   else if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() < rec.maxU))
							   {
								   points_down.push_back(&points[i][j]);//�±߽總�����ڲ���
							   }
							   else if ((points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() >= rec.maxU))
							   {
								   points_right.push_back(&points[i][j]);//�ұ߽總�����ڲ���
							   }
						   }
					   }
				   }
			   }
			   //�Բ�����߽��Ե��ڲ�������������£�
			   std::sort(points_down.begin(), points_down.end(), edge_compare1);//��point_down�еĵ㰴u��С��������
			   std::sort(points_right.begin(), points_right.end(), edge_compare4);//��point_right�еĵ㰴v�Ӵ�С����
			  

			   //�ж��ϱ߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[0][j].point.x() <= rec.maxU)//ȷ�����ڷֿ��ڣ�Ҫ�ߵ㣬������=//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_up.size(); ++k)
						   {
							   if ((points[0][j].point.x() == uvedge_up[k].point.x()) && (points[0][j].point.y() == uvedge_up[k].point.y()))
							   {
								   points[0][j].Id = uvedge_up[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�
			   for (int i = 0; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if (points[i][0].point.y() >= rec.minV)//ȷ�����ڷֿ��ڣ�Ҫ�ߵ㣬������=//*****************************��䶯
					   {
						   for (size_t k = 0; k < uvedge_left.size(); ++k)
						   {
							   if ((points[i][0].point.x() == uvedge_left[k].point.x()) && (points[i][0].point.y() == uvedge_left[k].point.y()))
							   {
								   points[i][0].Id = uvedge_left[k].Id;
								   break;
							   }
						   }
					   }
				   }
			   }
			   //�Ը÷ֿ�ı߽�����������***************************************************************************************
			   //����Ĺ�ϵ��nId��¼����********************************************************************************************

			   //double up_value = 2, down_value = 2;//��¼�߽����ڵ���߽��ľ���,2�ǳ�ʼֵ
			  // int up_index = 0, down_index = 0;//��¼�߽����ڵ�������ı߽��ı��
			   Member m;//���ڼ�¼����Ϣ
			   std::vector<Member> edge_m;//��¼������߽紦�ĸ˼�

			   //���ȴ���������½ǵ�ı߽磬��ʱ��Ϊ�ұ߽磬�������Ա�֤�ǵ�����
			   //�ұ߽紦��
			   edge_m = Edge_Member(uvEdge, uvedge_right, points_right, 4);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //�±߽紦��
			   edge_m = Edge_Member(uvEdge, uvedge_down, points_down, 1);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//����mDispalyVec��
			   //�ǵ㴦��
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());

			   //���±߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int j = 1; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if ((points[0][j].point.x() < rec.maxU) && (points[1][j].point.x()<rec.maxU))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //���ұ߽����ڲ�����������,�˴�i��1��ʼ����ȥ���½ǵ�
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
				   {
					   if ((points[i][0].point.y() > rec.minV) && (points[i][1].point.y()>rec.minV))//ȷ�����ڷֿ��ڣ���Ҫ�ߵ㣬����û��=//*****************************��䶯
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//����
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //�����ڲ��㣬���������±߽�ĵ㣬�˴�i,j��ʼΪ1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//�жϵ�points[i][j]�Ƿ񱻸���x��yֵ
					   {
						   if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() > rec.minV))//�жϵ��Ƿ�������ֿ���//*****************************��䶯
						   {
							   //Ѱ������߽�������ڲ���						
							   if (points[i][j + 1].point.x() >= rec.maxU)//�ж�j��u��vС����߽磬j+1��u���ڵ�����߽�//*****************************��䶯
							   {								  
								   if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y()>rec.minV) && (points[i][j].next.size() != 0))//�жϵ���·����Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() <= rec.minV)//Ѱ�����ϱ߽�������ڲ���//*****************************��䶯
							   {								   
								   if ((points[i][j + 1].point.x() < rec.maxU) && (points[i][j + 1].point.y()>rec.minV) && (points[i][j].next.size() != 0))//�жϵ���ҷ����Ƿ��ڷֿ���
								   {
									   m = Add_2Member(points[i][j], points[i ][j+1]);//���������γɸ˼�
									   mDisplayVec.push_back(m);//�˼�����mDisplayVec						
								   }
							   }
							   else
							   {
								   for (size_t k = 0; k < points[i][j].next.size(); ++k)
								   {
									   m = Add_2Member(points[i][j], points[points[i][j].next[k].x][points[i][j].next[k].y]);
									   mDisplayVec.push_back(m);
								   }
							   }
						   }
					   }
				   }
			   }

			   break;
	}
	default:
		break;
	}

	for (int i = 0; i < 1000; ++i)
	{
		delete[] points[i];
		/*delete[] nodeForDisplay[i];
		delete[] nodeForSave[i];*/
	}
	delete[] points;
	/*delete[] nodeForDisplay;
	delete[] nodeForSave;*/


	return mDisplayVec;
}


//��Ӹ˼���point_in��һ�㣨�ڲ��㣩��point_edge�ڶ��㣨������߽�㣩
Member NurbsSurface_qj::Add_Member(PointStruct point_in, PointStruct point_edge)
{
	Member m;
	Node *tempNode1;
	Node *tempNode2;
	m.setInode(point_in.Id);//д���ϵ�
	m.setJnode(point_edge.Id);
	tempNode1 = GetSurfacePointatUandV(point_in.point.x(), point_in.point.y());
	tempNode2 = GetSurfacePointatUandV(point_edge.point.x(),point_edge.point.y());
	m.setCutLength((*tempNode1 - *tempNode2).GetLength());
	m.set_ID(mId);
	mId++;
	delete tempNode1;
	delete tempNode2;
	return m;
}

//��Ӹ˼���point_in��һ�㣨�ڲ��㣩��point_edge�ڶ��㣨����߽磩
Member NurbsSurface_qj::Add_2Member(PointStruct point_in, PointStruct point_edge)
{
	Member m;
	Node *tempNode1;
	Node *tempNode2;
	m.setInode(point_in.Id);//д���ϵ�
	m.setJnode(point_edge.Id);
	tempNode1 = GetSurfacePointatUandV(point_in.point.x(), point_in.point.y());
	tempNode2 = GetSurfacePointatUandV(point_edge.point.x(), point_edge.point.y());
	m.setCutLength((*tempNode1 - *tempNode2).GetLength());
	m.set_ID(mId);
	mId++;
	delete tempNode1;
	delete tempNode2;
	return m;
}

//�жϿ����߽�ĵ���ڱ߽������䣬����up_index,down_index********(�߽�㣬����ֵ����ֵ��xy�����ϵ㣬�µ㣩
void NurbsSurface_qj::BetweenSearch(std::vector<PointStruct> uv_edge, double key1, double key2, int xy, int &up_index, int &down_index)
{
	double up_value = 2, down_value = 2;
	up_index = -1;
	down_index = -1;
	switch (xy)
	{
	case 0://��ʾѰ��y���������
	{
			   for (int k = 0; k != uv_edge.size(); ++k)
			   {
				   if (uv_edge[k].point.x() == key1)//ɸѡ��ĳ�߽��ϵĵ�
				   {
					   if (uv_edge[k].point.y() - key2 >= 0)//�˵�λ���Ϸ�
					   {
						   if (uv_edge[k].point.y() - key2 <= up_value)//���ϸ���ֵС
						   {
							   up_value = uv_edge[k].point.y() - key2;
							   up_index = k;//��¼���std::vector���
						   }
					   }
					   else//�˵�λ���·�
					   {
						   if (key2 - uv_edge[k].point.y() <= down_value)//���ϸ���ֵС
						   {
							   down_value =  key2 - uv_edge[k].point.y();
							   down_index = k;// ��¼���std::vector���
						   }
					   }
				   }
			   }
			   break;
	}
	case 1:
	{
			  for (int k = 0; k != uv_edge.size(); ++k)
			  {
				  if (uv_edge[k].point.y() == key1)//ɸѡ��ĳ�߽��ϵĵ�
				  {
					  if (uv_edge[k].point.x() - key2 >= 0)//�˵�λ���Ϸ�
					  {
						  if (uv_edge[k].point.x() - key2 <= up_value)//���ϸ���ֵС
						  {
							  up_value = uv_edge[k].point.x() - key2;
							  up_index = k;//��¼���std::vector���
						  }
					  }
					  else//�˵�λ���·�
					  {
						  if (key2-uv_edge[k].point.x()  <= down_value)//���ϸ���ֵС
						  {
							  down_value =  key2- uv_edge[k].point.x() ;
							  down_index = k;// ��¼���std::vector���
						  }
					  }
				  }
			  }
			  break;
	}
		default:
		break;
	}
}

//���طֿ���������ұ߽������������Ա߽����������ţ����߽����뵽������У�&�����޸�uvEdge��Id
void NurbsSurface_qj::Edge_udlr(std::vector<PointStruct> &uvEdge, rectangleRegion rec, std::vector<PointStruct> &uvedge_up, std::vector<PointStruct> &uvedge_down, std::vector<PointStruct> &uvedge_left, std::vector<PointStruct> &uvedge_right)
{
	//std::vector<PointStruct> uv_edge;
	//int indexld = -1, indexrd = -1, indexlu = -1, indexru = -1;//��Ƿֲ��ĸ��ǵ��Ƿ��Ѿ����������ţ�����leftdown��rightdown��leftup��rightup		
	Node addspace;//��������NodeForOut�����Ĵ�С
	//PointStruct indexAdd;//��ʱ���������ӽǵ���������uvEdge�Ĵ�С��������ű��������еķֲ��ǵ�
	//indexAdd.Id = -1;//Id��ʼ��
	for (size_t i = 0; i != uvEdge.size(); ++i)//ȡ��uvEdge[i]
	{
		//1�ж�uvEdge[i]��vֵ�Ƿ����ϱ߽����,���Ƿ���ڸò���,�����ҽǵ�,��������ǵ�
		if ((uvEdge[i].point.y() == rec.maxV) && (uvEdge[i].point.x()>rec.minU) && (uvEdge[i].point.x()<=rec.maxU))
		{
			if (uvEdge[i].Id == -1)//�жϸõ��Ƿ��Ѿ������˱�ţ���û��ʱ�Ž������²���
			{
				uvEdge[i].Id = nId;//���߽���id��ֵΪ����nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//��uvEdge��ֵ����addspace
				UVForOut.push_back(addspace);//���߽����뵽������
				++nId;//������������
			}
			uvedge_up.push_back(uvEdge[i]);//���÷ֿ��ϱ߽�ĵ����uvedge_up
			continue;
		}

		//2�ж�uvEdge[i]��vֵ�Ƿ����±߽����,���Ƿ���ڸò��֣�������ǵ㣬�������ҽǵ�
		if ((uvEdge[i].point.y() == rec.minV) && (uvEdge[i].point.x()>=rec.minU) && (uvEdge[i].point.x()<rec.maxU))
		{
			if (uvEdge[i].Id == -1)//�жϸõ��Ƿ��Ѿ������˱�ţ���û��ʱ�Ž������²���
			{
				uvEdge[i].Id = nId;//���߽���id��ֵΪ����nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//��uvEdge��ֵ����addspace
				UVForOut.push_back(addspace);//���߽����뵽������
				++nId;//������������
			}
			uvedge_down.push_back(uvEdge[i]);//���÷ֿ��ϱ߽�ĵ����uvedge_down
			continue;
		}

		//3�ж�uvEdge[i]��uֵ�Ƿ�����߽����,���Ƿ���ڸò��֣������Ͻǵ�,�������½ǵ�
		if ((uvEdge[i].point.x() == rec.minU) && (uvEdge[i].point.y()>rec.minV) && (uvEdge[i].point.y()<=rec.maxV))
		{
			if (uvEdge[i].Id == -1)//�жϸõ��Ƿ��Ѿ������˱�ţ���û��ʱ�Ž������²���
			{
				uvEdge[i].Id = nId;//���߽���id��ֵΪ����nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//��uvEdge��ֵ����addspace
				UVForOut.push_back(addspace);//���߽����뵽������
				++nId;//������������
			}
			uvedge_left.push_back(uvEdge[i]);//���÷ֿ��ϱ߽�ĵ����uvedge_left
			continue;
		}

		//4�ж�uvEdge[i]��uֵ�Ƿ����ұ߽����,���Ƿ���ڸò��֣������½ǵ㣬�������Ͻǵ�
		if ((uvEdge[i].point.x() == rec.maxU) && (uvEdge[i].point.y()>=rec.minV) && (uvEdge[i].point.y() < rec.maxV))
		{
			if (uvEdge[i].Id == -1)//�жϸõ��Ƿ��Ѿ������˱�ţ���û��ʱ�Ž������²���
			{
				uvEdge[i].Id = nId;//���߽���id��ֵΪ����nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//��uvEdge��ֵ����addspace
				UVForOut.push_back(addspace);//���߽����뵽������
				++nId;//������������
			}
			uvedge_right.push_back(uvEdge[i]);//���÷ֿ��ϱ߽�ĵ����uvedge_right
		}
	}
	
	
}


//�߽�˼�����(ȫ���߽�������������ӱ߽�㡾���ڲ���λ�����߽���м�����߽�ܽ�ʱ����ĳһ�߽�ĵ㣬�����߽���ڲ���,uv����
//ֻ������һ���߽������������������ڲ��㣬�����������ڲ�����������������������ܣ�
std::vector<Member> NurbsSurface_qj::Edge_Member(std::vector<PointStruct> &uvEdge, std::vector<PointStruct> uvedge, std::vector<PointStruct*> pointsedge, int xy, double length_differ)
{
	Member m;//��������edge_m
	std::vector<Member> edge_m;//�洢�߽�˼�
	Member point_m;//��¼�ڲ�������±߽���������
	std::vector<Member> edge_point_m;//�洢point_m
	int up_index = -1;//���µ��ID	
	int max_or_min = 0;//���point_corner�Ǽ����Ǽ�Сֵ
	PointStruct addpoint;
	Node addspace;
	
	switch (xy)
	{
	case 1://u���������
	{
			  
			   for (size_t i = 0; i < pointsedge.size(); i++)//Ѱ���ڲ�������±߽��
			   {	
				   up_index = -1;				   
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.x() - pointsedge[i]->point.x() > 0)//u���ڴ˵�ĵ�һ���߽��
					   {
						   up_index = j;//���ұ߽����±괫��
						   break;
					   }
					}

				   if (up_index >= 1 && up_index <= uvedge.size())//�ڲ��������ұ߽����֮�������������ж�����֮��ľ��룬�˴��趨Ϊ1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());					   
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ҵ�����
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index-1].point.x(), uvedge[up_index-1].point.y());//�������
					   tempNode4 = GetSurfacePointatUandV(pointsedge[i]->point.x(), uvedge[up_index].point.y());//ͶӰ����
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//���ҵ�˼�������С��2m��Ĭ�ϣ�
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x()) && (pointsedge[i - 1]->point.x() < uvedge[up_index].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id) )//�ж��ϸ��ڲ���������Ƿ��ǡ�index-1��index����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i-1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�����˼�������С��2m��Ĭ�ϣ�
					   {
						   pointsedge[i]->Id = uvedge[up_index-1].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index-1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index-1].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.x() > uvedge[up_index - 2].point.x()) && (pointsedge[i - 1]->point.x() < uvedge[up_index-1].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//�ж��ϸ��ڲ���������Ƿ�Ϊ��index-2��index-1����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i-1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //ͶӰ�����С��2m��Ĭ�ϣ�����u�����м�
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) && 
						   (fabs(pointsedge[i]->point.x() - (uvedge[up_index - 1].point.x() + uvedge[up_index].point.x()) / 2)<(uvedge[up_index ].point.x()- uvedge[up_index - 1].point.x()) / 5))
					   {//�˴������������ڲ�������
						   pointsedge[i]->Id = nId;//���õ�����Ϊ�߽��
						   /*pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   addpoint.Id = nId;
						   addpoint.point.y() = uvedge[up_index].point.y();
						   addpoint.point.x() = pointsedge[i]->point.x();
						   uvEdge.push_back(addpoint);//��ӵ�ȫ�ֱ߽�����
						   //uvedge.insert(uvedge.begin() + up_index, addpoint);//��ӵ��ֲ��߽�
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //һ�����
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index-1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;					   
				   }
				   else if (up_index == 0)//�ڲ���ֻ���ұ߽����֮����
				   {
					   Node *tempNode1;					   
					   Node *tempNode2;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());					   
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ҵ�����					   
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�µ�˼�������С��1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i+1]->point.x()< uvedge[up_index+1].point.x())//�������һ���ڲ��������Ϊ��0,1����
						   {
							   pointsedge[i]->next.clear();//���õ���ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
						   delete tempNode1;						  
						   delete tempNode2;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //����1ʱ
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(-1);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;

					   delete tempNode1;					  
					   delete tempNode2;					 
				   }
				   else if (up_index == -1 )//�ڲ�������߽����֮����
				   {
					   Node *tempNode1;
					   Node *tempNode3;	
					   up_index = uvedge.size() - 1;//ʹup_indexָ��߽������±�
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�������					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�ϵ�˼�������С��1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�������һ���ڲ�������䡾index-1��index����==��ʾ�ϸ��㱻�鲢
						   {
							   pointsedge[i-1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
						   delete tempNode1;
						   delete tempNode3;						  
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }					  
					   //����1ʱ
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;					   
				   }
			   }
			   //�ж��Ƿ��������ڲ������ͬ�������߽������
			   for (size_t i = 0; i < edge_point_m.size()-1; i++)
			   {
				   for (size_t j = i+1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//�������ͬ����ĳ��ı�������
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }

			   //�����ڲ�����߽�������ĸ˼�
			   for (size_t i = 0; i < edge_point_m.size(); i++)
			   {
				   if (edge_point_m[i].getInode() != -1)
				   {	
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getInode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
				   if (edge_point_m[i].getJnode() != -1)
				   {
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getJnode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
			   }		

			  break;
	}

	case 2://u������ҵ���
	{

			   for (size_t i = 0; i < pointsedge.size(); i++)//Ѱ���ڲ�������±߽��
			   {
				   up_index = -1;
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.x() - pointsedge[i]->point.x() < 0)//u���ڴ˵�ĵ�һ���߽�*************************
					   {
						   up_index = j;//���ұ߽����±괫��
						   break;
					   }
				   }

				   if (up_index >= 1 && up_index <= uvedge.size())//�ڲ��������ұ߽����֮�������������ж�����֮��ľ��룬�˴��趨Ϊ1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ҵ�����
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//�������
					   tempNode4 = GetSurfacePointatUandV(pointsedge[i]->point.x(), uvedge[up_index].point.y());//ͶӰ����
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//���ҵ�˼�������С��2m��Ĭ�ϣ�
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.x() < uvedge[up_index - 1].point.x()) && (pointsedge[i - 1]->point.x() > uvedge[up_index].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�ж��ϸ��ڲ���������Ƿ��ǡ�index-1��index����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�����˼�������С��2m��Ĭ�ϣ�
					   {
						   pointsedge[i]->Id = uvedge[up_index - 1].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index-1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index-1].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.x()<uvedge[up_index - 2].point.x()) && (pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//�ж��ϸ��ڲ���������Ƿ�Ϊ��index-2��index-1����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //ͶӰ�����С��2m��Ĭ�ϣ�����u�����м�
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						   (fabs(pointsedge[i]->point.x() - (uvedge[up_index - 1].point.x() + uvedge[up_index].point.x()) / 2)<(uvedge[up_index].point.x() - uvedge[up_index - 1].point.x()) / 5))
					   {//�˴������������ڲ�������
						   pointsedge[i]->Id = nId;//���õ�����Ϊ�߽��
						   /*pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   addpoint.Id = nId;
						   addpoint.point.y() = uvedge[up_index].point.y();
						   addpoint.point.x() = pointsedge[i]->point.x();
						   uvEdge.push_back(addpoint);//��ӵ�ȫ�ֱ߽�����
						   //uvedge.insert(uvedge.begin() + up_index, addpoint);//��ӵ��ֲ��߽�
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //һ�����
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index - 1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;
				   }
				   else if (up_index == 0)//�ڲ���ֻ���ұ߽����֮����
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ҵ�����					   
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�µ�˼�������С��1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i + 1]->point.x()> uvedge[up_index + 1].point.x())//�������һ���ڲ��������Ϊ��0,1����
						   {
							   pointsedge[i]->next.clear();//���õ���ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
						   delete tempNode1;
						   delete tempNode2;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //����1ʱ
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(-1);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;

					   delete tempNode1;
					   delete tempNode2;
				   }
				   else if (up_index == -1)//�ڲ�������߽����֮����*******�����
				   {
					   Node *tempNode1;
					   Node *tempNode3;
					   up_index = uvedge.size() - 1;//ʹup_indexָ��߽������±�
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�������					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�ϵ�˼�������С��1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.x()< uvedge[up_index - 1].point.x() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�������һ���ڲ�������䡾index-1��index����==��ʾ�ϸ��㱻�鲢
						   {//*****************�˴�����Ϊ��
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
						   delete tempNode1;
						   delete tempNode3;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //����1ʱ
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;
				   }
			   }
			   //�ж��Ƿ��������ڲ������ͬ�������߽������
			   for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			   {
				   for (size_t j = i + 1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//�������ͬ����ĳ��ı�������
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }

			   //�����ڲ�����߽�������ĸ˼�
			   for (size_t i = 0; i < edge_point_m.size(); i++)
			   {
				   if (edge_point_m[i].getInode() != -1)
				   {
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getInode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
				   if (edge_point_m[i].getJnode() != -1)
				   {
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getJnode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
			   }

			   break;
	}

	case 3://v�����������
	{
			  
			  for (size_t i = 0; i < pointsedge.size(); i++)//Ѱ���ڲ�������±߽��
			  {
				  up_index = -1;				  
				  for (size_t j = 0; j < uvedge.size(); j++)
				  {
					  if (uvedge[j].point.y() - pointsedge[i]->point.y() > 0)//v���ڴ˵�ĵ�һ���߽��
					  {
						  up_index = j;//���ϱ߽����±괫��
						  break;
					  }					 
				  }

				  if (up_index >= 1 && up_index <= uvedge.size())//�ж��ڲ����Ƿ����ϡ��±߽����֮�������������ж�����֮��ľ��룬�˴��趨Ϊ1m
				  {
					  Node *tempNode1;
					  Node *tempNode2;
					  Node *tempNode3;
					  Node *tempNode4;
					  tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					  tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ϵ�����
					  tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//�µ�����
					  tempNode4 = GetSurfacePointatUandV(uvedge[up_index].point.x(), pointsedge[i]->point.y());//ͶӰ����
					  
					  if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�ϵ�˼�������С��2m
					  {
						  pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						 /* pointsedge[i]->point.x() = uvedge[up_index].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						  if (i != 0)//��֤��Ϊ��һ���ڲ���
						  {
							  if ((pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y()) && (pointsedge[i - 1]->point.y() < uvedge[up_index].point.y())
								  || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id) )//�ж��ϸ��ڲ���������Ƿ��ǡ�index-1��index����==��ʾ�ϸ��㱻�鲢
							  {
								  pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							  }
						  }
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//�õ��Ѵ����꣬�����¸�i
					  }
					  if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�µ�˼�������С��2m
					  {
						  pointsedge[i]->Id = uvedge[up_index - 1].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						  /*pointsedge[i]->point.x() = uvedge[up_index - 1].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index - 1].point.y();*/
						  if (i != 0)//��֤��Ϊ��һ���ڲ���
						  {
							  if ((pointsedge[i - 1]->point.y() > uvedge[up_index - 2].point.y()) && (pointsedge[i - 1]->point.y() < uvedge[up_index-1].point.y())
								  || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id) )//�ж��ϸ��ڲ���������Ƿ��ǡ�index-2��index-1����==��ʾ�ϸ��㱻�鲢
							  {
								  pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							  }
						  }
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//�õ��Ѵ����꣬�����¸�i
					  }
					  if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						  (fabs(pointsedge[i]->point.y() - (uvedge[up_index - 1].point.y() + uvedge[up_index].point.y()) / 2)<(uvedge[up_index].point.y() - uvedge[up_index - 1].point.y()) / 5))
					  {//�˴������������ڲ�������
						  pointsedge[i]->Id = nId;//���õ�����Ϊ�߽��
						 /* pointsedge[i]->point.x() = uvedge[up_index].point.x();*/
						  addpoint.Id = nId;
						  addpoint.point.x() = uvedge[up_index].point.x();
						  addpoint.point.y() = pointsedge[i]->point.y();
						  uvEdge.push_back(*pointsedge[i]);//��ӵ�ȫ�ֱ߽�����
						  uvedge.insert(uvedge.begin() + up_index, addpoint);//��ӵ��ֲ��߽�
						  addspace.setX(addpoint.point.x());
						  addspace.setY(addpoint.point.y());
						  addspace.setID(nId);
						  UVForOut.push_back(addspace);
						  nId++;
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//�õ��Ѵ����꣬�����¸�i
					  }
					  //һ�����
					  pointsedge[i]->Id = nId;//�Ըõ����������
					  point_m.set_ID(nId);
					  point_m.setInode(uvedge[up_index - 1].Id);
					  point_m.setJnode(uvedge[up_index].Id);					 
					  edge_point_m.push_back(point_m);//������������ϵ��¼����
					  addspace.setX(pointsedge[i]->point.x());
					  addspace.setY(pointsedge[i]->point.y());
					  addspace.setID(nId);
					  UVForOut.push_back(addspace);
					  nId++;
					  if (i != 0)//��֤��Ϊ��һ���ڲ���
					  {
						  if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						  {
							  pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						  }
					  }
					  delete tempNode1;
					  delete tempNode2;
					  delete tempNode3;
					  delete tempNode4;
				  }
				  else if (up_index == 0)//�ڲ���ֻ���ϱ߽����֮����
				  {					  
						Node *tempNode1;
						Node *tempNode2;
						tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
						tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�µ�����					   
						if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�µ�˼�������С��1m
						{
							pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
							/*pointsedge[i]->point.x() = uvedge[up_index].point.x();
							pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

							if (pointsedge[i + 1]->point.y() < uvedge[up_index + 1].point.y())//�����һ���ڲ��������Ϊ��0,1����
							{
								pointsedge[i]->next.clear();//���õ���ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							}
							delete tempNode1;
							delete tempNode2;
							continue;//�õ��Ѵ����꣬�����¸�i
						}
						//������1ʱ
						pointsedge[i]->Id = nId;//�Ըõ����������
						point_m.set_ID(nId);
						point_m.setInode(-1);
						point_m.setJnode(uvedge[up_index].Id);
						edge_point_m.push_back(point_m);//������������ϵ��¼����
						addspace.setX(pointsedge[i]->point.x());
						addspace.setY(pointsedge[i]->point.y());
						addspace.setID(nId);
						UVForOut.push_back(addspace);
						nId++;
						delete tempNode1;
						delete tempNode2;					 
				  }
				  else if (up_index == -1)//�ڲ������±߽����֮����
				  {
					  Node *tempNode1;
					  Node *tempNode3;
					  up_index = uvedge.size() - 1;//ʹup_indexָ��߽������±�
					  tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					  tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ϵ�����					   
					  if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�ϵ�˼�������С��1m
					  {
						  pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						  /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						  if (pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�����һ���ڲ��������Ϊ��index-1,index��������==��ʾ�ϸ��㱻�鲢
						  {
							  pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						  }
						  delete tempNode1;
						  delete tempNode3;
						  continue;//�õ��Ѵ����꣬�����¸�i
					  }
					  //����1ʱ
					  pointsedge[i]->Id = nId;//�Ըõ����������
					  point_m.set_ID(nId);
					  point_m.setInode(uvedge[up_index].Id);
					  point_m.setJnode(-1);
					  edge_point_m.push_back(point_m);//������������ϵ��¼����
					  addspace.setX(pointsedge[i]->point.x());
					  addspace.setY(pointsedge[i]->point.y());
					  addspace.setID(nId);
					  UVForOut.push_back(addspace);
					  nId++;
					  if (i != 0)//��֤��Ϊ��һ���ڲ���
					  {
						  if (pointsedge[i - 1]->Id == uvedge[up_index ].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						  {
							  pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						  }
					  }
					  delete tempNode1;
					  delete tempNode3;
				  }
			  }
			  //�ж��Ƿ��������ڲ������ͬ�������߽������
			  for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			  {
				  for (size_t j = i + 1; j < edge_point_m.size(); j++)
				  {
					  if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					  {//�������ͬ����ĳ��ı�������
						  edge_point_m[i].setJnode(-1);
						  edge_point_m[j].setInode(-1);
					  }
				  }
			  }
			  //�����˼�
			  for (size_t i = 0; i < edge_point_m.size(); i++)
			  {
				  if (edge_point_m[i].getInode() != -1)
				  {
					  m.setInode(edge_point_m[i].getID());
					  m.setJnode(edge_point_m[i].getInode());
					  m.set_ID(mId);
					  edge_m.push_back(m);
					  mId++;
				  }
				  if (edge_point_m[i].getJnode() != -1)
				  {
					  m.setInode(edge_point_m[i].getID());
					  m.setJnode(edge_point_m[i].getJnode());
					  m.set_ID(mId);
					  edge_m.push_back(m);
					  mId++;
				  }
			  }
			  break;
	}

	case 4://v�����������
	{

			   for (size_t i = 0; i < pointsedge.size(); i++)//Ѱ���ڲ�������±߽��
			   {
				   up_index = -1;
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.y() - pointsedge[i]->point.y() < 0)//v���ڴ˵�ĵ�һ���߽��
					   {
						   up_index = j;//���ϱ߽����±괫��
						   break;
					   }
				   }

				   if (up_index >= 1 && up_index <= uvedge.size())//�ж��ڲ����Ƿ����ϡ��±߽����֮�������������ж�����֮��ľ��룬�˴��趨Ϊ1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ϵ�����
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//�µ�����
					   tempNode4 = GetSurfacePointatUandV(uvedge[up_index].point.x(), pointsedge[i]->point.y());//ͶӰ����

					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�ϵ�˼�������С��2m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /* pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.y() < uvedge[up_index - 1].point.y()) && (pointsedge[i - 1]->point.y() > uvedge[up_index].point.y())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�ж��ϸ��ڲ���������Ƿ��ǡ�index-1��index����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�µ�˼�������С��2m
					   {
						   pointsedge[i]->Id = uvedge[up_index - 1].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index - 1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index - 1].point.y();*/
						   if (i != 0)//��֤��Ϊ��һ���ڲ���
						   {
							   if ((pointsedge[i - 1]->point.y() < uvedge[up_index - 2].point.y()) && (pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//�ж��ϸ��ڲ���������Ƿ��ǡ�index-2��index-1����==��ʾ�ϸ��㱻�鲢
							   {
								   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						   (fabs(pointsedge[i]->point.y() - (uvedge[up_index - 1].point.y() + uvedge[up_index].point.y()) / 2)<(uvedge[up_index].point.y() - uvedge[up_index - 1].point.y()) / 5))
					   {//�˴������������ڲ�������
						   pointsedge[i]->Id = nId;//���õ�����Ϊ�߽��
						   /* pointsedge[i]->point.x() = uvedge[up_index].point.x();*/
						   addpoint.Id = nId;
						   addpoint.point.x() = uvedge[up_index].point.x();
						   addpoint.point.y() = pointsedge[i]->point.y();
						   uvEdge.push_back(*pointsedge[i]);//��ӵ�ȫ�ֱ߽�����
						   uvedge.insert(uvedge.begin() + up_index, addpoint);//��ӵ��ֲ��߽�
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //һ�����
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index - 1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;
				   }
				   else if (up_index == 0)//�ڲ���ֻ���ϱ߽����֮����
				   {					   
						Node *tempNode1;
						Node *tempNode2;
						tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
						tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�µ�����					   
						if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//�µ�˼�������С��1m
						{
							pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
							/*pointsedge[i]->point.x() = uvedge[up_index].point.x();
							pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

							if (pointsedge[i + 1]->point.y() > uvedge[up_index + 1].point.y())//�����һ���ڲ��������Ϊ��0,1����
							{
								pointsedge[i]->next.clear();//���õ���ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
							}
							delete tempNode1;
							delete tempNode2;
							continue;//�õ��Ѵ����꣬�����¸�i
						}
						//������1ʱ
						pointsedge[i]->Id = nId;//�Ըõ����������
						point_m.set_ID(nId);
						point_m.setInode(-1);
						point_m.setJnode(uvedge[up_index].Id);
						edge_point_m.push_back(point_m);//������������ϵ��¼����
						addspace.setX(pointsedge[i]->point.x());
						addspace.setY(pointsedge[i]->point.y());
						addspace.setID(nId);
						UVForOut.push_back(addspace);
						nId++;
						delete tempNode1;
						delete tempNode2;
					  
				   }
				   else if (up_index == -1)//�ڲ������±߽����֮����
				   {
					   Node *tempNode1;
					   Node *tempNode3;
					   up_index = uvedge.size() - 1;//ʹup_indexָ��߽������±�
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//�ϵ�����					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//�ϵ�˼�������С��1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//���ڲ��㸳ֵΪ�߽�㣬�����õ�ͱ��鲢���߽���
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.y() < uvedge[up_index - 1].point.y() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//�����һ���ڲ��������Ϊ��index-1,index��������==��ʾ�ϸ��㱻�鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
						   delete tempNode1;
						   delete tempNode3;
						   continue;//�õ��Ѵ����꣬�����¸�i
					   }
					   //����1ʱ
					   pointsedge[i]->Id = nId;//�Ըõ����������
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//������������ϵ��¼����
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//��֤��Ϊ��һ���ڲ���
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//�ж��ϸ��ڲ����Ƿ񱻹鲢
						   {
							   pointsedge[i - 1]->next.clear();//����һ����ڲ����������ϵɾ�����������Ա�֤û���ظ��˼�
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;
				   }
			   }
			   //�ж��Ƿ��������ڲ������ͬ�������߽������
			   for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			   {
				   for (size_t j = i + 1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//�������ͬ����ĳ��ı�������
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }
			   //�����˼�
			   for (size_t i = 0; i < edge_point_m.size(); i++)
			   {
				   if (edge_point_m[i].getInode() != -1)
				   {
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getInode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
				   if (edge_point_m[i].getJnode() != -1)
				   {
					   m.setInode(edge_point_m[i].getID());
					   m.setJnode(edge_point_m[i].getJnode());
					   m.set_ID(mId);
					   edge_m.push_back(m);
					   mId++;
				   }
			   }
			   break;
	}

	default:
		break;
	}
	return edge_m;
}

//����߽�����
std::vector<Member> NurbsSurface_qj::EdgeTotal_Member(std::vector<PointStruct> uvEdge, std::vector<PLib::Point2Dd> extrmumVec, int num, double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV)
{
	std::vector<PointStruct> *uv_eg = new std::vector<PointStruct>[2 * num + 4];//�洢���߽��
	std::vector<Member> mVecTemp;//�����洢�����ɵı߽�member
	//ż����u����ĵ㣬������v����ĵ�
	for (size_t i = 0; i < uvEdge.size(); ++i)//����uvEdge
	{
		if (uvEdge[i].Id != -1)//ѡ��ʹ�õ��ı߽��
		{
			if (uvEdge[i].point.y() == minV)//�±߽�
			{
				uv_eg[0].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.x() == maxU)//�ұ߽�
			{
				uv_eg[1].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.y() == maxV)//�ϱ߽�
			{
				uv_eg[2].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.x() == minU)//��߽�
			{
				uv_eg[3].push_back(uvEdge[i]);
			}

			for (size_t j = 4; j < 2 * num + 4; j +=2)//�ڲ��߽�u����
			{
				if (uvEdge[i].point.y() == extrmumVec[(j - 4) / 2].y())
				{
					uv_eg[j].push_back(uvEdge[i]);
				}
			}

			for (size_t j = 5; j < 2 * num + 4; j+=2)//�ڲ��߽�v����
			{
				if (uvEdge[i].point.x() == extrmumVec[(j - 5) / 2].x())
				{
					uv_eg[j].push_back(uvEdge[i]);
				}
			}

		}
	}
	//������

	std::sort(uv_eg[0].begin(), uv_eg[0].end(), compare1_x);//�±߽�
	std::sort(uv_eg[1].begin(), uv_eg[1].end(), compare3_y);//�ұ߽�
	std::sort(uv_eg[2].begin(), uv_eg[2].end(), compare1_x);//�ϱ߽�
	std::sort(uv_eg[3].begin(), uv_eg[3].end(), compare3_y);//��߽�
	for (size_t j = 4; j < 2 * num + 4; j +=2)
	{
		std::sort(uv_eg[j].begin(), uv_eg[j].end(), compare1_x);//�ڲ�u�߽�
	}
	for (size_t j = 5; j < 2 * num + 4; j +=2)
	{
		std::sort(uv_eg[j].begin(), uv_eg[j].end(), compare3_y);//�ڲ�v�߽�
	}
	for (size_t i = 0; i < 2 * num + 4; ++i)//����uv_eg
	{
		for (size_t j = 0; j < uv_eg[i].size() - 1; ++j)//��ȡuv_eg[i]��ǰsize-1��Ԫ��
		{
			Member m;
			Node *tempNode1;
			Node *tempNode2;
			m.setInode(uv_eg[i][j].Id);//д���ϵ�
			m.setJnode(uv_eg[i][j + 1].Id);
			tempNode1 = GetSurfacePointatUandV(uv_eg[i][j].point.x(), uv_eg[i][j].point.y());
			tempNode2 = GetSurfacePointatUandV(uv_eg[i][j + 1].point.x(), uv_eg[i][j + 1].point.y());
			m.setCutLength((*tempNode1 - *tempNode2).GetLength());
			m.set_ID(mId);
			mId++;
			mVecTemp.push_back(m);
			delete tempNode1;
			delete tempNode2;
		}
	}
	//����߽���������

	for (size_t i = 4; i < 2 * num + 4; i +=2)//u��������
	{
		for (size_t j = 1; j < 2 * num + 4; j +=2)//����v�����uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());//u�����ϵ�����Ӧ����ά����
			tempNode2 = GetSurfacePointatUandV(uv_eg[j].begin()->point.x(), uv_eg[i].begin()->point.y());//�������v����߽��ϵ�ͶӰ�����ά����
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthU)//�жϵ�����Ǹ�v����߽����
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//Ѱ�����v�߽��ϸ��������±߽��
				{
					if (uv_eg[i].begin()->point.y()< iter->point.y())//�߽�㶼�Ǵ�С��������ģ�������һ��ʱ���ҵ������±߽��
					{

						Member m;
						m.setInode(uv_eg[i].begin()->Id);//д���ϵ�
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode(uv_eg[i].begin()->Id);//д���µ�
						m.setJnode((iter - 1)->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV((iter - 1)->point.x(), (iter - 1)->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);						
						break;
					}
				}
			}
			delete tempNode1;
			delete tempNode2;
		}
	}
	for (size_t i = 4; i < 2 * num + 4; i +=2)//u������յ�
	{
		for (size_t j = 1; j < 2 * num + 4; j +=2)//����v�����uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());//u�����ϵ�����Ӧ����ά����
			tempNode2 = GetSurfacePointatUandV(uv_eg[j].begin()->point.x(), (uv_eg[i].end() - 1)->point.y());//�������v����߽��ϵ�ͶӰ�����ά����
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthU)//�жϵ�����Ǹ�v����߽����
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//Ѱ�����v�߽��ϸ��������±߽��
				{
					if ((uv_eg[i].end() - 1)->point.y()< iter->point.y())//�߽�㶼�Ǵ�С��������ģ�������һ��ʱ���ҵ������±߽��
					{
						Member m;
						m.setInode((uv_eg[i].end() - 1)->Id);//д���ϵ�
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode((uv_eg[i].end() - 1)->Id);//д���µ�
						m.setJnode((iter - 1)->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV((iter - 1)->point.x(), (iter - 1)->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						break;						
					}
				}
			}
			delete tempNode1;
			delete tempNode2;
		}
	}

	for (size_t i = 5; i < 2 * num + 4; i +=2)//v��������
	{
		for (size_t j = 0; j < 2 * num + 4; j +=2)//����u�����uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());//v�����ϵ�����Ӧ����ά����
			tempNode2 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[j].begin()->point.y());//�������u����߽��ϵ�ͶӰ�����ά����
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthV)//�жϵ�����Ǹ�u����߽����
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//Ѱ�����u�߽��ϸ��������±߽��
				{
					if (uv_eg[i].begin()->point.x()< iter->point.x())//�߽�㶼�Ǵ�С��������ģ�������һ��ʱ���ҵ������±߽��
					{
						Member m;
						m.setInode(uv_eg[i].begin()->Id);//д���ϵ�
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode(uv_eg[i].begin()->Id);//д���µ�
						m.setJnode((iter - 1)->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV((iter - 1)->point.x(), (iter - 1)->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						break;
					}
				}
			}
			delete tempNode1;
			delete tempNode2;
		}
	}
	for (size_t i = 5; i < 2 * num + 4; i +=2)//v������յ�
	{
		for (size_t j = 0; j < 2 * num + 4; j +=2)//����u�����uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());//v�����ϵ�����Ӧ����ά����
			tempNode2 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), uv_eg[j].begin()->point.y());//�������u����߽��ϵ�ͶӰ�����ά����
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthV)//�жϵ�����Ǹ�u����߽����
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//Ѱ�����u�߽��ϸ��������±߽��
				{
					if ((uv_eg[i].end() - 1)->point.x()< iter->point.x())//�߽�㶼�Ǵ�С��������ģ�������һ��ʱ���ҵ������±߽��
					{
						Member m;
						m.setInode((uv_eg[i].end() - 1)->Id);//д���ϵ�
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode((uv_eg[i].end() - 1)->Id);//д���µ�
						m.setJnode((iter - 1)->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV((iter - 1)->point.x(), (iter - 1)->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						break;
					}
				}
			}
			delete tempNode1;
			delete tempNode2;
		}
	}
	delete[] uv_eg;
	return mVecTemp;
}

//����߽�Ľǵ�
std::vector<Member> NurbsSurface_qj::EdgeCorner_Member(std::vector<PointStruct*> &point_corner, double length_differ)
{
	
	Member m;//��������edge_m
	std::vector<Member> edge_m;//�洢�߽�˼�	
	Node addspace;
	int index34 = -1, index56 = -1;//��ǵ�0���3456�Ĺ�ϵ
	int max3456 = -1;
	

	if (point_corner[0]->point.x()>point_corner[6]->point.x())
	{
		if (point_corner[0]->point.x()>point_corner[5]->point.x())//0����56
		{
			index56 = 6;//0��6
		}
		else
		{
			index56 = 56;//0��56
		}
	}
	else
	{
		if (point_corner[0]->point.x()>point_corner[5]->point.x())//0����56
		{
			index56 = 56;//0��56
		}
		else
		{
			index56 = 6;//0��6
		}
	}

	if (point_corner[0]->point.y()>point_corner[4]->point.y())
	{
		if (point_corner[0]->point.y()>point_corner[3]->point.y())//0����56
		{
			index34= 4;//0��4
		}
		else
		{
			index34 = 34;//0��34
		}
	}
	else
	{
		if (point_corner[0]->point.y()>point_corner[3]->point.y())//
		{
			index34 = 34;//0��34
		}
		else
		{
			index34 = 4;//0��4
		}
	}

	Node *tempNode0;
	Node *tempNode3;
	Node *tempNode4;
	Node *tempNode5;
	Node *tempNode6;
	Node *tempNode1;
	Node *tempNode2;

	tempNode0 = GetSurfacePointatUandV(point_corner[0]->point.x(), point_corner[0]->point.y());//�ǵ�
	tempNode3 = GetSurfacePointatUandV(point_corner[3]->point.x(), point_corner[3]->point.y());//3����
	tempNode4 = GetSurfacePointatUandV(point_corner[4]->point.x(), point_corner[4]->point.y());//4����
	tempNode5 = GetSurfacePointatUandV(point_corner[5]->point.x(), point_corner[5]->point.y());//5����
	tempNode6 = GetSurfacePointatUandV(point_corner[6]->point.x(), point_corner[6]->point.y());//6����
	tempNode1 = ((*tempNode0 - *tempNode3).GetLength() < (*tempNode0 - *tempNode4).GetLength()) ? tempNode3 :tempNode4;
	tempNode2 = ((*tempNode0 - *tempNode5).GetLength() < (*tempNode0 - *tempNode6).GetLength()) ? tempNode5 : tempNode6;
	max3456 = ((*tempNode0 - *tempNode1).GetLength() < (*tempNode0 - *tempNode2).GetLength()) ? 1 : 2;//1��ʾ34�㣬2��ʾ56��


	//�ǵ����ұ߽����
	if ((max3456==1) && ((*tempNode0- *tempNode1).GetLength() < length_differ))
	{	
		//���4����
		if (fabs((point_corner[0]->point.y() - point_corner[4]->point.y())) < fabs((point_corner[3]->point.y() - point_corner[4]->point.y()) /2))
		{
			point_corner[0]->Id = point_corner[4]->Id;//����0�鲢����4
			/*point_corner[0]->point.x() = point_corner[4]->point.x();
			point_corner[0]->point.y() = point_corner[4]->point.y();*/
			if ((point_corner[2]->Id == point_corner[3]->Id) || (point_corner[2]->Id == point_corner[4]->Id))//�жϵ�2yu3,4�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
		}
		//���3����
		if (fabs((point_corner[0]->point.y() - point_corner[3]->point.y())) < fabs((point_corner[3]->point.y() - point_corner[4]->point.y()) / 2))
		{
			point_corner[0]->Id = point_corner[3]->Id;//����0�鲢����3
			/*point_corner[0]->point.x() = point_corner[3]->point.x();
			point_corner[0]->point.y() = point_corner[3]->point.y();*/
			if (point_corner[2]->Id == point_corner[3]->Id)//�жϵ�2yu3�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
		}
		if (index56 == 56)//��0��5,6
		{
			if ((point_corner[1]->Id == point_corner[5]->Id) || (point_corner[2]->Id == point_corner[6]->Id))//�жϵ�1��5,6�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����1���Դ���ϵɾ��				
			}
		}
		if (index56 == 6)//��0��6
		{
			if ((point_corner[1]->Id == point_corner[6]->Id))//�жϵ�1��6�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����1���Դ���ϵɾ��				
			}
		}			
	}	

	//�ǵ����ϱ߽����
	if ((max3456==2) && ((*tempNode0 - *tempNode2).GetLength() < length_differ))
	{
		//���6����
		if (fabs((point_corner[0]->point.x() - point_corner[6]->point.x())) < fabs((point_corner[5]->point.x() - point_corner[6]->point.x()) / 2))
		{
			point_corner[0]->Id = point_corner[6]->Id;//����0�鲢����4
			/*point_corner[0]->point.x() = point_corner[6]->point.x();
			point_corner[0]->point.y() = point_corner[6]->point.y();*/
			if ((point_corner[1]->Id == point_corner[5]->Id) || (point_corner[1]->Id == point_corner[6]->Id))//�жϵ�1yu5,6�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����1���Դ���ϵɾ��				
			}
		}
		//���5����
		if (fabs((point_corner[0]->point.x() - point_corner[5]->point.x())) < fabs((point_corner[5]->point.x() - point_corner[6]->point.x()) / 2))
		{
			point_corner[0]->Id = point_corner[5]->Id;//����0�鲢����5
			/*point_corner[0]->point.x() = point_corner[5]->point.x();
			point_corner[0]->point.y() = point_corner[5]->point.y();*/
			if (point_corner[1]->Id == point_corner[5]->Id)//�жϵ�1yu5�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����2���Դ���ϵɾ��				
			}
		}
		if (index34 == 34)//��0��34
		{
			if ((point_corner[2]->Id == point_corner[3]->Id) || (point_corner[2]->Id == point_corner[4]->Id))//�жϵ�2��3,4�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
		}
		if (index34 == 4)//��0��4
		{
			if ((point_corner[2]->Id == point_corner[4]->Id))//�жϵ�2��4�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
		}
	}

	//�������򳤶ȴ���2m
	if (((*tempNode0 - *tempNode1).GetLength() > length_differ) && ((*tempNode0 - *tempNode2).GetLength() > length_differ))
	{
		point_corner[0]->Id = nId;//�Ըõ����������			
		addspace.setX(point_corner[0]->point.x());
		addspace.setY(point_corner[0]->point.y());
		addspace.setID(nId);
		UVForOut.push_back(addspace);
		nId++;
		if (index34 == 34)
		{
			if ((point_corner[2]->Id == point_corner[3]->Id))//�жϵ�2��3�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[4]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[3]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
		}
		else
		{
			if ((point_corner[2]->Id == point_corner[4]->Id))//�жϵ�2��4�Ƿ��غ�
			{
				point_corner[2]->next.clear();//����2���Դ���ϵɾ��				
			}
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[4]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
		}

		if (index56 == 56)
		{
			if ((point_corner[1]->Id == point_corner[5]->Id))//�жϵ�1��5�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����2���Դ���ϵɾ��				
			}
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[5]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[6]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
		}
		else
		{
			if ((point_corner[1]->Id == point_corner[6]->Id))//�жϵ�1��6�Ƿ��غ�
			{
				point_corner[1]->next.clear();//����1���Դ���ϵɾ��				
			}
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[6]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
		}
	}

	delete tempNode0;
	/*delete tempNode1;
	delete tempNode2;*/
	delete tempNode3;
	delete tempNode4;
	delete tempNode5;
	delete tempNode6;
	
	return edge_m;
}


//******************************************ȡz���ֵ�㣬ֻ���Ŀ�������滮��****************************************************************
std::vector< Member> NurbsSurface_qj::SurfaceAverageSpit_Rectangle(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec)
{
	//������z����ֵ����
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ�Сֵ��
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ���ֵ��
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//�洢��ֵ��
	int num = 0;
	for (std::vector<PLib::Point2Dd>::iterator itr = pMinExtrmums.begin(); itr != pMinExtrmums.end(); ++itr)//����Сֵ���뼫ֵ������
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	//pMaxExtrmums.push_back(PLib::Point2Dd(0.5, 0.5));
	for (std::vector<PLib::Point2Dd>::iterator itr = pMaxExtrmums.begin(); itr != pMaxExtrmums.end(); ++itr)//������ֵ���뼫ֵ������
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}

	CPoint2D maxPoint;//�洢z���ĵ� 
	double MaxzofPoint;
	maxPoint = extrmums[0];
	Node *tempNode1;
	tempNode1 = GetSurfacePointatUandV(extrmums[0].x, extrmums[0].y);
	MaxzofPoint = fabs(tempNode1->getZ());
	delete tempNode1;
	for (int i = 1; i<num; i++)//Ѱ������z���ֵ	
	{
		tempNode1 = GetSurfacePointatUandV(extrmums[i].x, extrmums[i].y);
		if (fabs(tempNode1->getZ()) > MaxzofPoint)
		{
			MaxzofPoint = fabs(tempNode1->getZ());
			maxPoint = extrmums[i];
		}
		delete tempNode1;
	}
	//��ֵ�ĸ��ֿ�ı߽�
	rectangleRegion rec[4];
	rec[0].maxU = maxU; rec[0].minU = maxPoint.x; rec[0].maxV = maxV; rec[0].minV = maxPoint.y; rec[0].type = 0; rec[0].type = 0;
	rec[1].maxU = maxPoint.x; rec[1].minU = minU; rec[1].maxV = maxV; rec[1].minV = maxPoint.y; rec[1].type = 1; rec[1].type = 1;
	rec[2].maxU = maxPoint.x; rec[2].minU = minU; rec[2].maxV = maxPoint.y; rec[2].minV = minV; rec[2].type = 2; rec[2].type = 2;
	rec[3].maxU = maxU; rec[3].minU = maxPoint.x; rec[3].maxV = maxPoint.y; rec[3].minV = minV; rec[3].type = 3; rec[3].type = 3;

	

	//�ֿ����ɸ˼�
	std::vector< Member> mVec;//���ڼ�¼�ܵĸ����ֿ��ڵ�member
	std::vector<Member> mVecTemp;//��ʱ�����ڼ�¼ѭ���зֿ��ڵ�member���β�	
	std::vector<IdRelationship> IdRel_Rectangle;//��¼�ı����ĵ��ϵ����ʱ��

	//�ȷ��ڲ��߽磬����uvEdge��
	std::vector<PointStruct> uvEdge;
	mVecTemp = InsideEdge_Rectangle(minU, maxU, minV, maxV, maxPoint, splitLengthU, splitLengthV, nVec, uvEdge);
	(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//��mVecTemp�ӵ�mVec��,�ĸ��˼�

	for (int i = 0; i < 4; ++i)//������������еȷ֣������ȷֵĵ�͸���Ϣ�ֱ�洢��nVec��mVec�У�4�ֿ���		
	{
		mVecTemp = MemberofRec_Rectangle(rec[i], splitLengthU, splitLengthU, nVec, uvEdge/*,Gui_index*/);//��i����еȷ�,���ɸ˼�
		(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//��mVecTemp�ӵ�mVec��	
	}	

	return mVec;	
}

//�ڲ��߽紦��
std::vector<Member> NurbsSurface_qj::InsideEdge_Rectangle(double minU, double maxU, double minV, double maxV, CPoint2D maxPoint, double splitLengthU, double splitLengthV, std::vector< Node> & nVec, std::vector<PointStruct> &uvEdge)
{
	std::vector<double> temp;
	PointStruct divPoint;
	Node AddSpace;
	std::vector<Member> mVecTemp;
	Member m;

	divPoint.point.x() = maxU;//��0
	divPoint.point.y() = maxPoint.y;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxU);
	AddSpace.setY(maxPoint.y);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = maxPoint.x;//��1
	divPoint.point.y() =maxV ;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxPoint.x);
	AddSpace.setY(maxV);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = minU;//��2
	divPoint.point.y() = maxPoint.y;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(minU);
	AddSpace.setY(maxPoint.y);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = maxPoint.x;//��3
	divPoint.point.y() = minV;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxPoint.x);
	AddSpace.setY(minV);
	nVec.push_back(AddSpace);
	nId++;

	temp = getUAverageLength(minU, maxPoint.x, maxPoint.y, splitLengthU, 1);//��������ȷ֣���Ӧtype=1��2,��¼�ȷֵ�
	for (int j = 0; j < temp.size(); ++j)
	{
		divPoint.point.x() = temp[j];//�����temp��u���PointStruct���������uvEdge��
		divPoint.point.y() = maxPoint.y;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//���ȷֱ߽�����uvEdge		
		AddSpace.setID(nId);
		AddSpace.setX(temp[j]);
		AddSpace.setY(maxPoint.y);
		nVec.push_back(AddSpace);//�����Ժ��ٸ�ֵ����ΪuvEdge�е��п��ܱ仯
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[2]);//������2����
	mVecTemp.push_back(m);


	temp = getUAverageLength(maxPoint.x, maxU, maxPoint.y, splitLengthU, 0);//�������ҵȷ֣���Ӧtype=0,3����
	temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ�����Ѿ����
	for (int j = 0; j < temp.size(); ++j)
	{
		divPoint.point.x() = temp[j];//���Ұ�temp��u���PointStruct���������uvEdge��			
		divPoint.point.y() = maxPoint.y;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//���ȷֱ߽�����uvEdge		
		AddSpace.setID(nId);
		AddSpace.setX(temp[j]);
		AddSpace.setY(maxPoint.y);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[0]);//������2����
	mVecTemp.push_back(m);


	temp = getVAverageLength(minV, maxPoint.y, maxPoint.x, splitLengthV, 2);//�������µȷ֣���Ӧtype=2��3,��¼�ȷֵ�
	temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ�����Ѿ����
	for (int j = 0; j != temp.size(); ++j)
	{
		divPoint.point.y() = temp[j];//���°�temp��v���PointStruct���������vVec��				
		divPoint.point.x() = maxPoint.x;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//���ȷֱ߽�����uvEdge		
		AddSpace.setID(nId);
		AddSpace.setY(temp[j]);
		AddSpace.setX(maxPoint.x);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[3]);//������2����
	mVecTemp.push_back(m);


	temp = getVAverageLength(maxPoint.y, maxV, maxPoint.x, splitLengthV, 0);//�������ϵȷ֣���Ӧtype=0,1��
	temp.erase(temp.begin());//��ȥ��ֵ�㣬��Ϊ�����Ѿ����
	for (int j = 0; j != temp.size(); ++j)
	{
		divPoint.point.y() = temp[j];//���ϰ�temp��v���PointStruct���������vVec��					
		divPoint.point.x() = maxPoint.x;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//���ȷֱ߽�����uvEdge		
		AddSpace.setID(nId);
		AddSpace.setY(temp[j]);
		AddSpace.setX(maxPoint.x);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[1]);//������2����
	mVecTemp.push_back(m);

	return mVecTemp;
}

//�ֿ黭�˼�
std::vector<Member> NurbsSurface_qj::MemberofRec_Rectangle(rectangleRegion rec, double splitLengthU, double splitLengthV, std::vector<Node> &nVec, std::vector<PointStruct> &uvEdge/*, std::vector<bool> &Gui_index*/)
{
	//��������
	double maxU = rec.maxU + 0.2;
	//***********�˴�0.3�д��Ż�***************************************************
	double maxV = rec.maxV + 0.2;
	double minU = rec.minU - 0.2;
	double minV = rec.minV - 0.2;

	int u_num;
	std::vector<std::vector<PointStruct> > points(1);//�����άvector
	PointStruct AddPoint;//��������points��size
	AddPoint.Id = -1;//��ʼ�����
	AddPoint.point = PLib::Point2Dd(-1, -1);//��ʼ��uv
	
	int row0 = 0;
	int col0 = 0;
	int row1 = 1;
	int col1 = 0;

	double nextU = 0;
	double nextV = 0;
	PLib::Point2Dd nextCrossP;
	Node *tempNode1;
	Node *tempNode2;
	Node *tempNode3;
	double cosValue = 0;
	double u = 0;
	double v = 0;

	//switch (rec.type)
	//{
	//case 0://�������������ҵ�
	//	u = rec.minU;
	//	v = rec.minV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u <= maxU)//����һ�е�������ҵ�	
	//	{
	//		nextU = nextUPoint(u, v, splitLengthU);
	//		points[row0].push_back(AddPoint);
	//		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
	//		col0 = col0 + 1;
	//		u = nextU;
	//	}	

	//	for (row0 = 0, row1 = 1; v <= maxV; row0++, row1++)
	//	{
	//		u = rec.minU;
	//		col0 = 0; col1 = 0;
	//		points.resize(row1 + 1);//����һ��
	//		nextV = nextVPoint(u, v, splitLengthV);//Ѱ�ҵڶ��еĵ�һ����v����
	//	
	//		points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����
	//		

	//		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV);//��1��2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1 )//����ҵ�������Ϊ-1 -1����������
	//			{
	//				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV);//��1��3
	//					//��ӵ�֮���ϵ
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
	//					
	//					col1 = col1 + 1;//�ڶ�������һ��
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV);//��1��4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//��ӵ�֮���ϵ
	//						points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0+1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//���col0+2������
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));	
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));							
	//						break;
	//					}
	//				}
	//			}
	//			//�������
	//			points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0+1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));				
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 1://������������
	//	u = rec.maxU;
	//	v = rec.minV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u >=minU)//����һ�е�������ҵ�	
	//	{
	//		nextU = nextUPoint(u, v, splitLengthU,1);
	//		points[row0].push_back(AddPoint);
	//		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
	//		col0 = col0 + 1;
	//		u = nextU;
	//	}

	//	for (row0 = 0, row1 = 1; v <= maxV; row0++, row1++)
	//	{
	//		u = rec.maxU;
	//		col0 = 0; col1 = 0;
	//		points.resize(row1 + 1);//����һ��
	//		nextV = nextVPoint(u, v, splitLengthV,1);//Ѱ�ҵڶ��еĵ�һ����v����

	//		points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


	//		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,1);//��1��2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
	//			{
	//				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,1);//��1��3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//��ӵ�֮���ϵ
	//					points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//�ڶ�������һ��
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,1);//��1��4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//��ӵ�֮���ϵ
	//						points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//���col0+2������
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//�������
	//			points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 2://��������������
	//	u = rec.maxU;
	//	v = rec.maxV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u >=minU)//����һ�е�������ҵ�	
	//	{
	//		nextU = nextUPoint(u, v, splitLengthU,2);
	//		points[row0].push_back(AddPoint);
	//		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
	//		col0 = col0 + 1;
	//		u = nextU;
	//	}

	//	for (row0 = 0, row1 = 1; v >=minV; row0++, row1++)
	//	{
	//		u = rec.maxU;
	//		col0 = 0; col1 = 0;
	//		points.resize(row1 + 1);//����һ��
	//		nextV = nextVPoint(u, v, splitLengthV,2);//Ѱ�ҵڶ��еĵ�һ����v����

	//		points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


	//		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,2);//��1��2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
	//			{
	//				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,2);//��1��3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//��ӵ�֮���ϵ
	//					points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//�ڶ�������һ��
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,2);//��1��4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//��ӵ�֮���ϵ
	//						points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//���col0+2������
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//�������
	//			points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 3://��������������
	//	u = rec.minU;
	//	v = rec.maxV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u <= maxU)//����һ�е�������ҵ�	
	//	{
	//		nextU = nextUPoint(u, v, splitLengthU,3);
	//		points[row0].push_back(AddPoint);
	//		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
	//		col0 = col0 + 1;
	//		u = nextU;
	//	}

	//	for (row0 = 0, row1 = 1; v >=minV ; row0++, row1++)
	//	{
	//		u = rec.minU;
	//		col0 = 0; col1 = 0;
	//		points.resize(row1 + 1);//����һ��
	//		nextV = nextVPoint(u, v, splitLengthV,3);//Ѱ�ҵڶ��еĵ�һ����v����

	//		points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


	//		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,3);//��1��2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
	//			{
	//				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,3);//��1��3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//��ӵ�֮���ϵ
	//					points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//�ڶ�������һ��
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,3);//��1��4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//��ӵ�֮���ϵ
	//						points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//���col0+2������
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0+1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//�������
	//			points[row1].push_back(AddPoint);//���ڶ��е�size����1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//default:
	//	break;
	//}


switch (rec.type)
{
case 0://�������������ҵ�
	u = rec.minU;
	v = rec.minV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u <= maxU)//����һ�е�������ҵ�	
	{
		nextU = nextUPoint(u, v, splitLengthU);
		points[row0].push_back(AddPoint);
		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
		col0 = col0 + 1;
		u = nextU;
	}

	for (row0 = 0, row1 = 1; v <= maxV; row0++, row1++)
	{
		u = rec.minU;
		col0 = 0; col1 = 0;
		points.resize(row1 + 1);//����һ��
		nextV = nextVPoint(u, v, splitLengthV);//Ѱ�ҵڶ��еĵ�һ����v����

		points[row1].push_back(AddPoint);//���ڶ��е�size����1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV);//��1��2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
			{
				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV);//��1��3
					//��ӵ�֮���ϵ
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					points[row1].push_back(AddPoint);//���ڶ��е�size����1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//�ڶ�������һ��
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV);//��1��4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//��ӵ�֮���ϵ
						points[row1].push_back(AddPoint);//���ڶ��е�size����1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//���col0+2������
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//�������
			points[row1].push_back(AddPoint);//���ڶ��е�size����1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 1://������������
	u = rec.maxU;
	v = rec.minV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u >= minU)//����һ�е�������ҵ�	
	{
		nextU = nextUPoint(u, v, splitLengthU, 1);
		points[row0].push_back(AddPoint);
		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
		col0 = col0 + 1;
		u = nextU;
	}

	for (row0 = 0, row1 = 1; v <= maxV; row0++, row1++)
	{
		u = rec.maxU;
		col0 = 0; col1 = 0;
		points.resize(row1 + 1);//����һ��
		nextV = nextVPoint(u, v, splitLengthV, 1);//Ѱ�ҵڶ��еĵ�һ����v����

		points[row1].push_back(AddPoint);//���ڶ��е�size����1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 1);//��1��2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
			{
				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 1);//��1��3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//��ӵ�֮���ϵ
					points[row1].push_back(AddPoint);//���ڶ��е�size����1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//�ڶ�������һ��
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 1);//��1��4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//��ӵ�֮���ϵ
						points[row1].push_back(AddPoint);//���ڶ��е�size����1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//���col0+2������
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//�������
			points[row1].push_back(AddPoint);//���ڶ��е�size����1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 2://��������������
	u = rec.maxU;
	v = rec.maxV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u >= minU)//����һ�е�������ҵ�	
	{
		nextU = nextUPoint(u, v, splitLengthU, 2);
		points[row0].push_back(AddPoint);
		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
		col0 = col0 + 1;
		u = nextU;
	}

	for (row0 = 0, row1 = 1; v >= minV; row0++, row1++)
	{
		u = rec.maxU;
		col0 = 0; col1 = 0;
		points.resize(row1 + 1);//����һ��
		nextV = nextVPoint(u, v, splitLengthV, 2);//Ѱ�ҵڶ��еĵ�һ����v����

		points[row1].push_back(AddPoint);//���ڶ��е�size����1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 2);//��1��2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
			{
				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 2);//��1��3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//��ӵ�֮���ϵ
					points[row1].push_back(AddPoint);//���ڶ��е�size����1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//�ڶ�������һ��
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 2);//��1��4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//��ӵ�֮���ϵ
						points[row1].push_back(AddPoint);//���ڶ��е�size����1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//���col0+2������
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//�������
			points[row1].push_back(AddPoint);//���ڶ��е�size����1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 3://��������������
	u = rec.minU;
	v = rec.maxV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u <= maxU)//����һ�е�������ҵ�	
	{
		nextU = nextUPoint(u, v, splitLengthU, 3);
		points[row0].push_back(AddPoint);
		points[row0][col0 + 1].point = PLib::Point2Dd(nextU, v);
		col0 = col0 + 1;
		u = nextU;
	}

	for (row0 = 0, row1 = 1; v >= minV; row0++, row1++)
	{
		u = rec.minU;
		col0 = 0; col1 = 0;
		points.resize(row1 + 1);//����һ��
		nextV = nextVPoint(u, v, splitLengthV, 3);//Ѱ�ҵڶ��еĵ�һ����v����

		points[row1].push_back(AddPoint);//���ڶ��е�size����1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//�ڶ��е�һ����


		for (; col0 < points[row0].size() - 1;)//col0ȡ����һ�еĵ����ڶ�����
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//�ڶ��е�һ����
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//��һ�еڶ�����
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//��һ�е�һ����

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************�˴�Ҫɾ��

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 3);//��1��2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//����ҵ�������Ϊ-1 -1����������
			{
				if (cosValue < 0)//�����㹹�ɵ������εĽǶȴ���90�ȣ�����������ε������������һ��
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 3);//��1��3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//��ӵ�֮���ϵ
					points[row1].push_back(AddPoint);//���ڶ��е�size����1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//ֻ����1��3�㣬2��3���´���
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//�ڶ�������һ��
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//�����һ��col0+2�����
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 3);//��1��4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//��ӵ�֮���ϵ
						points[row1].push_back(AddPoint);//���ڶ��е�size����1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//���col0+2������
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//�������
			points[row1].push_back(AddPoint);//���ڶ��е�size����1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

default:
	break;
}
	//�˼�����	
	std::vector<Member> mDisplayVec; //��Ϣ�����洢�˵��������˵�˳�����ID��ȫһ��	
	Node addspace;//��������NodeForOut�����Ĵ�С
	Member m;
	std::vector<PointStruct> EdgePoints;//�洢�߽��
	//bool guibingPoint = false;//��¼�ϵ��Ƿ�Ϊ�ڲ���鲢������Ϊtrue��Ĭ��Ϊfalse
	//double error_guibing = 0.2;//�鲢���
	PointStruct AddPoint2;
	
	
	switch (rec.type)
	{
	case 0:
	{
			  //��points���б��,�ֿ��ڲ��㣬��ȥ�߽�
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()<rec.maxU && points[i][j].point.y()<rec.maxV)//�жϵ��ڷֿ��ڲ�
					  {
						  points[i][j].Id = nId;
						  addspace.setX(points[i][j].point.x());
						  addspace.setY(points[i][j].point.y());
						  addspace.setID(nId);
						  nVec.push_back(addspace);
						  ++nId;
					  }
				  }
			  }

			  //�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[0][j].point.x() == uvEdge[k].point.x()) && (points[0][j].point.y() == uvEdge[k].point.y()))
						  {
							  points[0][j].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }
			  //�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//ȷ�����ڷֿ���//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[i][0].point.x() == uvEdge[k].point.x()) && (points[i][0].point.y() == uvEdge[k].point.y()))
						  {
							  points[i][0].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }

			  //���±߽翪ʼ���ӵ�һ�п�ʼ���˴�j��0��ʼ��j�������¹��򣬷ֿ�0��2��0��ʼ���ֿ�1��3��1��ʼ���������Է�ֹ��ֵ�㴦�ĸ˼��ظ�
			  EdgePoints.push_back(uvEdge[0]);//����0���뵽�߽��
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//�ж��Ƿ�Ϊ���ҵ�
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//�ж����ҵ�������ϵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.maxU, rec.maxV, 0);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  for (size_t k = 0; k < points[0][j].next.size(); k++)//�˴��������������
						  {
							  m = Add_Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //���ڲ��㽨���˼�����ȥ�±߽硢��߽�
			  for (size_t i = 1; i < points.size(); i++)//��������
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//��������
				  {
					  if (points[i][j].Id != -1)//�ҵ��ڲ���
					  {
						  if (points[i][j + 1].Id == -1)//��i��j���ҵ����ⲿ
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.maxU, rec.maxV, 0);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i][j + 1]);
							  mDisplayVec.push_back(m);
						  }

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//��i��j���ϵ����ⲿ���˴�û�п������������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.maxU, rec.maxV, 0);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i + 1][points[i][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //������߽�,�ӵڶ��п�ʼ

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//�ж��Ƿ�Ϊ���ϵ�
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//�ж����ϵ�������ҵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.maxU, rec.maxV, 0);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  for (size_t k = 0; k < points[i][0].next.size(); k++)//�˴��������������
						  {
							  m = Add_Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }
			  EdgePoints.push_back(uvEdge[1]);

			  for (size_t i = 0; i < EdgePoints.size() - 1; i++)
			  {
				  m = Add_Member(EdgePoints[i], EdgePoints[i + 1]);
				  mDisplayVec.push_back(m);
			  }

			  break;
	}

	case 1:
	{
			  //��points���б��,�ֿ��ڲ��㣬��ȥ�߽�
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()>rec.minU && points[i][j].point.y()<rec.maxV)//�жϵ��ڷֿ��ڲ�
					  {
						  points[i][j].Id = nId;
						  addspace.setX(points[i][j].point.x());
						  addspace.setY(points[i][j].point.y());
						  addspace.setID(nId);
						  nVec.push_back(addspace);
						  ++nId;
					  }
				  }
			  }

			  //�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[0][j].point.x() == uvEdge[k].point.x()) && (points[0][j].point.y() == uvEdge[k].point.y()))
						  {
							  points[0][j].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }
			  //�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//ȷ�����ڷֿ���//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[i][0].point.x() == uvEdge[k].point.x()) && (points[i][0].point.y() == uvEdge[k].point.y()))
						  {
							  points[i][0].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }

			  //���±߽翪ʼ���ӵ�һ�п�ʼ���˴�j��1��ʼ��j�������¹��򣬷ֿ�0��2��0��ʼ���ֿ�1��3��1��ʼ���������Է�ֹ��ֵ�㴦�ĸ˼��ظ�
			  EdgePoints.push_back(uvEdge[2]);//����0���뵽�߽��
			  for (int j = 1; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//�ж��Ƿ�Ϊ���ҵ�
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//�ж����ҵ�������ϵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.minU, rec.maxV, 1);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  
							  m = Add_Member(points[0][j], points[points[0][j].next[1].x][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  
					  }
				  }
			  }

			  //���ڲ��㽨���˼�����ȥ�±߽硢��߽�
			  for (size_t i = 1; i < points.size(); i++)//��������
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//��������
				  {
					  if (points[i][j].Id != -1)//�ҵ��ڲ���
					  {
						  if (points[i][j + 1].Id == -1)//��i��j���ҵ����ⲿ
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.minU, rec.maxV, 1);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i][j + 1]);
							  mDisplayVec.push_back(m);
						  }

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//��i��j���ϵ����ⲿ���˴�û�п������������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.minU, rec.maxV, 1);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i + 1][points[i][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //������߽�,�ӵڶ��п�ʼ

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//�ж��Ƿ�Ϊ���ϵ�
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//�ж����ϵ�������ҵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.minU, rec.maxV, 1);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  
							m = Add_Member(points[i][0], points[points[i][0].next[0].x][points[i][0].next[0].y]);//ֻ��next��0����next��1���ڷֲ�0���Ѿ�����
							mDisplayVec.push_back(m);
						  
					  }
				  }
			  }
			  EdgePoints.push_back(uvEdge[1]);

			  for (size_t i = 0; i < EdgePoints.size() - 1; i++)
			  {
				  m = Add_Member(EdgePoints[i], EdgePoints[i + 1]);
				  mDisplayVec.push_back(m);
			  }
			  break;
	}

	case 2:
	{
			  //��points���б��,�ֿ��ڲ��㣬��ȥ�߽�
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()>rec.minU && points[i][j].point.y()>rec.minV)//�жϵ��ڷֿ��ڲ�
					  {
						  points[i][j].Id = nId;
						  addspace.setX(points[i][j].point.x());
						  addspace.setY(points[i][j].point.y());
						  addspace.setID(nId);
						  nVec.push_back(addspace);
						  ++nId;
					  }
				  }
			  }

			  //�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[0][j].point.x() == uvEdge[k].point.x()) && (points[0][j].point.y() == uvEdge[k].point.y()))
						  {
							  points[0][j].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }
			  //�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//ȷ�����ڷֿ���//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[i][0].point.x() == uvEdge[k].point.x()) && (points[i][0].point.y() == uvEdge[k].point.y()))
						  {
							  points[i][0].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }

			  //���±߽翪ʼ���ӵ�һ�п�ʼ���˴�j��0��ʼ��j�������¹��򣬷ֿ�0��2��0��ʼ���ֿ�1��3��1��ʼ���������Է�ֹ��ֵ�㴦�ĸ˼��ظ�
			  EdgePoints.push_back(uvEdge[2]);//����0���뵽�߽��
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//�ж��Ƿ�Ϊ���ҵ�
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//�ж����ҵ�������ϵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.minU, rec.minV,2);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  for (size_t k = 0; k < points[0][j].next.size(); k++)//�˴��������������
						  {
							  m = Add_Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);
							  mDisplayVec.push_back(m);
						  }						  
					  }
				  }
			  }

			  //���ڲ��㽨���˼�����ȥ�±߽硢��߽�
			  for (size_t i = 1; i < points.size(); i++)//��������
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//��������
				  {
					  if (points[i][j].Id != -1)//�ҵ��ڲ���
					  {
						  if (points[i][j + 1].Id == -1)//��i��j���ҵ����ⲿ
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.minU, rec.minV, 2);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i][j + 1]);
							  mDisplayVec.push_back(m);
						  }

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//��i��j���ϵ����ⲿ���˴�û�п������������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.minU, rec.minV, 2);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i + 1][points[i][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //������߽�,�ӵڶ��п�ʼ

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//�ж��Ƿ�Ϊ���ϵ�
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//�ж����ϵ�������ҵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.minU, rec.minV, 2);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  for (size_t k = 0; k < points[i][0].next.size(); k++)//�˴��������������
						  {
							  m = Add_Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }
			  EdgePoints.push_back(uvEdge[3]);

			  for (size_t i = 0; i < EdgePoints.size() - 1; i++)
			  {
				  m = Add_Member(EdgePoints[i], EdgePoints[i + 1]);
				  mDisplayVec.push_back(m);
			  }
			  break;
	}

	case 3:
	{
			  //��points���б��,�ֿ��ڲ��㣬��ȥ�߽�
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()<rec.maxU && points[i][j].point.y()>rec.minV)//�жϵ��ڷֿ��ڲ�
					  {
						  points[i][j].Id = nId;
						  addspace.setX(points[i][j].point.x());
						  addspace.setY(points[i][j].point.y());
						  addspace.setID(nId);
						  nVec.push_back(addspace);
						  ++nId;
					  }
				  }
			  }

			  //�ж��±߽��ϵĵ�,��uvEdge�е������Ÿ����±߽��ϵĵ�
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[0][j].point.x() == uvEdge[k].point.x()) && (points[0][j].point.y() == uvEdge[k].point.y()))
						  {
							  points[0][j].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }
			  //�ж���߽��ϵĵ�,��uvEdge�е������Ÿ�����߽��ϵĵ�,�˴�i=1���������½ǵ��������Ѿ������
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//ȷ�����ڷֿ���//*****************************��䶯
				  {
					  for (size_t k = 0; k < uvEdge.size(); ++k)
					  {
						  if ((points[i][0].point.x() == uvEdge[k].point.x()) && (points[i][0].point.y() == uvEdge[k].point.y()))
						  {
							  points[i][0].Id = uvEdge[k].Id;
							  break;
						  }
					  }
				  }
			  }

			  //���±߽翪ʼ���ӵ�һ�п�ʼ���˴�j��1��ʼ��j�������¹��򣬷ֿ�0��2��0��ʼ���ֿ�1��3��1��ʼ���������Է�ֹ��ֵ�㴦�ĸ˼��ظ�
			  EdgePoints.push_back(uvEdge[0]);//����0���뵽�߽��
			  for (int j = 1; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//�ж��Ƿ�Ϊ���ҵ�
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//�ж����ҵ�������ϵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.maxU, rec.minV, 3);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  
							  m = Add_Member(points[0][j], points[points[0][j].next[1].x][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  
					  }
				  }
			  }

			  //���ڲ��㽨���˼�����ȥ�±߽硢��߽�
			  for (size_t i = 1; i < points.size(); i++)//��������
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//��������
				  {
					  if (points[i][j].Id != -1)//�ҵ��ڲ���
					  {
						  if (points[i][j + 1].Id == -1)//��i��j���ҵ����ⲿ
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.maxU, rec.minV, 3);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i][j + 1]);
							  mDisplayVec.push_back(m);
						  }

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//��i��j���ϵ����ⲿ���˴�û�п������������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.maxU, rec.minV, 3);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
						  else
						  {
							  m = Add_Member(points[i][j], points[i + 1][points[i][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //������߽�,�ӵڶ��п�ʼ

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//ȷ�����ڷֿ��ڣ�//*****************************��䶯
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//�ж��Ƿ�Ϊ���ϵ�
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//�ж����ϵ�������ҵ��Ƿ����ڲ�����Ϊ��
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//��ӱ߽�˼�
						  }
						  else//�ϵ�������
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.maxU, rec.minV, 3);//�����±߽��
							  EdgePoints.push_back(AddPoint2);//���õ���ӵ��߽��
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//�����ڲ���
					  {
						  
							  m = Add_Member(points[i][0], points[points[i][0].next[0].x][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);
						  
					  }
				  }
			  }
			  EdgePoints.push_back(uvEdge[3]);

			  for (size_t i = 0; i < EdgePoints.size() - 1; i++)
			  {
				  m = Add_Member(EdgePoints[i], EdgePoints[i + 1]);
				  mDisplayVec.push_back(m);
			  }
			  break;
	}
	default:
		break;
	}
	
	
	
	return mDisplayVec;
}

//��û�õ�
//�����ı��ζԽǵ�������������߱߳���ȷ����һ��������,type��ʾѰ����һ������ͣ� Ĭ��Ϊtype=0����ʾ������������Ѱ�ң�
//type=1��������������Ѱ�ң�type=2�������ϵ����£�type =3�������ϵ�����
PLib::Point2Dd NurbsSurface_qj::NextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type, double error)
//
{
	if ((point2.x() == -1 && point2.y() == -1) || (point1.x() == -1 && point1.y() == -1))
	{
		return PLib::Point2Dd(-1, -1);

	}
	//cout<<point1.x()<<"	"<<point1.y()<<"	"<<point2.x()<<"	"<<point2.y()<<endl;
	int Uflag = 0;
	int Vflag = 0;
	int preUFlag = 0;
	int preVFlag = 0;
	double u;
	double v;

	double deltaU = 0.005;
	double deltaV = 0.005;
	Node *tempNode1;
	Node *tempNode2;
	Node *tempNode3;
	double r1 = 0;
	double r2 = 0;
	tempNode1 = GetSurfacePointatUandV(point1.x(), point1.y());
	tempNode2 = GetSurfacePointatUandV(point2.x(), point2.y());
	unsigned int count = 0;
	switch (type)
	{
	case 0://�������������ҵ�
	{
			   u = (point1.x() + point2.x()) / 2 + (point1.y() - point2.y()) / 2;
			   v = (point1.y() + point2.y()) / 2 + (point2.x() - point1.x()) / 2;
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2�ϴ���Ҫ��С
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= 0)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < 1e-10 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 > 0)//r2��С����Ҫ����
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1�ϴ���Ҫ��С
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= 0)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < 1e-10 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 > 0)//r2��С����Ҫ����
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 1://������������Ѱ��
	{
			   u = (point1.x() + point2.x()) / 2 - (point1.y() - point2.y()) / 2;
			   v = (point1.y() + point2.y()) / 2 - (point2.x() - point1.x()) / 2;
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2�ϴ���Ҫ��С
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 > 0)//r2��С����Ҫ����
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1��С����Ҫ��Сuʹr1����
				   {
					   if (Uflag == 1)//˵����һ����Ҫ����u�����˵���Ҫ��Сu�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= -0.3)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < -0.3 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 < 0)//r1�ϴ���Ҫ����uʹr1��С
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 2://�����ϵ�����
	{
			   u = (point1.x() + point2.x()) / 2 + (point1.y() - point2.y()) / 2;
			   v = (point1.y() + point2.y()) / 2 + (point2.x() - point1.x()) / 2;
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2��С����Ҫ��Сvʹr2����
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 < 0)//r2�ϴ���Ҫ����vʹr2��С
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1��С����Ҫ��Сuʹr1����
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= -0.3)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < -0.3 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 < 0)//r1�ϴ���Ҫ����u��r1��С
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;

						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;
	case 3://�����ϵ�����
	{
			   u = (point1.x() + point2.x()) / 2 - (point1.y() - point2.y()) / 2;
			   v = (point1.y() + point2.y()) / 2 - (point2.x() - point1.x()) / 2;
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//���������ľ���
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2��С����Ҫ��Сvʹr2����
				   {

					   if (Vflag == 1)//˵����һ����Ҫ���󣬶��˲���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v - deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v + deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   if (v >= -0.3)
						   {
							   v = v - deltaV;
						   }
						   else
						   {
							   v = v + deltaV;
							   deltaV = deltaV / 2;
							   if (v < -0.3 && deltaV < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }
						   preVFlag = Vflag;
						   Vflag = -1;
					   }
				   }
				   else if (vChrodLength - r2 < 0)//r2�ϴ���Ҫ����vʹr2��С
				   {
					   if (Vflag == -1)//˵����һ����Ҫ��С�����˲���Ҫ���󣬿ɶ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   v = v + deltaV;//�˻ص���һ��
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//����
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1�ϴ���Ҫ��С
				   {
					   if (Uflag == 1)//˵����һ����Ҫ���󣬶��˵���Ҫ��С�����Զ϶�ҪѰ�ҵĵ�����һ��ʹ˵�֮��
					   {
						   u = u - deltaU;
						   deltaU = deltaU / 2;
						   u = u + deltaU;
						   Uflag = 0;
					   }
					   else
					   {

						   if (u >= 0)
						   {
							   u = u - deltaU;

						   }
						   else
						   {
							   u = u + deltaU;
							   deltaU = deltaU / 2;
							   if (u < 1e-10 && deltaU < 1e-10)
							   {
								   return PLib::Point2Dd(-1, -1);
							   }
						   }


						   Uflag = -1;
					   }

				   }
				   else if (uChrodLength - r1 > 0)//r2��С����Ҫ����
				   {
					   if (Uflag == -1)
					   {
						   u = u + deltaU;
						   deltaU = deltaU / 2;
						   u = u - deltaU;
						   Uflag = 0;
					   }
					   else
					   {
						   u = u + deltaU;
						   Uflag = 1;
					   }

				   }

			   } while (abs(r1 - vChrodLength) > error || abs(r2 - uChrodLength) > error);

	}
		break;

	default:
		break;
	}
	delete tempNode1;
	delete tempNode2;
	return PLib::Point2Dd(u, v);
}


//���������uvֵ�������֮��ı߽�㣬���õ��u��vΪ1,u_edgeΪ�߽�uֵ,01��ʾ�ϲ����飬23��ʾ�²�����
PointStruct NurbsSurface_qj::PointBtween_Rectangle(PointStruct p1, PointStruct p2, double u_edge, double v_edge,int type)
{
	PointStruct p3;

	switch (type)
	{
	case 0:
	case 1:
	{
			  if ((p1.point.x() - u_edge)*(p2.point.x() - u_edge) < 0)//�����u����
			  {
				  p3.point.x() = u_edge;
				  p3.point.y() = (u_edge - p1.point.x()) / (p2.point.x() - p1.point.x())*p2.point.y() +
					  (p2.point.x() - u_edge) / (p2.point.x() - p1.point.x())*p1.point.y();
				  if (p3.point.y() <= v_edge)//p3���ڲ�
				  {
					  p3.Id = nId;
					  nId++;
					  break;
				  }
			  }
			  if ((p1.point.y() - v_edge)*(p2.point.y() - v_edge) < 0)//�����v����
			  {
				  p3.point.y() = v_edge;
				  p3.point.x() = (v_edge - p1.point.y()) / (p2.point.y() - p1.point.y())*p2.point.x() +
					  (p2.point.y() - v_edge) / (p2.point.y() - p1.point.y())*p1.point.x();
				  p3.Id = nId;
				  nId++;
			  }
			  break;
	}
	case 2:
	case 3:
	{
			  if ((p1.point.x() - u_edge)*(p2.point.x() - u_edge) < 0)//�����u����
			  {
				  p3.point.x() = u_edge;
				  p3.point.y() = (u_edge - p1.point.x()) / (p2.point.x() - p1.point.x())*p2.point.y() +
					  (p2.point.x() - u_edge) / (p2.point.x() - p1.point.x())*p1.point.y();
				  if (p3.point.y() >= v_edge)//p3���ڲ�
				  {
					  p3.Id = nId;
					  nId++;
					  break;
				  }
			  }
			  if ((p1.point.y() - v_edge)*(p2.point.y() - v_edge) < 0)//�����v����
			  {
				  p3.point.y() = v_edge;
				  p3.point.x() = (v_edge - p1.point.y()) / (p2.point.y() - p1.point.y())*p2.point.x() +
					  (p2.point.y() - v_edge) / (p2.point.y() - p1.point.y())*p1.point.x();
				  p3.Id = nId;
				  nId++;
			  }
			  break;
	}
	
	default:
		break;
	}
	return p3;
}
















//******************************************ȡz���ֵ�㣬ֻ���Ŀ�������滮��****************************************************************

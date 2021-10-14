#include "NurbsSurface_qj.h"
//#include <nurbsS.h>

bool compare1_x( PointStruct p1,  PointStruct p2)//比较大小,从小到大
{
	return p1.point.x() < p2.point.x();
}
bool compare2_x(PointStruct p1, PointStruct p2)//比较大小,从大到小
{
	return p1.point.x() > p2.point.x();
}
bool compare3_y(PointStruct p1, PointStruct p2)//比较大小,从小到大
{
	return p1.point.y() < p2.point.y();
}
bool compare4_y(PointStruct p1, PointStruct p2)//比较大小,从大到小
{
	return p1.point.y() > p2.point.y();
}

bool edge_compare1(PointStruct *p1, PointStruct *p2)//用于不规则边界旁的内部点排序
{
	return p1->point.x() < p2->point.x();
}
bool edge_compare2(PointStruct *p1, PointStruct *p2)//用于不规则边界旁的内部点排序
{
	return p1->point.x() > p2->point.x();
}
bool edge_compare3(PointStruct *p1, PointStruct *p2)//用于不规则边界旁的内部点排序
{
	return p1->point.y() < p2->point.y();
}
bool edge_compare4(PointStruct *p1, PointStruct *p2)//用于不规则边界旁的内部点排序
{
	return p1->point.y() > p2->point.y();
}


bool compare_CPoint2D(CPoint2D p1, CPoint2D p2)//比较大小,从大到小
{
	return p1.y > p2.y;
}



NurbsSurface_qj::NurbsSurface_qj(void)
{
}


NurbsSurface_qj::~NurbsSurface_qj(void)
{
}



NurbsSurface_qj::NurbsSurface_qj(Vector_HPoint3Dd ccurve_points, int cdegree, Vector_HPoint3Dd tcurve_points, int tdegree)//给定曲线C的采样点和度，给定曲线T的采样点和度，求曲面
{

	this->CDegree = cdegree;//给度赋值
	CCurvepoint.resize(ccurve_points.size());

	//******************给曲线上的点赋值********************
	for (int i = 0; i < ccurve_points.size(); i++)
	{
		this->CCurvepoint[i] = ccurve_points[i];
	}
	//******************给曲线上的点赋值********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//构造曲线C

	this->TDegree = tdegree;//给度赋值
	TCurvepoint.resize(tcurve_points.size());

	//******************给曲线上的点赋值********************
	for (int i = 0; i < tcurve_points.size(); i++)
	{
		this->TCurvepoint[i] = tcurve_points[i];
	}
	//******************给曲线上的点赋值********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//构造曲线T

	nurbssurface.sweep(tnurbscurve, cnurbscurve, 5);
}

NurbsSurface_qj::NurbsSurface_qj(std::vector<Node*> Surface_points)//给定曲面的采样点，求曲面
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

NurbsSurface_qj::NurbsSurface_qj(std::vector<Node*> ccurve_points, int cdegree, std::vector<Node*> tcurve_points, int tdegree)//给定曲线C的采样点是Node形式和度，给定曲线T的采样点是Node形式和度，求曲面	
{
	Vector_HPoint3Dd theccurve_points;
	int CDegree;

	this->CDegree = cdegree;//给度赋值

	//********************给曲线上的采样点赋值***********************
	for (int i = 0; i <= ccurve_points.size(); i++)
	{

		(theccurve_points[i]).data[0] = (ccurve_points[0]->getX());
		(theccurve_points[i]).data[1] = (ccurve_points[0]->getY());
		(theccurve_points[i]).data[2] = (ccurve_points[0]->getZ());
		(theccurve_points[i]).data[3] = 1;

		this->CCurvepoint[i] = theccurve_points[i];
	}
	//********************给曲线上的采样点赋值***********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//构造曲线

	Vector_HPoint3Dd thetcurve_points;
	int TDegree;
	this->TDegree = tdegree;//给度赋值

	//********************给曲线上的采样点赋值***********************
	for (int i = 0; i <= tcurve_points.size(); i++)
	{

		(thetcurve_points[i]).data[0] = (tcurve_points[0]->getX());
		(thetcurve_points[i]).data[1] = (tcurve_points[0]->getY());
		(thetcurve_points[i]).data[2] = (tcurve_points[0]->getZ());
		(thetcurve_points[i]).data[3] = 1;

		this->TCurvepoint[i] = thetcurve_points[i];
	}
	//********************给曲线上的采样点赋值***********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//构造曲线
	nurbssurface.sweep(tnurbscurve, cnurbscurve, 5);
}


//读入参数，uv次数，节点矢量，控制点
NurbsSurface_qj::NurbsSurface_qj(int  uDeg, int vDeg, const Vector_DOUBLE& uKnots, const Vector_DOUBLE& vKnots, const Matrix_HPoint3Dd& controlPoints, unsigned long long int n_Id, unsigned long long int m_Id)
{
	nId = n_Id;
	mId = m_Id;
	PlNurbsSurfaced tempNurbsSurface(uDeg, vDeg, uKnots, vKnots, controlPoints);
	nurbssurface = tempNurbsSurface;
}


Node* NurbsSurface_qj::GetSurfacePointatUandV(double u, double v)//给定曲面和比例u，v,求曲面上的点
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
	this->CDegree = cdegree;//给度赋值
	CCurvepoint.resize(ccurve_points.size());

	//******************给曲线上的点赋值********************
	for (int i = 0; i < ccurve_points.size(); i++)
	{
		this->CCurvepoint[i] = ccurve_points[i];
	}
	//******************给曲线上的点赋值********************

	this->cnurbscurve.globalInterpH(CCurvepoint, CDegree);//构造曲线C
	this->TDegree = tdegree;//给度赋值
	TCurvepoint.resize(tcurve_points.size());

	//******************给曲线上的点赋值********************
	for (int i = 0; i < tcurve_points.size(); i++)
	{
		this->TCurvepoint[i] = tcurve_points[i];
	}
	//******************给曲线上的点赋值********************

	this->tnurbscurve.globalInterpH(TCurvepoint, TDegree);//构造曲线T

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


//获取曲面上一条曲线的极值点,v值固定,u 值可变，(u的最小值，u的最大值，v值）
std::vector<double> NurbsSurface_qj::getUExtremums(double min, double max, const double v)//获取曲面上一条曲线的极值点,v值固定,u 值可变
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//设定步长
	double u = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 1, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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


//获取曲面上一条曲线的极值点,u值固定,v 值可变,(v的最小值，v的最大值，u值）
std::vector<double> NurbsSurface_qj::getVExtremums(double min, double max, const double u)//获取曲面上一条曲线的极值点,u值固定,v 值可变
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//设定步长
	double v = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 1, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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


//获取u方向的拐点,(u的最小值，u的最大值，v值）
std::vector<double> NurbsSurface_qj::getUInflectionPoints(double min, double max, const double v)//获取u方向的拐点
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//设定步长
	double u = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 2, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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


//获取v方向的拐点,(v的最小值，v的最大值，u值）
std::vector<double> NurbsSurface_qj::getVInflectionPoints(double min, double max, const double u)
{
	std::vector<double> vec;
	double length = max - min;
	double deltaD = length / 20;//设定步长
	double v = min;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	nurbssurface.deriveAtH(u, v, 2, mat);

	Node *node1;
	Node *node2;
	node1 = GetSurfacePointatUandV(u, v);

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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
	Node node4 = (*node1 - *node3)*(*node2 - *node3);//叉积为(min, v)点，(max, v)点与((min + max)/2, v)点所构成平面的法向量l
	Node node5 = (*node1 - *node2)*node4;//*node5为所(min, v)点，(max, v)点和法向量l所构成平面的法向量
	double length = max - min;
	double deltaD = length / 20;


	Node tempNode1;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	double u = min;

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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
	Node node4 = (*node1 - *node3)*(*node2 - *node3);//叉积为(u, min)点，(u, max)点与(u, (min + max)/2)点所构成平面的法向量l
	Node node5 = (*node1 - *node2)*node4;//*node5为所(u, min)点，(u, max)点和法向量l所构成平面的法向量
	double length = max - min;
	double deltaD = length / 20;


	Node tempNode1;

	PLib::Matrix<PLib::HPoint3Dd> mat;
	double v = min;

	int preFlagZ = 0;//记录上次标志值，如果上次标志值的符号与该次的不同，而且曲线上两点之间的距离满足精度要求，则认为两点的中点是极值点
	int flagZ = 0;//记录此次标志

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


//v值固定，将u方向曲线按等弦长分割,
std::vector<double> NurbsSurface_qj::getUAveragePerSection(double min, double max, const double v, int splitNum)//v值固定，将u方向曲线按等弦长分割
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
	for (int i = 0; i < splitNum; ++i)//计算每段的弦长和总弦长，总弦长尾各段弦长之和；
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
	double R = totalStringLength / splitNum;//计算平均弦长,并以此做为初始的半径
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


//u值固定，将v方向曲线按等弦长分割
std::vector<double> NurbsSurface_qj::getVAveragePersection(double min, double max, const double u, int splitNum)//u值固定，将v方向曲线按等弦长分割
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
	for (int i = 0; i < splitNum; ++i)//计算每段的弦长和总弦长，总弦长尾各段弦长之和；
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
	double R = totalStringLength / splitNum;//计算平均弦长,并以此做为初始的半径
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
	std::vector<double> tmp1;//保存每一步的点的u值
	std::vector<double> tmp2;//保存每一步的点的u值
	tmp1 = getUExtremums(min, max, v);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//将极值点加入
	tmp1.clear();
	tmp1 = getUInflectionPoints(min, max, v);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//将拐点加入
	tmp1.clear();
	tmp2.push_back(min);//将起始点加入
	tmp2.push_back(max);//将终止点加入
	std::sort(tmp2.begin(), tmp2.end());//将起始点、终止点、拐点按从小到大排序
	result.insert(result.end(), tmp2.begin(), tmp2.end());
	for (size_t i = 0; i < result.size() - 1; ++i)
	{
		tmp1 = getUParallelPoints(result[i], result[i + 1], v);//获取与起点和终点平行的点
		tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());
		tmp1.clear();
	}
	std::sort(tmp2.begin(), tmp2.end());//将起始点、终止点、拐点以及平行点按从小到大排序


	int * splitNumPerSection = new int[tmp2.size()];//存储每段按比例分多少段
	Node *tempNode1;
	Node *tempNode2;
	double *l = new double[tmp2.size()];//记录每段的弦长
	double totalLength = 0;//记录总的长度

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		tempNode1 = GetSurfacePointatUandV(tmp2[i], v);
		tempNode2 = GetSurfacePointatUandV(tmp2[i + 1], v);
		l[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		//计算曲线上两点tempNode1和tempNode2的距离
		totalLength += l[i];
		delete tempNode1;
		delete tempNode2;
	}

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		splitNumPerSection[i] = int(l[i] / totalLength*splitNum + 0.5);//根据每段弦长占总长度的比例，确定分割多少次
	}

	result.clear();
	for (size_t i = 0; i < tmp2.size() - 1; ++i)
	{

		tmp1 = getUAveragePerSection(tmp2[i], tmp2[i + 1], v, splitNumPerSection[i]);
		//获取每段等分点的坐标

		result.insert(result.end(), tmp1.begin(), tmp1.end());
	}
	delete[] l;
	delete[] splitNumPerSection;
	return result;

}

//将v值固定，u方向上的曲线按等弦长划分；type=0表示从下到上，type = 1表示从上到下
//设定每份的长度，按一定长度将u方向曲线等分
std::vector<double> NurbsSurface_qj::getUAverageLength(double min, double max, const double v, double lengthPerSection, int type)//设定每份的长度，按一定长度将u方向曲线等分
{
	double length = max - min;
	std::vector<double> vec;
	double nextU = min;
	switch (type)
	{
	case 0://0，3都表示从左往右寻找等分点
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


//将u值固定，v方向上的曲线按等弦长划分；type=0表示从下到上，type = 1表示从上到下
//设定每份的长度，按一定长度将v方向曲线等分
std::vector<double> NurbsSurface_qj::getVAverageLength(double min, double max, const double u, double lengthPerSection, int type)//设定每份的长度，按一定长度将v方向曲线等分
{
	double length = max - min;
	std::vector<double> vec;
	double nextV;
	switch (type)
	{
	case 0://从下往上
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


//设定均分多少份，u值固定，将v方向曲线按等弦长分割
std::vector<double> NurbsSurface_qj::getVAverage(double min, double max, const double u, int splitNum)
{
	std::vector<double> result;
	std::vector<double> tmp1;//保存每一步的点的u值
	std::vector<double> tmp2;//保存每一步的点的u值
	tmp1 = getVExtremums(min, max, u);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//将极值点加入
	tmp1.clear();
	tmp1 = getVInflectionPoints(min, max, u);
	tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());//将拐点加入
	tmp1.clear();
	tmp2.push_back(min);//将起始点加入
	tmp2.push_back(max);//将终止点加入
	std::sort(tmp2.begin(), tmp2.end());//将起始点、终止点、拐点按从小到大排序
	result.insert(result.end(), tmp2.begin(), tmp2.end());
	for (size_t i = 0; i < result.size() - 1; ++i)
	{
		tmp1 = getVParallelPoints(result[i], result[i + 1], u);//获取与起点和终点平行的点
		tmp2.insert(tmp2.end(), tmp1.begin(), tmp1.end());
		tmp1.clear();
	}
	std::sort(tmp2.begin(), tmp2.end());//将起始点、终止点、拐点以及平行点按从小到大排序


	int * splitNumPerSection = new int[tmp2.size()];//存储每段按比例分多少段
	Node *tempNode1;
	Node *tempNode2;
	double *l = new double[tmp2.size()];//记录每段的弦长
	double totalLength = 0;//记录总的长度

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		tempNode1 = GetSurfacePointatUandV(u, tmp2[i]);
		tempNode2 = GetSurfacePointatUandV(u, tmp2[i + 1]);
		l[i] = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
			(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
			(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
		//计算曲线上两点tempNode1和tempNode2的距离
		totalLength += l[i];
		delete tempNode1;
		delete tempNode2;
	}

	for (int i = 0; i < tmp2.size() - 1; ++i)
	{
		splitNumPerSection[i] = int(l[i] / totalLength*splitNum + 0.5);//根据每段弦长占总长度的比例，确定分割多少次
	}

	result.clear();
	for (size_t i = 0; i < tmp2.size() - 1; ++i)
	{

		tmp1 = getVAveragePersection(tmp2[i], tmp2[i + 1], u, splitNumPerSection[i]);
		//获取每段等分点的坐标

		result.insert(result.end(), tmp1.begin(), tmp1.end());
	}
	delete[] l;
	delete[] splitNumPerSection;
	return result;
}


//寻找曲面的极大值点
std::vector<PLib::Point2Dd> NurbsSurface_qj::GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV)//寻找曲面的极大值点
{
	double lengthU = maxU - minU;
	double lenghtV = maxV - minV;

	double deltaU = lengthU / 50; //u方向初始步长
	double deltaV = lenghtV / 50; //v方向初始步长

	Node** record = new Node*[50];//记录u，v坐标，以及对应的z值；
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
				)//极大值点z值大于周围8个点的z值
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

	for (int i = 0; i < 50; ++i)//销毁动态分配的二维数组
	{
		delete[] record[i];
		record[i] = NULL;
	}
	delete[] record;
	record = NULL;
	return result;
}


//寻找曲面的极小值点
std::vector<PLib::Point2Dd> NurbsSurface_qj::GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV)//寻找曲面的极大值点
{
	double lengthU = maxU - minU;
	double lenghtV = maxV - minV;

	double deltaU = lengthU / 50; //u方向步长
	double deltaV = lenghtV / 50; //v方向步长

	Node** record = new Node*[50];//记录u，v坐标，以及对应的z值；
	for (int i = 0; i < 50; ++i)
	{
		record[i] = new Node[50];
	}

	double u = minU;
	double v = minV;
	Node * tempNode;
	//先将曲面粗略划分，确定极值点的大致范围，然后再在找到的范围内用二分法找极值点
	for (int i = 0; i < 50; ++i)//先将曲面划分，划分点的坐标存储起来，以方便确定极值点的范围
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
				)//极小值点的z值小于周围8个点的z值
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

	for (int i = 0; i < 50; ++i)//销毁动态分配的二维数组
	{
		delete[] record[i];
		record[i] = NULL;
	}
	delete[] record;
	record = NULL;
	return result;
}


//主程序：划分曲面，返回各个分块的内部点和边界点，等待边界处理
std::vector< Member> NurbsSurface_qj::surfaceAverageTotalSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> &Node_uv)
{
	
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极小值点
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极大值点
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//存储极值点

	int num = 0;//记录极值点数

	for (std::vector<PLib::Point2Dd>::iterator itr = pMinExtrmums.begin(); itr != pMinExtrmums.end(); ++itr)//将极小值加入极值点数组
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	
	for (std::vector<PLib::Point2Dd>::iterator itr = pMaxExtrmums.begin(); itr != pMaxExtrmums.end(); ++itr)//将极大值加入极值点数组
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	std::sort(extrmums, extrmums + num, compare_CPoint2D);//对极值点数组按v值从大到小排序，compare_CPoint2D是个函数，
	//std::sort自定义结构
	std::vector <PointStruct> uvEdge;//用于记录边界点的uv值
	//std::vector <PointStruct> vVec;//用于记录沿v方向的边界点，它们的标识符value记录u值大小
	//equalDivPoints divPoint;		//将极值点所在的纵横直线上的所有等分点存储下来,equalDivPoints有value（标识符）这个变量
	PointStruct divPoint;//用于临时记录极值点	
	divPoint.Id = -1;
	//Node AddSpace;//用于增加UVForOut
	std::vector<double> temp;//用于临时记录等分点
	//最外边界等分点
	//下边界
	temp = getUAverageLength(minU, maxU,minV, splitLengthU, 0);//从左往右
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
	if ((maxU - (uvEdge.end() - 1)->point.x()) < 0.001)//判断两点是否靠的太近
	{
		uvEdge.pop_back();
	}

	//右边界
	temp = getVAverageLength(minV, maxV, maxU, splitLengthV, 1);//从下往上	
	divPoint.point.x() = maxU;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.y() = temp[i];	
		uvEdge.push_back(divPoint);
	}
	if ((maxV - (uvEdge.end() - 1)->point.y()) < 0.001)//判断两点是否靠的太近
	{
		uvEdge.pop_back();
	}

	//上边界
	temp = getUAverageLength(minU,maxU,maxV, splitLengthU, 1);//从右往左	
	divPoint.point.y() = maxV;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.x() = temp[i];
		uvEdge.push_back(divPoint);
	}
	if (((uvEdge.end() - 1)->point.x()-minU ) < 0.001)//判断两点是否靠的太近
	{
		uvEdge.pop_back();
	}

	//左边界
	temp = getVAverageLength(minV, maxV, minU, splitLengthV, 2);//从上往下
	divPoint.point.x() =minU;
	for (int i = 0; i != temp.size(); ++i)
	{
		divPoint.point.y() = temp[i];
		uvEdge.push_back(divPoint);
	}
	if (((uvEdge.end() - 1)->point.y()-minV ) < 0.001)//判断两点是否靠的太近
	{
		uvEdge.pop_back();
	}

	for (int i = 0; i < num; ++i)//极值点的v值固定，将u方向的等分点全部记录下来，以方便划块后边界的处理；num表示极值点总数
	{
		divPoint.point.y() = extrmums[i].y;//记录极值点的v值
		temp = getUAverageLength(minU, extrmums[i].x, extrmums[i].y, splitLengthU, 1);//从右往左等分，对应type=1，2,记录等分点，是std::vector
		temp.erase(temp.begin());//除去极值点，因为下面还会添加一次
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.x() = temp[j];//将左半temp的u存成PointStruct，方便存入uvEdge中
			uvEdge.push_back(divPoint);//将等分边界点存入uvEdge
		}

		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());//将temp插入到divPoint中

		temp = getUAverageLength(extrmums[i].x, maxU, extrmums[i].y, splitLengthV, 0);//从左往右等分，对应type=0,3
		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.x() = temp[j];//将右半temp的u存成PointStruct，方便存入uvEdge中			
			uvEdge.push_back(divPoint);
		}
		//std::sort(uvEdge[i].points.begin(), uvEdge[i].points.end());//重新排序，u值从小到大
	}

	
	//此处compare在surface.cpp中有定义，因此改为compare1
	for (int i = 0; i < num; ++i)//极值点的u值固定，将v方向的等分点全部记录下来，以方便划块后边界的处理；num表示极值点总数
	{
		divPoint.point.x() = extrmums[i].x;//记录极值点的v值
		temp = getVAverageLength(minV, extrmums[i].y, extrmums[i].x, splitLengthU, 2);//从右往左等分，对应type=1，2,记录等分点，是std::vector
		temp.erase(temp.begin());//除去极值点，因为上面已经添加
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.y() = temp[j];//将下半temp的v存成PointStruct，方便存入vVec中
			uvEdge.push_back(divPoint);//将等分边界点存入vVec
		}

		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());//将temp插入到divPoint中

		temp = getVAverageLength(extrmums[i].y, maxV, extrmums[i].x, splitLengthV, 0);//从左往右等分，对应type=0,3
		temp.erase(temp.begin());//除去极值点，因为上面已经添加
		//(divPoint.points).insert(divPoint.points.end(), temp.begin(), temp.end());
		for (int j = 0; j != temp.size(); ++j)
		{
			divPoint.point.y() = temp[j];//将上半temp的v存成PointStruct，方便存入vVec中			
			uvEdge.push_back(divPoint);
		}

	}
	PLib::Point2Dd * extrmums2Dd = new PLib::Point2Dd[num];

	for (int i = 0; i < num; ++i)//将极值点的类型转换为PointDd存储
	{
		extrmums2Dd[i].x() = extrmums[i].x;
		extrmums2Dd[i].y() = extrmums[i].y;
	}
	//std::sort(extrmums2Dd, extrmums2Dd + num, compare_CPoint2D);//将极值点按照v从大到小排序，以方便划块后边界的处理；
	PBSTreeNode BSTree = NULL;
	BiTree_Create(&BSTree, extrmums2Dd, num);//根据极值点建立查找二叉树
	rectangleRegion * rec = new rectangleRegion[3 * num + 1];//num个极值点要划分3*num+1个区域
	PreOrder(BSTree, rec, 3 * num + 1, minU, maxU, minV, maxV);//按极值点将划块，划块的信息保存在rec数组中

	std::vector< Member> mVec;//用于记录总的各个分块内的member
	
	std::vector<Member> mVecTemp;//临时，用于记录循环中分块内的member，形参
	//std::vector<Node> nVecTemp;//临时，形参，Node

	for (int i =0; i < 3*num+1; ++i)//对所有区域进行等分，并将等分的点和杆信息分别存储在nVec和mVec中，3*num+1为分块数
		//***************************************************************************************************7待修改
	{
		mVecTemp = surfaceAverageSpit(rec[i], splitLengthU, splitLengthU, uvEdge);//将i块进行等分
		(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//将mVecTemp加到mVec中

		//mVec.push_back(mVecTemp);//这个会传递到main中
		//nVec.push_back(nVecTemp);//这个会传递到main中
		mVecTemp.clear();
		//nVecTemp.clear();
	}
	
	//边界进行杆件生成*********************************************
	std::vector<PLib::Point2Dd> extrmumVec;
	for (size_t i = 0; i < num; i++)
	{
		extrmumVec.push_back(extrmums2Dd[i]);
	}
	mVecTemp = EdgeTotal_Member(uvEdge, extrmumVec, num, minU, maxU, minV, maxV, splitLengthU, splitLengthV);
	(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//将mVecTemp加到mVec中
	//边界进行杆件生成*********************************************
	
	delete[] extrmums;
	delete[] rec;
	delete[] extrmums2Dd;
	
	Node_uv = UVForOut;
	return mVec;

}


//给定一点的u，v坐标，保持v不变，寻找u方向上弦长等于uChordLength的u值,type表示搜索方向
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
	case 0://从左下到右上
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r较大，需要减小
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}
			else if (uChordLength - r > error)//r较小，需要增加
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
	case 1://从右下到左上
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r较大，u需要增大
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}
			else if (uChordLength - r > error)//r较小，u需要减小
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
	case 2: //从右上到左下
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r较大，u需要增大
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				nextU = nextU + deltaU;
				flag = 1;
			}
			else if (uChordLength - r > error)//r较小，u需要减小
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
	case 3://从左上到右下
		do
		{

			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(nextU, v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - uChordLength > error)//r较大，需要减小
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				flag = -1;
			}
			else if (uChordLength - r > error)//r较小，需要增加
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


//给定一点的u，v坐标，寻找v方向上弦长等于vChordLength的u值,type表示搜索方向
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
	case 0://从左下往右上找点
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r较大，需要减小
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}
			else if (vChordLength - r > error)//r较小，需要增加
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
	case 1://从右下向左上寻找
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r较大，需要减小
			{
				if (flag == 1)
				{
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				flag = -1;
			}
			else if (vChordLength - r > error)//r较小，需要增加
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
	case 2://从右上到左下
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r较大，v需要增大
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}
			else if (vChordLength - r > error)//r较小，v需要减小
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

	case 3://从左上到右下
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(u, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - vChordLength > error)//r较大，v需要增大
			{
				if (flag == -1)
				{
					deltaV /= 2;
				}
				nextV = nextV + deltaV;
				flag = 1;
			}
			else if (vChordLength - r > error)//r较小，v需要减小
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


//给定两点，按照uChordLength和vChordLength求第三点（u，v），type表示搜索方向
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
	case 0://从左下往右上找点
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2较大，需要减小
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 > 0)//r2较小，需要增大
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1较大，需要减小
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 > 0)//r2较小，需要增大
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
	case 1://从右下向左上寻找
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2较大，需要减小
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 > 0)//r2较小，需要增大
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1较小，需要减小u使r1增大
				   {
					   if (Uflag == 1)//说明上一步需要增大u，而此点需要减小u，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 < 0)//r1较大，需要增大u使r1减小
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
	case 2://从右上到左下
	{

			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2较小，需要减小v使r2增大
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 < 0)//r2较大，需要增大v使r2减小
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1较小，需要减小u使r1增大
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 < 0)//r1较大，需要增大u是r1减小
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
	case 3://从左上到右下
	{
			   do
			   {
				   count++;
				   if (count >= 1e3)
				   {
					   return PLib::Point2Dd(-1, -1);;
				   }
				   tempNode3 = GetSurfacePointatUandV(u, v);
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2较小，需要减小v使r2增大
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 < 0)//r2较大，需要增大v使r2减小
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1较大，需要减小
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 > 0)//r2较小，需要增大
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
//extern std::vector<equalDivPoints> uvEdge;	//记录u方向均分点坐标
//extern std::vector<equalDivPoints> vVec;	//记录v方向均分点坐标



//分块内部等弦长划分，pDisplayVec用于存储要显示的点，传递到surfaceAverageTotalSpit中的nVecTemp，然后再传递到main中的nVec
//mDisplayVec用于存储要显示的杆，传递到surfaceAverageTotalSpit中的mVecTemp，然后再传递到main中的mVec
//返回mDisplayVec
std::vector<Member> NurbsSurface_qj::surfaceAverageSpit(rectangleRegion &rec, double splitLengthU, double splitLengthV, std::vector<PointStruct> &uvEdge)
{
	double maxU = rec.maxU + 0.3;
	//***********此处0.3有待优化***************************************************
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
	case 0://从左下往右上找点
		u = rec.minU;
		v = rec.minV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u <= maxU)//将第一列的坐标点找到
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
			nextV = nextVPoint(u, v, splitLengthV);//寻找每行的第一个点v坐标
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
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//如果1 2 点坐标正常，而找到的坐标为-1 -1，则画三角形
				{
					if (cosValue < -0.866)//根据1 2点坐标找第三点没有找到，若三点构成的三角形的角度大于150度，则根据三角形的左边两点找下一点
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

	case 1://从右下向左上寻找
		u = rec.maxU;
		v = rec.minV;

		points[0][0].point = PLib::Point2Dd(u, v);
		while (u >= minU)//将第一列的坐标点找到
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
			nextV = nextVPoint(u, v, splitLengthV, 1);//寻找每行的第一个点v坐标
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
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//如果1 2 点坐标正常，而找到的坐标为-1 -1，则画三角形
				{
					if (cosValue < -0.866)//根据1 2点坐标找第三点没有找到，若三点构成的三角形的角度大于150度，则根据三角形的左边两点找下一点
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

	case 2://从右上往左下找
		u = rec.maxU;
		v = rec.maxV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u >= minU)//将第一列的坐标点找到
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
			nextV = nextVPoint(u, v, splitLengthV, 2);//寻找每行的第一个点v坐标
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
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//如果1 2 点坐标正常，而找到的坐标为-1 -1，则画三角形
				{
					if (cosValue < -0.866)//根据1 2点坐标找第三点没有找到，若三点构成的三角形的角度大于150度，则根据三角形的左边两点找下一点
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
					//	&& points[row1][col1].x() != -1 && points[row1][col1].y() != -1))//画三角形仍找不到第三点，根据1点的坐标在其u方向上找下一点
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

	case 3://从左上往右下找
		u = rec.minU;
		v = rec.maxV;
		points[0][0].point = PLib::Point2Dd(u, v);
		while (u <= maxU)//将第一列的坐标点找到
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
			nextV = nextVPoint(u, v, splitLengthV, 3);//寻找每行的第一个点v坐标
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
					&& points[row0][col0 + 1].point.x() != -1 && points[row0][col0 + 1].point.y() != -1)//如果1 2 点坐标正常，而找到的坐标为-1 -1，则画三角形
				{
					if (cosValue < -0.866)//根据1 2点坐标找第三点没有找到，若三点构成的三角形的角度大于150度，则根据三角形的左边两点找下一点
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


	//Node ** nodeForDisplay = new Node *[1000];//存储显示的坐标点，记录u，v值，最终功能不需要，最后存在pDisplayVec
	////nodeForSave最后没用
	//Node ** nodeForSave = new Node *[1000];//存储坐标点，为了方便杆查找，点的ID顺序与点在容器中的顺序一致
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


	//unsigned long long int pID = 0;//点ID
	//unsigned long long int mID = 0;//杆ID
	//std::vector<Node> pSaveVec;//存储点的容器，点的顺序与ID完全一致，Node的x，y，z坐标表示实际值
	std::vector<Member> mDisplayVec; //信息流，存储杆的容器，杆的顺序与杆ID完全一致
	std::vector<Member> mSaveVec; //存储杆的容器，杆的顺序与杆ID完全一致
	Node addspace;//用来增加NodeForOut容器的大小
	//std::vector<double> uDirPoints;//记录u方向上不规则边界处的边界等分点
	//std::vector<double> vDirPoints;//记录v方向上不规则边界处的边界等分点
	//std::vector<PointStruct> uv_edge;//临时记录该分块的边界点,以便缩小收缩范围
	std::vector<PointStruct> uvedge_up, uvedge_down, uvedge_left, uvedge_right;//存储分块四边的边界点
	std::vector<PointStruct*> points_up, points_down,points_left,points_right;//存储最靠近不规则边界的内部点的指针
	std::vector<PointStruct*> point_corner;//记录角点


	switch (rec.type)
	{
	case 0://从左下到右上
	{			  
		//对该分块的边界点进行整体编号**********************************************************************************
		//将边界点的uv、编号都存在整体变量UVForOut中

		//返回分块的上下左右边界点的容器，并对边界点进行整体编号，将边界点加入到整体点
		Edge_udlr(uvEdge, rec,uvedge_up,uvedge_down,uvedge_left,uvedge_right);
		std::sort(uvedge_down.begin(), uvedge_down.end(), compare1_x);//对边界点从小到大排序
		std::sort(uvedge_up.begin(), uvedge_up.end(), compare1_x);
		std::sort(uvedge_left.begin(), uvedge_left.end(), compare3_y);
		std::sort(uvedge_right.begin(), uvedge_right.end(), compare3_y);
					
		//对分块内部的点进行整体编号,除离不规则边界最近的内部点,并返回与不规则边界（右上边界）最近的内部点容器
		for (int i = 1; i < 999; ++i)
		{
			for (int j = 1; j < 999; ++j)
			{
				if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				{
					if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() < rec.maxV))//*****************************需变动
						//将满足要求的区域内的点顺序存储	
					{
						//判断是否为不靠近不规则边界的点**********************************需改动
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
						//判断是否为靠近不规则上边界的点********************************需改动
						//与[i][j]相邻的两点都有可能是超越边界的,但是简单起见只考虑通常情况
						if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() >= rec.maxU))
						{
							point_corner.push_back ( &points[i][j]);//存入角点
							point_corner.push_back(&points[i][j-1]);//1
							point_corner.push_back(&points[i-1][j]);//2
							point_corner.push_back(&uvedge_right[uvedge_right.size() - 2]);//3
							point_corner.push_back(&uvedge_right.back());//4
							point_corner.push_back(&uvedge_up[uvedge_up.size() - 2]);//5
							point_corner.push_back(&uvedge_up.back());//6
						}
						else if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() < rec.maxU))
						{
							points_up.push_back(&points[i][j]);//上边界附近的内部点
						}
						else if ((points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() >= rec.maxU))
						{
							points_right.push_back(&points[i][j]);//右边界附近的内部点
						}
						//判断是否为靠近不规则右边界的点********************************需改动
						//与[i][j]相邻的两点都有可能是超越边界的
						
					}
				}				
			}
		}
		
		std::sort(points_up.begin(), points_up.end(), edge_compare1);//对point_up中的点按u从小到大排序
		std::sort(points_right.begin(), points_right.end(), edge_compare3);//对point_right中的点按v从小到大排序
		

		//判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
		for (int j = 0; j < 999; ++j)
		{
			if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
			{
				if (points[0][j].point.x() <=rec.maxU)//确保点在分块内，//*****************************需变动
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
		//判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点,此处i=1，表明左下角点在上面已经编号了
		for (int i = 1; i < 999; ++i)
		{
			if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
			{
				if (points[i][0].point.y() <= rec.maxV)//确保点在分块内，要边点，所以有=//*****************************需变动
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

		//对该分块的边界点进行整体编号***************************************************************************************
		//将点的关系用nId记录下来********************************************************************************************

		//double up_value = 2, down_value = 2;//记录边界相邻点与边界点的距离,2是初始值
		//int up_index = 0, down_index = 0;//记录边界相邻点与最近的边界点的编号,此编号为uv_edge的下标
		Member m;//用于记录杆信息
		std::vector<Member> edge_m;//记录不规则边界处的杆件

		//优先处理包含右上角点的边界，此时即为上边界，这样可以保证角点正常
		//上边界处理
		edge_m = Edge_Member(uvEdge,uvedge_up, points_up, 1);//x从左往右
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
		//右边界处理
		edge_m = Edge_Member(uvEdge, uvedge_right, points_right,  3);//y从下往上
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
		//角点处理
		edge_m = EdgeCorner_Member(point_corner);
		(mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
		
		//将下边界与内部点连接起来,此处i从1开始，抛去左下角点
		for (int j = 1; j < 999; ++j)
		{
			if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
			{
				if ((points[0][j].point.x() < rec.maxU) && (points[1][j].point.x()<rec.maxU))//确保点[0],[1]在分块内，//*****************************需变动
				{
					for (size_t k = 0; k < points[0][j].next.size(); ++k)
					{
						if (points[0][j].next[k].x == 1)
						{
							m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//两点相连形成杆件
							mDisplayVec.push_back(m);//杆件存入mDisplayVec						
						}
					}
				}
			}
		}

		//将左边界与内部点连接起来,此处i从1开始，抛去左下角点
		for (int i = 1; i < 999; ++i)
		{
			if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
			{
				if ((points[i][0].point.y() < rec.maxV)&&(points[i][1].point.y()<rec.maxV))//确保点[0],[1]在分块内，//*****************************需变动
				{
					for (size_t k = 0; k < points[i][0].next.size(); ++k)
					{
						if (points[i][0].next[k].y == 1)
						{
							m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//两点相连形成杆件
							mDisplayVec.push_back(m);//杆件存入mDisplayVec		
						}
					}
				}
			}
		}
		
		//处理内部点，不包含左下边界的点，不包含靠近右上边界的内部点（但其自身的联系需要处理），此处i,j初始为1
		for (int i = 1; i < 999; ++i)
		{
			for (int j = 1; j < 999; ++j)
			{
				if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				{
					if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() < rec.maxV))//判断点是否在这个分块中//*****************************需变动
					{
						//寻找与右边界最近的内部点						
						if (points[i][j + 1].point.x() >= rec.maxU)//判断j的u、v小于右边界，j+1的u大于等于右边界//*****************************需变动
						{							
							if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//判断点的上方点是否在分块内，判断有相连点
							{	
								m = Add_2Member(points[i][j], points[i+1][j]);//两点相连形成杆件
								mDisplayVec.push_back(m);//杆件存入mDisplayVec																						
							}
						}
						if (points[i + 1][j].point.y() >= rec.maxV)//寻找与上边界最近的内部点//*****************************需变动
						{							
							if ((points[i][j + 1].point.x() < rec.maxU) && (points[i][j + 1].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//判断点的右方点是否在分块内
							{
								m = Add_2Member(points[i][j], points[i][j + 1]);//两点相连形成杆件
								mDisplayVec.push_back(m);//杆件存入mDisplayVec								
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
	//将点的关系用nId记录下来********************************************************************************************



	case 1://从右下到左上
	{
			   //对该分块的边界点进行整体编号**********************************************************************************
			   //将边界点的uv、编号都存在整体变量UVForOut中
			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare2_x);//对边界点从大到小排序
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare2_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare3_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare3_y);

			   //对分块内部的点进行整体编号
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() >rec.minU) && (points[i][j].point.y() < rec.maxV))//*****************************需变动
							   //将满足要求的区域内的点顺序存储	
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
							   //判断是否为靠近不规则上边界的点********************************需改动
							   //与[i][j]相邻的两点都有可能是超越边界的,但是简单起见只考虑通常情况
							   //if (points[i + 1][j].point.y() >= rec.maxV)
							   //{
								  // points_up.push_back(&points[i][j]);
							   //}
							   ////判断是否为靠近不规则右边界的点********************************需改动
							   ////与[i][j]相邻的两点都有可能是超越边界的
							   //if (points[i][j + 1].point.x() <= rec.minU)
							   //{
								  // points_left.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   //point_corner = &points[i][j];//存入角点
								   point_corner.push_back(&points[i][j]);//存入角点
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_left[uvedge_left.size() - 2]);//3
								   point_corner.push_back(&uvedge_left.back());//4
								   point_corner.push_back(&uvedge_up[uvedge_up.size() - 2]);//5
								   point_corner.push_back(&uvedge_up.back());//6
							   }
							   else if ((points[i + 1][j].point.y() >= rec.maxV) && (points[i][j + 1].point.x() >rec.minU))
							   {
								   points_up.push_back(&points[i][j]);//上边界附近的内部点
							   }
							   else if ((points[i + 1][j].point.y() < rec.maxV) && (points[i][j + 1].point.x() < rec.minU))
							   {
								   points_left.push_back(&points[i][j]);//左边界附近的内部点
							   }
						   }
					   }
				   }
			   }
			   //对不规则边界旁的内部点进行排序（上左）
			   std::sort(points_up.begin(), points_up.end(),edge_compare2);//对point_up中的点按u从大到小排序
			   std::sort(points_left.begin(), points_left.end(), edge_compare3);//对point_right中的点按v从小到大排序
			   

			   //判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[0][j].point.x() >= rec.minU)//确保点在分块内，//*****************************需变动
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
			   //判断右边界上的点,将uvEdge中的整体编号赋予右边界上的点,此处i=1，表明左下角点在上面已经编号了
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[i][0].point.y() <= rec.maxV)//确保点在分块内//*****************************需变动
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
			   //对该分块的边界点进行整体编号***************************************************************************************



			   //将点的关系用nId记录下来********************************************************************************************

			  // double up_value = 2, down_value = 2;//记录边界相邻点与边界点的距离,2是初始值
			   //int up_index = 0, down_index = 0;//记录边界相邻点与最近的边界点的编号
			   Member m;//用于记录杆信息
				std::vector<Member> edge_m;//记录不规则边界处的杆件
			 
			   //优先处理包含左上角点的边界，此时即为左边界，这样可以保证角点正常
			   //左边界处理
				edge_m = Edge_Member(uvEdge, uvedge_left, points_left,3);//y从下往上
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //上边界处理
			   edge_m = Edge_Member(uvEdge, uvedge_up, points_up, 2);//x从右往左
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //角点处理
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
			 
			   //将下边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int j = 1; j < 999; ++j)
			   {				   
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if ((points[0][j].point.x() > rec.minU) && (points[1][j].point.x()>rec.minU))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //将右边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {

					   if ((points[i][0].point.y() < rec.maxV) && (points[i][1].point.y()<rec.maxV))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //处理内部点，不包含右下边界的点，此处i,j初始为1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() > rec.minU) && (points[i][j].point.y() < rec.maxV))//判断点是否在这个分块中//*****************************需变动
						   {
							   //寻找与左边界最近的内部点						
							   if (points[i][j + 1].point.x() <= rec.minU)//判断j的u、v小于左边界，j+1的u大于等于左边界//*****************************需变动
							   {
								   if ((points[i + 1][j].point.x() >rec.minU) && (points[i + 1][j].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//判断点的上方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() >= rec.maxV)//寻找与上边界最近的内部点//*****************************需变动
							   {								   
								   if ((points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y()<rec.maxV) && (points[i][j].next.size() != 0))//判断点的左方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i][j+1]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
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
		//将点的关系用nId记录下来********************************************************************************************			   

	case 2://从右上到左下
	{
			   //对该分块的边界点进行整体编号**********************************************************************************
			   //将边界点的uv、编号都存在整体变量UVForOut中

			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare2_x);//对边界点从小到大排序
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare2_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare4_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare4_y);

			   //对分块内部的点进行整体编号
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() >rec.minU) && (points[i][j].point.y() > rec.minV))//*****************************需变动
							   //将满足要求的区域内的点顺序存储	
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
							   //判断是否为靠近不规则上边界的点********************************需改动
							   //与[i][j]相邻的两点都有可能是超越边界的,但是简单起见只考虑通常情况
							   //if (points[i + 1][j].point.y() <= rec.minV)
							   //{
								  // points_down.push_back(&points[i][j]);
							   //}
							   ////判断是否为靠近不规则右边界的点********************************需改动
							   ////与[i][j]相邻的两点都有可能是超越边界的
							   //if (points[i][j + 1].point.x() <= rec.minU)
							   //{
								  // points_left.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   //point_corner = &points[i][j];//存入角点
								   point_corner.push_back(&points[i][j]);//存入角点
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_left[uvedge_left.size()-2]);//3
								   point_corner.push_back(&uvedge_left.back());//4
								   point_corner.push_back(&uvedge_down[uvedge_down.size()-2]);//5
								   point_corner.push_back(&uvedge_down.back());//6
							   }
							   else if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() > rec.minU))
							   {
								   points_down.push_back(&points[i][j]);//下边界附近的内部点
							   }
							   else if ((points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() <= rec.minU))
							   {
								   points_left.push_back(&points[i][j]);//左边界附近的内部点
							   }
						   }
					   }
				   }
			   }
			   //对不规则边界旁的内部点进行排序（左下）
			   std::sort(points_down.begin(), points_down.end(), edge_compare2);//对point_down中的点按u从大到小排序
			   std::sort(points_left.begin(), points_left.end(), edge_compare4);//对point_left中的点按v从大到下排序
			  

			   //判断上边界上的点,将uvEdge中的整体编号赋予下边界上的点
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[0][j].point.x() >= rec.minU)//确保点在分块内，要边点，所以有=//*****************************需变动
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
			   //判断右边界上的点,将uvEdge中的整体编号赋予左边界上的点
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[i][0].point.y() >= rec.minV)//确保点在分块内，要边点，所以有=//*****************************需变动
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
			   //对该分块的边界点进行整体编号***************************************************************************************



			   //将点的关系用nId记录下来********************************************************************************************

			   //double up_value = 2, down_value = 2;//记录边界相邻点与边界点的距离,2是初始值
			   //int up_index = 0, down_index = 0;//记录边界相邻点与最近的边界点的编号
			   Member m;//用于记录杆信息
			   std::vector<Member> edge_m;//记录不规则边界处的杆件

			   //优先处理包含左下角点的边界，此时即为下边界，这样可以保证角点正常
			   //下边界处理
			   edge_m = Edge_Member(uvEdge, uvedge_down, points_down, 2);//x从右往左
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //左边界处理
			   edge_m = Edge_Member(uvEdge, uvedge_left, points_left, 4);//y从上往下
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //角点处理
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());
			   

			   //将下边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int j = 1; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if ((points[0][j].point.x() > rec.minU) && (points[1][j].point.x()>rec.minU))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //将右边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if ((points[i][0].point.y() > rec.minV) && (points[i][1].point.y()> rec.minV))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //处理内部点，不包含右下边界的点，此处i,j初始为1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() > rec.minU) && (points[i][j].point.y() > rec.minV))//判断点是否在这个分块中//*****************************需变动
						   {
							   //寻找与左边界最近的内部点						
							   if (points[i][j + 1].point.x() <= rec.minU)//判断j的u、v小于左边界，j+1的u大于等于左边界//*****************************需变动
							   {								  
								   if ((points[i + 1][j].point.x() >rec.minU) && (points[i + 1][j].point.y()>rec.minV) && (points[i][j].next.size() != 0))//判断点的下方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() <= rec.minV)//寻找与上边界最近的内部点//*****************************需变动
							   {								   
								   if ((points[i][j + 1].point.x() >rec.minU) && (points[i][j + 1].point.y()>rec.minV) && (points[i][j].next.size() != 0))//判断点的左方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i][j+1]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
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

	case 3://从左上到右下
	{
			   //对该分块的边界点进行整体编号**********************************************************************************
			   //将边界点的uv、编号都存在整体变量UVForOut中
			   Edge_udlr(uvEdge, rec, uvedge_up, uvedge_down, uvedge_left, uvedge_right);
			   std::sort(uvedge_down.begin(), uvedge_down.end(), compare1_x);//对边界点从小到大排序
			   std::sort(uvedge_up.begin(), uvedge_up.end(), compare1_x);
			   std::sort(uvedge_left.begin(), uvedge_left.end(), compare4_y);
			   std::sort(uvedge_right.begin(), uvedge_right.end(), compare4_y);

			   //对分块内部的点进行整体编号
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() <rec.maxU) && (points[i][j].point.y() > rec.minV))//*****************************需变动
							   //将满足要求的区域内的点顺序存储	
						   {
							   //判断是否为不靠近不规则边界的点**********************************需改动
							   if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() <rec.maxU) && (points[i][j + 1].point.y() >rec.minV))
							   {
								   points[i][j].Id = nId;
								   addspace.setX(points[i][j].point.x());
								   addspace.setY(points[i][j].point.y());
								   addspace.setID(nId);
								   UVForOut.push_back(addspace);
								   ++nId;
							   }
							   //判断是否为靠近不规则上边界的点********************************需改动
							   //与[i][j]相邻的两点都有可能是超越边界的,但是简单起见只考虑通常情况
							   //if (points[i + 1][j].point.y() <= rec.minV)
							   //{
								  // points_down.push_back(&points[i][j]);
							   //}
							   ////判断是否为靠近不规则右边界的点********************************需改动
							   ////与[i][j]相邻的两点都有可能是超越边界的
							   //if (points[i][j + 1].point.x() >= rec.maxU)
							   //{
								  // points_right.push_back(&points[i][j]);
							   //}
							   if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() >= rec.maxU))
							   {
								  // point_corner = &points[i][j];//存入角点
								   point_corner.push_back(&points[i][j]);//存入角点
								   point_corner.push_back(&points[i][j - 1]);//1
								   point_corner.push_back(&points[i - 1][j]);//2
								   point_corner.push_back(&uvedge_right[uvedge_right.size()-2]);//3
								   point_corner.push_back(&uvedge_right.back());//4
								   point_corner.push_back(&uvedge_down[uvedge_down.size()-2]);//5
								   point_corner.push_back(&uvedge_down.back());//6
							   }
							   else if ((points[i + 1][j].point.y() <= rec.minV) && (points[i][j + 1].point.x() < rec.maxU))
							   {
								   points_down.push_back(&points[i][j]);//下边界附近的内部点
							   }
							   else if ((points[i + 1][j].point.y() > rec.minV) && (points[i][j + 1].point.x() >= rec.maxU))
							   {
								   points_right.push_back(&points[i][j]);//右边界附近的内部点
							   }
						   }
					   }
				   }
			   }
			   //对不规则边界旁的内部点进行排序（右下）
			   std::sort(points_down.begin(), points_down.end(), edge_compare1);//对point_down中的点按u从小到大排序
			   std::sort(points_right.begin(), points_right.end(), edge_compare4);//对point_right中的点按v从大到小排序
			  

			   //判断上边界上的点,将uvEdge中的整体编号赋予下边界上的点
			   for (int j = 0; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[0][j].point.x() <= rec.maxU)//确保点在分块内，要边点，所以有=//*****************************需变动
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
			   //判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点
			   for (int i = 0; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if (points[i][0].point.y() >= rec.minV)//确保点在分块内，要边点，所以有=//*****************************需变动
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
			   //对该分块的边界点进行整体编号***************************************************************************************
			   //将点的关系用nId记录下来********************************************************************************************

			   //double up_value = 2, down_value = 2;//记录边界相邻点与边界点的距离,2是初始值
			  // int up_index = 0, down_index = 0;//记录边界相邻点与最近的边界点的编号
			   Member m;//用于记录杆信息
			   std::vector<Member> edge_m;//记录不规则边界处的杆件

			   //优先处理包含右下角点的边界，此时即为右边界，这样可以保证角点正常
			   //右边界处理
			   edge_m = Edge_Member(uvEdge, uvedge_right, points_right, 4);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //下边界处理
			   edge_m = Edge_Member(uvEdge, uvedge_down, points_down, 1);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());//存入mDispalyVec中
			   //角点处理
			   edge_m = EdgeCorner_Member(point_corner);
			   (mDisplayVec).insert(mDisplayVec.end(), edge_m.begin(), edge_m.end());

			   //将下边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int j = 1; j < 999; ++j)
			   {
				   if ((points[0][j].point.x() != -1 || points[0][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if ((points[0][j].point.x() < rec.maxU) && (points[1][j].point.x()<rec.maxU))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[0][j].next.size(); ++k)
						   {
							   if (points[0][j].next[k].x == 1)
							   {
								   m = Add_2Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //将右边界与内部点连接起来,此处i从1开始，抛去左下角点
			   for (int i = 1; i < 999; ++i)
			   {
				   if ((points[i][0].point.x() != -1 || points[i][0].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
				   {
					   if ((points[i][0].point.y() > rec.minV) && (points[i][1].point.y()>rec.minV))//确保点在分块内，不要边点，所以没有=//*****************************需变动
					   {
						   for (size_t k = 0; k < points[i][0].next.size(); ++k)
						   {
							   if (points[i][0].next[k].y == 1)
							   {
								   m = Add_2Member(points[i][0], points[points[i][0].next[k].x][points[i][0].next[k].y]);//连杆
								   mDisplayVec.push_back(m);
							   }
						   }
					   }
				   }
			   }

			   //处理内部点，不包含右下边界的点，此处i,j初始为1
			   for (int i = 1; i < 999; ++i)
			   {
				   for (int j = 1; j < 999; ++j)
				   {
					   if ((points[i][j].point.x() != -1 || points[i][j].point.y() != -1))//判断点points[i][j]是否被赋予x，y值
					   {
						   if ((points[i][j].point.x() < rec.maxU) && (points[i][j].point.y() > rec.minV))//判断点是否在这个分块中//*****************************需变动
						   {
							   //寻找与左边界最近的内部点						
							   if (points[i][j + 1].point.x() >= rec.maxU)//判断j的u、v小于左边界，j+1的u大于等于左边界//*****************************需变动
							   {								  
								   if ((points[i + 1][j].point.x() < rec.maxU) && (points[i + 1][j].point.y()>rec.minV) && (points[i][j].next.size() != 0))//判断点的下方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i + 1][j]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
								   }
							   }
							   else if (points[i + 1][j].point.y() <= rec.minV)//寻找与上边界最近的内部点//*****************************需变动
							   {								   
								   if ((points[i][j + 1].point.x() < rec.maxU) && (points[i][j + 1].point.y()>rec.minV) && (points[i][j].next.size() != 0))//判断点的右方点是否在分块内
								   {
									   m = Add_2Member(points[i][j], points[i ][j+1]);//两点相连形成杆件
									   mDisplayVec.push_back(m);//杆件存入mDisplayVec						
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


//添加杆件，point_in第一点（内部点），point_edge第二点（不规则边界点）
Member NurbsSurface_qj::Add_Member(PointStruct point_in, PointStruct point_edge)
{
	Member m;
	Node *tempNode1;
	Node *tempNode2;
	m.setInode(point_in.Id);//写入上点
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

//添加杆件，point_in第一点（内部点），point_edge第二点（规则边界）
Member NurbsSurface_qj::Add_2Member(PointStruct point_in, PointStruct point_edge)
{
	Member m;
	Node *tempNode1;
	Node *tempNode2;
	m.setInode(point_in.Id);//写入上点
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

//判断靠近边界的点界于边界点的区间，返回up_index,down_index********(边界点，方向值，点值，xy方向，上点，下点）
void NurbsSurface_qj::BetweenSearch(std::vector<PointStruct> uv_edge, double key1, double key2, int xy, int &up_index, int &down_index)
{
	double up_value = 2, down_value = 2;
	up_index = -1;
	down_index = -1;
	switch (xy)
	{
	case 0://表示寻找y方向的区间
	{
			   for (int k = 0; k != uv_edge.size(); ++k)
			   {
				   if (uv_edge[k].point.x() == key1)//筛选出某边界上的点
				   {
					   if (uv_edge[k].point.y() - key2 >= 0)//此点位于上方
					   {
						   if (uv_edge[k].point.y() - key2 <= up_value)//比上个差值小
						   {
							   up_value = uv_edge[k].point.y() - key2;
							   up_index = k;//记录点的std::vector编号
						   }
					   }
					   else//此点位于下方
					   {
						   if (key2 - uv_edge[k].point.y() <= down_value)//比上个差值小
						   {
							   down_value =  key2 - uv_edge[k].point.y();
							   down_index = k;// 记录点的std::vector编号
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
				  if (uv_edge[k].point.y() == key1)//筛选出某边界上的点
				  {
					  if (uv_edge[k].point.x() - key2 >= 0)//此点位于上方
					  {
						  if (uv_edge[k].point.x() - key2 <= up_value)//比上个差值小
						  {
							  up_value = uv_edge[k].point.x() - key2;
							  up_index = k;//记录点的std::vector编号
						  }
					  }
					  else//此点位于下方
					  {
						  if (key2-uv_edge[k].point.x()  <= down_value)//比上个差值小
						  {
							  down_value =  key2- uv_edge[k].point.x() ;
							  down_index = k;// 记录点的std::vector编号
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

//返回分块的上下左右边界点的容器，并对边界点进行整体编号，将边界点加入到整体点中，&用于修改uvEdge的Id
void NurbsSurface_qj::Edge_udlr(std::vector<PointStruct> &uvEdge, rectangleRegion rec, std::vector<PointStruct> &uvedge_up, std::vector<PointStruct> &uvedge_down, std::vector<PointStruct> &uvedge_left, std::vector<PointStruct> &uvedge_right)
{
	//std::vector<PointStruct> uv_edge;
	//int indexld = -1, indexrd = -1, indexlu = -1, indexru = -1;//标记分部四个角点是否已经进行整体标号，按照leftdown，rightdown，leftup，rightup		
	Node addspace;//用来增加NodeForOut容器的大小
	//PointStruct indexAdd;//临时，用于增加角点用来增加uvEdge的大小，用来存放本身不在其中的分部角点
	//indexAdd.Id = -1;//Id初始化
	for (size_t i = 0; i != uvEdge.size(); ++i)//取遍uvEdge[i]
	{
		//1判断uvEdge[i]的v值是否与上边界相等,且是否介于该部分,包含右角点,不包含左角点
		if ((uvEdge[i].point.y() == rec.maxV) && (uvEdge[i].point.x()>rec.minU) && (uvEdge[i].point.x()<=rec.maxU))
		{
			if (uvEdge[i].Id == -1)//判断该点是否已经进行了编号，当没有时才进行以下操作
			{
				uvEdge[i].Id = nId;//将边界点的id赋值为整体nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//将uvEdge的值赋给addspace
				UVForOut.push_back(addspace);//将边界点加入到整体中
				++nId;//整体坐标增加
			}
			uvedge_up.push_back(uvEdge[i]);//将该分块上边界的点存入uvedge_up
			continue;
		}

		//2判断uvEdge[i]的v值是否与下边界相等,且是否介于该部分，包含左角点，不包含右角点
		if ((uvEdge[i].point.y() == rec.minV) && (uvEdge[i].point.x()>=rec.minU) && (uvEdge[i].point.x()<rec.maxU))
		{
			if (uvEdge[i].Id == -1)//判断该点是否已经进行了编号，当没有时才进行以下操作
			{
				uvEdge[i].Id = nId;//将边界点的id赋值为整体nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//将uvEdge的值赋给addspace
				UVForOut.push_back(addspace);//将边界点加入到整体中
				++nId;//整体坐标增加
			}
			uvedge_down.push_back(uvEdge[i]);//将该分块上边界的点存入uvedge_down
			continue;
		}

		//3判断uvEdge[i]的u值是否与左边界相等,且是否介于该部分，包含上角点,不包含下角点
		if ((uvEdge[i].point.x() == rec.minU) && (uvEdge[i].point.y()>rec.minV) && (uvEdge[i].point.y()<=rec.maxV))
		{
			if (uvEdge[i].Id == -1)//判断该点是否已经进行了编号，当没有时才进行以下操作
			{
				uvEdge[i].Id = nId;//将边界点的id赋值为整体nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//将uvEdge的值赋给addspace
				UVForOut.push_back(addspace);//将边界点加入到整体中
				++nId;//整体坐标增加
			}
			uvedge_left.push_back(uvEdge[i]);//将该分块上边界的点存入uvedge_left
			continue;
		}

		//4判断uvEdge[i]的u值是否与右边界相等,且是否介于该部分，包含下角点，不包含上角点
		if ((uvEdge[i].point.x() == rec.maxU) && (uvEdge[i].point.y()>=rec.minV) && (uvEdge[i].point.y() < rec.maxV))
		{
			if (uvEdge[i].Id == -1)//判断该点是否已经进行了编号，当没有时才进行以下操作
			{
				uvEdge[i].Id = nId;//将边界点的id赋值为整体nId
				addspace.setX(uvEdge[i].point.x());
				addspace.setY(uvEdge[i].point.y());
				addspace.setID(nId);//将uvEdge的值赋给addspace
				UVForOut.push_back(addspace);//将边界点加入到整体中
				++nId;//整体坐标增加
			}
			uvedge_right.push_back(uvEdge[i]);//将该分块上边界的点存入uvedge_right
		}
	}
	
	
}


//边界杆件处理(全部边界点容器用于增加边界点【当内部点位于两边界点中间且与边界很近时】，某一边界的点，靠近边界的内部点,uv方向）
//只考虑在一个边界区间内至多有两个内部点，不考虑三个内部点的情况（此情况基本不可能）
std::vector<Member> NurbsSurface_qj::Edge_Member(std::vector<PointStruct> &uvEdge, std::vector<PointStruct> uvedge, std::vector<PointStruct*> pointsedge, int xy, double length_differ)
{
	Member m;//用于增加edge_m
	std::vector<Member> edge_m;//存储边界杆件
	Member point_m;//记录内部点的上下边界点的整体编号
	std::vector<Member> edge_point_m;//存储point_m
	int up_index = -1;//上下点的ID	
	int max_or_min = 0;//标记point_corner是极大还是极小值
	PointStruct addpoint;
	Node addspace;
	
	switch (xy)
	{
	case 1://u方向从左到右
	{
			  
			   for (size_t i = 0; i < pointsedge.size(); i++)//寻找内部点的上下边界点
			   {	
				   up_index = -1;				   
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.x() - pointsedge[i]->point.x() > 0)//u大于此点的第一个边界点
					   {
						   up_index = j;//将右边界点的下标传出
						   break;
					   }
					}

				   if (up_index >= 1 && up_index <= uvedge.size())//内部点有左、右边界点与之相连，若有则判断两点之间的距离，此处设定为1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());					   
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//右点坐标
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index-1].point.x(), uvedge[up_index-1].point.y());//左点坐标
					   tempNode4 = GetSurfacePointatUandV(pointsedge[i]->point.x(), uvedge[up_index].point.y());//投影坐标
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//与右点杆件长度若小于2m（默认）
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x()) && (pointsedge[i - 1]->point.x() < uvedge[up_index].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id) )//判断上个内部点的区间是否是【index-1，index】，==表示上个点被归并
							   {
								   pointsedge[i-1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//与左点杆件长度若小于2m（默认）
					   {
						   pointsedge[i]->Id = uvedge[up_index-1].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index-1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index-1].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.x() > uvedge[up_index - 2].point.x()) && (pointsedge[i - 1]->point.x() < uvedge[up_index-1].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//判断上个内部点的区间是否为【index-2，index-1】，==表示上个点被归并
							   {
								   pointsedge[i-1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   //投影点距离小于2m（默认），且u靠近中间
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) && 
						   (fabs(pointsedge[i]->point.x() - (uvedge[up_index - 1].point.x() + uvedge[up_index].point.x()) / 2)<(uvedge[up_index ].point.x()- uvedge[up_index - 1].point.x()) / 5))
					   {//此处不考虑三个内部点的情况
						   pointsedge[i]->Id = nId;//将该点升级为边界点
						   /*pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   addpoint.Id = nId;
						   addpoint.point.y() = uvedge[up_index].point.y();
						   addpoint.point.x() = pointsedge[i]->point.x();
						   uvEdge.push_back(addpoint);//添加到全局边界容器
						   //uvedge.insert(uvedge.begin() + up_index, addpoint);//添加到局部边界
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   //一般情况
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index-1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;					   
				   }
				   else if (up_index == 0)//内部点只有右边界点与之相连
				   {
					   Node *tempNode1;					   
					   Node *tempNode2;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());					   
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//右点坐标					   
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//下点杆件长度若小于1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i+1]->point.x()< uvedge[up_index+1].point.x())//如果与下一个内部点的区间为【0,1】，
						   {
							   pointsedge[i]->next.clear();//将该点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
						   delete tempNode1;						  
						   delete tempNode2;
						   continue;//该点已处理完，处理下个i
					   }
					   //大于1时
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(-1);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;

					   delete tempNode1;					  
					   delete tempNode2;					 
				   }
				   else if (up_index == -1 )//内部点有左边界点与之相连
				   {
					   Node *tempNode1;
					   Node *tempNode3;	
					   up_index = uvedge.size() - 1;//使up_index指向边界左点的下标
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//左点坐标					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//上点杆件长度若小于1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//如果与上一个内部点的区间【index-1，index】，==表示上个点被归并
						   {
							   pointsedge[i-1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
						   delete tempNode1;
						   delete tempNode3;						  
						   continue;//该点已处理完，处理下个i
					   }					  
					   //大于1时
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;					   
				   }
			   }
			   //判断是否有两个内部点和相同的两个边界点相连
			   for (size_t i = 0; i < edge_point_m.size()-1; i++)
			   {
				   for (size_t j = i+1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//如果有相同的则改成四边形相连
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }

			   //建立内部点与边界点相连的杆件
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

	case 2://u方向从右到左
	{

			   for (size_t i = 0; i < pointsedge.size(); i++)//寻找内部点的上下边界点
			   {
				   up_index = -1;
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.x() - pointsedge[i]->point.x() < 0)//u大于此点的第一个边界*************************
					   {
						   up_index = j;//将右边界点的下标传出
						   break;
					   }
				   }

				   if (up_index >= 1 && up_index <= uvedge.size())//内部点有左、右边界点与之相连，若有则判断两点之间的距离，此处设定为1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//右点坐标
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//左点坐标
					   tempNode4 = GetSurfacePointatUandV(pointsedge[i]->point.x(), uvedge[up_index].point.y());//投影坐标
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//与右点杆件长度若小于2m（默认）
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.x() < uvedge[up_index - 1].point.x()) && (pointsedge[i - 1]->point.x() > uvedge[up_index].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//判断上个内部点的区间是否是【index-1，index】，==表示上个点被归并
							   {
								   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//与左点杆件长度若小于2m（默认）
					   {
						   pointsedge[i]->Id = uvedge[up_index - 1].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index-1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index-1].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.x()<uvedge[up_index - 2].point.x()) && (pointsedge[i - 1]->point.x() > uvedge[up_index - 1].point.x())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//判断上个内部点的区间是否为【index-2，index-1】，==表示上个点被归并
							   {
								   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   //投影点距离小于2m（默认），且u靠近中间
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						   (fabs(pointsedge[i]->point.x() - (uvedge[up_index - 1].point.x() + uvedge[up_index].point.x()) / 2)<(uvedge[up_index].point.x() - uvedge[up_index - 1].point.x()) / 5))
					   {//此处不考虑三个内部点的情况
						   pointsedge[i]->Id = nId;//将该点升级为边界点
						   /*pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   addpoint.Id = nId;
						   addpoint.point.y() = uvedge[up_index].point.y();
						   addpoint.point.x() = pointsedge[i]->point.x();
						   uvEdge.push_back(addpoint);//添加到全局边界容器
						   //uvedge.insert(uvedge.begin() + up_index, addpoint);//添加到局部边界
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   //一般情况
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index - 1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;
				   }
				   else if (up_index == 0)//内部点只有右边界点与之相连
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//右点坐标					   
					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//下点杆件长度若小于1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i + 1]->point.x()> uvedge[up_index + 1].point.x())//如果与下一个内部点的区间为【0,1】，
						   {
							   pointsedge[i]->next.clear();//将该点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
						   delete tempNode1;
						   delete tempNode2;
						   continue;//该点已处理完，处理下个i
					   }
					   //大于1时
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(-1);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;

					   delete tempNode1;
					   delete tempNode2;
				   }
				   else if (up_index == -1)//内部点有左边界点与之相连*******最外点
				   {
					   Node *tempNode1;
					   Node *tempNode3;
					   up_index = uvedge.size() - 1;//使up_index指向边界左点的下标
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//左点坐标					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//上点杆件长度若小于1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.x()< uvedge[up_index - 1].point.x() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//如果与上一个内部点的区间【index-1，index】，==表示上个点被归并
						   {//*****************此处》改为《
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
						   delete tempNode1;
						   delete tempNode3;
						   continue;//该点已处理完，处理下个i
					   }
					   //大于1时
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;
				   }
			   }
			   //判断是否有两个内部点和相同的两个边界点相连
			   for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			   {
				   for (size_t j = i + 1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//如果有相同的则改成四边形相连
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }

			   //建立内部点与边界点相连的杆件
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

	case 3://v方向从下往上
	{
			  
			  for (size_t i = 0; i < pointsedge.size(); i++)//寻找内部点的上下边界点
			  {
				  up_index = -1;				  
				  for (size_t j = 0; j < uvedge.size(); j++)
				  {
					  if (uvedge[j].point.y() - pointsedge[i]->point.y() > 0)//v大于此点的第一个边界点
					  {
						  up_index = j;//将上边界点的下标传出
						  break;
					  }					 
				  }

				  if (up_index >= 1 && up_index <= uvedge.size())//判断内部点是否有上、下边界点与之相连，若有则判断两点之间的距离，此处设定为1m
				  {
					  Node *tempNode1;
					  Node *tempNode2;
					  Node *tempNode3;
					  Node *tempNode4;
					  tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					  tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//上点坐标
					  tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//下点坐标
					  tempNode4 = GetSurfacePointatUandV(uvedge[up_index].point.x(), pointsedge[i]->point.y());//投影坐标
					  
					  if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//上点杆件长度若小于2m
					  {
						  pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						 /* pointsedge[i]->point.x() = uvedge[up_index].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						  if (i != 0)//保证不为第一个内部点
						  {
							  if ((pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y()) && (pointsedge[i - 1]->point.y() < uvedge[up_index].point.y())
								  || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id) )//判断上个内部点的区间是否是【index-1，index】，==表示上个点被归并
							  {
								  pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							  }
						  }
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//该点已处理完，处理下个i
					  }
					  if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//下点杆件长度若小于2m
					  {
						  pointsedge[i]->Id = uvedge[up_index - 1].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						  /*pointsedge[i]->point.x() = uvedge[up_index - 1].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index - 1].point.y();*/
						  if (i != 0)//保证不为第一个内部点
						  {
							  if ((pointsedge[i - 1]->point.y() > uvedge[up_index - 2].point.y()) && (pointsedge[i - 1]->point.y() < uvedge[up_index-1].point.y())
								  || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id) )//判断上个内部点的区间是否是【index-2，index-1】，==表示上个点被归并
							  {
								  pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							  }
						  }
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//该点已处理完，处理下个i
					  }
					  if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						  (fabs(pointsedge[i]->point.y() - (uvedge[up_index - 1].point.y() + uvedge[up_index].point.y()) / 2)<(uvedge[up_index].point.y() - uvedge[up_index - 1].point.y()) / 5))
					  {//此处不考虑三个内部点的情况
						  pointsedge[i]->Id = nId;//将该点升级为边界点
						 /* pointsedge[i]->point.x() = uvedge[up_index].point.x();*/
						  addpoint.Id = nId;
						  addpoint.point.x() = uvedge[up_index].point.x();
						  addpoint.point.y() = pointsedge[i]->point.y();
						  uvEdge.push_back(*pointsedge[i]);//添加到全局边界容器
						  uvedge.insert(uvedge.begin() + up_index, addpoint);//添加到局部边界
						  addspace.setX(addpoint.point.x());
						  addspace.setY(addpoint.point.y());
						  addspace.setID(nId);
						  UVForOut.push_back(addspace);
						  nId++;
						  delete tempNode1;
						  delete tempNode2;
						  delete tempNode3;
						  delete tempNode4;
						  continue;//该点已处理完，处理下个i
					  }
					  //一般情况
					  pointsedge[i]->Id = nId;//对该点进行整体编号
					  point_m.set_ID(nId);
					  point_m.setInode(uvedge[up_index - 1].Id);
					  point_m.setJnode(uvedge[up_index].Id);					 
					  edge_point_m.push_back(point_m);//将这种相连关系记录下来
					  addspace.setX(pointsedge[i]->point.x());
					  addspace.setY(pointsedge[i]->point.y());
					  addspace.setID(nId);
					  UVForOut.push_back(addspace);
					  nId++;
					  if (i != 0)//保证不为第一个内部点
					  {
						  if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//判断上个内部点是否被归并
						  {
							  pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						  }
					  }
					  delete tempNode1;
					  delete tempNode2;
					  delete tempNode3;
					  delete tempNode4;
				  }
				  else if (up_index == 0)//内部点只有上边界点与之相连
				  {					  
						Node *tempNode1;
						Node *tempNode2;
						tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
						tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//下点坐标					   
						if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//下点杆件长度若小于1m
						{
							pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
							/*pointsedge[i]->point.x() = uvedge[up_index].point.x();
							pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

							if (pointsedge[i + 1]->point.y() < uvedge[up_index + 1].point.y())//如果下一个内部点的区间为【0,1】，
							{
								pointsedge[i]->next.clear();//将该点的内部点自身的联系删除，这样可以保证没有重复杆件
							}
							delete tempNode1;
							delete tempNode2;
							continue;//该点已处理完，处理下个i
						}
						//都大于1时
						pointsedge[i]->Id = nId;//对该点进行整体编号
						point_m.set_ID(nId);
						point_m.setInode(-1);
						point_m.setJnode(uvedge[up_index].Id);
						edge_point_m.push_back(point_m);//将这种相连关系记录下来
						addspace.setX(pointsedge[i]->point.x());
						addspace.setY(pointsedge[i]->point.y());
						addspace.setID(nId);
						UVForOut.push_back(addspace);
						nId++;
						delete tempNode1;
						delete tempNode2;					 
				  }
				  else if (up_index == -1)//内部点有下边界点与之相连
				  {
					  Node *tempNode1;
					  Node *tempNode3;
					  up_index = uvedge.size() - 1;//使up_index指向边界左点的下标
					  tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					  tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//上点坐标					   
					  if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//上点杆件长度若小于1m
					  {
						  pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						  /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						  pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						  if (pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//如果上一个内部点的区间为【index-1,index】，，，==表示上个点被归并
						  {
							  pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						  }
						  delete tempNode1;
						  delete tempNode3;
						  continue;//该点已处理完，处理下个i
					  }
					  //大于1时
					  pointsedge[i]->Id = nId;//对该点进行整体编号
					  point_m.set_ID(nId);
					  point_m.setInode(uvedge[up_index].Id);
					  point_m.setJnode(-1);
					  edge_point_m.push_back(point_m);//将这种相连关系记录下来
					  addspace.setX(pointsedge[i]->point.x());
					  addspace.setY(pointsedge[i]->point.y());
					  addspace.setID(nId);
					  UVForOut.push_back(addspace);
					  nId++;
					  if (i != 0)//保证不为第一个内部点
					  {
						  if (pointsedge[i - 1]->Id == uvedge[up_index ].Id)//判断上个内部点是否被归并
						  {
							  pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						  }
					  }
					  delete tempNode1;
					  delete tempNode3;
				  }
			  }
			  //判断是否有两个内部点和相同的两个边界点相连
			  for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			  {
				  for (size_t j = i + 1; j < edge_point_m.size(); j++)
				  {
					  if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					  {//如果有相同的则改成四边形相连
						  edge_point_m[i].setJnode(-1);
						  edge_point_m[j].setInode(-1);
					  }
				  }
			  }
			  //建立杆件
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

	case 4://v方向从上往下
	{

			   for (size_t i = 0; i < pointsedge.size(); i++)//寻找内部点的上下边界点
			   {
				   up_index = -1;
				   for (size_t j = 0; j < uvedge.size(); j++)
				   {
					   if (uvedge[j].point.y() - pointsedge[i]->point.y() < 0)//v大于此点的第一个边界点
					   {
						   up_index = j;//将上边界点的下标传出
						   break;
					   }
				   }

				   if (up_index >= 1 && up_index <= uvedge.size())//判断内部点是否有上、下边界点与之相连，若有则判断两点之间的距离，此处设定为1m
				   {
					   Node *tempNode1;
					   Node *tempNode2;
					   Node *tempNode3;
					   Node *tempNode4;
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//上点坐标
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index - 1].point.x(), uvedge[up_index - 1].point.y());//下点坐标
					   tempNode4 = GetSurfacePointatUandV(uvedge[up_index].point.x(), pointsedge[i]->point.y());//投影坐标

					   if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//上点杆件长度若小于2m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /* pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.y() < uvedge[up_index - 1].point.y()) && (pointsedge[i - 1]->point.y() > uvedge[up_index].point.y())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//判断上个内部点的区间是否是【index-1，index】，==表示上个点被归并
							   {
								   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//下点杆件长度若小于2m
					   {
						   pointsedge[i]->Id = uvedge[up_index - 1].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index - 1].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index - 1].point.y();*/
						   if (i != 0)//保证不为第一个内部点
						   {
							   if ((pointsedge[i - 1]->point.y() < uvedge[up_index - 2].point.y()) && (pointsedge[i - 1]->point.y() > uvedge[up_index - 1].point.y())
								   || (pointsedge[i - 1]->Id == uvedge[up_index - 2].Id))//判断上个内部点的区间是否是【index-2，index-1】，==表示上个点被归并
							   {
								   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
							   }
						   }
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   if (((*tempNode1 - *tempNode4).GetLength() < length_differ) &&
						   (fabs(pointsedge[i]->point.y() - (uvedge[up_index - 1].point.y() + uvedge[up_index].point.y()) / 2)<(uvedge[up_index].point.y() - uvedge[up_index - 1].point.y()) / 5))
					   {//此处不考虑三个内部点的情况
						   pointsedge[i]->Id = nId;//将该点升级为边界点
						   /* pointsedge[i]->point.x() = uvedge[up_index].point.x();*/
						   addpoint.Id = nId;
						   addpoint.point.x() = uvedge[up_index].point.x();
						   addpoint.point.y() = pointsedge[i]->point.y();
						   uvEdge.push_back(*pointsedge[i]);//添加到全局边界容器
						   uvedge.insert(uvedge.begin() + up_index, addpoint);//添加到局部边界
						   addspace.setX(addpoint.point.x());
						   addspace.setY(addpoint.point.y());
						   addspace.setID(nId);
						   UVForOut.push_back(addspace);
						   nId++;
						   delete tempNode1;
						   delete tempNode2;
						   delete tempNode3;
						   delete tempNode4;
						   continue;//该点已处理完，处理下个i
					   }
					   //一般情况
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index - 1].Id);
					   point_m.setJnode(uvedge[up_index].Id);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode2;
					   delete tempNode3;
					   delete tempNode4;
				   }
				   else if (up_index == 0)//内部点只有上边界点与之相连
				   {					   
						Node *tempNode1;
						Node *tempNode2;
						tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
						tempNode2 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//下点坐标					   
						if ((*tempNode1 - *tempNode2).GetLength() < length_differ)//下点杆件长度若小于1m
						{
							pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
							/*pointsedge[i]->point.x() = uvedge[up_index].point.x();
							pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

							if (pointsedge[i + 1]->point.y() > uvedge[up_index + 1].point.y())//如果下一个内部点的区间为【0,1】，
							{
								pointsedge[i]->next.clear();//将该点的内部点自身的联系删除，这样可以保证没有重复杆件
							}
							delete tempNode1;
							delete tempNode2;
							continue;//该点已处理完，处理下个i
						}
						//都大于1时
						pointsedge[i]->Id = nId;//对该点进行整体编号
						point_m.set_ID(nId);
						point_m.setInode(-1);
						point_m.setJnode(uvedge[up_index].Id);
						edge_point_m.push_back(point_m);//将这种相连关系记录下来
						addspace.setX(pointsedge[i]->point.x());
						addspace.setY(pointsedge[i]->point.y());
						addspace.setID(nId);
						UVForOut.push_back(addspace);
						nId++;
						delete tempNode1;
						delete tempNode2;
					  
				   }
				   else if (up_index == -1)//内部点有下边界点与之相连
				   {
					   Node *tempNode1;
					   Node *tempNode3;
					   up_index = uvedge.size() - 1;//使up_index指向边界左点的下标
					   tempNode1 = GetSurfacePointatUandV(pointsedge[i]->point.x(), pointsedge[i]->point.y());
					   tempNode3 = GetSurfacePointatUandV(uvedge[up_index].point.x(), uvedge[up_index].point.y());//上点坐标					   
					   if ((*tempNode1 - *tempNode3).GetLength() < length_differ)//上点杆件长度若小于1m
					   {
						   pointsedge[i]->Id = uvedge[up_index].Id;//将内部点赋值为边界点，这样该点就被归并到边界上
						   /*pointsedge[i]->point.x() = uvedge[up_index].point.x();
						   pointsedge[i]->point.y() = uvedge[up_index].point.y();*/

						   if (pointsedge[i - 1]->point.y() < uvedge[up_index - 1].point.y() || (pointsedge[i - 1]->Id == uvedge[up_index - 1].Id))//如果上一个内部点的区间为【index-1,index】，，，==表示上个点被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
						   delete tempNode1;
						   delete tempNode3;
						   continue;//该点已处理完，处理下个i
					   }
					   //大于1时
					   pointsedge[i]->Id = nId;//对该点进行整体编号
					   point_m.set_ID(nId);
					   point_m.setInode(uvedge[up_index].Id);
					   point_m.setJnode(-1);
					   edge_point_m.push_back(point_m);//将这种相连关系记录下来
					   addspace.setX(pointsedge[i]->point.x());
					   addspace.setY(pointsedge[i]->point.y());
					   addspace.setID(nId);
					   UVForOut.push_back(addspace);
					   nId++;
					   if (i != 0)//保证不为第一个内部点
					   {
						   if (pointsedge[i - 1]->Id == uvedge[up_index].Id)//判断上个内部点是否被归并
						   {
							   pointsedge[i - 1]->next.clear();//将上一点的内部点自身的联系删除，这样可以保证没有重复杆件
						   }
					   }
					   delete tempNode1;
					   delete tempNode3;
				   }
			   }
			   //判断是否有两个内部点和相同的两个边界点相连
			   for (size_t i = 0; i < edge_point_m.size() - 1; i++)
			   {
				   for (size_t j = i + 1; j < edge_point_m.size(); j++)
				   {
					   if ((edge_point_m[i].getInode() == edge_point_m[j].getInode()) && (edge_point_m[i].getJnode() == edge_point_m[j].getJnode()))
					   {//如果有相同的则改成四边形相连
						   edge_point_m[i].setJnode(-1);
						   edge_point_m[j].setInode(-1);
					   }
				   }
			   }
			   //建立杆件
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

//处理边界自身
std::vector<Member> NurbsSurface_qj::EdgeTotal_Member(std::vector<PointStruct> uvEdge, std::vector<PLib::Point2Dd> extrmumVec, int num, double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV)
{
	std::vector<PointStruct> *uv_eg = new std::vector<PointStruct>[2 * num + 4];//存储各边界点
	std::vector<Member> mVecTemp;//用来存储新生成的边界member
	//偶数存u方向的点，奇数存v方向的点
	for (size_t i = 0; i < uvEdge.size(); ++i)//遍历uvEdge
	{
		if (uvEdge[i].Id != -1)//选择使用到的边界点
		{
			if (uvEdge[i].point.y() == minV)//下边界
			{
				uv_eg[0].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.x() == maxU)//右边界
			{
				uv_eg[1].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.y() == maxV)//上边界
			{
				uv_eg[2].push_back(uvEdge[i]);
			}
			if (uvEdge[i].point.x() == minU)//左边界
			{
				uv_eg[3].push_back(uvEdge[i]);
			}

			for (size_t j = 4; j < 2 * num + 4; j +=2)//内部边界u方向
			{
				if (uvEdge[i].point.y() == extrmumVec[(j - 4) / 2].y())
				{
					uv_eg[j].push_back(uvEdge[i]);
				}
			}

			for (size_t j = 5; j < 2 * num + 4; j+=2)//内部边界v方向
			{
				if (uvEdge[i].point.x() == extrmumVec[(j - 5) / 2].x())
				{
					uv_eg[j].push_back(uvEdge[i]);
				}
			}

		}
	}
	//点排序

	std::sort(uv_eg[0].begin(), uv_eg[0].end(), compare1_x);//下边界
	std::sort(uv_eg[1].begin(), uv_eg[1].end(), compare3_y);//右边界
	std::sort(uv_eg[2].begin(), uv_eg[2].end(), compare1_x);//上边界
	std::sort(uv_eg[3].begin(), uv_eg[3].end(), compare3_y);//左边界
	for (size_t j = 4; j < 2 * num + 4; j +=2)
	{
		std::sort(uv_eg[j].begin(), uv_eg[j].end(), compare1_x);//内部u边界
	}
	for (size_t j = 5; j < 2 * num + 4; j +=2)
	{
		std::sort(uv_eg[j].begin(), uv_eg[j].end(), compare3_y);//内部v边界
	}
	for (size_t i = 0; i < 2 * num + 4; ++i)//遍历uv_eg
	{
		for (size_t j = 0; j < uv_eg[i].size() - 1; ++j)//读取uv_eg[i]的前size-1个元素
		{
			Member m;
			Node *tempNode1;
			Node *tempNode2;
			m.setInode(uv_eg[i][j].Id);//写入上点
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
	//处理边界三角区域

	for (size_t i = 4; i < 2 * num + 4; i +=2)//u方向的起点
	{
		for (size_t j = 1; j < 2 * num + 4; j +=2)//遍历v方向的uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());//u方向上的起点对应的三维坐标
			tempNode2 = GetSurfacePointatUandV(uv_eg[j].begin()->point.x(), uv_eg[i].begin()->point.y());//该起点在v方向边界上的投影点的三维坐标
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthU)//判断点距离那个v方向边界最近
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//寻找最近v边界上该起点的上下边界点
				{
					if (uv_eg[i].begin()->point.y()< iter->point.y())//边界点都是从小到大排序的，当增大到一定时就找到了上下边界点
					{

						Member m;
						m.setInode(uv_eg[i].begin()->Id);//写入上点
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode(uv_eg[i].begin()->Id);//写入下点
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
	for (size_t i = 4; i < 2 * num + 4; i +=2)//u方向的终点
	{
		for (size_t j = 1; j < 2 * num + 4; j +=2)//遍历v方向的uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());//u方向上的起点对应的三维坐标
			tempNode2 = GetSurfacePointatUandV(uv_eg[j].begin()->point.x(), (uv_eg[i].end() - 1)->point.y());//该起点在v方向边界上的投影点的三维坐标
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthU)//判断点距离那个v方向边界最近
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//寻找最近v边界上该起点的上下边界点
				{
					if ((uv_eg[i].end() - 1)->point.y()< iter->point.y())//边界点都是从小到大排序的，当增大到一定时就找到了上下边界点
					{
						Member m;
						m.setInode((uv_eg[i].end() - 1)->Id);//写入上点
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode((uv_eg[i].end() - 1)->Id);//写入下点
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

	for (size_t i = 5; i < 2 * num + 4; i +=2)//v方向的起点
	{
		for (size_t j = 0; j < 2 * num + 4; j +=2)//遍历u方向的uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());//v方向上的起点对应的三维坐标
			tempNode2 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[j].begin()->point.y());//该起点在u方向边界上的投影点的三维坐标
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthV)//判断点距离那个u方向边界最近
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//寻找最近u边界上该起点的上下边界点
				{
					if (uv_eg[i].begin()->point.x()< iter->point.x())//边界点都是从小到大排序的，当增大到一定时就找到了上下边界点
					{
						Member m;
						m.setInode(uv_eg[i].begin()->Id);//写入上点
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV(uv_eg[i].begin()->point.x(), uv_eg[i].begin()->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode(uv_eg[i].begin()->Id);//写入下点
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
	for (size_t i = 5; i < 2 * num + 4; i +=2)//v方向的终点
	{
		for (size_t j = 0; j < 2 * num + 4; j +=2)//遍历u方向的uv_eg
		{
			Node *tempNode1;
			Node *tempNode2;
			tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());//v方向上的起点对应的三维坐标
			tempNode2 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), uv_eg[j].begin()->point.y());//该起点在u方向边界上的投影点的三维坐标
			if ((*tempNode1 - *tempNode2).GetLength() <= splitLengthV)//判断点距离那个u方向边界最近
			{
				for (std::vector<PointStruct>::iterator iter = uv_eg[j].begin(); iter != uv_eg[j].end(); ++iter)//寻找最近u边界上该起点的上下边界点
				{
					if ((uv_eg[i].end() - 1)->point.x()< iter->point.x())//边界点都是从小到大排序的，当增大到一定时就找到了上下边界点
					{
						Member m;
						m.setInode((uv_eg[i].end() - 1)->Id);//写入上点
						m.setJnode(iter->Id);
						tempNode1 = GetSurfacePointatUandV((uv_eg[i].end() - 1)->point.x(), (uv_eg[i].end() - 1)->point.y());
						tempNode2 = GetSurfacePointatUandV(iter->point.x(), iter->point.y());
						m.setCutLength((*tempNode1 - *tempNode2).GetLength());
						m.set_ID(mId);
						mId++;
						mVecTemp.push_back(m);
						m.setInode((uv_eg[i].end() - 1)->Id);//写入下点
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

//处理边界的角点
std::vector<Member> NurbsSurface_qj::EdgeCorner_Member(std::vector<PointStruct*> &point_corner, double length_differ)
{
	
	Member m;//用于增加edge_m
	std::vector<Member> edge_m;//存储边界杆件	
	Node addspace;
	int index34 = -1, index56 = -1;//标记点0与点3456的关系
	int max3456 = -1;
	

	if (point_corner[0]->point.x()>point_corner[6]->point.x())
	{
		if (point_corner[0]->point.x()>point_corner[5]->point.x())//0大于56
		{
			index56 = 6;//0连6
		}
		else
		{
			index56 = 56;//0连56
		}
	}
	else
	{
		if (point_corner[0]->point.x()>point_corner[5]->point.x())//0大于56
		{
			index56 = 56;//0连56
		}
		else
		{
			index56 = 6;//0连6
		}
	}

	if (point_corner[0]->point.y()>point_corner[4]->point.y())
	{
		if (point_corner[0]->point.y()>point_corner[3]->point.y())//0大于56
		{
			index34= 4;//0连4
		}
		else
		{
			index34 = 34;//0连34
		}
	}
	else
	{
		if (point_corner[0]->point.y()>point_corner[3]->point.y())//
		{
			index34 = 34;//0连34
		}
		else
		{
			index34 = 4;//0连4
		}
	}

	Node *tempNode0;
	Node *tempNode3;
	Node *tempNode4;
	Node *tempNode5;
	Node *tempNode6;
	Node *tempNode1;
	Node *tempNode2;

	tempNode0 = GetSurfacePointatUandV(point_corner[0]->point.x(), point_corner[0]->point.y());//角点
	tempNode3 = GetSurfacePointatUandV(point_corner[3]->point.x(), point_corner[3]->point.y());//3坐标
	tempNode4 = GetSurfacePointatUandV(point_corner[4]->point.x(), point_corner[4]->point.y());//4坐标
	tempNode5 = GetSurfacePointatUandV(point_corner[5]->point.x(), point_corner[5]->point.y());//5坐标
	tempNode6 = GetSurfacePointatUandV(point_corner[6]->point.x(), point_corner[6]->point.y());//6坐标
	tempNode1 = ((*tempNode0 - *tempNode3).GetLength() < (*tempNode0 - *tempNode4).GetLength()) ? tempNode3 :tempNode4;
	tempNode2 = ((*tempNode0 - *tempNode5).GetLength() < (*tempNode0 - *tempNode6).GetLength()) ? tempNode5 : tempNode6;
	max3456 = ((*tempNode0 - *tempNode1).GetLength() < (*tempNode0 - *tempNode2).GetLength()) ? 1 : 2;//1表示34点，2表示56点


	//角点与右边界更近
	if ((max3456==1) && ((*tempNode0- *tempNode1).GetLength() < length_differ))
	{	
		//与点4更近
		if (fabs((point_corner[0]->point.y() - point_corner[4]->point.y())) < fabs((point_corner[3]->point.y() - point_corner[4]->point.y()) /2))
		{
			point_corner[0]->Id = point_corner[4]->Id;//将点0归并到点4
			/*point_corner[0]->point.x() = point_corner[4]->point.x();
			point_corner[0]->point.y() = point_corner[4]->point.y();*/
			if ((point_corner[2]->Id == point_corner[3]->Id) || (point_corner[2]->Id == point_corner[4]->Id))//判断点2yu3,4是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
			}
		}
		//与点3更近
		if (fabs((point_corner[0]->point.y() - point_corner[3]->point.y())) < fabs((point_corner[3]->point.y() - point_corner[4]->point.y()) / 2))
		{
			point_corner[0]->Id = point_corner[3]->Id;//将点0归并到点3
			/*point_corner[0]->point.x() = point_corner[3]->point.x();
			point_corner[0]->point.y() = point_corner[3]->point.y();*/
			if (point_corner[2]->Id == point_corner[3]->Id)//判断点2yu3是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
			}
		}
		if (index56 == 56)//点0连5,6
		{
			if ((point_corner[1]->Id == point_corner[5]->Id) || (point_corner[2]->Id == point_corner[6]->Id))//判断点1与5,6是否重合
			{
				point_corner[1]->next.clear();//将点1的自带关系删除				
			}
		}
		if (index56 == 6)//点0连6
		{
			if ((point_corner[1]->Id == point_corner[6]->Id))//判断点1与6是否重合
			{
				point_corner[1]->next.clear();//将点1的自带关系删除				
			}
		}			
	}	

	//角点与上边界更近
	if ((max3456==2) && ((*tempNode0 - *tempNode2).GetLength() < length_differ))
	{
		//与点6更近
		if (fabs((point_corner[0]->point.x() - point_corner[6]->point.x())) < fabs((point_corner[5]->point.x() - point_corner[6]->point.x()) / 2))
		{
			point_corner[0]->Id = point_corner[6]->Id;//将点0归并到点4
			/*point_corner[0]->point.x() = point_corner[6]->point.x();
			point_corner[0]->point.y() = point_corner[6]->point.y();*/
			if ((point_corner[1]->Id == point_corner[5]->Id) || (point_corner[1]->Id == point_corner[6]->Id))//判断点1yu5,6是否重合
			{
				point_corner[1]->next.clear();//将点1的自带关系删除				
			}
		}
		//与点5更近
		if (fabs((point_corner[0]->point.x() - point_corner[5]->point.x())) < fabs((point_corner[5]->point.x() - point_corner[6]->point.x()) / 2))
		{
			point_corner[0]->Id = point_corner[5]->Id;//将点0归并到点5
			/*point_corner[0]->point.x() = point_corner[5]->point.x();
			point_corner[0]->point.y() = point_corner[5]->point.y();*/
			if (point_corner[1]->Id == point_corner[5]->Id)//判断点1yu5是否重合
			{
				point_corner[1]->next.clear();//将点2的自带关系删除				
			}
		}
		if (index34 == 34)//点0连34
		{
			if ((point_corner[2]->Id == point_corner[3]->Id) || (point_corner[2]->Id == point_corner[4]->Id))//判断点2与3,4是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
			}
		}
		if (index34 == 4)//点0连4
		{
			if ((point_corner[2]->Id == point_corner[4]->Id))//判断点2与4是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
			}
		}
	}

	//两个方向长度大于2m
	if (((*tempNode0 - *tempNode1).GetLength() > length_differ) && ((*tempNode0 - *tempNode2).GetLength() > length_differ))
	{
		point_corner[0]->Id = nId;//对该点进行整体编号			
		addspace.setX(point_corner[0]->point.x());
		addspace.setY(point_corner[0]->point.y());
		addspace.setID(nId);
		UVForOut.push_back(addspace);
		nId++;
		if (index34 == 34)
		{
			if ((point_corner[2]->Id == point_corner[3]->Id))//判断点2与3是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
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
			if ((point_corner[2]->Id == point_corner[4]->Id))//判断点2与4是否重合
			{
				point_corner[2]->next.clear();//将点2的自带关系删除				
			}
			m.setInode(point_corner[0]->Id);
			m.setJnode(point_corner[4]->Id);
			m.set_ID(mId);
			edge_m.push_back(m);
			mId++;
		}

		if (index56 == 56)
		{
			if ((point_corner[1]->Id == point_corner[5]->Id))//判断点1与5是否重合
			{
				point_corner[1]->next.clear();//将点2的自带关系删除				
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
			if ((point_corner[1]->Id == point_corner[6]->Id))//判断点1与6是否重合
			{
				point_corner[1]->next.clear();//将点1的自带关系删除				
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


//******************************************取z最大值点，只分四块进行曲面划分****************************************************************
std::vector< Member> NurbsSurface_qj::SurfaceAverageSpit_Rectangle(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec)
{
	//找曲面z绝对值最大点
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极小值点
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极大值点
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//存储极值点
	int num = 0;
	for (std::vector<PLib::Point2Dd>::iterator itr = pMinExtrmums.begin(); itr != pMinExtrmums.end(); ++itr)//将极小值加入极值点数组
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}
	//pMaxExtrmums.push_back(PLib::Point2Dd(0.5, 0.5));
	for (std::vector<PLib::Point2Dd>::iterator itr = pMaxExtrmums.begin(); itr != pMaxExtrmums.end(); ++itr)//将极大值加入极值点数组
	{
		extrmums[num].x = (*itr).x();
		extrmums[num].y = (*itr).y();
		num++;
	}

	CPoint2D maxPoint;//存储z最大的点 
	double MaxzofPoint;
	maxPoint = extrmums[0];
	Node *tempNode1;
	tempNode1 = GetSurfacePointatUandV(extrmums[0].x, extrmums[0].y);
	MaxzofPoint = fabs(tempNode1->getZ());
	delete tempNode1;
	for (int i = 1; i<num; i++)//寻找曲面z最大值	
	{
		tempNode1 = GetSurfacePointatUandV(extrmums[i].x, extrmums[i].y);
		if (fabs(tempNode1->getZ()) > MaxzofPoint)
		{
			MaxzofPoint = fabs(tempNode1->getZ());
			maxPoint = extrmums[i];
		}
		delete tempNode1;
	}
	//赋值四个分块的边界
	rectangleRegion rec[4];
	rec[0].maxU = maxU; rec[0].minU = maxPoint.x; rec[0].maxV = maxV; rec[0].minV = maxPoint.y; rec[0].type = 0; rec[0].type = 0;
	rec[1].maxU = maxPoint.x; rec[1].minU = minU; rec[1].maxV = maxV; rec[1].minV = maxPoint.y; rec[1].type = 1; rec[1].type = 1;
	rec[2].maxU = maxPoint.x; rec[2].minU = minU; rec[2].maxV = maxPoint.y; rec[2].minV = minV; rec[2].type = 2; rec[2].type = 2;
	rec[3].maxU = maxU; rec[3].minU = maxPoint.x; rec[3].maxV = maxPoint.y; rec[3].minV = minV; rec[3].type = 3; rec[3].type = 3;

	

	//分块生成杆件
	std::vector< Member> mVec;//用于记录总的各个分块内的member
	std::vector<Member> mVecTemp;//临时，用于记录循环中分块内的member，形参	
	std::vector<IdRelationship> IdRel_Rectangle;//记录四边形四点关系，逆时针

	//等分内部边界，存入uvEdge中
	std::vector<PointStruct> uvEdge;
	mVecTemp = InsideEdge_Rectangle(minU, maxU, minV, maxV, maxPoint, splitLengthU, splitLengthV, nVec, uvEdge);
	(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//将mVecTemp加到mVec中,四根杆件

	for (int i = 0; i < 4; ++i)//对所有区域进行等分，并将等分的点和杆信息分别存储在nVec和mVec中，4分块数		
	{
		mVecTemp = MemberofRec_Rectangle(rec[i], splitLengthU, splitLengthU, nVec, uvEdge/*,Gui_index*/);//将i块进行等分,生成杆件
		(mVec).insert(mVec.end(), mVecTemp.begin(), mVecTemp.end());//将mVecTemp加到mVec中	
	}	

	return mVec;	
}

//内部边界处理
std::vector<Member> NurbsSurface_qj::InsideEdge_Rectangle(double minU, double maxU, double minV, double maxV, CPoint2D maxPoint, double splitLengthU, double splitLengthV, std::vector< Node> & nVec, std::vector<PointStruct> &uvEdge)
{
	std::vector<double> temp;
	PointStruct divPoint;
	Node AddSpace;
	std::vector<Member> mVecTemp;
	Member m;

	divPoint.point.x() = maxU;//点0
	divPoint.point.y() = maxPoint.y;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxU);
	AddSpace.setY(maxPoint.y);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = maxPoint.x;//点1
	divPoint.point.y() =maxV ;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxPoint.x);
	AddSpace.setY(maxV);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = minU;//点2
	divPoint.point.y() = maxPoint.y;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(minU);
	AddSpace.setY(maxPoint.y);
	nVec.push_back(AddSpace);
	nId++;

	divPoint.point.x() = maxPoint.x;//点3
	divPoint.point.y() = minV;
	divPoint.Id = nId;
	uvEdge.push_back(divPoint);
	AddSpace.setID(nId);
	AddSpace.setX(maxPoint.x);
	AddSpace.setY(minV);
	nVec.push_back(AddSpace);
	nId++;

	temp = getUAverageLength(minU, maxPoint.x, maxPoint.y, splitLengthU, 1);//从右往左等分，对应type=1，2,记录等分点
	for (int j = 0; j < temp.size(); ++j)
	{
		divPoint.point.x() = temp[j];//将左半temp的u存成PointStruct，方便存入uvEdge中
		divPoint.point.y() = maxPoint.y;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//将等分边界点存入uvEdge		
		AddSpace.setID(nId);
		AddSpace.setX(temp[j]);
		AddSpace.setY(maxPoint.y);
		nVec.push_back(AddSpace);//坐标以后再赋值，因为uvEdge中点有可能变化
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[2]);//左点与点2相连
	mVecTemp.push_back(m);


	temp = getUAverageLength(maxPoint.x, maxU, maxPoint.y, splitLengthU, 0);//从左往右等分，对应type=0,3，右
	temp.erase(temp.begin());//除去极值点，因为上面已经添加
	for (int j = 0; j < temp.size(); ++j)
	{
		divPoint.point.x() = temp[j];//将右半temp的u存成PointStruct，方便存入uvEdge中			
		divPoint.point.y() = maxPoint.y;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//将等分边界点存入uvEdge		
		AddSpace.setID(nId);
		AddSpace.setX(temp[j]);
		AddSpace.setY(maxPoint.y);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[0]);//左点与点2相连
	mVecTemp.push_back(m);


	temp = getVAverageLength(minV, maxPoint.y, maxPoint.x, splitLengthV, 2);//从上往下等分，对应type=2，3,记录等分点
	temp.erase(temp.begin());//除去极值点，因为上面已经添加
	for (int j = 0; j != temp.size(); ++j)
	{
		divPoint.point.y() = temp[j];//将下半temp的v存成PointStruct，方便存入vVec中				
		divPoint.point.x() = maxPoint.x;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//将等分边界点存入uvEdge		
		AddSpace.setID(nId);
		AddSpace.setY(temp[j]);
		AddSpace.setX(maxPoint.x);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[3]);//左点与点2相连
	mVecTemp.push_back(m);


	temp = getVAverageLength(maxPoint.y, maxV, maxPoint.x, splitLengthV, 0);//从下往上等分，对应type=0,1，
	temp.erase(temp.begin());//除去极值点，因为上面已经添加
	for (int j = 0; j != temp.size(); ++j)
	{
		divPoint.point.y() = temp[j];//将上半temp的v存成PointStruct，方便存入vVec中					
		divPoint.point.x() = maxPoint.x;
		divPoint.Id = nId;
		uvEdge.push_back(divPoint);//将等分边界点存入uvEdge		
		AddSpace.setID(nId);
		AddSpace.setY(temp[j]);
		AddSpace.setX(maxPoint.x);
		nVec.push_back(AddSpace);
		nId++;
	}
	m = Add_Member(uvEdge.back(), uvEdge[1]);//左点与点2相连
	mVecTemp.push_back(m);

	return mVecTemp;
}

//分块画杆件
std::vector<Member> NurbsSurface_qj::MemberofRec_Rectangle(rectangleRegion rec, double splitLengthU, double splitLengthV, std::vector<Node> &nVec, std::vector<PointStruct> &uvEdge/*, std::vector<bool> &Gui_index*/)
{
	//划分曲面
	double maxU = rec.maxU + 0.2;
	//***********此处0.3有待优化***************************************************
	double maxV = rec.maxV + 0.2;
	double minU = rec.minU - 0.2;
	double minV = rec.minV - 0.2;

	int u_num;
	std::vector<std::vector<PointStruct> > points(1);//定义二维vector
	PointStruct AddPoint;//用于增加points的size
	AddPoint.Id = -1;//初始化编号
	AddPoint.point = PLib::Point2Dd(-1, -1);//初始化uv
	
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
	//case 0://从左下往右上找点
	//	u = rec.minU;
	//	v = rec.minV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u <= maxU)//将第一行的坐标点找到	
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
	//		points.resize(row1 + 1);//增加一行
	//		nextV = nextVPoint(u, v, splitLengthV);//寻找第二行的第一个点v坐标
	//	
	//		points[row1].push_back(AddPoint);//将第二行的size增加1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点
	//		

	//		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************此处要删吗

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV);//点1与2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1 )//如果找到的坐标为-1 -1，则画三角形
	//			{
	//				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV);//点1与3
	//					//添加点之间关系
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					points[row1].push_back(AddPoint);//将第二行的size增加1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));
	//					
	//					col1 = col1 + 1;//第二行增加一点
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV);//点1与4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//添加点之间关系
	//						points[row1].push_back(AddPoint);//将第二行的size增加1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0+1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//如果col0+2不存在
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));	
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));							
	//						break;
	//					}
	//				}
	//			}
	//			//正常情况
	//			points[row1].push_back(AddPoint);//将第二行的size增加1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0+1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));				
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 1://从右下往左上
	//	u = rec.maxU;
	//	v = rec.minV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u >=minU)//将第一行的坐标点找到	
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
	//		points.resize(row1 + 1);//增加一行
	//		nextV = nextVPoint(u, v, splitLengthV,1);//寻找第二行的第一个点v坐标

	//		points[row1].push_back(AddPoint);//将第二行的size增加1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


	//		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************此处要删吗

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,1);//点1与2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
	//			{
	//				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,1);//点1与3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//添加点之间关系
	//					points[row1].push_back(AddPoint);//将第二行的size增加1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//第二行增加一点
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,1);//点1与4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//添加点之间关系
	//						points[row1].push_back(AddPoint);//将第二行的size增加1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//如果col0+2不存在
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//正常情况
	//			points[row1].push_back(AddPoint);//将第二行的size增加1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 2://从右上往左下找
	//	u = rec.maxU;
	//	v = rec.maxV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u >=minU)//将第一行的坐标点找到	
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
	//		points.resize(row1 + 1);//增加一行
	//		nextV = nextVPoint(u, v, splitLengthV,2);//寻找第二行的第一个点v坐标

	//		points[row1].push_back(AddPoint);//将第二行的size增加1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


	//		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************此处要删吗

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,2);//点1与2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
	//			{
	//				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,2);//点1与3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//添加点之间关系
	//					points[row1].push_back(AddPoint);//将第二行的size增加1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//第二行增加一点
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,2);//点1与4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//添加点之间关系
	//						points[row1].push_back(AddPoint);//将第二行的size增加1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//如果col0+2不存在
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//正常情况
	//			points[row1].push_back(AddPoint);//将第二行的size增加1
	//			points[row1][col1 + 1].point = nextCrossP;
	//			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//			col0 = col0 + 1;
	//			col1 = col1 + 1;
	//		}
	//		v = nextV;
	//	}
	//	break;

	//case 3://从左上往右下找
	//	u = rec.minU;
	//	v = rec.maxV;
	//	points[0].push_back(AddPoint);
	//	points[0][0].point = PLib::Point2Dd(u, v);

	//	while (u <= maxU)//将第一行的坐标点找到	
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
	//		points.resize(row1 + 1);//增加一行
	//		nextV = nextVPoint(u, v, splitLengthV,3);//寻找第二行的第一个点v坐标

	//		points[row1].push_back(AddPoint);//将第二行的size增加1
	//		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


	//		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
	//		{
	//			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
	//			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
	//			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

	//			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
	//			delete tempNode1;
	//			delete tempNode2;
	//			delete tempNode3;//*******************************************************************************此处要删吗

	//			nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV,3);//点1与2

	//			
	//			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
	//			{
	//				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
	//				{
	//					nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV,3);//点1与3
	//					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//					{
	//						break;
	//					}
	//					//添加点之间关系
	//					points[row1].push_back(AddPoint);//将第二行的size增加1
	//					points[row1][col1 + 1].point = nextCrossP;
	//					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
	//					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
	//					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

	//					col1 = col1 + 1;//第二行增加一点
	//					continue;
	//				}
	//				else
	//				{
	//					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
	//					{
	//						nextCrossP = nextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV,3);//点1与4
	//						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
	//						{
	//							break;
	//						}
	//						//添加点之间关系
	//						points[row1].push_back(AddPoint);//将第二行的size增加1
	//						points[row1][col1 + 1].point = nextCrossP;
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

	//						col1 = col1 + 1;
	//						col0 = col0 + 2;
	//						continue;
	//					}
	//					else//如果col0+2不存在
	//					{
	//						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
	//						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
	//						points[row0][col0+1].next.push_back(pointCoordinate(-1, -1));
	//						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
	//						break;
	//					}
	//				}
	//			}
	//			//正常情况
	//			points[row1].push_back(AddPoint);//将第二行的size增加1
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
case 0://从左下往右上找点
	u = rec.minU;
	v = rec.minV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u <= maxU)//将第一行的坐标点找到	
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
		points.resize(row1 + 1);//增加一行
		nextV = nextVPoint(u, v, splitLengthV);//寻找第二行的第一个点v坐标

		points[row1].push_back(AddPoint);//将第二行的size增加1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************此处要删吗

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV);//点1与2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
			{
				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV);//点1与3
					//添加点之间关系
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					points[row1].push_back(AddPoint);//将第二行的size增加1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//第二行增加一点
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV);//点1与4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//添加点之间关系
						points[row1].push_back(AddPoint);//将第二行的size增加1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//如果col0+2不存在
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//正常情况
			points[row1].push_back(AddPoint);//将第二行的size增加1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 1://从右下往左上
	u = rec.maxU;
	v = rec.minV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u >= minU)//将第一行的坐标点找到	
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
		points.resize(row1 + 1);//增加一行
		nextV = nextVPoint(u, v, splitLengthV, 1);//寻找第二行的第一个点v坐标

		points[row1].push_back(AddPoint);//将第二行的size增加1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************此处要删吗

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 1);//点1与2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
			{
				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 1);//点1与3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//添加点之间关系
					points[row1].push_back(AddPoint);//将第二行的size增加1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//第二行增加一点
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 1);//点1与4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//添加点之间关系
						points[row1].push_back(AddPoint);//将第二行的size增加1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//如果col0+2不存在
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//正常情况
			points[row1].push_back(AddPoint);//将第二行的size增加1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 2://从右上往左下找
	u = rec.maxU;
	v = rec.maxV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u >= minU)//将第一行的坐标点找到	
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
		points.resize(row1 + 1);//增加一行
		nextV = nextVPoint(u, v, splitLengthV, 2);//寻找第二行的第一个点v坐标

		points[row1].push_back(AddPoint);//将第二行的size增加1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************此处要删吗

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 2);//点1与2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
			{
				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 2);//点1与3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//添加点之间关系
					points[row1].push_back(AddPoint);//将第二行的size增加1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//第二行增加一点
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 2);//点1与4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//添加点之间关系
						points[row1].push_back(AddPoint);//将第二行的size增加1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//如果col0+2不存在
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//正常情况
			points[row1].push_back(AddPoint);//将第二行的size增加1
			points[row1][col1 + 1].point = nextCrossP;
			points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
			points[row0][col0].next.push_back(pointCoordinate(row1, col1));
			col0 = col0 + 1;
			col1 = col1 + 1;
		}
		v = nextV;
	}
	break;

case 3://从左上往右下找
	u = rec.minU;
	v = rec.maxV;
	points[0].push_back(AddPoint);
	points[0][0].point = PLib::Point2Dd(u, v);

	while (u <= maxU)//将第一行的坐标点找到	
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
		points.resize(row1 + 1);//增加一行
		nextV = nextVPoint(u, v, splitLengthV, 3);//寻找第二行的第一个点v坐标

		points[row1].push_back(AddPoint);//将第二行的size增加1
		points[row1][col1].point = PLib::Point2Dd(u, nextV);//第二行第一个点


		for (; col0 < points[row0].size() - 1;)//col0取到第一行的倒数第二个点
		{
			tempNode1 = GetSurfacePointatUandV(points[row1][col1].point.x(), points[row1][col1].point.y());//第二行第一个点
			tempNode2 = GetSurfacePointatUandV(points[row0][col0 + 1].point.x(), points[row0][col0 + 1].point.y());//第一行第二个点
			tempNode3 = GetSurfacePointatUandV(points[row0][col0].point.x(), points[row0][col0].point.y());//第一行第一个点

			cosValue = ((*tempNode1 - *tempNode3).GetNormal()) | ((*tempNode2 - *tempNode3).GetNormal());
			delete tempNode1;
			delete tempNode2;
			delete tempNode3;//*******************************************************************************此处要删吗

			nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 1].point, splitLengthU, splitLengthV, 3);//点1与2


			if (nextCrossP.x() == -1 && nextCrossP.y() == -1)//如果找到的坐标为-1 -1，则画三角形
			{
				if (cosValue < 0)//若三点构成的三角形的角度大于90度，则根据三角形的左边两点找下一点
				{
					nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0].point, splitLengthU, splitLengthV, 3);//点1与3
					if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
					{
						break;
					}
					//添加点之间关系
					points[row1].push_back(AddPoint);//将第二行的size增加1
					points[row1][col1 + 1].point = nextCrossP;
					points[row0][col0].next.push_back(pointCoordinate(-1, -1));
					points[row0][col0].next.push_back(pointCoordinate(row1, col1));//只需连1、3点，2、3点下次连
					//points[row0][col0].next.push_back(pointCoordinate(row1, col1 + 1));

					col1 = col1 + 1;//第二行增加一点
					continue;
				}
				else
				{
					if (col0 + 2 < points[row0].size())//如果第一行col0+2点存在
					{
						nextCrossP = NextCrossPoint(points[row1][col1].point, points[row0][col0 + 2].point, splitLengthU, splitLengthV, 3);//点1与4
						if (nextCrossP.x() == -1 && nextCrossP.y() == -1)
						{
							break;
						}
						//添加点之间关系
						points[row1].push_back(AddPoint);//将第二行的size增加1
						points[row1][col1 + 1].point = nextCrossP;
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row0, col0 + 2));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));

						col1 = col1 + 1;
						col0 = col0 + 2;
						continue;
					}
					else//如果col0+2不存在
					{
						points[row0][col0].next.push_back(pointCoordinate(row0, col0 + 1));
						points[row0][col0].next.push_back(pointCoordinate(row1, col1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(-1, -1));
						points[row0][col0 + 1].next.push_back(pointCoordinate(row1, col1));
						break;
					}
				}
			}
			//正常情况
			points[row1].push_back(AddPoint);//将第二行的size增加1
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
	//杆件生成	
	std::vector<Member> mDisplayVec; //信息流，存储杆的容器，杆的顺序与杆ID完全一致	
	Node addspace;//用来增加NodeForOut容器的大小
	Member m;
	std::vector<PointStruct> EdgePoints;//存储边界点
	//bool guibingPoint = false;//记录上点是否为内部点归并，是则为true，默认为false
	//double error_guibing = 0.2;//归并误差
	PointStruct AddPoint2;
	
	
	switch (rec.type)
	{
	case 0:
	{
			  //对points进行编号,分块内部点，除去边界
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()<rec.maxU && points[i][j].point.y()<rec.maxV)//判断点在分块内部
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

			  //判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//确保点在分块内，//*****************************需变动
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
			  //判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点,此处i=1，表明左下角点在上面已经编号了
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//确保点在分块内//*****************************需变动
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

			  //从下边界开始，从第一列开始，此处j从0开始，j按照如下规则，分块0、2从0开始，分块1、3从1开始，这样可以防止极值点处的杆件重复
			  EdgePoints.push_back(uvEdge[0]);//将点0加入到边界点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//确保点在分块内，//*****************************需变动
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//判断是否为最右点
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//判断最右点的相连上点是否在内部，在为真
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.maxU, rec.maxV, 0);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  for (size_t k = 0; k < points[0][j].next.size(); k++)//此处考虑三角形情况
						  {
							  m = Add_Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);
							  mDisplayVec.push_back(m);
						  }
					  }
				  }
			  }

			  //对内部点建立杆件，除去下边界、左边界
			  for (size_t i = 1; i < points.size(); i++)//从下往上
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//从右往左
				  {
					  if (points[i][j].Id != -1)//找到内部点
					  {
						  if (points[i][j + 1].Id == -1)//点i，j的右点在外部
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.maxU, rec.maxV, 0);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//点i，j的上点在外部，此处没有考虑三角形情况
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.maxU, rec.maxV, 0);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

			  //处理左边界,从第二行开始

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//确保点在分块内，//*****************************需变动
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//判断是否为最上点
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//判断最上点的相连右点是否在内部，在为真
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.maxU, rec.maxV, 0);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  for (size_t k = 0; k < points[i][0].next.size(); k++)//此处考虑三角形情况
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
			  //对points进行编号,分块内部点，除去边界
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()>rec.minU && points[i][j].point.y()<rec.maxV)//判断点在分块内部
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

			  //判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//确保点在分块内，//*****************************需变动
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
			  //判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点,此处i=1，表明左下角点在上面已经编号了
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//确保点在分块内//*****************************需变动
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

			  //从下边界开始，从第一列开始，此处j从1开始，j按照如下规则，分块0、2从0开始，分块1、3从1开始，这样可以防止极值点处的杆件重复
			  EdgePoints.push_back(uvEdge[2]);//将点0加入到边界点
			  for (int j = 1; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//确保点在分块内，//*****************************需变动
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//判断是否为最右点
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//判断最右点的相连上点是否在内部，在为真
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.minU, rec.maxV, 1);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  
							  m = Add_Member(points[0][j], points[points[0][j].next[1].x][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  
					  }
				  }
			  }

			  //对内部点建立杆件，除去下边界、左边界
			  for (size_t i = 1; i < points.size(); i++)//从下往上
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//从右往左
				  {
					  if (points[i][j].Id != -1)//找到内部点
					  {
						  if (points[i][j + 1].Id == -1)//点i，j的右点在外部
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.minU, rec.maxV, 1);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//点i，j的上点在外部，此处没有考虑三角形情况
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.minU, rec.maxV, 1);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

			  //处理左边界,从第二行开始

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() < rec.maxV)//确保点在分块内，//*****************************需变动
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//判断是否为最上点
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//判断最上点的相连右点是否在内部，在为真
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.minU, rec.maxV, 1);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  
							m = Add_Member(points[i][0], points[points[i][0].next[0].x][points[i][0].next[0].y]);//只连next【0】，next【1】在分部0中已经建立
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
			  //对points进行编号,分块内部点，除去边界
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()>rec.minU && points[i][j].point.y()>rec.minV)//判断点在分块内部
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

			  //判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//确保点在分块内，//*****************************需变动
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
			  //判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点,此处i=1，表明左下角点在上面已经编号了
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//确保点在分块内//*****************************需变动
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

			  //从下边界开始，从第一列开始，此处j从0开始，j按照如下规则，分块0、2从0开始，分块1、3从1开始，这样可以防止极值点处的杆件重复
			  EdgePoints.push_back(uvEdge[2]);//将点0加入到边界点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() > rec.minU)//确保点在分块内，//*****************************需变动
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//判断是否为最右点
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//判断最右点的相连上点是否在内部，在为真
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.minU, rec.minV,2);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  for (size_t k = 0; k < points[0][j].next.size(); k++)//此处考虑三角形情况
						  {
							  m = Add_Member(points[0][j], points[points[0][j].next[k].x][points[0][j].next[k].y]);
							  mDisplayVec.push_back(m);
						  }						  
					  }
				  }
			  }

			  //对内部点建立杆件，除去下边界、左边界
			  for (size_t i = 1; i < points.size(); i++)//从下往上
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//从右往左
				  {
					  if (points[i][j].Id != -1)//找到内部点
					  {
						  if (points[i][j + 1].Id == -1)//点i，j的右点在外部
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.minU, rec.minV, 2);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//点i，j的上点在外部，此处没有考虑三角形情况
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.minU, rec.minV, 2);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

			  //处理左边界,从第二行开始

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//确保点在分块内，//*****************************需变动
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//判断是否为最上点
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//判断最上点的相连右点是否在内部，在为真
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.minU, rec.minV, 2);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  for (size_t k = 0; k < points[i][0].next.size(); k++)//此处考虑三角形情况
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
			  //对points进行编号,分块内部点，除去边界
			  for (size_t i = 1; i < points.size(); i++)
			  {
				  for (size_t j = 1; j < points[i].size(); j++)
				  {
					  if (points[i][j].point.x()<rec.maxU && points[i][j].point.y()>rec.minV)//判断点在分块内部
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

			  //判断下边界上的点,将uvEdge中的整体编号赋予下边界上的点
			  for (int j = 0; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//确保点在分块内，//*****************************需变动
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
			  //判断左边界上的点,将uvEdge中的整体编号赋予左边界上的点,此处i=1，表明左下角点在上面已经编号了
			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//确保点在分块内//*****************************需变动
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

			  //从下边界开始，从第一列开始，此处j从1开始，j按照如下规则，分块0、2从0开始，分块1、3从1开始，这样可以防止极值点处的杆件重复
			  EdgePoints.push_back(uvEdge[0]);//将点0加入到边界点
			  for (int j = 1; j < points[0].size(); ++j)
			  {
				  if (points[0][j].point.x() < rec.maxU)//确保点在分块内，//*****************************需变动
				  {
					  if (points[0][points[0][j].next[0].y].Id == -1)//判断是否为最右点
					  {
						  if (points[1][points[0][j].next[1].y].Id != -1)//判断最右点的相连上点是否在内部，在为真
						  {
							  m = Add_Member(points[0][j], points[1][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[0][j], points[1][points[0][j].next[1].y], rec.maxU, rec.minV, 3);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[0][j], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
					  {
						  
							  m = Add_Member(points[0][j], points[points[0][j].next[1].x][points[0][j].next[1].y]);
							  mDisplayVec.push_back(m);
						  
					  }
				  }
			  }

			  //对内部点建立杆件，除去下边界、左边界
			  for (size_t i = 1; i < points.size(); i++)//从下往上
			  {
				  for (size_t j = points[i].size() - 1; j >= 1; j--)//从右往左
				  {
					  if (points[i][j].Id != -1)//找到内部点
					  {
						  if (points[i][j + 1].Id == -1)//点i，j的右点在外部
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i][j + 1], rec.maxU, rec.minV, 3);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

						  if (points[i + 1][points[i][j].next[1].y].Id == -1)//点i，j的上点在外部，此处没有考虑三角形情况
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][j], points[i + 1][points[i][j].next[1].y], rec.maxU, rec.minV, 3);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
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

			  //处理左边界,从第二行开始

			  for (int i = 1; i < points.size(); ++i)
			  {
				  if (points[i][0].point.y() > rec.minV)//确保点在分块内，//*****************************需变动
				  {
					  if (points[points[i][0].next[1].x][0].Id == -1)//判断是否为最上点
					  {
						  if (points[i][points[i][0].next[0].y].Id != -1)//判断最上点的相连右点是否在内部，在为真
						  {
							  m = Add_Member(points[i][0], points[i][points[i][0].next[0].y]);
							  mDisplayVec.push_back(m);//添加边界杆件
						  }
						  else//上点在外面
						  {
							  AddPoint2 = PointBtween_Rectangle(points[i][0], points[i][points[i][0].next[0].y], rec.maxU, rec.minV, 3);//生成新边界点
							  EdgePoints.push_back(AddPoint2);//将该点添加到边界点
							  addspace.setX(AddPoint2.point.x());
							  addspace.setY(AddPoint2.point.y());
							  addspace.setID(AddPoint2.Id);
							  nVec.push_back(addspace);
							  m = Add_Member(points[i][0], AddPoint2);
							  mDisplayVec.push_back(m);
						  }
					  }
					  else//其余内部点
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

//还没用到
//给出四边形对角的两点坐标和两边边长，确定另一顶点坐标,type表示寻找下一点的类型， 默认为type=0，表示从左下向右上寻找；
//type=1：从右下向左上寻找，type=2：从右上到左下，type =3：从左上到右下
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
	case 0://从左下往右上找点
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
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2较大，需要减小
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 > 0)//r2较小，需要增大
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1较大，需要减小
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 > 0)//r2较小，需要增大
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
	case 1://从右下向左上寻找
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
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength > 0)//r2较大，需要减小
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 > 0)//r2较小，需要增大
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1较小，需要减小u使r1增大
				   {
					   if (Uflag == 1)//说明上一步需要增大u，而此点需要减小u，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 < 0)//r1较大，需要增大u使r1减小
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
	case 2://从右上到左下
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
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2较小，需要减小v使r2增大
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 < 0)//r2较大，需要增大v使r2减小
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength < 0)//r1较小，需要减小u使r1增大
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 < 0)//r1较大，需要增大u是r1减小
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
	case 3://从左上到右下
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
				   r1 = (*tempNode1 - *tempNode3).GetLength();//计算两点间的距离
				   r2 = (*tempNode2 - *tempNode3).GetLength();
				   delete tempNode3;

				   if (r2 - vChrodLength < 0)//r2较小，需要减小v使r2增大
				   {

					   if (Vflag == 1)//说明上一步需要增大，而此步需要减小，可以断定要寻找的点在上一点和此点之间
					   {
						   v = v - deltaV;//退回到上一点
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
				   else if (vChrodLength - r2 < 0)//r2较大，需要增大v使r2减小
				   {
					   if (Vflag == -1)//说明上一步需要减小，而此步需要增大，可断定要寻找的点在上一点和此点之间
					   {
						   v = v + deltaV;//退回到上一点
						   deltaV = deltaV / 2;
						   v = v - deltaV;
						   Vflag = 0;
					   }
					   else
					   {
						   v = v + deltaV;//增大
						   Vflag = 1;
					   }

				   }

				   if (r1 - uChrodLength > 0)//r1较大，需要减小
				   {
					   if (Uflag == 1)//说明上一步需要增大，而此点需要减小，可以断定要寻找的点在上一点和此点之间
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
				   else if (uChrodLength - r1 > 0)//r2较小，需要增大
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


//根据两点的uv值求出两点之间的边界点，即该点的u或v为1,u_edge为边界u值,01表示上部两块，23表示下部两块
PointStruct NurbsSurface_qj::PointBtween_Rectangle(PointStruct p1, PointStruct p2, double u_edge, double v_edge,int type)
{
	PointStruct p3;

	switch (type)
	{
	case 0:
	case 1:
	{
			  if ((p1.point.x() - u_edge)*(p2.point.x() - u_edge) < 0)//如果是u超限
			  {
				  p3.point.x() = u_edge;
				  p3.point.y() = (u_edge - p1.point.x()) / (p2.point.x() - p1.point.x())*p2.point.y() +
					  (p2.point.x() - u_edge) / (p2.point.x() - p1.point.x())*p1.point.y();
				  if (p3.point.y() <= v_edge)//p3在内部
				  {
					  p3.Id = nId;
					  nId++;
					  break;
				  }
			  }
			  if ((p1.point.y() - v_edge)*(p2.point.y() - v_edge) < 0)//如果是v超限
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
			  if ((p1.point.x() - u_edge)*(p2.point.x() - u_edge) < 0)//如果是u超限
			  {
				  p3.point.x() = u_edge;
				  p3.point.y() = (u_edge - p1.point.x()) / (p2.point.x() - p1.point.x())*p2.point.y() +
					  (p2.point.x() - u_edge) / (p2.point.x() - p1.point.x())*p1.point.y();
				  if (p3.point.y() >= v_edge)//p3在内部
				  {
					  p3.Id = nId;
					  nId++;
					  break;
				  }
			  }
			  if ((p1.point.y() - v_edge)*(p2.point.y() - v_edge) < 0)//如果是v超限
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
















//******************************************取z最大值点，只分四块进行曲面划分****************************************************************

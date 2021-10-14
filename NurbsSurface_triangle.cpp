#include"NurbsSurface_triangle.h"

//比较函数
//bool compare_Node_u(Node p1, Node p2)
//{
//	return p1.getX() < p2.getX();
//}
//bool compare_Point2Dd_u(PLib::Point2Dd p1, PLib::Point2Dd p2)
//{
//	return p1.x() < p2.x();
//}

//构造函数
NurbsSurface_triangle::NurbsSurface_triangle()
{}
NurbsSurface_triangle::~NurbsSurface_triangle()
{}
NurbsSurface_triangle::NurbsSurface_triangle(int uDeg, int vDeg, const Vector_DOUBLE& uKnots, const Vector_DOUBLE& vKnots, const Matrix_HPoint3Dd& controlPoints, unsigned long long int n_Id, unsigned long long int m_Id)
{
	nId = n_Id;
	mId = m_Id;
	PlNurbsSurfaced tempNurbsSurface(uDeg, vDeg, uKnots, vKnots, controlPoints);
	nurbssurface = tempNurbsSurface;
}


//设定每份的长度，按此长度将u方向曲线等分
std::vector<Node> NurbsSurface_triangle::GetUPoints(double min, double max,double u0, double v0, double splitLengthU)
{
	std::vector<Node> vecleft;	
	std::vector<Node> vecright;
	std::vector<Node> vec;
	Node pointuv;
	pointuv.setX(u0) ;
	pointuv.setY(v0) ;
	while (pointuv.getX() <= max)//向右找等分点
	{
		vecright.push_back(pointuv);
		pointuv = GetnextUPoint(pointuv.getX(), v0, splitLengthU,0);
	}

	pointuv =GetnextUPoint(u0,v0,splitLengthU,1);	
	while (pointuv.getX() >= min)//向左找等分点
	{
		vecleft.push_back(pointuv);
		pointuv = GetnextUPoint(pointuv.getX(), v0, splitLengthU, 1);
	}

	if (!vecleft.empty())//当vecletf非空时
	{
		std::vector<Node>::iterator iter = vecleft.end();
		iter--;	
		for (;; iter--)
		{
			vec.push_back(*iter);
			if (iter == vecleft.begin())
			{
				break;
			}
		}		
	}
	(vec).insert(vec.end(), vecright.begin(), vecright.end());
	return vec;
}

//给定一点的u，v坐标，寻找u方向上弦长等于splitLengthU的u值
Node NurbsSurface_triangle::GetnextUPoint(double u, const double v, double splitLengthU,int type ,double error0)
{
	double deltaU = 0.025;//此处取值保证nextU在u+-deltaU*2附近振动
	Node pointuv;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	pointuv.setY(v);
	switch (type)
	{
	case 0://从左到右
		pointuv.setX( u + deltaU * 2);
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(pointuv.getX(), v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - splitLengthU > error0)//r较大，需要减小
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX()-  deltaU) ;
				flag = -1;
			}
			else if (splitLengthU - r > error0)//r较小，需要增加
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX() + deltaU);
				flag = 1;
			}

		} while (abs(r - splitLengthU) > error0);
		break;

	case 1://从右到左
		pointuv.setX( u - deltaU * 2);
		do
		{
			tempNode1 = GetSurfacePointatUandV(u, v);
			tempNode2 = GetSurfacePointatUandV(pointuv.getX(), v);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - splitLengthU > error0)//r较大，u需要增大
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX() + deltaU);
				flag = 1;
			}
			else if (splitLengthU - r > error0)//r较小，u需要减小
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX() - deltaU);
				flag = -1;
			}

		} while (abs(r - splitLengthU) > error0);
		break;
	default:
		break;
	}

	return pointuv;
}

//设定每段长度，按此长度将沿任意方向的曲线等分（除去平行于u方向的曲线），以v为自变量，u随之变化
std::vector<Node> NurbsSurface_triangle::GetPoints(double minU, double maxU, double minV, double maxV, double u0, double v0, double k,  double PointLength)
{
	std::vector<Node> Nodeuv;//存等分点的uv值
	std::vector<Node> Nodeup;
	std::vector<Node> Nodedown;
	Node nextNode;
	nextNode.setX(u0);
	nextNode.setY(v0);

	while (nextNode.getY() <= maxV && nextNode.getX() <= maxU && nextNode.getX() >= minU)//u0上方区域
	{
		Nodeup.push_back(nextNode);
		nextNode = GetNextPoint(nextNode,k, PointLength, 0);
	}

	nextNode.setX(u0);
	nextNode.setY(v0);
	nextNode = GetNextPoint(nextNode, k, PointLength, 1);
	while (nextNode.getY() >= minV && nextNode.getX() <= maxU && nextNode.getX() >= minU)//u0下方区域
	{
		Nodedown.push_back(nextNode);
		nextNode = GetNextPoint(nextNode, k,  PointLength, 1);
	}

	if (!Nodedown.empty())//当Nodedown非空
	{
		std::vector<Node>::iterator iter = Nodedown.end();
		iter--;
		if (iter != Nodedown.begin())//保证Nodedown不为空
		{
			for (;; iter--)
			{
				Nodeuv.push_back(*iter);
				if (iter == Nodedown.begin())
				{
					break;
				}
			}
		}		
	}
	(Nodeuv).insert(Nodeuv.end(), Nodeup.begin(), Nodeup.end());
	return Nodeuv;
}

//设定每段长度，按此长度将沿任意方向的曲线求下一点，以v为自变量，u随之变化（除去平行于u方向的曲线）,以v的正负分界
Node NurbsSurface_triangle::GetNextPoint(Node nextNode, double k, double PointLength, int type, double error0)
{
	double deltaV = 0.025;//v是自变量，u随其变化
	double deltaU = deltaV / k;	
	double nextU=0;
	double nextV=0;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	switch (type)
	{
	case 0://v为正
		nextU=nextNode.getX() + deltaU * 2;
		nextV=nextNode.getY() + deltaV * 2;
		do
		{
			tempNode1 = GetSurfacePointatUandV(nextNode.getX(), nextNode.getY());
			tempNode2 = GetSurfacePointatUandV(nextU, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - PointLength > error0)//r较大，v需要减小
			{
				if (flag == 1)
				{
					deltaV /= 2;
					deltaU /= 2;
				}
				nextU = nextU - deltaU;
				nextV = nextV - deltaV;
				flag = -1;
			}
			else if (PointLength - r > error0)//r较小，v需要增加
			{
				if (flag == -1)
				{
					deltaU /= 2;
					deltaV /= 2;
				}
				nextU = nextU + deltaU;
				nextV = nextV + deltaV;
				flag = 1;
			}

		} while (abs(r - PointLength) > error0);
		break;
	case 1://v为负
		nextU = nextNode.getX() - deltaU * 2;
		nextV = nextNode.getY() - deltaV * 2;
		do
		{
			tempNode1 = GetSurfacePointatUandV(nextNode.getX(), nextNode.getY());
			tempNode2 = GetSurfacePointatUandV(nextU, nextV);

			r = sqrt((tempNode1->getX() - tempNode2->getX()) *  (tempNode1->getX() - tempNode2->getX()) +
				(tempNode1->getY() - tempNode2->getY()) *  (tempNode1->getY() - tempNode2->getY()) +
				(tempNode1->getZ() - tempNode2->getZ()) *  (tempNode1->getZ() - tempNode2->getZ()));
			delete tempNode1;
			delete tempNode2;
			if (r - PointLength > error0)//r较大，v需要减小
			{
				if (flag == 1)
				{
					deltaV /= 2;
					deltaU /= 2;
				}
				nextV = nextV + deltaV;
				nextU = nextU + deltaU;				
				flag = -1;
			}
			else if (PointLength - r > error0)//r较小，v需要增加
			{
				if (flag == -1)
				{
					deltaU /= 2;
					deltaV /= 2;
				}
				nextV = nextV - deltaV;
				nextU = nextU - deltaU;				
				flag = 1;
			}
		} while (abs(r - PointLength) > error0);
		break;
	default:
		break;
	}
	nextNode.setX(nextU);
	nextNode.setY(nextV);
	return nextNode;
}

//根据uv得到三维坐标
Node* NurbsSurface_triangle::GetSurfacePointatUandV(double u, double v)//给定曲面和比例u，v,求曲面上的点
{
	PLib::HPoint3Dd hpoint = this->nurbssurface.pointAt(u, v);
	PLib::Point3Dd point = project(hpoint);

	Node* thepoint = new Node();
	thepoint->setX((point.data[0]));
	thepoint->setY((point.data[1]));
	thepoint->setZ((point.data[2]));

	return thepoint;
}

//寻找曲面的极大值点
std::vector<PLib::Point2Dd> NurbsSurface_triangle::GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV)//寻找曲面的极大值点
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
		u = minU + i * deltaU;

		for (int j = 0; j < 50; ++j)
		{
			v = minV + j * deltaV;
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
std::vector<PLib::Point2Dd> NurbsSurface_triangle::GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV)//寻找曲面的极大值点
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
		u = minU + i * deltaU;

		for (int j = 0; j < 50; ++j)
		{
			v = minV + j * deltaV;
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

//主程序，等分曲面
std::vector<Member> NurbsSurface_triangle::SurfaceAverageSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double PointLength, std::vector< Node> & nVec)
{
	//找曲面z绝对值最大点
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极小值点
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//寻找曲面内的极大值点
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//存储极值点
	int num=0;
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
	for (int i=1; i<num; i++)//寻找曲面z最大值	
	{
		tempNode1 = GetSurfacePointatUandV(extrmums[i].x, extrmums[i].y);
		if (fabs(tempNode1->getZ()) > MaxzofPoint)
		{
			MaxzofPoint = fabs(tempNode1->getZ());
			maxPoint = extrmums[i];
		}
		delete tempNode1;
	}

	std::vector<Node> LineofSixty;//60度角的直线上的点（u，v）
	LineofSixty = GetPoints(minU, maxU, minV, maxV, maxPoint.x, maxPoint.y, 2, PointLength);
	//std::sort(LineofSixty.begin(), LineofSixty.end(), compare_Node_u);//使60度线上点按u从小到大排序

	std::vector<Node> LineOne;//存储u方向直线上的等分点
	std::vector<Node> LineTwo;
	std::vector<Member> mVec;//存储全部杆件
	std::vector<Member> mofLine;//一条线的杆件
	Node addspace;
	Member m;

	LineOne = GetUPoints(minU, maxU, LineofSixty.begin()->getX(), LineofSixty.begin()->getY(), splitLengthU);
	//std::sort(LineOne.begin(), LineOne.end(),compare_Node_u);//从小到大排序
	for (std::vector<Node>::iterator iter = LineOne.begin(); iter != LineOne.end(); iter++)//对lingOne整体编号
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineOne.begin(); iter != LineOne.end() - 1; iter++) //连第一行杆
	{
		m.setInode(iter->getID());
		m.setJnode((iter+1)->getID());
		m.set_ID(mId);
		mVec.push_back(m);
		mId++;
	}

	for (std::vector<Node>::iterator iter=LineofSixty.begin()+1; iter!=LineofSixty.end(); iter++)
	{
		LineTwo = GetUPoints(minU, maxU,iter->getX(), iter->getY(), splitLengthU);
		//std::sort(LineTwo.begin(), LineTwo.end(), compare_Node_u);
		//连杆************
		mofLine = MemberofLine2(LineOne, LineTwo, nVec);
		(mVec).insert(mVec.end(), mofLine.begin(), mofLine.end());
		//连杆************	
		if (iter == LineofSixty.end() - 1)//最后一条线
		{
			break;
		}
		LineOne = LineTwo;
		LineTwo = GetUPoints(minU, maxU, (iter+1)->getX(), (iter+1)->getY(), splitLengthU);
	}
	
	return mVec;
}

//连杆，对点、杆进行整体编号,判断区间
std::vector<Member> NurbsSurface_triangle::MemberofLine1(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec)
{
	Member m;//用于增加edge_m
	std::vector<Member> mofLine;//存储边界杆件
	int up_index = -1;//lineOne中某点在LineTwo中的右区间点的Id		
	Node addspace;

	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end(); iter++)//整体编号
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end() - 1; iter++) //连第二行杆
	{
		m.setInode(iter->getID());
		m.setJnode((iter+1)->getID());
		m.set_ID(mId);
		mofLine.push_back(m);
		mId++;
	}

	for (size_t i = 0; i < LineOne.size(); i++)//遍历LineOne
	{
		up_index = -1;
		for (size_t j = 0; j < LineTwo.size(); j++)
		{
			if (LineTwo[j].getX()-LineOne[i].getX() > 0)//lineOne中某点在LineTwo中的右区间点
			{
				up_index = j;//将右边界点的下标传出
				break;
			}
		}
		if (up_index >= 1 && up_index <= LineTwo.size())//内部点有左、右边界点与之相连
		{
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[up_index - 1].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;
			m.setJnode(LineTwo[up_index].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;			
			continue;
		}
		if (up_index == 0)//内部点只有右边界点与之相连
		{
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[up_index].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;
			continue;
		}
		if (up_index == -1)//内部点有左边界点与之相连
		{	
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[LineTwo.size() - 1].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;			
		}
	}

	return mofLine;
}

//连杆，对点、杆进行整体编号，不判断区间
std::vector<Member> NurbsSurface_triangle::MemberofLine2(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec)
{
	Member m;//用于增加edge_m
	std::vector<Member> mofLine;//存储边界杆件
	int up_index = -1;//lineOne中某点在LineTwo中的右区间点的Id		
	Node addspace;

	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end(); iter++)//整体编号
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end() - 1; iter++) //连第二行杆
	{
		m.setInode(iter->getID());
		m.setJnode((iter + 1)->getID());
		m.set_ID(mId);
		mofLine.push_back(m);
		mId++;
	}

	
	for (size_t i = 0; i < LineTwo.size(); i++)//遍历LineOne,up记录LineTwo的右点
	{
		if (LineOne[0].getX() < LineTwo[i].getX())
		{
			up_index = i;
			break;
		}
	}	

	for (size_t i = 0; i < LineOne.size(); i++)
	{
		if ((i+up_index-1>=0) && (i+up_index-1<LineTwo.size()))
		{
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[i+up_index  - 1].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;
		}
		if ((i + up_index >= 0) && (i + up_index  < LineTwo.size()))
		{
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[i+up_index].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;
		}
	}
	

	/*if (up_index >= 1)
	{
		for (size_t i = 0; i < LineOne.size(); i++)
		{			
			if (i +up_index-1< LineTwo.size())
			{
				m.setInode(LineOne[i].getID());
				m.setJnode(LineTwo[i+up_index-1].getID());
				m.set_ID(mId);
				mofLine.push_back(m);
				mId++;				
			}
			if (i+up_index<LineTwo.size())
			{
				m.setInode(LineOne[i].getID());
				m.setJnode(LineTwo[i + up_index].getID());
				m.set_ID(mId);
				mofLine.push_back(m);
				mId++;
			}
		}
	}*/
	return mofLine;
}
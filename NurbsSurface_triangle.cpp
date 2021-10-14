#include"NurbsSurface_triangle.h"

//�ȽϺ���
//bool compare_Node_u(Node p1, Node p2)
//{
//	return p1.getX() < p2.getX();
//}
//bool compare_Point2Dd_u(PLib::Point2Dd p1, PLib::Point2Dd p2)
//{
//	return p1.x() < p2.x();
//}

//���캯��
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


//�趨ÿ�ݵĳ��ȣ����˳��Ƚ�u�������ߵȷ�
std::vector<Node> NurbsSurface_triangle::GetUPoints(double min, double max,double u0, double v0, double splitLengthU)
{
	std::vector<Node> vecleft;	
	std::vector<Node> vecright;
	std::vector<Node> vec;
	Node pointuv;
	pointuv.setX(u0) ;
	pointuv.setY(v0) ;
	while (pointuv.getX() <= max)//�����ҵȷֵ�
	{
		vecright.push_back(pointuv);
		pointuv = GetnextUPoint(pointuv.getX(), v0, splitLengthU,0);
	}

	pointuv =GetnextUPoint(u0,v0,splitLengthU,1);	
	while (pointuv.getX() >= min)//�����ҵȷֵ�
	{
		vecleft.push_back(pointuv);
		pointuv = GetnextUPoint(pointuv.getX(), v0, splitLengthU, 1);
	}

	if (!vecleft.empty())//��vecletf�ǿ�ʱ
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

//����һ���u��v���꣬Ѱ��u�������ҳ�����splitLengthU��uֵ
Node NurbsSurface_triangle::GetnextUPoint(double u, const double v, double splitLengthU,int type ,double error0)
{
	double deltaU = 0.025;//�˴�ȡֵ��֤nextU��u+-deltaU*2������
	Node pointuv;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	pointuv.setY(v);
	switch (type)
	{
	case 0://������
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
			if (r - splitLengthU > error0)//r�ϴ���Ҫ��С
			{
				if (flag == 1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX()-  deltaU) ;
				flag = -1;
			}
			else if (splitLengthU - r > error0)//r��С����Ҫ����
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

	case 1://���ҵ���
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
			if (r - splitLengthU > error0)//r�ϴ�u��Ҫ����
			{
				if (flag == -1)
				{
					deltaU /= 2;
				}
				pointuv.setX(pointuv.getX() + deltaU);
				flag = 1;
			}
			else if (splitLengthU - r > error0)//r��С��u��Ҫ��С
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

//�趨ÿ�γ��ȣ����˳��Ƚ������ⷽ������ߵȷ֣���ȥƽ����u��������ߣ�����vΪ�Ա�����u��֮�仯
std::vector<Node> NurbsSurface_triangle::GetPoints(double minU, double maxU, double minV, double maxV, double u0, double v0, double k,  double PointLength)
{
	std::vector<Node> Nodeuv;//��ȷֵ��uvֵ
	std::vector<Node> Nodeup;
	std::vector<Node> Nodedown;
	Node nextNode;
	nextNode.setX(u0);
	nextNode.setY(v0);

	while (nextNode.getY() <= maxV && nextNode.getX() <= maxU && nextNode.getX() >= minU)//u0�Ϸ�����
	{
		Nodeup.push_back(nextNode);
		nextNode = GetNextPoint(nextNode,k, PointLength, 0);
	}

	nextNode.setX(u0);
	nextNode.setY(v0);
	nextNode = GetNextPoint(nextNode, k, PointLength, 1);
	while (nextNode.getY() >= minV && nextNode.getX() <= maxU && nextNode.getX() >= minU)//u0�·�����
	{
		Nodedown.push_back(nextNode);
		nextNode = GetNextPoint(nextNode, k,  PointLength, 1);
	}

	if (!Nodedown.empty())//��Nodedown�ǿ�
	{
		std::vector<Node>::iterator iter = Nodedown.end();
		iter--;
		if (iter != Nodedown.begin())//��֤Nodedown��Ϊ��
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

//�趨ÿ�γ��ȣ����˳��Ƚ������ⷽ�����������һ�㣬��vΪ�Ա�����u��֮�仯����ȥƽ����u��������ߣ�,��v�������ֽ�
Node NurbsSurface_triangle::GetNextPoint(Node nextNode, double k, double PointLength, int type, double error0)
{
	double deltaV = 0.025;//v���Ա�����u����仯
	double deltaU = deltaV / k;	
	double nextU=0;
	double nextV=0;
	double r = 0;
	int flag = 0;
	Node * tempNode1;
	Node * tempNode2;
	switch (type)
	{
	case 0://vΪ��
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
			if (r - PointLength > error0)//r�ϴ�v��Ҫ��С
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
			else if (PointLength - r > error0)//r��С��v��Ҫ����
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
	case 1://vΪ��
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
			if (r - PointLength > error0)//r�ϴ�v��Ҫ��С
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
			else if (PointLength - r > error0)//r��С��v��Ҫ����
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

//����uv�õ���ά����
Node* NurbsSurface_triangle::GetSurfacePointatUandV(double u, double v)//��������ͱ���u��v,�������ϵĵ�
{
	PLib::HPoint3Dd hpoint = this->nurbssurface.pointAt(u, v);
	PLib::Point3Dd point = project(hpoint);

	Node* thepoint = new Node();
	thepoint->setX((point.data[0]));
	thepoint->setY((point.data[1]));
	thepoint->setZ((point.data[2]));

	return thepoint;
}

//Ѱ������ļ���ֵ��
std::vector<PLib::Point2Dd> NurbsSurface_triangle::GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV)//Ѱ������ļ���ֵ��
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
std::vector<PLib::Point2Dd> NurbsSurface_triangle::GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV)//Ѱ������ļ���ֵ��
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

//�����򣬵ȷ�����
std::vector<Member> NurbsSurface_triangle::SurfaceAverageSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double PointLength, std::vector< Node> & nVec)
{
	//������z����ֵ����
	std::vector<PLib::Point2Dd> pMinExtrmums = GetSurfaceMinExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ�Сֵ��
	std::vector<PLib::Point2Dd> pMaxExtrmums = GetSurfaceMaxExtrmums(minU, maxU, minV, maxV);//Ѱ�������ڵļ���ֵ��
	CPoint2D * extrmums = new CPoint2D[pMaxExtrmums.size() + pMinExtrmums.size()];//�洢��ֵ��
	int num=0;
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
	for (int i=1; i<num; i++)//Ѱ������z���ֵ	
	{
		tempNode1 = GetSurfacePointatUandV(extrmums[i].x, extrmums[i].y);
		if (fabs(tempNode1->getZ()) > MaxzofPoint)
		{
			MaxzofPoint = fabs(tempNode1->getZ());
			maxPoint = extrmums[i];
		}
		delete tempNode1;
	}

	std::vector<Node> LineofSixty;//60�Ƚǵ�ֱ���ϵĵ㣨u��v��
	LineofSixty = GetPoints(minU, maxU, minV, maxV, maxPoint.x, maxPoint.y, 2, PointLength);
	//std::sort(LineofSixty.begin(), LineofSixty.end(), compare_Node_u);//ʹ60�����ϵ㰴u��С��������

	std::vector<Node> LineOne;//�洢u����ֱ���ϵĵȷֵ�
	std::vector<Node> LineTwo;
	std::vector<Member> mVec;//�洢ȫ���˼�
	std::vector<Member> mofLine;//һ���ߵĸ˼�
	Node addspace;
	Member m;

	LineOne = GetUPoints(minU, maxU, LineofSixty.begin()->getX(), LineofSixty.begin()->getY(), splitLengthU);
	//std::sort(LineOne.begin(), LineOne.end(),compare_Node_u);//��С��������
	for (std::vector<Node>::iterator iter = LineOne.begin(); iter != LineOne.end(); iter++)//��lingOne������
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineOne.begin(); iter != LineOne.end() - 1; iter++) //����һ�и�
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
		//����************
		mofLine = MemberofLine2(LineOne, LineTwo, nVec);
		(mVec).insert(mVec.end(), mofLine.begin(), mofLine.end());
		//����************	
		if (iter == LineofSixty.end() - 1)//���һ����
		{
			break;
		}
		LineOne = LineTwo;
		LineTwo = GetUPoints(minU, maxU, (iter+1)->getX(), (iter+1)->getY(), splitLengthU);
	}
	
	return mVec;
}

//���ˣ��Ե㡢�˽���������,�ж�����
std::vector<Member> NurbsSurface_triangle::MemberofLine1(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec)
{
	Member m;//��������edge_m
	std::vector<Member> mofLine;//�洢�߽�˼�
	int up_index = -1;//lineOne��ĳ����LineTwo�е���������Id		
	Node addspace;

	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end(); iter++)//������
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end() - 1; iter++) //���ڶ��и�
	{
		m.setInode(iter->getID());
		m.setJnode((iter+1)->getID());
		m.set_ID(mId);
		mofLine.push_back(m);
		mId++;
	}

	for (size_t i = 0; i < LineOne.size(); i++)//����LineOne
	{
		up_index = -1;
		for (size_t j = 0; j < LineTwo.size(); j++)
		{
			if (LineTwo[j].getX()-LineOne[i].getX() > 0)//lineOne��ĳ����LineTwo�е��������
			{
				up_index = j;//���ұ߽����±괫��
				break;
			}
		}
		if (up_index >= 1 && up_index <= LineTwo.size())//�ڲ��������ұ߽����֮����
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
		if (up_index == 0)//�ڲ���ֻ���ұ߽����֮����
		{
			m.setInode(LineOne[i].getID());
			m.setJnode(LineTwo[up_index].getID());
			m.set_ID(mId);
			mofLine.push_back(m);
			mId++;
			continue;
		}
		if (up_index == -1)//�ڲ�������߽����֮����
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

//���ˣ��Ե㡢�˽��������ţ����ж�����
std::vector<Member> NurbsSurface_triangle::MemberofLine2(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec)
{
	Member m;//��������edge_m
	std::vector<Member> mofLine;//�洢�߽�˼�
	int up_index = -1;//lineOne��ĳ����LineTwo�е���������Id		
	Node addspace;

	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end(); iter++)//������
	{
		addspace.setX(iter->getX());
		addspace.setY(iter->getY());
		addspace.setID(nId);
		iter->setID(nId);
		nVec.push_back(addspace);
		nId++;
	}
	for (std::vector<Node>::iterator iter = LineTwo.begin(); iter != LineTwo.end() - 1; iter++) //���ڶ��и�
	{
		m.setInode(iter->getID());
		m.setJnode((iter + 1)->getID());
		m.set_ID(mId);
		mofLine.push_back(m);
		mId++;
	}

	
	for (size_t i = 0; i < LineTwo.size(); i++)//����LineOne,up��¼LineTwo���ҵ�
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
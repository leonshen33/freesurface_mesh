#pragma once

#include <iostream>
#include <assert.h>
#include "nurbs.h"// for nurbs++
#include <gl/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <vector>
#include "Node.h"
#include "nurbsS.h"
#include <fstream>
#include <algorithm>
#include "Member.h"
//#include<algorithm>
//#include<functional>

#include <Windows.h>
//using namespace std;

//#pragma comment(linker, "/subsystem:console") 
#pragma comment( lib, "glew32.lib")
#pragma comment( lib, "glu32.lib")
#pragma comment( lib, "glut32.lib")
#pragma comment( lib, "opengl32.lib")
//using namespace PLib;


class NurbsSurface_triangle
{
public:
	NurbsSurface_triangle(void);
	NurbsSurface_triangle(int  uDeg, int vDeg, const Vector_DOUBLE& uKnots, const Vector_DOUBLE& vKnots, const Matrix_HPoint3Dd& controlPoints, unsigned long long int n_Id = 0, unsigned long long int m_Id = 0);
	~NurbsSurface_triangle(void);

	std::vector<Node> GetUPoints(double minU, double maxU, double u0, double v0, double splitLengthU);//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�u�������ߵȷ�	
	Node GetnextUPoint(double u, const double v, double splitLengthU, int type=0,double error0 = 1e-4);//����һ���u��v���꣬Ѱ��u�������ҳ�����splitLengthU��uֵ

	//�趨ÿ�γ��ȣ����˳��Ƚ������ⷽ������ߵȷ֣���ȥƽ����u��������ߣ�����vΪ�Ա�����u��֮�仯
	//(u0,v0),k�ǹ����ֱ�ߵ�б�ʣ�au-v=au0-v0
	std::vector<Node> GetPoints(double minU, double maxU, double minV, double maxV, double u0, double v0, double k, double PointLength);
	//�趨ÿ�γ��ȣ����˳��Ƚ������ⷽ�����������һ�㣬��vΪ�Ա�����u��֮�仯����ȥƽ����u��������ߣ�
	Node GetNextPoint(Node nextNode,  double k,  double PointLength,int type=0,double error0=1e-4);
	
	std::vector<PLib::Point2Dd> GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV);//Ѱ������ļ���ֵ��
	std::vector<PLib::Point2Dd> GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV);//Ѱ������ļ�Сֵ��
	
	//�õ�(U,V)��Ӧ�����ά����
	Node* GetSurfacePointatUandV(double u, double v);

	//�����򣬵ȷ�����
	std::vector< Member> SurfaceAverageSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double PointLength, std::vector< Node> & nVec);

	//���ˣ��Ե㡢�˽��������ţ��ж�����
	std::vector<Member> MemberofLine1(std::vector<Node> &LineOne,std::vector<Node> &LineTwo,std::vector<Node> &nVec);
	//���ˣ��Ե㡢�˽��������ţ����ж�����
	std::vector<Member> MemberofLine2(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec);
	


private:
	unsigned long long int nId;
	unsigned long long int mId;

	PlNurbsSurfaced nurbssurface;

};
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

	std::vector<Node> GetUPoints(double minU, double maxU, double u0, double v0, double splitLengthU);//设定每份的长度，按一定长度将u方向曲线等分	
	Node GetnextUPoint(double u, const double v, double splitLengthU, int type=0,double error0 = 1e-4);//给定一点的u，v坐标，寻找u方向上弦长等于splitLengthU的u值

	//设定每段长度，按此长度将沿任意方向的曲线等分（除去平行于u方向的曲线），以v为自变量，u随之变化
	//(u0,v0),k是过点的直线的斜率，au-v=au0-v0
	std::vector<Node> GetPoints(double minU, double maxU, double minV, double maxV, double u0, double v0, double k, double PointLength);
	//设定每段长度，按此长度将沿任意方向的曲线求下一点，以v为自变量，u随之变化（除去平行于u方向的曲线）
	Node GetNextPoint(Node nextNode,  double k,  double PointLength,int type=0,double error0=1e-4);
	
	std::vector<PLib::Point2Dd> GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV);//寻找曲面的极大值点
	std::vector<PLib::Point2Dd> GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV);//寻找曲面的极小值点
	
	//得到(U,V)对应点的三维坐标
	Node* GetSurfacePointatUandV(double u, double v);

	//主程序，等分曲面
	std::vector< Member> SurfaceAverageSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double PointLength, std::vector< Node> & nVec);

	//连杆，对点、杆进行整体编号，判断区间
	std::vector<Member> MemberofLine1(std::vector<Node> &LineOne,std::vector<Node> &LineTwo,std::vector<Node> &nVec);
	//连杆，对点、杆进行整体编号，不判断区间
	std::vector<Member> MemberofLine2(std::vector<Node> &LineOne, std::vector<Node> &LineTwo, std::vector<Node> &nVec);
	


private:
	unsigned long long int nId;
	unsigned long long int mId;

	PlNurbsSurfaced nurbssurface;

};
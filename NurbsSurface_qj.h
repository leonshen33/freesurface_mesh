#pragma once

#include <iostream>
//#include <assert.h>
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
#include "Rectangle.h"

#include<functional>

#include <Windows.h>
//using namespace std;

#pragma comment(linker, "/subsystem:console") 
#pragma comment( lib, "glew32.lib")
#pragma comment( lib, "glu32.lib")
#pragma comment( lib, "glut32.lib")
#pragma comment( lib, "opengl32.lib")
//using namespace PLib;

struct equalDivPoints//记录边界和标识符
{
	std::vector<double> points;
	double value; 
	
};

struct pointCoordinate//定义
{
	int x;
	int y;
	pointCoordinate(int x0, int y0) :x(x0), y(y0){};
};
struct PointStruct
{
	PLib::Point2Dd point;//记录u，v坐标
	std::vector<pointCoordinate> next;//记录该点所指向的点坐标
	 long long int Id;//************************************************************与pId相同的局部编号,初始为-1
};

struct pointWithNodeId
{
	PLib::Point2Dd point;
	unsigned long long int nodeId;
};


class NurbsSurface_qj
{
public:
	NurbsSurface_qj(void);
	//NurbsSurface_qj(Vector_HPoint3Dd ccontrol_points,Vector_DOUBLE cknots,int cdegree,Vector_HPoint3Dd tcontrol_points,Vector_DOUBLE tknots,int tdegree);//正算
	NurbsSurface_qj(Vector_HPoint3Dd ccurve_points,int cdegree,Vector_HPoint3Dd tcurve_points,int tdegree);
	NurbsSurface_qj(std::vector<Node*> Surface_points/*,int cdegree,int tdegree*/);
	/*NurbsSurface_qj(std::vector<Node*> ccontrol_points,std::vector<double > cknots,int cdegree,std::vector<Node*> tcontrol_points,std::vector<double > tknots,int tdegree);*/
	NurbsSurface_qj(std::vector<Node*> ccurve_points,int cdegree,std::vector<Node*> tcurve_points,int tdegree);
	NurbsSurface_qj(int  uDeg, int vDeg, const Vector_DOUBLE& uKnots, const Vector_DOUBLE& vKnots, const Matrix_HPoint3Dd& controlPoints, unsigned long long int n_Id, unsigned long long int m_Id);
	//NurbsSurface(int DegU, int DegV, const Vector<T>& Uk, const Vector<T>& Vk, const Matrix< HPoint_nD<T,N> >& Cp) ;
	
	void reset(Vector_HPoint3Dd ccurve_points,int cdegree,Vector_HPoint3Dd tcurve_points,int tdegree);

	void reset(Vector_HPoint3Dd Surface_points);
	
	Node* GetSurfacePointatUandV(double u,double v);

	int getControlPoint(std::vector<Node*>& controlpoint);
	
	
	int getKnotUPoint(double* Knotpoint);
	int getKnotVPoint(double* Knotpoint);


	/****************getExtremums*************************
	***函数功能：获取一条曲线上所有极值点坐标
	***参数：min为曲线对应的最小的v坐标，max为曲线对应的最大v坐标；
	***返回值：所有曲线上极值点所对应的u坐标
	*****************getExtremums*************************/
	std::vector<double> getVExtremums(double min, double max, const double u);//u值固定，沿v方向寻找极值点，
	std::vector<double> getUExtremums(double min, double max, const double u);//v值固定，沿u方向寻找极值点，


	/****************getUInflectionPoints*************************
	***函数功能：给定v方向的值，获取U方向曲线上的拐点；
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：曲线上拐点所对应的u值坐标
	***思路：拐点两侧的z方向的二阶导异号，设定较大步长，逐个遍历，若前一个点的z方向
	***      二阶导的符号与该点的符号异号，且该点与前一点的距离很短，即可认为该点为拐点
	*****************getInflectionPoints*************************/
	std::vector<double> getUInflectionPoints(double min, double max, const double v);

	/****************getVInflectionPoints*************************
	***函数功能：给定u方向值，获取V方向曲线上的拐点；
	***参数：min为曲线对应的最小的v坐标，max为曲线对应的最大v坐标；
	***返回值：曲线上拐点所对应的v值坐标
	***思路：拐点两侧的z方向的二阶导异号，设定较大步长，逐个遍历，若前一个点的z方向
	***      二阶导的符号与该点的符号异号，且该点与前一点的距离很短，即可认为该点为拐点
	*****************getInflectionPoints*************************/
	std::vector<double> getVInflectionPoints(double min, double max, const double u);

	/***********************getUParallelPoints********************
	***函数功能：给定v方向值，获取u方向曲线上的导数向量与始终点连线平行的点的u值；
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：曲线上拐点所对应的u值坐标
	***思路：先寻找过起点和终点以及起始点与终点的中点的平面l, 找平面l的法向量n，然后沿曲线寻找与n
	***      垂直的点对应的u值
	***********************getUParallelPoints********************/
	std::vector<double> getUParallelPoints(double min, double max, const double v);

	/***********************getVParallelPoints********************
	***函数功能：量与始给定u方向值，获取v方向曲线上的导数向终点连线平行的点的v值；
	***参数：min为曲线对应的最小的v坐标，max为曲线对应的最大v坐标；
	***返回值：曲线上拐点所对应的u值坐标
	***思路：先寻找过起点和终点以及起始点与终点的中点的平面l, 找平面l的法向量n，然后沿曲线寻找与n
	***      垂直的点对应的v值
	***********************getVParallelPoints********************/
	std::vector<double> getVParallelPoints(double min, double max, const double u);

	std::vector<double> getUAveragePerSection(double min, double max,const double v, int splitNum);//v值固定，将u方向曲线按等弦长分割
	std::vector<double> getVAveragePersection(double min, double max,const double u, int splitNum);//u值固定，将u方向曲线按等弦长分割



	
	std::vector<double> getUAverage(double min, double max, const double v, int splitNum);//设定均分多少份，v值固定，将u方向曲线按等弦长分割
	std::vector<double> getVAverage(double min, double max, const double u, int splitNum);//设定均分多少份，u值固定，将v方向曲线按等弦长分割

	std::vector<double> getUAverageLength(double min, double max, const double v, double lengthPerSection, int type = 0);//设定每份的长度，按一定长度将u方向曲线等分
	std::vector<double> getVAverageLength(double min, double max, const double u, double lengthPerSection, int type = 0);//设定每份的长度，按一定长度将v方向曲线等分
	
	std::vector<PLib::Point2Dd> GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV);//寻找曲面的极大值点
	std::vector<PLib::Point2Dd> GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV);//寻找曲面的极小值点

	
	//将曲面按弦长等分
	//type表示寻找下一点的类型， 默认为type = 0，表示从左下向右上寻找；type=1：从右下向左上寻找，type=2：从右上到左下，type =3：从左上到右下		
	//std::vector<Member> surfaceAverageSpit(double minU, double realMaxU, double minV, double realMaxV, double splitLengthU, double splitLengthV, std::vector<Node>& pDisplayVec, int type=0);
	std::vector<Member> surfaceAverageSpit(rectangleRegion &rec, double splitLengthU, double splitLengthV, std::vector<PointStruct> &uvEdge);
	std::vector< Member> surfaceAverageTotalSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec);
	//将nurbs曲面三角等分，并将等分构成的点 和杆分别存储在nVec和返回值中
	
	//std::vector<Member> surfaceTriSplit(double minU, double maxU, double minV, double maxV, int uSplitNum, int vSplitNum, std::vector<Node> & nVec);
	
	double nextUPoint(double u, const double v, double uChordLength, int type = 0, double error = 1e-4);//给定一点的u，v坐标，寻找u方向上弦长等于chordLength的u值
	double nextVPoint(const double u, double v, double vChordLength, int type = 0, double error = 1e-4);//给定一点的u，v坐标，寻找v方向上弦长等于chordLength的v值
	//给出四边形对角的两点坐标和两边边长，确定另一顶点坐标,type表示寻找下一点的类型， 默认为type=0，表示从左下向右上寻找；
	//type=1：从右下向左上寻找，type=2：从右上到左下，type =3：从左上到右下
	PLib::Point2Dd nextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type = 0, double error = 1e-4);
	void drawNurbsSurface();

	~NurbsSurface_qj(void);
	//*********************************新添加的public成员**********************************
	std::vector< Member> SurfaceAverageSpit_Rectangle(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec);
	std::vector<Member> MemberofRec_Rectangle(rectangleRegion rec, double splitLengthU, double splitLengthV, std::vector<Node> &nVec,std::vector<PointStruct> &uvEdge/*,std::vector<bool> &Gui_index*/);//分块画杆件
	
	//内部边界处理
	std::vector<Member> InsideEdge_Rectangle(double minU, double maxU, double minV, double maxV, CPoint2D maxPoint, double splitLengthU, double splitLengthV, std::vector< Node> & nVec, std::vector<PointStruct> &uvEdge);
	//根据两点的uv值求出两点之间的边界点，即该点的u或v为1,u_edge为边界u值,01表示上部两块，23表示下部两块
	PointStruct PointBtween_Rectangle(PointStruct p1, PointStruct p2,double u_edge,double v_edge,int type);


	//给出四边形对角的两点坐标和两边边长，确定另一顶点坐标,type表示寻找下一点的类型， 默认为type=0，表示从左下向右上寻找；
	//type=1：从右下向左上寻找，type=2：从右上到左下，type =3：从左上到右下
	PLib::Point2Dd NextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type = 0, double error = 1e-4);

	struct IdRelationship//记录四边形四点Id
	{
		unsigned long long int IdOne;//1ID
		unsigned long long int IdTwo;//2ID
		unsigned long long int IdThree;//3ID
		unsigned long long int IdFour;//4ID
	};

	std::vector<Node> UVForOut;//记录点的uv值，整体编号

	/*map(NodeForOut的下标，NodeForOut的ID）*/

	
	//std::vector<IdRelationship> IdRelShip;//点相互关系std::vector
	void BetweenSearch(std::vector<PointStruct> uv_edge, double key1, double key2,int xy,int &up_index, int &down_index);//判断靠近边界的点界于边界点的区间，返回up_index,down_index
	Member Add_Member(PointStruct point_in,PointStruct point_edge);//添加杆件，point_in第一点（内部点），point_edge第二点（不规则边界点）
	Member Add_2Member(PointStruct point_in, PointStruct point_edge);//添加杆件，point_in第一点（内部点），point_edge第二点（规则边界点）

	//返回分块的上下左右边界点的容器，并对边界点进行整体编号，将边界点加入到整体点中，&用于修改uvEdge的Id
	void Edge_udlr(std::vector<PointStruct> &uvEdge, rectangleRegion rec, std::vector<PointStruct> &uvedge_up, std::vector<PointStruct> &uvedge_down,std::vector<PointStruct> &uvedge_left,std::vector<PointStruct> &uvedge_right);
	
	//边界杆件处理(全部边界点容器用于增加边界点【当内部点位于两边界点中间且与边界很近时】，某一边界的点，靠近边界的内部点,uv方向,）
	//只考虑在一个边界区间内至多有两个内部点，不考虑三个内部点的情况（此情况基本不可能）
	std::vector<Member> Edge_Member(std::vector<PointStruct> &uvEdge, std::vector<PointStruct> uvedge, std::vector<PointStruct*> pointsedge,  int xy,double length_differ = 2);

	//处理边界自身
	std::vector<Member> EdgeTotal_Member(std::vector<PointStruct> uvEdge, std::vector<PLib::Point2Dd> extrmumVec, int num, double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV);

	//处理边界的角点
	std::vector<Member> EdgeCorner_Member(std::vector<PointStruct*> &point_corner, double length_differ = 2);

	//*********************************新添加的public成员**********************************


private:
	Vector_HPoint3Dd TControlpoint;//nurbs曲线T的控制点
	Vector_HPoint3Dd CControlpoint;//nurbs曲线C的控制点
	
	Vector_HPoint3Dd TCurvepoint;//nurbs曲线T的采样点
	Vector_HPoint3Dd CCurvepoint;//nurbs曲线C的采样点

	Vector_DOUBLE TKnotpoint;//nurbs曲线T的节点
	Vector_DOUBLE CKnotpoint;//nurbs曲线C的节点
	
	int TDegree;//nurbs曲线T的度
	int CDegree;//nurbs曲线C的度
	
	PlNurbsCurved tnurbscurve ;//nurbs曲线T
	PlNurbsCurved cnurbscurve ;//nurbs曲线C

	PlNurbsSurfaced nurbssurface;


	Vector_HPoint3Dd Surfacepoint;//nurbs曲面的采样点

	//*********************************新添加的private成员**********************************
	unsigned long long int nId ;//整体节点编号
	unsigned long long int mId ;//整体杆件编号
	//*********************************新添加的private成员**********************************
	


};


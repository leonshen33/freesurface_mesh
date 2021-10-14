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

struct equalDivPoints//��¼�߽�ͱ�ʶ��
{
	std::vector<double> points;
	double value; 
	
};

struct pointCoordinate//����
{
	int x;
	int y;
	pointCoordinate(int x0, int y0) :x(x0), y(y0){};
};
struct PointStruct
{
	PLib::Point2Dd point;//��¼u��v����
	std::vector<pointCoordinate> next;//��¼�õ���ָ��ĵ�����
	 long long int Id;//************************************************************��pId��ͬ�ľֲ����,��ʼΪ-1
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
	//NurbsSurface_qj(Vector_HPoint3Dd ccontrol_points,Vector_DOUBLE cknots,int cdegree,Vector_HPoint3Dd tcontrol_points,Vector_DOUBLE tknots,int tdegree);//����
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
	***�������ܣ���ȡһ�����������м�ֵ������
	***������minΪ���߶�Ӧ����С��v���꣬maxΪ���߶�Ӧ�����v���ꣻ
	***����ֵ�����������ϼ�ֵ������Ӧ��u����
	*****************getExtremums*************************/
	std::vector<double> getVExtremums(double min, double max, const double u);//uֵ�̶�����v����Ѱ�Ҽ�ֵ�㣬
	std::vector<double> getUExtremums(double min, double max, const double u);//vֵ�̶�����u����Ѱ�Ҽ�ֵ�㣬


	/****************getUInflectionPoints*************************
	***�������ܣ�����v�����ֵ����ȡU���������ϵĹյ㣻
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�������Ϲյ�����Ӧ��uֵ����
	***˼·���յ������z����Ķ��׵���ţ��趨�ϴ󲽳��������������ǰһ�����z����
	***      ���׵��ķ�����õ�ķ�����ţ��Ҹõ���ǰһ��ľ���̣ܶ�������Ϊ�õ�Ϊ�յ�
	*****************getInflectionPoints*************************/
	std::vector<double> getUInflectionPoints(double min, double max, const double v);

	/****************getVInflectionPoints*************************
	***�������ܣ�����u����ֵ����ȡV���������ϵĹյ㣻
	***������minΪ���߶�Ӧ����С��v���꣬maxΪ���߶�Ӧ�����v���ꣻ
	***����ֵ�������Ϲյ�����Ӧ��vֵ����
	***˼·���յ������z����Ķ��׵���ţ��趨�ϴ󲽳��������������ǰһ�����z����
	***      ���׵��ķ�����õ�ķ�����ţ��Ҹõ���ǰһ��ľ���̣ܶ�������Ϊ�õ�Ϊ�յ�
	*****************getInflectionPoints*************************/
	std::vector<double> getVInflectionPoints(double min, double max, const double u);

	/***********************getUParallelPoints********************
	***�������ܣ�����v����ֵ����ȡu���������ϵĵ���������ʼ�յ�����ƽ�еĵ��uֵ��
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�������Ϲյ�����Ӧ��uֵ����
	***˼·����Ѱ�ҹ������յ��Լ���ʼ�����յ���е��ƽ��l, ��ƽ��l�ķ�����n��Ȼ��������Ѱ����n
	***      ��ֱ�ĵ��Ӧ��uֵ
	***********************getUParallelPoints********************/
	std::vector<double> getUParallelPoints(double min, double max, const double v);

	/***********************getVParallelPoints********************
	***�������ܣ�����ʼ����u����ֵ����ȡv���������ϵĵ������յ�����ƽ�еĵ��vֵ��
	***������minΪ���߶�Ӧ����С��v���꣬maxΪ���߶�Ӧ�����v���ꣻ
	***����ֵ�������Ϲյ�����Ӧ��uֵ����
	***˼·����Ѱ�ҹ������յ��Լ���ʼ�����յ���е��ƽ��l, ��ƽ��l�ķ�����n��Ȼ��������Ѱ����n
	***      ��ֱ�ĵ��Ӧ��vֵ
	***********************getVParallelPoints********************/
	std::vector<double> getVParallelPoints(double min, double max, const double u);

	std::vector<double> getUAveragePerSection(double min, double max,const double v, int splitNum);//vֵ�̶�����u�������߰����ҳ��ָ�
	std::vector<double> getVAveragePersection(double min, double max,const double u, int splitNum);//uֵ�̶�����u�������߰����ҳ��ָ�



	
	std::vector<double> getUAverage(double min, double max, const double v, int splitNum);//�趨���ֶ��ٷݣ�vֵ�̶�����u�������߰����ҳ��ָ�
	std::vector<double> getVAverage(double min, double max, const double u, int splitNum);//�趨���ֶ��ٷݣ�uֵ�̶�����v�������߰����ҳ��ָ�

	std::vector<double> getUAverageLength(double min, double max, const double v, double lengthPerSection, int type = 0);//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�u�������ߵȷ�
	std::vector<double> getVAverageLength(double min, double max, const double u, double lengthPerSection, int type = 0);//�趨ÿ�ݵĳ��ȣ���һ�����Ƚ�v�������ߵȷ�
	
	std::vector<PLib::Point2Dd> GetSurfaceMaxExtrmums(double minU, double maxU, double minV, double maxV);//Ѱ������ļ���ֵ��
	std::vector<PLib::Point2Dd> GetSurfaceMinExtrmums(double minU, double maxU, double minV, double maxV);//Ѱ������ļ�Сֵ��

	
	//�����水�ҳ��ȷ�
	//type��ʾѰ����һ������ͣ� Ĭ��Ϊtype = 0����ʾ������������Ѱ�ң�type=1��������������Ѱ�ң�type=2�������ϵ����£�type =3�������ϵ�����		
	//std::vector<Member> surfaceAverageSpit(double minU, double realMaxU, double minV, double realMaxV, double splitLengthU, double splitLengthV, std::vector<Node>& pDisplayVec, int type=0);
	std::vector<Member> surfaceAverageSpit(rectangleRegion &rec, double splitLengthU, double splitLengthV, std::vector<PointStruct> &uvEdge);
	std::vector< Member> surfaceAverageTotalSpit(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec);
	//��nurbs�������ǵȷ֣������ȷֹ��ɵĵ� �͸˷ֱ�洢��nVec�ͷ���ֵ��
	
	//std::vector<Member> surfaceTriSplit(double minU, double maxU, double minV, double maxV, int uSplitNum, int vSplitNum, std::vector<Node> & nVec);
	
	double nextUPoint(double u, const double v, double uChordLength, int type = 0, double error = 1e-4);//����һ���u��v���꣬Ѱ��u�������ҳ�����chordLength��uֵ
	double nextVPoint(const double u, double v, double vChordLength, int type = 0, double error = 1e-4);//����һ���u��v���꣬Ѱ��v�������ҳ�����chordLength��vֵ
	//�����ı��ζԽǵ�������������߱߳���ȷ����һ��������,type��ʾѰ����һ������ͣ� Ĭ��Ϊtype=0����ʾ������������Ѱ�ң�
	//type=1��������������Ѱ�ң�type=2�������ϵ����£�type =3�������ϵ�����
	PLib::Point2Dd nextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type = 0, double error = 1e-4);
	void drawNurbsSurface();

	~NurbsSurface_qj(void);
	//*********************************����ӵ�public��Ա**********************************
	std::vector< Member> SurfaceAverageSpit_Rectangle(double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV, std::vector< Node> & nVec);
	std::vector<Member> MemberofRec_Rectangle(rectangleRegion rec, double splitLengthU, double splitLengthV, std::vector<Node> &nVec,std::vector<PointStruct> &uvEdge/*,std::vector<bool> &Gui_index*/);//�ֿ黭�˼�
	
	//�ڲ��߽紦��
	std::vector<Member> InsideEdge_Rectangle(double minU, double maxU, double minV, double maxV, CPoint2D maxPoint, double splitLengthU, double splitLengthV, std::vector< Node> & nVec, std::vector<PointStruct> &uvEdge);
	//���������uvֵ�������֮��ı߽�㣬���õ��u��vΪ1,u_edgeΪ�߽�uֵ,01��ʾ�ϲ����飬23��ʾ�²�����
	PointStruct PointBtween_Rectangle(PointStruct p1, PointStruct p2,double u_edge,double v_edge,int type);


	//�����ı��ζԽǵ�������������߱߳���ȷ����һ��������,type��ʾѰ����һ������ͣ� Ĭ��Ϊtype=0����ʾ������������Ѱ�ң�
	//type=1��������������Ѱ�ң�type=2�������ϵ����£�type =3�������ϵ�����
	PLib::Point2Dd NextCrossPoint(PLib::Point2Dd &point1, PLib::Point2Dd &point2, double uChrodLength, double vChrodLength, int type = 0, double error = 1e-4);

	struct IdRelationship//��¼�ı����ĵ�Id
	{
		unsigned long long int IdOne;//1ID
		unsigned long long int IdTwo;//2ID
		unsigned long long int IdThree;//3ID
		unsigned long long int IdFour;//4ID
	};

	std::vector<Node> UVForOut;//��¼���uvֵ��������

	/*map(NodeForOut���±꣬NodeForOut��ID��*/

	
	//std::vector<IdRelationship> IdRelShip;//���໥��ϵstd::vector
	void BetweenSearch(std::vector<PointStruct> uv_edge, double key1, double key2,int xy,int &up_index, int &down_index);//�жϿ����߽�ĵ���ڱ߽������䣬����up_index,down_index
	Member Add_Member(PointStruct point_in,PointStruct point_edge);//��Ӹ˼���point_in��һ�㣨�ڲ��㣩��point_edge�ڶ��㣨������߽�㣩
	Member Add_2Member(PointStruct point_in, PointStruct point_edge);//��Ӹ˼���point_in��һ�㣨�ڲ��㣩��point_edge�ڶ��㣨����߽�㣩

	//���طֿ���������ұ߽������������Ա߽����������ţ����߽����뵽������У�&�����޸�uvEdge��Id
	void Edge_udlr(std::vector<PointStruct> &uvEdge, rectangleRegion rec, std::vector<PointStruct> &uvedge_up, std::vector<PointStruct> &uvedge_down,std::vector<PointStruct> &uvedge_left,std::vector<PointStruct> &uvedge_right);
	
	//�߽�˼�����(ȫ���߽�������������ӱ߽�㡾���ڲ���λ�����߽���м�����߽�ܽ�ʱ����ĳһ�߽�ĵ㣬�����߽���ڲ���,uv����,��
	//ֻ������һ���߽������������������ڲ��㣬�����������ڲ�����������������������ܣ�
	std::vector<Member> Edge_Member(std::vector<PointStruct> &uvEdge, std::vector<PointStruct> uvedge, std::vector<PointStruct*> pointsedge,  int xy,double length_differ = 2);

	//����߽�����
	std::vector<Member> EdgeTotal_Member(std::vector<PointStruct> uvEdge, std::vector<PLib::Point2Dd> extrmumVec, int num, double minU, double maxU, double minV, double maxV, double splitLengthU, double splitLengthV);

	//����߽�Ľǵ�
	std::vector<Member> EdgeCorner_Member(std::vector<PointStruct*> &point_corner, double length_differ = 2);

	//*********************************����ӵ�public��Ա**********************************


private:
	Vector_HPoint3Dd TControlpoint;//nurbs����T�Ŀ��Ƶ�
	Vector_HPoint3Dd CControlpoint;//nurbs����C�Ŀ��Ƶ�
	
	Vector_HPoint3Dd TCurvepoint;//nurbs����T�Ĳ�����
	Vector_HPoint3Dd CCurvepoint;//nurbs����C�Ĳ�����

	Vector_DOUBLE TKnotpoint;//nurbs����T�Ľڵ�
	Vector_DOUBLE CKnotpoint;//nurbs����C�Ľڵ�
	
	int TDegree;//nurbs����T�Ķ�
	int CDegree;//nurbs����C�Ķ�
	
	PlNurbsCurved tnurbscurve ;//nurbs����T
	PlNurbsCurved cnurbscurve ;//nurbs����C

	PlNurbsSurfaced nurbssurface;


	Vector_HPoint3Dd Surfacepoint;//nurbs����Ĳ�����

	//*********************************����ӵ�private��Ա**********************************
	unsigned long long int nId ;//����ڵ���
	unsigned long long int mId ;//����˼����
	//*********************************����ӵ�private��Ա**********************************
	


};


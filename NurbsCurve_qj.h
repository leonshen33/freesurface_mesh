#pragma once

#include <iostream>
#include <fstream>
#include "nurbs.h"// for nurbs++
#include <gl/glew.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include "Node.h"

#include <Windows.h>
//using namespace std;

#pragma comment(linker, "/subsystem:console") 
#pragma comment( lib, "glew32.lib")
#pragma comment( lib, "glu32.lib")
#pragma comment( lib, "glut32.lib")
#pragma comment( lib, "opengl32.lib")
//using namespace PLib;

class NurbsCurve_qj
{
public:
	NurbsCurve_qj(void);
	NurbsCurve_qj(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree);//����
	NurbsCurve_qj(Vector_HPoint3Dd curve_points,int degree);//����
	NurbsCurve_qj(std::vector<Node*> control_points,std::vector<double > knots,int degree);
	NurbsCurve_qj(std::vector<Node*> curve_points,int degree);
	void reset(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree);
	Node* getCruvePointatU(double u);
	int getControlPoint(std::vector<Node*>& controlpoint);
	int getKnotPoint(double* knotpoint);
	bool isParallel(const Node &v0,const Node &v1);//����û��Node���ж����_IsParallel��

	//std::vector<Node> getAverageCurvePointsPerSection (int splitNum, double min, double max);//���ݷָ����splitNum��ȡÿ�������ϵȳ��ҵķָ��,�ָ��洢��nodes��
	

	/****************getAverageCurvePoints*************************
	***�������ܣ���ȡһ�����߰��ҳ��ȷֵ����е������
	***�����������м�ֵ�㡢�յ��һЩ����ĵ���Ҫ�������������ߵĵȷֲ����Ǿ��������ϵĵȷ֣����Ǳ�֤�����֮��������ǵȷֵģ�
	***������splitNumΪ�ȷֶ��ٷݣ�minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�����������ϵȷֵ������ֵ
	*****************getAverageCurvePoints:*************************/
	std::vector<Node*> getAverageCurvePoints(int splitNum, double min, double max);//���ݼ�ֵ��ָ�����䣬����getAverageCurvePointsPerSection���зָ�
	
	/****************getExtremums*************************
	***�������ܣ���ȡһ�����������м�ֵ������
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�����������ϼ�ֵ������Ӧ��u����
	*****************getExtremums*************************/
	std::vector<double> getExtremums(double min, double max);

	
	/****************getKnees*************************
	***�������ܣ���ȡһ��������һЩ��Ҫ����������㣬��Щ��������루������յ㣩���ߵķ���һ�£�
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�����������������Ӧ��uֵ����
	*****************getKnees*************************/
	std::vector<double> getKnees(double min,  double max);//��ȡ���������߷�����min��max����ƽ�еĵ�

	/****************getInflectionPoints*************************
	***�������ܣ���ȡһ�������ϵĹյ㣻
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�������Ϲյ�����Ӧ��uֵ����
	*****************getInflectionPoints*************************/
	std::vector<double> getInflectionPoints(double min, double max);

	/****************getExtremumsAndKnees*************************
	***�������ܣ���ȡһ�����������м�ֵ��͹յ��Լ�һЩ��Ҫ�����������
	***������minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�����������ϼ�ֵ�㡢�յ��Լ����������Ӧ��u����
	*****************getExtremumsAndKnees*************************/
	std::vector<double> getExtremumsAndKnees(double min, double max);


	/****************getAvergreArcPointsPerSection*************************
	***�������ܣ���ȡһ�������Ͻ����߰��ҳ��ȷֵĵȷֵ�
	***������splitNumΪ�ö�Ҫ�ȷֵĴ�����minΪ���߶�Ӧ����С��u���꣬maxΪ���߶�Ӧ�����u���ꣻ
	***����ֵ�����ö����߰��ҳ��ȷֵĵȷֵ�
	*****************getAvergreArcPointsPerSection*************************/
	std::vector<Node*> getAvergreArcPointsPerSection(int splitNum, double min, double max);//��ȡÿ�������ϵȷָ����ߵĵ�����

    void drawNurbsCurve();

	//void display();
	
	~NurbsCurve_qj(void);

private:
	Vector_HPoint3Dd Controlpoint;
	Vector_HPoint3Dd Curvepoint;
	Vector_DOUBLE Knotpoint;
	int Degree;
	PlNurbsCurved nurbscurve ;


};


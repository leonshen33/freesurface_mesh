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
	NurbsCurve_qj(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree);//正算
	NurbsCurve_qj(Vector_HPoint3Dd curve_points,int degree);//反算
	NurbsCurve_qj(std::vector<Node*> control_points,std::vector<double > knots,int degree);
	NurbsCurve_qj(std::vector<Node*> curve_points,int degree);
	void reset(Vector_HPoint3Dd control_points,Vector_DOUBLE knots,int degree);
	Node* getCruvePointatU(double u);
	int getControlPoint(std::vector<Node*>& controlpoint);
	int getKnotPoint(double* knotpoint);
	bool isParallel(const Node &v0,const Node &v1);//精度没有Node类中定义的_IsParallel高

	//std::vector<Node> getAverageCurvePointsPerSection (int splitNum, double min, double max);//根据分割次数splitNum获取每段曲线上等长弦的分割点,分割点存储在nodes中
	

	/****************getAverageCurvePoints*************************
	***函数功能：获取一条曲线按弦长等分的所有点的坐标
	***由于曲线上有极值点、拐点和一些特殊的点需要保留，所以曲线的等分并不是绝对意义上的等分，而是保证特殊点之间的曲线是等分的；
	***参数：splitNum为等分多少份，min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：所有曲线上等分点的坐标值
	*****************getAverageCurvePoints:*************************/
	std::vector<Node*> getAverageCurvePoints(int splitNum, double min, double max);//根据极值点分割的区间，调用getAverageCurvePointsPerSection进行分割
	
	/****************getExtremums*************************
	***函数功能：获取一条曲线上所有极值点坐标
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：所有曲线上极值点所对应的u坐标
	*****************getExtremums*************************/
	std::vector<double> getExtremums(double min, double max);

	
	/****************getKnees*************************
	***函数功能：获取一条曲线上一些需要保留的特殊点，这些点的切线与（起点与终点）连线的方向一致；
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：曲线上特殊点所对应的u值坐标
	*****************getKnees*************************/
	std::vector<double> getKnees(double min,  double max);//获取曲线上切线方向与min和max连线平行的点

	/****************getInflectionPoints*************************
	***函数功能：获取一条曲线上的拐点；
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：曲线上拐点所对应的u值坐标
	*****************getInflectionPoints*************************/
	std::vector<double> getInflectionPoints(double min, double max);

	/****************getExtremumsAndKnees*************************
	***函数功能：获取一条曲线上所有极值点和拐点以及一些需要保留的特殊点
	***参数：min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：所有曲线上极值点、拐点以及特殊点所对应的u坐标
	*****************getExtremumsAndKnees*************************/
	std::vector<double> getExtremumsAndKnees(double min, double max);


	/****************getAvergreArcPointsPerSection*************************
	***函数功能：获取一段曲线上将曲线按弦长等分的等分点
	***参数：splitNum为该段要等分的次数，min为曲线对应的最小的u坐标，max为曲线对应的最大u坐标；
	***返回值：将该段曲线按弦长等分的等分点
	*****************getAvergreArcPointsPerSection*************************/
	std::vector<Node*> getAvergreArcPointsPerSection(int splitNum, double min, double max);//获取每段曲线上等分该曲线的点坐标

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


#pragma once
#ifndef _MATRIX_H
#define _MATRIX_H
#include<iostream>
#include "Node.h"
#include "CadBase.h"
//using namespace std;

class  CMatrix3D;

class CMatrix
{
    public:

        CMatrix(int row=1, int col=1);             // 带默认参数的构造函数
        CMatrix(int row, int col, double mtx[]);    // 用数组创建一个矩阵
        CMatrix(const CMatrix &obj);                 // copy构造函数
		//CMatrix(Matrix3d &thematrix3d);
		CMatrix(CMatrix3D &thematrix3d);

       ~CMatrix() { delete[] this->mtx; }   


		int getRow()const ;    // 访问矩阵行数
        int getCol()const;     // 访问矩阵列数
        int getN()const ;      // 访问矩阵元素个数
        double* getMtx()const; // 获取该矩阵的数组
        double get(const int i, const int j)const; // 用下标访问矩阵元素
        void set(const int i, const int j, const double replacedata); // 用下标修改矩阵元素值
		void output()const;    // 利用输出函数 输出矩阵
		friend std::ostream &operator <<(std::ostream & out,const CMatrix &copymtx);//重载运算符<<，实现矩阵按照行列的格式输出
	
     // 重载了一些常用操作符，包括 +，-，x，=，负号，正号，
		CMatrix &operator= (const CMatrix &obj);// A = B
		CMatrix  operator+ ()const ;// +A
        CMatrix  operator- ()const;// -A
        friend  CMatrix  operator+ (const CMatrix &A, const CMatrix &B);// A + B
        friend  CMatrix  operator- (const CMatrix &A, const CMatrix &B);// A - B
        friend  CMatrix  operator* (const CMatrix &A, const CMatrix &B);// A * B 两矩阵相乘
        friend  CMatrix  operator* (const double &a, const CMatrix &B);// a * B 实数与矩阵相乘 
        CMatrix  transpose(const CMatrix &A); // A 的转置 
		double  det(CMatrix A); // 求A 的行列式值，采用列主元消去法,将矩阵化为三角阵,此处为了防止修改原矩阵，采用传值调用    
		CMatrix  inv(CMatrix A);// A 的逆矩阵，采用高斯-若当列主元消去法
      
private:	

		int row;    // 矩阵的行数
		int col;    // 矩阵的列数
		int n;      // 矩阵元素个数
		double*  mtx; // 动态分配用来存放数组的空间
};
#endif        


 

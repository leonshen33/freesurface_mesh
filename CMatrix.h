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

        CMatrix(int row=1, int col=1);             // ��Ĭ�ϲ����Ĺ��캯��
        CMatrix(int row, int col, double mtx[]);    // �����鴴��һ������
        CMatrix(const CMatrix &obj);                 // copy���캯��
		//CMatrix(Matrix3d &thematrix3d);
		CMatrix(CMatrix3D &thematrix3d);

       ~CMatrix() { delete[] this->mtx; }   


		int getRow()const ;    // ���ʾ�������
        int getCol()const;     // ���ʾ�������
        int getN()const ;      // ���ʾ���Ԫ�ظ���
        double* getMtx()const; // ��ȡ�þ��������
        double get(const int i, const int j)const; // ���±���ʾ���Ԫ��
        void set(const int i, const int j, const double replacedata); // ���±��޸ľ���Ԫ��ֵ
		void output()const;    // ����������� �������
		friend std::ostream &operator <<(std::ostream & out,const CMatrix &copymtx);//���������<<��ʵ�־��������еĸ�ʽ���
	
     // ������һЩ���ò����������� +��-��x��=�����ţ����ţ�
		CMatrix &operator= (const CMatrix &obj);// A = B
		CMatrix  operator+ ()const ;// +A
        CMatrix  operator- ()const;// -A
        friend  CMatrix  operator+ (const CMatrix &A, const CMatrix &B);// A + B
        friend  CMatrix  operator- (const CMatrix &A, const CMatrix &B);// A - B
        friend  CMatrix  operator* (const CMatrix &A, const CMatrix &B);// A * B ���������
        friend  CMatrix  operator* (const double &a, const CMatrix &B);// a * B ʵ���������� 
        CMatrix  transpose(const CMatrix &A); // A ��ת�� 
		double  det(CMatrix A); // ��A ������ʽֵ����������Ԫ��ȥ��,������Ϊ������,�˴�Ϊ�˷�ֹ�޸�ԭ���󣬲��ô�ֵ����    
		CMatrix  inv(CMatrix A);// A ������󣬲��ø�˹-��������Ԫ��ȥ��
      
private:	

		int row;    // ���������
		int col;    // ���������
		int n;      // ����Ԫ�ظ���
		double*  mtx; // ��̬���������������Ŀռ�
};
#endif        


 

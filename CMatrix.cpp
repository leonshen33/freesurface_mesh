//#include "StdAfx.h"
#include "CMatrix.h"
#include<iostream>
#include<math.h>                
#include<stdlib.h>                
#include<iomanip>    
#include "Node.h"

//using namespace std;

// ��Ĭ�ϲ���ֵ�Ĺ��캯��,����һ��row��col�е������
CMatrix::CMatrix(int row, int col)
{                
    this->row = row;
    this->col = col;                                                         
    this->n = row * col;
    this->mtx = new double[n];
                                               
    for(int i=0; i<n; i++)
        this->mtx[i] = 0.0;
}    

// ��һ�������ʼ������                                                        
CMatrix::CMatrix(int row, int col, double mtx[])                                  
{
    this->row = row;
    this->col = col;                                                         
    this->n   = row * col;    
    this->mtx = new double[n];                                               
    
	for(int i=0; i<n; i++)                                                     
		this->mtx[i] = mtx[i];                                                
}    

// �������캯������Ϊ��Ա�������ж�̬�ռ䣬��ֹ���ݲ��� �Ȳ�����������
CMatrix::CMatrix(const CMatrix &obj)
{            
    this->row = obj.getRow();                                         
    this->col = obj.getCol();                                          
    this->n  = obj.getN();                                            
    this->mtx = new double[n];
                                            
    for(int i=0; i<n; i++)
       this->mtx[i] = obj.getMtx()[i];
}  

CMatrix::CMatrix(CMatrix3D &thematrix3d)
{
	this->row = 4;
	this->col = 4;                                                         
	this->n   = 16;    
	this->mtx = new double[16];                                               

	for(int i=0; i<4; i++)  
		for (int j = 0;j<4; j++)
		{
			this->mtx[i*4+j] = thematrix3d.A[i][j];     
		}
		
		

}


int CMatrix::getRow()const // ���ʾ�������
{
	return this->row;
}
int CMatrix::getCol()const// ���ʾ������� 
{ 
	return this->col;
}
int CMatrix::getN()const // ���ʾ���Ԫ�ظ���  
{ 
	return this->n;  
}
double*  CMatrix::getMtx()const // ��ȡ�þ��������
{ 
	return this->mtx;
}
                          
double CMatrix::get(const int i, const int j)const// ��ȡ����Ԫ��,ע����������±��(0,0)��ʼ
{                
    return this->mtx[i*this->col + j];                           
}
 
void CMatrix::set(const int i, const int j, const double replacedata)// �޸ľ���Ԫ�� 
{                  
    this->mtx[i*this->col + j] = replacedata;                                 
} 

// ��������� �����������Ԫ��
void CMatrix::output()const
{                
    for(int i=0; i<this->row; i++)
	{                                         
       for(int j=0; j<this->col; j++)
          if(fabs(this->get(i,j)) <= 1.0e-10)
              std::cout<< std::setiosflags(std::ios::left) << std::setw(4) << 0.0 << ' ';
          else
              std::cout<< std::setiosflags(std::ios::left) << std::setw(4) << this->get(i,j) << ' ';
              std::cout <<std::endl;      
    }          
}

//���������<<��ʵ�־��������еĸ�ʽ���
std::ostream &operator <<(std::ostream & out,const CMatrix &copymtx)
{
	out<<std::endl;
	for(int i=0;i<copymtx.row;i++)
	{
		for(int j=0;j<copymtx.col;j++)
			 if(fabs(copymtx.get(i,j)) <= 1.0e-10)
				 out<< std::setiosflags(std::ios::left) << std::setw(4) << 0.0 << ' ';
			 else
				 out<<std::setiosflags(std::ios::left) <<std::setw(4)<<*(copymtx.mtx +i*(copymtx.col)+j);
				 out<<std::endl;
	}
	return out;
}  

// ���ظ�ֵ�����������ڳ�Ա�����к��ж�̬����
CMatrix &CMatrix::operator= (const CMatrix &obj)
{                
    if(this == &obj)    // ��һ�����󸳸����Լ�ʱ�������ؼ���        
       return *this;
    delete[] this->mtx; // ����ɾ��Ŀ�Ĳ������Ķ�̬�ռ�
    this->row = obj.getRow();
    this->col = obj.getCol();
    this->n   = obj.getN();
    this->mtx = new double[n]; // ���·���ռ䱣��obj������

    for(int i=0; i<n; i++)
       this->mtx[i] = obj.getMtx()[i];                              
    return *this;
}   
 
// �ӺŲ�����������ֵΪ�þ���
CMatrix  CMatrix::operator+ ()const
{ 
	return *this;
}
// ���Ų�����������ֵΪ�þ���ĸ�����ԭ���󲻱�
CMatrix CMatrix::operator- ()const
{                
    CMatrix newmtx(this->row, this->col);// Ϊ�˲��ı�ԭ���ľ��󣬴˴����¹���һ������

    for(int i=0; i<newmtx.n; i++)
       newmtx.mtx[i] = -(this->mtx[i]);
    return newmtx;       
}    

// ������ͣ���Ӧλ��Ԫ�����
CMatrix operator+ (const CMatrix &A, const CMatrix &B)
{                
    CMatrix plusAB(A.row, A.col);

    if(A.row!=B.row || A.col!=B.col)
	{
       std::cout << "Can't do A+B\n"; // �������A��B��������һ���򲻿����
       exit(0);        
    }

    for(int i=0; i<plusAB.n; i++)
       plusAB.mtx[i] = A.mtx[i] + B.mtx[i];
    
        return plusAB;       
}     

// ����������ü���һ����������ʵ��
CMatrix operator- (const CMatrix &A, const CMatrix &B)
{                
    return (A + (-B));
}  

// ����˷�
CMatrix operator* (const CMatrix &A, const CMatrix &B)
{                
    if(A.col != B.row){    // A�����������B������һ��
       std::cout << "Can't multiply\n";
       exit(0);          
    }                                                                          

    CMatrix AB(A.row, B.col); // AB���ڱ���˻�

    for(int i=0; i<AB.row; i++)
       for(int j=0; j<AB.col; j++)
          for(int k=0; k<A.col; k++)
              AB.set(i, j, AB.get(i,j) + A.get(i,k)*B.get(k,j));
    return AB;       
}

// ������ʵ�����
CMatrix operator* (const double &a, const CMatrix &B)
{
    CMatrix aB(B);
    for(int i=0; i<aB.row; i++)
       for(int j=0; j<aB.col; j++)
          aB.set(i,j, a*B.get(i,j));    
    return aB;
}

// �����ת�� ��(i,j)��(j,i)����  
// �˺�������һ�������ת�þ��󣬲����ı�ԭ���ľ���  
CMatrix CMatrix::transpose(const CMatrix &A)
{  
	CMatrix AT(A.col, A.row);  
	for(int i=0; i<AT.row; i++)  
	    for(int j=0; j<AT.col; j++)  
	        AT.set(i, j, A.get(j,i));  
	return AT;         
}

// ��������ʽֵ����������Ԫ��ȥ��
double CMatrix::det(CMatrix A)
{
    if(A.row != A.col) 
	{    // �������Ϊn*n�Ĳſɽ�������ʽ��ֵ       
        std::cout << "error" << std::endl;
        return 0.0;     //�����������������ȷ���0.0    
	}

    double detValue = 1.0;  // ���ڱ�������ʽֵ
    for(int i=0; i<A.getRow()-1; i++)
	{ // ��Ҫn-1���л������

         //------------------ ѡ��Ԫ ---------------------------------    
       double max = fabs(A.get(i,i));    // ��Ԫ��ʼĬ��Ϊ���·������׸�Ԫ��     
        int    ind = i;  
          // ��Ԫ�к�Ĭ��Ϊ���·���������      
    
       for(int j=i+1; j<A.getRow(); j++)
	   {  // ѡ������Ԫ
        if(fabs(A.get(j,i)) > max){     // ��������ֵ�����Ԫ��
           max = fabs(A.get(j,i));      // ������Ԫֵ
           ind = j;                       // ������Ԫ�к�
        }                  
      }//loop j    
    
       //------------------- �ƶ���Ԫ�� -----------------------------
       if(max <= 1.0e-10) return 0.0;    // ���·���������Ϊ�㣬��Ȼ����ʽֵΪ��      
        if(ind != i)
		{// ��Ԫ�з����·���������                                  
			 for(int k=i; k<A.getRow(); k++)
			 {     // ����Ԫ�������·��������л���          
				double temp = A.get(i,k);
				A.set(i,k,A.get(ind,k));
				A.set(ind,k,temp);
			}
          detValue = -detValue;             // ��������ʽ���У�����ʽֵ����    
        }                            
    
       //------------------- ��Ԫ ----------------------------------                                  
        for(int j=i+1; j<A.getRow(); j++)
		{     // ������  
           double temp = A.get(j,i)/A.get(i,i);        
           for(int k=i; k<A.getRow(); k++) // ��������ÿ��Ԫ�أ�������0
           A.set(j, k, A.get(j,k)-A.get(i,k)*temp);     
		}                    
        detValue *= A.get(i,i); // ÿ����Ԫ�������һ���Խ�����Ԫ�أ������۳�
       }// loop i     
      
      // ע��������һ��Ԫ������Ԫ�Ĺ�����û�б��۳˵�
      return detValue * A.get(A.getRow()-1, A.getRow()-1);
}


// A������� ��˹-������ȥ��������ѡ��Ԫ
CMatrix CMatrix::inv(CMatrix A)
{                
	if(A.row != A.col)
	{ // ֻ������������󣬼���������ͬ
		std::cout << "CMatrix should be N x N\n";
		exit(0);       
	}                                                                                                
	// ����һ����A������ͬ�ĵ�λ��B
	CMatrix B(A.row,A.col);

	for(int r=0; r<A.row; r++)
		for(int c=0; c<A.col; c++)
			if(r == c) B.set(r,c,1.0); 

	// �Ծ���A����A.row����Ԫ���㣬ÿ�α�֤��K��ֻ�жԽ����Ϸ���
	// ͬʱ��ͬ���Ĳ���ʩ�����B�����A��Ϊ��λ��BΪ��������
	for(int k=0; k<A.row; k++){
		//------------------ ѡ��Ԫ --------------------------------------    
		double max = fabs(A.get(k,k));    // ��Ԫ��ʼĬ��Ϊ���·������׸�Ԫ��     
		int    ind = k;            // ��Ԫ�к�Ĭ��Ϊ���·���������
		// �����ind��Ϊ����Ԫ��
		for(int n=k+1; n<A.getRow(); n++){
			if(fabs(A.get(n,k)) > max){// ��������ֵ�����Ԫ��             
				max = fabs(A.get(n,k));    // ������Ԫֵ
				ind = n;// ������Ԫ�к�
			}
		}                   
		//------------------- �ƶ���Ԫ�� --------------------------------
		if(ind != k)
		{// ��Ԫ�в������·���������
			for(int m=k; m<A.row; m++)
			{// ����Ԫ�������·��������л���         
				double tempa = A.get(k,m);
				A.set(k, m, A.get(ind,m));
				A.set(ind, m, tempa);            
			}
			for(int m=0; m<B.row; m++)
			{
				double tempb = B.get(k,m); // �Ծ���Bʩ����ͬ����                  
				B.set(k, m, B.get(ind,m)); // B��A������ͬ������һ��ѭ����               
				B.set(ind, m, tempb);    
			}    
		}        
		//--------------------- ��Ԫ -----------------------------------
		// ��k����Ԫ�������Ե�k����Ϊ��Ԫ�У��������¸��еĵ�k��Ԫ�ػ�Ϊ��
		// ͬʱ��ͬ���Ĳ�����Bʩ��ͬ���Ĳ�������ʱ���Խ�B����A�����һ����          
		for(int i=0; i<A.col; i++)
		{                                 
			if(i != k)
			{
				double Mik = -A.get(i,k)/A.get(k,k);                         
				for(int j=k+1; j<A.row; j++)
					A.set(i, j, A.get(i,j) + Mik*A.get(k,j));                     
				for(int j=0; j<B.row; j++)
					B.set(i, j, B.get(i,j) + Mik*B.get(k,j));                   
			}//end if 
		}//loop i  
		double Mkk = 1.0/A.get(k,k);

		for(int j=0; j<A.row; j++)
			A.set(k, j, A.get(k,j) * Mkk);
		for(int j=0; j<B.row; j++)
			B.set(k, j, B.get(k,j) * Mkk);
	}//loop k      
	return B;
}
                                                               
                


                                                                                  
//#include "StdAfx.h"
#include "CMatrix.h"
#include<iostream>
#include<math.h>                
#include<stdlib.h>                
#include<iomanip>    
#include "Node.h"

//using namespace std;

// 带默认参数值的构造函数,构造一个row行col列的零矩阵
CMatrix::CMatrix(int row, int col)
{                
    this->row = row;
    this->col = col;                                                         
    this->n = row * col;
    this->mtx = new double[n];
                                               
    for(int i=0; i<n; i++)
        this->mtx[i] = 0.0;
}    

// 用一个数组初始化矩阵                                                        
CMatrix::CMatrix(int row, int col, double mtx[])                                  
{
    this->row = row;
    this->col = col;                                                         
    this->n   = row * col;    
    this->mtx = new double[n];                                               
    
	for(int i=0; i<n; i++)                                                     
		this->mtx[i] = mtx[i];                                                
}    

// 拷贝构造函数，因为成员变量含有动态空间，防止传递参数 等操作发生错误
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


int CMatrix::getRow()const // 访问矩阵行数
{
	return this->row;
}
int CMatrix::getCol()const// 访问矩阵列数 
{ 
	return this->col;
}
int CMatrix::getN()const // 访问矩阵元素个数  
{ 
	return this->n;  
}
double*  CMatrix::getMtx()const // 获取该矩阵的数组
{ 
	return this->mtx;
}
                          
double CMatrix::get(const int i, const int j)const// 获取矩阵元素,注意这里矩阵下标从(0,0)开始
{                
    return this->mtx[i*this->col + j];                           
}
 
void CMatrix::set(const int i, const int j, const double replacedata)// 修改矩阵元素 
{                  
    this->mtx[i*this->col + j] = replacedata;                                 
} 

// 用输出函数 输出矩阵所有元素
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

//重载运算符<<，实现矩阵按照行列的格式输出
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

// 重载赋值操作符，由于成员变量中含有动态分配
CMatrix &CMatrix::operator= (const CMatrix &obj)
{                
    if(this == &obj)    // 将一个矩阵赋给它自己时简单做返回即可        
       return *this;
    delete[] this->mtx; // 首先删除目的操作数的动态空间
    this->row = obj.getRow();
    this->col = obj.getCol();
    this->n   = obj.getN();
    this->mtx = new double[n]; // 重新分配空间保存obj的数组

    for(int i=0; i<n; i++)
       this->mtx[i] = obj.getMtx()[i];                              
    return *this;
}   
 
// 加号操作符，返回值为该矩阵
CMatrix  CMatrix::operator+ ()const
{ 
	return *this;
}
// 负号操作符，返回值为该矩阵的负矩阵，原矩阵不变
CMatrix CMatrix::operator- ()const
{                
    CMatrix newmtx(this->row, this->col);// 为了不改变原来的矩阵，此处从新构造一个矩阵

    for(int i=0; i<newmtx.n; i++)
       newmtx.mtx[i] = -(this->mtx[i]);
    return newmtx;       
}    

// 矩阵求和，对应位置元素相加
CMatrix operator+ (const CMatrix &A, const CMatrix &B)
{                
    CMatrix plusAB(A.row, A.col);

    if(A.row!=B.row || A.col!=B.col)
	{
       std::cout << "Can't do A+B\n"; // 如果矩阵A和B行列数不一致则不可相加
       exit(0);        
    }

    for(int i=0; i<plusAB.n; i++)
       plusAB.mtx[i] = A.mtx[i] + B.mtx[i];
    
        return plusAB;       
}     

// 矩阵减法，用加上一个负矩阵来实现
CMatrix operator- (const CMatrix &A, const CMatrix &B)
{                
    return (A + (-B));
}  

// 矩阵乘法
CMatrix operator* (const CMatrix &A, const CMatrix &B)
{                
    if(A.col != B.row){    // A的列数必须和B的行数一致
       std::cout << "Can't multiply\n";
       exit(0);          
    }                                                                          

    CMatrix AB(A.row, B.col); // AB用于保存乘积

    for(int i=0; i<AB.row; i++)
       for(int j=0; j<AB.col; j++)
          for(int k=0; k<A.col; k++)
              AB.set(i, j, AB.get(i,j) + A.get(i,k)*B.get(k,j));
    return AB;       
}

// 矩阵与实数相乘
CMatrix operator* (const double &a, const CMatrix &B)
{
    CMatrix aB(B);
    for(int i=0; i<aB.row; i++)
       for(int j=0; j<aB.col; j++)
          aB.set(i,j, a*B.get(i,j));    
    return aB;
}

// 矩阵的转置 将(i,j)与(j,i)互换  
// 此函数返回一个矩阵的转置矩阵，并不改变原来的矩阵  
CMatrix CMatrix::transpose(const CMatrix &A)
{  
	CMatrix AT(A.col, A.row);  
	for(int i=0; i<AT.row; i++)  
	    for(int j=0; j<AT.col; j++)  
	        AT.set(i, j, A.get(j,i));  
	return AT;         
}

// 矩阵行列式值，采用列主元消去法
double CMatrix::det(CMatrix A)
{
    if(A.row != A.col) 
	{    // 矩阵必须为n*n的才可进行行列式求值       
        std::cout << "error" << std::endl;
        return 0.0;     //如果不满足行列数相等返回0.0    
	}

    double detValue = 1.0;  // 用于保存行列式值
    for(int i=0; i<A.getRow()-1; i++)
	{ // 需要n-1步列化零操作

         //------------------ 选主元 ---------------------------------    
       double max = fabs(A.get(i,i));    // 主元初始默认为右下方矩阵首个元素     
        int    ind = i;  
          // 主元行号默认为右下方矩阵首行      
    
       for(int j=i+1; j<A.getRow(); j++)
	   {  // 选择列主元
        if(fabs(A.get(j,i)) > max){     // 遇到绝对值更大的元素
           max = fabs(A.get(j,i));      // 更新主元值
           ind = j;                       // 更新主元行号
        }                  
      }//loop j    
    
       //------------------- 移动主元行 -----------------------------
       if(max <= 1.0e-10) return 0.0;    // 右下方矩阵首行为零，显然行列式值为零      
        if(ind != i)
		{// 主元行非右下方矩阵首行                                  
			 for(int k=i; k<A.getRow(); k++)
			 {     // 将主元行与右下方矩阵首行互换          
				double temp = A.get(i,k);
				A.set(i,k,A.get(ind,k));
				A.set(ind,k,temp);
			}
          detValue = -detValue;             // 互换行列式两行，行列式值反号    
        }                            
    
       //------------------- 消元 ----------------------------------                                  
        for(int j=i+1; j<A.getRow(); j++)
		{     // 遍历行  
           double temp = A.get(j,i)/A.get(i,i);        
           for(int k=i; k<A.getRow(); k++) // 遍历行中每个元素，行首置0
           A.set(j, k, A.get(j,k)-A.get(i,k)*temp);     
		}                    
        detValue *= A.get(i,i); // 每步消元都会产生一个对角线上元素，将其累乘
       }// loop i     
      
      // 注意矩阵最后一个元素在消元的过程中没有被累乘到
      return detValue * A.get(A.getRow()-1, A.getRow()-1);
}


// A的逆矩阵 高斯-若当消去法，按列选主元
CMatrix CMatrix::inv(CMatrix A)
{                
	if(A.row != A.col)
	{ // 只可求狭义逆矩阵，即行列数相同
		std::cout << "CMatrix should be N x N\n";
		exit(0);       
	}                                                                                                
	// 构造一个与A行列相同的单位阵B
	CMatrix B(A.row,A.col);

	for(int r=0; r<A.row; r++)
		for(int c=0; c<A.col; c++)
			if(r == c) B.set(r,c,1.0); 

	// 对矩阵A进行A.row次消元运算，每次保证第K列只有对角线上非零
	// 同时以同样的操作施与矩阵B，结果A变为单位阵B为所求逆阵
	for(int k=0; k<A.row; k++){
		//------------------ 选主元 --------------------------------------    
		double max = fabs(A.get(k,k));    // 主元初始默认为右下方矩阵首个元素     
		int    ind = k;            // 主元行号默认为右下方矩阵首行
		// 结果第ind行为列主元行
		for(int n=k+1; n<A.getRow(); n++){
			if(fabs(A.get(n,k)) > max){// 遇到绝对值更大的元素             
				max = fabs(A.get(n,k));    // 更新主元值
				ind = n;// 更新主元行号
			}
		}                   
		//------------------- 移动主元行 --------------------------------
		if(ind != k)
		{// 主元行不是右下方矩阵首行
			for(int m=k; m<A.row; m++)
			{// 将主元行与右下方矩阵首行互换         
				double tempa = A.get(k,m);
				A.set(k, m, A.get(ind,m));
				A.set(ind, m, tempa);            
			}
			for(int m=0; m<B.row; m++)
			{
				double tempb = B.get(k,m); // 对矩阵B施以相同操作                  
				B.set(k, m, B.get(ind,m)); // B与A阶数相同，可在一个循环中               
				B.set(ind, m, tempb);    
			}    
		}        
		//--------------------- 消元 -----------------------------------
		// 第k次消元操作，以第k行作为主元行，将其上下各行的第k列元素化为零
		// 同时以同样的参数对B施以同样的操作，此时可以将B看作A矩阵的一部分          
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
                                                               
                


                                                                                  
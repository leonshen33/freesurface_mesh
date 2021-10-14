#include "Node.h"
//#include "stdafx.h"
#include "math.h"
#include <iostream>
//Node::Node(void)
//{
//	this->No = 0;
//	this->X = 0;
//	this->Y = 0;
//	this->Z = 0;
//	this->Color = 0;
//}

//operators

Node Node::operator*(const CMatrix3D &matrix ) const
{
	double rx,ry,rz,sc;
	rx = X* matrix.A[0][0] + Y * matrix.A[1][0] + Z * matrix.A[2][0] + matrix.A[3][0];
	ry = X * matrix.A[0][1] + Y * matrix.A[1][1] + Z * matrix.A[2][1] + matrix.A[3][1];
	rz = X * matrix.A[0][2] + Y * matrix.A[1][2] + Z * matrix.A[2][2] + matrix.A[3][2];
	sc = X * matrix.A[0][3] + Y* matrix.A[1][3] + Z * matrix.A[2][3] + matrix.A[3][3];
	rx /= sc;
	ry /= sc;
	rz /= sc;
	return Node(rx,ry,rz);
}

 void  Node::operator*=(const CMatrix3D& matrix)
 {
	 (*this) = (*this)*matrix;
 }

 // offsetting with vector
 Node Node::operator+(Node v) const
 {
	 return Node(X+v.X,Y+v.Y,Z+v.Z);
 }

 void Node::operator+=(Node v)
 {
	 X+=v.X;
	 Y+=v.Y;
	 Z+=v.Z;
 }

 Node Node::operator-(Node v) const
 {
	 return Node(X-v.X,Y-v.Y,Z-v.Z);
 }

 void Node::operator-=(Node v)
 {
	 X-=v.X;
	 Y-=v.Y;
	 Z-=v.Z;
 }

 bool Node::operator==(Node pos) const
 {
	 Node  vect(X-pos.X,Y-pos.Y,Z-pos.Z);

	 if( IS_ZERO( vect.GetLength( ) ) ) return true;
	 else return false;

 }
 bool Node::operator!=(Node pos) const
 {
	 Node vect(X-pos.X,Y-pos.Y,Z-pos.Z);

	 if( IS_ZERO( vect.GetLength( ) ) ) return false;
	 else return true;
 }

 double Node::GetLength() const
 {
	 return sqrt(X*X+Y*Y+Z*Z);
 }
Node Node::operator*(double d) const
 {
	 return Node(X*d,Y*d,Z*d);
 }
void Node::operator*=(double d)
{
	X *= d;
	Y *= d;
	Z *= d;
}
Node Node::operator/(double d) const
{
	return Node(X/d,Y/d,Z/d);
}

void Node::operator/=(double d)
{
	X /= d;
	Y /= d;
	Z/= d;
}
// cross product
Node Node::operator*(Node v) const
{
	return Node(Y*v.Z-v.Y*Z,v.X*Z-X*v.Z,X*v.Y-v.X*Y);
}

// dot product
double Node::operator|(Node v) const
{
	return X*v.X+Y*v.Y+Z*v.Z;
}
double Node::GetLengthXY() const
{
	return sqrt(X*X+Y*Y);
}

double Node::GetLengthYZ() const
{
	return sqrt(Y*Y+Z*Z);
}

double Node::GetLengthZX() const
{
	return sqrt(X*X+Z*Z);
}
Node Node::GetNormal() const
{
	double len = GetLength();
	return Node(X/len,Y/len,Z/len);
}

void Node::Normalize()
{
	double len = GetLength();
	X /= len;
	Y /= len;
	Z /= len;
}

bool Node::IsZeroLength() const
{
	return IS_ZERO(GetLength());
}





//Node  Node::operator =(Node temp)//将一个点赋给另外一个点
//{
//	Node ret;
//	ret.setX(temp.getX());
//	ret.setY(temp.getY());
//	ret.setZ(temp.getZ());
//	return ret;
//}

/************************************************************************************/

//	//construction&&destruction
//	CMatrix3D::Matrix3d()
//{
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++)
//		{
//			A[i][j] = (i==j)?1.0:0.0;
//		}
//}
//
//CMatrix3D::CMatrix3D(const CMatrix3D& matrix)
//{
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++)
//		{
//			A[i][j] = matrix.A[i][j];
//		}
//}
//
//CMatrix3D::Matrix3d(const double *matrix)
//{
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++)
//		{
//			A[i][j] = matrix[i*4+j];
//		}
//}
//
//CMatrix3D::~CMatrix3D()
//{
//}
////operators
//CMatrix3D CMatrix3D::operator*(const CMatrix3D& matrix2) const
//{ 
//	cCc matrix;
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++){
//			matrix.A[i][j] = A[i][0]*matrix2.A[0][j]
//			+ A[i][1]*matrix2.A[1][j]
//			+ A[i][2]*matrix2.A[2][j]
//			+ A[i][3]*matrix2.A[3][j];
//		}
//		return matrix;
//}
//
//void CMatrix3D::operator*=(const CMatrix3D& matrix)
//{
//	(*this) = (*this)*matrix;
//}
////methods
//void Matrix3d::IdenticalMatrix()
//{
//	for(int i=0;i<4;i++)
//		for(int j=0;j<4;j++)
//		{
//			A[i][j] = (i==j)?1.0:0.0;
//		}
//}
//
//double Matrix3d::GetValue() const
//{
//	return  A[0][0]*A[1][1]*A[2][2] +
//		A[0][1]*A[1][2]*A[2][0] +
//		A[0][2]*A[1][0]*A[2][1] -
//		A[0][2]*A[1][1]*A[2][0] -
//		A[0][1]*A[1][0]*A[2][2] -
//		A[0][0]*A[1][2]*A[2][1];
//}
//
//// static member functions
//double Matrix3d::GetValue(double a00, double a01, double a02,
//	double a10, double a11, double a12,
//	double a20, double a21, double a22)
//{
//	return  a00*a11*a22 +
//		a01*a12*a20 +
//		a02*a10*a21 -
//		a02*a11*a20 -
//		a01*a10*a22 -
//		a00*a12*a21;
//}
//
//
//
//Matrix3d Matrix3d::CreateMirrorMatrix(Node v)
//{
//	/*double len=((Node)v).GetLength();
//	Matrix3d matrix;
//	matrix.A[0][0]= (v.getX()*v.getX() -1.0)*2.0/len/len;
//	matrix.A[1][1]= (v.getY()*v.getY() -1.0)*2.0/len/len;
//	matrix.A[2][2]= (v.getZ()*v.getZ() -1.0)*2.0/len/len;
//	matrix.A[0][1]=matrix.A[1][0]= v.getX()*v.getY()*2.0/len/len;
//	matrix.A[0][2]=matrix.A[2][0]= v.getX()*v.getZ()*2.0/len/len;
//	matrix.A[1][2]=matrix.A[2][1]= v.getZ()*v.getY()*2.0/len/len;
//	return matrix;*/
//
//
//	
//	double len=((Node)v).GetLength();
//	Matrix3d matrix;
//	matrix.A[0][0]= (1.0-2*v.getX()*v.getX() )/len/len;
//	matrix.A[1][1]= (1.0-2*v.getY()*v.getY())/len/len;
//	matrix.A[2][2]= (1.0-2*v.getZ()*v.getZ())/len/len;
//	matrix.A[0][1]=matrix.A[1][0]= -2*v.getX()*v.getY()/len/len;
//	matrix.A[0][2]=matrix.A[2][0]= -2*v.getX()*v.getZ()/len/len;
//	matrix.A[1][2]=matrix.A[2][1]= -2*v.getZ()*v.getY()/len/len;
//	return matrix;
//
//}
//
//
//
//
//
////Matrix3d Matrix3d::CreateMirrorMatrix(Node v,double d)
////{
////	
////
////	double len=v.GetLength();
////	Matrix3d matrix;
////	matrix.A[0][0]= (1.0 -2*v.getX()*v.getX());
////	matrix.A[1][1]= (1.0 -2*v.getY()*v.getY());
////	matrix.A[2][2]= (1.0 -2*v.getZ()*v.getZ());
////	matrix.A[3][3]= 1;
////	matrix.A[0][1]=matrix.A[1][0]= -v.getX()*v.getY()*2.0;
////	matrix.A[0][2]=matrix.A[2][0]= -v.getX()*v.getZ()*2.0;
////	matrix.A[1][2]=matrix.A[2][1]= -v.getZ()*v.getY()*2.0;
////	matrix.A[3][0]=-v.getX()*d*2.0;
////	matrix.A[3][1]=-v.getY()*d*2.0;
////	matrix.A[3][2]=-v.getZ()*d*2.0;
////	return matrix;
////}
//
//CMatrix3D CMatrix3D::CreateRotateMatrix(double da,Node v)
//{
//	CMatrix3D R;
//	Node bv(v);
//
//	if(IS_ZERO(da))	return R;
//
//	//ASSERT(!bv.IsZeroLength());
//
//	double lxy=bv.GetLengthXY();
//	if(IS_ZERO(lxy))
//	{
//		if(bv.getZ() < 0.0) da *= -1.0;
//		R.A[0][0]=R.A[1][1]=cos(da);
//		R.A[0][1]=sin(da);R.A[1][0]=-sin(da);
//		return R;
//	}
//	double lyz=bv.GetLengthYZ();
//	if(IS_ZERO(lyz))
//	{
//		if(bv.getX() < 0.0) da *= -1.0;
//		R.A[2][2]=R.A[1][1]=cos(da);
//		R.A[1][2]=sin(da);R.A[2][1]= -sin(da);
//		return R;
//	}
//	double lxz=bv.GetLengthZX();
//	if(IS_ZERO(lxz))
//	{
//		if(bv.getY() < 0.0) da *= -1.0;
//		R.A[0][0]=R.A[2][2]=cos(da);
//		R.A[0][2]=-sin(da);R.A[2][0]=sin(da);
//		return R;
//	}
//
//	CMatrix3D Rz;
//	Rz.A[0][0]=Rz.A[1][1]=bv.getY()/lxy;
//	Rz.A[0][1]=bv.getX()/lxy;Rz.A[1][0]= -bv.getX()/lxy;
//
//	double len=bv.GetLength();
//	CMatrix3D Rx;
//	Rx.A[2][2]=Rx.A[1][1]=bv.getZ()/len;
//	Rx.A[1][2]=lxy/len;Rx.A[2][1]= -lxy/len;
//
//	R.A[0][0]=R.A[1][1]=cos(da);
//	R.A[0][1]=sin(da);R.A[1][0]= -sin(da);
//
//	Matrix3d Rxn;
//	Rxn.A[2][2]=Rxn.A[1][1]=bv.getZ()/len;
//	Rxn.A[2][1]=lxy/len;Rxn.A[1][2]= -lxy/len;
//
//	Matrix3d Rzn;
//	Rzn.A[0][0]=Rzn.A[1][1]=bv.getY()/lxy;
//	Rzn.A[1][0]=bv.getX()/lxy;Rzn.A[0][1]= - bv.getX()/lxy;
//
//	return Rz*Rx*R*Rxn*Rzn;
//}
//
//Matrix3d Matrix3d::CreateScaleMatrix(double d)
//{
//	Matrix3d m;
//	m.A[0][0]=m.A[1][1]=m.A[2][2]=d;
//	return m;
//}
//
//Matrix3d Matrix3d::CreateTransfMatrix(Node vec)
//{
//	Matrix3d m;
//	m.A[3][0]=vec.getX();
//	m.A[3][1]=vec.getY();
//	m.A[3][2]=vec.getZ();
//	return m;
//}
//
double Node::_DistOf(Node pt)
{
	Node vec(pt.getX()-this->getX(),pt.getY()-this->getY(),pt.getZ()-this->getZ());
	return vec.GetLength();
}
bool Node::_IsParallel(Node v0,Node v1)
{
	Node cv0(v0),cv1(v1);
	return IS_ZERO((cv0*cv1).GetLength());
}

bool Node::_IsOrthogonal(Node v0,Node v1)
{
	Node cv0(v0),cv1(v1);
	return IS_ZERO(cv0|cv1);
}
double Node::_AngleBetween(Node vec)
{
	if(_IsParallel(this,vec))	
		return 0;
	Node cv1(this),cv2(vec);
	return acos((cv1|cv2.GetNormal())/cv1.GetLength());
}
//
//
//
//
//
//Node Node::CalPlaneLineIntersectPoint(Node planenode,Node planenormalvector,Node linenode,Node linedirectvector)
//{  
//	Node returnResult;  
//	double vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t,vpt;  
//	vp1 = planenormalvector.getX();
//	vp2 = planenormalvector.getY();
//	vp3 = planenormalvector.getZ();
//	n1 = planenode.getX();
//	n2 = planenode.getY();
//	n3 = planenode.getZ();
//	v1 = linedirectvector.getX();  
//	v2 = linedirectvector.getY();  
//	v3 = linedirectvector.getZ();  
//	m1 = linenode.getX();  
//	m2 = linenode.getY();  
//	m3 = linenode.getZ();
//
//
//	vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;  
//	//首先判断直线是否与平面平行  
//	if (vpt == 0)  
//	{  
//		/*returnResult = NULL;  */
//	}  
//	else  
//	{  
//		t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;  
//		returnResult.setX(m1 + v1 * t);  
//		returnResult.setY(m2 + v2 * t);
//		returnResult.setZ(m3 + v3 * t);  
//	}  
//	return returnResult;  
//
//
//}
//
//
////Node multiplycross(Node A,Node B)
////{
////	Node ret(A.getY()*B.getZ()-B.getY()*A.getZ(),A.getZ()*B.getX()-A.getX()*B.getZ(),A.getX()*B.getY()-A.getY()*B.getX());
////	return ret;
////}
////
////Node subt(Node A,Node B)
////{
////	Node ret(A.getX()-B.getX(),A.getY()-B.getY(),A.getZ()-B.getZ());
////	return ret;
////}
////
////Node normalvector(Node a,Node b,Node c)
////{
////	if (this->isNodeinline(a,b,c))//判断三点是否共线
////		exit(-1);
////	else
////		return  multiplycross(subt(a,b),subt(b,c));	
////}
////
////double vectorlen(Node p)
////{
////	return sqrt(p.getX()*p.getX()+p.getY()*p.getY()+p.getZ()*p.getZ());
////}
////
////double multiplydot(Node A,Node B)
////{
////	return A.getX()*B.getX()+ A.getY()*B.getY()+ A.getZ()*B.getZ();
////}
////
////double Nodetoplane(Node p,Node a,Node b,Node c)
////{
////	if (this->isNodeinline(a,b,c))//判断三点是否共线
////		exit(-1);
////	else
////		return fabs(multiplydot(normalvector(a,b,c),subt(p,a)))/vectorlen(normalvector(a,b,c));
////
////}
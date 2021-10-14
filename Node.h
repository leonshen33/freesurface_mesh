#pragma once
#include "CadBase.h"
#include "CMatrix.h"
#include <math.h>

#define CAD_ZERO		1.0E-6
#define NC_ZERO		1.0E-3
#define IS_ZERO(x)		(fabs(x)<=CAD_ZERO)
#define IS_NCZERO(x)		(fabs(x)<=NC_ZERO)
#define IS_BETWEEN(x,min,max) (x<=max && x>=min)

#define PI	3.1415926535
class CMatrix3D;
class Node;
//class Matrix3d
//{
//public:
//	Matrix3d();
//    Matrix3d(const Matrix3d&);
//	Matrix3d(const double *);
//	virtual ~Matrix3d();
//public:
//	double A[4][4];
//	//operators
//	Matrix3d operator*(const Matrix3d& matrix)const;
//	void operator*=(const Matrix3d& matrix);
//
//	//methods
//	void   IdenticalMatrix();
//	double GetValue() const;
//	// static member functions
//	static double GetValue(double a00, double a01, double a02,
//		double a10, double a11, double a12,
//		double a20, double a21, double a22);
//	/*static Matrix3d CreateMirrorMatrix(Node plnNorm,double d);*/
//	static Matrix3d CreateMirrorMatrix(Node plnNorm);
//	static Matrix3d CreateRotateMatrix(double da,Node bv);
//	static Matrix3d CreateScaleMatrix(double);
//	static Matrix3d CreateTransfMatrix(Node vec);
//
//} ; 

	



class Node
{
public:
	
	Node(double x,double y,double z){X=x;Y=y;Z=z;};
	Node(double *p){X=p[0];Y=p[1];Z=p[2];};
	//Node(Node &p){X=p.getX();Y=p.getY();Z=p.getY();};


	Node(){
		this->ID = 0;
		this->No = 0;
		this->FEANo = 0;
		this->Layer = 0;
		this->X = 0;
		this->Y = 0;
		this->Z = 0;
		this->Color = 0;
		this->c_group = 0;
		this->Dia = 0;
		this->Th = 0;
		this->Spec = 0;
		this->c_constraint = 0;
		this->Type = 0;
		this->Classify = 0;
		this->TimeStep = 0;
		this->MadeStep =0;
		this->CreatTimeStep = 0;
		this->Face = 0;
		this->ModelNo =0;
		this->StandardModel = 0;
		this->OriX = 0;
		this->OriY = 0;
		this->OriZ = 0;
		this->MainMember1 =0;
		this->MainMember2 = 0;
		this->Strength = 0;
		this->ST = 0;
		this->Weight = 0;
		this->MasterJoint = 0;
		this->Selected = 0;
		this->Displayed = 0;
		this->Future1 = 0;
		this->Future2 = 0;
		this->Future3 = 0;
		this->Future4 = 0;
		this->Future5 = 0;
		this->Future6 = 0;
		this->Future7 = 0;
		this->Future8 = 0;
		this->Future9 = 0;
		this->Future10 = 0;
	}
	Node(long int ID,int No,int FEANo,int Layer,double X,double Y,double Z,int Color,int c_group,double Dia,
		double Th,int Spec,int c_constraint,int Type,int Classify,int TimeStep,int MadeStep,int CreatTimeStep,int Face,int ModelNo,
		int StandardModel,double OriX,double OriY,double OriZ,int MainMember1,int MainMember2,double Strength,int ST,
		double Weight,int MasterJoint,bool Selected,bool Displayed,int Future1,int Future2,int Future3,int Future4,int Future5,double Future6,
		double Future7,double Future8,double Future9,double Future10){
			this->ID = ID;
			this->No = No;
			this->FEANo=FEANo;
			this->Layer = Layer;
			this->X = X;
			this->Y = Y;
			this->Z = Z;
			this->Color = Color;
			this->c_group = c_group;
			this->Dia = Dia;
			this->Th = Th;
			this->Spec = Spec;
			this->c_constraint = c_constraint;
			this->Type = Type;
			this->Classify = Classify;
			this->TimeStep = TimeStep;
			this->MadeStep = MadeStep;
			this->CreatTimeStep = CreatTimeStep;
			this->Face = Face;
			this->ModelNo =ModelNo;
			this->StandardModel = StandardModel;
			this->OriX = OriX;
			this->OriY = OriY;
			this->OriZ = OriZ;
			this->MainMember1 = MainMember1;
			this->MainMember2 = MainMember2;
			this->Strength = Strength;
			this->ST = ST;
			this->Weight = Weight;
			this->MasterJoint = MasterJoint;
			this->Selected = Selected;
			this->Displayed = Displayed;
			this->Future1 = Future1;
			this->Future2 = Future2;
			this->Future3 = Future3;
			this->Future4 = Future4;
			this->Future5 = Future5;
			this->Future6 = Future6;
			this->Future7 = Future7;
			this->Future8 = Future8;
			this->Future9 = Future9;
			this->Future10 = Future10;
	}
	Node(Node* node){
		this->ID = node-> ID;
		this->No = node-> No;
		this->FEANo= node->FEANo;
		this->Layer = node-> Layer;
		this->X = node-> X;
		this->Y = node-> Y;
		this->Z = node-> Z;
		this->Color = node-> Color;
		this->c_group = node-> c_group;
		this->Dia = node-> Dia;
		this->Th = node-> Th;
		this->Spec = node-> Spec;
		this->c_constraint = node-> c_constraint;
		this->Type = node-> Type;
		this->Classify = node-> Classify;
		this->TimeStep = node-> TimeStep;
		this->MadeStep = node-> MadeStep;
		this->CreatTimeStep = node-> CreatTimeStep;
		this->Face = node-> Face;
		this->ModelNo = node->ModelNo;
		this->StandardModel = node-> StandardModel;
		this->OriX = node-> OriX;
		this->OriY = node-> OriY;
		this->OriZ = node-> OriZ;
		this->MainMember1 = node-> MainMember1;
		this->MainMember2 = node-> MainMember2;
		this->Strength = node-> Strength;
		this->ST = node-> ST;
		this->Weight = node-> Weight;
		this->MasterJoint = node-> MasterJoint;
		this->Selected = node-> Selected;
		this->Displayed = node-> Displayed;
		this->Future1 = node-> Future1;
		this->Future2 = node-> Future2;
		this->Future3 = node-> Future3;
		this->Future4 = node-> Future4;
		this->Future5 = node-> Future5;
		this->Future6 = node-> Future6;
		this->Future7 = node-> Future7;
		this->Future8 = node-> Future8;
		this->Future9 = node-> Future9;
		this->Future10 = node-> Future10;
	}
	~Node(){};

	void setID(unsigned long int id) {ID=id;}
	void setNumber(long int no) { No=no;}
	void setFEANo(int feano) { FEANo=feano;}
	void setLayer(int layer) { Layer=layer;}
	void setX(double x) { X=x;}
	void setY(double y) { Y=y;}
	void setZ(double z) { Z=z;}
	void setColor(int color) { Color=color;}
	void setc_group(int c_group1) { c_group=c_group1;}
	void setDia(double dia) { Dia=dia;}
	void setTh(double th) { Th=th;}
	void setSpec(int spec) { Spec=spec;}
	void setc_constraint(int c_constraint1) { c_constraint=c_constraint1;}
	void setType(int type) { Type=type;}
	void setClassify(int classify) { Classify=classify;}
	void setTimeStep(int timedtep) { TimeStep=timedtep;}
	void setMadeStep(int madestep) { MadeStep=madestep;}
	void setCreatTimeStep(int craettimestep) { CreatTimeStep=craettimestep;}
	void setFace(int face) { Face=face;}
	void setModelNo(int modelno) { ModelNo=modelno;}
	void setStandardModel(int standardmodel) { StandardModel=standardmodel; }
	void setOriX(double orix) { OriX=orix;}
	void setOriY(double oriy) { OriY=oriy;}
	void setOriZ(double oriz) { OriZ=oriz;}
	void setMainMember1(int mainmember1) { MainMember1=mainmember1;}
	void setMainMember2(int mainmember2) { MainMember2=mainmember2;}
	void setStrength(double strength) { Strength=strength;}
	void setST(int st) { ST=st;}
	void setWeight(double weight) { Weight=weight;}
	void getMasterJoint(int masterjoint) { MasterJoint=masterjoint;}
	void setSelected(bool selected) { Selected=selected;}
	void setDisplayed(bool displayed) { Displayed=displayed;}
	void setFuture1(int future1) { Future1=future1;}
	void setFuture2(int future2) { Future2=future2;}
	void setFuture3(int future3) { Future3=future3;}
	void setFuture4(int future4) { Future4=future4;}
	void setFuture5(int future5) { Future5=future5;}
	void setFuture6(double future6) { Future6=future6;}
	void setFuture7(double future7) { Future7=future7;}
	void setFuture8(double future8) { Future8=future8;}
	void setFuture9(double future9) { Future9=future9;}
	void setFuture10(double future10) { Future10=future10;}




	long int getID() {return ID;}
	int getNumber() {return No;}
	int getFEANo() {return FEANo;}
	int getLayer() {return Layer;}
	double getX() {return X;}
	double getY() {return Y;}
	double getZ() {return Z;}
	int getColor() {return Color;}
	int getc_group() {return c_group;}
	double getDia() {return Dia;}
	double getTh() {return Th;}
	int getSpec() {return Spec;}
	int getc_constraint() {return c_constraint;}
	int getType() {return Type;}
	int getClassify() {return Classify;}
	int getTimeStep() {return TimeStep;}
	int getMadeStep() {return MadeStep;}
	int getCreatTimeStep() {return CreatTimeStep;}
	int getFace() {return Face;}
	int getModelNo() {return ModelNo;}
	int getStandardModel() {return StandardModel; }
	double getOriX() {return OriX;}
	double getOriY() {return OriY;}
	double getOriZ() {return OriZ;}
	int getMainMember1() {return MainMember1;}
	int getMainMember2() {return MainMember2;}
	double getStrength() {return Strength;}
	int getST() {return ST;}
	double getWeight() {return Weight;}
	int getMasterJoint() {return MasterJoint;}
	bool getSelected() {return Selected;}
	bool getDisplayed() {return Displayed;}
	int getFuture1() {return Future1;}
	int getFuture2() {return Future2;}
	int getFuture3() {return Future3;}
	int getFuture4() {return Future4;}
	int getFuture5() {return Future5;}
	double getFuture6() {return Future6;}
	double getFuture7() {return Future7;}
	double getFuture8() {return Future8;}
	double getFuture9() {return Future9;}
	double getFuture10() {return Future10;}


	//Node operator =(Node temp);//将一个点赋给另外一个点

	//operators
	 Node operator*(const CMatrix3D& matrix) const;
	void  operator*=(const CMatrix3D& matrix);
	Node operator*(double d) const;
	void operator*=(double d);

	Node operator/(double d) const;
	void operator/=(double d);
	//offsetting with vector
	Node operator+(Node v) const;
	void operator+=(Node v);
	Node operator-(Node v) const;
	void operator-=(Node v);
	bool operator==(Node  pos) const;
	bool operator!=(Node  pos) const;
	//cross product
	Node operator*(Node v) const;
	//dot product
	double operator|(Node v) const;
	//length
	double GetLength() const;
	double GetLengthXY() const;
	double GetLengthYZ() const;
	double GetLengthZX() const;

	Node	GetNormal() const;
	void	Normalize();
	bool	IsZeroLength() const;
	void   IdenticalMatrix();
	double GetValue() const;
	bool	/*AFX_EXT_API*/ _IsParallel(Node v0,Node v1);
	double	/*AFX_EXT_API*/ _AngleBetween(Node vec);
	double	/*AFX_EXT_API*/ _DistOf(Node pt);
	bool	/*AFX_EXT_API*/ _IsOrthogonal(Node v0,Node v1);


	Node CalPlaneLineIntersectPoint(Node planenode,Node planenormalvector,Node linenode,Node linedirectvector);

private:
	unsigned long int ID;
	long int No;
	int FEANo;
	int Layer;
	double X;
	double Y;
	double Z;
	int Color;
	int c_group;
	double Dia;
	double Th;
	int Spec;
	int c_constraint;
	int Type;
	int Classify;
	int TimeStep;
	int MadeStep;
	int CreatTimeStep;;
	int Face;
	int ModelNo;
	int StandardModel;
	double OriX;
	double OriY;
	double OriZ;
	int MainMember1;
	int MainMember2;
	double Strength;
	int ST;
	double Weight;
	int MasterJoint;
	bool Selected;
	bool Displayed;
	int Future1;
	int Future2;
	int Future3;
	int Future4;
	int Future5;
	double Future6;
	double Future7;
	double Future8;
	double Future9;
	double Future10;
};




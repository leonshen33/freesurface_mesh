#ifndef Member_h
#define Member_h
#include<string>
#include"Node.h"
class Node;
//using namespace std;
class Member{
private:
	unsigned long long int ID;
	int Number;
	int Color;
	int Layer;
	int Group;
	int Type;
	unsigned long long int Inode;
	unsigned long long int Jnode;
	unsigned long long int Knode;
	int Surface;
	int TrussNo;
	int Chord;
	int Cone;
	int DeadTimeStep;
	int MadeStep;
	int CreateTimeStep;
	int ModelNo;
	int StandardModel;
	double L2;
	double L3;
	double AllowLmd;
	double Lmd;
	double UnitW;
	int SectionType;
	int SectionIndex;
	int SectionOptIndex;
	int SectionOptMin;
	int SectionOptMax;
	bool Selected;
	bool Displayed;
	double AllowStressRatio;
	int Material;
	int ReleaseStart;
	int ReleaseEnd;
	int WrapStart;
	int WrapEnd;
	double RigidLengS;
	double RigidLengE;
	int StiffStart;
	int StiffEnd;
	double MaxTension;
	double MaxCompression;
	double DesignStressRatio;
	int LeftNutNo;
	int RightNutNo;
	int LeftBoltNo;
	int RightBoltNo;
	int LeftConeNo;
	int RightConeNo;
	int MainMember1;
	int MainMember2;
	double CutLength;
	double WeldLength;
	int MaterialLabel;
	int TrussMatLabel;
	int ParallelMember;
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
	unsigned long long int EleXSymID;
	unsigned long long int EleYSymID;
	unsigned long long int LeftChordEleID;
	unsigned long long int RightChordEleID; 
	double PreTension;
	Node*  inode;
	Node* jnode;
	//DBConn* dbconn;
	//SavePoint* savepoint;

public:
	Member(){
		this->ID = 1;
		this->Number = 0;
		this->Color = 0;
		this->Layer = 0;
		this->Group = 0;
		this->Type = 0;
		this->Inode = 0;
		this->Jnode = 0;
		this->Knode = 0;
		this->Surface = 0;
		this->TrussNo = 0;
		this->Chord = 0;
		this->Cone = 0;
		this->DeadTimeStep = 0;
		this->MadeStep = 0;
		this->CreateTimeStep = 0;
		this->ModelNo = 0;
		this->StandardModel = 0;
		this->L2 = 0.0;
		this->L3 = 0.0;
		this->AllowLmd = 0.0;
		this->Lmd = 0.0;
		this->UnitW = 0.0;
		this->SectionType = 0;
		this->SectionIndex = 0;
		this->SectionOptIndex = 0;
		this->SectionOptMin = 0;
		this->SectionOptMax = 0;
		this->Selected = false;
		this->Displayed = true;
		this->AllowStressRatio = 0.0;
		this->Material = 0;
		this->ReleaseStart = 0;
		this->ReleaseEnd = 0;
		this->WrapStart = 0;
		this->WrapEnd = 0;
		this->RigidLengS = 0.0;
		this->RigidLengE = 0.0;
		this->StiffStart = 0;
		this->StiffEnd = 0;
		this->MaxTension = 0.0;
		this->MaxCompression = 0.0;
		this->DesignStressRatio = 0.0;
		this->LeftNutNo = 0;
		this->RightNutNo = 0;
		this->LeftBoltNo = 0;
		this->RightBoltNo = 0;
		this->LeftConeNo = 0;
		this->RightConeNo = 0;
		this->MainMember1 = 0;
		this->MainMember2 = 0;
		this->CutLength = 0.0;
		this->WeldLength = 0.0;
		this->MaterialLabel = 0;
		this->TrussMatLabel = 0;
		this->ParallelMember = 0;
		this->Future1 = 0;
		this->Future2 = 0;
		this->Future3 = 0;
		this->Future4 = 0;
		this->Future5 = 0;
		this->Future6 = 0.0;
		this->Future7 = 0.0;
		this->Future8 = 0.0;
		this->Future9 = 0.0;
		this->Future10 = 0.0;
		this-> EleXSymID=0;
		this-> EleYSymID=0;
		this-> LeftChordEleID=0;
		this-> RightChordEleID=0;
		this-> PreTension=0;
		this->inode = nullptr;
		this->jnode = nullptr;
		//this->dbconn = NULL;
		//this->savepoint = NULL;
	}
	~Member()
	{
		 inode=NULL;
		 jnode=NULL;
		 //dbconn = NULL;
		// savepoint = NULL;
	}
	Member(unsigned long long int ID, int Number, int Color, int Layer, int Group, int Type, unsigned long long int Inode, unsigned long long int Jnode, unsigned long long int Knode, int Surface, int TrussNo,
		int Chord, int Cone, int DeadTimeStep, int MadeStep, int CreateTimeStep, int ModelNo, int StandardModel, double L2,
		double L3, double AllowLmd, double Lmd, double UnitW, int SectionType, int SectionIndex, int SectionOptIndex, int SectionOptMin,
		int SectionOptMax, bool Selected, bool Displayed, double AllowStressRatio, int Material, int ReleaseStart, int ReleaseEnd,
		int WrapStart, int WrapEnd, double RigidLengS, double RigidLengE, int StiffStart, int StiffEnd, double MaxTension,
		double MaxCompression, double DesignStressRatio, int LeftNutNo, int RightNutNo, int LeftBoltNo, int RightBoltNo, int LeftConeNo,
		int RightConeNo, int MainMember1, int MainMember2, double CutLength, double WeldLength, int MaterialLabel, int TrussMatLabel, int ParallelMember,
		int Future1, int Future2, int Future3, int Future4, int Future5, double Future6, double Future7, double Future8, double Future9, double Future10, unsigned long long int EleXSymID,
	    unsigned long long int EleYSymID,unsigned long long int LeftChordEleID,unsigned long long int RightChordEleID,double PreTension)
	{
		this->ID = ID;
		this->Number = Number;
		this->Color = Color;
		this->Layer = Layer;
		this->Group = Group;
		this->Type = Type;
		this->Inode = Inode;
		this->Jnode = Jnode;
		this->Knode = Knode;
		this->Surface = Surface;
		this->TrussNo = TrussNo;
		this->Chord = Chord;
		this->Cone = Cone;
		this->DeadTimeStep = DeadTimeStep;
		this->MadeStep = MadeStep;
		this->CreateTimeStep = CreateTimeStep;
		this->ModelNo = ModelNo;
		this->StandardModel = StandardModel;
		this->L2 = L2;
		this->L3 = L3;
		this->AllowLmd = AllowLmd;
		this->Lmd = Lmd;
		this->UnitW = UnitW;
		this->SectionType = SectionType;
		this->SectionIndex = SectionIndex;
		this->SectionOptIndex = SectionOptIndex;
		this->SectionOptMin = SectionOptMin;
		this->SectionOptMax = SectionOptMax;
		this->Selected = Selected;
		this->Displayed = Displayed;
		this->AllowStressRatio = AllowStressRatio;
		this->Material = Material;
		this->ReleaseStart = ReleaseStart;
		this->ReleaseEnd = ReleaseEnd;
		this->WrapStart = WrapStart;
		this->WrapEnd = WrapEnd;
		this->RigidLengS = RigidLengS;
		this->RigidLengE = RigidLengE;
		this->StiffStart = StiffStart;
		this->StiffEnd = StiffEnd;
		this->MaxTension = MaxTension;
		this->MaxCompression = MaxCompression;
		this->DesignStressRatio = DesignStressRatio;
		this->LeftNutNo = LeftNutNo;
		this->RightNutNo = RightNutNo;
		this->LeftBoltNo = LeftBoltNo;
		this->RightBoltNo = RightBoltNo;
		this->LeftConeNo = LeftConeNo;
		this->RightConeNo = RightConeNo;
		this->MainMember1 = MainMember1;
		this->MainMember2 = MainMember2;
		this->CutLength = CutLength;
		this->WeldLength = WeldLength;
		this->MaterialLabel = MaterialLabel;
		this->TrussMatLabel = TrussMatLabel;
		this->ParallelMember = ParallelMember;
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
		this->EleXSymID = EleXSymID;
		this->EleYSymID = EleYSymID;
		this->LeftChordEleID = LeftChordEleID;
		this->RightChordEleID = RightChordEleID;
		this->PreTension = PreTension;
	}
	Member(unsigned long long int ID, int Color, int Layer, int Group, unsigned long long int Inode, unsigned long long int Jnode)//获取绘图信息时所用的构造函数
	{
		this->ID = ID;
		this->Number = 0;
		this->Color = Color;
		this->Layer = Layer;
		this->Group = Group;
		this->Type = 0;
		this->Inode = Inode;
		this->Jnode = Jnode;
		this->Knode = 0;
		this->Surface = 0;
		this->TrussNo = 0;
		this->Chord = 0;
		this->Cone = 0;
		this->DeadTimeStep = 0;
		this->MadeStep = 0;
		this->CreateTimeStep = 0;
		this->ModelNo = 0;
		this->StandardModel = 0;
		this->L2 = 0.0;
		this->L3 = 0.0;
		this->AllowLmd = 0.0;
		this->Lmd = 0.0;
		this->UnitW = 0.0;
		this->SectionType = 0;
		this->SectionIndex = 0;
		this->SectionOptIndex = 0;
		this->SectionOptMin = 0;
		this->SectionOptMax = 0;
		this->Selected = false;
		this->Displayed = true;
		this->AllowStressRatio = 0.0;
		this->Material = 0;
		this->ReleaseStart = 0;
		this->ReleaseEnd = 0;
		this->WrapStart = 0;
		this->WrapEnd = 0;
		this->RigidLengS = 0.0;
		this->RigidLengE = 0.0;
		this->StiffStart = 0;
		this->StiffEnd = 0;
		this->MaxTension = 0.0;
		this->MaxCompression = 0.0;
		this->DesignStressRatio = 0.0;
		this->LeftNutNo = 0;
		this->RightNutNo = 0;
		this->LeftBoltNo = 0;
		this->RightBoltNo = 0;
		this->LeftConeNo = 0;
		this->RightConeNo = 0;
		this->MainMember1 = 0;
		this->MainMember2 = 0;
		this->CutLength = 0.0;
		this->WeldLength = 0.0;
		this->MaterialLabel = 0;
		this->TrussMatLabel = 0;
		this->ParallelMember = 0;
		this->Future1 = 0;
		this->Future2 = 0;
		this->Future3 = 0;
		this->Future4 = 0;
		this->Future5 = 0;
		this->Future6 = 0.0;
		this->Future7 = 0.0;
		this->Future8 = 0.0;
		this->Future9 = 0.0;
		this->Future10 = 0.0;
		this->EleXSymID = 0;
		this->EleYSymID = 0;
		this->LeftChordEleID = 0;
		this->RightChordEleID = 0;
		this->PreTension = 0;
	}
	Member(int Color, int Layer, int Group, unsigned long long int Inode, unsigned long long int Jnode)//新建杆件时所用的构造函数
	{
		this->ID = 0;
		this->Number = 0;
		this->Color = Color;
		this->Layer = Layer;
		this->Group = Group;
		this->Type = 0;
		this->Inode = Inode;
		this->Jnode = Jnode;
		this->Knode = 0;
		this->Surface = 0;
		this->TrussNo = 0;
		this->Chord = 0;
		this->Cone = 0;
		this->DeadTimeStep = 0;
		this->MadeStep = 0;
		this->CreateTimeStep = 0;
		this->ModelNo = 0;
		this->StandardModel = 0;
		this->L2 = 0.0;
		this->L3 = 0.0;
		this->AllowLmd = 0.0;
		this->Lmd = 0.0;
		this->UnitW = 0.0;
		this->SectionType = 0;
		this->SectionIndex = 0;
		this->SectionOptIndex = 0;
		this->SectionOptMin = 0;
		this->SectionOptMax = 0;
		this->Selected = false;
		this->Displayed = true;
		this->AllowStressRatio = 0.0;
		this->Material = 0;
		this->ReleaseStart = 0;
		this->ReleaseEnd = 0;
		this->WrapStart = 0;
		this->WrapEnd = 0;
		this->RigidLengS = 0.0;
		this->RigidLengE = 0.0;
		this->StiffStart = 0;
		this->StiffEnd = 0;
		this->MaxTension = 0.0;
		this->MaxCompression = 0.0;
		this->DesignStressRatio = 0.0;
		this->LeftNutNo = 0;
		this->RightNutNo = 0;
		this->LeftBoltNo = 0;
		this->RightBoltNo = 0;
		this->LeftConeNo = 0;
		this->RightConeNo = 0;
		this->MainMember1 = 0;
		this->MainMember2 = 0;
		this->CutLength = 0.0;
		this->WeldLength = 0.0;
		this->MaterialLabel = 0;
		this->TrussMatLabel = 0;
		this->ParallelMember = 0;
		this->Future1 = 0;
		this->Future2 = 0;
		this->Future3 = 0;
		this->Future4 = 0;
		this->Future5 = 0;
		this->Future6 = 0.0;
		this->Future7 = 0.0;
		this->Future8 = 0.0;
		this->Future9 = 0.0;
		this->Future10 = 0.0;
		this->EleXSymID = 0;
		this->EleYSymID = 0;
		this->LeftChordEleID = 0;
		this->RightChordEleID = 0;
		this->PreTension = 0;

	}

	Member(Member *member){
		this->inode = member->inode;
		this->jnode = member->jnode;
		/*this->inode = nullptr;
		this->jnode = nullptr;*/

		this->ID = member->ID;
		this->Number = member->Number;
		this->Color = member->Color;
		this->Layer = member->Layer;
		this->Group = member->Group;
		this->Type = member->Type;
		this->Inode = member->Inode;
		this->Jnode = member->Jnode;
		this->Knode = member->Knode;
		this->Surface = member->Surface;
		this->TrussNo = member->TrussNo;
		this->Chord = member->Chord;
		this->Cone = member->Cone;
		this->DeadTimeStep = member->DeadTimeStep;
		this->MadeStep = member->MadeStep;
		this->CreateTimeStep = member->CreateTimeStep;
		this->ModelNo = member->ModelNo;
		this->StandardModel = member->StandardModel;
		this->L2 = member->L2;
		this->L3 = member->L3;
		this->AllowLmd = member->AllowLmd;
		this->Lmd = member->Lmd;
		this->UnitW = member->UnitW;
		this->SectionType = member->SectionType;
		this->SectionIndex = member->SectionIndex;
		this->SectionOptIndex = member->SectionOptIndex;
		this->SectionOptMin = member->SectionOptMin;
		this->SectionOptMax = member->SectionOptMax;
		this->Selected = member->Selected;
		this->Displayed = member->Displayed;
		this->AllowStressRatio = member->AllowStressRatio;
		this->Material = member->Material;
		this->ReleaseStart = member->ReleaseStart;
		this->ReleaseEnd = member->ReleaseEnd;
		this->WrapStart = member->WrapStart;
		this->WrapEnd = member->WrapEnd;
		this->RigidLengS = member->RigidLengS;
		this->RigidLengE = member->RigidLengE;
		this->StiffStart = member->StiffStart;
		this->StiffEnd = member->StiffEnd;
		this->MaxTension = member->MaxTension;
		this->MaxCompression = member->MaxCompression;
		this->DesignStressRatio = member->DesignStressRatio;
		this->LeftNutNo = member->LeftNutNo;
		this->RightNutNo = member->RightNutNo;
		this->LeftBoltNo = member->LeftBoltNo;
		this->RightBoltNo = member->RightBoltNo;
		this->LeftConeNo = member->LeftConeNo;
		this->RightConeNo = member->RightConeNo;
		this->MainMember1 = member->MainMember1;
		this->MainMember2 = member->MainMember2;
		this->CutLength = member->CutLength;
		this->WeldLength = member->WeldLength;
		this->MaterialLabel = member->MaterialLabel;
		this->TrussMatLabel = member->TrussMatLabel;
		this->ParallelMember = member->ParallelMember;
		this->Future1 = member->Future1;
		this->Future2 = member->Future2;
		this->Future3 = member->Future3;
		this->Future4 = member->Future4;
		this->Future5 = member->Future5;
		this->Future6 = member->Future6;
		this->Future7 = member->Future7;
		this->Future8 = member->Future8;
		this->Future9 = member->Future9;
		this->Future10 = member->Future10;
		this->EleXSymID = member->EleXSymID;
		this->EleYSymID = member->EleYSymID;
		this->LeftChordEleID = member->LeftChordEleID;
		this->RightChordEleID = member->RightChordEleID;
		this->PreTension = member->PreTension;
	}

	Member(Node *theinode, Node *thejnode)
	{
		this->inode = theinode;
		this->jnode = thejnode;
		this->Inode = theinode->getID();
		this->Jnode = thejnode->getID();

		this->ID = 1;
		this->Number = 0;
		this->Color = 0;
		this->Layer = 0;
		this->Group = 0;
		this->Type = 0;
		this->Knode = 0;
		this->Surface = 0;
		this->TrussNo = 0;
		this->Chord = 0;
		this->Cone = 0;
		this->DeadTimeStep = 0;
		this->MadeStep = 0;
		this->CreateTimeStep = 0;
		this->ModelNo = 0;
		this->StandardModel = 0;
		this->L2 = 0.0;
		this->L3 = 0.0;
		this->AllowLmd = 0.0;
		this->Lmd = 0.0;
		this->UnitW = 0.0;
		this->SectionType = 0;
		this->SectionIndex = 0;
		this->SectionOptIndex = 0;
		this->SectionOptMin = 0;
		this->SectionOptMax = 0;
		this->Selected = false;
		this->Displayed = true;
		this->AllowStressRatio = 0.0;
		this->Material = 0;
		this->ReleaseStart = 0;
		this->ReleaseEnd = 0;
		this->WrapStart = 0;
		this->WrapEnd = 0;
		this->RigidLengS = 0.0;
		this->RigidLengE = 0.0;
		this->StiffStart = 0;
		this->StiffEnd = 0;
		this->MaxTension = 0.0;
		this->MaxCompression = 0.0;
		this->DesignStressRatio = 0.0;
		this->LeftNutNo = 0;
		this->RightNutNo = 0;
		this->LeftBoltNo = 0;
		this->RightBoltNo = 0;
		this->LeftConeNo = 0;
		this->RightConeNo = 0;
		this->MainMember1 = 0;
		this->MainMember2 = 0;
		this->CutLength = 0.0;
		this->WeldLength = 0.0;
		this->MaterialLabel = 0;
		this->TrussMatLabel = 0;
		this->ParallelMember = 0;
		this->Future1 = 0;
		this->Future2 = 0;
		this->Future3 = 0;
		this->Future4 = 0;
		this->Future5 = 0;
		this->Future6 = 0.0;
		this->Future7 = 0.0;
		this->Future8 = 0.0;
		this->Future9 = 0.0;
		this->Future10 = 0.0;
		this->EleXSymID = 0;
		this->EleYSymID = 0;
		this->LeftChordEleID = 0;
		this->RightChordEleID = 0;
		this->PreTension = 0;
		//this->dbconn = NULL;
		//this->savepoint = NULL;
		
	}

	unsigned long long int getID();


	int getNumber();

	int getColor();

	int getLayer();

	int getGroup();

	unsigned long long int getInode();

	int getType();

	unsigned long long int getJnode();

	unsigned long long int getKnode();

	int getSurface();

	int getTrussNo();

	int getChord();

	int getCone();

	int getDeadTimeStep();

	int getMadeStep();

	int getCreateTimeStep();

	int getModelNo();

	int getStandardModel();

	double getL2();

	double getL3();

	double  getAllowLmd();

	double getLmd();

	double getUnitW();

	int getSectionType();

	int getSectionIndex();

	int getSectionOptIndex();

	int getSectionOptMin();

	int getSectionOptMax();

	bool getSelected();

	bool getDisplayed();

	double getAllowStressRatio();

	int getMaterial();

	int getReleaseStart();

	int getReleaseEnd();

	int getWrapStart();

	int getWrapEnd();

	double getRigidLengS();

	double getRigidLengE();

	int getStiffStart();

	int getStiffEnd();

	double getMaxTension();

	double getMaxCompression();

	double getDesignStressRatio();

	int getLeftNutNo();

	int getRightNutNo();

	int getLeftBoltNo();

	int getRightBoltNo();

	int getLeftConeNo();

	int getRightConeNo();

	int getMainMember1();

	int getMainMember2();

	double getCutLength();

	double getWeldLength();

	int getMaterialLabel();

	int getTrussMatLabel();

	int getParallelMember();

	int getFutre1();

	int getFutre2();

	int getFutre3();

	int getFutre4();

	int getFutre5();

	double getFuture6();

	double getFuture7();

	double getFuture8();

	double getFuture9();

	double getFuture10();
	
	unsigned long long int getEleXSymID();

	unsigned long long int getEleYSymID();
	
	unsigned long long int getLeftChordEleID();
	
	unsigned long long int getRightChordEleID();

	double getPreTension();



	Node* getinode();

	Node* getjnode();

	//void setID(unsigned long long int ID);

	void setNumber(int Number);

	void setColor(int Color);

	void setLayer(int Layer);

	void setGroup(int Group);

	void setType(int Type);

	void setInode(unsigned long long int Inode);

	void setJnode(unsigned long long int Jnode);

	void setKnode(unsigned long long int Knode);

	void setSurface(int Surface);

	void setTrussNo(int TrussNo);

	void setChord(int Chord);

	void setCone(int Cone);

	void setDeadTimeStep(int DeadTimeStep);

	void setMadeStep(int MadeStep);

	void setCreateTimeStep(int CreateTimeStep);

	void setModelNo(int ModelNo);

	void setStandardModel(int StandardModel);

	void setL2(double L2);

	void setL3(double L3);

	void setAllowLmd(double AllowLmd);

	void setLmd(double Lmd);

	void setUnitW(double UnitW);

	void setSectionType(int SectionType);

	void setSectionIndex(int SectionIndex);

	void setSectionOptIndex(int SectionOptIndex);

	void setSectionOptMin(int SectionOptMin);

	void setSectionOptMax(int SectionOptMax);

	void setSelected(bool Selected);

	void setDisplayed(bool Displayed);

	void setAllowStressRatio(double AllowStressRatio);

	void setMaterial(int Material);

	void setReleaseStart(int ReleaseStart);

	void setReleaseEnd(int ReleaseEnd);

	void setWrapStart(int WrapStart);

	void setWrapEnd(int WrapEnd);

	void setRigidLengS(double RigidLengS);

	void setRigidLengE(double RigidLengE);

	void setStiffStart(int StiffStart);

	void setStiffEnd(int StiffEnd);

	void setMaxTension(double MaxTension);

	void setMaxCompression(double MaxCompression);

	void setDesignStressRatio(double DesignStressRatio);

	void setLeftNutNo(int LeftNutNo);

	void setRightNutNo(int RightNutNo);

	void setLeftBoltNo(int LeftBoltNo);

	void setRightBoltNo(int RightBoltNo);

	void setLeftConeNo(int LeftConeNo);

	void setRightConeNo(int RightConeNo);

	void setMainMember1(int MainMember1);

	void setMainMember2(int MainMember2);

	void setCutLength(double CutLength);

	void setWeldLength(double WeldLengt);

	void setMaterialLabel(int MaterialLabel);

	void setTrussMatLabel(int TrussMatLabel);

	void setParallelMember(int ParallelMember);

	void setFuture1(int Future1);

	void setFuture2(int Future2);

	void setFuture3(int Future3);

	void setFuture4(int Future4);

	void setFuture5(int Future5);

	void setFuture6(double Future6);

	void setFuture7(double Future7);

	void setFuture8(double Future8);

	void setFuture9(double Future9);

	void setFuture10(double Future10);

	void setEleXSymID(unsigned long long int EleXSymID);

	void setEleYSymID(unsigned long long int EleYSymID);

	void setLeftChordEleID(unsigned long long int LeftChordEleID);

	void setRightChordEleID(unsigned long long int RightChordEleID);

	void setPreTension(double PreTension);


	void set_ID(unsigned long long int ID);

	void set_Number(int Number);

	void set_Color(int Color);

	void set_Layer(int Layer);

	void set_Group(int Group);

	void set_Type(int Type);

	void set_Inode(unsigned long long int Inode);

	void set_Jnode(unsigned long long int Jnode);

	void set_Knode(unsigned long long int Knode);

	void set_Surface(int Surface);

	void set_TrussNo(int TrussNo);

	void set_Chord(int Chord);

	void set_Cone(int Cone);

	void set_DeadTimeStep(int DeadTimeStep);

	void set_MadeStep(int MadeStep);

	void set_CreateTimeStep(int CreateTimeStep);

	void set_ModelNo(int ModelNo);

	void set_StandardModel(int StandardModel);

	void set_L2(double L2);

	void set_L3(double L3);

	void set_AllowLmd(double AllowLmd);

	void set_Lmd(double Lmd);

	void set_UnitW(double UnitW);

	void set_SectionType(int SectionType);

	void set_SectionIndex(int SectionIndex);

	void set_SectionOptIndex(int SectionOptIndex);

	void set_SectionOptMin(int SectionOptMin);

	void set_SectionOptMax(int SectionOptMax);

	void set_Selected(bool Selected);

	void set_Displayed(bool Displayed);

	void set_AllowStressRatio(double AllowStressRatio);

	void set_Material(int Material);

	void set_ReleaseStart(int ReleaseStart);

	void set_ReleaseEnd(int ReleaseEnd);

	void set_WrapStart(int WrapStart);

	void set_WrapEnd(int WrapEnd);

	void set_RigidLengS(double RigidLengS);

	void set_RigidLengE(double RigidLengE);

	void set_StiffStart(int StiffStart);

	void set_StiffEnd(int StiffEnd);

	void set_MaxTension(double MaxTension);

	void set_MaxCompression(double MaxCompression);

	void set_DesignStressRatio(double DesignStressRatio);

	void set_LeftNutNo(int LeftNutNo);

	void set_RightNutNo(int RightNutNo);

	void set_LeftBoltNo(int LeftBoltNo);

	void set_RightBoltNo(int RightBoltNo);

	void set_LeftConeNo(int LeftConeNo);

	void set_RightConeNo(int RightConeNo);

	void set_MainMember1(int MainMember1);

	void set_MainMember2(int MainMember2);

	void set_CutLength(double CutLength);

	void set_WeldLength(double WeldLengt);

	void set_MaterialLabel(int MaterialLabel);

	void set_TrussMatLabel(int TrussMatLabel);

	void set_ParallelMember(int ParallelMember);

	void set_Future1(int Future1);

	void set_Future2(int Future2);

	void set_Future3(int Future3);

	void set_Future4(int Future4);

	void set_Future5(int Future5);

	void set_Future6(double Future6);

	void set_Future7(double Future7);

	void set_Future8(double Future8);

	void set_Future9(double Future9);

	void set_Future10(double Future10);
	
	void set_EleXSymID(unsigned long long int EleXSymID);

	void set_EleYSymID(unsigned long long int EleYSymID);

	void set_LeftChordEleID(unsigned long long int LeftChordEleID);

	void set_RightChordEleID(unsigned long long int RightChordEleID);

	void set_PreTension(double PreTension);

	void set_inode(Node* node);

	void set_jnode(Node* node);

	//void set_dbconn(DBConn* dbconn);

	//void set_savepoint(SavePoint* savepoint);

	std::string intTostring(int num);
	std::string doubleTostring(double num);
	std::string boolTostring(bool num);
	std::string longintTostring(int num);
	std::string unsignedlonglongintTostring(unsigned long long int num);

};
#endif

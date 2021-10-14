#include"Member.h"
#include <sstream>
unsigned long long int Member::getID()
{
	return ID;
}
int Member::getNumber()
{
	return Number;
}
int Member::getColor()
{
	return Color;
}
int Member::getLayer()
{
	return Layer;
}
int Member::getGroup()
{
	return Group;
}
int Member::getType()
{
	return Type;
}
unsigned long long int Member::getInode()
{
	return Inode;
}
unsigned long long int Member::getJnode()
{
	return Jnode;
}
unsigned long long int Member::getKnode()
{
	return Knode;
}
int Member::getSurface()
{
	return Surface;
}
int Member::getTrussNo()
{
	return TrussNo;
}
int Member::getChord()
{
	return Chord;
}
int Member::getCone()
{
	return Cone;
}
int Member::getDeadTimeStep()
{
	return DeadTimeStep;
}
int Member::getMadeStep()
{
	return MadeStep;
}
int Member::getCreateTimeStep()
{
	return CreateTimeStep;
}
int Member::getModelNo()
{
	return ModelNo;
}
int Member::getStandardModel()
{
	return StandardModel;
}
double Member::getL2()
{
	return L2;
}
double Member::getL3()
{
	return L3;
}
double Member::getAllowLmd()
{
	return AllowLmd;
}
double Member::getLmd()
{
	return Lmd;
}
double Member::getUnitW()
{
	return UnitW;
}
int Member::getSectionType()
{
	return SectionType;
}
int Member::getSectionIndex()
{
	return SectionIndex;
}
int Member::getSectionOptIndex()
{
	return SectionOptIndex;
}
int Member::getSectionOptMin()
{
	return SectionOptMin;
}
int Member::getSectionOptMax()
{
	return SectionOptMax;
}
bool Member::getSelected()
{
	return Selected;
}
bool Member::getDisplayed()
{
	return  Displayed;
}
double Member::getAllowStressRatio()
{
	return AllowStressRatio;
}
int Member::getMaterial()
{
	return Material;
}
int Member::getReleaseStart()
{
	return ReleaseStart;
}
int Member::getReleaseEnd()
{
	return ReleaseEnd;
}
int Member::getWrapStart()
{
	return WrapStart;
}
int Member::getWrapEnd()
{
	return WrapEnd;
}
double Member::getRigidLengS()
{
	return RigidLengS;
}
double Member::getRigidLengE()
{
	return RigidLengE;
}
int Member::getStiffStart()
{
	return StiffStart;
}
int Member::getStiffEnd()
{
	return StiffEnd;
}
double Member::getMaxTension()
{
	return MaxTension;
}
double Member::getMaxCompression()
{
	return MaxCompression;
}
double Member::getDesignStressRatio()
{
	return DesignStressRatio;
}
int Member::getLeftNutNo()
{
	return LeftNutNo;
}
int Member::getRightNutNo()
{
	return RightNutNo;
}
int Member::getLeftBoltNo()
{
	return LeftBoltNo;
}
int Member::getRightBoltNo()
{
	return RightBoltNo;
}
int Member::getLeftConeNo()
{
	return LeftConeNo;
}
int Member::getRightConeNo()
{
	return RightConeNo;
}
int Member::getMainMember1()
{
	return MainMember1;
}
int Member::getMainMember2()
{
	return MainMember2;
}
double Member::getCutLength()
{
	return CutLength;
}
double Member::getWeldLength()
{
	return WeldLength;
}
int Member::getMaterialLabel()
{
	return MaterialLabel;
}
int Member::getTrussMatLabel()
{
	return TrussMatLabel;
}
int Member::getParallelMember()
{
	return ParallelMember;
}
int Member::getFutre1()
{
	return Future1;
}
int Member::getFutre2()
{
	return Future2;
}
int Member::getFutre3()
{
	return Future3;
}
int Member::getFutre4()
{
	return Future4;
}
int Member::getFutre5()
{
	return Future5;
}
double Member::getFuture6()
{
	return Future6;
}
double Member::getFuture7()
{
	return Future7;
}
double Member::getFuture8()
{
	return Future8;
}
double Member::getFuture9()
{
	return Future9;
}
double Member::getFuture10()
{
	return Future10;
}
Node* Member::getinode()
{
	return inode;
}
Node* Member::getjnode()
{
	return jnode;
}
//void Member::setID(unsigned long long int ID)
//{
//	this->ID = ID;
//	std::string sql = "UPDATE member SET ID =" + unsignedlonglongintTostring(ID);
//	//->SaveSql(sql);
//	//->executeUpdate(sql);
//	
//}
void Member::setNumber(int Number)
{
	this->Number = Number;
	std::string sql = "UPDATE member SET No =" + longintTostring(Number) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setColor(int Color)
{
	this->Color = Color;
	std::string sql = "UPDATE member SET Color =" + intTostring(Color) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setLayer(int Layer)
{
	this->Layer = Layer;
	std::string sql = "UPDATE member SET Layer =" + intTostring(Layer) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setGroup(int Group)
{
	this->Group = Group;
	std::string sql = "UPDATE member SET c_Group =" + intTostring(Group) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setType(int Type)
{
	this->Type = Type;
	std::string sql = "UPDATE member SET Type =" + intTostring(Type) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setInode(unsigned long long int Inode)
{
	this->Inode = Inode;
	std::string sql = "UPDATE member SET Inode =" + unsignedlonglongintTostring(Inode) + " WHERE ID = " + unsignedlonglongintTostring(this->ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setJnode(unsigned long long int Jnode)
{
	this->Jnode = Jnode;
	std::string sql = "UPDATE member SET Jnode =" + unsignedlonglongintTostring(Jnode) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setKnode(unsigned long long int Knode)
{
	this->Knode = Knode;
	std::string sql = "UPDATE member SET Knode =" + unsignedlonglongintTostring(Knode) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setSurface(int Surface)
{
	this->Surface = Surface;
	std::string sql = "UPDATE member SET Surface =" + intTostring(Surface) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setTrussNo(int TrussNo)
{
	this->TrussNo = TrussNo;
	std::string sql = "UPDATE member SET TrussNo =" + intTostring(TrussNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	////->SaveSql(sql);
	////->executeUpdate(sql);

}
void Member::setChord(int Chord)
{
	this->Chord = Chord;
	std::string sql = "UPDATE member SET Chord  =" + intTostring(Chord) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setCone(int Cone)
{
	this->Cone = Cone;
	std::string sql = "UPDATE member SET Cone  =" + intTostring(Cone) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setDeadTimeStep(int DeadTimeStep)
{
	this->DeadTimeStep = DeadTimeStep;
	std::string sql = "UPDATE member SET DeadTimeStep  =" + intTostring(DeadTimeStep) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMadeStep(int MadeStep)
{
	this->MadeStep = MadeStep;
	std::string sql = "UPDATE member SET MadeStep  =" + intTostring(MadeStep) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setCreateTimeStep(int CreateTimeStep)
{
	this->CreateTimeStep = CreateTimeStep;
	std::string sql = "UPDATE member SET CreateTimeStep  =" + intTostring(CreateTimeStep) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setModelNo(int ModelNo)
{
	this->ModelNo = ModelNo;
	std::string sql = "UPDATE member SET ModelNo  =" + intTostring(ModelNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setStandardModel(int StandardModel)
{
	this->StandardModel = StandardModel;
	std::string sql = "UPDATE member SET StandardModel  =" + intTostring(StandardModel) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setL2(double L2)
{
	this->L2 = L2;
	std::string sql = "UPDATE member SET L2 =" + doubleTostring(L2) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setL3(double L3)
{
	this->L3 = L3;
	std::string sql = "UPDATE member SET L3 =" + doubleTostring(L3) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setAllowLmd(double AllowLmd)
{
	this->AllowLmd = AllowLmd;
	std::string sql = "UPDATE member SET AllowLmd =" + doubleTostring(AllowLmd) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setLmd(double Lmd)
{
	this->Lmd = Lmd;
	std::string sql = "UPDATE member SET Lmd =" + doubleTostring(Lmd) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setUnitW(double UnitW)
{
	this->UnitW = UnitW;
	std::string sql = "UPDATE member SET  UnitW =" + doubleTostring(UnitW) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setSectionType(int SectionType)
{
	this->SectionType = SectionType;
	std::string sql = "UPDATE member SET  SectionType =" + intTostring(SectionType) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setSectionIndex(int SectionIndex)
{
	this->SectionIndex = SectionIndex;
	std::string sql = "UPDATE member SET  SectionIndex =" + intTostring(SectionIndex) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setSectionOptIndex(int SectionOptIndex)
{
	this->SectionOptIndex = SectionOptIndex;
	std::string sql = "UPDATE member SET  SectionOptIndex=" + intTostring(SectionOptIndex) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setSectionOptMin(int SectionOptMin)
{
	this->SectionOptMin = SectionOptMin;
	std::string sql = "UPDATE member SET  SectionOptMin=" + intTostring(SectionOptMin) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setSectionOptMax(int SectionOptMax)
{
	this->SectionOptMax = SectionOptMax;
	std::string sql = "UPDATE member SET  SectionOptMax=" + intTostring(SectionOptMax) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setSelected(bool Selected)
{
	this->Selected = Selected;
	std::string sql = "UPDATE member SET  Selected=" + boolTostring(Selected) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setDisplayed(bool Displayed)
{
	this->Displayed = Displayed;
	std::string sql = "UPDATE member SET  Displayed=" + boolTostring(Displayed) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setAllowStressRatio(double AllowStressRatio)
{
	this->AllowStressRatio = AllowStressRatio;
	std::string sql = "UPDATE member SET  AllowStressRatio=" + doubleTostring(AllowStressRatio) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMaterial(int Material)
{
	this->Material = Material;
	std::string sql = "UPDATE member SET  Material=" + intTostring(Material) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setReleaseStart(int ReleaseStart)
{
	this->ReleaseStart = ReleaseStart;
	std::string sql = "UPDATE member SET  ReleaseStart=" + intTostring(ReleaseStart) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setReleaseEnd(int ReleaseEnd)
{
	this->ReleaseEnd = ReleaseEnd;
	std::string sql = "UPDATE member SET  ReleaseEnd=" + intTostring(ReleaseEnd) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setWrapStart(int WrapStart)
{
	this->WrapStart = WrapStart;
	std::string sql = "UPDATE member SET  WrapStart=" + intTostring(WrapStart) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setWrapEnd(int WrapEnd)
{
	this->WrapEnd = WrapEnd;
	std::string sql = "UPDATE member SET  WrapEnd=" + intTostring(WrapEnd) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setRigidLengS(double RigidLengS)
{
	this->RigidLengS = RigidLengS;
	std::string sql = "UPDATE member SET  RigidLengS=" + doubleTostring(RigidLengS) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setRigidLengE(double RigidLengE)
{
	this->RigidLengE = RigidLengE;
	std::string sql = "UPDATE member SET RigidLengE=" + doubleTostring(RigidLengE) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setStiffStart(int StiffStart)
{
	this->StiffStart = StiffStart;
	std::string sql = "UPDATE member SET StiffStart =" + intTostring(StiffStart) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setStiffEnd(int StiffEnd)
{
	this->StiffEnd = StiffEnd;
	std::string sql = "UPDATE member SET StiffEnd =" + intTostring(StiffEnd) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMaxTension(double MaxTension)
{
	this->MaxTension = MaxTension;
	std::string sql = "UPDATE member SET  MaxTension =" + doubleTostring(MaxTension) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMaxCompression(double MaxCompression)
{
	this->MaxCompression = MaxCompression;
	std::string sql = "UPDATE member SET MaxCompression =" + doubleTostring(MaxCompression) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setDesignStressRatio(double DesignStressRatio)
{
	this->DesignStressRatio = DesignStressRatio;
	std::string sql = "UPDATE member SET DesignStressRatio =" + doubleTostring(DesignStressRatio) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setLeftNutNo(int LeftNutNo)
{
	this->LeftNutNo = LeftNutNo;
	std::string sql = "UPDATE member SET  LeftNutNo =" + intTostring(LeftNutNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setRightNutNo(int RightNutNo)
{
	this->RightNutNo = RightNutNo;
	std::string sql = "UPDATE member SET RightNutNo =" + intTostring(RightNutNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setLeftBoltNo(int LeftBoltNo)
{
	this->LeftBoltNo = LeftBoltNo;
	std::string sql = "UPDATE member SET  LeftBoltNo =" + intTostring(LeftBoltNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setRightBoltNo(int RightBoltNo)
{
	this->RightBoltNo = RightBoltNo;
	std::string sql = "UPDATE member SET  RightBoltNo =" + intTostring(RightBoltNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setLeftConeNo(int LeftConeNo)
{
	this->LeftConeNo = LeftConeNo;
	std::string sql = "UPDATE member SET  LeftConeNo =" + intTostring(LeftConeNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setRightConeNo(int RightConeNo)
{
	this->RightConeNo = RightConeNo;
	std::string sql = "UPDATE member SET  RightConeNo =" + intTostring(RightConeNo) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMainMember1(int MainMember1)
{
	this->MainMember1 = MainMember1;
	std::string sql = "UPDATE member SET  MainMember1 =" + intTostring(MainMember1) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setMainMember2(int MainMember2)
{
	this->MainMember2 = MainMember2;
	std::string sql = "UPDATE member SET  MainMember2 =" + intTostring(MainMember2) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setCutLength(double CutLength)
{
	this->CutLength = CutLength;
	std::string sql = "UPDATE member SET CutLength =" + doubleTostring(CutLength) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setWeldLength(double WeldLengt)
{
	this->WeldLength = WeldLength;
	std::string sql = "UPDATE member SET WeldLengt =" + doubleTostring(WeldLength) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setMaterialLabel(int MaterialLabel)
{
	this->MaterialLabel = MaterialLabel;
	std::string sql = "UPDATE member SET MaterialLabel =" + intTostring(MaterialLabel) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setTrussMatLabel(int TrussMatLabel)
{
	this->TrussMatLabel = TrussMatLabel;
	std::string sql = "UPDATE member SET TrussMatLabel =" + intTostring(TrussMatLabel) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);


}
void Member::setParallelMember(int ParallelMember)
{
	this->ParallelMember = ParallelMember;
	std::string sql = "UPDATE member SET ParallelMember =" + longintTostring(ParallelMember) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture1(int Future1)
{
	this->Future1 = Future1;
	std::string sql = "UPDATE member SET Future1 =" + intTostring(Future1) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture2(int Future2)
{
	this->Future2 = Future2;
	std::string sql = "UPDATE member SET Future2 =" + intTostring(Future2) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture3(int Future3)
{
	this->Future3 = Future3;
	std::string sql = "UPDATE member SET Future3 =" + intTostring(Future3) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture4(int Future4)
{
	this->Future4 = Future4;
	std::string sql = "UPDATE member SET Future4 =" + intTostring(Future4) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture5(int Future5)
{
	this->Future5 = Future5;
	std::string sql = "UPDATE member SET Future5 =" + intTostring(Future5) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture6(double Future6)
{
	this->Future6 = Future6;
	std::string sql = "UPDATE member SET Future6 =" + doubleTostring(Future6) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture7(double Future7)
{
	this->Future7 = Future7;
	std::string sql = "UPDATE member SET Future7 =" + doubleTostring(Future7) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture8(double Future8)
{
	this->Future8 = Future8;
	std::string sql = "UPDATE member SET Future8 =" + doubleTostring(Future8) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture9(double Future9)
{
	this->Future9 = Future9;
	std::string sql = "UPDATE member SET Future9 =" + doubleTostring(Future9) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::setFuture10(double Future10)
{
	this->Future10 = Future10;
	std::string sql = "UPDATE member SET Future10 =" + doubleTostring(Future10) + " WHERE ID = " + unsignedlonglongintTostring(ID);
	//->SaveSql(sql);
	//->executeUpdate(sql);

}
void Member::set_inode(Node* node)
{
	this->inode = node;
}
void Member::set_jnode(Node* node)
{
	this->jnode = node;
}

void Member::set_ID(unsigned long long int ID)
{
	this->ID = ID;
}
void Member::set_Number(int Number)
{
	this->Number = Number;
}
void Member::set_Color(int Color)
{
	this->Color = Color;
}
void Member::set_Layer(int Layer)
{
	this->Layer = Layer;
}
void Member::set_Group(int Group)
{
	this->Group = Group;
}
void Member::set_Type(int Type)
{
	this->Type = Type;
}
void Member::set_Inode(unsigned long long int Inode)
{
	this->Inode = Inode;
}
void Member::set_Jnode(unsigned long long int Jnode)
{
	this->Jnode = Jnode;
}
void Member::set_Knode(unsigned long long int Knode)
{
	this->Knode = Knode;
}
void Member::set_Surface(int Surface)
{
	this->Surface = Surface;
}
void Member::set_TrussNo(int TrussNo)
{
	this->TrussNo = TrussNo;
}
void Member::set_Chord(int Chord)
{
	this->Chord = Chord;
}
void Member::set_Cone(int Cone)
{
	this->Cone = Cone;
}
void Member::set_DeadTimeStep(int DeadTimeStep)
{
	this->DeadTimeStep = DeadTimeStep;
}
void Member::set_MadeStep(int MadeStep)
{
	this->MadeStep = MadeStep;
}
void Member::set_CreateTimeStep(int CreateTimeStep)
{
	this->CreateTimeStep = CreateTimeStep;
}
void Member::set_ModelNo(int ModelNo)
{
	this->ModelNo = ModelNo;
}
void Member::set_StandardModel(int StandardModel)
{
	this->StandardModel = StandardModel;
}
void Member::set_L2(double L2)
{
	this->L2 = L2;
}
void Member::set_L3(double L3)
{
	this->L3 = L3;
}
void Member::set_AllowLmd(double AllowLmd)
{
	this->AllowLmd = AllowLmd;
}
void Member::set_Lmd(double Lmd)
{
	this->Lmd = Lmd;
}
void Member::set_UnitW(double UnitW)
{
	this->UnitW = UnitW;
}
void Member::set_SectionType(int SectionType)
{
	this->SectionType = SectionType;
}
void Member::set_SectionIndex(int SectionIndex)
{
	this->SectionIndex = SectionIndex;
}
void Member::set_SectionOptIndex(int SectionOptIndex)
{
	this->SectionOptIndex = SectionOptIndex;
}
void Member::set_SectionOptMin(int SectionOptMin)
{
	this->SectionOptMin = SectionOptMin;
}
void Member::set_SectionOptMax(int SectionOptMax)
{
	this->SectionOptMax = SectionOptMax;
}
void Member::set_Selected(bool Selected)
{
	this->Selected = Selected;
}
void Member::set_Displayed(bool Displayed)
{
	this->Displayed = Displayed;
}
void Member::set_AllowStressRatio(double AllowStressRatio)
{
	this->AllowStressRatio = AllowStressRatio;
}
void Member::set_Material(int Material)
{
	this->Material = Material;
}
void Member::set_ReleaseStart(int ReleaseStart)
{
	this->ReleaseStart = ReleaseStart;
}
void Member::set_ReleaseEnd(int ReleaseEnd)
{
	this->ReleaseEnd = ReleaseEnd;
}
void Member::set_WrapStart(int WrapStart)
{
	this->WrapStart = WrapStart;
}
void Member::set_WrapEnd(int WrapEnd)
{
	this->WrapEnd = WrapEnd;

}
void Member::set_RigidLengS(double RigidLengS)
{
	this->RigidLengS = RigidLengS;
}
void Member::set_RigidLengE(double RigidLengE)
{
	this->RigidLengE = RigidLengE;

}
void Member::set_StiffStart(int StiffStart)
{
	this->StiffStart = StiffStart;
}
void Member::set_StiffEnd(int StiffEnd)
{
	this->StiffEnd = StiffEnd;
}
void Member::set_MaxTension(double MaxTension)
{
	this->MaxTension = MaxTension;
}
void Member::set_MaxCompression(double MaxCompression)
{
	this->MaxCompression = MaxCompression;
}
void Member::set_DesignStressRatio(double DesignStressRatio)
{
	this->DesignStressRatio = DesignStressRatio;
}
void Member::set_LeftNutNo(int LeftNutNo)
{
	this->LeftNutNo = LeftNutNo;
}
void Member::set_RightNutNo(int RightNutNo)
{
	this->RightNutNo = RightNutNo;

}
void Member::set_LeftBoltNo(int LeftBoltNo)
{
	this->LeftBoltNo = LeftBoltNo;

}
void Member::set_RightBoltNo(int RightBoltNo)
{
	this->RightBoltNo = RightBoltNo;

}
void Member::set_LeftConeNo(int LeftConeNo)
{
	this->LeftConeNo = LeftConeNo;

}
void Member::set_RightConeNo(int RightConeNo)
{
	this->RightConeNo = RightConeNo;

}
void Member::set_MainMember1(int MainMember1)
{
	this->MainMember1 = MainMember1;


}
void Member::set_MainMember2(int MainMember2)
{
	this->MainMember2 = MainMember2;

}
void Member::set_CutLength(double CutLength)
{
	this->CutLength = CutLength;

}
void Member::set_WeldLength(double WeldLengt)
{
	this->WeldLength = WeldLength;

}
void Member::set_MaterialLabel(int MaterialLabel)
{
	this->MaterialLabel = MaterialLabel;


}
void Member::set_TrussMatLabel(int TrussMatLabel)
{
	this->TrussMatLabel = TrussMatLabel;


}
void Member::set_ParallelMember(int ParallelMember)
{
	this->ParallelMember = ParallelMember;

}
void Member::set_Future1(int Future1)
{
	this->Future1 = Future1;

}
void Member::set_Future2(int Future2)
{
	this->Future2 = Future2;

}
void Member::set_Future3(int Future3)
{
	this->Future3 = Future3;

}
void Member::set_Future4(int Future4)
{
	this->Future4 = Future4;

}
void Member::set_Future5(int Future5)
{
	this->Future5 = Future5;

}
void Member::set_Future6(double Future6)
{
	this->Future6 = Future6;

}
void Member::set_Future7(double Future7)
{
	this->Future7 = Future7;

}
void Member::set_Future8(double Future8)
{
	this->Future8 = Future8;

}
void Member::set_Future9(double Future9)
{
	this->Future9 = Future9;

}
void Member::set_Future10(double Future10)
{
	this->Future10 = Future10;

}
//void Member::set_dbconn(//* //)
//{
	//this->// = //;
//}
//void Member::set_savepoint(//* //)
//{
	//this->// = //;
//}
std::string Member::intTostring(int num)//用于int到std::string的转换
{
	std::stringstream stream;
	std::string temp;

	stream << num;
	stream >> temp;

	return temp;

}
std::string Member::doubleTostring(double num)//用于double到std::string的转换
{
	std::stringstream stream;
	std::string temp;
	stream << num;
	stream >> temp;

	return temp;
}
std::string Member::boolTostring(bool num)//用于bool到std::string的转换
{
	std::stringstream stream;
	std::string temp;
	stream << num;
	stream >> temp;

	return temp;
}
std::string Member::longintTostring(int num)
{
	std::stringstream stream;
	std::string temp;
	stream << num;
	stream >> temp;

	return temp;
}
std::string Member::unsignedlonglongintTostring(unsigned long long int num)
{
	std::stringstream stream;
	std::string temp;
	stream << num;
	stream >> temp;
	return temp;

}
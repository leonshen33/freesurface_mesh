#include <vector>
#include <stack>
#include "CadBase.h"
#include "point_nd.h"
//using namespace PLib;

//一块矩形区域，存储u v方向的最大值和最小值，以及划分网格时的类型
//type= 0，从左下到右上；
//type = 1，从右下到左上
//type = 2，从右上到左下
//type = 3，从左上到右下



typedef struct rectangleRegion
{
	double minU;	//矩形区域内最小U 值；
	double minV;	//矩形区域内最小V值；
	double maxU;	//矩形区域内最大U值；
	double maxV;	//矩形区域内最大V值；

	double minUOfEdge;	//矩形的延长边界的U方向最大值， 设置此值的目的是为了方便边界处理
	double maxUOfEdege;
	double minVOfEdge;
	double maxVOfEdge;
	int type;//划分网格时的类型
};

typedef struct BSTreeNode
{
	PLib::Point2Dd key;          //关键字  
	struct BSTreeNode *lChild;   //左孩子指针  
	struct BSTreeNode *rChild;  //右孩子指针  
	struct BSTreeNode *parent; //指向父节点指针  
	BSTreeNode(PLib::Point2Dd key0 = PLib::Point2Dd(-1, -1)) :key(key0), lChild(NULL), rChild(NULL), parent(NULL){};//初始化
}BSTreeNode, *PBSTreeNode;




//往二叉查找树中插入结点  
//插入的话，可能要改变根结点的地址，所以传的是二级指针  
void insert(PBSTreeNode * root, PLib::Point2Dd key);

void BiTree_Create(PBSTreeNode* root, PLib::Point2Dd *keyArray, int length);//构建二叉查找树

PBSTreeNode search(PBSTreeNode root, PLib::Point2Dd key);//在二叉树中查找值为key.x的点

PBSTreeNode searchMin(PBSTreeNode root);

int deleteNode(PBSTreeNode* root, PLib::Point2Dd key);//


PBSTreeNode searchMax(PBSTreeNode root);//查找最大关键字,空树时返回NULL  
 
PBSTreeNode searchPredecessor(PBSTreeNode p);//查找某个结点的前驱 

PBSTreeNode searchSuccessor(PBSTreeNode p);//查找某个结点的后继  

void PreOrder(PBSTreeNode T, rectangleRegion *rec, int num, double minUOfRegion, double maxUOfRegion, double minVOfRegion, double maxVOfRegion	);//根据二叉搜索树的前序遍历来确定划分网格的区域

int binarySearch(std::vector<double> &array, double key);//使用二分法寻找key所在数组中的位置i， 满足array[i] <= key < array[i+1];


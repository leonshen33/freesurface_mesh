#include <vector>
#include <stack>
#include "CadBase.h"
#include "point_nd.h"
//using namespace PLib;

//һ��������򣬴洢u v��������ֵ����Сֵ���Լ���������ʱ������
//type= 0�������µ����ϣ�
//type = 1�������µ�����
//type = 2�������ϵ�����
//type = 3�������ϵ�����



typedef struct rectangleRegion
{
	double minU;	//������������СU ֵ��
	double minV;	//������������СVֵ��
	double maxU;	//�������������Uֵ��
	double maxV;	//�������������Vֵ��

	double minUOfEdge;	//���ε��ӳ��߽��U�������ֵ�� ���ô�ֵ��Ŀ����Ϊ�˷���߽紦��
	double maxUOfEdege;
	double minVOfEdge;
	double maxVOfEdge;
	int type;//��������ʱ������
};

typedef struct BSTreeNode
{
	PLib::Point2Dd key;          //�ؼ���  
	struct BSTreeNode *lChild;   //����ָ��  
	struct BSTreeNode *rChild;  //�Һ���ָ��  
	struct BSTreeNode *parent; //ָ�򸸽ڵ�ָ��  
	BSTreeNode(PLib::Point2Dd key0 = PLib::Point2Dd(-1, -1)) :key(key0), lChild(NULL), rChild(NULL), parent(NULL){};//��ʼ��
}BSTreeNode, *PBSTreeNode;




//������������в�����  
//����Ļ�������Ҫ�ı�����ĵ�ַ�����Դ����Ƕ���ָ��  
void insert(PBSTreeNode * root, PLib::Point2Dd key);

void BiTree_Create(PBSTreeNode* root, PLib::Point2Dd *keyArray, int length);//�������������

PBSTreeNode search(PBSTreeNode root, PLib::Point2Dd key);//�ڶ������в���ֵΪkey.x�ĵ�

PBSTreeNode searchMin(PBSTreeNode root);

int deleteNode(PBSTreeNode* root, PLib::Point2Dd key);//


PBSTreeNode searchMax(PBSTreeNode root);//�������ؼ���,����ʱ����NULL  
 
PBSTreeNode searchPredecessor(PBSTreeNode p);//����ĳ������ǰ�� 

PBSTreeNode searchSuccessor(PBSTreeNode p);//����ĳ�����ĺ��  

void PreOrder(PBSTreeNode T, rectangleRegion *rec, int num, double minUOfRegion, double maxUOfRegion, double minVOfRegion, double maxVOfRegion	);//���ݶ�����������ǰ�������ȷ���������������

int binarySearch(std::vector<double> &array, double key);//ʹ�ö��ַ�Ѱ��key���������е�λ��i�� ����array[i] <= key < array[i+1];


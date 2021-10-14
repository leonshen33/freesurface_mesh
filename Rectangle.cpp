#include "Rectangle.h"





//������������в�����  
//����Ļ�������Ҫ�ı�����ĵ�ַ�����Դ����Ƕ���ָ��  
void insert(PBSTreeNode * root, PLib::Point2Dd key)
{
	//��ʼ��������  
	PBSTreeNode p = (PBSTreeNode)malloc(sizeof(BSTreeNode));
	p->key = key;
	p->lChild = p->rChild = p->parent = NULL;
	//����ʱ��ֱ����Ϊ�����  
	if ((*root) == NULL)
	{
		*root = p;
		return;
	}
	//���뵽��ǰ��㣨*root��������  
	if ((*root)->lChild == NULL && (*root)->key.x() > key.x()){
		p->parent = (*root);
		(*root)->lChild = p;
		return;
	}
	//���뵽��ǰ��㣨*root�����Һ���  
	if ((*root)->rChild == NULL && (*root)->key.x() < key.x()){
		p->parent = (*root);
		(*root)->rChild = p;
		return;
	}
	if ((*root)->key.x() > key.x())
		insert(&(*root)->lChild, key);
	else if ((*root)->key.x() < key.x())
		insert(&(*root)->rChild, key);
	else
		return;
}

//����Ԫ��,�ҵ����عؼ��ֵĽ��ָ�룬û�ҵ�����NULL  
PBSTreeNode search(PBSTreeNode root, PLib::Point2Dd key)
{
	if (root == NULL)
		return NULL;
	if (key.x() > root->key.x()) //����������  
		return search(root->rChild, key);
	else if (key.x() < root->key.x()) //����������  
		return search(root->lChild, key);
	else
		return root;
}

//������С�ؼ���,����ʱ����NULL  
PBSTreeNode searchMin(PBSTreeNode root)
{
	if (root == NULL)
		return NULL;
	if (root->lChild == NULL)
		return root;
	else  //һֱ�������ң�ֱ��û�����ӵĽ��  
		return searchMin(root->lChild);
}

//�������ؼ���,����ʱ����NULL  
PBSTreeNode searchMax(PBSTreeNode root)
{
	if (root == NULL)
		return NULL;
	if (root->rChild == NULL)
		return root;
	else  //һֱ���Һ����ң�ֱ��û���Һ��ӵĽ��  
		return searchMax(root->rChild);
}

//����ĳ������ǰ��  
PBSTreeNode searchPredecessor(PBSTreeNode p)
{
	//����  
	if (p == NULL)
		return p;
	//�����������������������Ǹ�  
	if (p->lChild)
		return searchMax(p->lChild);
	//��������,����ĳ��������������������  
	else{
		if (p->parent == NULL)
			return NULL;
		//����Ѱ��ǰ��  
		while (p){
			if (p->parent->rChild == p)
				break;
			p = p->parent;
		}
		return p->parent;
	}
}

//����ĳ�����ĺ��  
PBSTreeNode searchSuccessor(PBSTreeNode p)
{
	//����  
	if (p == NULL)
		return p;
	//��������������������С���Ǹ�  
	if (p->rChild)
		return searchMin(p->rChild);
	//��������,����ĳ��������������������  
	else{
		if (p->parent == NULL)
			return NULL;
		//����Ѱ�Һ��  
		while (p){
			if (p->parent->lChild == p)
				break;
			p = p->parent;
		}
		return p->parent;
	}
}

//���ݹؼ���ɾ��ĳ�����,ɾ���ɹ�����1,���򷵻�0  
//����Ѹ����ɾ������ôҪ�ı�����ĵ�ַ�����Դ�����ָ��  
int deleteNode(PBSTreeNode* root, PLib::Point2Dd key)
{
	PBSTreeNode q;
	//���ҵ�Ҫɾ���Ľ��  
	PBSTreeNode p = search(*root, key);
	PLib::Point2Dd temp;    //�ݴ��̽���ֵ  
	//û�鵽�˹ؼ���  
	if (!p)
		return 0;
	//1.��ɾ�����Ҷ�ӽ�㣬ֱ��ɾ��  
	if (p->lChild == NULL && p->rChild == NULL){
		//ֻ��һ��Ԫ�أ�ɾ��֮����һ�ſ���  
		if (p->parent == NULL){
			free(p);
			(*root) = NULL;
		}
		else{
			//ɾ���Ľ���Ǹ��ڵ������  
			if (p->parent->lChild == p)
				p->parent->lChild = NULL;
			else  //ɾ���Ľ���Ǹ��ڵ���Һ���  
				p->parent->rChild = NULL;
			free(p);
		}
	}

	//2.��ɾ���ֻ��������  
	else if (p->lChild && !(p->rChild)){
		p->lChild->parent = p->parent;
		//���ɾ���Ǹ���㣬Ҫ�ı丸�ڵ�ָ��  
		if (p->parent == NULL)
			*root = p->lChild;
		//ɾ���Ľ���Ǹ��ڵ������  
		else if (p->parent->lChild == p)
			p->parent->lChild = p->lChild;
		else //ɾ���Ľ���Ǹ��ڵ���Һ���  
			p->parent->rChild = p->lChild;
		free(p);
	}
	//3.��ɾ���ֻ���Һ���  
	else if (p->rChild && !(p->lChild)){
		p->rChild->parent = p->parent;
		//���ɾ���Ǹ���㣬Ҫ�ı丸�ڵ�ָ��  
		if (p->parent == NULL)
			*root = p->rChild;
		//ɾ���Ľ���Ǹ��ڵ������  
		else if (p->parent->lChild == p)
			p->parent->lChild = p->rChild;
		else //ɾ���Ľ���Ǹ��ڵ���Һ���  
			p->parent->rChild = p->rChild;
		free(p);
	}
	//4.��ɾ���Ľ��������ӣ������Һ���  
	//�ý��ĺ�̽��϶���������(�ο�������Һ�̽�㺯��)  
	//ɾ����̽��,��̽���ֵ����ý��  
	else{
		//�ҵ�Ҫɾ�����ĺ��  
		q = searchSuccessor(p);
		temp = q->key;
		//ɾ����̽��  
		deleteNode(root, q->key);
		p->key = temp;
	}
	return 1;
}

//����һ�ö��������  
void BiTree_Create(PBSTreeNode* root, PLib::Point2Dd *keyArray, int length)
{
	int i;
	//����������������  
	for (i = 0; i < length; i++)
		insert(root, keyArray[i]);
}


/* �������(�ǵݹ�)
˼·������T->data�󣬽�T��ջ������������������������������ʱ��ջ��Ԫ��ӦΪT����ջ�����������T����������
	�ڱ����Ĺ����н��ָ�����򻮷ֲ����浽rec������
*/
void PreOrder(PBSTreeNode T, rectangleRegion *rec, int num,double minUOfRegion, double maxUOfRegion, double minVOfRegion, double maxVOfRegion)
{
	std::stack<PBSTreeNode> stack;
	//p�Ǳ���ָ��  
	PBSTreeNode p = T;
	//ջ���ջ���p����ʱѭ��  
	int count = 0;//������Ԫ�ص��±�
	BSTreeNode *q;
	while (p || !stack.empty()){
		if (p != NULL){
			//����ջ��  
			stack.push(p);
			//���ʸ��ڵ�  
			cout << p->key; 
			rec[count].minUOfEdge = minUOfRegion;
			rec[count].maxUOfEdege = maxUOfRegion;
			rec[count].minVOfEdge = minVOfRegion;
			rec[count].maxVOfEdge = maxVOfRegion;
			if (p->parent == NULL)//��һ����
			{
				rec[count].minU = p->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = maxUOfRegion;
				rec[count].maxV = maxVOfRegion;
				rec[count].type = 0;			//�����µ�����Ѱ��
				count++;

				rec[count].minU = minUOfRegion;
				rec[count].minV = p->key.y();
				rec[count].maxU = p->key.x();
				rec[count].maxV = maxVOfRegion;
				rec[count].type = 1;			//�����µ�����Ѱ��
				count++;

				if (p->lChild == NULL)
				{
					rec[count].minU = minUOfRegion;
					rec[count].minV = minVOfRegion;
					rec[count].maxU = p->key.x();
					rec[count].maxV = p->key.y();
					rec[count].type = 2;			//�����ϵ�����Ѱ��
					count++;
				}
				if (p->rChild == NULL)
				{
					rec[count].minU = p->key.x();
					rec[count].minV = minVOfRegion;
					rec[count].maxU = maxUOfRegion;
					rec[count].maxV = p->key.y();
					rec[count].type = 3;			//�����µ�����Ѱ��
					count++;
				}
			}
			else if (p->parent->lChild == p)	//pΪ�������ϵĽڵ�
			{
				rec[count].minU = p->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = p->parent->key.x();
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 0;			//�����µ�����Ѱ��
				count++;

				q = p->parent;
				while (q != NULL)
				{
					if (q->key.x() < p->key.x())
					{
						break;
					}
					q = q->parent;
				}
				if (q == NULL)
				{
					rec[count].minU = minUOfRegion;
				}
				else
				{
					rec[count].minU = q->key.x();
				}
				rec[count].minV = p->key.y();
				rec[count].maxU = p->key.x();
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 1;			//�����µ�����Ѱ��
				count++;


				if (p->lChild == NULL)
				{
					/*q = p->parent;
					while (q != NULL)
					{
						if (q->key.x() < p->key.x())
						{
							break;
						}
						q = q->parent;
					}		*/
					if (q==NULL)
					{
						rec[count].minU = minUOfRegion;
					}
					else
					{
						rec[count].minU = q->key.x();
					}
					//rec[count].minV = minVOfRegion;
					rec[count].maxU = p->key.x();
					rec[count].maxV = p->key.y();
					rec[count].type = 2;			//�����µ�����Ѱ��
					count++;
				}
				if (p->rChild == NULL)
				{
					rec[count].minU = p->key.x();
					rec[count].minV = minVOfRegion;
					q = p->parent;
					while (q != NULL)
					{
						if (q->key.x() > p->key.x())
						{
							break;
						}
						q = q->parent;
					}
					if (q == NULL)
					{
						rec[count].maxU = maxUOfRegion;
					}
					else
					{
						rec[count].maxU = q->key.x();
					}
					//rec[count].maxU = p->parent->key.x();
					rec[count].maxV = p->key.y();
					rec[count].type = 3;			//�����µ�����Ѱ��
					count++;					
				}
			}
			else if (p->parent->rChild == p)	//pΪ�������ϵĽڵ�
			{
				//�����µ�����Ѱ��
				rec[count].minU = p->key.x();
				rec[count].minV = p->key.y();
				q = p->parent;
				while (q != NULL)
				{
					if (q->key.x() > p->key.x())
					{
						break;
					}
					q = q->parent;
				}
				if (q == NULL)
				{
					rec[count].maxU = maxUOfRegion;
				}
				else
				{
					rec[count].maxU = q->key.x();
				}
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 0;			
				count++;

				//�����µ�����Ѱ��
				rec[count].minU = p->parent->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = p->key.x();
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 1;			
				count++;

				//�����ϵ�����Ѱ��
				if (p->lChild == NULL)
				{				
					q = p->parent;
					while (q != NULL)
					{
						if (q->key.x() < p->key.x())
						{
							break;
						}
						q = q->parent;
					}
					if (q == NULL)
					{
						rec[count].minU = minUOfRegion;
					}
					else
					{
						rec[count].minU = q->key.x();
					}
					rec[count].minU = p->parent->key.x();
					rec[count].minV = minVOfRegion;
					rec[count].maxU = p->key.x();
					rec[count].maxV = p->key.y();
					rec[count].type = 2;			
					count++;
				}
				if (p->rChild == NULL)
				{
					//�����ϵ�����Ѱ��
					rec[count].minU = p->key.x();
					rec[count].minV = minVOfRegion;
					q = p->parent;
					while (q != NULL)
					{
						if (q->key.x() > p->key.x())
						{
							break;
						}
						q = q->parent;
					}
					if (q == NULL)
					{
						rec[count].maxU = maxUOfRegion;
					}
					else
					{
						rec[count].maxU = q->key.x();
					}
					rec[count].maxV = p->key.y();
					rec[count].type = 3;			//�����µ�����Ѱ��
					count++;
									
				}

				
			}
			//����������  
			p = p->lChild;
		}
		else{
			//��ջ  
			p = stack.top();
			stack.pop();
			//����������  
			p = p->rChild;
		}
	}//while  
}


int binarySearch(std::vector<double> &array, double key)	//ʹ�ö��ַ�Ѱ��key���������е�λ��i�� ����array[i] <= key < array[i+1]
{
	int length = array.size();
	if (length <= 0)
	{
		std::cerr << "����Ϊ�գ��޷����Ҷ�Ӧ��keyֵ��" << endl;
	}
	if (array[length-1] < key)
	{
		return  length - 1;
	}
	if (array[0] > key)
	{
		return  0;
	}
	
	int low = 0, high = length;
	int mid = (low + high) / 2;
	while (low <= high)
	{
		if (array[mid] ==  key)
		{
			return mid;
		}
		if (array[mid] <  key)
		{
			low = mid+1;
			mid = (low + high) / 2;
		}
		else
		{
			high = mid -1;
			mid = (low + high) / 2;
		}

	}
	return mid;
}

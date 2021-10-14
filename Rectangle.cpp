#include "Rectangle.h"





//往二叉查找树中插入结点  
//插入的话，可能要改变根结点的地址，所以传的是二级指针  
void insert(PBSTreeNode * root, PLib::Point2Dd key)
{
	//初始化插入结点  
	PBSTreeNode p = (PBSTreeNode)malloc(sizeof(BSTreeNode));
	p->key = key;
	p->lChild = p->rChild = p->parent = NULL;
	//空树时，直接作为根结点  
	if ((*root) == NULL)
	{
		*root = p;
		return;
	}
	//插入到当前结点（*root）的左孩子  
	if ((*root)->lChild == NULL && (*root)->key.x() > key.x()){
		p->parent = (*root);
		(*root)->lChild = p;
		return;
	}
	//插入到当前结点（*root）的右孩子  
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

//查找元素,找到返回关键字的结点指针，没找到返回NULL  
PBSTreeNode search(PBSTreeNode root, PLib::Point2Dd key)
{
	if (root == NULL)
		return NULL;
	if (key.x() > root->key.x()) //查找右子树  
		return search(root->rChild, key);
	else if (key.x() < root->key.x()) //查找左子树  
		return search(root->lChild, key);
	else
		return root;
}

//查找最小关键字,空树时返回NULL  
PBSTreeNode searchMin(PBSTreeNode root)
{
	if (root == NULL)
		return NULL;
	if (root->lChild == NULL)
		return root;
	else  //一直往左孩子找，直到没有左孩子的结点  
		return searchMin(root->lChild);
}

//查找最大关键字,空树时返回NULL  
PBSTreeNode searchMax(PBSTreeNode root)
{
	if (root == NULL)
		return NULL;
	if (root->rChild == NULL)
		return root;
	else  //一直往右孩子找，直到没有右孩子的结点  
		return searchMax(root->rChild);
}

//查找某个结点的前驱  
PBSTreeNode searchPredecessor(PBSTreeNode p)
{
	//空树  
	if (p == NULL)
		return p;
	//有左子树、左子树中最大的那个  
	if (p->lChild)
		return searchMax(p->lChild);
	//无左子树,查找某个结点的右子树遍历完了  
	else{
		if (p->parent == NULL)
			return NULL;
		//向上寻找前驱  
		while (p){
			if (p->parent->rChild == p)
				break;
			p = p->parent;
		}
		return p->parent;
	}
}

//查找某个结点的后继  
PBSTreeNode searchSuccessor(PBSTreeNode p)
{
	//空树  
	if (p == NULL)
		return p;
	//有右子树、右子树中最小的那个  
	if (p->rChild)
		return searchMin(p->rChild);
	//无右子树,查找某个结点的左子树遍历完了  
	else{
		if (p->parent == NULL)
			return NULL;
		//向上寻找后继  
		while (p){
			if (p->parent->lChild == p)
				break;
			p = p->parent;
		}
		return p->parent;
	}
}

//根据关键字删除某个结点,删除成功返回1,否则返回0  
//如果把根结点删掉，那么要改变根结点的地址，所以传二级指针  
int deleteNode(PBSTreeNode* root, PLib::Point2Dd key)
{
	PBSTreeNode q;
	//查找到要删除的结点  
	PBSTreeNode p = search(*root, key);
	PLib::Point2Dd temp;    //暂存后继结点的值  
	//没查到此关键字  
	if (!p)
		return 0;
	//1.被删结点是叶子结点，直接删除  
	if (p->lChild == NULL && p->rChild == NULL){
		//只有一个元素，删完之后变成一颗空树  
		if (p->parent == NULL){
			free(p);
			(*root) = NULL;
		}
		else{
			//删除的结点是父节点的左孩子  
			if (p->parent->lChild == p)
				p->parent->lChild = NULL;
			else  //删除的结点是父节点的右孩子  
				p->parent->rChild = NULL;
			free(p);
		}
	}

	//2.被删结点只有左子树  
	else if (p->lChild && !(p->rChild)){
		p->lChild->parent = p->parent;
		//如果删除是父结点，要改变父节点指针  
		if (p->parent == NULL)
			*root = p->lChild;
		//删除的结点是父节点的左孩子  
		else if (p->parent->lChild == p)
			p->parent->lChild = p->lChild;
		else //删除的结点是父节点的右孩子  
			p->parent->rChild = p->lChild;
		free(p);
	}
	//3.被删结点只有右孩子  
	else if (p->rChild && !(p->lChild)){
		p->rChild->parent = p->parent;
		//如果删除是父结点，要改变父节点指针  
		if (p->parent == NULL)
			*root = p->rChild;
		//删除的结点是父节点的左孩子  
		else if (p->parent->lChild == p)
			p->parent->lChild = p->rChild;
		else //删除的结点是父节点的右孩子  
			p->parent->rChild = p->rChild;
		free(p);
	}
	//4.被删除的结点既有左孩子，又有右孩子  
	//该结点的后继结点肯定无左子树(参考上面查找后继结点函数)  
	//删掉后继结点,后继结点的值代替该结点  
	else{
		//找到要删除结点的后继  
		q = searchSuccessor(p);
		temp = q->key;
		//删除后继结点  
		deleteNode(root, q->key);
		p->key = temp;
	}
	return 1;
}

//创建一棵二叉查找树  
void BiTree_Create(PBSTreeNode* root, PLib::Point2Dd *keyArray, int length)
{
	int i;
	//逐个结点插入二叉树中  
	for (i = 0; i < length; i++)
		insert(root, keyArray[i]);
}


/* 先序遍历(非递归)
思路：访问T->data后，将T入栈，遍历左子树；遍历完左子树返回时，栈顶元素应为T，出栈，再先序遍历T的右子树。
	在遍历的过程中将分割的区域划分并保存到rec数组中
*/
void PreOrder(PBSTreeNode T, rectangleRegion *rec, int num,double minUOfRegion, double maxUOfRegion, double minVOfRegion, double maxVOfRegion)
{
	std::stack<PBSTreeNode> stack;
	//p是遍历指针  
	PBSTreeNode p = T;
	//栈不空或者p不空时循环  
	int count = 0;//数组中元素的下标
	BSTreeNode *q;
	while (p || !stack.empty()){
		if (p != NULL){
			//存入栈中  
			stack.push(p);
			//访问根节点  
			cout << p->key; 
			rec[count].minUOfEdge = minUOfRegion;
			rec[count].maxUOfEdege = maxUOfRegion;
			rec[count].minVOfEdge = minVOfRegion;
			rec[count].maxVOfEdge = maxVOfRegion;
			if (p->parent == NULL)//第一个点
			{
				rec[count].minU = p->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = maxUOfRegion;
				rec[count].maxV = maxVOfRegion;
				rec[count].type = 0;			//从左下到右上寻找
				count++;

				rec[count].minU = minUOfRegion;
				rec[count].minV = p->key.y();
				rec[count].maxU = p->key.x();
				rec[count].maxV = maxVOfRegion;
				rec[count].type = 1;			//从右下到左上寻找
				count++;

				if (p->lChild == NULL)
				{
					rec[count].minU = minUOfRegion;
					rec[count].minV = minVOfRegion;
					rec[count].maxU = p->key.x();
					rec[count].maxV = p->key.y();
					rec[count].type = 2;			//从右上到左下寻找
					count++;
				}
				if (p->rChild == NULL)
				{
					rec[count].minU = p->key.x();
					rec[count].minV = minVOfRegion;
					rec[count].maxU = maxUOfRegion;
					rec[count].maxV = p->key.y();
					rec[count].type = 3;			//从右下到左上寻找
					count++;
				}
			}
			else if (p->parent->lChild == p)	//p为左子树上的节点
			{
				rec[count].minU = p->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = p->parent->key.x();
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 0;			//从左下到右上寻找
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
				rec[count].type = 1;			//从左下到右上寻找
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
					rec[count].type = 2;			//从左下到右上寻找
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
					rec[count].type = 3;			//从左下到右上寻找
					count++;					
				}
			}
			else if (p->parent->rChild == p)	//p为右子树上的节点
			{
				//从左下到右上寻找
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

				//从左下到右上寻找
				rec[count].minU = p->parent->key.x();
				rec[count].minV = p->key.y();
				rec[count].maxU = p->key.x();
				rec[count].maxV = p->parent->key.y();
				rec[count].type = 1;			
				count++;

				//从右上到左下寻找
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
					//从左上到右下寻找
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
					rec[count].type = 3;			//从左下到右上寻找
					count++;
									
				}

				
			}
			//遍历左子树  
			p = p->lChild;
		}
		else{
			//退栈  
			p = stack.top();
			stack.pop();
			//访问右子树  
			p = p->rChild;
		}
	}//while  
}


int binarySearch(std::vector<double> &array, double key)	//使用二分法寻找key所在数组中的位置i， 满足array[i] <= key < array[i+1]
{
	int length = array.size();
	if (length <= 0)
	{
		std::cerr << "容器为空，无法查找对应的key值！" << endl;
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

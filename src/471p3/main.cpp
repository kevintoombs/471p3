#include "includes.h"
#include "genomics.h"
#include "dynamicTable.h"
#include "suffixTree.h"

using namespace std;

int PREPARE_NEXTINDEX;
int MINLENGTH = -1;

void PrepareST(McSuffixTree sT);
void DFS_PrepareST(Node* T, int A[]);
Node* lastSibling(Node* n);

int main(int argc, char *argv[])
{
	//The below function shows the original functionality from program 1.
	//DP_table::demo(argv[1], argv[3], argv[2]);

	//The two functions below are as follows
	//1: creates a dp table from arguments
	//2: demos that table like in program 1
	//DP_table t(argv[1], argv[3], atoi(argv[2]));
	//t.demoTable();

	//The below function shows the original functionality from program 2.
	//McSuffixTree::demo(argv[4], argv[5]);

	//The two functions below are as follows
	//1: creates a st from arguments
	//2: demos that s2 like in program 2
	McSuffixTree mST(argv[4], argv[5]);
	//mST.demoTree();
	PrepareST(mST);

	cin.ignore();
	return 0;
	
}


void PrepareST(McSuffixTree sT)
{
	//1. Create an array A of size n(string length of input + 1 for $), and initialize content with - 1.
	int ci = sT.s.length() + 1;
	int* A = new int[ci];
	fill_n(A, ci, -1);

	//2. Initialize a global integer variable nextIndex to the start of the A array(i.e., in all my pseudocodes I have been using start indices as 1. Please initialize this to 0 in your real code.) This variable represents the next index in A which has to be populated.
	PREPARE_NEXTINDEX = 0;

	//2. Call DFS_PrepareST(root, A); // where root is the root of the suffix tree of the reference genome
	DFS_PrepareST(sT.root, A);
	
	return;
	//delete[] A;
}

void DFS_PrepareST(Node* T, int A[])
{
	if (T == NULL) return;
	
	if (T->child == NULL) {         // case: T is a leaf node
		A[PREPARE_NEXTINDEX] = T->suffixID; //suffix ID of this leaf node;
		if (McSuffixTree::deep(T) >= MINLENGTH) {
			T->start_leaf_index = PREPARE_NEXTINDEX;
			T->end_leaf_index = PREPARE_NEXTINDEX;
		}
		PREPARE_NEXTINDEX++;
		DFS_PrepareST(T->sibling, A);
		return;

	}
	else //case: T is an internal node
	{
		DFS_PrepareST(T->child, A);
		DFS_PrepareST(T->sibling, A);
		if (McSuffixTree::deep(T) >= MINLENGTH) {
			Node* u_left = T->child;
			Node* u_right = lastSibling(T->child);
			T->start_leaf_index = u_left->start_leaf_index;
			T->end_leaf_index = u_right->end_leaf_index;
		}
	}
}

Node* lastSibling(Node* n)
{
	if (n->sibling == NULL)
		return n;
	else
		return lastSibling(n->sibling);
}
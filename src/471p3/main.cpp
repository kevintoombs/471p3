#include "includes.h"
#include "genomics.h"
#include "dynamicTable.h"
#include "suffixTree.h"

using namespace std;

int PREPARE_NEXTINDEX;
int MINLENGTH = -1;

string READS[1000000];
int A[6000000];
McSuffixTree* ST;

void PrepareST();
void mapReads();
void DFS_PrepareST(Node* T);
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
	ST = new McSuffixTree(argv[4], argv[5]);
	//mST.demoTree();
	PrepareST();

	for (int i = 0; i < 10; i++)
	{
		cout << A[i] << endl;
	}

	cin.ignore();
	return 0;
	
}


void PrepareST()
{
	//1. Create an array A of size n(string length of input + 1 for $), and initialize content with - 1.
	int ci = ST->s.length() + 1;

	fill_n(A, ci, -1);

	//2. Initialize a global integer variable nextIndex to the start of the A array(i.e., in all my pseudocodes I have been using start indices as 1. Please initialize this to 0 in your real code.) This variable represents the next index in A which has to be populated.
	PREPARE_NEXTINDEX = 0;

	//2. Call DFS_PrepareST(root, A); // where root is the root of the suffix tree of the reference genome
	DFS_PrepareST(ST->root);

	mapReads();
	
	return;
	//delete[] A;
}

void DFS_PrepareST(Node* T)
{
	if (T == NULL) return;
	
	if (T->child == NULL) {         // case: T is a leaf node
		A[PREPARE_NEXTINDEX] = T->suffixID; //suffix ID of this leaf node;
		if (McSuffixTree::deep(T) >= MINLENGTH) {
			T->start_leaf_index = PREPARE_NEXTINDEX;
			T->end_leaf_index = PREPARE_NEXTINDEX;
		}
		PREPARE_NEXTINDEX++;
		DFS_PrepareST(T->sibling);
		return;

	}
	else //case: T is an internal node
	{
		DFS_PrepareST(T->child);
		DFS_PrepareST(T->sibling);
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

void mapReads()
{

	/*
	Step 3) (MapReads) 
    For i=1 to m do {

    Step 3a) Let ri be the ith read. Let l denote the length of ri.

    Step 3b) (FindLoc) 
                  Find the set Li of all locations on the genome G that share a longest common substring of length >=x characters with the read ri. For example, if the reference genome sequence is the string "accgaccgtact" and the read is "tacaccg", then the longest common substring between the read and the reference is "accg", which occurs starting from positions 1 and 5 along the reference genome. So the Li for this read should be output as {1,5}, assuming x is 2 or 3.
                   The implementation of this step will use the suffix tree constructed in Step 1.  x is a parameter to the code, and in practice it can be calculated as a function of both the read length and the estimated error rate of the sequencing process that led to the generation of the reads. For this project, however, just use x=25 in all your experiments. 
                  The proposed algorithm to do this "FindLoc" is elaborated below.
                  Note that the set Li represents a candidate list of all indices j's along the reference genome G which are starting positions for the longest common exact match of length >=x characters between ri and G. (This also implies, 1≤j≤N-x+1). 
                  Let the number of identified indices along G for ri be |Li|.

    Step 3c) (Align) 
                   For each j Î Li  {
                            i) Extract substring G[j-l... j+l], where l is the length of the read.  (Of course, make sure you handle boundary cases here - i.e., if j is at the beginning or ending parts of G, then you should correspondingly retrieve as many characters that exist in G without going out of bounds.).

                            ii) Perform a local alignment computation (using Smith-Waterman) between read ri and G[j-l... j+l]. For the alignment, you can use following parameters: {m_a =+1, m_i=-2, h=-5, g=-1}.  After computing the DP table, record two pieces of information corresponding to the computed Optimal Local Alignment: a) the number of matches (#matches), and b) the alignment length (#alignlen), which is nothing but the number of aligned columns in your final alignment (should be equal to #matches + #mismatches + #gaps). One simple way to calculate these two values (#matches, #alignlen) will be to simply call the optimal path traceback function (without doing a display of the path) and calculate  the numbers from there. (PS: There is actually second, more efficient way that will allow you to calculate these numbers during the forward computation of the DP table itself.)

                            iii) Let PercentIdentity = #matches/#alignlen.

                            iv) Let LengthCoverage = #alignlen / l .

                            v) If (PercentIdentity ≥ X% AND LengthCoverage ≥ Y%) {
                                        // that means, the alignment is of "good" quality. All we need to do now is to keep track of and output the best quality alignment. 
                                        // Note: Here, X and Y are user supplied parameters. By default, set X=90%, Y=80%.
                                        Is this value of LengthCoverage > previously seen best value for LengthCoverage for this read from any other value of j seen so far? 
                                        If so, update the new value of LengthCoverage and also the record this as the "best hit". Basically, "best hit" = (j0,j1), which are the start and end indices of the substring on the reference genome sequence corresponding to this optimal local alignment.

                                }

                    } // end for j

    Step 4) (Output)
                Output the best hit calculated from Step 3c.v as the hit for read ri. Your output can be simply: <Read_name> <Start index of hit> <End index of hit>.

                If no hit was identified for this read, output <Read_name>  "No hit found".

          } // end for i               

********************** ReadMapping: MAIN end ***********************
	*/
}
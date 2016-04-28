#include "includes.h"
#include "genomics.h"
#include "dynamicTable.h"
#include "suffixTree.h"

using namespace std;

int PREPARE_NEXTINDEX;
int MINLENGTH = 15; //25
Sequence READS[1834925];
int NUMREADS;
int NUMALIGNS;
int A[6000000];
McSuffixTree* ST;
Node *DEEPESTNODE;
int X = 90;
int Y = 80;
ofstream OUT;

int main(int argc, char * argv[]);

void PrepareST();
void mapReads();
void DFS_PrepareST(Node* T);
Node* lastSibling(Node* n);
void findLoc(int i);
void align(int i);
void output(int i, int bestStart, int bestEnd);

int newFindPath(Node *&N, int i, int read_ptr);

int main(int argc, char *argv[])
{
	OUT.open("MappingResults_fasta.txt", ios::app);

	unsigned int a = clock();
	ST = new McSuffixTree(argv[1], argv[2]);
	cout << "The tree took: " << clock() - a << "ms to build including all memory allocation for the program." << endl;
	

	unsigned int b = clock();
	PrepareST();
	cout << "The tree took: " << clock() - b << "ms to prepare" << endl;

	unsigned int d = clock();
	Sequence::parseFastaIntoReads(argv[3]);
	cout << "It took " << clock() - d << "ms to parse " << NUMREADS << " reads." << endl;

	
	unsigned int c = clock();
	mapReads();
	OUT.close();
	cout << "It took: " << clock() - c << "ms to map " << NUMREADS << " reads." << endl;
	cout << "Final execution time: " << clock() - a << " ms." << endl;
	cout << endl;

	cout << "Reference genome length: " << ST->s.length() << endl;
	cout << "Reads: " << NUMREADS << endl;
	cout << "ms/Read: " << (float)c / (float)NUMREADS << endl;
	cout << "Average alignments per read: " << (float)NUMALIGNS / (float)NUMREADS << endl;
	cout << endl;

	cout << "Press enter to exit." << endl;
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
	
	return;
	//delete[] A;
}

void DFS_PrepareST(Node* T)
{
	if (T == NULL) return;
	T->deep = McSuffixTree::deep(T);
	if (T->child == NULL) {         // case: T is a leaf node
		A[PREPARE_NEXTINDEX] = T->suffixID; //suffix ID of this leaf node;
		if (T->deep >= MINLENGTH) {
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
		if (T->deep >= MINLENGTH) {
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
	for (int i = 0; i < 10000 /*NUMREADS*/; i++)
	{
		if (i % 100 == 0)cout << "!";
		//cout << i << "/" << NUMREADS << endl;
		//unsigned int a = clock();
		findLoc(i);
		//cout << "find loc: " << clock() - a << "ms" << endl;
		//a = clock();
		align(i);
		//cout << "align loc: " << clock() - a << "ms" << endl;
	}
}

void findLoc(int i)
{
	Node* N = ST->root;
	int read_ptr = 0;

	int deepest = 0;

	//unsigned int a = clock();
	while (read_ptr < READS[i].seq.length())
	{
		//printf("new find path about to be called with %i, %i\n", i, read_ptr);
		int new_read_ptr = newFindPath(N, i, read_ptr);
		//printf("new find path returned  with %i\n",new_read_ptr);
		if (read_ptr == new_read_ptr)
		{
			read_ptr++;
		}
		else
			read_ptr = new_read_ptr;
		

		int deep = McSuffixTree::deep(N);
		if (deep >= MINLENGTH && deep >= deepest)
		{
			deepest = deep;
			DEEPESTNODE = N;
		}
		N = N->sL;
	}
	//cout << "newFindPath loop inside find loc took: " << clock() - a << "ms" << endl;
	return;

}

void align(int i)
{
	if (!DEEPESTNODE)
	{
		//cout << READS[i].header << " No hit found.\n" << endl;
		OUT << READS[i].header + ": No hit found.\n";
		return;
	}
	config c;
	c.matchScore;

	float bestCoverage = 0;
	int bestStart;
	int bestEnd;
	//cout << "INDEXES: ";
	int a = DEEPESTNODE->start_leaf_index;
	int b = DEEPESTNODE->end_leaf_index;
	for (int j = a; j <= b; j++)
	{
		int l = READS[i].seq.length();
		int start = fmax(A[j] - l, 0);
		int length = fmin(l * 2 + 1, ST->s.length()-1-start);
		//cout << A[j] << ST->s[A[j]] << j << "[" << ST->s.substr(start, length) << "]" << start << "/" << length << ", ";
		NUMALIGNS++;
		report r = DP_table::align(ST->s.substr(start, length), READS[i].seq);
		
		//printf("REPORT: %f %f\n", r.lengthCoverage, r.PercentIdentity);

		if (r.lengthCoverage >= bestCoverage)
		{
			bestCoverage = r.lengthCoverage;
			bestStart = start;
			bestEnd = start + length;
		}
		//cout << bestStart << " - " << bestEnd << endl;
	}
	//cout << endl;
	output(i, bestStart, bestEnd);
	DEEPESTNODE = NULL;
	return;
}
	
void output(int i, int bestStart, int bestEnd)
{
	//cout << READS[i].header << " mapped to : G[" << bestStart << " ... " << bestEnd << "]" << endl;
	string s = READS[i].header;
	s += " mapped to : G[";
	s += to_string(bestStart);
	s += " ... ";
	s += to_string(bestEnd);
	s += "]\n";
	OUT << s;
	return;
}

int newFindPath(Node*& T, int I, int read_ptr)
{
	//printf("new find path called with %i, %i\n", I, read_ptr);
	//check the children of n
	Node *u = T->child;
	string s = ST->s;
 	int sumI = 0;

	//if n has any children or hasn't run out of possiblities
	while (u != NULL)
	{
		//if the first character of the child matches matches...
		if (s[u->startIndex - 1] == READS[I].seq[read_ptr + sumI])
		{
			//set matches to 1 
			int i = 1;
			//check for more matches
			while ((s[u->startIndex - 1 + i] == READS[I].seq[read_ptr + sumI + i]) && (i < u->stringSize))
			{
				i++;
				if (read_ptr + sumI + i == READS[I].seq.length())
				{
					//matches sumI + i
					//at end of read
					//printf("Total matches this time: %i", sumI);
					return read_ptr + sumI + i;
				}
			}
			sumI += i;
			//if all characters matched
			if (i == u->stringSize)
			{
				//set this child to the new parent
				T = u;
				//set it's first child to the child
				u = T->child;
				//note how many matches have been made total.
			}
			else //otherwise break the edge that far down.
			{
				return sumI + read_ptr;
			}
		}
		//if not go to the next child
		else
		{
			u = u->sibling;
		}
	}
	//if there is no matching child, move on!
	//printf("Total matches this time: %i", sumI);
	return sumI + read_ptr;
}

void Sequence::parseFastaIntoReads(char * fastaFile)
{
	ifstream file(fastaFile);
	int i = 0;
	string line;

	while (file.good())
	{
		NUMREADS++;
		file >> READS[i].header;
		//file >> line;

		file >> READS[i].seq;
		//cout << READS[i].seq << endl;
		i++;
	}
	return;
}

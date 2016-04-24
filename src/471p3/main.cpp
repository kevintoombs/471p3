#include "includes.h"
#include "genomics.h"
#include "dynamicTable.h"
#include "suffixTree.h"

using namespace std;
int P2main(char* fasta, char* alphabet);

int main(int argc, char *argv[])
{
	//The below function shows the original functionality from program 1.
	DP_table::demo(argv[1], argv[3], argv[2]);

	//The two functions below are as follows
	//1: creates a dp table from arguments
	//2: demos that table like in program 1
	DP_table t(argv[1], argv[3], atoi(argv[2]));
	t.demoTable();

	//The below function shows the original functionality from program 2.
	McSuffixTree::demo(argv[4], argv[5]);

	//The two functions below are as follows
	//1: creates a st from arguments
	//2: demos that s2 like in program 2
	McSuffixTree mST(argv[4], argv[5]);
	mST.demoTree();

	cin.ignore();
	return 0;
	
}

int P2main(char* fasta, char* alphabet)
{
	string inputString = Sequence::parseFasta(fasta); //input argument one should be a fasta input file
	Alphabet a = Alphabet::parseAlphabet(alphabet, '$'); //input argument two should be an alphabet file

														//Tree is constructed from the given alphabet and string.
	unsigned int start = clock();
	McSuffixTree mST(inputString, a);
	cout << "The tree took: " << clock() - start << "ms to build" << endl;
	cout << endl;

	cout << "TREE STATS" << endl;
	cout << "Number of internal nodes: " << mST.nodes - mST.s.length() + 1 << endl; //cheating
	cout << "Number of leaves: " << mST.s.length() << endl; //cheating
	cout << "Number of total nodes: " << mST.nodes + 1 << endl;
	cout << "sizeof() tree in bytes: " << sizeof(mST) << endl;
	cout << endl;

	srand((int)time(NULL));
	cout << "PRINT SOME RANDOM CHILDREN" << endl;
	cout << "==Showing children of node:" << 0 << endl;
	mST.printNode(0);
	for (int i = 1; i < 5; i++)
	{
		int randN = rand() % mST.nodes;
		cout << "==Showing children of node:" << randN << endl;
		mST.printNode(randN);
	}
	cout << endl;

	mST.findLCS();
	cout << "longest matching repeat:\n" << mST.printString(mST.LCS) << endl;
	cout << endl;
	cout << "coords of longest matching repeat: " << mST.index1 << "," << mST.index2 << endl;
	cout << "string-depth of deepest eternal node: " << mST.LCSdepth << endl;
	float depthSum = (float)mST.depthSum;
	float internalNodes = (float)(mST.nodes - mST.s.length() + 1);
	float avgSDoIN = depthSum / internalNodes;
	cout << "average string depth of internal nodes:" << avgSDoIN << endl;
	cout << endl;


	if (mST.s.length() < 100000)
	{
		mST.printDFST();
		cout << "DFST written to out file." << endl;
		mST.BWT();
		cout << "BWT written to out file" << endl;
	}
	else
	{
		cout << "Input string is at least 100000 characters" << endl;
		cout << "Not outputting DFTS or BWT" << endl;
	}
	cout << endl;

	cout << "Press enter to exit." << endl;
	cin.ignore();
	return 0;


}
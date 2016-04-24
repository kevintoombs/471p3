#pragma once

#include "includes.h"
#include "genomics.h"


class Node {
public:
	int startIndex = 0; //where the string starts [1.....n]
	int stringSize = 0; //where the string ends when added to  startIndex [1.....n]
	int nodeDepth = 0;

	int nodeNumber = 0;
	int suffixID = 0;

	Node *parent = NULL;
	Node *sL = NULL;
	Node *child = NULL;
	Node *sibling = NULL;
	int start_leaf_index = -1; 
	int end_leaf_index = -1;


	Node();
	Node(int sI, int sD, int nD);
	Node(int sI, int sS, int nD, Node *pP, int nN);
};

class McSuffixTree
{
public:
	Node *root;
	Node *u;
	Node *LongestInternalNode;
	int ExactMatchLength;
	std::string s;
	Alphabet sigma;
	Printer printer;

	bool DEBUG = 0;
	int nodes = 0;

	Node* LCS;
	int LCSdepth;
	int index1;
	int index2;

	int depthSum = 0;
	int DFSTcounter;
	int suffix = 1;

	McSuffixTree(std::string sIn, Alphabet aIn);

	McSuffixTree(char * fasta, char * alphabet);

	void demoTree();

	void init();

	void displayAllChildren(Node *u);

	void displayChildren(Node *u);

	void printDFST();

	void DFSTHelper(Node *node);

	void findLCS();

	void findLCSHelper(Node *node);

	void findPath(Node *v, int t);

	std::string printString(Node* n);

	std::string printStringHelper(Node* n);
	std::string pathLabel(Node *u);

	void insertSuffix(int i);

	//(SL(u)) && (u != root)
	void case1a(Node *u, int i);

	//(SL(u)) && (u == root)
	void case1b(Node *u, int i);

	//!(SL(u)) && (u' != root)
	void case2a(Node *u, int i);

	//!(SL(u)) && (u' == root)
	void case2b(Node *u, int i);

	static int deep(Node* n);

	void nodeHops(Node* vPrime, int betaStart, int betaLength);

	//called to insert a string under a node by adding a new terminal node.
	void insertNode(Node *parent, int stringStart);

	void fixOrder(Node* parent);

	Node* edgeBreak(Node* v, Node* vChild, Node* vSibling, int correctComparisons);

	void increaseDepth(Node* n);

	void BWT();

	void BWTHelper(Node* node, std::vector<int>* bwtArray);

	void ExactMatchingRepeat();

	void ExactMatchHelper(int currentDepth, Node* node);

	void printNode(int n);

	void printNodeHelper(Node *node, int n);
	static void demo(char * fasta, char * alphabet);
};
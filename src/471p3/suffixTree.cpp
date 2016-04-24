#include "suffixTree.h"

using namespace std;

Node::Node() {};
Node::Node(int sI, int sD, int nD)
{
	startIndex = sI;
	stringSize = sD;
	nodeDepth = nD;
}
Node::Node(int sI, int sS, int nD, Node *pP, int nN)
{
	startIndex = sI;
	stringSize = sS;
	nodeDepth = nD;
	parent = pP;
	nodeNumber = nN;
}



McSuffixTree::McSuffixTree(string sIn, Alphabet aIn)
{
	sigma = aIn;
	s = sIn + sigma.rootDelimiter->sym;

	printer = Printer();

	init();
}
McSuffixTree::McSuffixTree(char* fasta, char* alphabet)
{
	string inputString = Sequence::parseFasta(fasta); //input argument one should be a fasta input file
	Alphabet a = Alphabet::parseAlphabet(alphabet, '$'); //input argument two should be an alphabet file
	s = inputString + a.rootDelimiter->sym;
	sigma = a;
	init();

														 //Tree is constructed from the given alphabet and string.
	unsigned int start = clock();
	McSuffixTree mST(inputString, a);
	//cout << "The tree took: " << clock() - start << "ms to build" << endl;
	//cout << endl;

	

}

void McSuffixTree::demoTree()
{
	cout << "TREE STATS" << endl;
	cout << "Number of internal nodes: " << this->nodes - this->s.length() + 1 << endl; //cheating
	cout << "Number of leaves: " << this->s.length() << endl; //cheating
	cout << "Number of total nodes: " << this->nodes + 1 << endl;
	cout << "sizeof() tree in bytes: " << sizeof(this) << endl;
	cout << endl;

	srand((int)time(NULL));
	cout << "PRINT SOME RANDOM CHILDREN" << endl;
	cout << "==Showing children of node:" << 0 << endl;
	this->printNode(0);
	for (int i = 1; i < 5; i++)
	{
		int randN = rand() % this->nodes;
		cout << "==Showing children of node:" << randN << endl;
		this->printNode(randN);
	}
	cout << endl;

	this->findLCS();
	cout << "longest matching repeat:\n" << this->printString(this->LCS) << endl;
	cout << endl;
	cout << "coords of longest matching repeat: " << this->index1 << "," << this->index2 << endl;
	cout << "string-depth of deepest eternal node: " << this->LCSdepth << endl;
	float depthSum = (float)this->depthSum;
	float internalNodes = (float)(this->nodes - this->s.length() + 1);
	float avgSDoIN = depthSum / internalNodes;
	cout << "average string depth of internal nodes:" << avgSDoIN << endl;
	cout << endl;


	if (this->s.length() < 100000)
	{
		this->printDFST();
		cout << "DFST written to out file." << endl;
		this->BWT();
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
}

void McSuffixTree::init()
{
	root = new Node(0, 0, 0);
	root->parent = root;
	root->sL = root;
	u = root->parent;

	for (unsigned int i = 1; i <= s.length(); i++)
	{
		insertSuffix(i);
		//displayAllChildren(root);
		__noop;
	}
}

void McSuffixTree::displayAllChildren(Node *u)
{
	Node *child = u->child;

	if (DEBUG == 1) cout << u->nodeDepth << ":";
	for (int i = 0; i < u->nodeDepth; i++)
	{
		if (DEBUG == 1) cout << " ";
	}
	if (DEBUG == 1) cout << "[";
	while (child != NULL)
	{
		cout << child->startIndex << "," << child->stringSize << "|" << deep(child) << endl;
		if (DEBUG == 1) cout << "(" << child->startIndex << "," << child->stringSize;
		if (child->nodeDepth != u->nodeDepth + 1)
		{
			if (DEBUG == 1) cout << "{!!!" << child->nodeDepth << ": not correct node depth}";
		}

		if (child->sL != NULL && false)
		{
			if (DEBUG == 1) cout << "|!!!" << child->sL->startIndex << "," << child->sL->stringSize << "suffix link|";
		}
		if (DEBUG == 1) cout << ")";
		child = child->sibling;
	}
	if (DEBUG == 1) cout << "]" << endl;

	child = u->child;
	while (child != NULL)
	{
		displayAllChildren(child);
		child = child->sibling;
	}

}

void McSuffixTree::displayChildren(Node *u)
{
	Node *child = u->child;
	while (child != NULL)
	{
		//cout << child->startIndex << "," << child->stringSize << "|" << deep(child) << endl;
		string s = printString(child);
		string l = pathLabel(child);
		int sLNN = -1;
		if (child->sL)
		{
			sLNN = child->sL->nodeNumber;
		}
		printf("Node:%7i   Index,Length:(%7i,%7i)   StringDepth:%7i   SL:%7i\n", child->nodeNumber, child->startIndex, child->stringSize, deep(child), sLNN);
		//cout << l << endl;
		child = child->sibling;
	}
}

void McSuffixTree::printDFST()
{
	cout << endl;
	DFSTHelper(root);
	cout << endl;
}

void McSuffixTree::DFSTHelper(Node *node)
{
	if (node == NULL)
	{
		return;
	}
	if (node->nodeNumber % 100 == 0)
	{
		cout << "\rTraversing through the tree.";
		printf(" On node (%i/%i).", (DFSTcounter++) * 100, nodes);
	}
	int d = deep(node);

	printer.printToFile(to_string(d), "DFST");
	printer.printToFile(" | ", "DFST");

	DFSTHelper(node->child);
	DFSTHelper(node->sibling);
}

void McSuffixTree::findLCS()
{
	LCS = root;
	LCSdepth = deep(root);
	findLCSHelper(root);
}

void McSuffixTree::findLCSHelper(Node *node)
{

	if (node == NULL)
	{
		return;
	}


	if (node->child != NULL)
	{
		int depth = deep(node);
		depthSum += depth;
		if (depth > LCSdepth)
		{
			if (DEBUG == 1) cout << "New LCS node: (" << node->startIndex << "," << node->stringSize << ")." << endl;
			LCSdepth = depth;

			index1 = node->child->startIndex - depth;
			index2 = node->child->sibling->startIndex - depth;
			LCS = node;
		}
	}
	findLCSHelper(node->child);
	findLCSHelper(node->sibling);
}

void McSuffixTree::findPath(Node *v, int t)
{
	//check the children of n
	Node *n = v->child;
	Node *lastN = NULL;


	int sumI = 0;
	//if n has any children or hasn't run out of possiblities
	while (n != NULL)
	{

		//if the first character of the child matches matches...
		if (s[n->startIndex - 1] == s[t - 1 + sumI])
		{
			//set matches to 1 
			int i = 1;
			//check for more matches
			while ((s[n->startIndex - 1 + i] == s[t - 1 + i + sumI]) && (i < n->stringSize))
			{
				i++;
				if (i == n->stringSize) break;
			}

			//if all characters matched
			if (i == n->stringSize)
			{
				//set this child to the new parent
				v = n;
				//set it's first child to the child
				n = v->child;
				lastN = NULL;
				//note how many matches have been made total.
				sumI += i;
			}
			//otherwise break the edge that far down.
			else
			{
				if (DEBUG == 1) printf("===!!!!===break called with i = %i\n", i);
				Node* temp = edgeBreak(v, n, lastN, i);
				insertNode(temp, t + sumI + i);

				return;
			}
		}
		//if not go to the next child
		else
		{
			lastN = n;
			n = n->sibling;
		}
	}

	//if there is no matching child, insert one!
	insertNode(v, t + sumI);
	return;
}

string McSuffixTree::printString(Node* n)
{
	return printStringHelper(n);
}

string McSuffixTree::printStringHelper(Node* n)
{
	if (n == root)
	{
		return "";
	}

	string s = printStringHelper(n->parent);
	return s + pathLabel(n);
}
string McSuffixTree::pathLabel(Node *u)
{
	if (u == root)
	{
		return "";
	}
	return s.substr(u->startIndex - 1, u->stringSize);
}

void McSuffixTree::insertSuffix(int i)
{
	//DEBUG
	//if (DEBUG == 0) cout << i << "/" << s.length() << endl;
	if (DEBUG == 1) cout << "Inserting by index: " << s.substr(i - 1) << endl;
	//END DEBUG
	if (u->sL && (u != root)) //case1a (SL(u)) && (u != root) 
	{
		if (DEBUG == 1) cout << "case 1a called for suffix " << i << "." << endl;
		case1a(u, i);
	}
	else if (u->sL && (u == root)) //case1b (SL(u)) && (u == root)
	{
		if (DEBUG == 1) cout << "case 1b called for suffix " << i << "." << endl;
		case1b(u, i);
	}
	else if (!(u->sL) && (u->parent != root)) //case2a !(SL(u)) && (u != root)
	{
		if (DEBUG == 1) cout << "case 2a called for suffix " << i << "." << endl;
		case2a(u, i);
	}
	else if (!(u->sL) && (u->parent == root)) //case2b !(SL(u)) && (u == root)
	{
		if (DEBUG == 1) cout << "case 2b called for suffix " << i << "." << endl;
		case2b(u, i);
	}

}

//(SL(u)) && (u != root)
void McSuffixTree::case1a(Node *u, int i)
{
	Node *v = u->sL;
	findPath(v, i + deep(u) - 1);
}

//(SL(u)) && (u == root)
void McSuffixTree::case1b(Node *u, int i)
{
	Node *v = u->sL;
	findPath(v, i);
}

//!(SL(u)) && (u' != root)
void McSuffixTree::case2a(Node *u, int i)
{
	int sizeOfAlphaAndC = deep(u);
	Node *uPrime = u->parent;
	Node *vPrime = uPrime->sL;
	nodeHops(vPrime, u->startIndex, u->stringSize);
	Node *v = u->sL;
	findPath(v, i + sizeOfAlphaAndC - 1);
}

//!(SL(u)) && (u' == root)
void McSuffixTree::case2b(Node *u, int i)
{
	int sizeOfBeta = u->stringSize;
	Node *uPrime = u->parent;
	Node *vPrime = uPrime->sL;
	nodeHops(vPrime, u->startIndex + 1, u->stringSize - 1);
	Node *v = u->sL;
	findPath(v, i + sizeOfBeta - 1);
}

int McSuffixTree::deep(Node* n)
{
	int i = 0;
	while (n->stringSize != 0)
	{
		i += n->stringSize;
		n = n->parent;
	}

	return i;
}

void McSuffixTree::nodeHops(Node* vPrime, int betaStart, int betaLength)
{
	if (betaLength == 0)
	{
		u->sL = vPrime;
		return;
	}
	Node* child = vPrime->child;
	Node* parent = child->parent;
	Node* prevSibling = NULL;
	bool b = true;
	//begin node hop, stopping when we are at the endpoint, or will insert on an edge
	while (/*betaLength >= (child->stringSize)*/b && betaLength > 0)
	{
		parent = child->parent;
		prevSibling = NULL;
		while (s[child->startIndex - 1] != s[betaStart - 1])
		{
			prevSibling = child;
			child = child->sibling;
		}


		if (betaLength < (child->stringSize))
		{
			b = false;
		}
		else
		{
			betaLength -= child->stringSize;
			betaStart += child->stringSize;
			//at this point we have the right child

			if (betaLength > 0) child = child->child;
		}
	}

	//mid-edge insertion
	if (betaLength > 0)
	{
		//NEED TO INSERT FIND PATH HERE WITH SPECIAL
		Node* v = edgeBreak(parent, child, prevSibling, betaLength);

		u->sL = v; //TODO
					//u->sL = child;
	}
	//endpoint insertion
	else
	{
		u->sL = child;
		//findPath(child, betaStart);
	}
}

//called to insert a string under a node by adding a new terminal node.
void McSuffixTree::insertNode(Node *parent, int stringStart)
{
	if (DEBUG == 1) cout << "Making node for index " << stringStart << ".\n";

	Node *child = parent->child;
	Node *prevChild = parent;

	//if parent has no children...
	if (child == NULL)
	{
		parent->child = new Node(stringStart, s.length() - stringStart + 1, parent->nodeDepth + 1, parent, ++nodes);
		u = parent;
		parent->child->suffixID = suffix;
		suffix++;
		//fixOrder(parent);
		return;
	}
	//otherwise loop through children
	else
	{

		while (1)
		{
			if (prevChild == parent && s[stringStart - 1] < s[child->startIndex - 1])
			{
				Node* n = new Node(stringStart, s.length() - stringStart + 1, parent->nodeDepth + 1, parent, ++nodes);
				parent->child = n;
				n->sibling = child;
				u = parent;
				n->suffixID = suffix;
				suffix++;
				return;
			}
			else if (s[stringStart - 1] < s[child->startIndex - 1])
			{
				Node* n = new Node(stringStart, s.length() - stringStart + 1, parent->nodeDepth + 1, parent, ++nodes);
				n->sibling = child;
				prevChild->sibling = n;
				u = parent;
				n->suffixID = suffix;
				suffix++;
				return;
			}
			else
			{
				prevChild = child;
				child = child->sibling;
			}

			if (child == NULL)
			{
				Node* n = new Node(stringStart, s.length() - stringStart + 1, parent->nodeDepth + 1, parent, ++nodes);
				prevChild->sibling = n;
				u = parent;
				n->suffixID = suffix;
				suffix++;
				return;
			}
		}
	}

	//make new node
	//find correct spot in children chain for it lexographically .. 
	//we have children 1...i,j...n
	//child j is the one we want to insert, child i is the previous.
	//when we insert child j, child j's sibling pointer points to child i, and child i's sibling pointer is what j's used to be 


	//this is always the final function call, now we set the next u and repeat.
	exit(1);
	return;
}

void McSuffixTree::fixOrder(Node* parent)
{
	Node* moveMeDown = parent->child;
	while (moveMeDown != NULL)
	{
		Node* next = moveMeDown->sibling;
		while (next != NULL)
		{
			if (s[moveMeDown->startIndex - 1] > s[next->startIndex - 1])
			{
				std::swap(moveMeDown, next);
				next = moveMeDown->sibling;
			}
			else
			{
				next = next->sibling;
			}
		}
		moveMeDown = moveMeDown->sibling;
	}
}

Node* McSuffixTree::edgeBreak(Node* v, Node* vChild, Node* vSibling, int correctComparisons)
{
	if (DEBUG == 1) cout << "===!!!!===break called\n";
	Node* Ui = new Node(vChild->startIndex, correctComparisons, vChild->nodeDepth, v, ++nodes);
	Ui->child = vChild;
	//if vChild is not v's first child
	if (vSibling != NULL)
	{
		vSibling->sibling = Ui;
	}

	//if vChild IS v's first child
	else
		v->child = Ui;
	vChild->parent = Ui;

	//handle the case where the child being replaced has a sibling pointer
	if (vChild->sibling != NULL)
	{
		Ui->sibling = vChild->sibling;
		vChild->sibling = NULL;
	}
	//replaced child's node depth increases by 1
	increaseDepth(vChild);
	vChild->stringSize -= correctComparisons;
	//set previous child's start index to the last character of Ui + 1

	vChild->startIndex = Ui->startIndex + Ui->stringSize;
	//if (u->sL == NULL) u->sL = Ui;
	return Ui;
}

void McSuffixTree::increaseDepth(Node* n)
{
	n->nodeDepth++;
	if (n->child != NULL) increaseDepth(n->child);
	if (n->sibling != NULL) increaseDepth(n->sibling);
}

void McSuffixTree::BWT()
{
	vector<int> bwtArray;
	BWTHelper(root, &bwtArray);
	map<int, int> lookup;
	vector<int> sortedArray = bwtArray;
	sort(sortedArray.begin(), sortedArray.end());
	vector<char> bwtResult;
	unsigned int i;
	for (i = 0; i < sortedArray.size(); i++)
	{
		lookup.insert(pair<int, int>(sortedArray[i], i + 1));
	}
	for (i = 0; i < sortedArray.size(); i++)
	{
		bwtArray[i] = lookup[bwtArray[i]];
	}
	if (DEBUG == 1) cout << "[";
	for (i = 0; i < bwtArray.size(); i++)
	{
		if (DEBUG == 1) cout << bwtArray[i] << ", ";
	}
	if (DEBUG == 1) cout << "]" << endl;
	for (i = 0; i < bwtArray.size(); i++)
	{
		if (bwtArray[i] == 1)
			bwtResult.push_back('$');
		else
			bwtResult.push_back(s[bwtArray[i] - 2]);
	}
	ofstream outfile;
	outfile.open("outfile.txt");
	for (i = 0; i < bwtResult.size(); i++)
		outfile << bwtResult[i] << endl;
	outfile.close();
}

void McSuffixTree::BWTHelper(Node* node, vector<int>* bwtArray)
{
	if (node == NULL)
		return;
	if (node->child == NULL)
	{
		bwtArray->push_back(node->nodeNumber);
		BWTHelper(node->sibling, bwtArray);
		return;
	}
	else
	{
		BWTHelper(node->child, bwtArray);
	}
	BWTHelper(node->sibling, bwtArray);
	return;
}

void McSuffixTree::ExactMatchingRepeat()
{
	//find longest internal node
	ExactMatchLength = 0;
	ExactMatchHelper(0, root);
	cout << "----------------------" << endl;
	cout << "Exact Matching Repeat:" << endl;
	cout << "Length: " << ExactMatchLength << endl;
	cout << "Starting Coordinates: " << LongestInternalNode->startIndex << " and " << abs(LongestInternalNode->child->stringSize - LongestInternalNode->child->sibling->stringSize) + LongestInternalNode->startIndex << endl;
	cout << "----------------------" << endl;
}

void McSuffixTree::ExactMatchHelper(int currentDepth, Node* node)
{
	if (node == NULL)
		return;
	if (node->child == NULL)
		return;
	currentDepth += node->stringSize;
	if (currentDepth > ExactMatchLength)
	{
		ExactMatchLength = currentDepth;
		LongestInternalNode = node;
	}
	ExactMatchHelper(currentDepth, node->child);
	ExactMatchHelper(currentDepth - node->stringSize, node->sibling);
	return;
}

void McSuffixTree::printNode(int n)
{
	printNodeHelper(root, n);
}

void McSuffixTree::printNodeHelper(Node *node, int n)
{
	if (node == NULL)
	{
		return;
	}
	if (node->nodeNumber == n)
	{
		displayChildren(node);
		return;
	}

	printNodeHelper(node->child, n);
	printNodeHelper(node->sibling, n);
}

void McSuffixTree::demo(char* fasta, char* alphabet)
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

}
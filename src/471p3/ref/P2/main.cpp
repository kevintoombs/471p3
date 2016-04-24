#include "suffixTree.cpp"

bool FORCEPRINT = 0;

int main(int argc, char *argv[])
{
	string inputString = Sequence::parseFasta(argv[1]); //input argument one should be a fasta input file
	Alphabet a = Alphabet::parseAlphabet(argv[2], '$'); //input argument two should be an alphabet file

	//Tree is constructed from the given alphabet and string.
	unsigned int start = clock();
	McSuffixTree mST(inputString, a);
	cout << "The tree took: " << clock() - start << "ms to build" << endl;
	cout << endl;

	cout << "TREE STATS" << endl;
	cout << "Number of internal nodes: " << mST.nodes - mST.s.length() << endl; //cheating
	cout << "Number of leaves: " << mST.s.length() << endl; //cheating
	cout << "Number of total nodes: " << mST.nodes << endl;
	cout << "sizeof() tree in bytes: " << sizeof(mST) << endl;
	cout << endl;

	srand(time(NULL));
	cout << "print some random children:" << "[[TODO]]" << endl;
	for (int i = 1; i < 5; i++)
	{
		int randN = rand() % inputString.length() + 1;

	}
	cout << endl;

	mST.findLCS();
	cout << "longest matching repeat: " << mST.printString(mST.LCS) << endl;
	
	cout << endl;
	cout << "coords of longest matching repeat: " << mST.index1 << "," << mST.index2 << endl;
	cout << "string-depth of deepest eternal node: " << mST.LCSdepth << endl;
	cout << "average string depth of internal nodes:" << "[[TODO]]" << endl;
	cout << endl;


	if (mST.s.length() < 100000 || FORCEPRINT)
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

int oldMain(int argc, char *argv[])
{

	string inputString = Sequence::parseFasta(argv[1]);
	//cout << inputString << endl;

	Alphabet a = Alphabet::parseAlphabet(argv[2], '$');
	//a.display();


	unsigned int start = clock();
	McSuffixTree mST(inputString, a);
	std::cout << "Time taken in millisecs: " << clock() - start << endl;

	//cout << mST.s << endl;
	//mST.sigma.display();
	//mST.BWT();


	mST.printDFST();

	//cout << mST.s << endl;
	//mST.sigma.display();

	//mST.ExactMatchingRepeat();

	cin.ignore();
	return 0;


}
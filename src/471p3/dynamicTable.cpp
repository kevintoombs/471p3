#include "dynamicTable.h"

using namespace std;



DP_table::DP_table()
{
}

DP_table::DP_table(std::string fasta, std::string config, int type)
{
	this->alightmentType = type;
	this->c = config::getConfig(config);
	parseFasta(fasta);

	buildTable();
	calcTable();
	//retrace();
}

bool DP_table::parseFasta(string fileName)
{
	ifstream fasta;
	if (true)
	{
		fasta.open(fileName, ios::in);
		if (!fasta.good())
		{
			return false;
		}

		char ch;
		int state = 0;
		int getId = 0;

		//probably overly complicated state based method to parse the fasta file
		while (fasta.get(ch))
		{

			//should always evaluate true at first 
			if (state == 0 && getId == 1)
			{
				if (isspace(ch))
				{
					//cout << "not alnum" << endl;
					getId = 0;
				}
				else if (!isspace(ch))
				{
					this->id1 += ch;
					//cout << this->id1 << endl;
				}
			}

			if (ch == '>' && state == 0 && getId == 0)
			{
				getId = 1;
			}



			//after the first line
			if (ch == '\n' && state == 0)
			{
				state = 1;
			}

			if (state == 2 && getId == 1)
			{
				if (isspace(ch))
				{
					//cout << "not alnum" << endl;
					getId = 0;
				}
				else if (!isspace(ch))
				{
					this->id2 += ch;
					//cout << this->id2 << endl;
				}
			}

			if (state == 1)
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					this->sequence1 += toupper(ch);
					//cout << this->sequence1 << endl;
					//cout << this->sequence1.length() << endl;
				}
				if (ch == '>')
				{
					state = 2;
					getId = 1;
				}
			}




			if (ch == '\n' && state == 2)
			{
				state = 3;
			}

			if (state == 3)
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					this->sequence2 += toupper(ch);
					//cout << this->sequence2 << endl;
					//cout << this->sequence2.length() << endl;
				}
				//cout << ch;
				if (ch == '>') {
					state = 4;
				}
			}


		}

		//cout << this->id1 << ',' << this->id2 << endl;
		fasta.close();
		return true;
	}

	else
	{
		cout << endl << "No input file specified. Correct usage is  $ <executable name> "
			<< "<input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>"
			<< endl;
		return false;
	}
}

void DP_table::setAlignmentType(char* arg)
{
	//argv[2]
	this->alightmentType = atoi(arg);
}

void DP_table::setAlignmentType(int type)
{
	this->alightmentType = type;
}

void DP_table::buildTable()
{
	//printf("Building Table.");
	DP_cell c;
	this->t.resize(this->sequence1.length() + 1, vector<DP_cell>(this->sequence2.length() + 1, c));
	//cout << endl;
	initTable();

	return;
}

void DP_table::initTable()
{
	int h = this->c.startGapScore;
	int g = this->c.continueGapScore;

	this->t[0][0].S = 0;
	this->t[0][0].D = -1147483648;
	this->t[0][0].I = -1147483648;

	for (size_t i = 1; i <= this->sequence1.length(); i++)
	{
		this->t[i][0].S = -1147483648;
		this->t[i][0].D = h + (int)i * g;
		this->t[i][0].I = -1147483648;
		this->t[i][0].dDir = 2;
	}
	for (size_t j = 1; j <= this->sequence2.length(); j++)
	{
		this->t[0][j].S = -1147483648;
		this->t[0][j].D = -1147483648;
		this->t[0][j].I = h + (int)j * g;
		this->t[0][j].iDir = 3;
	}
}

void DP_table::calcTable()
{

	int maxValue = 0;
	tuple<int, int> maxPair;

	//printf("Calculating Table[rows remaining]:");
	for (size_t i = 1; i <= this->sequence1.length(); i++)
	{
		//if (i % ((this->sequence1.length() / 100) + 1) == 0) printf("[%i]", (int)this->sequence1.length() - (int)i);
		for (size_t j = 1; j <= this->sequence2.length(); j++)
		{
			int sSub = subFunction(this->sequence1[i - 1], this->sequence2[j - 1], this->c);
			DP_cell subCell = this->t[i - 1][j - 1];
			DP_cell deleteCell = this->t[i - 1][j];
			DP_cell insertCell = this->t[i][j - 1];


			this->t[i][j].S = maximum(subCell.S + sSub, subCell.D + sSub, subCell.I + sSub, this->alightmentType);
			if (this->t[i][j].S == subCell.S + sSub) this->t[i][j].sDir = 1;
			else if (this->t[i][j].S == subCell.D + sSub) this->t[i][j].sDir = 2;
			else if (this->t[i][j].S == subCell.I + sSub) this->t[i][j].sDir = 3;
			this->t[i][j].D = maximum
				(deleteCell.S + this->c.startGapScore + this->c.continueGapScore,
					deleteCell.D + this->c.continueGapScore,
					deleteCell.I + this->c.startGapScore + this->c.continueGapScore,
					this->alightmentType);
			if (this->t[i][j].D == deleteCell.S + this->c.startGapScore + this->c.continueGapScore) this->t[i][j].dDir = 1;
			else if (this->t[i][j].D == deleteCell.D + this->c.continueGapScore) this->t[i][j].dDir = 2;
			else if (this->t[i][j].D == deleteCell.I + this->c.startGapScore + this->c.continueGapScore) this->t[i][j].dDir = 3; //the bug is not here though. (still never called)
			this->t[i][j].I = maximum
				(insertCell.S + this->c.startGapScore + this->c.continueGapScore,
					insertCell.D + this->c.startGapScore + this->c.continueGapScore, //FIXED//this is a bug lol, but it is like an impossible condition.
					insertCell.I + this->c.continueGapScore,
					this->alightmentType);
			if (this->t[i][j].I == insertCell.S + this->c.startGapScore + this->c.continueGapScore) this->t[i][j].iDir = 1;
			else if (this->t[i][j].I == insertCell.D + this->c.startGapScore + this->c.continueGapScore) this->t[i][j].iDir = 2; //FIXED//saaaame bug.
			else if (this->t[i][j].I == insertCell.I + this->c.continueGapScore) this->t[i][j].iDir = 3;
			// s.i = I + G is above, was I + H before. This bug was so hard to track down for a couple reasons. I'll detail what I think they are.
			// First. We were taught in class that you put the shorter string as your s2 so that your space complexity is not quadratic.
			//		My implementation takes almost 4GB's in debug after I implemented retrace (it was 2 before, like Ananth said).
			//		This means that usually the first input string is longer than the first. Because of that any insertions on 
			//		String 1 (if it is longer) will have to be matched by AT LEAST 1 MORE DELETION from string 2 in a global alignment.
			//		A shorter strings score still has to account for those empty, deleted, character... unless start gap penalties are nyah.
			// Second: It's the last of 3 in a shitty chain of horrible to read if else statements.
			// Third: my naming convention is just bad! I also should have kept ordering consitant. 
			// Fourth: It really is just a rare call for the alignment. I can't imagine many cases where we wouldn't just be better off doing a mismatch in the 
			//		first place. This is shown by the fact that it only caused the final number to be off by %10 in all of those calculations.
			if (this->alightmentType == 1)
			{
				int thisMax = this->t[i][j].cellMax();
				if (thisMax > maxValue)
				{
					maxValue = thisMax;
					maxPair = make_tuple(i, j);
					this->maxPair = maxPair;
				}
			}
		}
	}
	//cout << endl << endl;

}

int DP_cell::cellMax()
{
	int max = this->S;
	if (this->D > max)
	{
		max = this->D;
	}
	if (this->I > max)
	{
		max = this->I;
	}
	return max;
}

int DP_table::maximum(int S, int D, int I, int alignmentType)
{
	int max = S;
	if (D > max)
	{
		max = D;
	}
	if (I > max)
	{
		max = I;
	}
	if (alignmentType == 1 && max < 0) max = 0;
	return max;
}

int DP_table::subFunction(char a, char b, config c)
{
	if (a == b)
	{
		return c.matchScore;
	}
	else
	{
		return c.mismatchScore;
	}
}

void DP_table::printTable()
{
	int m = 0;
	int a = this->alightmentType;
	int s1length = (int)this->sequence1.length();
	int s2length = (int)this->sequence2.length();

	cout << " X";
	for (int j = 0; j <= /*20/**/s2length; j++) //do min(100, slength)
	{
		printf("%12d", j);
	}
	cout << endl;
	for (int i = 0; i <= /*20/**/s1length; i++) //do min(100, slength)
	{
		printf("%2d", i);
		for (int j = 0; j <= /*20/**/s2length; j++) //do min(100, slength)
		{
			DP_cell c = this->t[i][j];
			if (c.S >= -100) { printf("|%3d", c.S); }
			else printf("| - ");
			if (c.D >= -100) { printf(" %3d", c.D); }
			else printf("  - ");
			if (c.I >= -100) { printf(" %3d", c.I); }
			else printf("  - ");
		}
		cout << endl;
	}
}

void DP_table::retrace()
{
	//cout << "retrace\n";
	int matches = 0, mismatches = 0, gaps = 0, openingGaps = 0;
	int lastValue = -1;
	int lastDir = -1;
	int counter = 0;
	stack<char> s1, s2, r;
	int i = (int)this->sequence1.length();
	int j = (int)this->sequence2.length();

	if (this->alightmentType == 1)
	{
		i = get<0>(this->maxPair);
		j = get<1>(this->maxPair);
	}


	cout << "reverse dirs" << endl;
	DP_cell c = this->t[i][j];
	int dirSDI = 0;
	int max = c.cellMax();
	if (max == c.S) dirSDI = 1;
	if (max == c.D) dirSDI = 2;
	if (max == c.I) dirSDI = 3;
	int moveDir = c.cellMax2(max, dirSDI);
	while (i > 0 || j > 0)
	{
		if (this->alightmentType == 1 && this->t[i][j].cellMax() == 0)
		{
			DP_cell subCell = this->t[i - 1][j - 1];
			DP_cell deleteCell = this->t[i - 1][j];
			DP_cell insertCell = this->t[i][j - 1];
			break;
		}

		if (moveDir == 1) //s 
		{
			j--;
			i--;

			if (i >= 0 && j >= 0)
			{
				s1.push(this->sequence1[i]);
				s2.push(this->sequence2[j]);
				if (DP_table::subFunction(this->sequence1[i], this->sequence2[j], this->c) > 0)
				{
					matches++;
					r.push('|');
					//cout << "match" << endl;
				}
				else
				{
					mismatches++;
					r.push(' ');
					//cout << "mismatch" << endl;
				}
			}
		}
		else if (moveDir == 2) //d
		{
			i--;
			s2.push('-');
			r.push(' ');
			s1.push(this->sequence1[i]);
			if (lastDir == moveDir)
			{
				gaps++;
				//cout << "del" << endl;
			}
			else
			{
				openingGaps++;
				gaps++;
				//cout << "start gap" << endl;
			}
		}
		else if (moveDir == 3) //i
		{
			j--;
			s1.push('-');
			s2.push(this->sequence2[j]);
			r.push(' ');
			if (lastDir == moveDir)
			{
				//cout << "in" << endl;
				gaps++;
			}
			else
			{
				openingGaps++;
				gaps++;
				//cout << "start gap" << endl;
			}
		}
		lastDir = moveDir;

		//c is the cell we moved into based on where (S/D/I) the cell got its value from
		c = this->t[i][j];

		//find the appropriate value in this cell based in the dir arrow from the last S/D/I value
		int find = 0;
		if (dirSDI == 1) find = c.S;
		if (dirSDI == 2) find = c.D;
		if (dirSDI == 3) find = c.I;
		moveDir = c.cellMax2(find, dirSDI);

		printf("Moving to %i, for %i\n", moveDir, dirSDI);

	}

	cout << endl << endl;

	while (!s1.empty())
	{
		int to60 = 0;

		// printf("s1 %5d ", );
		while (to60 != 60)
		{
			if (!s1.empty()) cout << s1.top();
			if (!s1.empty()) s1.pop();
			to60++;
		}
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!r.empty()) cout << r.top();
			if (!r.empty()) r.pop();
			to60++;
		}
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!s2.empty()) cout << s2.top();
			if (!s2.empty()) s2.pop();
			to60++;
		}
		to60 = 0;
		cout << endl << endl;
	}
	cout << endl;

	c = this->t[i][j];

	cout << "Report:\n\n";
	if (this->alightmentType == 0)/**/ printf("Global optimal score = %i.\n", this->t[this->sequence1.size()][this->sequence2.size()].cellMax());
	if (this->alightmentType == 1)/**/ printf("Local optimal score = %i found at %i,%i.\n", this->t[get<0>(this->maxPair)][get<1>(this->maxPair)].cellMax(), get<0>(this->maxPair), get<1>(this->maxPair));
	printf("Number of:  matches = %i, mismatches = %i, gaps = %i, opening gaps = %i\n", matches, mismatches, gaps, openingGaps);
	float total = (float)gaps + (float)matches + (float)mismatches;
	float identities = (float)matches + (float)mismatches;
	float p1 = 100 * identities / total;
	float p2 = 100 * gaps / total;
	printf("Identities = %.0f/%.0f (%2.1f percent), Gaps = %i/%.0f (%2.1f percent)\n", identities, total, p1, gaps, total, p2);
	printf("Sanity Check: %i = ", matches * this->c.matchScore + mismatches * this->c.mismatchScore + gaps * this->c.continueGapScore + openingGaps * this->c.startGapScore);
	printf("%i * %i + %i * %i + %i * %i + %i * %i\n\n", matches, this->c.matchScore, mismatches, this->c.mismatchScore, gaps, this->c.continueGapScore, openingGaps, this->c.startGapScore);
	return;
}

int DP_table::direction(int i, int j)
{
	//1 = S, 2 = D, 3 = I 
	if (j == 0) return 3;
	if (i == 0) return 1;

	DP_cell c = this->t[i][j];
	DP_cell subCell = this->t[i - 1][j - 1];
	DP_cell deleteCell = this->t[i - 1][j];
	DP_cell insertCell = this->t[i][j - 1];

	int tv1 = subCell.cellMax();
	int tv2 = deleteCell.cellMax();
	int tv3 = insertCell.cellMax();


	DP_cell fromCell;
	int max = c.cellMax();

	int dir = 0;

	if (i != 0 && j != 0) //s 
	{
		int sSub = DP_table::subFunction(this->sequence1[i - 1], this->sequence2[j - 1], this->c);

		if (max == tv1 + sSub)
		{
			fromCell = subCell;
			dir = 1;
		}
	}
	if (i != 0) //d 
	{
		int testValue1 = deleteCell.S + this->c.startGapScore + this->c.continueGapScore;
		int testValue2 = deleteCell.D + this->c.continueGapScore;
		int testValue3 = deleteCell.I + this->c.startGapScore + this->c.continueGapScore;

		if (max == testValue1 || max == testValue2 /*|| max == testValue3*/)
		{
			fromCell = deleteCell;
			dir = 2;
		}
	}
	if (j != 0) //i
	{
		int testValue1 = insertCell.S + this->c.startGapScore + this->c.continueGapScore;
		int testValue2 = insertCell.D + this->c.startGapScore + this->c.continueGapScore;
		int testValue3 = insertCell.I + this->c.continueGapScore;

		if (max == testValue1 || max == testValue3 /*|| max == testValue2*/)
		{
			fromCell = insertCell;
			dir = 3;
		}
	}

	if (dir == 0)
	{
		exit(2);
	}

	//testDirection(max, c, this->c, dir, i, j, t);
	return dir;
}

//DONE: above
///NEED TO TRANSFER TO DP_TABLE:: below

int DP_table::demo(char* fastaFile, char* configFile, char* typeText)
{

	DP_table t;
	t.setAlignmentType(typeText);
	t.c = config::getConfig(configFile);
	if (!t.parseFasta(fastaFile))
	{
		cin.ignore();
		return 1;
	}

	printf("Scores: Match = %i, Mismatch = %i, h = %i, g = %i\n", t.c.matchScore, t.c.mismatchScore, t.c.startGapScore, t.c.continueGapScore);
	cout << endl;

	cout << "Sequence 1 = \"" << t.id1;
	printf("\", length = %i characters\n", (int)t.sequence1.length());
	cout << "Sequence 2 = \"" << t.id2;
	printf("\", length = %i characters\n", (int)t.sequence2.length());
	cout << endl;

	t.buildTable();
	t.calcTable();

	t.retrace();

	//printTable(t);
	cout << "Press enter to exit." << endl;
	cin.ignore();
	return 0;
}

void DP_table::demoTable()
{

	printf("Scores: Match = %i, Mismatch = %i, h = %i, g = %i\n", this->c.matchScore, this->c.mismatchScore, this->c.startGapScore, this->c.continueGapScore);
	cout << endl;

	cout << "Sequence 1 = \"" << this->id1;
	printf("\", length = %i characters\n", (int)this->sequence1.length());
	cout << "Sequence 2 = \"" << this->id2;
	printf("\", length = %i characters\n", (int)this->sequence2.length());
	cout << endl;


	this->retrace();

	//printTable(t);
	cout << "Press enter to exit." << endl;
	cin.ignore();
}


//a + b = c
void DP_table::testDirection(int lastValue, DP_cell to, int dir, int i, int j)
{
	//1 = S, 2 = D, 3 = I 
	bool y = true;
	char a = this->sequence1[i];
	char b = this->sequence2[j];

	if (dir == 2)
	{
		int sSub = DP_table::subFunction(this->sequence1[i - 1], this->sequence2[j - 1], this->c);
		y = to.cellMax() == lastValue;
	}
	else if (dir == 1)
	{
		y = true;

	}
	else if (dir == 3)
	{
		y = true;
	}

	cout << "[" << lastValue << " dir(" << 0 << ")" << to.cellMax() << "]";
	if (!y)
	{

		printf("\n@@@@X@@@@\n");
	}
}

int DP_cell::cellMax2(int find, int &mDir)
{
	int SDI = mDir;
	//printf("cout lol %i", "\n");
	if (find == this->S && SDI == 1)
	{
		cout << "(found " << find << ") ";
		mDir = this->sDir;
		return 1;
	}
	if (find == this->D && SDI == 2)
	{
		cout << "(found " << find << ") ";
		mDir = this->dDir;
		return 2;
	}
	if (find == this->I && SDI == 3)
	{
		cout << "(found " << find << ") ";
		mDir = this->iDir;
		return 3;
	}
	int r = 12;
	printf("cout lol %s", "\n");
	cout << "WHOOPS";
	cin.ignore();
	//exit(3);
	return -100;
}

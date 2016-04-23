#pragma once

#include "includes.h"
#include "genomics.h"

class DP_cell
{
public:
	int S = -1147483648;// numeric_limits<int>::min();
	int D = -1147483648;// numeric_limits<int>::min();
	int I = -1147483648;// numeric_limits<int>::min();
	int sDir = -1, dDir = -1, iDir = -1;
};

class DP_table
{
public:
	std::string sequence1, sequence2;
	std::string id1, id2;
	std::vector< std::vector<DP_cell> > t;
	config c;
	int alightmentType;
	std::tuple<int, int> maxPair;

	bool parseFasta(std::string filename);

};



int getAlignmentType(int argc, char *argv[]);
int demo(int argc, char * argv[]);
void buildTable(DP_table &t);
void initTable(DP_table &t);
void calcTable(DP_table &t);
int maximum(int S, int D, int I, int alignmentType);
int subFunction(char a, char b, config c);
void printTable(DP_table &t);
void retrace(DP_table &t);
void recursivelyPrintChildren(DP_table &t, int i, int j);
int direction(DP_table &t, int i, int j);
int maximum2(DP_cell &c, int &mDir);
int cellMax(DP_cell c);
void testDirection(int lastValue, DP_cell to, config c, int dir, int i, int j, DP_table &t);
int cellMax2(DP_cell &c, int find, int &mDir);
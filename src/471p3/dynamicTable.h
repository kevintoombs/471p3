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

	int cellMax();
	int cellMax2(int find, int &mDir);

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

	DP_table();
	DP_table(std::string fasta, std::string config, int type);

	bool parseFasta(std::string filename);
	void setAlignmentType(char* arg);
	void setAlignmentType(int type);
	void buildTable();
	void initTable();
	void calcTable();
	void printTable();
	void retrace();
	int direction(int i, int j);
	void testDirection(int lastValue, DP_cell to, int dir, int i, int j);
	void demoTable();

	static int maximum(int S, int D, int I, int alignmentType);
	static int subFunction(char a, char b, config c);
	static int demo(char* fastaFile, char* configFile, char* typeText);
	
};
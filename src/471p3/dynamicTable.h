#pragma once

#include "includes.h"
#include "genomics.h"

struct DP_cell
{
	int S = -1147483648;// numeric_limits<int>::min();
	int D = -1147483648;// numeric_limits<int>::min();
	int I = -1147483648;// numeric_limits<int>::min();
	int sDir = -1, dDir = -1, iDir = -1;
};

struct DP_table
{
	string sequence1, sequence2;
	string id1, id2;
	vector< vector<DP_cell> > t;
	config c;
	int alightmentType;
	tuple<int, int> maxPair;
};
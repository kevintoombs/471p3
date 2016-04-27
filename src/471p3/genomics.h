#pragma once

#include "includes.h"

class report
{
public:
	float PercentIdentity;
	float lengthCoverage;
	int start;
	int end;

	report::report(int matches, int alignlen, int l);
};

class config
{
public:
	int matchScore = 1;
	int mismatchScore = -2;
	int startGapScore = -5; //h
	int continueGapScore = -1; //g

	static config getConfig(int argc, char * argv[]);
	static config getConfig(std::string fileName);
};

class Alphabet
{
public:
	class Symbol
	{
	public:
		char sym;
		Symbol *next = NULL;
		Symbol::Symbol(char s);
	};

	Symbol *rootDelimiter;
	Symbol *currentSymbol;
	std::string aString;

	bool DEBUG = 0;

	Alphabet::Alphabet();
	Alphabet::Alphabet(char delimiter);
	void Alphabet::addSymbol(char newSymbol);
	void Alphabet::display();
	void Alphabet::display(Symbol *sym);
	static Alphabet Alphabet::parseAlphabet(char *alphabetFile, char delimiter);
};

class Sequence
{
public:
	std::string seq;
	std::string header;

	static std::string Sequence::parseFasta(char *fastaFile);
	static void Sequence::parseFastaIntoReads(char *fastaFile);
};

class Printer
{
public:
	int time = clock();
	int counter = 0;

	void Printer::printToFile(std::string stringToPrint, std::string relativeFilename);
	static void printP3(std::string stringToPrint, std::string relativeFilename);
};
#pragma once

#include "includes.h"

struct config
{
	int matchScore = 0;
	int mismatchScore = 0;
	int startGapScore = 0; //h
	int continueGapScore = 0; //g
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
};

class Printer
{
public:
	int time = clock();
	int counter = 0;

	void Printer::printToFile(std::string stringToPrint, std::string relativeFilename);
};
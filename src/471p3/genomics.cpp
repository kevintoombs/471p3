#include "genomics.h"

using namespace std;

Alphabet::Symbol::Symbol(char s)
{
	this->sym = s;
}

Alphabet::Alphabet()
	{

	}
Alphabet::Alphabet(char delimiter)
	{
		rootDelimiter = new Symbol(delimiter);
		currentSymbol = rootDelimiter;
	}
	void Alphabet::addSymbol(char newSymbol)
	{
		Symbol *sym = new Symbol(newSymbol);
		this->currentSymbol->next = sym;
		this->currentSymbol = sym;
	}
	void Alphabet::display()
	{
		display(rootDelimiter);
		cout << endl;
		cout << this->aString << endl;
	}
	void Alphabet::display(Symbol *sym)
	{
		cout << sym->sym << ",";
		if (sym->next)
			display(sym->next);
	}

	Alphabet Alphabet::parseAlphabet(char *alphabetFile, char delimiter)
	{
		Alphabet a = Alphabet(delimiter);

		ifstream fin;
		fin.open(alphabetFile);

		string sym, line;
		stringstream iss;
		getline(fin, line, '\n');
		iss << line;
		a.aString = delimiter;
		while (getline(iss, sym, ' '))
		{
			a.aString.append(sym);
			a.addSymbol(sym[0]);
		}

		fin.close();
		return a;
	}


	std::string Sequence::parseFasta(char *fastaFile)
	{
		string out;
		ifstream fin;
		fin.open(fastaFile);

		string seqHeader;
		getline(fin, seqHeader, '\n');

		while (fin.good())
		{
			string line;
			getline(fin, line, '\n');
			out += line;
		}

		fin.close();
		return out;
	}


	void Printer::printToFile(string stringToPrint, std::string relativeFilename)
	{
		relativeFilename += ".";
		relativeFilename += to_string(time);
		relativeFilename += ".out.txt";

		ofstream outFile;
		outFile.open(relativeFilename, ios::app);
		outFile << stringToPrint;

		counter++;
		if (counter == 20)
		{
			counter = 0;
			outFile << "\n";
		}

		outFile.close();
	}
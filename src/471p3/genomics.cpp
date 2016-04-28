#include "genomics.h"

using namespace std;

config config::getConfig(int argc, char *argv[])
{
	
	string filename;
	

	if (argc == 4)
	{
		cout << argv[3];
		filename = argv[3];
	}
	else
	{
		filename = "parameters.config";
		cout << "Default file: " << filename << " opened." << endl;
	}

	return getConfig(filename);
	
}

config config::getConfig(string fileName)
{
	ifstream conFile;
	config c;

	conFile.open(fileName, ios::in);
	if (!conFile.good())
	{
		cout << endl << fileName << " not found. All paremters default to 0." << endl;
		return c;
	}

	string line;
	while (getline(conFile, line))
	{
		stringstream s(line);
		string tmp1;
		string tmp2;
		while (!s.eof()) {
			s >> tmp1;
			s >> tmp2;
		} //cout << tmp1 << tmp2 << endl;
		if (tmp1 == "match")
		{
			c.matchScore = stoi(tmp2);
		}
		if (tmp1 == "mismatch")
		{
			c.mismatchScore = stoi(tmp2);
		}
		if (tmp1 == "h")
		{
			c.startGapScore = stoi(tmp2);
		}
		if (tmp1 == "g")
		{
			c.continueGapScore = stoi(tmp2);
		}
	}

	return c;
}

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

	report::report(int matches, int alignlen, int l)
	{
		this->PercentIdentity = (float)matches / (float)alignlen;
		this->lengthCoverage = (float)alignlen / (float)l;
	}

	void Printer::printP3(string stringToPrint, std::string relativeFilename)
	{
		std::string fName = "MappingResults_";
		relativeFilename = fName + "fasta";
		relativeFilename += ".txt";

		ofstream outFile;
		outFile.open(relativeFilename, ios::app);

		outFile << stringToPrint;
		outFile.close();
	}

	
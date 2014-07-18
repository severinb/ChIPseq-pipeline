#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>
using namespace std;
using namespace boost::iostreams;



/* COMPILE dynamic:
 * /import/bc2/soft/app/gcc/4.6.3/Linux/bin/g++ -I/import/bc2/soft/app/boost/1.42.0/Linux/include -L/import/bc2/soft/app/boost/1.42.0/Linux/lib -lboost_program_options -lboost_iostreams -o binReads binReads.cpp
 * RUN:
 * LD_LIBRARY_PATH=/import/bc2/soft/app/gcc/4.6.3/Linux/lib64/:/import/bc2/soft/app/boost/1.42.0/Linux/lib ./binReads
 *
 * COMPILE static:
 * /import/bc2/soft/app/gcc/4.6.3/Linux/bin/g++ -I/import/bc2/soft/app/boost/1.42.0/Linux/include -L/import/bc2/soft/app/boost/1.42.0/Linux/lib -o binReads binReads.cpp  -lboost_program_options -lboost_iostreams -lz -static
 */

class windowT
{
	int beg;
	int end;
	string chr;
};
class chrInfoT
{
	public:
	int len;		// length (the max beg position)
	int lastRead;	// keeps the head of the list of all reads at a given chromosome
	chrInfoT() : len(0), lastRead(0) {}
};
typedef boost::unordered_map<string, chrInfoT> string2chrInfoT;
namespace po = boost::program_options;
using namespace boost;
int main(int argc, char* argv[])
{
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
		    ("help", "produce help message")
		    ("input-file,i", po::value<string>(), "input file")
		;
		
		po::positional_options_description p;
		p.add("input-file", -1);
		
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
		po::notify(vm);
	
		if (vm.count("help") || ! vm.count("input-file"))
		{
		    cerr << desc << "\n";
		    return 1;
		}
		//if (vm.count("input-file"))
		//{
		//	cout << "Input file: " << vm["input-file"].as<string>() << "\n";
		//}
		ifstream file(vm["input-file"].as<string>().c_str(), ios_base::in | ios_base::binary);
		filtering_istream in;
		in.push(gzip_decompressor());
		in.push(file);
		int beg, end;
		string id;
		char str;
		double score;
		windowT currWindow;
		for(string chr; in >> chr >> beg >> end >> id >> score >> str; )
		{
			cout << chr << "\t" << str << '\n';
			if(beg < currWindow.end)
			{

			}
			else
			{

			}
		}
		return 0;
	}
	catch(std::exception& e)
    {
        cout << e.what() << "\n";
        return 1;
    }
}

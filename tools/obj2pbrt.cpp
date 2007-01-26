// obj2pbrt.cpp : Defines the entry point for the console application.
//

#include <list>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

struct tuple 
{
	float x;
	float y;
	float z;
};

int main(int argc, char *argv[])
{
	char* buffer = (char*)calloc(256, sizeof(char));

	if ( argc < 2 )
	{
		printf("filename not specified");
		exit(1);
	}

	string filename(argv[1]);

	// Create files names
	string infile = filename + ".obj";
	string outfile = filename + ".pbrt";

	ifstream input( infile.c_str() );
	ofstream output( outfile.c_str() );

	list<tuple> vertices;
	list<tuple> tris;

	tuple t;

	while ( ! input.eof() )
	{ 
		input >> buffer;

		// Skip comments
		if (buffer[0] == '#')
		{
			input.ignore(256, '\n');
			continue;
		}
		
		// Vertices
		if (buffer[0] == 'v')
		{

			// get the three vertices
			input >> t.x >> t.y >> t.z;
			vertices.push_back( t );
		}

		// Faces
		if (buffer[0] == 'f')
		{

			// get the three face indices
			input >> t.x >> t.y >> t.z;
			tris.push_back( t );
		}
	}


	/************************************************************************/
	/* WRITE TO PBRT FILE                                                   */
	/************************************************************************/
	output << "AttributeBegin \n Shape \"trianglemesh\" \n\"point P\" [\n";

	// Write vertices
	list<tuple>::iterator vertIt = vertices.begin();
	for (; vertIt != vertices.end(); vertIt++)
	{
		tuple vertex = *vertIt;
		output << " " << vertex.x << " " << vertex.y << " " << vertex.z << endl;
	}

	output << "]\n\n\t\"integer indices\" [";

	// Write faces
	list<tuple>::iterator trisIt = tris.begin();
	for (; trisIt != tris.end(); trisIt++)
	{
		tuple tri = *trisIt;
		output << (tri.x - 1) << " " << (tri.y - 1) << " " << (tri.z - 1) << endl;
	}

	output << "] \nAttributeEnd";

	return 0;
}


/* This file is part of Cloudy and is copyright (C)1978-2023 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <stdint.h>
#include "vectorhash.h"

using namespace std;

int main(int argc, char** argv)
{
	if( argc < 2 )
	{
		cout << "usage: " << argv[0] << " <file>..." << endl;
		return 1;
	}
	for( int i=1; i < argc; ++i )
	{
		FILE* io = fopen( argv[i], "r" );
		if( io == 0 )
		{
			cout << argv[0] << ": " << argv[i] << ": No such file or directory" << endl;
			return 1;
		}
		string vh128sum = VHstream( io );
		fclose( io );
		cout << vh128sum << "  " << argv[i] << endl;
	}
	return 0;
}

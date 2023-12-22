#include "cddefines.h"

void test(const multi_arr<long,3,ARPA_TYPE,false>& arr)
{
	multi_arr<double,3,ARPA_TYPE,false>::const_iterator p = arr.ptr(0,0,0);
}

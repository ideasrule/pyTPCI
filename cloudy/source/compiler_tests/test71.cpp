#include "cddefines.h"

void test(const flex_arr<long,false>& arr)
{
	flex_arr<double,false>::const_iterator p = arr.ptr(0);
}

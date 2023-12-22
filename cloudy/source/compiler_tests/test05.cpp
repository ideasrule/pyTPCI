#include "cddefines.h"

void test(multi_arr<long,3>& arr)
{
	ml3ci p = arr.ptr(0,0,0);
	*(--p) = 1;
}

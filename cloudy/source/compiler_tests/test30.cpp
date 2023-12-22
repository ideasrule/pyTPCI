#include "cddefines.h"

struct k
{
	long n;
};

void test(multi_arr<k,3>& arr)
{
	multi_arr<k,3>::const_iterator p = arr.ptr(0,0,0);
	(2+p)->n = 1;
}

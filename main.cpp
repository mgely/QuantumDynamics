#include "header.h"
using namespace std;

int main()
{
	simulation s;
	for (int i = 0; i < 20; ++i)
	{
		cout << i << endl;
		s.step();
	}
	s.step();
	s.print_to_file(NORM);
	return 0;
}
#include "Matrix.hpp"
#include "Complex.h"
#include <iostream>
#include <vector>

using namespace std;
int main()
{
	vector<int> v = {1,2,3,4};
	Matrix<int> m(2,2,v);

	vector<int> v1 = {5,6,7,8};
	Matrix<int> m1(2,2,v1);


	cout << m << endl;
	cout << m1 << endl;	
	cout << m + m1 << endl;
	
}

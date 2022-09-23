#include <iostream>

using namespace std;

int main() {
	// C-style null-terminated string
	char w[5] = {'H', 'i', '!', '\n', '\0'};
	// C++ std string
	string x = "Hello, World!\n";
	cout << w;
	cout << x;
	return 0;
}

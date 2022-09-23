#include <iostream>

int main() {
	// C-style null-terminated string
	char w[5] = {'H', 'i', '!', '\n', '\0'};
	// C++ std string
	std::string x = "Hello, World!\n";
	std::cout << w;
	std::cout << x;
	return 0;
}

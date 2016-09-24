#include "stage.hpp"

#include <iostream>

using std::cout;
using std::endl;
using std::string;

void stage::exec(string name, const std::vector<string> &steps)
{
    for (auto& step: steps) {
        if (_executable.find(step) != _executable.end()) {
            cout << "Execute: " << name << "/" << step << endl;
            _executable[step]();
        } else {
            cout << "Skip unknown step: " << name << "/" << step << endl;
        }
    }
}


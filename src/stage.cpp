#include "stage.hpp"

#include <chrono>
#include <iostream>

using std::cout;
using std::endl;
using std::string;

void stage::exec(string name, const std::vector<string> &steps)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    for (auto& step: steps) {
        if (_executable.find(step) != _executable.end()) {
            cout << "Execute: " << name << "/" << step << endl;
            start = std::chrono::system_clock::now();
            _executable[step]();
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "Time " << name << "/" << step << ": " << time.count()
                 << endl;
        } else {
            cout << "Skip unknown step: " << name << "/" << step << endl;
        }
    }
}

/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

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

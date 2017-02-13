#include "input.hpp"

#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

void input_queue::init()
{
    std::ifstream ifs(_conf.input_file);
    if (!ifs.good()) {
        cout << "Invalid input file" << endl;
        exit(1);
    }

    for (string file; std::getline(ifs, file);) {
        sm_chunk chunk;
        chunk.file = file;
        chunk.begin = -1;
        chunk.end = -1;
        chunk.kind = NORMAL_READ;
        std::size_t found = file.find("_T_");
        if (found != std::string::npos) {
            chunk.kind = CANCER_READ;
        }

        _queue.enqueue(chunk);
        len++;
    }
}

bool input_queue::try_dequeue(sm_chunk &chunk)
{
    return _queue.try_dequeue(chunk);
}

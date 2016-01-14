#include <getopt.h>
#include <iostream>
#include <chrono>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <kseq.h>
#include <google/sparse_hash_map>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

KSEQ_INIT(int, read);

typedef google::sparse_hash_map<string, bool> seen_table;
seen_table seen;

void load_fq(string file)
{
    int len;
    int nreads = 0;
    FILE* in = fopen(file.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(in));
    while ((len = kseq_read(seq)) >= 0) {
        nreads++;
        string id(seq->name.s);
        string s(seq->seq.s);
        string q(seq->qual.s);
        if (!seen[id]) {
            cout << "@" << id << endl << s << endl << "+" << endl << q << endl;
        }
        seen[id] = true;
    }
    kseq_destroy(seq);
    fclose(in);
}

void display_usage()
{
    cout << "Usage: sm-join-fq [OPTIONS] -i INPUT_PATH" << endl;
    cout << "Options:" << endl;
    cout << " -i, --input INPUT_PATH" << endl;
    cout << " -h, --help" << endl;
}

int main(int argc, char *argv[])
{
    string input_path;

    static const char *opts = "i:h";
    static const struct option opts_long[] = {
        { "input", required_argument, NULL, 'i' },
        { "help", no_argument, NULL, 'h' },
        { NULL, no_argument, NULL, 0 },
    };

    int opt = 0;
    int opt_index;
    while (opt != -1) {
        opt = getopt_long(argc, argv, opts, opts_long, &opt_index);
        switch (opt) {
            case 'i': input_path = string(optarg); break;
            case 'h':
                display_usage();
                return 0;
        }
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    seen.resize(2000000000);
    for (int i = 0; i < 8; i++) {
        start = std::chrono::system_clock::now();
        load_fq(input_path + "p" + std::to_string(i) + "/filter-tn.fastq");
        end = std::chrono::system_clock::now();
        time = end - start;
        cerr << "FASTQ load time: " << time.count() << endl;
    }
}

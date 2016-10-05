#include <stdio.h>

#include <chrono>
#include <fstream>
#include <iostream>

#include "kseq.h"

#include "db.hpp"
#include "common.hpp"
#include "merge.hpp"

using std::cout;
using std::endl;
using std::string;

KSEQ_INIT(int, read);

void load_i2p(rocksdb::DB* db, string path, string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = path + "/filter-" + set + "." + std::to_string(i) + ".i2p";
    cout << "Merge: " << file << endl;
    std::ifstream in(file);
    if (!in.good()) {
        cout << "Failed to open: " << file << endl;
        return;
    }

    int n = 0;
    string sid;
    sm_pos_bitmap p;

    while (in >> sid >> p.a[0] >> p.a[1] >> p.b[0] >> p.b[1]) {
        string serialized;
        encode_pos(serialized, p);
        db->Merge(rocksdb::WriteOptions(), sid, serialized);
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void load_k2i(rocksdb::DB* db, string path, string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = path + "/filter-" + set + "." + std::to_string(i) + ".k2i";
    cout << "Merge: " << file << endl;
    std::ifstream in(file);
    if (!in.good()) {
        cout << "Failed to open: " << file << endl;
        return;
    }

    int n = 0;
    string kmer;
    int len = 0;

    while (in >> kmer >> len) {
        std::stringstream s;
        for (int i = 0; i < len; i++) {
            string sid;
            in >> sid;
            s << sid << " ";
        }
        db->Merge(rocksdb::WriteOptions(), kmer, s.str());
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
}

void load_seq(rocksdb::DB* db, string path, string set, int i)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    string file = path + "/filter-" + set + "." + std::to_string(i) + ".fastq";
    cout << "Merge: " << file << endl;
    FILE* in = fopen(file.c_str(), "r");
    if (in == NULL) {
        cout << "Failed to open: " << file << " (" << errno << ")" << endl;
        return;
    }

    int n = 0;
    int len;
    kseq_t *seq = kseq_init(fileno(in));
    rocksdb::WriteBatch batch;
    while ((len = kseq_read(seq)) >= 0) {
        string id(seq->name.s);
        string s(seq->seq.s);
        batch.Put(id, s);
        if (n % 10000 == 0) {
            db->Write(rocksdb::WriteOptions(), &batch);
            batch.Clear();
        }
        if (n % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "M: " << i << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
        n++;
    }
    db->Write(rocksdb::WriteOptions(), &batch);
    kseq_destroy(seq);
    fclose(in);
}

#include "prune.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "util.hpp"

using std::cout;
using std::endl;
using std::string;

prune::prune(const sm_config &conf) : stage(conf)
{
    std::ifstream ifs(_conf.input_file);
    if (!ifs.good()) {
        cout << "Invalid input file" << endl;
        exit(1);
    }

    for (string line; std::getline(ifs, line);) {
        _input_queue.enqueue(line);
        _input_len++;
    }

    _all_size = _conf.all_size / _conf.num_partitions / _conf.num_storers;
    _allowed_size = _conf.allowed_size / _conf.num_partitions /
                    _conf.num_storers;

    _executable["run"] = std::bind(&prune::run, this);
}

void prune::run()
{
    if (_conf.num_loaders > MAX_LOADERS) {
        cout << "Number of loaders is larger than MAX_LOADERS" << endl;
        exit(1);
    }

    if (_conf.num_storers > MAX_STORERS) {
        cout << "Number of storers is larger than MAX_STORERS" << endl;
        exit(1);
    }

    // Initialize message queues
    for (int i = 0; i < _conf.num_storers; i++) {
        for (int j = 0; j < _conf.num_loaders; j++) {
            _queues[i][j] = new sm_prune_queue(PRUNE_QUEUE_LEN);
        }
    }

    std::vector<std::thread> loaders;
    for (int i = 0; i < _conf.num_loaders; i++)
        loaders.push_back(std::thread(&prune::load, this, i));
    cout << "Spawned " << loaders.size() << " prune loader threads" << endl;

    std::vector<std::thread> storers;
    for (int i = 0; i < _conf.num_storers; i++)
        storers.push_back(std::thread(&prune::add, this, i));
    cout << "Spawned " << storers.size() << " prune storer threads" << endl;

    for (auto& loader: loaders)
        loader.join();
    _done = true;
    for (auto& storer: storers)
        storer.join();
}

void prune::load(int lid)
{
    string file;
    while (_input_len > 0) {
        while (_input_queue.try_dequeue(file)) {
            load_file(lid, file);
            _input_len--;
        }
    }
}

void prune::load_file(int lid, string file)
{
    int len;
    int nreads = 0;
    sm_bulk_key bulks[MAX_STORERS];
    gzFile in = gzopen(file.c_str(), "rb");

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    kseq_t *seq = kseq_init(in);
    while ((len = kseq_read(seq)) >= 0) {
        nreads++;

        if (lq_count(seq->qual.s, seq->qual.l) > len/10)
            continue;

        int p = 0;
        int l = seq->seq.l;
        int n = seq->seq.l;
        char *ps;

        while ((ps = (char*) memchr(&seq->seq.s[p], 'N', l - p)) != NULL) {
            n = ps - &seq->seq.s[p];
            if (n > 0) {
                load_sub(lid, &seq->seq.s[p], n, bulks);
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            load_sub(lid, &seq->seq.s[p], n, bulks);
        }

        if (nreads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << lid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
    }

    for (int sid = 0; sid < _conf.num_storers; sid++) {
        while (!_queues[sid][lid]->write(bulks[sid])) {
            continue;
        }
        bulks[sid].num = 0;
    }

    kseq_destroy(seq);
    gzclose(in);
}

inline void prune::load_sub(int lid, const char* sub, int len,
                            sm_bulk_key* bulks)
{
    if (len < _conf.k)
        return;

    int stem_len = _conf.k - 2;
    char stem[stem_len + 1];
    for (int i = 0; i <= len - _conf.k; i++) {
        strncpy(stem, &sub[i + 1], stem_len);
        stem[stem_len] = '\0';

        uint64_t m = 0;
        memcpy(&m, stem, MAP_LEN);
        hash_5mer(m);

        if (map_l1[m] != _conf.pid)
            continue;
        int sid = map_l2[m];
        sm_key key = strtob4(stem);

        bulks[sid].array[bulks[sid].num] = key;
        bulks[sid].num++;

        if (bulks[sid].num == BULK_KEY_LEN) {
            while (!_queues[sid][lid]->write(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void prune::add(int sid)
{
    double fp = _conf.false_positive_rate;
    _all[sid] = new bf::basic_bloom_filter(fp, _all_size);
    _allowed[sid] = new bf::basic_bloom_filter(fp, _allowed_size);

    sm_bulk_key* pmsg;
    while (!_done) {
        for (int lid = 0; lid < _conf.num_loaders; lid++) {
            pmsg = _queues[sid][lid]->frontPtr();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    add_key(sid, pmsg->array[i]);
                }
                _queues[sid][lid]->popFront();
                pmsg = _queues[sid][lid]->frontPtr();
            }
        }
    }

    for (int lid = 0; lid < _conf.num_loaders; lid++) {
        pmsg = _queues[sid][lid]->frontPtr();
        while (pmsg) {
            for (int i = 0; i < pmsg->num; i++) {
                add_key(sid, pmsg->array[i]);
            }
            _queues[sid][lid]->popFront();
            pmsg = _queues[sid][lid]->frontPtr();
        }
    }

    delete _all[sid];
}

void prune::add_key(int sid, sm_key key)
{
    if (!_all[sid]->lookup(key)) {
        _all[sid]->add(key);
    } else if (!_allowed[sid]->lookup(key)) {
        _allowed[sid]->add(key);
    }
}

#include "count.hpp"

#include <errno.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <thread>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

count::count(const sm_config &conf) : stage(conf)
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

    // Initialize tables and message queues
    for (int i = 0; i < NUM_STORERS; i++) {
        for (int j = 0; j < MAX_LOADERS; j++) {
            _queues[i][j] = new sm_queue(QMSG_LEN);
        }
    }

    _executable["run"] = std::bind(&count::run, this);
    _executable["dump"] = std::bind(&count::dump, this);
    _executable["restore"] = std::bind(&count::restore, this);
    _executable["stats"] = std::bind(&count::stats, this);
}

void count::run()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    std::vector<std::thread> loaders;
    for (int i = 0; i < _conf.num_loaders; i++)
        loaders.push_back(std::thread(&count::load, this, i));
    cout << "Spawned " << loaders.size() << " loader threads" << endl;

    std::vector<std::thread> storers;
    for (int i = 0; i < _conf.num_storers; i++)
        storers.push_back(std::thread(&count::incr, this, i));
    cout << "Spawned " << storers.size() << " storer threads" << endl;

    for (auto& loader: loaders)
        loader.join();
    _done = true;
    for (auto& storer: storers)
        storer.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Count run time: " << time.count() << endl;
}

void count::load(int lid)
{
    string file;
    while (_input_len > 0) {
        while (_input_queue.try_dequeue(file)) {
            load_file(lid, file);
            _input_len--;
        }
    }
}

void count::load_file(int lid, string file)
{
    int len;
    int nreads = 0;
    sm_bulk bulks[NUM_STORERS];
    gzFile in = gzopen(file.c_str(), "rb");

    // Identify read kind from file name.
    sm_read_kind kind = NORMAL_READ;
    std::size_t found = file.find("_T_");
    if (found != std::string::npos) {
        kind = CANCER_READ;
    }

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
                load_sub(lid, &seq->seq.s[p], n, kind, bulks);
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            load_sub(lid, &seq->seq.s[p], n, kind, bulks);
        }

        if (nreads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "C: " << lid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
    }

    for (int sid = 0; sid < NUM_STORERS; sid++) {
        while (!_queues[sid][lid]->write(bulks[sid])) {
            continue;
        }
        bulks[sid].num = 0;
    }

    kseq_destroy(seq);
    gzclose(in);
}

inline void count::load_sub(int lid, const char* sub, int len,
                            sm_read_kind kind, sm_bulk* bulks)
{
    if (len < KMER_LEN)
        return;

    char imer[IMER_LEN + 1];
    for (int i = 0; i <= len - KMER_LEN; i++) {
        strncpy(imer, &sub[i + 1], IMER_LEN);
        imer[IMER_LEN] = '\0';

        uint64_t m = 0;
        memcpy(&m, imer, MAP_LEN);
        hash_5c_map(m);

        if (map_l1[m] != _conf.pid)
            continue;
        int sid = map_l2[m];
        sm_key key = strtob4(imer);

        sm_value_offset off;
        off.first = code[sub[i]] - '0';
        off.last = code[sub[i + KMER_LEN - 1]] - '0';
        off.kind = kind;

        bulks[sid].array[bulks[sid].num] = sm_msg(key, off);
        bulks[sid].num++;

        if (bulks[sid].num == BULK_LEN) {
            while (!_queues[sid][lid]->write(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void count::incr(int sid)
{
    _tables[sid] = new sm_table();
    _tables[sid]->resize(TABLE_LEN);
    sm_cache cache = sm_cache();
    cache.resize(CACHE_LEN);

    sm_bulk* pmsg;
    while (!_done) {
        for (int lid = 0; lid < _conf.num_loaders; lid++) {
            pmsg = _queues[sid][lid]->frontPtr();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    incr_key(sid, &cache, pmsg->array[i].first, pmsg->array[i].second);
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
                incr_key(sid, &cache, pmsg->array[i].first, pmsg->array[i].second);
            }
            _queues[sid][lid]->popFront();
            pmsg = _queues[sid][lid]->frontPtr();
        }
    }

    cout << "Cache " << sid << ": " << cache.size() << endl;
    cout << "Table " << sid << ": " << _tables[sid]->size() << endl;
}

inline void count::incr_key(int sid, sm_cache* cache, sm_key key,
                            sm_value_offset off)
{
    // Use sm_cache to hold keys with a single appearance; as soon as a key in
    // increased more than once, it is placed into sm_table. The steps are as
    // follows:
    //
    // - Find key in cache.
    //   - Key doesn't exist in cache: insert in cache.
    //   - Key exists in cache: find key in table.
    //     - Key doesn't exist in table: insert key and cache in table.
    //     - Key exists in table: update entry if there's no overflow.

    sm_cache::const_iterator cit = cache->find(key);
    if (cit == cache->end()) {
        uint8_t val = (off.first << 6) | (off.last << 4) | (off.kind << 2);
        cache->insert(std::pair<sm_key, uint8_t>(key, val));
        return;
    }

    sm_table::const_iterator it = _tables[sid]->find(key);
    if (it == _tables[sid]->end()) {
        uint8_t cache_value = cit->second;
        sm_value_offset coff;
        coff.first = cache_value >> 6;
        coff.last = (cache_value >> 4) & 0x03;
        coff.kind = (sm_read_kind) ((cache_value >> 2) & 0x03);
        sm_value val;
        val.v[coff.first][coff.last][coff.kind] = 1;
        val.v[off.first][off.last][off.kind]++;
        _tables[sid]->insert(std::pair<sm_key, sm_value>(key, val));
    } else {
        uint32_t inc = it->second.v[off.first][off.last][off.kind] + 1;
        uint16_t over = inc >> 16;
        uint16_t count = inc & 0x0000FFFF;
        if (over == 0)
            (*_tables[sid])[key].v[off.first][off.last][off.kind] = count;
    }
}

void count::dump()
{
    spawn("dump", std::bind(&count::dump_table, this, std::placeholders::_1),
          NUM_STORERS);
}

void count::dump_table(int sid)
{
    string file = string("table-") + std::to_string(_conf.pid) + string("-") +
                  std::to_string(sid) + string(".sht");

    char buf[PATH_MAX] = "";
    if (getcwd(buf, PATH_MAX) != NULL) {
        cout << "Serialize " << string(buf) << "/" << file << endl;
    }

    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(1);
    }

    if (!_tables[sid]->serialize(sm_table::NopointerSerializer(), fp)) {
        cout << "Failed to serialize table " << _conf.pid << "-" << sid << endl;
        exit(1);
    }

    fclose(fp);
}

void count::restore()
{
    spawn("restore", std::bind(&count::restore_table, this, std::placeholders::_1),
          NUM_STORERS);
}

void count::restore_table(int sid)
{
    _tables[sid] = new sm_table();
    _tables[sid]->resize(TABLE_LEN);
    string file = string("table-") + std::to_string(_conf.pid) + string("-") +
                  std::to_string(sid) + string(".sht");

    char buf[PATH_MAX] = "";
    if (getcwd(buf, PATH_MAX) != NULL) {
        cout << "Unserialize " << string(buf) << "/" << file << endl;
    }

    FILE* fp = fopen(file.c_str(), "r");
    if (fp == NULL) {
        cout << "Failed to open " << file << " (" << errno << ")" << endl;
        exit(0);
    }

    _tables[sid]->unserialize(sm_table::NopointerSerializer(), fp);
    fclose(fp);
}

void count::stats()
{
    std::map<uint64_t, uint64_t> hist_n, hist_t;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    start = std::chrono::system_clock::now();

    uint64_t total_stems = 0;
    uint64_t total_kmers = 0;
    uint64_t total_sum = 0;

    for (int i = 0; i < _conf.num_storers; i++) {
        uint64_t num_stems = _tables[i]->size();
        uint64_t num_kmers = 0;
        uint64_t sum = 0;
        for (sm_table::const_iterator it = _tables[i]->begin();
             it != _tables[i]->end(); ++it) {
            uint64_t sum_t = 0;
            uint64_t sum_n = 0;
            for (int f = 0; f < 4; f++) {
                for (int l = 0; l < 4; l++) {
                    uint16_t nc = it->second.v[f][l][NORMAL_READ];
                    uint16_t tc = it->second.v[f][l][CANCER_READ];
                    sum_n += nc;
                    sum_t += tc;
                    num_kmers += (nc + tc > 0) ? 1 : 0;
                }
            }
            if (sum_n > 0)
                hist_n[int(log2(sum_n))]++;
            if (sum_t > 0)
                hist_t[int(log2(sum_t))]++;
            sum += sum_n + sum_t;
        }
        cout << "Table " << i << ": " << num_stems << " " << num_kmers
             << " " << sum << endl;
        total_stems += num_stems;
        total_kmers += num_kmers;
        total_sum += sum;
    }

    for (auto& h: hist_n) {
        cout << "Histo N: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    for (auto& h: hist_t) {
        cout << "Histo T: " << h.first << " " << exp2(h.first)
             << " " << h.second << endl;
    }

    cout << "Number of stems: " << total_stems << endl;
    cout << "Number of kmers: " << total_kmers << endl;
    cout << "Sum of counters: " << total_sum << endl;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Stats time: " << time.count() << endl;
}

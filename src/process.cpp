#include <common.hpp>
#include <process.hpp>

#include <errno.h>
#include <string>
#include <iostream>
#include <kseq.h>
#include <boost/algorithm/string.hpp>
#include <google/sparse_hash_map>

using std::cout;
using std::endl;
using std::string;

KSEQ_INIT(gzFile, gzread);

void process_load(int pid, int lid)
{
    string file;
    while (input_count > 0) {
        while (input_queue.try_dequeue(file)) {
            process_load_file(pid, lid, file);
            input_count--;
        }
    }
}

void process_load_file(int pid, int lid, string file)
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
                process_load_sub(pid, lid, &seq->seq.s[p], n, kind, bulks);
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            process_load_sub(pid, lid, &seq->seq.s[p], n, kind, bulks);
        }

        if (nreads % 100000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << lid << " " << time.count() << endl;
            start = std::chrono::system_clock::now();
        }
    }

    for (int sid = 0; sid < NUM_STORERS; sid++) {
        while (!queues[sid][lid]->write(bulks[sid])) {
            continue;
        }
        bulks[sid].num = 0;
    }

    kseq_destroy(seq);
    gzclose(in);
}

inline void process_load_sub(int pid, int lid, const char* sub, int len,
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

        if (map_l1[m] != pid)
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
            while (!queues[sid][lid]->write(bulks[sid])) {
                continue;
            }
            bulks[sid].num = 0;
        }
    }
}

void process_incr(int pid, int sid, int num_loaders)
{
    sm_table table = sm_table();
    table.resize(TABLE_LEN);
    sm_cache cache = sm_cache();
    cache.resize(CACHE_LEN);

    sm_bulk* pmsg;
    while (!process_done) {
        for (int lid = 0; lid < num_loaders; lid++) {
            pmsg = queues[sid][lid]->frontPtr();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    process_incr_key(&table, &cache, pmsg->array[i].first, pmsg->array[i].second);
                }
                queues[sid][lid]->popFront();
                pmsg = queues[sid][lid]->frontPtr();
            }
        }
    }

    for (int lid = 0; lid < num_loaders; lid++) {
        pmsg = queues[sid][lid]->frontPtr();
        while (pmsg) {
            for (int i = 0; i < pmsg->num; i++) {
                process_incr_key(&table, &cache, pmsg->array[i].first, pmsg->array[i].second);
            }
            queues[sid][lid]->popFront();
            pmsg = queues[sid][lid]->frontPtr();
        }
    }

    cout << "Cache " << sid << ": " << cache.size() << endl;

    std::this_thread::sleep_for(std::chrono::seconds(sid * 7));
    string file = string("table-") + std::to_string(pid) + string("-") +
                  std::to_string(sid) + string(".data");

    char buf[PATH_MAX] = "";
    if (getcwd(buf, PATH_MAX) != NULL) {
        cout << "Serialize " << string(buf) << "/" << file << endl;
    }

    FILE* fp = fopen(file.c_str(), "w");
    if (fp == NULL) {
        cout << "Failed to open: " << file << " (" << errno << ")" << endl;
        exit(1);
    }

    if (!table.serialize(sm_table::NopointerSerializer(), fp)) {
        cout << "Failed to serialize table " << pid << "-" << sid << endl;
        exit(1);
    }

    fclose(fp);
}

inline void process_incr_key(sm_table* table, sm_cache* cache, sm_key key, sm_value_offset off)
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

    sm_table::const_iterator it = table->find(key);
    if (it == table->end()) {
        uint8_t cache_value = cit->second;
        sm_value_offset coff;
        coff.first = cache_value >> 6;
        coff.last = (cache_value >> 4) & 0x03;
        coff.kind = (sm_read_kind) ((cache_value >> 2) & 0x03);
        sm_value val;
        val.v[coff.first][coff.last][coff.kind] = 1;
        val.v[off.first][off.last][off.kind]++;
        table->insert(std::pair<sm_key, sm_value>(key, val));
    } else {
        uint32_t inc = it->second.v[off.first][off.last][off.kind] + 1;
        uint16_t over = inc >> 16;
        uint16_t count = inc & 0x0000FFFF;
        if (over == 0)
            (*table)[key].v[off.first][off.last][off.kind] = count;
    }
}

#include <sm_common.hpp>
#include <sm_process.hpp>

#include <string>
#include <iostream>
#include <kseq.h>
#include <boost/algorithm/string.hpp>
#include <simdcomp.h>

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

        if (lq_count(seq->qual.s) > 8) {
            continue;
        }

        int p = 0;
        int n = 80;
        char *ps;

        while ((ps = (char*) memchr(&seq->seq.s[p], 'N', 80 - p)) != NULL) {
            n = ps - &seq->seq.s[p];
            if (n > 0) {
                process_load_sub(pid, lid, &seq->seq.s[p], n, kind, bulks);
                p += n;
            }
            p++;
        }

        n = 80 - p;
        if (n > 0) {
            process_load_sub(pid, lid, &seq->seq.s[p], n, kind, bulks);
        }

        if (nreads % 50000 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << nreads << " reads (" << time.count() << ")" << endl;
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

        uint32_t m = 0;
        memcpy(&m, imer, MAP_LEN);
        hash_4c_map(m);

        if (map_l1[m] != pid)
            continue;
        int sid = map_l2[m];
        sm_key key = strtob4(imer);

        sm_value_offset off;
        off.first = code[sub[i]] - '0';
        off.last = code[sub[i + KMER_LEN - 1]] - '0';
        off.kind = kind;

        // cout << "P: " << &imer[0] << " " << key << endl;

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

void process_incr(int sid, int num_loaders)
{
    sm_bulk* pmsg;
    while (!process_done) {
        for (int lid = 0; lid < num_loaders; lid++) {
            pmsg = queues[sid][lid]->frontPtr();
            while (pmsg) {
                for (int i = 0; i < pmsg->num; i++) {
                    process_incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
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
                process_incr_key(sid, pmsg->array[i].first, pmsg->array[i].second);
            }
            queues[sid][lid]->popFront();
            pmsg = queues[sid][lid]->frontPtr();
        }
    }
}

inline void process_incr_key(int sid, sm_key key, sm_value_offset off)
{
    sm_tally tally;
    sm_value *value;
    uint32_t b = 0;

    sm_table::iterator it = tables[sid].find(key);
    if (it != tables[sid].end()) {
        value = &it->second;
        simdunpack_length((const __m128i *) value->v, 32, (uint32_t *) &tally.v, value->b);
    }

    tally.v[off.first][off.last][off.kind]++;
    b = maxbits_length((const uint32_t *) &tally.v, 32);
    simdpack_length((const uint32_t *) &tally.v, 32, (__m128i *) &tables[sid][key].v, b);
    tables[sid][key].b = b;
}

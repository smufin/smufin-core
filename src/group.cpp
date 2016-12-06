#include "group.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <rocksdb/db.h>

#include "db.hpp"

using std::cout;
using std::endl;
using std::string;

group::group(const sm_config &conf) : stage(conf)
{
    _executable["run"] = std::bind(&group::run, this);
}

void group::run()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    for (int i = 0; i < MAX_GROUPERS; i++) {
        _l2p[i] = new l2p_table(100000);
        _l2k[i] = new l2k_table(100000);
        _l2i[i] = new l2i_table(100000);
        _l2r[i] = new l2r_table(100000);
    }

    for (int i = 0; i < 2; i++) {
        _i2r[i] = new i2r_table();
        _k2i[i] = new kmer_table();
    }

    _i2r[0]->resize(420000000);
    _i2r[1]->resize(500000000);
    _k2i[0]->resize(30000000);
    _k2i[1]->resize(60000000);

    string dir;
    rocksdb::Status status;

    rocksdb::DB* i2p_db;
    rocksdb::Options options;
    options.WAL_ttl_seconds = 0;
    options.WAL_size_limit_MB = 0;
    options.merge_operator.reset(new PositionsMapOperator());
    dir = _conf.output_path + "/i2p-" + sm::sets[TM] + ".rdb";
    cout << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options, dir, &i2p_db);
    assert(status.ok());

    rocksdb::DB* i2r_tm_db;
    rocksdb::Options options2;
    options2.WAL_ttl_seconds = 0;
    options2.WAL_size_limit_MB = 0;
    dir = _conf.output_path + "/seq-" + sm::sets[TM] + ".rdb";
    cout << "Open RocksDB: " << dir << endl;
    status = rocksdb::DB::OpenForReadOnly(options2, dir, &i2r_tm_db);
    assert(status.ok());

    string sid;
    sm_pos_bitmap p;

    int num_all = 0;
    int num_exist = 0;
    int num_map = 0;
    int num_match_a = 0;
    int num_match_b = 0;

    std::chrono::time_point<std::chrono::system_clock> istart, iend;
    std::chrono::duration<double> itime;
    istart = std::chrono::system_clock::now();

    // 1. Iterate through candidates.

    rocksdb::Iterator* it = i2p_db->NewIterator(rocksdb::ReadOptions());
    for (it->SeekToFirst(); it->Valid(); it->Next()) {
        num_all++;

        sid = it->key().ToString();
        p = decode_pos(it->value().data());

        std::vector<int> a_pos;
        std::vector<int> b_pos;

        string read;
        rocksdb::Status status;
        status = i2r_tm_db->Get(rocksdb::ReadOptions(), sid, &read);
        if (!status.ok())
            continue;

        if (num_all % 100000 == 0) {
            iend = std::chrono::system_clock::now();
            itime = iend - istart;
            cout << "S: " << num_all << " " << itime.count() << endl;
            istart = std::chrono::system_clock::now();
        }

        num_exist++;

        int read_length = read.size();

        if (read_length == 0)
            continue;

        if (read.find("N") != std::string::npos)
            continue;

        uint64_t m = 0;
        string sub = read.substr(0, MAP_LEN);
        memcpy(&m, sub.c_str(), MAP_LEN);
        hash_5c_map(m);
        if (map_l1[m] != _conf.pid)
            continue;
        int gid = map_l2[m];

        num_map++;

        get_positions_a(p.a, &a_pos);
        get_positions_b(p.b, &b_pos, read_length);
        std::reverse(b_pos.begin(), b_pos.end());

        int a_len = a_pos.size();
        int b_len = b_pos.size();

        if (a_len >= KMIN && a_len <= KMAX && match_window(a_pos)) {
            select_candidate(gid, sid, read, read, a_pos, 0);
            num_match_a++;
        }

        if (b_len >= KMIN && b_len <= KMAX && match_window(b_pos)) {
            char buf[RMAX + 1];
            copy(read.begin(), read.end(), buf);
            buf[read_length] = '\0';
            revcomp(buf, read_length);
            string directed_read = string(buf);
            select_candidate(gid, sid, read, directed_read, b_pos, 1);
            num_match_b++;
        }
    }

    delete it;
    delete i2p_db;
    delete i2r_tm_db;

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Candidate selection time: " << time.count() << endl;

    cout << "Number of iterated I2P: " << num_all << endl;
    cout << "Number of existing I2P: " << num_exist << endl;
    cout << "Number of mapped I2P: " << num_map << endl;
    cout << "Number of A-matching I2P: " << num_match_a << endl;
    cout << "Number of B-matching I2P: " << num_match_b << endl;
    for (int i = 0; i < MAX_GROUPERS; i++) {
        cout << "Number of candidates (" << std::to_string(i) << "): "
             << _l2k[i]->size() << endl;
    }

    cout << "Number of kmers (0): " << _k2i[0]->size() << endl;
    cout << "Number of kmers (1): " << _k2i[1]->size() << endl;

    // 2. Iterate K2I retrieving all k-mers seen in candidate positions.

    for (int i = 0; i < 2; i++) {
        start = std::chrono::system_clock::now();

        rocksdb::DB* k2i_db;
        rocksdb::Options options;
        options.WAL_ttl_seconds = 0;
        options.WAL_size_limit_MB = 0;
        options.merge_operator.reset(new IDListOperator());
        dir = _conf.output_path + "/k2i-" + sm::sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(options, dir, &k2i_db);
        assert(status.ok());

        int num_seen = 0;
        int num_kmer = 0;
        istart = std::chrono::system_clock::now();
        it = k2i_db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_kmer++;

            string kmer = it->key().ToString();
            string list = it->value().ToString();
            kmer_table::const_iterator kit = _k2i[i]->find(kmer);
            if (kit != _k2i[i]->end()) {
                (*_k2i[i])[kmer] = list;
                num_seen++;
            }

            if (num_kmer % 1000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "K: " << i << " " << num_kmer << " " << num_seen << " "
                     << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        cout << "Number of kmers seen (" << sm::sets[i] << "): "
             << num_seen << endl;
        cout << "Number of IDs seen (" << sm::sets[i] << "): "
             << _i2r[i]->size() << endl;


        delete it;
        delete k2i_db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Kmer iteration time (" << sm::sets[i] << "): "
             << time.count() << endl;
    }

    // 3. Iterate & collect reads.

    for (int i = 0; i < 2; i++) {
        start = std::chrono::system_clock::now();

        rocksdb::DB* i2r_db;
        dir = _conf.output_path + "/seq-" + sm::sets[i] + ".rdb";
        cout << "Open RocksDB: " << dir << endl;
        status = rocksdb::DB::OpenForReadOnly(rocksdb::Options(), dir, &i2r_db);
        assert(status.ok());

        int num_seen = 0;
        int num_read = 0;
        istart = std::chrono::system_clock::now();
        it = i2r_db->NewIterator(rocksdb::ReadOptions());
        for (it->SeekToFirst(); it->Valid(); it->Next()) {
            num_read++;

            string sid = it->key().ToString();
            string seq = it->value().ToString();
            sm_read read;
            encode_read(seq, read);
            (*_i2r[i])[sid] = read;

            if (num_read % 10000000 == 0) {
                iend = std::chrono::system_clock::now();
                itime = iend - istart;
                cout << "R: " << i << " " << num_read << " " << num_seen << " "
                     << itime.count() << endl;
                istart = std::chrono::system_clock::now();
            }
        }

        delete it;
        delete i2r_db;

        end = std::chrono::system_clock::now();
        time = end - start;
        cout << "Read iteration time (" << sm::sets[i] << "): "
             << time.count() << endl;
    }

    // 4. Populate candidate groups.

    start = std::chrono::system_clock::now();

    std::vector<std::thread> populators;
    for (int i = 0; i < _conf.num_groupers; i++)
        populators.push_back(std::thread(&group::populate, this, i));
    cout << "Spawned " << populators.size() << " populator threads" << endl;
    for (auto& populator: populators)
        populator.join();

    end = std::chrono::system_clock::now();
    time = end - start;
    cout << "Populate candidates time: " << time.count() << endl;
}

void group::encode_read(std::string& str, sm_read& read)
{
    read.len = str.size();
    for (int i = 0; i < read.len; i += 32) {
        read.seq[i/32] = strtob4(str.substr(i, 32).c_str());
    }
}

void group::decode_read(sm_read& read, std::string& str)
{
    char decode[4] = {'A', 'C', 'G', 'T'};
    str = "";

    int len_bits = read.len * 2;
    int pos = len_bits / 64;
    int shift = 64 - (len_bits % 64);
    read.seq[pos] = read.seq[pos] << shift;

    int l = 0;
    int b = 0;
    while (l < read.len) {
        int pos = b / 64;
        int shift = 62 - (b % 64);
        int i = ((read.seq[pos] >> shift) & 3);
        char c = decode[i];
        str += c;
        l++;
        b += 2;
    }
}

void group::get_positions_a(uint64_t bitmap[2], std::vector<int> *pos)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(p + offset);
        }
    }
}

void group::get_positions_b(uint64_t bitmap[2], std::vector<int> *pos, int len)
{
    for (int i = 0; i < 2; i++) {
        unsigned long tmp = bitmap[i];
        int offset = i * 64;
        while (tmp > 0) {
            int p = __builtin_ffsl(tmp) - 1;
            tmp &= (tmp - 1);
            pos->push_back(len - _conf.k - (p + offset));
        }
    }
}

bool group::match_window(std::vector<int> pos)
{
    if (pos.size() < _conf.window_min) {
        return false;
    }

    for (int i = 0; i <= pos.size() - _conf.window_min; i++) {
        std::vector<int> sub(pos.begin() + i, pos.begin() + i + _conf.window_min);
        if (sub.back() - sub.front() < _conf.window_len) {
            return true;
        }
    }

    return false;
}

void group::select_candidate(int gid, string& sid, string& seq, string& dseq,
                             std::vector<int>& pos, int dir)
{
    std::vector<string> kmers;
    for (int p: pos) {
        string kmer = dseq.substr(p, _conf.k);
        kmers.push_back(kmer);
        kmer_table::const_iterator it = _k2i[0]->find(kmer);
        if (it == _k2i[0]->end()) {
            (*_k2i[0])[kmer] = string();
            (*_k2i[1])[kmer] = string();
        }
    }
    (*_l2r[gid])[sid] = seq;
    (*_l2p[gid])[sid][dir] = pos;
    (*_l2k[gid])[sid][dir] = kmers;
}

void group::populate(int gid)
{
    const char comp_code[] = "ab";
    const char kind_code[] = "nt";

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;
    start = std::chrono::system_clock::now();

    std::ofstream ofs;
    string file = _conf.output_path + "/group." + std::to_string(_conf.pid) +
                  "-" + std::to_string(gid) + ".json";
    ofs.open(file);
    ofs << "{";

    int num_leads = 0;
    bool first_group = true;
    for (l2k_table::const_iterator it = _l2k[gid]->begin(); it != _l2k[gid]->end(); ++it) {
        num_leads++;
        string lid = it->first;
        index_count keep;
        index_count drop;

        for (int i = 0; i < 2; i++) {
            for (string kmer: it->second[i]) {
                keep[i][kmer] = 0;
                drop[i][kmer] = 0;
            }
        }

        populate_index(gid, lid, it->second[0], NN, keep, drop);
        populate_index(gid, lid, it->second[0], TN, keep, drop);
        populate_index(gid, lid, it->second[1], NN, keep, drop);
        populate_index(gid, lid, it->second[1], TN, keep, drop);

        l2r_table::const_iterator lit = _l2r[gid]->find(lid);
        if (lit == _l2r[gid]->end())
            continue;
        string seq = lit->second;

        if (!first_group)
            ofs << ",";
        first_group = false;

        ofs << "\"" << lid << "\":{";
        ofs  << "\"lead\":["
             << "\"" << lid << "\","
             << "\"" << seq << "\""
             << "],";

        for (int i = 0; i < 2; i++) {
            ofs << "\"pos-" << comp_code[i] << "\":[";
            bool first_pos = true;
            for (int p: (*_l2p[gid])[lid][i]) {
                if (!first_pos)
                    ofs << ",";
                first_pos = false;
                ofs << p;
            }
            ofs << "],";

            ofs << "\"kmers-" << comp_code[i] << "\":[";
            bool first_kmer = true;
            for (string kmer: it->second[i]) {
                int kept_n = keep[0][kmer];
                int dropped_n = drop[0][kmer];
                int kept_t = keep[1][kmer];
                int dropped_t = drop[1][kmer];
                if (!first_kmer)
                    ofs << ",";
                first_kmer = false;
                ofs << "[\"" << kmer << "\"," << kept_n << "," << kept_t << ","
                     << dropped_n << "," << dropped_t << "]";
            }
            ofs << "],";
        }

        for (int i = 0; i < 2; i++) {
            bool first_read = true;
            ofs << "\"reads-" << kind_code[i] << "\":[";
            for (string sid: (*_l2i[gid])[lid][i]) {
                i2r_table::const_iterator sit = _i2r[i]->find(sid);
                if (sit == _i2r[i]->end())
                    continue;
                if (!first_read)
                    ofs << ",";
                first_read = false;
                sm_read read = sit->second;
                seq = "";
                decode_read(read, seq);
                ofs << "[\"" << sid << "\",\"" << seq << "\"]";
            }
            ofs << ( i == 1 ? "]" : "]," );
        }

        ofs << "}";

        if (num_leads % 100 == 0) {
            end = std::chrono::system_clock::now();
            time = end - start;
            cout << "P: " << gid << " " << num_leads << " "
                 << time.count() << "\n";
            start = std::chrono::system_clock::now();
        }
    }

    ofs << "}";
}

void group::populate_index(int gid, string& lid,
                           const std::vector<string>& kmers, int kind,
                           index_count& keep, index_count& drop)
{
    for (string kmer: kmers) {
        kmer_table::const_iterator it = _k2i[kind]->find(kmer);
        if (it == _k2i[kind]->end()) {
            continue;
        }

        string list = it->second;
        boost::trim_if(list, boost::is_any_of(" "));
        if (list.size() == 0) {
            continue;
        }

        std::unordered_set<string> sids;
        boost::split(sids, list, boost::is_any_of(" "));

        if (sids.size() > DROP) {
            drop[kind][kmer] += sids.size();
            continue;
        }

        (*_l2i[gid])[lid][kind].insert(sids.begin(), sids.end());
        keep[kind][kmer] += sids.size();
    }
}

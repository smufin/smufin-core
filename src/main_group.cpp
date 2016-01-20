#include <getopt.h>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <chrono>
#include <zlib.h>
#include <kseq.h>
#include <google/sparse_hash_map>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using namespace boost::iostreams;

#define RLEN 100
#define KLEN 30
#define KMIN 0
#define KMAX 100
#define WMIN 7
#define WLEN 10
#define DROP 500

KSEQ_INIT(gzFile, gzread);

typedef std::pair<string, string> read_value;
typedef std::array<read_value, 2> i2r_value;
typedef google::sparse_hash_map<string, i2r_value> i2r_table;
typedef std::array<std::vector<string>, 2> k2i_value;
typedef google::sparse_hash_map<string, k2i_value> k2i_table;
typedef std::array<std::vector<int>, 2> pos_value;
typedef google::sparse_hash_map<string, pos_value> l2p_table;
typedef std::array<std::vector<string>, 2> kmer_value;
typedef google::sparse_hash_map<string, kmer_value> l2k_table;
typedef std::array<std::unordered_set<string>, 2> id_value;
typedef google::sparse_hash_map<string, id_value> l2i_table;

typedef std::array<std::unordered_map<string, int>, 2> index_count;

i2r_table i2r;
k2i_table k2i;
l2p_table l2p;
l2k_table l2k;
l2i_table l2i;

const char comp_code[] = "ab";
const char kind_code[] = "nt";

void rrevcomp(char read[])
{
    const char comp[] = "-------------------------------------------"
                        "----------------------T-G---C------N-----A";
    int c, i, j;
    for (i = 0, j = RLEN - 1; i < j; i++, j--) {
        c = read[i];
        read[i] = comp[read[j]];
        read[j] = comp[c];
    }
}

void load_fq(string file, int kind)
{
    int len;
    int nreads = 0;
    gzFile in = gzopen(file.c_str(), "rb");
    kseq_t *seq = kseq_init(in);
    while ((len = kseq_read(seq)) >= 0) {
        nreads++;
        string id(seq->name.s);
        string s(seq->seq.s);
        string q(seq->qual.s);
        i2r[id][kind] = read_value(s, q);
    }
    kseq_destroy(seq);
    gzclose(in);
}

void load_k2i(string file, int kind)
{
    std::ifstream ifs(file, std::ios_base::in | std::ios_base::binary);
    filtering_istream in;
    in.push(gzip_decompressor());
    in.push(ifs);
    string kmer;
    string sid;
    while (in >> kmer) {
        in >> sid;
        k2i[kmer][kind].push_back(sid);
    }
}

bool match_window(std::vector<int> pos)
{
    if (pos.size() < WMIN) {
        return false;
    }

    for (int i = 0; i <= pos.size() - WMIN; i++) {
        std::vector<int> sub(pos.begin() + i, pos.begin() + i + WMIN);
        if (sub.back() - sub.front() < WLEN) {
            return true;
        }
    }

    return false;
}

void select_candidate(string sid, string seq, std::vector<int>& pos, int dir)
{
    l2p[sid][dir] = pos;
    // TODO: Switch to unordered_set instead of vector?
    std::vector<string> kmers;
    for (int p: pos) {
        kmers.push_back(seq.substr(p, KLEN));
    }
    l2k[sid][dir] = kmers;
}

void populate_index(string& lid, const std::vector<string>& kmers, int kind,
                    index_count& keep, index_count& drop)
{
    for (string kmer: kmers) {
        k2i_table::const_iterator kit = k2i.find(kmer);
        if (kit == k2i.end()) {
            continue;
        }
        std::unordered_set<string> sids = std::unordered_set<string>();
        for (string sid: kit->second[kind]) {
            sids.insert(sid);
        }
        if (sids.size() > DROP) {
            drop[kind][kmer] += sids.size();
            continue;
        }
        l2i[lid][kind].insert(sids.begin(), sids.end());
        keep[kind][kmer] += sids.size();
    }
}

void display_usage()
{
    cout << "Usage: sm-group [OPTIONS] -i INPUT_PATH" << endl;
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

    std::ios_base::sync_with_stdio(false);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> time;

    i2r.resize(500000000);
    start = std::chrono::system_clock::now();
    load_fq(input_path + "filter-nn.fastq.gz", 0);
    load_fq(input_path + "filter-tn.fastq.gz", 1);
    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "FASTQ read time: " << time.count() << endl;

    k2i.resize(4000000);
    start = std::chrono::system_clock::now();
    load_k2i(input_path + "filter-nn.k2i.gz", 0);
    load_k2i(input_path + "filter-tn.k2i.gz", 1);
    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "K2I read time: " << time.count() << endl;

    start = std::chrono::system_clock::now();
    std::ifstream in(input_path + "filter-tm.i2p");
    string sid;
    int flen, rlen;
    while (in >> sid >> flen >> rlen) {
        int pos;
        std::vector<int> fpos;
        std::vector<int> rpos;

        for (int i = 0; i < flen; i++) {
            in >> pos;
            fpos.push_back(pos);
        }
        for (int i = 0; i < rlen; i++) {
            in >> pos;
            rpos.push_back(RLEN - KLEN - pos);
        }

        i2r_table::const_iterator it = i2r.find(sid);
        if (it == i2r.end())
            continue;

        read_value read = it->second[1];

        if (read.first.size() == 0)
            continue;

        if (read.first.find("N") != std::string::npos)
            continue;

        if (flen >= KMIN && flen <= KMAX && rlen > 0 && match_window(fpos)) {
            select_candidate(sid, read.first, fpos, 0);
        }

        std::reverse(rpos.begin(), rpos.end());

        if (rlen >= KMIN && rlen <= KMAX && flen > 0 && match_window(rpos)) {
            char buf[RLEN + 1];
            copy(read.first.begin(), read.first.end(), buf);
            buf[RLEN] = '\0';
            rrevcomp(buf);
            select_candidate(sid, string(buf), rpos, 1);
        }
    }
    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "Candidate selection time: " << time.count() << endl;

    cout << "{";
    start = std::chrono::system_clock::now();
    bool first_group = true;
    for (l2k_table::const_iterator it = l2k.begin(); it != l2k.end(); ++it) {
        string lid = it->first;
        index_count keep;
        index_count drop;

        for (int i = 0; i < 2; i++) {
            for (string kmer: it->second[i]) {
                keep[i][kmer] = 0;
                drop[i][kmer] = 0;
            }
        }

        populate_index(lid, it->second[0], 0, keep, drop);
        populate_index(lid, it->second[0], 1, keep, drop);
        populate_index(lid, it->second[1], 0, keep, drop);
        populate_index(lid, it->second[1], 1, keep, drop);

        if (!first_group)
            cout << ",";
        first_group = false;

        cout << "\"" << lid << "\":{";
        cout << "\"lead\":["
             << "\"" << lid << "\","
             << "\"" << i2r[lid][1].first << "\","
             << "\"" << i2r[lid][1].second << "\""
             << "],";

        for (int i = 0; i < 2; i++) {
            cout << "\"pos-" << comp_code[i] << "\":[";
            bool first_pos = true;
            for (int p: l2p[lid][i]) {
                if (!first_pos)
                    cout << ",";
                first_pos = false;
                cout << p;
            }
            cout << "],";

            cout << "\"kmers-" << comp_code[i] << "\":[";
            bool first_kmer = true;
            for (string kmer: it->second[i]) {
                int kept_n = keep[0][kmer];
                int dropped_n = drop[0][kmer];
                int kept_t = keep[1][kmer];
                int dropped_t = drop[1][kmer];
                if (!first_kmer)
                    cout << ",";
                first_kmer = false;
                cout << "[\"" << kmer << "\"," << kept_n << "," << kept_t << ","
                     << dropped_n << "," << dropped_t << "]";
            }
            cout << "],";
        }

        for (int i = 0; i < 2; i++) {
            bool first_read = true;
            cout << "\"reads-" << kind_code[i] << "\":[";
            for (string read: l2i[lid][i]) {
                if (!first_read)
                    cout << ",";
                first_read = false;
                cout << "[\"" << read << "\",\"" << i2r[read][i].first
                     << "\",\"" << i2r[read][i].second << "\"]";
            }
            cout << ( i == 1 ? "]" : "]," );
        }

        cout << "}";
    }
    cout << "}";
    end = std::chrono::system_clock::now();
    time = end - start;
    cerr << "Populate candidates time: " << time.count() << endl;
}

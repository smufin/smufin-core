#include "input.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include <sys/stat.h>

#include <htslib/sam.h>
#include <htslib/bgzf.h>

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

    ifs.close();
}

bool input_queue::try_dequeue(sm_chunk &chunk)
{
    return _queue.try_dequeue(chunk);
}

void input_queue_bam_chunks::init()
{
    std::ifstream ifs(_conf.input_file);
    if (!ifs.good()) {
        cout << "Invalid input file" << endl;
        exit(1);
    }

    for (string file; std::getline(ifs, file);) {
        sm_read_kind kind = NORMAL_READ;
        std::size_t found = file.find("_T_");
        if (found != std::string::npos) {
            kind = CANCER_READ;
        }

        std::vector<uint64_t> offsets(_conf.num_loaders);
        if (!chunk_bam(file, _conf.num_loaders, offsets)) {
            cout << "Failed to chunk BAM file " << file << endl;
            exit(1);
        }

        for (int i = 0; i < offsets.size(); i++) {
            sm_chunk chunk;
            chunk.file = file;
            chunk.begin = offsets[i];
            chunk.end = offsets[i + 1];
            chunk.kind = kind;
            _queue.enqueue(chunk);
            len++;
        }
    }

    ifs.close();
}

bool input_queue_bam_chunks::chunk_bam(const string bam_file,
                                       const int num_chunks,
                                       std::vector<uint64_t> &offsets)
{
    struct stat st;
    stat(bam_file.c_str(), &st);

    const size_t bam_file_size = st.st_size;
    const uint64_t chunk_size = bam_file_size / num_chunks;
    const string bai_file = string(bam_file + ".bai");

    for (auto& o: offsets)
        o = (chunk_size * num_chunks) << 16;

    BGZF *fp = bgzf_open(bai_file.c_str(), "r");
    uint8_t magic[4];
    uint32_t n_ref;

    if (fp == NULL)
        return false;

    if (bgzf_read(fp, magic, 4) != 4)
        return false;

    if (memcmp(magic, "BAI\1", 4) != 0)
        return false;

    if (bgzf_read(fp, &n_ref, 4) != 4)
        return false;

    for (int i = 0; i < n_ref; i++) {
        uint32_t n_bin;
        if (bgzf_read(fp, &n_bin, 4) != 4)
            return false;

        for (int j = 0; j < n_bin; j++) {
            uint32_t bin, n_chunk;

            if (bgzf_read(fp, &bin, 4) != 4)
                return false;

            if (bgzf_read(fp, &n_chunk, 4) != 4)
                return false;

            size_t off_size = n_chunk << 4;
            uint64_t* off_list = (uint64_t*) malloc(off_size);
            if (bgzf_read(fp, off_list, off_size) != off_size)
                return false;
        }

        uint32_t n_intv;
        if (bgzf_read(fp, &n_intv, 4) != 4)
            return false;

        size_t off_size = n_intv * sizeof(uint64_t);
        uint64_t* off_list = (uint64_t*) malloc(off_size);
        if (bgzf_read(fp, off_list, off_size) != off_size)
            return false;

        for (int j = 0; j < n_intv; j++) {
            uint64_t new_offset = off_list[j] >> 16;
            int bin = new_offset / chunk_size;
            uint64_t old_offset = offsets[bin] >> 16;
            if (new_offset < old_offset)
                offsets[bin] = off_list[j];
        }
    }

    bgzf_close(fp);

    offsets.push_back(bam_file_size << 16);
    return true;
}

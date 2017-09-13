#include "input.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/stat.h>

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <boost/algorithm/string.hpp>

using std::cout;
using std::endl;
using std::string;

void input_queue::init()
{
    // Interleave normal and tumoral files.
    std::vector<std::pair<string, sm_read_kind>> files;
    int min = std::min(_conf.list_normal.size(), _conf.list_tumor.size());
    for (int i = 0; i < min; i++) {
        files.push_back({_conf.list_normal[i], NORMAL_READ});
        files.push_back({_conf.list_tumor[i], CANCER_READ});
    }

    for (int i = min; i < _conf.list_normal.size(); i++)
        files.push_back({_conf.list_normal[i], NORMAL_READ});
    for (int i = min; i < _conf.list_tumor.size(); i++)
        files.push_back({_conf.list_tumor[i], CANCER_READ});

    for (auto& file: files) {
        sm_chunk chunk;
        chunk.file = file.first;
        chunk.begin = -1;
        chunk.end = -1;
        chunk.kind = file.second;
        _queue.enqueue(chunk);
        len++;
    }
}

bool input_queue::try_dequeue(sm_chunk &chunk)
{
    return _queue.try_dequeue(chunk);
}

void input_queue_bam_chunks::init()
{
    std::vector<std::pair<string, sm_read_kind>> files;
    for (auto& file: _conf.list_normal)
        files.push_back({file, NORMAL_READ});
    for (auto& file: _conf.list_tumor)
        files.push_back({file, CANCER_READ});

    for (auto& file: files) {
        std::vector<uint64_t> offsets;
        if (!chunk_bam(file.first, _conf.num_loaders, offsets)) {
            cout << "Failed to chunk BAM file " << file.first << endl;

            // Default to addressing the entire BAM file without chunking.
            sm_chunk chunk;
            chunk.file = file.first;
            chunk.begin = -1;
            chunk.end = -1;
            chunk.kind = file.second;
            _queue.enqueue(chunk);
            len++;
            continue;
        }

        for (int i = 0; i < offsets.size() - 1; i++) {
            sm_chunk chunk;
            chunk.file = file.first;
            chunk.begin = offsets[i];
            chunk.end = offsets[i + 1];
            chunk.kind = file.second;
            _queue.enqueue(chunk);
            len++;
        }
    }
}

// Manually reads the BAI index of «bam_file», jumping to its linear index and
// trying to find «num_chunks» addresses (offsets) at approximately the same
// distance from one another. This ensures that all BGZF blocks are read,
// including unmapped alignments.
bool input_queue_bam_chunks::chunk_bam(const string bam_file,
                                       const int num_chunks,
                                       std::vector<uint64_t> &offsets)
{
    struct stat st;
    stat(bam_file.c_str(), &st);

    const size_t bam_file_size = st.st_size;
    const uint64_t chunk_size = bam_file_size / num_chunks;
    const string bai_file = string(bam_file + ".bai");

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

    for (int i = 0; i < num_chunks; i++)
        offsets.push_back((chunk_size * num_chunks) << 16);

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

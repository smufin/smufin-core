#include "input_iterator_bam.hpp"

#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

input_iterator_bam::input_iterator_bam(const sm_config &conf,
                                       const sm_chunk &chunk)
    : _conf(conf), _chunk(chunk)
{
    _in = sam_open(chunk.file.c_str(), "r");
    _record = bam_init1();

    BGZF *fp = _in->fp.bgzf;
    bam_hdr_t *header = bam_hdr_read(fp);

    uint64_t begin_address = chunk.begin >> 16;
    int begin_offset = chunk.begin & 0xFFFF;
    uint64_t end_address = chunk.end >> 16;
    int end_offset = chunk.end & 0xFFFF;

    if (chunk.begin > 0) {
        uint64_t seek = bgzf_seek(fp, chunk.begin, SEEK_SET);
    }
}

bool input_iterator_bam::next(sm_read *read)
{
    int len;
    BGZF *fp = _in->fp.bgzf;

    while((len = bam_read1(fp, _record)) >= 0 && bgzf_tell(fp) <= _chunk.end) {
        const int read_len = _record->core.l_qseq;
        const uint8_t *s = bam_get_seq(_record);
        const uint8_t *q = bam_get_qual(_record);

        for (int i = 0; i < read_len; i++) {
            _seq[i] = seq_nt16_str[bam_seqi(s, i)];
            _qual[i] = 33 + q[i];
        }

        if (bam_is_rev(_record)) {
            revcomp(_seq, read_len);
            rev(_qual, read_len);
        }

        _seq[read_len] = '\0';
        _qual[read_len] = '\0';

        string dir;
        if (_record->core.flag & BAM_FREAD1)
            dir = "/1";
        else if (_record->core.flag & BAM_FREAD2)
            dir = "/2";

        if (lq_count(_qual, read_len) > read_len / 10)
            continue;

        sprintf(_id, "%s%s", bam_get_qname(_record), dir.c_str());
        read->id = _id;
        read->seq = _seq;
        read->qual = _qual;
        read->len = read_len;

        read->num_splits = 0;
        int p = 0;
        int l = read_len;
        int n = read_len;
        char *ps;

        while ((ps = (char*) memchr(&_seq[p], 'N', l - p)) != NULL) {
            n = ps - &_seq[p];
            if (n > 0) {
                read->splits[read->num_splits][0] = p;
                read->splits[read->num_splits][1] = n;
                read->num_splits++;
                p += n;
            }
            p++;
        }

        n = l - p;
        if (n > 0) {
            read->splits[read->num_splits][0] = p;
            read->splits[read->num_splits][1] = n;
            read->num_splits++;
        }

        return true;
    }

    return false;
}

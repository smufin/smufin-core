/*
 * Copyright © 2015-2019 Barcelona Supercomputing Center (BSC)
 *
 * This file is part of SMUFIN Core. SMUFIN Core is released under the SMUFIN
 * Public License, and may not be used except in compliance with it. This file
 * is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; see the SMUFIN Public License for more details. You should have
 * received a copy of the SMUFIN Public License along with this file. If not,
 * see <https://github.com/smufin/smufin-core/blob/master/COPYING>.
 *
 * Jordà Polo <jorda.polo@bsc.es>, 2015-2018
 */

#include "input_iterator_bam.hpp"

#include <fstream>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

input_iterator_bam::input_iterator_bam(const sm_config &conf,
                                       const sm_chunk &chunk)
    : input_iterator(conf, chunk)
{
    check = conf.check_quality;

    _in = sam_open(chunk.file.c_str(), "r");
    _header = sam_hdr_read(_in);
    _record = bam_init1();

    uint64_t begin_address = chunk.begin >> 16;
    int begin_offset = chunk.begin & 0xFFFF;
    uint64_t end_address = chunk.end >> 16;
    int end_offset = chunk.end & 0xFFFF;

    if (chunk.begin > 0) {
        uint64_t seek = bgzf_seek(_in->fp.bgzf, chunk.begin, SEEK_SET);
    }
}

bool input_iterator_bam::next(sm_read *read)
{
    int len;
    while((len = sam_read1(_in, _header, _record)) >= 0 &&
          bgzf_tell(_in->fp.bgzf) <= _chunk.end) {
        const int read_len = _record->core.l_qseq;
        assert(read_len <= MAX_READ_LEN);

        const uint8_t *s = bam_get_seq(_record);
        const uint8_t *q = bam_get_qual(_record);

        // Skip non-primary and supplementary alignments.
        if (_record->core.flag & 256 || _record->core.flag & 2048)
            continue;

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

        if (check && lq_count(_qual, read_len) > (read_len / 10))
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
            if (n >= _conf.k) {
                assert(read->num_splits < MAX_SPLITS);
                read->splits[read->num_splits][0] = p;
                read->splits[read->num_splits][1] = n;
                read->num_splits++;
            }
            p += n + 1;
        }

        n = l - p;
        if (n >= _conf.k) {
            read->splits[read->num_splits][0] = p;
            read->splits[read->num_splits][1] = n;
            read->num_splits++;
        }

        return true;
    }

    return false;
}

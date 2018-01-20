/*
 * Copyright © 2015-2018 Barcelona Supercomputing Center (BSC)
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

#ifndef __SM_STAGE_H__
#define __SM_STAGE_H__

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "common.hpp"

class stage {
public:
    stage(const sm_config &conf) : _conf(conf) {};
    virtual void chain(const stage* prev) {};
    void exec(std::string name, const std::vector<std::string> &steps);

    template<typename T> static stage* create(const sm_config &conf)
    {
        return new T(conf);
    };

protected:
    const sm_config &_conf;
    std::map<std::string, std::function<void(void)>> _executable;
};

typedef std::function<stage*(const sm_config &conf)> stage_s;

#endif

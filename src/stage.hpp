#ifndef __SM_STAGE_H__
#define __SM_STAGE_H__

#include <map>
#include <string>
#include <vector>

#include "common.hpp"

class stage {
public:
    stage(const sm_config &conf) : _conf(conf) {};
    virtual void chain(const stage* prev) {};
    void exec(std::string name, const std::vector<std::string> &steps);

protected:
    const sm_config &_conf;
    std::map<std::string, std::function<void(void)>> _executable;
};

#endif

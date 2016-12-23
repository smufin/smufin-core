#include "db.hpp"

#include <sstream>
#include <string>

#include "common.hpp"

rocksdb::Options get_rocks_options(std::string type)
{
    rocksdb::Options options;
    options.create_if_missing = true;
    options.IncreaseParallelism(4);

    if (type == "k2i") {
        options.merge_operator.reset(new IDListOperator());
    }

    if (type == "i2p") {
        options.merge_operator.reset(new PositionsMapOperator());
    }

    return options;
}

void encode_pos(std::string &s, sm_pos_bitmap &p)
{
    std::stringstream e;
    e << std::hex << p.a[0] << " " << std::hex << p.a[1] << " "
      << std::hex << p.b[0] << " " << std::hex << p.b[1];
    s = e.str();
}

sm_pos_bitmap decode_pos(char const *s)
{
   sm_pos_bitmap p;
   std::istringstream(s)
       >> std::hex >> p.a[0] >> std::hex >> p.a[1]
       >> std::hex >> p.b[0] >> std::hex >> p.b[1];
   return p;
}

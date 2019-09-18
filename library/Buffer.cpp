#include "Buffer.h"

#include <cassert>

void Buffer::flush()
{
    fields.clear();
}

void Buffer::add(const vecDbl &vec)
{
    fields.push_back(vec);
}

void Buffer::write(const vecDbl &vec, std::ofstream * &out)
{
    if (vec.size() > 0u) {
        for (auto x : vec) {
            //assert(x); // nullptr access violation around here
            out->write((char *) &x, sizeof(x));
        }
    }
}

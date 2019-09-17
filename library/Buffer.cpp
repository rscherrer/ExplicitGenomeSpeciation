#include "Buffer.h"

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
    if (vec.size() > 0u)
        for (auto x : vec)
            out->write((char *) &x, sizeof(x));
}

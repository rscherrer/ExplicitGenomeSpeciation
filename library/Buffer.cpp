#include "Buffer.h"

void Buffer::flush()
{
    fields.clear();
}

void Buffer::add(const double &x)
{
    fields.push_back(x);
}

void Buffer::write(std::ofstream * &out, const double &value)
{
    out->write((char *) &value, sizeof(value));
}

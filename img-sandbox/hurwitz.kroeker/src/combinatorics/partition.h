#include <vector>
#include <iostream>
#include <cstdio>
#include <assert.h>
#include <iterator>


///todo: templatize std::vector<T> and at the same time shrink possible parameters (integers)
inline bool next_partition_desc(std::vector<int> *p_ptr)
{
    assert(p_ptr != NULL);
    std::vector<int> &p = *p_ptr;
    if (p.size() < 2)
        return false;

    int sum = p.back();
    int i;
    for (i = p.size() - 2; i > 0 && p[i] + 1 > p[i - 1]; --i)
    {
        sum += p[i];
    }
    --sum;
    p[i]++;
    p.resize(i + 1 + sum);
    std::fill(p.begin() + i + 1, p.begin() + i + 1 + sum, 1);
    
    return true;
}

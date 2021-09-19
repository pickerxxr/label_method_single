#pragma once

#include <utility>
#include <chrono>
#include <stdint.h>
#include <sys/stat.h>
#include <iostream>
#include <sstream>


typedef uint32_t vid_t;
const vid_t INVALID_VID = -1;

struct edge_t {
    vid_t first, second;
    edge_t() : first(0), second(0) {}
    edge_t(vid_t first, vid_t second) : first(first), second(second) {} // 构造函数
    const bool valid() { return first != INVALID_VID; } // 查看边是否符合规则 根据规则 如果边不对应该将第一个设定为-1
    void remove() { first = INVALID_VID; } // 如上所述，将某个边删除只需要将其第一个节点设定为 -1
};

struct label_mirror{
    vid_t select_label;
    int max_mirror;
    int num_adj_labeled;
    int num_adj_adj_labeled;
};


inline bool is_exists(const std::string &name)
{
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

inline std::string degree_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".degree";
    return ss.str();
}

inline std::string binedgelist_name(const std::string &basefilename)
{
    std::stringstream ss;
    ss << basefilename << ".binedgelist";
    return ss.str();
}

class Timer
{
  private:
    std::chrono::system_clock::time_point t1, t2;
    double total;

  public:
    Timer() : total(0) {}
    void reset() { total = 0; }
    void start() { t1 = std::chrono::system_clock::now(); }
    void stop()
    {
        t2 = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = t2 - t1;
        total += diff.count();
    }
    double get_time() { return total; }
};

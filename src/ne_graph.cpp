#include "ne_graph.hpp"

void ne_graph_t::build(const std::vector<edge_t> &edges) // 构建节点、边、邻接列表的数据结构
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].first]++;

    vdata[0] = ne_adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = ne_adjlist_t(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].first].push_back(i);
}

void ne_graph_t::build_reverse(const std::vector<edge_t> &edges) // 对于有向图来说需要反向构建一次
{
    if (edges.size() > nedges)
        neighbors = (uint40_t *)realloc(neighbors, sizeof(uint40_t) * edges.size());
    nedges = edges.size();

    std::vector<size_t> count(num_vertices, 0);
    for (size_t i = 0; i < nedges; i++)
        count[edges[i].second]++;

    vdata[0] = ne_adjlist_t(neighbors);
    for (vid_t v = 1; v < num_vertices; v++) {
        count[v] += count[v-1];
        vdata[v] = ne_adjlist_t(neighbors + count[v-1]);
    }
    for (size_t i = 0; i < edges.size(); i++)
        vdata[edges[i].second].push_back(i);
}
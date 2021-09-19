#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>

#include "util.hpp"
#include "label_max_heap.hpp"
#include "dense_bitset.hpp"
#include "partitioner.hpp"
#include "ne_graph.hpp"

/* Neighbor Expansion (NE) */
class LabelAssign : public Partitioner
        {
        private:
            const double BALANCE_RATIO = 1.00;

            std::string basefilename;

            vid_t num_vertices;
            size_t num_edges, assigned_edges;
            int p, bucket;
            int batch;
            double average_degree;
            size_t capacity;
            size_t batch_capacity;

            std::vector<edge_t> edges;
            ne_graph_t adj_out, adj_in;
            LabelMaxHeap<vid_t, vid_t> max_heap;
            std::vector<size_t> occupied;
            std::vector<vid_t> degrees;
            std::vector<std::vector<vid_t>> assigned_edges_set;
            std::vector<vid_t> v_label;
            std::vector<std::vector<vid_t> > v_label_set;
            std::vector<int8_t> master;
            std::vector<bool> adj_to_labeled;
            std::vector<bool> is_labeled;
            std::vector<vid_t> labeled_neighbor_num;
//            LabelMaxHeap<vid_t, vid_t> to_be_label_heap;
//            // TODO: find the label range
            std::vector<LabelMaxHeap<vid_t, vid_t>> to_be_labeled;
            int batch_number_bound = 20000;
            double lambda;
            std::random_device rd;
            std::mt19937 gen;
            std::uniform_int_distribution<vid_t> dis;

            void add_label_set(vid_t vid, vid_t add_vid_label){
                if (!in_set(add_vid_label, v_label_set[vid])){
                    v_label_set[vid].push_back(add_vid_label);
                }
            }

            std::vector<vid_t> find_available_label_set(vid_t vid){
                std::vector<vid_t> res;
                for (int direction = 0; direction < 2; ++direction){
                    ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                    for (size_t i = 0; i < neighbors.size(); i++){
                        vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                        for (auto iter = v_label_set[u].begin(); iter != v_label_set[u].end(); ++iter){
                            if (!in_set(*iter, res)){
                                res.push_back(*iter);
                            }
                        }
                    }
                }
                return res;
            }

            static long get_vidt_index(std::vector<vid_t> v, vid_t K)
            {
                auto it = find(v.begin(), v.end(), K);
                // If element was found
                if (it != v.end())
                {
                    // calculating the index
                    long index = it - v.begin();
                    return index;
                }
                else {
                    return -1;
                }
            }

            static long get_int_index(std::vector<int> v, int K)
            {
                auto it = find(v.begin(), v.end(), K);
                // If element was found
                if (it != v.end())
                {
                    // calculating the index
                    long index = it - v.begin();
                    return index;
                }
                else {
                    return -1;
                }
            }



            label_mirror min_rf_label(vid_t vid){
                std::vector<vid_t> available_set = find_available_label_set(vid);
                std::vector<int> count_set;
                int adj_labeled = 0;
                int adj_adj_labeled = 0;
                int max_element;
                count_set.resize(available_set.size());
                count_set.assign(available_set.size(), 0);
                vid_t res;
                for (auto iter = available_set.begin(); iter != available_set.end(); ++iter){
                    for (int direction = 0; direction < 2; ++direction) {
                        ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                        for (size_t i = 0; i < neighbors.size(); i++) {
                            vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                            // find the most popular one
                            if (in_set(*iter, v_label_set[u])) {
                                long idx = get_vidt_index(available_set, *iter);
                                count_set[idx]++;
                            }
                            if (is_labeled[u]){
                                adj_labeled ++;
                            }
                            if (adj_to_labeled[u]){
                                adj_adj_labeled ++;
                            }
                        }
                    }
                }
                max_element = get_max_element(count_set);
                res = available_set[get_int_index(count_set, max_element)];
                label_mirror max_label_mirror = {res, max_element, adj_labeled, adj_adj_labeled};
                return max_label_mirror;
            }


            int get_max_element(std::vector<int> s){
                int max_res = *s.begin();
                for (auto iter = s.begin(); iter != s.end(); ++iter){
                    if (*iter > max_res){
                        max_res = *iter;
                    }
                }
                return max_res;
            }

            void assign_labeled(vid_t vid, vid_t d, vid_t _label){
                if (d == 0) return;
                v_label[vid] = _label;
                for (int direction = 0; direction < 2; direction++){
                    ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];

                    for (size_t i = 0; i < neighbors.size(); i++){
                        vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                        if (is_labeled[u]){
                            drop_bucket(direction ? vid : u, direction ? u : vid, _label);
                        }
                    }
                }
            }


            bool get_free_vertex(vid_t &vid)
            {
                vid = dis(gen);
                vid_t count = 0;
                while ((count < num_vertices &&
                (adj_out[vid].size() + adj_in[vid].size() == 0 ||
                adj_out[vid].size() + adj_in[vid].size() >
                2 * average_degree ||
                is_labeled[vid])) ||
                adj_to_labeled[vid]){
                    vid = (vid + ++count) % num_vertices;
                }
                if (count == num_vertices)
                    return false;
                return true;
            }

            bool in_set(vid_t tar, std::vector<vid_t> _set){
                if (std::find(_set.begin(), _set.end(), tar) != _set.end()){
                    return true;
                }else{
                    return false;
                }
            }

            void drop_bucket(vid_t from, vid_t to, vid_t _bucket){
                occupied[_bucket]++;
                degrees[from]--;
                degrees[to]--;
                assigned_edges++;
                if (!in_set(from, assigned_edges_set[_bucket])){ assigned_edges_set[_bucket].push_back(from); }
                if (!in_set(to, assigned_edges_set[_bucket])){ assigned_edges_set[_bucket].push_back(to); }
            }

            std::vector<std::vector<vid_t>> greedy_assign(std::vector<std::vector<vid_t>> finish_edge_set, int final_partition_num){
                std::vector<std::vector<vid_t>> res;
                res.resize(final_partition_num);
                std::vector<int> length_set;
                length_set.resize(finish_edge_set.size());
                for (int i = 0; i < finish_edge_set.size(); i ++){
                    length_set.push_back(finish_edge_set[i].size());
                }
                std::vector<int> desc_index_set;
                while (!length_set.empty()){
                    int max_l = 0;
                    int max_index = 0;
                    for (int i = 0; i < length_set.size(); i++) {
                        if (length_set[i] > max_l){
                            max_l = length_set[i];
                            max_index = i;
                        }
                    }
                    desc_index_set.push_back(max_index);
                    length_set.erase (length_set.begin() + max_index);
                }
                // keep a record of final length for greedy
                std::vector<int> length_keep;
                length_keep.assign(final_partition_num, 0);
                for (int i = 0; i < desc_index_set.size(); i ++){
                    int edge_set_index = desc_index_set[i];
                    // following: greedy assign each bucket to final partition...
                    int min_partition;
                    min_partition = get_min_length_index(length_keep);
                    // TODO The problem is in the func above.
                    length_keep[min_partition] += desc_index_set[edge_set_index];
                    for (int j = 0; j < finish_edge_set[edge_set_index].size(); j ++){
                        if (!in_set(finish_edge_set[edge_set_index][j], res[min_partition])){
                            res[min_partition].push_back(finish_edge_set[edge_set_index][j]);
                        }
                    }
                }
                return res;
            }


            static int get_min_length_index(std::vector<int> _set){
                int min_index = 0;
                int min_value = 99999999;
                for (int i = 0; i < _set.size(); i ++){
                    if (_set[i] < min_value){
                        min_index = i;
                        min_value = _set[i];
                    }
                }
                return min_index;
            }

            double compute_rf(std::vector<std::vector<vid_t>> _res) const{
                double all_mirrors_num = 0;
                for (auto iter = _res.begin(); iter != _res.end(); iter++){
                    all_mirrors_num += iter->size();
                }
                double rf;
                rf = (1.0 * all_mirrors_num) / num_vertices;
                return rf;
            }

        public:
            LabelAssign(std::string basefilename);
            void split();
        };
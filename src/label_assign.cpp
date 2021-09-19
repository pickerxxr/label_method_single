#include "label_assign.hpp"
#include "conversions.hpp"


LabelAssign::LabelAssign(std::string basefilename)
: basefilename(basefilename), rd(), gen(rd())
{
    Timer convert_timer;
    convert_timer.start();
    Converter *converter = new Converter(basefilename);
    convert(basefilename, converter);
    delete converter;
    convert_timer.stop();

    total_time.start();
    std::cout << "initializing partitioner" << std::endl;
    std::string filename;

    filename = binedgelist_name(basefilename);

    std::ifstream fin(filename,
                      std::ios::binary | std::ios::ate);
    auto filesize = fin.tellg();
    std::cout << "file size: " << filesize;
    fin.seekg(0, std::ios::beg);

    fin.read((char *)&num_vertices, sizeof(num_vertices));
    fin.read((char *)&num_edges, sizeof(num_edges));

    std::cout << "num_vertices: " << num_vertices
    << ", num_edges: " << num_edges;

    p = 32;

    lambda = 1.0;
    average_degree = (double)num_edges * 2 / num_vertices;
    batch_capacity = num_edges / p;
    is_labeled.resize(num_vertices);
    occupied.assign(batch_number_bound, 0);
    adj_out.resize(num_vertices);
    adj_in.resize(num_vertices);
    v_label_set.resize(num_vertices);
    assigned_edges_set.resize(batch_number_bound * 2);
    dis.param(std::uniform_int_distribution<vid_t>::param_type(0, num_vertices - 1));

    Timer read_timer;
    read_timer.start();
    std::cout << "loading...";
    edges.resize(num_edges);
    fin.read((char *)&edges[0], sizeof(edge_t) * num_edges);

    std::cout << "constructing...";
    adj_out.build(edges);
    adj_in.build_reverse(edges);
    adj_to_labeled.resize(num_vertices);
    degrees.resize(num_vertices);

    labeled_neighbor_num.resize(num_vertices);
    labeled_neighbor_num.assign(num_vertices, 0);
    std::ifstream degree_file(degree_name(basefilename), std::ios::binary);


    degree_file.read((char *)&degrees[0], num_vertices * sizeof(vid_t));
    degree_file.close();
    read_timer.stop();
    std::cout << "time used for graph input and construction: " << read_timer.get_time();
};

/*
 * main label algorithm implementation
 */
void LabelAssign::split()
{
    std::cout << "starting the label process..." << std::endl;
    Timer computing_timer;
    max_heap.reserve(num_vertices);
    std::cout << "give random label first... " << std::endl;
    // 初始化
    v_label.assign(num_vertices, 0);
    vid_t max_label = 1;
    v_label[0] = max_label;
    v_label_set[0].push_back(max_label);
    is_labeled[0] = true;
    for (int direction = 0; direction < 2; ++direction){
        ne_adjlist_t &neighbors = direction ? adj_out[0] : adj_in[0];
        for (size_t i = 0; i < neighbors.size(); i++){
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                adj_to_labeled[u] = true;
                labeled_neighbor_num[u]++;
                v_label_set[u].push_back(1);

        }
        for (size_t i = 0; i < neighbors.size(); i++){
                vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                max_heap.insert(labeled_neighbor_num[u], u);
        }
    }
    int mirrors_record = 0;
    int assigned_vertex = 0;
    while (assigned_edges < num_edges)
    {
        vid_t vid, d;
        label_mirror min_rf_label_mirror;
        if (!max_heap.get_max(d, vid)){
            std::cout << " a random vertex born... " << std::endl;
            max_label++;
            if (!get_free_vertex(vid)){
                std::cout << "No free vertex... somewhere wrong... " << std::endl;
                break;
            }
            v_label[vid] = max_label;
            is_labeled[vid] = true;
            v_label_set[vid].push_back(max_label);
            d = adj_out[vid].size() + adj_in[vid].size();
            for (int direction = 0; direction < 2; ++direction){
                ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size(); i++){
                        vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                        adj_to_labeled[u] = true;
                        labeled_neighbor_num[u]++;
                        v_label_set[u].push_back(max_label);
                }
                for (size_t i = 0; i < neighbors.size(); i++){
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (!is_labeled[u] && !adj_to_labeled[u]) { max_heap.insert(labeled_neighbor_num[u], u); }
                }
            }
            }else{
            d = adj_out[vid].size() + adj_in[vid].size();

            // TODO: need to find the theta (partial what????)

            vid_t picked_label;
            min_rf_label_mirror = min_rf_label(vid);
            int mirror_tmp = min_rf_label_mirror.num_adj_labeled - min_rf_label_mirror.max_mirror; //real new mirror
            int mirror_if_new = min_rf_label_mirror.num_adj_adj_labeled; // if a new label given, the RF
            double theta = 1.0 * mirror_if_new / ((mirrors_record + 1) / max_label);
            if ( theta < 0.001 ) {
                picked_label = ++max_label;
                mirrors_record += mirror_if_new;
                std::cout << "max label ++ " << max_label <<std::endl;
            }else{
                picked_label = min_rf_label_mirror.select_label;
                mirrors_record += mirror_tmp;
            }


            is_labeled[vid] = true;
            if(!in_set(picked_label, v_label_set[vid])){ v_label_set[vid].push_back(picked_label); }
            for (int direction = 0; direction < 2; ++direction){
                ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size(); i++ ){
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    add_label_set(u, picked_label);
                }
            }
            assign_labeled(vid, d, picked_label);
            max_heap.remove(vid);
        }
        // 下面进行boundary的维护
        if (adj_out[vid].size() + adj_in[vid].size() > 0){
            for (int direction = 0; direction < 2; ++direction) {
                ne_adjlist_t &neighbors = direction ? adj_out[vid] : adj_in[vid];
                for (size_t i = 0; i < neighbors.size(); i++) {
                    vid_t &u = direction ? edges[neighbors[i].v].second : edges[neighbors[i].v].first;
                    if (!is_labeled[u] && adj_to_labeled[u]) {
                        labeled_neighbor_num[u]++;
                        max_heap.increase_key(u);
                        if (!in_set(v_label[vid], v_label_set[u])) { v_label_set[u].push_back(v_label[vid]); }
                    } else if (!is_labeled[u] && !adj_to_labeled[u]) {
                        max_heap.insert(1, u);
                        adj_to_labeled[u] = true;
                    }
                }
            }
        }
        assigned_vertex++;
        std::cout << "one vertex processed ... remaining: " << (num_vertices - assigned_vertex) << std::endl;
    }
    std::cout << "------ all info stats ------" << std::endl;
    std::cout << "maximum label: " << max_label << std::endl;
    std::cout << "mirrors: " << mirrors_record << std::endl;
    std::cout << "num vertices: " << num_vertices << std::endl;
    std::cout << "Replication Factor: " << (double)(1.0 * mirrors_record / num_vertices) << std::endl;

    std::cout << "********** STARTING GREEDY ASSIGN *********" << std::endl;
    std::vector<std::vector<vid_t>> partition_result;
    partition_result = greedy_assign(assigned_edges_set, 128);
    std::cout<< "Greedy assign  finished ... " << std::endl;

    // RF after greedy assignment
    int greedy_mirrors = 0;
    for (int i = 0; i < assigned_edges_set.size(); i++){
        greedy_mirrors += assigned_edges_set[i].size();
    }
    double greedy_rf = 1.0 * greedy_mirrors / num_vertices;

    // balance factor
    int greedy_max_size = 0;
    int greedy_min_size = 9999999;
    for (int i = 0; i < partition_result.size(); i++){
        if (partition_result[i].size() > greedy_max_size){
            greedy_max_size = partition_result[i].size();
        }
        if (partition_result[i].size() < greedy_min_size){
            greedy_min_size = partition_result[i].size();
        }
    }
    double balance_factor;
    balance_factor = 1.0 * greedy_max_size / greedy_min_size;
    std::cout << "Final balance: " << balance_factor << std::endl;
    std::cout << "Final replication factor: " << greedy_rf << std::endl;
}


#include <string>
#include <cstdio>
#include <sys/stat.h>
#include <iostream>
#include "util.hpp"
#include "label_assign.hpp"
#include "util.hpp"

//DECLARE_bool(help);
//DECLARE_bool(helpshort);
//
//DEFINE_int32(p, 32, "number of partitions");
//DEFINE_string(filename, "", "the file name of the input graph");
//DEFINE_string(filetype, "edgelist",
//              "the type of input file (HEP supports only 'edgelist')");
//DEFINE_string(method, "hep", "partition method: hep or ne. v2e is a special case that just transforms a given vertex partitioning to an edge partitioning.");
//DEFINE_bool(write_results, false, "Should the result be written to the output file or not. Take care, writing is in ASCII format is is really really slow in the current implementation.");
//DEFINE_bool(write_low_degree_edgelist, false, "Should the list of edges incident to a low-degree vertex be written out to a file?");
//DEFINE_double(hdf, 100, "High-degree factor: hdf * average_degree = high-degree threshold (hdth). Called \tau in the paper. Vertices with than hdth neighbors are treated specially in fast NE");
//DEFINE_double(lambda, 1.1, "Lambda value to weigh in balancing score in streaming partitioning via HDRF");
//DEFINE_bool(extended_metrics, false, "Display extended metrics in the result");
//DEFINE_bool(random_streaming, false, "Use random streaming instead of HDRF in the second phase of HEP.");
//DEFINE_bool(hybrid_NE, false, "Perform hybrid partitioning in HEP-style, but use NE instead of NE++ for the first phase.");

int main()
{

    Timer timer;
    timer.start();

    Partitioner *partitioner = NULL;

    partitioner = new LabelAssign("./../data/com-amazon.ungraph.txt");
//    partitioner = new LabelAssign("./../data/test.txt");
    partitioner->split();
    timer.stop();
    std::cout << "total time: " << timer.get_time();

}

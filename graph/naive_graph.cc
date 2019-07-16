#include "graph/naive_graph.h"
#include <set>

namespace hitnlls {
namespace graph {

void NaiveGraph::SymbolicAnalysis() {
    ::std::set<int> node_ids;
    ::std::set<int> factor_node_ids;
    for (auto iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
        if (!(iter->second)->GetFixed()) {
            node_ids.insert(iter->first);
        }
    }
    for (size_t idx = 0; idx < factors_.size(); ++idx) {
        for (int nidx = 0; nidx < factors_[idx]->GetNnodes(); ++nidx) {
            int nid = factors_[idx]->GetNodeId(nidx);
            factor_node_ids.insert(nid);
        }
    }
    ::std::set<int> effe_node_ids;
    for (auto iter = node_ids.begin(); iter != node_ids.end(); ++iter) {
        if (1 == factor_node_ids.count(*iter)) {
            effe_node_ids.insert(*iter);
        }
    }
    ordering_ = ::hitnlls::matrix::Matrixxi(effe_node_ids.size(), 1);
    int ridx = 0;
    for (auto iter = effe_node_ids.begin(); iter != effe_node_ids.end(); ++iter) {
        ordering_[ridx++] = *iter;
    }
}

HITNLLS_REGISTER_GRAPH(NaiveGraph);

} // namespace graph
} // namespace hitnlls
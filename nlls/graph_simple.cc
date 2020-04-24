#include "nlls/graph_simple.h"

namespace hitnlls {

void GraphSimple::SymbolicAnalysis() {
    std::set<int> node_ids;
    std::set<int> factor_node_ids;
    for (auto iter = id2v_.begin(); iter != id2v_.end(); ++iter) {
        if ((iter->second)->GetStatus() == VertexBase::Active) {
            node_ids.insert(iter->first);
        }
    }
    for (size_t idx = 0; idx < factors_.size(); ++idx) {
        for (int nidx = 0; nidx < factors_[idx]->NumVertices(); ++nidx) {
            int nid = factors_[idx]->GetVertexIdAt(nidx);
            factor_node_ids.insert(nid);
        }
    }
    std::set<int> effe_node_ids;
    for (auto iter = node_ids.begin(); iter != node_ids.end(); ++iter) {
        if (1 == factor_node_ids.count(*iter)) {
            effe_node_ids.insert(*iter);
        }
    }
    ordering_ = SparseMatrixi(effe_node_ids.size(), 1);
    int ridx = 0;
    for (auto iter = effe_node_ids.begin(); iter != effe_node_ids.end(); ++iter) {
        ordering_(ridx++, 0) = *iter;
    }
}

HITNLLS_REGISTER_GRAPH(GraphSimple)

} // namespace hitnlls
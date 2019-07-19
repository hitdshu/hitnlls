#include "graph/minheap_graph.h"
#include <set>
#include <map>
#include <algorithm>
#include <vector>

namespace hitnlls {
namespace graph {

namespace {

struct HeapElement {
    int nid;
    int idx;
    int ncn;
};

struct HeapCompare {
    bool operator()(const HeapElement &e1, const HeapElement &e2) {
        return e1.ncn > e2.ncn;
    }
};

}

void MinheapGraph::SymbolicAnalysis() {
    NaiveGraph::SymbolicAnalysis();
    ::std::set<int> node_ids;
    ::std::map<int, int> nid2idx;
    ::hitnlls::matrix::Matrixsi symb_smat(ordering_.Rows(), ordering_.Rows());
    for (size_t idx =0; idx < ordering_.Rows(); ++idx) {
        node_ids.insert(ordering_[idx]);
        nid2idx[ordering_[idx]] = idx;
    }
    for (size_t idx = 0; idx < factors_.size(); ++idx) {
        for (int nidx1 = 0; nidx1 < factors_[idx]->GetNnodes(); ++nidx1) {
            int nid1 = factors_[idx]->GetNodeId(nidx1);
            if (1 == node_ids.count(nid1)) {
                for (int nidx2 = nidx1; nidx2 < factors_[idx]->GetNnodes(); ++nidx2) {
                    int nid2 = factors_[idx]->GetNodeId(nidx2);
                    if (1 == node_ids.count(nid2)) {
                        symb_smat(nid2idx[nid1], nid2idx[nid2]) = 1;
                        symb_smat(nid2idx[nid2], nid2idx[nid1]) = 1;
                    }
                }
            }
        }
    }
    ::std::vector<HeapElement> reord_heap;
    for (int idx = 0; idx < symb_smat.Rows(); ++idx) {
        HeapElement tmp_ele;
        tmp_ele.nid = ordering_[idx];
        tmp_ele.idx = idx;
        tmp_ele.ncn = symb_smat.NumElementsInRow(idx) - 1;
        reord_heap.push_back(tmp_ele);
    }
    HeapCompare comp;
    ::std::make_heap(reord_heap.begin(), reord_heap.end(), comp);
    ::std::sort_heap(reord_heap.begin(), reord_heap.end(), comp);
    for (int idx = 0; idx < symb_smat.Rows(); ++idx) {
        ordering_[idx] = reord_heap[idx].nid;
    }
}

HITNLLS_REGISTER_GRAPH(MinheapGraph);

} // namespace graph
} // namespace hitnlls
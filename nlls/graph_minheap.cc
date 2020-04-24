#include <set>
#include <map>
#include <algorithm>
#include <vector>
#include <list>
#include <limits>

#include "nlls/graph_minheap.h"

namespace hitnlls {
namespace {
struct Element {
    int nid;
    int idx;
    int ncn;
    int oid;
};
}

void GraphMinheap::SymbolicAnalysis() {
    GraphSimple::SymbolicAnalysis();
    std::set<int> node_ids;
    std::map<int, int> nid2idx;
    SparseMatrixi symb_smat(ordering_.Rows(), ordering_.Rows());
    for (size_t idx =0; idx < ordering_.Rows(); ++idx) {
        node_ids.insert(ordering_(idx, 0));
        nid2idx[ordering_(idx, 0)] = idx;
    }
    for (size_t idx = 0; idx < factors_.size(); ++idx) {
        for (int nidx1 = 0; nidx1 < factors_[idx]->NumVertices(); ++nidx1) {
            int nid1 = factors_[idx]->GetVertexIdAt(nidx1);
            if (1 == node_ids.count(nid1)) {
                for (int nidx2 = nidx1; nidx2 < factors_[idx]->NumVertices(); ++nidx2) {
                    int nid2 = factors_[idx]->GetVertexIdAt(nidx2);
                    if (1 == node_ids.count(nid2)) {
                        symb_smat(nid2idx[nid1], nid2idx[nid2]) = 1;
                        symb_smat(nid2idx[nid2], nid2idx[nid1]) = 1;
                    }
                }
            }
        }
    }
    std::vector<Element> elements;
    std::map<int, std::list<Element *>>  ele_heap;
    for (int idx = 0; idx < symb_smat.Rows(); ++idx) {
        Element tmp_ele;
        tmp_ele.nid = ordering_(idx, 0);
        tmp_ele.ncn = symb_smat.NumElementsAtRow(idx) - 1;
        tmp_ele.oid = symb_smat.Rows() - 1;
        tmp_ele.idx = idx;
        elements.push_back(tmp_ele);
        ele_heap[tmp_ele.ncn].push_back(&tmp_ele);
    }
    for (int idx = 0; idx < symb_smat.Rows() - 1; ++idx) {
        Element *lst_ncn = nullptr;
        while (1) {
            lst_ncn = ele_heap.begin()->second.front();
            if (ele_heap.begin()->first != lst_ncn->ncn) {
                ele_heap.begin()->second.pop_front();
                ele_heap[lst_ncn->ncn].push_back(lst_ncn);
                if (0 == ele_heap.begin()->second.size()) {
                    ele_heap.erase(ele_heap.begin());
                }
            } else {
                ele_heap.begin()->second.pop_front();
                if (0 == ele_heap.begin()->second.size()) {
                    ele_heap.erase(ele_heap.begin());
                }
            }
        }
        lst_ncn->oid = idx;
        std::vector<int> nn = symb_smat.IndicesAtRow(lst_ncn->idx);
        for (auto nn_iter = nn.begin(); nn_iter != nn.end(); ++nn_iter) {
            if (*nn_iter == lst_ncn->nid) {
                continue;
            }
            --(elements[*nn_iter].ncn);
            for (auto nn_iter_n = ++nn_iter; nn_iter_n != nn.end(); ++nn_iter_n) {
                if (*nn_iter_n == lst_ncn->nid) {
                    continue;
                }
                if (!symb_smat.CountIndex(*nn_iter, *nn_iter_n)) {
                    symb_smat(*nn_iter, *nn_iter_n) = 1;
                    symb_smat(*nn_iter_n, *nn_iter) = 1;
                    ++(elements[*nn_iter].ncn);
                    ++(elements[*nn_iter_n].ncn);
                }
            }
        }
    }
    for (auto iter = elements.begin(); iter != elements.end(); ++iter) {
        ordering_(iter->oid, 0) = iter->nid;
    }
    return;
}

HITNLLS_REGISTER_GRAPH(GraphMinheap)

} // namespace hitnlls
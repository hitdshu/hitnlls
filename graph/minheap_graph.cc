#include "graph/minheap_graph.h"
#include <set>

namespace hitnlls {
namespace graph {

void MinheapGraph::SymbolicAnalysis() {
    NaiveGraph::SymbolicAnalysis();
    
}

HITNLLS_REGISTER_GRAPH(MinheapGraph);

} // namespace graph
} // namespace hitnlls
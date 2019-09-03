//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_GOC_H
#define GOC_GOC_H

#include "goc/base/maybe.h"

#include "goc/collection/bitset_utils.h"
#include "goc/collection/collection_utils.h"
#include "goc/collection/matrix.h"
#include "goc/collection/vector_map.h"

#include "goc/exception/exception_utils.h"

#include "goc/graph/arc.h"
#include "goc/graph/digraph.h"
#include "goc/graph/edge.h"
#include "goc/graph/graph.h"
#include "goc/graph/graph_path.h"
#include "goc/graph/maxflow_mincut.h"
#include "goc/graph/path_finding.h"
#include "goc/graph/vertex.h"

#include "goc/json/json_utils.h"

#include "goc/lib/json.hpp"

#include "goc/linear_programming/cuts/separation_routine.h"
#include "goc/linear_programming/cuts/separation_strategy.h"
#include "goc/linear_programming/model/branch_priority.h"
#include "goc/linear_programming/model/constraint.h"
#include "goc/linear_programming/model/expression.h"
#include "goc/linear_programming/model/formulation.h"
#include "goc/linear_programming/model/valuation.h"
#include "goc/linear_programming/model/variable.h"
#include "goc/linear_programming/solver/bc_solver.h"
#include "goc/linear_programming/solver/cg_solver.h"
#include "goc/linear_programming/solver/lp_solver.h"

#include "goc/log/bcp_execution_log.h"
#include "goc/log/blb_execution_log.h"
#include "goc/log/cg_execution_log.h"
#include "goc/log/log.h"
#include "goc/log/lp_execution_log.h"
#include "goc/log/mlb_execution_log.h"

#include "goc/math/interval.h"
#include "goc/math/linear_function.h"
#include "goc/math/number_utils.h"
#include "goc/math/point_2d.h"
#include "goc/math/pwl_function.h"

#include "goc/print/print_utils.h"
#include "goc/print/printable.h"
#include "goc/print/table_stream.h"

#include "goc/runner/runner_utils.h"

#include "goc/string/string_utils.h"

#include "goc/time/date.h"
#include "goc/time/duration.h"
#include "goc/time/point_in_time.h"
#include "goc/time/stopwatch.h"
#include "goc/time/watch.h"

#include "goc/vrp/route.h"
#include "goc/vrp/vrp_solution.h"

#endif // GOC_GOC_H
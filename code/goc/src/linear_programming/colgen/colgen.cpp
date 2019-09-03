//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/colgen/colgen.h"

#include "goc/lib/json.hpp"
#include "goc/time/duration.h"
#include "goc/time/stopwatch.h"
#include "goc/string/string_utils.h"
#include "goc/math/number_utils.h"
#include "goc/print/table_stream.h"

using namespace std;
using namespace nlohmann;

namespace goc
{
namespace
{
unordered_map<LPStatus, CGStatus> mapper = {{LPStatus::DidNotStart, CGStatus::DidNotStart},
											{LPStatus::Infeasible, CGStatus::Infeasible},
											{LPStatus::Unbounded, CGStatus::Unbounded},
											{LPStatus::TimeLimitReached, CGStatus::TimeLimitReached}, {LPStatus::MemoryLimitReached, CGStatus::MemoryLimitReached},
											{LPStatus::Optimum, CGStatus::Optimum}};
CGStatus parse_lp_status(LPStatus status)
{
	return mapper[status];
}
}

CGExecutionLog solve_colgen(Formulation* formulation,
				   ostream* screen_output,
				   Duration time_limit,
				   const PricingFunction& pricing_function,
				   LPSolver* lp_solver,
				   const unordered_set<CGOption>& option)
{
	Stopwatch rolex(true);
	
	// Keep track of the execution.
	CGExecutionLog execution_log;
	if (!execution_log.lp_time.IsSet()) execution_log.lp_time = 0.0_sec;
	if (!execution_log.pricing_time.IsSet()) execution_log.pricing_time = 0.0_sec;
	if (!execution_log.iterations.IsSet()) execution_log.iterations = vector<json>{};
	if (!execution_log.iteration_count.IsSet()) execution_log.iteration_count = 0;
	
	execution_log.status = CGStatus::Optimum;
	
	TableStream output(screen_output, 1.0);
	output.AddColumn("time", 10).AddColumn("#", 5).AddColumn("value", 10).AddColumn("#cols", 7);
	
	// While the pricing problem finds new columns, keep iterating.
	int variable_count = -1, initial_variable_count = formulation->VariableCount();
	int row_count = -1;
	output.WriteHeader();
	double objective_value = 0.0;
	while (variable_count < formulation->VariableCount() || row_count < formulation->ConstraintCount())
	{
		// Check if time limit was exceeded.
		if (rolex.Peek() >= time_limit) {execution_log.status = CGStatus::TimeLimitReached; break; }
		
		// Solve LP relaxation to get dual variables.
		lp_solver->time_limit = time_limit - rolex.Peek();
		auto lp_log = lp_solver->Solve(formulation, {LPOption::Duals, LPOption::Incumbent});
		*execution_log.lp_time += *lp_log.time;
		
		if (*lp_log.status != LPStatus::Optimum) { execution_log.status = parse_lp_status(*lp_log.status); break; }
		objective_value = lp_log.incumbent_value;
		if (output.RegisterAttempt()) output.WriteRow({STR(rolex.Peek()), STR(execution_log.iteration_count++), STR(objective_value), STR(formulation->VariableCount())});
		
		// Update variable count before solving the pricing problem.
		variable_count = formulation->VariableCount();
		row_count = formulation->ConstraintCount();
		
		// Solve the pricing problem (i.e. add new variables to the formulation).
		Stopwatch pricing_rolex(true);
		pricing_function(*lp_log.duals, *lp_log.incumbent_value, time_limit - rolex.Peek(), &execution_log);
		*execution_log.pricing_time += pricing_rolex.Pause();
	}
	output.WriteRow({STR(rolex.Peek()), STR(execution_log.iteration_count), STR(objective_value), STR(formulation->VariableCount())});
	if (screen_output) *screen_output << endl;
	
	// If the column generation was solved to optimality, get the actual solution.
	if (*execution_log.status == CGStatus::Optimum)
	{
		lp_solver->time_limit = Duration::Max();
		auto lp_log = lp_solver->Solve(formulation, {LPOption::Incumbent});
		execution_log.incumbent_value = lp_log.incumbent_value;
		execution_log.incumbent = lp_log.incumbent;
	}
	execution_log.columns_added = formulation->VariableCount() - initial_variable_count;
	execution_log.time = rolex.Peek();
	
	return execution_log;
}
} // namespace goc
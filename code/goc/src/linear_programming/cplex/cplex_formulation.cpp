//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/cplex/cplex_formulation.h"

#include <map>

#include "goc/collection/collection_utils.h"
#include "goc/exception/exception_utils.h"
#include "goc/math/number_utils.h"

using namespace std;

namespace goc
{
namespace
{
// The internal structure of CPLEX to store a constraint.
struct CplexRow
{
	int nzcnt, nzind;
	double rhs;
	char sense;
	vector<int> rmatbeg, rmatind;
	vector<double> rmatval;
};

// Returns: a CPLEX row representing the constraint.
// Precondition: constraint must be normalized.
CplexRow constraint_to_cplex_row(const Constraint& constraint)
{
	CplexRow row;
	
	map<enum Constraint::Sense, char> cplex_senses = {{Constraint::LessEqual, 'L'}, {Constraint::GreaterEqual, 'G'},
													  {Constraint::Equality, 'E'}};
	row.nzcnt = constraint.LeftSide().NonZeroVariableTermCount();
	row.rhs = constraint.RightSide();
	row.sense = cplex_senses[constraint.Sense()];
	row.rmatbeg.push_back(0);
	row.nzind = 0;
	for (auto& term: constraint.LeftSide().Terms())
	{
		row.rmatind.push_back(term.first.Index());
		row.rmatval.push_back(term.second);
		row.nzind++;
	}
	return row;
}
}

CplexFormulation::CplexFormulation()
{
	env_memory_handler_ = shared_ptr<cpxenv>(cplex::openCPLEX(), [&] (CPXENVptr env_p)
	{
		cplex::closeCPLEX(&env_p);
	});
	env_ = env_memory_handler_.get();
	problem_ = cplex::createprob(env_, "formulation");
}

CplexFormulation::~CplexFormulation()
{
	cplex::freeprob(env_, &problem_);
	for (int* i: variable_indices_) delete i;
}

int CplexFormulation::AddConstraint(const Constraint& constraint)
{
	// Add constraint to CPLEX formulation.
	auto cplex_row = constraint_to_cplex_row(constraint);
	cplex::addrows(env_, problem_, 0, 1, cplex_row.nzind, &cplex_row.rhs, &cplex_row.sense, &cplex_row.rmatbeg[0], &cplex_row.rmatind[0], &cplex_row.rmatval[0], nullptr, nullptr);
	
	return ConstraintCount()-1;
}

void CplexFormulation::RemoveConstraint(int constraint_index)
{
	// Remove constraint from CPLEX.
	cplex::delrows(env_, problem_, constraint_index, constraint_index);
	
	// Reduce all constraints indices following 'constraint_id' by one.
	for (int i = constraint_index; i < constraint_indices_.size()-1; ++i)
	{
		swap(constraint_indices_[i], constraint_indices_[i+1]);
		*constraint_indices_[i] = i;
	}
	delete constraint_indices_.back();
	constraint_indices_.pop_back();
}

void CplexFormulation::AddLazyConstraint(SeparationRoutine* lazy_constraint)
{
	if (!lazy_constraint) return;
	lazy_constraints_.push_back(lazy_constraint);
}

void CplexFormulation::RemoveLazyConstraint(SeparationRoutine* lazy_constraint)
{
	if (!lazy_constraint) return;
	lazy_constraints_.erase(remove(lazy_constraints_.begin(), lazy_constraints_.end(), lazy_constraint));
}

Variable CplexFormulation::AddVariable(const string& name, VariableDomain domain, double lower_bound, double upper_bound)
{
	// Add variable to internal structure.
	variable_indices_.push_back(new int(variable_indices_.size()));
	variable_names_.push_back(name);
	
	// Add variable to CPLEX.
	char* colname[] = {(char*)variable_names_.back().c_str()};
	cplex::newcols(env_, problem_, 1, nullptr, nullptr, nullptr, nullptr, colname);
	
	Variable var(name, variable_indices_.back());
	
	// Set variable domain.
	SetVariableDomain(var, domain);
	
	// Set bounds.
	SetVariableBound(var, lower_bound, upper_bound);
	
	return var;
}

void CplexFormulation::RemoveVariable(const Variable& variable)
{
	// Remove variable from CPLEX.
	cplex::delcols(env_, problem_, variable.Index(), variable.Index());
	
	// Reduce all variable indices following the erased variable by one and delete its name from the vector of names.
	for (int i = variable.Index(); i < variable_indices_.size()-1; ++i)
	{
		swap(variable_names_[i], variable_names_[i+1]);
		swap(variable_indices_[i], variable_indices_[i+1]);
		*variable_indices_[i] = i;
	}
	variable_names_.pop_back();
	delete variable_indices_.back();
	variable_indices_.pop_back();
}

void CplexFormulation::SetVariableDomain(const Variable& variable, VariableDomain domain)
{
	std::map<VariableDomain, char> cplex_domains = {{VariableDomain::Real, 'C'}, {VariableDomain::Integer, 'I'},
													{VariableDomain::Binary, 'B'}};
	int indices[] = {variable.Index()};
	char xctype[] = {cplex_domains[domain]};
	cplex::chgctype(env_, problem_, 1, indices, xctype);
}

void CplexFormulation::SetVariableBound(const Variable& v, double lower_bound, double upper_bound)
{
	if (upper_bound == INFTY) upper_bound = CPX_INFBOUND;
	if (lower_bound == -INFTY) lower_bound = -CPX_INFBOUND;
	int indices[] = {v.Index(), v.Index()};
	double bd[] = {lower_bound, upper_bound};
	char type[] = {'L', 'U'};
	cplex::chgbds(env_, problem_, 2, indices, type, bd);
}

void CplexFormulation::SetVariableLowerBound(const Variable& v, double lower_bound)
{
	if (lower_bound == -INFTY) lower_bound = -CPX_INFBOUND;
	int indices[] = {v.Index()};
	double bd[] = {lower_bound};
	char type[] = {'L'};
	cplex::chgbds(env_, problem_, 1, indices, type, bd);
}

void CplexFormulation::SetVariableUpperBound(const Variable& v, double upper_bound)
{
	if (upper_bound == INFTY) upper_bound = CPX_INFBOUND;
	int indices[] = {v.Index()};
	double bd[] = {upper_bound};
	char type[] = {'U'};
	cplex::chgbds(env_, problem_, 1, indices, type, bd);
}

void CplexFormulation::Minimize(const Expression& objective_function)
{
	vector<int> indices = range(0, VariableCount());
	vector<double> values(VariableCount(), 0.0);
	int i = 0;
	for (auto& term: objective_function.Terms())
	{
		values[term.first.Index()] = term.second;
		++i;
	}
	cplex::chgobj(env_, problem_, VariableCount(), &indices[0], &values[0]);
	cplex::chgobjsen(env_, problem_, CPX_MIN);
}

void CplexFormulation::Maximize(const Expression& objective_function)
{
	vector<int> indices = range(0, VariableCount());
	vector<double> values(VariableCount(), 0.0);
	int i = 0;
	for (auto& term: objective_function.Terms())
	{
		values[term.first.Index()] = term.second;
		++i;
	}
	cplex::chgobj(env_, problem_, VariableCount(), &indices[0], &values[0]);
	cplex::chgobjsen(env_, problem_, CPX_MAX);
}

void CplexFormulation::SetConstraintRightHandSide(int constraint_index, double value)
{
	cplex::chgrhs(env_, problem_, 1, &constraint_index, &value);
}

void CplexFormulation::SetConstraintCoefficient(int constraint_index, const Variable& variable, double coefficient)
{
	cplex::chgcoef(env_, problem_, constraint_index, variable.Index(), coefficient);
}

void CplexFormulation::SetObjectiveCoefficient(const Variable& variable, double coefficient)
{
	int indices[] = {variable.Index()};
	cplex::chgobj(env_, problem_, 1, indices, &coefficient);
}

Formulation::ObjectiveSense CplexFormulation::GetObjectiveSense() const
{
	return cplex::getobjsen(env_, problem_) == 1 ? ObjectiveSense::Minimization : ObjectiveSense::Maximization;
}

double CplexFormulation::GetObjectiveCoefficient(const Variable& variable) const
{
	double values[] = {0};
	cplex::getobj(env_, problem_, values, variable.Index(), variable.Index());
	return values[0];
}

double CplexFormulation::GetConstraintRightHandSide(int constraint_index) const
{
	double rhs;
	cplex::getrhs(env_, problem_, &rhs, constraint_index, constraint_index);
	return rhs;
}

double CplexFormulation::GetConstraintCoefficient(int constraint_index, const Variable& variable)
{
	double coef;
	cplex::getcoef(env_, problem_, constraint_index, variable.Index(), &coef);
	return coef;
}

VariableDomain CplexFormulation::GetVariableDomain(const Variable& variable) const
{
	char type;
	cplex::getctype(env_, problem_, &type, variable.Index(), variable.Index());
	
	std::map<char, VariableDomain> cplex_domains = {{'C', VariableDomain::Real}, {'I', VariableDomain::Integer},
													{'B', VariableDomain::Binary}};
	return cplex_domains[type];
}

pair<double, double> CplexFormulation::GetVariableBound(const Variable& variable) const
{
	double lb, ub;
	cplex::getlb(env_, problem_, &lb, variable.Index(), variable.Index());
	cplex::getub(env_, problem_, &ub, variable.Index(), variable.Index());
	if (lb == -CPX_INFBOUND) lb = -INFTY;
	if (ub == CPX_INFBOUND) ub = INFTY;
	return {lb, ub};
}

Expression CplexFormulation::ObjectiveFunction() const
{
	vector<double> coefficients(VariableCount());
	cplex::getobj(env_, problem_, &coefficients[0], 0, VariableCount()-1);
	Expression obj;
	for (int i = 0; i < VariableCount(); ++i) obj += coefficients[i] * Variable(variable_names_[i], variable_indices_[i]);
	return obj;
}

vector<Variable> CplexFormulation::Variables() const
{
	vector<Variable> variables;
	for (int i = 0; i < VariableCount(); ++i) variables.push_back({variable_names_[i], variable_indices_[i]});
	return variables;
}

vector<Constraint> CplexFormulation::Constraints() const
{
	vector<Constraint> constraints;
	int variable_count = VariableCount();
	
	// Set buffer sizes.
	CplexRow row;
	row.rmatbeg.assign(1,0);
	row.rmatind.assign(variable_count, 0);
	row.rmatval.assign(variable_count, 0.0);
	int surplus = 0;
	for (int i = 0; i < ConstraintCount(); ++i)
	{
		// Get i-th row from CPLEX.
		cplex::getrows(env_, problem_, &row.nzcnt, &row.rmatbeg[0], &row.rmatind[0], &row.rmatval[0], variable_count, &surplus, i, i);
		
		// Create constraint with row values.
		Expression left;
		for (int j = 0; j < row.nzcnt; ++j) left += row.rmatval[j] * VariableAtIndex(row.rmatind[j]);
		double right = GetConstraintRightHandSide(i);
		char sense;
		cplex::getsense(env_, problem_, &sense, i, i);
		if (sense == 'L') constraints.push_back(left.LEQ(right));
		else if (sense == 'G') constraints.push_back(left.GEQ(right));
		else if (sense == 'E') constraints.push_back(left.EQ(right));
	}
	return constraints;
}

const vector<SeparationRoutine*>& CplexFormulation::LazyConstraints() const
{
	return lazy_constraints_;
}

int CplexFormulation::VariableCount() const
{
	return cplex::getnumcols(env_, problem_);
}

int CplexFormulation::ConstraintCount() const
{
	return cplex::getnumrows(env_, problem_);
}

Variable CplexFormulation::VariableAtIndex(int variable_index) const
{
	return Variable(variable_names_[variable_index], variable_indices_[variable_index]);
}

double CplexFormulation::EvaluateValuation(const Valuation& valuation) const
{
	return ObjectiveFunction().Value(valuation);
}

bool CplexFormulation::IsFeasibleValuation(const Valuation& v, bool verbose) const
{
	// Check if all constraints hold for valuation v.
	auto constraints = Constraints();
	for (auto& c: constraints)
	{
		if (!c.Holds(v))
		{
			if (verbose) clog << c << endl;
			return false;
		}
	}
	return true;
}

Formulation* CplexFormulation::Copy() const
{
	CplexFormulation* copy = new CplexFormulation(env_memory_handler_, cplex::cloneprob(env_, problem_));
	copy->variable_names_ = variable_names_;
	for (int i = 0; i < VariableCount(); ++i) copy->variable_indices_.push_back(new int(i));
	for (int i = 0; i < ConstraintCount(); ++i) copy->constraint_indices_.push_back(new int(i));
	copy->lazy_constraints_ = lazy_constraints_;
	return copy;
}

void CplexFormulation::Print(ostream& os) const
{
	os << ObjectiveFunction() << endl;
	os << "s.t.";
	for (auto& c: Constraints()) os << endl << c;
	map<VariableDomain, string> domain_to_str = {{VariableDomain::Real, "R"}, {VariableDomain::Binary, "{0,1}"}, {VariableDomain::Integer, "Z"}};
	for (auto& v: Variables()) os << endl << GetVariableBound(v).first << " <= " << v << " <= " << GetVariableBound(v).second;
	for (auto& v: Variables()) os << endl << v << " \\in " << domain_to_str[GetVariableDomain(v)];
}

CPXENVptr CplexFormulation::Environment() const
{
	return env_;
}

CPXLPptr CplexFormulation::Problem() const
{
	return problem_;
}

CplexFormulation::CplexFormulation(const std::shared_ptr<cpxenv>& env_memory_handler, CPXLPptr problem)
	: env_memory_handler_(env_memory_handler), env_(env_memory_handler.get()), problem_(problem)
{

}
} // namespace goc
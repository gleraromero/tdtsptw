//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#include "goc/linear_programming/cplex/cplex_wrapper.h"

#include <string>

#include "goc/exception/exception_utils.h"
#include "goc/string/string_utils.h"

using namespace std;

namespace goc
{
namespace cplex
{
namespace
{
void fail_with_error_message(CPXENVptr env, int status, const string& function_name)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	CPXgeterrorstring(env, status, errmsg);
	fail("Failed at function " + function_name + " with status: " + STR(status) + ", with message: " + string(errmsg));
}

void fail_with_error_message(CPXCENVptr env, int status, const string& function_name)
{
	char errmsg[CPXMESSAGEBUFSIZE];
	CPXgeterrorstring(env, status, errmsg);
	fail("Failed at function " + function_name + " with status: " + STR(status) + ", with message: " + string(errmsg));
}

void fail_with_error_message(int status, const string& function_name)
{
	fail("Failed at function " + function_name + " with status: " + STR(status));
}
}

CPXENVptr openCPLEX()
{
	int status = 0;
	CPXENVptr env = CPXopenCPLEX(&status);
	if (status == 0)
	{
		return env;
	}
	else
	{
		fail_with_error_message(env, status, "CPXopenCPLEX");
		return nullptr;
	}
}

CPXLPptr createprob(CPXENVptr env, char const* probname_str)
{
	int status = 0;
	CPXLPptr problem = CPXcreateprob(env, &status, probname_str);
	if (status == 0)
	{
		return problem;
	}
	else
	{
		fail_with_error_message(env, status, "CPXcreateprob");
		return nullptr;
	}
}

void closeCPLEX(CPXENVptr* env_p)
{
	int status = CPXcloseCPLEX(env_p);
	if (status != 0)
	{
		fail_with_error_message(*env_p, status, "CPXcloseCPLEX");
	}
}

void freeprob(CPXENVptr env, CPXLPptr* lp_p)
{
	int status = CPXfreeprob(env, lp_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXfreeprob");
	}
}

void
newcols(CPXENVptr env, CPXLPptr lp, int ccnt, double const* obj, double const* lb, double const* ub, char const* xctype,
		char** colname)
{
	int status = CPXnewcols(env, lp, ccnt, obj, lb, ub, xctype, colname);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXnewcols");
	}
}

void delcols(CPXCENVptr env, CPXLPptr lp, int begin, int end)
{
	int status = CPXdelcols(env, lp, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXdelcols");
	}
}

void chgctype(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, char const* xctype)
{
	int status = CPXchgctype(env, lp, cnt, indices, xctype);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgctype");
	}
}

void getctype(CPXENVptr env, CPXCLPptr lp, char* xctype, int begin, int end)
{
	int status = CPXgetctype(env, lp, xctype, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetctype");
	}
}

void addrows(CPXENVptr env, CPXLPptr lp, int ccnt, int rcnt, int nzcnt, double const* rhs, char const* sense,
			 int const* rmatbeg, int const* rmatind, double const* rmatval, char** colname, char** rowname)
{
	int status = CPXaddrows(env, lp, ccnt, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, colname, rowname);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXaddrows");
	}
}

void getlb(CPXENVptr env, CPXCLPptr lp, double* lb, int begin, int end)
{
	int status = CPXgetlb(env, lp, lb, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetlb");
	}
}

void getub(CPXENVptr env, CPXCLPptr lp, double* ub, int begin, int end)
{
	int status = CPXgetub(env, lp, ub, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetub");
	}
}

void getpi(CPXENVptr env, CPXCLPptr lp, double* pi, int begin, int end)
{
	int status = CPXgetpi(env, lp, pi, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetpi");
	}
}

void getdj(CPXCENVptr env, CPXCLPptr lp, double* dj, int begin, int end)
{
	int status = CPXgetdj(env, lp, dj, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetdj");
	}
}

void chgbds(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, char const* lu, double const* bd)
{
	int status = CPXchgbds(env, lp, cnt, indices, lu, bd);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgbds");
	}
}

void chgobj(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, double const* values)
{
	int status = CPXchgobj(env, lp, cnt, indices, values);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgobj");
	}
}

int getobjsen(CPXCENVptr env, CPXCLPptr lp)
{
	return CPXgetobjsen(env, lp);
}

void getobj(CPXENVptr env, CPXLPptr lp, double* obj, int begin, int end)
{
	int status = CPXgetobj(env, lp, obj, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetobj");
	}
}

void getrhs(CPXCENVptr env, CPXCLPptr lp, double* rhs, int begin, int end)
{
	int status = CPXgetrhs(env, lp, rhs, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetrhs");
	}
}

void chgrhs(CPXCENVptr env, CPXLPptr lp, int cnt, int const* indices, double const* values)
{
	int status = CPXchgrhs(env, lp, cnt, indices, values);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgrhs");
	}
}

void getcoef(CPXCENVptr env, CPXCLPptr lp, int i, int j, double* coef_p)
{
	int status = CPXgetcoef(env, lp, i, j, coef_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetcoef");
	}
}

void chgcoef(CPXENVptr env, CPXLPptr lp, int i, int j, double newvalue)
{
	int status = CPXchgcoef(env, lp, i, j, newvalue);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgcoef");
	}
}

void chgobjsen(CPXENVptr env, CPXLPptr lp, int maxormin)
{
	int status = CPXchgobjsen(env, lp, maxormin);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgobjsen");
	}
}

void delrows(CPXENVptr env, CPXLPptr lp, int begin, int end)
{
	int status = CPXdelrows(env, lp, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXdelrows");
	}
}

CPXLPptr cloneprob(CPXENVptr env, CPXCLPptr lp)
{
	int status;
	CPXLPptr problem = CPXcloneprob(env, lp, &status);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXcloneprob");
	}
	return problem;
}

void mipopt(CPXENVptr env, CPXLPptr lp)
{
	int status = CPXmipopt(env, lp);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXmipopt");
	}
}

void lpopt(CPXENVptr env, CPXLPptr lp)
{
	int status = CPXlpopt(env, lp);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXlpopt");
	}
}

void
solution(CPXENVptr env, CPXCLPptr lp, int* lpstat_p, double* objval_p, double* x, double* pi, double* slack, double* dj)
{
	int status = CPXsolution(env, lp, lpstat_p, objval_p, x, pi, slack, dj);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsolution");
	}
}

void getbestobjval(CPXENVptr env, CPXCLPptr lp, double* objval_p)
{
	int status = CPXgetbestobjval(env, lp, objval_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetbestobjval");
	}
}

int getnodecnt(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetnodecnt(env, lp);
}

int getnodeleftcnt(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetnodeleftcnt(env, lp);
}

bool getobjval(CPXENVptr env, CPXCLPptr lp, double* objval_p)
{
	return CPXgetobjval(env, lp, objval_p) == 0;
}

void getmiprelgap(CPXENVptr env, CPXCLPptr lp, double* gap_p)
{
	int status = CPXgetmiprelgap(env, lp, gap_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetmiprelgap");
	}
}

int getmipitcnt(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetmipitcnt(env, lp);
}

int getitcnt(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetitcnt(env, lp);
}

int getstat(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetstat(env, lp);
}

void getsense(CPXCENVptr env, CPXCLPptr lp, char* sense, int begin, int end)
{
	int status = CPXgetsense(env, lp, sense, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetsense");
	}
}

void getrows(CPXCENVptr env, CPXCLPptr lp, int* nzcnt_p, int* rmatbeg, int* rmatind, double* rmatval, int rmatspace,
			 int* surplus_p, int begin, int end)
{
	int status = CPXgetrows(env, lp, nzcnt_p, rmatbeg, rmatind, rmatval, rmatspace, surplus_p, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetrows");
	}
}

void setintparam(CPXENVptr env, int whichparam, CPXINT newvalue)
{
	int status = CPXsetintparam(env, whichparam, newvalue);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetintparam");
	}
}

void setdblparam(CPXENVptr env, int whichparam, double newvalue)
{
	int status = CPXsetdblparam(env, whichparam, newvalue);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetdblparam");
	}
}

void setlongparam(CPXENVptr env, int whichparam, long newvalue)
{
	int status = CPXsetlongparam(env, whichparam, newvalue);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetlongparam");
	}
}

void setstrparam(CPXENVptr env, int whichparam, char const* newvalue)
{
	int status = CPXsetstrparam(env, whichparam, newvalue);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetstrparam");
	}
}

void getintparam(CPXENVptr env, int whichparam, CPXINT* value_p)
{
	int status = CPXgetintparam(env, whichparam, value_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetintparam");
	}
}

void getlongparam(CPXENVptr env, int whichparam, CPXLONG* value_p)
{
	int status = CPXgetlongparam(env, whichparam, value_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetlongparam");
	}
}

void getdblparam(CPXENVptr env, int whichparam, double* value_p)
{
	int status = CPXgetdblparam(env, whichparam, value_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetdblparam");
	}
}

void chgprobtype(CPXENVptr env, CPXLPptr lp, int type)
{
	int status = CPXchgprobtype(env, lp, type);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXchgprobtype");
	}
}

int getprobtype(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetprobtype(env, lp);
}

void getparamname(CPXENVptr env, int whichparam, char* name_str)
{
	int status = CPXgetparamname(env, whichparam, name_str);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetparamname");
	}
}

void getparamtype(CPXENVptr env, int whichparam, int* paramtype)
{
	int status = CPXgetparamtype(env, whichparam, paramtype);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetparamtype");
	}
}

void getparamnum(CPXENVptr env, char const* name_str, int* whichparam_p)
{
	int status = CPXgetparamnum(env, name_str, whichparam_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetparamnum");
	}
}

void setlazyconstraintcallbackfunc(CPXENVptr env, int(CPXPUBLIC* lazyconcallback)(CALLBACK_CUT_ARGS),
												 void* cbhandle)
{
	int status = CPXsetlazyconstraintcallbackfunc(env, lazyconcallback, cbhandle);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetlazyconstraintcallbackfunc");
	}
}

void setusercutcallbackfunc(CPXENVptr env, int(CPXPUBLIC* cutcallback)(CALLBACK_CUT_ARGS),
										  void* cbhandle)
{
	int status = CPXsetusercutcallbackfunc(env, cutcallback, cbhandle);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetusercutcallbackfunc");
	}
}

void setbranchcallbackfunc(CPXENVptr env, int(CPXPUBLIC* branchcallback)(CALLBACK_BRANCH_ARGS),
										 void* cbhandle)
{
	int status = CPXsetbranchcallbackfunc(env, branchcallback, cbhandle);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetbranchcallbackfunc");
	}
}

void branchcallbackbranchbds(CPXENVptr env, void* cbdata, int wherefrom, int cnt, int const* indices,
										   char const* lu, double const* bd, double nodeest, void* userhandle,
										   int* seqnum_p)
{
	int status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, cnt, indices, lu, bd, nodeest,
											userhandle, seqnum_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXbranchcallbackbranchbds");
	}
}

void branchcallbackbranchconstraints(CPXENVptr env, void* cbdata, int wherefrom, int rcnt, int nzcnt,
												   double const* rhs, char const* sense, int const* rmatbeg,
												   int const* rmatind, double const* rmatval, double nodeest, void* userhandle,
												   int* seqnum_p)
{
	int status = CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, rcnt, nzcnt, rhs, sense,
													rmatbeg, rmatind, rmatval, nodeest, userhandle, seqnum_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXbranchcallbackbranchconstraints");
	}
}

void cutcallbackadd(CPXENVptr env, void* cbdata, int wherefrom, int nzcnt, double rhs, int sense,
								  int const* cutind, double const* cutval, int purgeable)
{
	int status = CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs, sense, cutind, cutval,
								   purgeable);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXcutcallbackadd");
	}
}

void getcallbacknodex(CPXENVptr env, void* cbdata, int wherefrom, double* x, int begin, int end)
{
	int status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetcallbacknodex");
	}
}

void getcallbackinfo(CPXCENVptr env, void* cbdata, int wherefrom, int whichinfo, void* result_p)
{
	int status = CPXgetcallbackinfo(env, cbdata, wherefrom, whichinfo, result_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetcallbackinfo");
	}
}

void addmipstarts(CPXENVptr env, CPXLPptr lp, int mcnt, int nzcnt, int const* beg, int const* varindices,
								double const* values, int const* effortlevel, char** mipstartname)
{
	int status = CPXaddmipstarts(env, lp, mcnt, nzcnt, beg, varindices, values, effortlevel,
								 mipstartname);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXaddmipstarts");
	}
}

int getnummipstarts(CPXCENVptr env, CPXCLPptr lp)
{
	return CPXgetnummipstarts(env, lp);
}

void delmipstarts(CPXCENVptr env, CPXLPptr lp, int begin, int end)
{
	int status = CPXdelmipstarts(env, lp, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXdelmipstarts");
	}
}

void getnumcuts(CPXENVptr env, CPXCLPptr lp, int cuttype, int* num_p)
{
	int status = CPXgetnumcuts(env, lp, cuttype, num_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetnumcuts");
	}
}

int getnumrows(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetnumrows(env, lp);
}

int getnumcols(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetnumcols(env, lp);
}

void copyorder(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, int const* priority,
							 int const* direction)
{
	int status = CPXcopyorder(env, lp, cnt, indices, priority, direction);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXcopyorder");
	}
}

void populate(CPXENVptr env, CPXLPptr lp)
{
	int status = CPXpopulate(env, lp);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXpopulate");
	}
}

int getsolnpoolnumsolns(CPXENVptr env, CPXCLPptr lp)
{
	return CPXgetsolnpoolnumsolns(env, lp);
}

void getsolnpoolx(CPXENVptr env, CPXCLPptr lp, int soln, double* x, int begin, int end)
{
	int status = CPXgetsolnpoolx(env, lp, soln, x, begin, end);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetsolnpoolx");
	}
}

void setdefaults(CPXENVptr env)
{
	int status = CPXsetdefaults(env);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetdefaults");
	}
}

void addfuncdest(CPXCENVptr env, CPXCHANNELptr channel, void* handle,
							   void(CPXPUBLIC* msgfunction)(void*, const char*))
{
	int status = CPXaddfuncdest(env, channel, handle, msgfunction);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXaddfuncdest");
	}
}

void delfuncdest(CPXCENVptr env, CPXCHANNELptr channel, void* handle,
							   void(CPXPUBLIC* msgfunction)(void*, const char*))
{
	int status = CPXdelfuncdest(env, channel, handle, msgfunction);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXdelfuncdest");
	}
}

void getchannels(CPXCENVptr env, CPXCHANNELptr* cpxresults_p, CPXCHANNELptr* cpxwarning_p,
							   CPXCHANNELptr* cpxerror_p, CPXCHANNELptr* cpxlog_p)
{
	int status = CPXgetchannels(env, cpxresults_p, cpxwarning_p, cpxerror_p, cpxlog_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetchannels");
	}
}

void setinfocallbackfunc(CPXENVptr env,
									   int(CPXPUBLIC* callback)(CPXCENVptr, void*, int, void*),
									   void* cbhandle)
{
	int status = CPXsetinfocallbackfunc(env, callback, cbhandle);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXsetinfocallbackfunc");
	}
}

void getcallbacknodeinfo(CPXCENVptr env, void* cbdata, int wherefrom, int nodeindex, int whichinfo,
									   void* result_p)
{
	int status = CPXgetcallbacknodeinfo(env, cbdata, wherefrom, nodeindex, whichinfo, result_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetcallbacknodeinfo");
	}
}

void getcallbacknodeobjval(CPXCENVptr env, void* cbdata, int wherefrom, double* objval_p)
{
	int status = CPXgetcallbacknodeobjval(env, cbdata, wherefrom, objval_p);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetcallbacknodeobjval");
	}
}

void callbacksetfunc(CPXENVptr env, CPXLPptr lp, CPXLONG contextmask, CPXCALLBACKFUNC callback,
								   void* userhandle)
{
	int status = CPXcallbacksetfunc(env, lp, contextmask, callback, userhandle);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXcallbacksetfunc");
	}
}

void callbackcandidateispoint(CPXCALLBACKCONTEXTptr context, int* ispoint_p)
{
	int status = CPXcallbackcandidateispoint(context, ispoint_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackcandidateispoint");
	}
}

void callbackgetcandidatepoint(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p)
{
	int status = CPXcallbackgetcandidatepoint(context, x, begin, end, obj_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackgetcandidatepoint");
	}
}

void callbackgetrelaxationpoint(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p)
{
	int status = CPXcallbackgetrelaxationpoint(context, x, begin, end, obj_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackgetrelaxationpoint");
	}
}

void callbackrejectcandidate(CPXCALLBACKCONTEXTptr context, int rcnt, int nzcnt, double const* rhs,
										   char const* sense, int const* rmatbeg, int const* rmatind,
										   double const* rmatval)
{
	int status = CPXcallbackrejectcandidate(context, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackrejectcandidate");
	}
}

void callbackaddusercuts(CPXCALLBACKCONTEXTptr context, int rcnt, int nzcnt, double const* rhs,
									   char const* sense, int const* rmatbeg, int const* rmatind, double const* rmatval,
									   int const* purgeable, int const* local)
{
	int status = CPXcallbackaddusercuts(context, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, purgeable, local);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackaddusercuts");
	}
}

void callbackgetinfoint(CPXCALLBACKCONTEXTptr context, CPXCALLBACKINFO what, CPXINT* data_p)
{
	int status = CPXcallbackgetinfoint(context, what, data_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackgetinfoint");
	}
}

void callbackgetinfodbl(CPXCALLBACKCONTEXTptr context, CPXCALLBACKINFO what, double* data_p)
{
	int status = CPXcallbackgetinfodbl(context, what, data_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackgetinfodbl");
	}
}

void callbackpostheursoln(CPXCALLBACKCONTEXTptr context, int cnt, int const* ind, double const* val,
										double obj, CPXCALLBACKSOLUTIONSTRATEGY strat)
{
	int status = CPXcallbackpostheursoln(context, cnt, ind, val, obj, strat);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackpostheursoln");
	}
}

void callbackgetincumbent(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p)
{
	int status = CPXcallbackgetincumbent(context, x, begin, end, obj_p);
	if (status != 0)
	{
		fail_with_error_message(status, "CPXcallbackgetincumbent");
	}
}

void copybase(CPXCENVptr env, CPXLPptr lp, int const* cstat, int const* rstat)
{
	int status = CPXcopybase(env, lp, cstat, rstat);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXcopybase");
	}
}

void getbase(CPXCENVptr env, CPXCLPptr lp, int* cstat, int* rstat)
{
	int status = CPXgetbase(env, lp, cstat, rstat);
	if (status != 0)
	{
		fail_with_error_message(env, status, "CPXgetbase");
	}
}
} // namespace cplex
} // namespace goc
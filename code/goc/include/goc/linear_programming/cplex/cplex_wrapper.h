//
// Created by Gonzalo Lera Romero.
// Grupo de Optimizacion Combinatoria (GOC).
// Departamento de Computacion - Universidad de Buenos Aires.
//

#ifndef GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_WRAPPER_H
#define GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_WRAPPER_H

#include "ilcplex/cplex.h"

namespace goc
{
namespace cplex
{
CPXENVptr openCPLEX();

CPXLPptr createprob(CPXENVptr env, char const* probname_str);

void closeCPLEX(CPXENVptr* env_p);

void freeprob(CPXENVptr env, CPXLPptr* lp_p);

void newcols(CPXENVptr env, CPXLPptr lp, int ccnt, double const* obj, double const* lb,
					double const* ub, char const* xctype, char** colname);

void delcols(CPXCENVptr env, CPXLPptr lp, int begin, int end);

void chgctype(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, char const* xctype);

void getctype(CPXENVptr env, CPXCLPptr lp, char* xctype, int begin, int end);


void addrows(CPXENVptr env, CPXLPptr lp, int ccnt, int rcnt, int nzcnt, double const* rhs,
					char const* sense, int const* rmatbeg, int const* rmatind, double const* rmatval, char** colname,
					char** rowname);

void getlb(CPXENVptr env, CPXCLPptr lp, double* lb, int begin, int end);

void getub(CPXENVptr env, CPXCLPptr lp, double* ub, int begin, int end);

void getpi(CPXENVptr env, CPXCLPptr lp, double* pi, int begin, int end);

void getdj(CPXCENVptr env, CPXCLPptr lp, double* dj, int begin, int end);

void chgbds(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, char const* lu,
				   double const* bd);

void chgobj(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, double const* values);

int getobjsen(CPXCENVptr env, CPXCLPptr lp);

void getobj(CPXENVptr env, CPXLPptr lp, double* obj, int begin, int end);

void getrhs(CPXCENVptr env, CPXCLPptr lp, double* rhs, int begin, int end);

void chgrhs(CPXCENVptr env, CPXLPptr lp, int cnt, int const* indices, double const* values);

void getcoef(CPXCENVptr env, CPXCLPptr lp, int i, int j, double* coef_p);

void chgcoef(CPXENVptr env, CPXLPptr lp, int i, int j, double newvalue);

void chgobjsen(CPXENVptr env, CPXLPptr lp, int maxormin);

void delrows(CPXENVptr env, CPXLPptr lp, int begin, int end);

CPXLPptr cloneprob(CPXENVptr env, CPXCLPptr lp);

void mipopt(CPXENVptr env, CPXLPptr lp);

void lpopt(CPXENVptr env, CPXLPptr lp);

void solution(CPXENVptr env, CPXCLPptr lp, int* lpstat_p, double* objval_p, double* x, double* pi,
					 double* slack, double* dj);

void getbestobjval(CPXENVptr env, CPXCLPptr lp, double* objval_p);

int getnodecnt(CPXENVptr env, CPXCLPptr lp);

int getnodeleftcnt(CPXENVptr env, CPXCLPptr lp);

bool getobjval(CPXENVptr env, CPXCLPptr lp, double* objval_p);

void getmiprelgap(CPXENVptr env, CPXCLPptr lp, double* gap_p);

int getmipitcnt(CPXENVptr env, CPXCLPptr lp);

int getitcnt(CPXENVptr env, CPXCLPptr lp);

int getstat(CPXENVptr env, CPXCLPptr lp);

void getsense(CPXCENVptr env, CPXCLPptr lp, char* sense, int begin, int end);

void getrows(CPXCENVptr env, CPXCLPptr lp, int* nzcnt_p, int* rmatbeg, int* rmatind, double* rmatval,
			 int rmatspace, int* surplus_p, int begin, int end);

void setintparam(CPXENVptr env, int whichparam, CPXINT newvalue);

void setdblparam(CPXENVptr env, int whichparam, double newvalue);

void setlongparam(CPXENVptr env, int whichparam, long newvalue);

void setstrparam(CPXENVptr env, int whichparam, char const* newvalue);

void getintparam(CPXENVptr env, int whichparam, CPXINT* value_p);

void getlongparam(CPXENVptr env, int whichparam, CPXLONG* value_p);

void getdblparam(CPXENVptr env, int whichparam, double* value_p);

void chgprobtype(CPXENVptr env, CPXLPptr lp, int type);

int getprobtype(CPXENVptr env, CPXCLPptr lp);

void getparamname(CPXENVptr env, int whichparam, char* name_str);

void getparamtype(CPXENVptr env, int whichparam, int* paramtype);

void getparamnum(CPXENVptr env, char const* name_str, int* whichparam_p);

void setlazyconstraintcallbackfunc(CPXENVptr env, int(CPXPUBLIC* lazyconcallback)(CALLBACK_CUT_ARGS),
void* cbhandle);

void setusercutcallbackfunc(CPXENVptr env, int(CPXPUBLIC* cutcallback)(CALLBACK_CUT_ARGS),
void* cbhandle);

void setbranchcallbackfunc(CPXENVptr env, int(CPXPUBLIC* branchcallback)(CALLBACK_BRANCH_ARGS),
void* cbhandle);

void branchcallbackbranchbds(CPXENVptr env, void* cbdata, int wherefrom, int cnt, int const* indices,
									char const* lu, double const* bd, double nodeest, void* userhandle,
									int* seqnum_p);

void branchcallbackbranchconstraints(CPXENVptr env, void* cbdata, int wherefrom, int rcnt, int nzcnt,
											double const* rhs, char const* sense, int const* rmatbeg, int const* rmatind,
											double const* rmatval, double nodeest, void* userhandle, int* seqnum_p);

void cutcallbackadd(CPXENVptr env, void* cbdata, int wherefrom, int nzcnt, double rhs, int sense,
						   int const* cutind, double const* cutval, int purgeable);

void getcallbacknodex(CPXENVptr env, void* cbdata, int wherefrom, double* x, int begin,
							 int end);

void getcallbackinfo(CPXCENVptr env, void* cbdata, int wherefrom, int whichinfo, void* result_p);

void addmipstarts(CPXENVptr env, CPXLPptr lp, int mcnt, int nzcnt, int const* beg,
						 int const* varindices, double const* values, int const* effortlevel, char** mipstartname);

int getnummipstarts(CPXCENVptr env, CPXCLPptr lp);

void delmipstarts(CPXCENVptr env, CPXLPptr lp, int begin, int end);

void getnumcuts(CPXENVptr env, CPXCLPptr lp, int cuttype, int* num_p);

int getnumrows(CPXENVptr env, CPXCLPptr lp);

int getnumcols(CPXENVptr env, CPXCLPptr lp);

void copyorder(CPXENVptr env, CPXLPptr lp, int cnt, int const* indices, int const* priority,
					  int const* direction);

void populate(CPXENVptr env, CPXLPptr lp);

int getsolnpoolnumsolns(CPXENVptr env, CPXCLPptr lp);

void getsolnpoolx(CPXENVptr env, CPXCLPptr lp, int soln, double* x, int begin, int end);

void setdefaults(CPXENVptr env);

void addfuncdest(CPXCENVptr env, CPXCHANNELptr channel, void* handle,
						void(CPXPUBLIC* msgfunction)(void*, const char*));

void delfuncdest(CPXCENVptr env, CPXCHANNELptr channel, void* handle,
						void(CPXPUBLIC* msgfunction)(void*, const char*));

void getchannels(CPXCENVptr env, CPXCHANNELptr* cpxresults_p, CPXCHANNELptr* cpxwarning_p,
						CPXCHANNELptr* cpxerror_p, CPXCHANNELptr* cpxlog_p);

void setinfocallbackfunc(CPXENVptr env,
								int(CPXPUBLIC* callback)(CPXCENVptr, void*, int, void*),
void* cbhandle);

void getcallbacknodeinfo(CPXCENVptr env, void* cbdata, int wherefrom, int nodeindex, int whichinfo,
								void* result_p);

void getcallbacknodeobjval (CPXCENVptr env, void* cbdata, int wherefrom, double* objval_p);

void callbacksetfunc(CPXENVptr env, CPXLPptr lp, CPXLONG contextmask, CPXCALLBACKFUNC callback,
							void* userhandle);

void callbackcandidateispoint(CPXCALLBACKCONTEXTptr context, int* ispoint_p);

void callbackgetcandidatepoint(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p);

void callbackgetrelaxationpoint(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p);

void callbackrejectcandidate(CPXCALLBACKCONTEXTptr context, int rcnt, int nzcnt, double const* rhs,
									char const* sense, int const* rmatbeg, int const* rmatind,
									double const* rmatval);

void callbackaddusercuts(CPXCALLBACKCONTEXTptr context, int rcnt, int nzcnt, double const* rhs, char const* sense,
								int const* rmatbeg, int const* rmatind, double const* rmatval, int const* purgeable,
								int const* local);

void callbackgetinfoint(CPXCALLBACKCONTEXTptr context, CPXCALLBACKINFO what, CPXINT* data_p);

void callbackgetinfodbl(CPXCALLBACKCONTEXTptr context, CPXCALLBACKINFO what, double* data_p);

void callbackpostheursoln(CPXCALLBACKCONTEXTptr context, int cnt, int const* ind, double const* val,
								 double obj, CPXCALLBACKSOLUTIONSTRATEGY strat);

void callbackgetincumbent(CPXCALLBACKCONTEXTptr context, double* x, int begin, int end, double* obj_p);

void callbackgetinfodbl(CPXCALLBACKCONTEXTptr context, CPXCALLBACKINFO what, double* data_p);

void copybase(CPXCENVptr env, CPXLPptr lp, int const* cstat, int const* rstat);

void getbase(CPXCENVptr env, CPXCLPptr lp, int* cstat, int* rstat);
} // namespace cplex
} // namespace goc

#endif //GOC_LINEAR_PROGRAMMING_CPLEX_CPLEX_WRAPPER_H

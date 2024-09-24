/* File: surfcal76module.c
 * This file is auto-generated with f2py (version:1.26.4).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Wed May  8 01:34:44 2024
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifndef PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN
#endif /* PY_SSIZE_T_CLEAN */

/* Unconditionally included */
#include <Python.h>
#include <numpy/npy_os.h>

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "fortranobject.h"
/*need_includes0*/

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *surfcal76_error;
static PyObject *surfcal76_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
typedef char * string;

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif


#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif


#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif


/* See fortranobject.h for definitions. The macros here are provided for BC. */
#define rank f2py_rank
#define shape f2py_shape
#define fshape f2py_shape
#define len f2py_len
#define flen f2py_flen
#define slen f2py_slen
#define size f2py_size


/************************ See f2py2e/cfuncs.py: cfuncs ************************/

static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = surfcal76_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}


static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj)) {
        PyErr_Clear();
        tmp = PyObject_GetAttrString(obj,"real");
    }
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj)) {
        /*pass*/;
    }
    else if (PySequence_Check(obj)) {
        PyErr_Clear();
        tmp = PySequence_GetItem(obj, 0);
    }

    if (tmp) {
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }

    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = surfcal76_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}


static int
float_from_pyobj(float* v, PyObject *obj, const char *errmess)
{
    double d=0.0;
    if (double_from_pyobj(&d,obj,errmess)) {
        *v = (float)d;
        return 1;
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(surfcal,SURFCAL)(void);
extern void F_FUNC(asorder,ASORDER)(int*,int*,float*,float*,int*);
extern void F_FUNC(bsurf,BSURF)(int*,string*,float*,int*,string*,size_t,size_t);
extern void F_FUNC(intsect,INTSECT)(float*,float*,float*,float*,float*,float*,int*,int*,float*,float*);
extern void F_FUNC(assvdw,ASSVDW)(string*,float*,int*,string*,size_t,size_t);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** surfcal **********************************/
static char doc_f2py_rout_surfcal76_surfcal[] = "\
surfcal()\n\nWrapper for ``surfcal``.\
\n";
/* extern void F_FUNC(surfcal,SURFCAL)(void); */
static PyObject *f2py_rout_surfcal76_surfcal(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(void)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    static char *capi_kwlist[] = {NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "|:surfcal76.surfcal",\
        capi_kwlist))
        return NULL;
/*frompyobj*/
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)();
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of surfcal *******************************/

/********************************** asorder **********************************/
static char doc_f2py_rout_surfcal76_asorder[] = "\
asorder(nuniq,ncount,parea,dtmin,iuni)\n\nWrapper for ``asorder``.\
\n\nParameters\n----------\n"
"nuniq : input rank-1 array('i') with bounds (1200)\n"
"ncount : input rank-1 array('i') with bounds (1200)\n"
"parea : input rank-1 array('f') with bounds (1200)\n"
"dtmin : input rank-1 array('f') with bounds (1200)\n"
"iuni : input int";
/* extern void F_FUNC(asorder,ASORDER)(int*,int*,float*,float*,int*); */
static PyObject *f2py_rout_surfcal76_asorder(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,int*,float*,float*,int*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int *nuniq = NULL;
    npy_intp nuniq_Dims[1] = {-1};
    const int nuniq_Rank = 1;
    PyArrayObject *capi_nuniq_as_array = NULL;
    int capi_nuniq_intent = 0;
    PyObject *nuniq_capi = Py_None;
    int *ncount = NULL;
    npy_intp ncount_Dims[1] = {-1};
    const int ncount_Rank = 1;
    PyArrayObject *capi_ncount_as_array = NULL;
    int capi_ncount_intent = 0;
    PyObject *ncount_capi = Py_None;
    float *parea = NULL;
    npy_intp parea_Dims[1] = {-1};
    const int parea_Rank = 1;
    PyArrayObject *capi_parea_as_array = NULL;
    int capi_parea_intent = 0;
    PyObject *parea_capi = Py_None;
    float *dtmin = NULL;
    npy_intp dtmin_Dims[1] = {-1};
    const int dtmin_Rank = 1;
    PyArrayObject *capi_dtmin_as_array = NULL;
    int capi_dtmin_intent = 0;
    PyObject *dtmin_capi = Py_None;
    int iuni = 0;
    PyObject *iuni_capi = Py_None;
    static char *capi_kwlist[] = {"nuniq","ncount","parea","dtmin","iuni",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOO|:surfcal76.asorder",\
        capi_kwlist,&nuniq_capi,&ncount_capi,&parea_capi,&dtmin_capi,&iuni_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable nuniq */
    nuniq_Dims[0]=1200;
    capi_nuniq_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.asorder: failed to create array from the 1st argument `nuniq`";
    capi_nuniq_as_array = ndarray_from_pyobj(  NPY_INT,1,nuniq_Dims,nuniq_Rank,  capi_nuniq_intent,nuniq_capi,capi_errmess);
    if (capi_nuniq_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        nuniq = (int *)(PyArray_DATA(capi_nuniq_as_array));

    /* Processing variable ncount */
    ncount_Dims[0]=1200;
    capi_ncount_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.asorder: failed to create array from the 2nd argument `ncount`";
    capi_ncount_as_array = ndarray_from_pyobj(  NPY_INT,1,ncount_Dims,ncount_Rank,  capi_ncount_intent,ncount_capi,capi_errmess);
    if (capi_ncount_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        ncount = (int *)(PyArray_DATA(capi_ncount_as_array));

    /* Processing variable parea */
    parea_Dims[0]=1200;
    capi_parea_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.asorder: failed to create array from the 3rd argument `parea`";
    capi_parea_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,parea_Dims,parea_Rank,  capi_parea_intent,parea_capi,capi_errmess);
    if (capi_parea_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        parea = (float *)(PyArray_DATA(capi_parea_as_array));

    /* Processing variable dtmin */
    dtmin_Dims[0]=1200;
    capi_dtmin_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.asorder: failed to create array from the 4th argument `dtmin`";
    capi_dtmin_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,dtmin_Dims,dtmin_Rank,  capi_dtmin_intent,dtmin_capi,capi_errmess);
    if (capi_dtmin_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        dtmin = (float *)(PyArray_DATA(capi_dtmin_as_array));

    /* Processing variable iuni */
        f2py_success = int_from_pyobj(&iuni,iuni_capi,"surfcal76.asorder() 5th argument (iuni) can't be converted to int");
    if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(nuniq,ncount,parea,dtmin,&iuni);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    } /*if (f2py_success) of iuni*/
    /* End of cleaning variable iuni */
    if((PyObject *)capi_dtmin_as_array!=dtmin_capi) {
        Py_XDECREF(capi_dtmin_as_array); }
    }  /* if (capi_dtmin_as_array == NULL) ... else of dtmin */
    /* End of cleaning variable dtmin */
    if((PyObject *)capi_parea_as_array!=parea_capi) {
        Py_XDECREF(capi_parea_as_array); }
    }  /* if (capi_parea_as_array == NULL) ... else of parea */
    /* End of cleaning variable parea */
    if((PyObject *)capi_ncount_as_array!=ncount_capi) {
        Py_XDECREF(capi_ncount_as_array); }
    }  /* if (capi_ncount_as_array == NULL) ... else of ncount */
    /* End of cleaning variable ncount */
    if((PyObject *)capi_nuniq_as_array!=nuniq_capi) {
        Py_XDECREF(capi_nuniq_as_array); }
    }  /* if (capi_nuniq_as_array == NULL) ... else of nuniq */
    /* End of cleaning variable nuniq */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of asorder *******************************/

/*********************************** bsurf ***********************************/
static char doc_f2py_rout_surfcal76_bsurf[] = "\
bsurf(rayflag,attype,vdwr,nvdwt,restyp)\n\nWrapper for ``bsurf``.\
\n\nParameters\n----------\n"
"rayflag : input int\n"
"attype : input rank-1 array('S') with bounds (50)\n"
"vdwr : input rank-1 array('f') with bounds (50)\n"
"nvdwt : input int\n"
"restyp : input rank-1 array('S') with bounds (50)";
/* extern void F_FUNC(bsurf,BSURF)(int*,string*,float*,int*,string*,size_t,size_t); */
static PyObject *f2py_rout_surfcal76_bsurf(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(int*,string*,float*,int*,string*,size_t,size_t)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    int rayflag = 0;
    PyObject *rayflag_capi = Py_None;
    string *attype = NULL;
    npy_intp attype_Dims[1] = {-1};
    const int attype_Rank = 1;
    PyArrayObject *capi_attype_as_array = NULL;
    int capi_attype_intent = 0;
    int slen(attype) = 0;
    PyObject *attype_capi = Py_None;
    float *vdwr = NULL;
    npy_intp vdwr_Dims[1] = {-1};
    const int vdwr_Rank = 1;
    PyArrayObject *capi_vdwr_as_array = NULL;
    int capi_vdwr_intent = 0;
    PyObject *vdwr_capi = Py_None;
    int nvdwt = 0;
    PyObject *nvdwt_capi = Py_None;
    string *restyp = NULL;
    npy_intp restyp_Dims[1] = {-1};
    const int restyp_Rank = 1;
    PyArrayObject *capi_restyp_as_array = NULL;
    int capi_restyp_intent = 0;
    int slen(restyp) = 0;
    PyObject *restyp_capi = Py_None;
    static char *capi_kwlist[] = {"rayflag","attype","vdwr","nvdwt","restyp",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOO|:surfcal76.bsurf",\
        capi_kwlist,&rayflag_capi,&attype_capi,&vdwr_capi,&nvdwt_capi,&restyp_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable rayflag */
        f2py_success = int_from_pyobj(&rayflag,rayflag_capi,"surfcal76.bsurf() 1st argument (rayflag) can't be converted to int");
    if (f2py_success) {
    /* Processing variable attype */
    attype_Dims[0]=50;
    capi_attype_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.bsurf: failed to create array from the 2nd argument `attype`";
    capi_attype_as_array = ndarray_from_pyobj(  NPY_STRING,3,attype_Dims,attype_Rank,  capi_attype_intent,attype_capi,capi_errmess);
    if (capi_attype_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        attype = (string *)(PyArray_DATA(capi_attype_as_array));

    slen(attype) = f2py_itemsize(attype);
    /* Processing variable vdwr */
    vdwr_Dims[0]=50;
    capi_vdwr_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.bsurf: failed to create array from the 3rd argument `vdwr`";
    capi_vdwr_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,vdwr_Dims,vdwr_Rank,  capi_vdwr_intent,vdwr_capi,capi_errmess);
    if (capi_vdwr_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        vdwr = (float *)(PyArray_DATA(capi_vdwr_as_array));

    /* Processing variable nvdwt */
        f2py_success = int_from_pyobj(&nvdwt,nvdwt_capi,"surfcal76.bsurf() 4th argument (nvdwt) can't be converted to int");
    if (f2py_success) {
    /* Processing variable restyp */
    restyp_Dims[0]=50;
    capi_restyp_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.bsurf: failed to create array from the 5th argument `restyp`";
    capi_restyp_as_array = ndarray_from_pyobj(  NPY_STRING,3,restyp_Dims,restyp_Rank,  capi_restyp_intent,restyp_capi,capi_errmess);
    if (capi_restyp_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        restyp = (string *)(PyArray_DATA(capi_restyp_as_array));

    slen(restyp) = f2py_itemsize(restyp);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(&rayflag,attype,vdwr,&nvdwt,restyp,slen(attype),slen(restyp));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    if((PyObject *)capi_restyp_as_array!=restyp_capi) {
        Py_XDECREF(capi_restyp_as_array); }
    }  /* if (capi_restyp_as_array == NULL) ... else of restyp */
    /* End of cleaning variable restyp */
    } /*if (f2py_success) of nvdwt*/
    /* End of cleaning variable nvdwt */
    if((PyObject *)capi_vdwr_as_array!=vdwr_capi) {
        Py_XDECREF(capi_vdwr_as_array); }
    }  /* if (capi_vdwr_as_array == NULL) ... else of vdwr */
    /* End of cleaning variable vdwr */
    if((PyObject *)capi_attype_as_array!=attype_capi) {
        Py_XDECREF(capi_attype_as_array); }
    }  /* if (capi_attype_as_array == NULL) ... else of attype */
    /* End of cleaning variable attype */
    } /*if (f2py_success) of rayflag*/
    /* End of cleaning variable rayflag */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************** end of bsurf ********************************/

/********************************** intsect **********************************/
static char doc_f2py_rout_surfcal76_intsect[] = "\
intsect(spc,dcl,dcm,dcn,ptc,sphrad,iretro,iflag,dmin,dist6)\n\nWrapper for ``intsect``.\
\n\nParameters\n----------\n"
"spc : input rank-1 array('f') with bounds (3)\n"
"dcl : input float\n"
"dcm : input float\n"
"dcn : input float\n"
"ptc : input rank-1 array('f') with bounds (3)\n"
"sphrad : input float\n"
"iretro : input int\n"
"iflag : input int\n"
"dmin : input float\n"
"dist6 : input float";
/* extern void F_FUNC(intsect,INTSECT)(float*,float*,float*,float*,float*,float*,int*,int*,float*,float*); */
static PyObject *f2py_rout_surfcal76_intsect(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(float*,float*,float*,float*,float*,float*,int*,int*,float*,float*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    float *spc = NULL;
    npy_intp spc_Dims[1] = {-1};
    const int spc_Rank = 1;
    PyArrayObject *capi_spc_as_array = NULL;
    int capi_spc_intent = 0;
    PyObject *spc_capi = Py_None;
    float dcl = 0;
    PyObject *dcl_capi = Py_None;
    float dcm = 0;
    PyObject *dcm_capi = Py_None;
    float dcn = 0;
    PyObject *dcn_capi = Py_None;
    float *ptc = NULL;
    npy_intp ptc_Dims[1] = {-1};
    const int ptc_Rank = 1;
    PyArrayObject *capi_ptc_as_array = NULL;
    int capi_ptc_intent = 0;
    PyObject *ptc_capi = Py_None;
    float sphrad = 0;
    PyObject *sphrad_capi = Py_None;
    int iretro = 0;
    PyObject *iretro_capi = Py_None;
    int iflag = 0;
    PyObject *iflag_capi = Py_None;
    float dmin = 0;
    PyObject *dmin_capi = Py_None;
    float dist6 = 0;
    PyObject *dist6_capi = Py_None;
    static char *capi_kwlist[] = {"spc","dcl","dcm","dcn","ptc","sphrad","iretro","iflag","dmin","dist6",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOOOOOOOO|:surfcal76.intsect",\
        capi_kwlist,&spc_capi,&dcl_capi,&dcm_capi,&dcn_capi,&ptc_capi,&sphrad_capi,&iretro_capi,&iflag_capi,&dmin_capi,&dist6_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable spc */
    spc_Dims[0]=3;
    capi_spc_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.intsect: failed to create array from the 1st argument `spc`";
    capi_spc_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,spc_Dims,spc_Rank,  capi_spc_intent,spc_capi,capi_errmess);
    if (capi_spc_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        spc = (float *)(PyArray_DATA(capi_spc_as_array));

    /* Processing variable dcl */
        f2py_success = float_from_pyobj(&dcl,dcl_capi,"surfcal76.intsect() 2nd argument (dcl) can't be converted to float");
    if (f2py_success) {
    /* Processing variable dcm */
        f2py_success = float_from_pyobj(&dcm,dcm_capi,"surfcal76.intsect() 3rd argument (dcm) can't be converted to float");
    if (f2py_success) {
    /* Processing variable dcn */
        f2py_success = float_from_pyobj(&dcn,dcn_capi,"surfcal76.intsect() 4th argument (dcn) can't be converted to float");
    if (f2py_success) {
    /* Processing variable ptc */
    ptc_Dims[0]=3;
    capi_ptc_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.intsect: failed to create array from the 5th argument `ptc`";
    capi_ptc_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,ptc_Dims,ptc_Rank,  capi_ptc_intent,ptc_capi,capi_errmess);
    if (capi_ptc_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        ptc = (float *)(PyArray_DATA(capi_ptc_as_array));

    /* Processing variable sphrad */
        f2py_success = float_from_pyobj(&sphrad,sphrad_capi,"surfcal76.intsect() 6th argument (sphrad) can't be converted to float");
    if (f2py_success) {
    /* Processing variable iretro */
        f2py_success = int_from_pyobj(&iretro,iretro_capi,"surfcal76.intsect() 7th argument (iretro) can't be converted to int");
    if (f2py_success) {
    /* Processing variable iflag */
        f2py_success = int_from_pyobj(&iflag,iflag_capi,"surfcal76.intsect() 8th argument (iflag) can't be converted to int");
    if (f2py_success) {
    /* Processing variable dmin */
        f2py_success = float_from_pyobj(&dmin,dmin_capi,"surfcal76.intsect() 9th argument (dmin) can't be converted to float");
    if (f2py_success) {
    /* Processing variable dist6 */
        f2py_success = float_from_pyobj(&dist6,dist6_capi,"surfcal76.intsect() 10th argument (dist6) can't be converted to float");
    if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(spc,&dcl,&dcm,&dcn,ptc,&sphrad,&iretro,&iflag,&dmin,&dist6);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    } /*if (f2py_success) of dist6*/
    /* End of cleaning variable dist6 */
    } /*if (f2py_success) of dmin*/
    /* End of cleaning variable dmin */
    } /*if (f2py_success) of iflag*/
    /* End of cleaning variable iflag */
    } /*if (f2py_success) of iretro*/
    /* End of cleaning variable iretro */
    } /*if (f2py_success) of sphrad*/
    /* End of cleaning variable sphrad */
    if((PyObject *)capi_ptc_as_array!=ptc_capi) {
        Py_XDECREF(capi_ptc_as_array); }
    }  /* if (capi_ptc_as_array == NULL) ... else of ptc */
    /* End of cleaning variable ptc */
    } /*if (f2py_success) of dcn*/
    /* End of cleaning variable dcn */
    } /*if (f2py_success) of dcm*/
    /* End of cleaning variable dcm */
    } /*if (f2py_success) of dcl*/
    /* End of cleaning variable dcl */
    if((PyObject *)capi_spc_as_array!=spc_capi) {
        Py_XDECREF(capi_spc_as_array); }
    }  /* if (capi_spc_as_array == NULL) ... else of spc */
    /* End of cleaning variable spc */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of intsect *******************************/

/*********************************** assvdw ***********************************/
static char doc_f2py_rout_surfcal76_assvdw[] = "\
assvdw(attype,vdwr,nvdwt,restyp)\n\nWrapper for ``assvdw``.\
\n\nParameters\n----------\n"
"attype : input rank-1 array('S') with bounds (50)\n"
"vdwr : input rank-1 array('f') with bounds (50)\n"
"nvdwt : input int\n"
"restyp : input rank-1 array('S') with bounds (50)";
/* extern void F_FUNC(assvdw,ASSVDW)(string*,float*,int*,string*,size_t,size_t); */
static PyObject *f2py_rout_surfcal76_assvdw(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string*,float*,int*,string*,size_t,size_t)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

    string *attype = NULL;
    npy_intp attype_Dims[1] = {-1};
    const int attype_Rank = 1;
    PyArrayObject *capi_attype_as_array = NULL;
    int capi_attype_intent = 0;
    int slen(attype) = 0;
    PyObject *attype_capi = Py_None;
    float *vdwr = NULL;
    npy_intp vdwr_Dims[1] = {-1};
    const int vdwr_Rank = 1;
    PyArrayObject *capi_vdwr_as_array = NULL;
    int capi_vdwr_intent = 0;
    PyObject *vdwr_capi = Py_None;
    int nvdwt = 0;
    PyObject *nvdwt_capi = Py_None;
    string *restyp = NULL;
    npy_intp restyp_Dims[1] = {-1};
    const int restyp_Rank = 1;
    PyArrayObject *capi_restyp_as_array = NULL;
    int capi_restyp_intent = 0;
    int slen(restyp) = 0;
    PyObject *restyp_capi = Py_None;
    static char *capi_kwlist[] = {"attype","vdwr","nvdwt","restyp",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOO|:surfcal76.assvdw",\
        capi_kwlist,&attype_capi,&vdwr_capi,&nvdwt_capi,&restyp_capi))
        return NULL;
/*frompyobj*/
    /* Processing variable attype */
    attype_Dims[0]=50;
    capi_attype_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.assvdw: failed to create array from the 1st argument `attype`";
    capi_attype_as_array = ndarray_from_pyobj(  NPY_STRING,3,attype_Dims,attype_Rank,  capi_attype_intent,attype_capi,capi_errmess);
    if (capi_attype_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        attype = (string *)(PyArray_DATA(capi_attype_as_array));

    slen(attype) = f2py_itemsize(attype);
    /* Processing variable vdwr */
    vdwr_Dims[0]=50;
    capi_vdwr_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.assvdw: failed to create array from the 2nd argument `vdwr`";
    capi_vdwr_as_array = ndarray_from_pyobj(  NPY_FLOAT,1,vdwr_Dims,vdwr_Rank,  capi_vdwr_intent,vdwr_capi,capi_errmess);
    if (capi_vdwr_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        vdwr = (float *)(PyArray_DATA(capi_vdwr_as_array));

    /* Processing variable nvdwt */
        f2py_success = int_from_pyobj(&nvdwt,nvdwt_capi,"surfcal76.assvdw() 3rd argument (nvdwt) can't be converted to int");
    if (f2py_success) {
    /* Processing variable restyp */
    restyp_Dims[0]=50;
    capi_restyp_intent |= F2PY_INTENT_IN;
    const char * capi_errmess = "surfcal76.surfcal76.assvdw: failed to create array from the 4th argument `restyp`";
    capi_restyp_as_array = ndarray_from_pyobj(  NPY_STRING,3,restyp_Dims,restyp_Rank,  capi_restyp_intent,restyp_capi,capi_errmess);
    if (capi_restyp_as_array == NULL) {
        PyObject* capi_err = PyErr_Occurred();
        if (capi_err == NULL) {
            capi_err = surfcal76_error;
            PyErr_SetString(capi_err, capi_errmess);
        }
    } else {
        restyp = (string *)(PyArray_DATA(capi_restyp_as_array));

    slen(restyp) = f2py_itemsize(restyp);
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
                (*f2py_func)(attype,vdwr,&nvdwt,restyp,slen(attype),slen(restyp));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
    if((PyObject *)capi_restyp_as_array!=restyp_capi) {
        Py_XDECREF(capi_restyp_as_array); }
    }  /* if (capi_restyp_as_array == NULL) ... else of restyp */
    /* End of cleaning variable restyp */
    } /*if (f2py_success) of nvdwt*/
    /* End of cleaning variable nvdwt */
    if((PyObject *)capi_vdwr_as_array!=vdwr_capi) {
        Py_XDECREF(capi_vdwr_as_array); }
    }  /* if (capi_vdwr_as_array == NULL) ... else of vdwr */
    /* End of cleaning variable vdwr */
    if((PyObject *)capi_attype_as_array!=attype_capi) {
        Py_XDECREF(capi_attype_as_array); }
    }  /* if (capi_attype_as_array == NULL) ... else of attype */
    /* End of cleaning variable attype */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/******************************* end of assvdw *******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

static FortranDataDef f2py_resids_def[] = {
  {"nseg",0,{{-1}},NPY_INT, 1},
  {"secnm",1,{{2000}},NPY_STRING, 3},
  {"secsq",1,{{2000}},NPY_INT, 1},
  {"atnm",1,{{20000}},NPY_STRING, 4},
  {"iats",1,{{2000}},NPY_INT, 1},
  {"iate",1,{{2000}},NPY_INT, 1},
  {"seqcd",1,{{2000}},NPY_INT, 1},
  {"natom",0,{{-1}},NPY_INT, 1},
  {NULL}
};
static void f2py_setup_resids(char *nseg,char *secnm,char *secsq,char *atnm,char *iats,char *iate,char *seqcd,char *natom) {
  int i_f2py=0;
  f2py_resids_def[i_f2py++].data = nseg;
  f2py_resids_def[i_f2py++].data = secnm;
  f2py_resids_def[i_f2py++].data = secsq;
  f2py_resids_def[i_f2py++].data = atnm;
  f2py_resids_def[i_f2py++].data = iats;
  f2py_resids_def[i_f2py++].data = iate;
  f2py_resids_def[i_f2py++].data = seqcd;
  f2py_resids_def[i_f2py++].data = natom;
}
extern void F_FUNC(f2pyinitresids,F2PYINITRESIDS)(void(*)(char*,char*,char*,char*,char*,char*,char*,char*));
static void f2py_init_resids(void) {
  F_FUNC(f2pyinitresids,F2PYINITRESIDS)(f2py_setup_resids);
}

static FortranDataDef f2py_coord_def[] = {
  {"cr",2,{{20000,3}},NPY_FLOAT, 1},
  {"secsqe",1,{{2000}},NPY_STRING, 1},
  {"rescid",1,{{2000}},NPY_STRING, 1},
  {NULL}
};
static void f2py_setup_coord(char *cr,char *secsqe,char *rescid) {
  int i_f2py=0;
  f2py_coord_def[i_f2py++].data = cr;
  f2py_coord_def[i_f2py++].data = secsqe;
  f2py_coord_def[i_f2py++].data = rescid;
}
extern void F_FUNC(f2pyinitcoord,F2PYINITCOORD)(void(*)(char*,char*,char*));
static void f2py_init_coord(void) {
  F_FUNC(f2pyinitcoord,F2PYINITCOORD)(f2py_setup_coord);
}

static FortranDataDef f2py_anput_def[] = {
  {"resinf",1,{{5}},NPY_STRING, 13},
  {"residen",0,{{-1}},NPY_STRING, 13},
  {NULL}
};
static void f2py_setup_anput(char *resinf,char *residen) {
  int i_f2py=0;
  f2py_anput_def[i_f2py++].data = resinf;
  f2py_anput_def[i_f2py++].data = residen;
}
extern void F_FUNC(f2pyinitanput,F2PYINITANPUT)(void(*)(char*,char*));
static void f2py_init_anput(void) {
  F_FUNC(f2pyinitanput,F2PYINITANPUT)(f2py_setup_anput);
}

static FortranDataDef f2py_vdw_def[] = {
  {"vdwrad",1,{{20000}},NPY_FLOAT, 1},
  {"inter",0,{{-1}},NPY_INT, 1},
  {"icter",0,{{-1}},NPY_INT, 1},
  {NULL}
};
static void f2py_setup_vdw(char *vdwrad,char *inter,char *icter) {
  int i_f2py=0;
  f2py_vdw_def[i_f2py++].data = vdwrad;
  f2py_vdw_def[i_f2py++].data = inter;
  f2py_vdw_def[i_f2py++].data = icter;
}
extern void F_FUNC(f2pyinitvdw,F2PYINITVDW)(void(*)(char*,char*,char*));
static void f2py_init_vdw(void) {
  F_FUNC(f2pyinitvdw,F2PYINITVDW)(f2py_setup_vdw);
}

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
    {"surfcal",-1,{{-1}},0,0,(char *)  F_FUNC(surfcal,SURFCAL),  (f2py_init_func)f2py_rout_surfcal76_surfcal,doc_f2py_rout_surfcal76_surfcal},
    {"asorder",-1,{{-1}},0,0,(char *)  F_FUNC(asorder,ASORDER),  (f2py_init_func)f2py_rout_surfcal76_asorder,doc_f2py_rout_surfcal76_asorder},
    {"bsurf",-1,{{-1}},0,0,(char *)  F_FUNC(bsurf,BSURF),  (f2py_init_func)f2py_rout_surfcal76_bsurf,doc_f2py_rout_surfcal76_bsurf},
    {"intsect",-1,{{-1}},0,0,(char *)  F_FUNC(intsect,INTSECT),  (f2py_init_func)f2py_rout_surfcal76_intsect,doc_f2py_rout_surfcal76_intsect},
    {"assvdw",-1,{{-1}},0,0,(char *)  F_FUNC(assvdw,ASSVDW),  (f2py_init_func)f2py_rout_surfcal76_assvdw,doc_f2py_rout_surfcal76_assvdw},

/*eof routine_defs*/
    {NULL}
};

static PyMethodDef f2py_module_methods[] = {

    {NULL,NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "surfcal76",
    NULL,
    -1,
    f2py_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit_surfcal76(void) {
    int i;
    PyObject *m,*d, *s, *tmp;
    m = surfcal76_module = PyModule_Create(&moduledef);
    Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
    import_array();
    if (PyErr_Occurred())
        {PyErr_SetString(PyExc_ImportError, "can't initialize module surfcal76 (failed to import numpy)"); return m;}
    d = PyModule_GetDict(m);
    s = PyUnicode_FromString("1.26.4");
    PyDict_SetItemString(d, "__version__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString(
        "This module 'surfcal76' is auto-generated with f2py (version:1.26.4).\nFunctions:\n"
"    surfcal()\n"
"    asorder(nuniq,ncount,parea,dtmin,iuni)\n"
"    bsurf(rayflag,attype,vdwr,nvdwt,restyp)\n"
"    intsect(spc,dcl,dcm,dcn,ptc,sphrad,iretro,iflag,dmin,dist6)\n"
"    assvdw(attype,vdwr,nvdwt,restyp)\n"
"COMMON blocks:\n""  /resids/ nseg,secnm(2000),secsq(2000),atnm(20000),iats(2000),iate(2000),seqcd(2000),natom\n""  /coord/ cr(20000,3),secsqe(2000),rescid(2000)\n""  /anput/ resinf(5),residen\n""  /vdw/ vdwrad(20000),inter,icter\n"".");
    PyDict_SetItemString(d, "__doc__", s);
    Py_DECREF(s);
    s = PyUnicode_FromString("1.26.4");
    PyDict_SetItemString(d, "__f2py_numpy_version__", s);
    Py_DECREF(s);
    surfcal76_error = PyErr_NewException ("surfcal76.error", NULL, NULL);
    /*
     * Store the error object inside the dict, so that it could get deallocated.
     * (in practice, this is a module, so it likely will not and cannot.)
     */
    PyDict_SetItemString(d, "_surfcal76_error", surfcal76_error);
    Py_DECREF(surfcal76_error);
    for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
        tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
        PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
        Py_DECREF(tmp);
    }





/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

  tmp = PyFortranObject_New(f2py_resids_def,f2py_init_resids);
  if (tmp == NULL) return NULL;
  if (F2PyDict_SetItemString(d, "resids", tmp) == -1) return NULL;
  Py_DECREF(tmp);
  tmp = PyFortranObject_New(f2py_coord_def,f2py_init_coord);
  if (tmp == NULL) return NULL;
  if (F2PyDict_SetItemString(d, "coord", tmp) == -1) return NULL;
  Py_DECREF(tmp);
  tmp = PyFortranObject_New(f2py_anput_def,f2py_init_anput);
  if (tmp == NULL) return NULL;
  if (F2PyDict_SetItemString(d, "anput", tmp) == -1) return NULL;
  Py_DECREF(tmp);
  tmp = PyFortranObject_New(f2py_vdw_def,f2py_init_vdw);
  if (tmp == NULL) return NULL;
  if (F2PyDict_SetItemString(d, "vdw", tmp) == -1) return NULL;
  Py_DECREF(tmp);
/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
    if (! PyErr_Occurred())
        on_exit(f2py_report_on_exit,(void*)"surfcal76");
#endif
    return m;
}
#ifdef __cplusplus
}
#endif

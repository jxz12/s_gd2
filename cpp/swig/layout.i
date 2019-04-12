%module layout
%{
    #define SWIG_FILE_WITH_INIT
    extern void sgd_unweighted(int n, double* X, int m, double* I, double* J, int t_max, double eps, bool horizontal);
    extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta, bool horizontal);
%}

%include "numpy.i"
%init %{
    import_array();
%}

// vertex positions
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}
// edge indices
%apply (double* IN_ARRAY1, int DIM1){(double* E, int m)}
// direct MDS
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}


extern void sgd_unweighted(int n, double* X, int m, double* I, double* J, int t_max, double eps, bool horizontal);
extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta, bool horizontal);

%rename (sgd_unweighted) np_sgd_unweighted;
%exception np_sgd_unweighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (sgd_direct) np_sgd_direct;
%exception np_sgd_direct {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
    void np_sgd_direct(double* X, int n, int kd,
                       double* d, int len_d,
                       double* w, int len_w,
                       double* eta, int len_eta) {

        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "only 2D positions are currently supported");
            return;
        }
        int nC2 = (n*(n-1))/2;
        if (len_d != nC2 || len_w != nC2) {
            PyErr_Format(PyExc_ValueError, "d or w not right length for condensed distance matrix");
            return;
        }
        sgd_direct(n, X, d, w, len_eta, eta, false);
    }
    void np_sgd_unweighted(double* X, int n, int kd,
                           double* I, int m,
                           double* J, int m2,
                           int t_max, double eps) {

        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "only 2D positions are currently supported");
            return;
        }
        if (m != m2) {
            PyErr_Format(PyExc_ValueError, "arrays of indices do not have same length");
            return;
        }
        sgd_unweighted(n, X, m, I, J, t_max, eps, false);
    }
%}


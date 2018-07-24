%module s_gd2
%{
    #define SWIG_FILE_WITH_INIT
    extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
    extern void sgd_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double mu_min);
%}

%include "numpy.i"
%init %{
    import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}
%apply (int* IN_ARRAY1, int DIM1){(int* I, int len_I),
                                  (int* J, int len_J)}

extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
extern void sgd_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double mu_min);

%rename (sgd_direct) my_sgd_direct;
%exception my_sgd_direct {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (sgd_unweighted) my_sgd_unweighted;
%exception my_sgd_unweighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
    void my_sgd_direct(double* X, int n, int kd,
                       double* d, int len_d,
                       double* w, int len_w,
                       double* eta, int len_eta) {

        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "only 2D positions are currently supported");
            return;
        }
        int nC2 = (n*(n-1))/2;
        if (len_d != nC2 || len_w != nC2) {
            PyErr_Format(PyExc_ValueError, "d or w not right length for condensed distnce matrix");
            return;
        }
        sgd_direct(n, X, d, w, len_eta, eta);
    }
    void my_sgd_unweighted(double* X, int n, int kd,
                           int* I, int len_I,
                           int* J, int len_J,
                           int t_max, double mu_min) {

        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "only 2D positions are currently supported");
            return;
        }
        if (len_I != len_J) {
            PyErr_Format(PyExc_ValueError, "number of I indices does not equal number of J");
            return;
        }
        int m = len_I;
        sgd_unweighted(n, X, m, I, J, t_max, mu_min);
    }
%}


%module s_gd2
%{
    #define SWIG_FILE_WITH_INIT
    extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
    extern void sgd_direct_horizontal(int n, double* X, double* d, double* w, int t_max, double* eta);
%}

%include "numpy.i"
%init %{
    import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}

extern void sgd_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
extern void sgd_direct_horizontal(int n, double* X, double* d, double* w, int t_max, double* eta);

%rename (sgd_direct) my_sgd_direct;
%exception my_sgd_direct {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (sgd_direct_horizontal) my_sgd_direct_horizontal;
%exception my_sgd_direct_horizontal {
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
            PyErr_Format(PyExc_ValueError, "d or w not right length for condensed distance matrix");
            return;
        }
        sgd_direct(n, X, d, w, len_eta, eta);
    }
    void my_sgd_direct_horizontal(double* X, int n, int kd,
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
        sgd_direct_horizontal(n, X, d, w, len_eta, eta);
    }
%}


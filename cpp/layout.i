%module s_gd2
%{
    #define SWIG_FILE_WITH_INIT
    #include "layout.hpp"
%}

%include "numpy.i"
%init %{
    import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}

%include "layout.hpp"

%rename (sgd) my_layout;
%exception my_layout {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
    void my_layout(double* X, int n, int kd,
                   double* d, int len_d,
                   double* w, int len_w,
                   double* eta, int len_eta) {

        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "positions not 2D");
            return;
        }
        int nC2 = (n*(n-1))/2;
        if (len_d != nC2 || len_w != nC2) {
            PyErr_Format(PyExc_ValueError, "Arrays not nC2");
            return;
        }
        sgd(n, X, d, w, len_eta, eta);
    }
%}


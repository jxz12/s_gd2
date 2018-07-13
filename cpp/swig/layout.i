%module s_gd2
%{
    #define SWIG_FILE_WITH_INIT
    extern void sgd(int n, double* X, double* d, double* w, int t_max, double* eta);
%}

%include "numpy.i"
%init %{
    import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}
// %apply (double* IN_ARRAY2, int DIM1, int DIM2){(double* d, int n_d, int m_d),
//                                                (double* w, int n_w, int m_w)}
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}

extern void sgd(int n, double* X, double* d, double* w, int t_max, double* eta);

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
            PyErr_Format(PyExc_ValueError, "d or w not right length for condensed distnce matrix");
            return;
        }
        sgd(n, X, d, w, len_eta, eta);
    }
%}


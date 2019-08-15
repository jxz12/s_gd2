%module layout
%{
    #define SWIG_FILE_WITH_INIT
    #include "layout.hpp"
%}

%include "numpy.i"
%init %{
    import_array();
%}



// vertex positions
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}

// edge indices
%apply (int* IN_ARRAY1, int DIM1){(int* I, int len_I),
                                  (int* J, int len_J)}
%apply (double* IN_ARRAY1, int DIM1){(double* V, int len_V)}

// for direct MDS with weights given
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}
#include "layout.hpp"

%rename (layout_unweighted) np_layout_unweighted;
%exception np_layout_unweighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_weighted) np_layout_weighted;
%exception np_layout_weighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_unweighted_convergent) np_layout_unweighted_convergent;
%exception np_layout_unweighted_convergent {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_weighted_convergent) np_layout_weighted_convergent;
%exception np_layout_weighted_convergent {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%rename (layout_sparse_unweighted) np_layout_sparse_unweighted;
%exception np_layout_sparse_unweighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_sparse_weighted) np_layout_sparse_weighted;
%exception np_layout_sparse_weighted {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%rename (mds_direct) np_mds_direct;
%exception np_mds_direct {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
    void dimension_check(int kd) {
        if (kd != 2) {
            PyErr_Format(PyExc_ValueError, "only 2D positions are currently supported");
            return;
        }
    }
    void unweighted_edge_check(int len_I, int len_J) {
        if (len_I != len_J) {
            PyErr_Format(PyExc_ValueError, "arrays of indices do not have same length");
            return;
        }
    }
    void weighted_edge_check(int len_I, int len_J, int len_V) {
        if (len_I != len_J || len_J != len_V) {
            PyErr_Format(PyExc_ValueError, "arrays of indices do not have same length");
            return;
        }
    }
    void np_layout_unweighted(double* X, int n, int kd,
                              int* I, int len_I,
                              int* J, int len_J,
                              int t_max, double eps) {

        dimension_check(kd);
        unweighted_edge_check(len_I, len_J);
        layout_unweighted(n, X, len_I, I, J, t_max, eps);
    }
    void np_layout_weighted(double* X, int n, int kd,
                            int* I, int len_I,
                            int* J, int len_J,
                            double* V, int len_V,
                            int t_max, double eps) {

        dimension_check(kd);
        weighted_edge_check(len_I, len_J, len_V);
        layout_weighted(n, X, len_I, I, J, V, t_max, eps);
    }
    void np_layout_unweighted_convergent(double* X, int n, int kd,
                                         int* I, int len_I,
                                         int* J, int len_J,
                                         int t_max, double eps, double delta, int t_maxmax) {

        dimension_check(kd);
        unweighted_edge_check(len_I, len_J);
        layout_unweighted_convergent(n, X, len_I, I, J, t_max, eps, delta, t_maxmax);
    }
    void np_layout_weighted_convergent(double* X, int n, int kd,
                                       int* I, int len_I,
                                       int* J, int len_J,
                                       double* V, int len_V,
                                       int t_max, double eps, double delta, int t_maxmax) {

        dimension_check(kd);
        weighted_edge_check(len_I, len_J, len_V);
        layout_weighted_convergent(n, X, len_I, I, J, V, t_max, eps, delta, t_maxmax);
    }
    void np_layout_sparse_unweighted(double* X, int n, int kd,
                                     int* I, int len_I,
                                     int* J, int len_J,
                                     int p, int t_max, double eps) {

        dimension_check(kd);
        unweighted_edge_check(len_I, len_J);
        layout_sparse_unweighted(n, X, len_I, I, J, p, t_max, eps);
    }
    void np_layout_sparse_weighted(double* X, int n, int kd,
                                     int* I, int len_I,
                                     int* J, int len_J,
                                     double* V, int len_V,
                                     int p, int t_max, double eps) {

        dimension_check(kd);
        weighted_edge_check(len_I, len_J, len_V);
        layout_sparse_weighted(n, X, len_I, I, J, V, p, t_max, eps);
    }

    void np_mds_direct(double* X, int n, int kd,
                       double* d, int len_d,
                       double* w, int len_w,
                       double* eta, int len_eta) {

        dimension_check(kd);
        int nC2 = (n*(n-1))/2;
        if (len_d != nC2 || len_w != nC2) {
            PyErr_Format(PyExc_ValueError, "d or w not right length for condensed distance matrix");
            return;
        }
        mds_direct(n, X, d, w, len_eta, eta);
    }
%}


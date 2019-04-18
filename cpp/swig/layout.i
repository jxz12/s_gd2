%module layout
%{
    #define SWIG_FILE_WITH_INIT
    extern void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps);
    extern void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps);
    extern void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
    extern void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);
    extern void layout_unweighted_focus(int n, double* X, int m, int* I, int* J, int f, int t_max, double eps, int t_maxmax);
    extern void layout_weighted_focus(int n, double* X, int m, int* I, int* J, int f, double* V, int t_max, double eps, int t_maxmax);
    extern void layout_unweighted_horizontal(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
    extern void layout_weighted_horizontal(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);
    extern void mds_direct(int n, double* X, double* d, double* w, int t_max, double* eta);
%}

%include "numpy.i"
%init %{
    import_array();
%}





// vertex positions
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2){(double* X, int n, int kd)}

// edge indices
%apply (int* IN_ARRAY1, int DIM1){(int* I, int len_I),
                                  (int* J, int len_J)} // edge weights
%apply (double* IN_ARRAY1, int DIM1){(double* V, int len_V)}

// direct MDS
%apply (double* IN_ARRAY1, int DIM1){(double* d, int len_d),
                                     (double* w, int len_w),
                                     (double* eta, int len_eta)}



extern void layout_unweighted(int n, double* X, int m, int* I, int* J, int t_max, double eps);
extern void layout_weighted(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps);
extern void layout_unweighted_convergent(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
extern void layout_weighted_convergent(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);
extern void layout_unweighted_focus(int n, double* X, int m, int* I, int* J, int f, int t_max, double eps, int t_maxmax);
extern void layout_weighted_focus(int n, double* X, int m, int* I, int* J, int f, double* V, int t_max, double eps, int t_maxmax);
extern void layout_unweighted_horizontal(int n, double* X, int m, int* I, int* J, int t_max, double eps, double delta, int t_maxmax);
extern void layout_weighted_horizontal(int n, double* X, int m, int* I, int* J, double* V, int t_max, double eps, double delta, int t_maxmax);
extern void mds_direct(int n, double* X, double* d, double* w, int t_max, double* eta);

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
%rename (layout_unweighted_focus) np_layout_unweighted_focus;
%exception np_layout_unweighted_focus {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_weighted_focus) np_layout_weighted_focus;
%exception np_layout_weighted_focus {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_unweighted_horizontal) np_layout_unweighted_horizontal;
%exception np_layout_unweighted_horizontal {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%rename (layout_weighted_horizontal) np_layout_weighted_horizontal;
%exception np_layout_weighted_horizontal {
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
    void np_layout_unweighted_focus(double* X, int n, int kd,
                                    int* I, int len_I,
                                    int* J, int len_J,
                                    int f,
                                    int t_max, double eps, int t_maxmax) {

        dimension_check(kd);
        unweighted_edge_check(len_I, len_J);
        layout_unweighted_focus(n, X, len_I, I, J, f, t_max, eps, t_maxmax);
    }
    void np_layout_weighted_focus(double* X, int n, int kd,
                                  int* I, int len_I,
                                  int* J, int len_J,
                                  int f,
                                  double* V, int len_V,
                                  int t_max, double eps, int t_maxmax) {

        dimension_check(kd);
        weighted_edge_check(len_I, len_J, len_V);
        layout_weighted_focus(n, X, len_I, I, J, f, V, t_max, eps, t_maxmax);
    }
    void np_layout_unweighted_horizontal(double* X, int n, int kd,
                                         int* I, int len_I,
                                         int* J, int len_J,
                                         int t_max, double eps, double delta, int t_maxmax) {

        dimension_check(kd);
        unweighted_edge_check(len_I, len_J);
        layout_unweighted_horizontal(n, X, len_I, I, J, t_max, eps, delta, t_maxmax);
    }
    void np_layout_weighted_horizontal(double* X, int n, int kd,
                                       int* I, int len_I,
                                       int* J, int len_J,
                                       double* V, int len_V,
                                       int t_max, double eps, double delta, int t_maxmax) {

        dimension_check(kd);
        weighted_edge_check(len_I, len_J, len_V);
        layout_weighted_horizontal(n, X, len_I, I, J, V, t_max, eps, delta, t_maxmax);
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


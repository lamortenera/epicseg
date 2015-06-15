/* This file gets automatically copied from your personal library to multiple projects.
To modify it change the version in /project/ale/home/data/cppstubs/array.hpp,
do not modify this version */
//timestamp: 2015-06-15 14:29:50.967004

#ifndef ARRAY_HPP
#define ARRAY_HPP

#include <Rcpp.h>
#include <R_ext/BLAS.h>

/* vector, matrix and cube classes.
 * All these classes do not own the memory (no alloc/dealloc), they
 * just wrap pointers to memory allocated/deallocated somewhere else */

//vector
template<typename T>
struct Vec {
    T* ptr;
    int len;
    
    Vec(){}
    
    Vec(T* _ptr, int _len):
        ptr(_ptr), len(_len){}
    
    inline T operator[] (int i) const {return ptr[i];}
    inline T& operator[] (int i){return ptr[i];}
    inline T* begin() {return ptr; }
    inline T* end() {return ptr + len;}
    
    Vec subset(int start, int end){
        return Vec(ptr + start, end - start);
    }
    
    Vec(SEXP x){
        const int RType = Rcpp::traits::r_sexptype_traits<T>::rtype;
        //no memory allocation
        if (TYPEOF(x) != RType) Rcpp::stop("incompatible types");

        Rcpp::Vector< RType > rcpp_vec(x);
        ptr = rcpp_vec.begin();
        len = rcpp_vec.length();
    }
};

//row of a matrix
template<typename T>
struct MatRow {
    T* ptr;
    int len;
    int step;
    
    MatRow(T* _ptr, int _len, int _step):
        ptr(_ptr), len(_len), step(_step){}
    
    inline int el(int i) { return i*step; }
    inline T operator[] (int i) const {return ptr[el(i)];}
    inline T& operator[] (int i){return ptr[el(i)];}
};


//normal matrix
template<typename T>
struct Mat {
    T* ptr;
    int nrow;
    int ncol;
    
    Mat(){}
    
    Mat(T* _ptr, int _nrow, int _ncol):
        ptr(_ptr), nrow(_nrow), ncol(_ncol){}
    
    inline int el(int row, int col){
        return row + nrow*col;
    }
    
    inline T operator[] (int i) const {return ptr[i];}
    inline T& operator[] (int i){return ptr[i];}
    
    inline T operator() (int row, int col) const {return ptr[el(row, col)];}
    inline T& operator() (int row, int col) {return ptr[el(row, col)];}
    
    inline T* colptr(int col){
        return ptr + el(0, col);
    }
    
    Mat subsetCol(int colStart, int colEnd){
        return Mat(colptr(colStart), nrow, colEnd-colStart);
    }
    
    inline Vec<T> getCol(int col){
        return Vec<T>(colptr(col), nrow);
    }
    
    inline MatRow<T> getRow(int row){
        return MatRow<T>(ptr+row, ncol, nrow);
    }
    
    Mat(SEXP x){
        const int RType = Rcpp::traits::r_sexptype_traits<T>::rtype;
        
        //no memory allocation
        if (TYPEOF(x) != RType) Rcpp::stop("incompatible types");
        Rcpp::Matrix< RType > rcpp_mat(x);
        ptr = rcpp_mat.begin();
        nrow = rcpp_mat.nrow();
        ncol = rcpp_mat.ncol();
    }
};

//3d array
template<typename T>
struct Cube {
    T* ptr;
    int nrow;
    int ncol;
    int nslice;
    
    Cube(){}
    
    Cube(T* _ptr, int _nrow, int _ncol, int _nslice):
        ptr(_ptr), nrow(_nrow), ncol(_ncol), nslice(_nslice){}
    
    inline int el(int row, int col, int slice){
        return row + nrow*(col + ncol*slice);
    }
    
    inline T operator[] (int i) const {return ptr[i];}
    inline T& operator[] (int i){return ptr[i];}
    
    inline T operator() (int row, int col, int slice) const {return ptr[el(row, col, slice)];}
    inline T& operator() (int row, int col, int slice) {return ptr[el(row, col, slice)];}
    
    inline T* sliceptr(int slice){
        return ptr + el(0, 0, slice);
    }
    
    Cube subsetSlice(int sliceStart, int sliceEnd){
        return Cube(sliceptr(sliceStart), nrow, ncol, sliceEnd-sliceStart);
    }
    
    inline Mat<T> getSlice(int slice){
        return Mat<T>(sliceptr(slice), nrow, ncol);
    }
    
    Cube(SEXP x){
        const int RType = Rcpp::traits::r_sexptype_traits<T>::rtype;
        
        //no memory allocation
        if (TYPEOF(x) != RType) Rcpp::stop("incompatible types");
        //need three dimensions
        Rcpp::IntegerVector dims = Rf_getAttrib(x, R_DimSymbol);
        if (dims.length() != 3) Rcpp::stop("not an array");
        Rcpp::Vector< RType > rcpp_vec(x);
        
        ptr = rcpp_vec.begin();
        nrow = dims[0];
        ncol = dims[1];
        nslice = dims[2];
    }
};



/* from std and rcpp datastructures extract pointers and wrap them */

//so far I am using it only to translate at compile-time INTSXP to int and REALSXP to double
#define CType(RType) typename Rcpp::traits::storage_type<RType>::type
//the opposite is:  Rcpp::traits::r_sexptype_traits<T>::rtype

//this is not going to work with vector<bool>,
//because it doesn't contain a pointer to bools
template<typename T>
inline Vec<T> asVec(std::vector<T>& v){
    return Vec<T>(v.data(), v.size());
}

template<int RType>
inline Vec< CType(RType) > asVec(Rcpp::Vector<RType>& v){
    return Vec< CType(RType) >(v.begin(), v.length());
}

template<typename T>
inline Mat<T> asMat(std::vector<T>& v, int ncol){
    if (v.size() % ncol != 0) throw std::invalid_argument("number of columns must be a divisor of vector length");
    return Mat<T>(v.data(), v.size()/ncol, ncol);
}

template<int RType>
inline Mat< CType(RType) > asMat(Rcpp::Matrix<RType>& m){
    return Mat< CType(RType) >(m.begin(), m.nrow(), m.ncol());
}

/* colsums and rowsums */
template<typename TNumMat, typename TNumVec, template <typename> class TMat>
static void colSums(TMat<TNumMat> mat, Vec<TNumVec> vec, int nthreads){
    if (mat.ncol != vec.len) throw std::invalid_argument("provided vector has invalid length");

    TNumVec*  cs = vec.ptr;
    int nrow = mat.nrow;
    int ncol = mat.ncol;
    
    #pragma omp parallel for schedule(static) num_threads(std::max(1, nthreads))
    for (int col = 0; col < ncol; ++col){
        TNumMat* ptr = mat.colptr(col);
        TNumMat tmp = 0;
        for (int row = 0; row < nrow; ++row){
            tmp += *ptr++;
        }
        cs[col] = tmp;
    }
}


//here vec must be clean at the beginning
template<typename TNumMat, typename TNumVec, template <typename> class TMat>
static void rowSums(TMat<TNumMat> mat, Vec<TNumVec> vec, int nthreads){
    if (mat.nrow != vec.len) throw std::invalid_argument("provided vector has invalid length");

    int nrow = mat.nrow;
    int ncol = mat.ncol;
    
    #pragma omp parallel num_threads(std::max(1, nthreads))
    {
        std::vector<TNumVec> acc(nrow, 0);
        TNumVec* accBegin = acc.data();
        #pragma omp for schedule(static) nowait
        for (int col = 0; col < ncol; ++col){
            TNumMat* matCol = mat.colptr(col);
            TNumVec* accIter = accBegin;
            for (int row = 0; row < nrow; ++row){//this loop should be unrolled...
                *accIter++ += *matCol++;
            }
        }
        #pragma omp critical
        {
            for (int row = 0; row < nrow; ++row){
                vec[row] += acc[row];
            }
        }
    }
}


inline void setDim(Mat<double> mat, bool trans, int* nrow, int* ncol, char* T){
    if (trans){
        *nrow = mat.ncol;
        *ncol = mat.nrow;
        T[0] = 'T';
        
    } else {
        *nrow = mat.nrow;
        *ncol = mat.ncol;
        T[0] = 'N';
    }
}

//C := alpha*op( A )*op( B ) + beta*C
inline void dgemm(Mat<double>& A, Mat<double>& B, Mat<double>& C, double alpha=1, double beta=0,
    bool transA=false, bool transB=false){
    
    int opA_nrow, opA_ncol; char TA[1];
    setDim(A, transA, &opA_nrow, &opA_ncol, TA);
    int opB_nrow, opB_ncol; char TB[1];
    setDim(B, transB, &opB_nrow, &opB_ncol, TB);
    
    int M = opA_nrow; 
    int N = opB_ncol;
    int K = opA_ncol;
    
    if (M != C.nrow || N != C.ncol || K != opB_nrow) {
        Rcpp::stop("dgemm: non conformable matrices");
    }
    
    F77_CALL(dgemm)(TA, TB, &M, &N, &K, &alpha, A.ptr, &A.nrow, 
                B.ptr, &B.nrow, &beta, C.ptr, &C.nrow);
}

//y := alpha*op(A)*x + beta*y
inline void dgemv(Mat<double> A, Vec<double> X, Vec<double> Y, double alpha=1, double beta=0,
    bool trans=false){
    
    //M = A.nrow, N = A.ncol
    char T[] = "N"; if (trans) T[0] = 'T';
    int one = 1;
    
    if ((!trans && (A.ncol != X.len || A.nrow != Y.len)) || 
        ( trans && (A.nrow != X.len || A.ncol != Y.len))){
        Rcpp::stop("dgemv: non conformable matrices");
    }
    
    F77_CALL(dgemv)(T, &A.nrow, &A.ncol, &alpha, A.ptr, &A.nrow, X.ptr, &one, 
            &beta, Y.ptr, &one);
}


#endif // ARRAY_HPP

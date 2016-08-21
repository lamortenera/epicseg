#include <Rcpp.h>
#include <random>
#include "tabbing.h"
#include "array.h"
#include "clist.h"



//t1 and t2 are vector<int>-like classes
template< class T1, class T2 >
inline void sortCounts(T1& source, T2& dest){
    if (source.len != dest.len) Rcpp::stop("non-matching arrays...");
    int len = source.len;
    if (len == 0) return;
    //get the maximum and minimum element
    int smin = source[0];
    int smax = source[0];
    for (int i = 1; i < len; ++i){
        if (source[i] > smax) smax = source[i];
        else if (source[i] < smin) smin = source[i];
    }
    if (smin < 0) Rcpp::stop("negative elements are not allowed");
    
    std::vector<int> tab(10);
    if (smax > 2*len) {//use std::sort
        //dest does not have iterators and cannot be used directly with std::sort
        //do everything with tab and copy back on dest
        tab.resize(len);
        std::copy(source.begin(), source.end(), tab.begin());
        std::sort(tab.begin(), tab.end());
        for (int i = 0; i < len; ++i){ dest[i] = tab[i]; }
        return;
    }
    //standard bucket-sort
    //tabulate all integers in source
    tabFast_impl(source.begin(), source.end(), tab);
    //do cumsum on the table
    for (int i = 0, e = tab.size(), acc = 0; i < e; ++i){
        int tmp = acc + tab[i];
        tab[i] = acc;
        acc = tmp;
    }
    //write destination
    for (int i = 0; i < len; ++i){
        int c = source[i];
        int pos = tab[c]++;
        dest[pos] = c;
    }
}


//note that vector v will be overwritten!!!
inline int median(Vec<int> v){
    int len = v.len;
    if (len == 1) return v[0];
    if (len == 2) return (v[0] + v[1])/2;
    
    //general case. I am not using partial sorting here, but I don't 
    //expect more than 100 elements, so it shouldn't make a big difference...
    std::sort(v.begin(), v.end());
    int mid = (len-1) >> 1;
    if (len & 1){//odd number
        return v[mid];
    } else {//even number
        return (v[mid] + v[mid+1]) >> 1;
    }
}

inline int mean(Vec<int> v){
    return std::accumulate(v.begin(), v.end(), 0)/v.len;
}
//not checking for 0-length vectors
inline int min(Vec<int> v){
    int res = v[0];
    for (int i = 1; i < v.len; ++i) if (v[i] < res) res = v[i];
    return res;
}

template < class Tfun >
inline void colSummary_loop(Mat<int> mat, Vec<int> ref, Tfun fun, int nthreads){
    int ncol = mat.ncol;
    #pragma omp parallel for num_threads(nthreads)
    for (int col = 0; col < ncol; ++col){
        ref[col] = fun(mat.getCol(col));
    }
}

// [[Rcpp::export]]
Rcpp::IntegerVector colSummary(Rcpp::IntegerMatrix mat, std::string type, int nthreads=1){
    //allocate final array
    Rcpp::IntegerVector ref(mat.ncol());
    Mat<int> smat = asMat(mat);
    Vec<int> sref = asVec(ref);
    
    //compute function 'type' for each column
    if (type == "median"){      colSummary_loop(smat, sref, median, nthreads);
    } else if (type == "mean"){ colSummary_loop(smat, sref, mean, nthreads);
    } else if (type == "min"){  colSummary_loop(smat, sref, min, nthreads);
    } else Rcpp::stop("invalid type");

    return ref;
}

// [[Rcpp::export]]
Rcpp::IntegerVector getRef(Rcpp::IntegerMatrix mat, std::string type, int nthreads=1){
    //allocate another matrix with transpose dimensions as mat
    int ncol = mat.ncol();
    int nrow = mat.nrow();
    if (ncol*nrow == 0) Rcpp::stop("empty input is invalid");
    Rcpp::IntegerMatrix mem_smat(ncol, nrow); 
    Mat<int> smat = asMat(mem_smat);
    Mat<int> omat = asMat(mat);
    //sort every column
    #pragma omp parallel for num_threads(nthreads)
    for (int col = 0; col < ncol; ++col){
        Vec<int> ovec = omat.getCol(col);
        MatRow<int> svec = smat.getRow(col);
        sortCounts(ovec, svec);
    }
    
    return colSummary(mem_smat, type, nthreads);
}

template < class Tgenerator >
inline void qtlnorm(Vec<int> source, Vec<int> ref, Vec<int> dest, 
    std::vector<std::pair<int, int> >& storage, Tgenerator& g){
    
    int len = source.len;
    if (len != ref.len || len != dest.len || len != storage.size()) {
        Rcpp::stop("incompatible vectors...");
    }
    //set the pairs
    for (int i = 0; i < len; ++i){
        storage[i].first = source[i];
        storage[i].second = i;
    }
    //sort
    std::sort(storage.begin(), storage.end());
    //shuffle draws
    //find intervals [start, end) such that all values have the same score
    int start = 0;
    while (start < len-1){
        if (storage[start].first != storage[start+1].first){
            ++start; continue;
        }
        //beginning of a draw at index start, find end
        double draw = storage[start].first;
        int end = start + 2;
        for (; end < len; ++end){
            if (storage[end].first != draw) break;
        }
        std::shuffle(storage.begin() + start, storage.begin() + end, g);
        start = end;
    }
    
    //assign reference values
    for (int i = 0; i < len; ++i){
        dest[storage[i].second] = ref[i];
    }
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix quantileNorm(Rcpp::IntegerMatrix mat, Rcpp::IntegerVector ref, int nthreads=1, int seed=13){
    if (mat.nrow() != ref.length()) Rcpp::stop("incompatible arrays...");
    if (!std::is_sorted(ref.begin(), ref.end())) Rcpp::stop("ref must be sorted");
    int ncol = mat.ncol();
    int nrow = mat.nrow();
    //allocate new matrix
    Rcpp::IntegerMatrix res(nrow, ncol);
    Mat<int> oldmat = asMat(mat); 
    Mat<int> newmat = asMat(res);
    Vec<int> ref2 = asVec(ref);
    //allocate a seed for each column
    std::seed_seq sseq{seed};
    std::vector<std::uint32_t> seeds(ncol);
    sseq.generate(seeds.begin(), seeds.end());
    
    #pragma omp parallel num_threads(nthreads)
    {
        std::vector<std::pair<int, int> > storage(nrow);//pairs <value, index>
        #pragma omp for 
        for (int col = 0; col < ncol; ++col){
            std::mt19937 gen(seeds[col]);
            qtlnorm(oldmat.getCol(col), ref2, newmat.getCol(col), storage, gen);
        }
    }
    
    res.attr("dimnames") = mat.attr("dimnames");
    return res;
}

// [[Rcpp::export]]
Rcpp::List clist2mlist(Rcpp::List clist, int nthreads=1){
    if (clist.length()==0) Rcpp::stop("empty list is invalid");
    int ncounts, nmarks, nbins = -1;
    std::vector<std::string> rnames;
    listcubedim(clist, &nmarks, &nbins, &ncounts, rnames);
    
    //allocate storage
    Rcpp::List mlist(nmarks);
    for (int mark = 0; mark < nmarks; ++mark){
        mlist[mark] = Rcpp::IntegerMatrix(nbins, ncounts);
    }
    if (rnames.size() > 0) mlist.attr("names") = rnames;
    
    //copy data
    #pragma omp parallel for num_threads(nthreads) collapse(2) 
    for (int c = 0; c < ncounts; ++c){
        for (int mark = 0; mark < nmarks; ++mark){
            MatRow<int> row = Mat<int>((SEXP)clist[c]).getRow(mark);
            Vec<int> col = Mat<int>((SEXP)mlist[mark]).getCol(c);
            for (int bin = 0; bin < nbins; ++bin){
                col[bin] = row[bin];
            }
        }
    }
    
    return mlist;
}

// [[Rcpp::export]]
Rcpp::List mlist2clist(Rcpp::List mlist, int nthreads=1){
    if (mlist.length()==0) Rcpp::stop("empty list is invalid");
    int ncounts, nmarks, nbins = -1;
    std::vector<std::string> foo;
    listcubedim(mlist, &nbins, &ncounts, &nmarks, foo);
    Rcpp::List newdnames(2); newdnames[0] = mlist.attr("names");
    
    //allocate storage
    Rcpp::List clist(ncounts);
    for (int c = 0; c < ncounts; ++c){
        Rcpp::IntegerMatrix mat(nmarks, nbins);
        if (!Rf_isNull(newdnames[0])) mat.attr("dimnames") = newdnames;
        clist[c] = mat;
    }
    
    //copy data
    #pragma omp parallel for num_threads(nthreads) collapse(2) 
    for (int c = 0; c < ncounts; ++c){
        for (int mark = 0; mark < nmarks; ++mark){
            Vec<int> col = Mat<int>((SEXP)mlist[mark]).getCol(c);
            MatRow<int> row = Mat<int>((SEXP)clist[c]).getRow(mark);
            for (int bin = 0; bin < nbins; ++bin){
                row[bin] = col[bin];
            }
        }
    }
    
    return clist;
}
//just for testing

// [[Rcpp::export]]
Rcpp::IntegerVector testSortCounts(Rcpp::IntegerVector v){
    Rcpp::IntegerVector res(v.length());
    Vec<int> v2 = asVec(v);
    Vec<int> v3 = asVec(res);
    sortCounts(v2, v3);
    return res;
}

// [[Rcpp::export]]
int testMeanAndMedian(Rcpp::IntegerVector v, std::string type){
    if (type=="mean") return mean(asVec(v));
    if (type=="median") return median(asVec(v));
    Rcpp::stop("invalid type");
    return -1;
}

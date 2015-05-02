#include <Rcpp.h>
#include "array.hpp"

//it also checks that s is a matrix
inline std::vector<std::string> getRownames(Rcpp::IntegerMatrix s){
    SEXP dm = s.attr("dimnames");
    if (Rf_isNull(dm)) return std::vector<std::string>();
    Rcpp::List dnames = Rcpp::as<Rcpp::List>(dm);
    if (Rf_isNull(dnames[0])) return std::vector<std::string>();
    return dnames[0];
}

inline void listcubedim(Rcpp::List clist, int* nrow, int* ncol, int* nslice, 
                                            std::vector<std::string>& rnames){
    if (clist.length()==0) Rcpp::stop("empty list is invalid");
    *nslice = clist.length();
    //check the dimension and the rownames of the first matrix
    {
        Mat<int> counts((SEXP)clist[0]);
        *nrow = counts.nrow;
        *ncol = counts.ncol;
        rnames = getRownames(clist[0]);
    }
    //check strict compatibility with the other matrices
    for (int i = 1, e = *nslice; i < e; ++i){
        Mat<int> counts((SEXP)clist[0]);
        if (*nrow != counts.nrow || *ncol != counts.ncol) {
            Rcpp::stop("matrices have incompatible dimensions!"); }
        std::vector<std::string> tmp(getRownames(clist[0]));
        if (tmp.size() != rnames.size()) Rcpp::stop("inconsistent naming");
        for (unsigned j = 0; j < rnames.size(); ++j){
            if (tmp[j] != rnames[j]) {
                Rcpp::stop("rownames don't match!");}
        }
    }
    
}


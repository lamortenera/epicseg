#include <Rcpp.h>
#include "array.h"

//it also checks that s is a matrix
inline std::vector<std::string> getRownames(Rcpp::RObject s){
    SEXP dm = s.attr("dimnames");
    if (Rf_isNull(dm)) return std::vector<std::string>();
    Rcpp::List dnames = Rcpp::as<Rcpp::List>(dm);
    if (Rf_isNull(dnames[0])) return std::vector<std::string>();
    return dnames[0];
}

inline std::vector<int> getDim(Rcpp::RObject s){
    return Rcpp::as<std::vector<int> >(s.attr("dim"));
}

inline void listcubedim(Rcpp::List clist, int* nrow, int* ncol, int* nslice, 
                                            std::vector<std::string>& rnames){
    if (clist.length()==0) Rcpp::stop("empty list is invalid");
    *nslice = clist.length();
    //check the dimension and the rownames of the first matrix
    {
        std::vector<int> dims = getDim(clist[0]);
        if (dims.size() != 2) Rcpp::stop("not a matrix!");
        *nrow = dims[0];
        *ncol = dims[1];
        rnames = getRownames(clist[0]);
    }
    //check strict compatibility with the other matrices
    for (int i = 1, e = *nslice; i < e; ++i){
        std::vector<int> dims = getDim(clist[0]);
        if (dims.size() != 2) Rcpp::stop("not a matrix!");
        if (*nrow != dims[0] || *ncol != dims[1]) {
            Rcpp::stop("matrices have incompatible dimensions!"); }
        std::vector<std::string> tmp(getRownames(clist[0]));
        if (tmp.size() != rnames.size()) Rcpp::stop("inconsistent naming");
        for (unsigned j = 0; j < rnames.size(); ++j){
            if (tmp[j] != rnames[j]) {
                Rcpp::stop("rownames don't match!");}
        }
    }
    
}


/* This file gets automatically copied from your personal library to multiple projects.
To modify it change the version in /project/ale/home/data/cppstubs/tabbing.hpp,
do not modify this version */
//timestamp: 2015-06-15 10:46:21.128639

//functions for tabulating the values of an array of integers
inline void shrink(std::vector<int>& v){
    int realsize = v.size();
    for (; realsize > 0 && v[realsize-1] == 0; --realsize){}
    v.resize(realsize);
}


template<typename TIter>
void tabFast_impl(TIter C, TIter E, std::vector<int>& tab, bool wantshrink=true){
    int csize = tab.size();
    for (; C < E; ++C){
        int c = *C;
        if (c < 0) Rcpp::stop("negative counts are not allowed");
        if (c >= csize) { csize += c; tab.resize(csize); }
        ++tab[c];
    }
    if (wantshrink) shrink(tab);
}

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "clist.h"
#include <R_ext/Applic.h>

struct Edge {
    double weight;
    int i;
    int j;
};

static inline bool edgeComparator(const Edge& a, const Edge& b){
    return b.weight > a.weight;
}

struct DisjointSets {
    struct Node {
        Node* parent;
        int rank;
    };
    
    std::vector<Node> nodes;
    
    DisjointSets(int n) : nodes(n){
        for (int i = 0; i < n; ++i){
            nodes[i].parent = &nodes[i];
            nodes[i].rank = 0;
        }
    }
    
    Node* findSet(Node* x){
        if (x->parent != x){
            x->parent = findSet(x->parent);
        }
        return x->parent;
    }
    
    void link(Node* x, Node* y){
        if (x->rank > y->rank){
            y->parent = x;
        } else {
            x->parent = y;
            if (x->rank == y->rank){
                ++(y->rank);
            }
        }
    }
    
    void merge(int i, int j){
        link(findSet(&nodes[i]), findSet(&nodes[j]));
    }
    
    bool sameSet(int i, int j){
        return findSet(&nodes[i]) == findSet(&nodes[j]);
    }
    
};

//[[Rcpp::export]]
Rcpp::IntegerVector smallWeightHamiltonianPath(Rcpp::NumericMatrix dmat){
    if (dmat.ncol() != dmat.nrow()) Rcpp::stop("'dmat' must be a square matrix");
    int n = dmat.ncol();
    if (n <= 1) return Rcpp::IntegerVector(n, n);
    //adjacency list representation of the graph
    std::vector<std::vector<int> > graph(n);
    //connected components of the graph
    DisjointSets ccs(n);
    //all possible edges
    int nedges = n*(n-1)/2;
    std::vector<Edge> edges(nedges);
    {
        int p = 0;
        for (int i = 0; i < n-1; ++i){
            for (int j = i+1; j < n; ++j){
                edges[p].weight = dmat(i,j);
                edges[p].i = i;
                edges[p].j = j;
                ++p;
            }
        }
    }
    //sort the edges
    std::sort(edges.begin(), edges.end(), edgeComparator);
    //add the edges from the smallest to the largest with
    //the constraint that:
    //1. an edge must connect two disconnected components
    //2. the nodes being connected must have degree less than 2
    int m = 0;
    for (int p = 0; p < nedges && m < n-1 ; ++p){
        Edge& e = edges[p];
        if (!ccs.sameSet(e.i, e.j) && graph[e.i].size() < 2 && graph[e.j].size() < 2){
            ccs.merge(e.i, e.j);
            graph[e.i].push_back(e.j);
            graph[e.j].push_back(e.i);
            ++m;
        }
    }
    if (m != n-1) Rcpp::stop("didn't manage to find n-1 edges...");
    //find one of the two ends of the graph
    int init = -1;
    for (int i = 0; i < n && init < 0; ++i){
        if (graph[i].size() == 1) init = i;
    }
    if (init < 0) Rcpp::stop("didn't manage to find a node with degree 1...");
    //reconstruct the path
    Rcpp::IntegerVector path(n);
    int last = init;
    int curr = graph[last][0];
    int nextPos = 2;
    path[0] = last;
    path[1] = curr;
    while (graph[curr].size() > 1 && nextPos < n){
        int next = graph[curr][0];
        if (next == last) next = graph[curr][1];
        path[nextPos++] = next;
        last = curr;
        curr = next;
    }
    if (nextPos < n || graph[curr].size() != 1) Rcpp::stop("something wrong with the other end of the path...");
    //convert the C indices to R indices
    for (int i = 0; i < n; ++i) ++path[i];
    
    return path;
}

inline void writeVectors(Rcpp::IntegerVector bigv, Rcpp::List vlist, int nthreads){
    int nv = vlist.length();
    std::vector<long int> breaks(nv + 1);
    for (long int i = 0, acc = 0; i < nv; ++i){
        acc += Rcpp::as<Rcpp::IntegerVector>(vlist[i]).length();
        breaks[i+1] = acc;
    }
    if (breaks[nv] != bigv.length()) Rcpp::stop("invalid length");
    
    #pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < nv; ++i){
        //if nthreads >> vlist.length() you might want to parallelize better...
        Rcpp::IntegerVector v = vlist[i];
        memcpy(bigv.begin() + breaks[i], v.begin(), sizeof(int)*v.length());
    }
}

//these 2 functions are unsafe for the following reasons:
//1. if obj is referenced by more than 1 variable, all variables
//will be modified (something that never happens in R)
//2. in setDim, if obj has some dimnames and the new dims are the same as the 
// old ones the dimnames will be lost for no reason
//3. there might be other reasons that I don't see
//of course I wrote them because obj can be a huge vector that I don't want
//to copy

// [[Rcpp::export]]
void setDim_unsafe(Rcpp::RObject obj, Rcpp::IntegerVector dims){
    obj.attr("dim") = dims;
}

// [[Rcpp::export]]
void setDimnames_unsafe(Rcpp::RObject obj, Rcpp::List dimnames){
    obj.attr("dimnames") = dimnames;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix bindCols(Rcpp::List vlist, int nthreads=1){
    int ncol = vlist.length();
    int nrow;
    {
        Rcpp::IntegerVector v = vlist[0];
        nrow = v.length();
    }
    for (int i = 1; i < ncol; ++i){
        Rcpp::IntegerVector v = vlist[i];
        if (v.length() != nrow) Rcpp::stop("The vectors in the list must have equal length");
    }
    
    Rcpp::IntegerMatrix mat(nrow, ncol);
    writeVectors(mat, vlist, nthreads);
    
    Rcpp::List dimnames(2);
    dimnames[1] = vlist.attr("names");
    mat.attr("dimnames") = dimnames;
    
    return mat;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix bindCList(Rcpp::List clist, int nthreads=1){
    if (clist.length()==0) Rcpp::stop("empty list is invalid");
    int ncounts, nmarks, nbins = -1;
    std::vector<std::string> rnames;
    listcubedim(clist, &nmarks, &nbins, &ncounts, rnames);
    
    Rcpp::IntegerMatrix bigmat(nmarks, ((long int)nbins)*ncounts);
    writeVectors(bigmat, clist, nthreads);
    
    Rcpp::List dnames = Rcpp::List(2); 
    dnames[0] = rnames;
    setDimnames_unsafe(bigmat, dnames);
    
    
    return bigmat;
}

// [[Rcpp::export]]
void writeCountsTXT(Rcpp::IntegerMatrix counts, std::vector<std::string> marks, std::string path){
    if (counts.nrow() != ((int)marks.size())) Rcpp::stop("rownames don't match with counts");
    if (counts.length() == 0) Rcpp::stop("empty count matrix");
    std::ofstream outstrm; outstrm.open(path);
    
    int nmarks = marks.size();
    int ncol = counts.ncol();
    //write marks
    outstrm << marks[0];
    for (int i = 1; i < nmarks; ++i) outstrm << "\t" << marks[i];
    outstrm << "\n";
    
    //write counts
    for (int col = 0; col < ncol; ++col){
        Rcpp::MatrixColumn<INTSXP> C = counts.column(col);
        
        outstrm << C[0];
        for (int i = 1; i < nmarks; ++i) outstrm << "\t" << C[i];
        outstrm << "\n";
    }
    
    outstrm.close();
}


// [[Rcpp::export]]
Rcpp::NumericMatrix avgCountsPerClust(Rcpp::IntegerMatrix counts, Rcpp::IntegerVector clusts){
    if (counts.ncol() != clusts.length()) Rcpp::stop("invalid input");
    int len = clusts.length();
    int nmarks = counts.nrow();
    int nclusts = 0;
    for (int i = 0; i < len; ++i){
        if (clusts[i] <= 0) Rcpp::stop("cluster labels must be from 1 to nclust (included)");
        if (clusts[i] > nclusts) nclusts = clusts[i];
    }
    
    Rcpp::NumericMatrix avgs(nmarks, nclusts);
    std::vector<int> denoms(nclusts);
    for (int i = 0; i < len; ++i){
        int c = clusts[i]-1;
        ++denoms[c];
        Rcpp::MatrixColumn<INTSXP> C = counts.column(i);
        Rcpp::MatrixColumn<REALSXP> A = avgs.column(c);
        for (int j = 0; j < nmarks; ++j){
            A[j] += C[j];
        }
    }
    
    for (int c = 0; c < nclusts; ++c){
        double denom = denoms[c];
        Rcpp::MatrixColumn<REALSXP> A = avgs.column(c);
        for (int j = 0; j < nmarks; ++j){
            A[j] /= denom;
        }
    }
    
    return avgs;
}

//functions for tabulating one or two factors at a decent speed

// [[Rcpp::export]]
Rcpp::IntegerVector tabf(Rcpp::IntegerVector v, bool naRm=true){
    if (!Rf_inherits(v, "factor")) Rcpp::stop("expecting a factor");
    Rcpp::CharacterVector levels = v.attr("levels");
    int nlevels = levels.length();
    
    Rcpp::IntegerVector ans(nlevels);
    bool outofrange = false;
    for (int i = 0, e = v.length(); i < e; ++i){
        if (v[i] <= 0 || v[i] > nlevels) {
            outofrange = true;
        } else ++ans[v[i]-1];
    }
    
    if (outofrange && !naRm) Rcpp::stop("there were NAs or values out of range");
    
    ans.attr("names") = levels;
    return ans;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix tabf2(Rcpp::IntegerVector v1, Rcpp::IntegerVector v2, bool naRm=true){
    if (!Rf_inherits(v1, "factor") || !Rf_inherits(v2, "factor")) Rcpp::stop("expecting factors");
    if (v1.length() != v2.length()) Rcpp::stop("the two factors must have the same length");
    Rcpp::CharacterVector levels1 = v1.attr("levels");
    Rcpp::CharacterVector levels2 = v2.attr("levels");
    int nlevels1 = levels1.length();
    int nlevels2 = levels2.length();
    
    Rcpp::IntegerMatrix ans(nlevels1, nlevels2);
    bool outofrange = false;
    for (int i = 0, e = v1.length(); i < e; ++i){
        if (v1[i] <= 0 || v1[i] > nlevels1 || 
            v2[i] <= 0 || v2[i] > nlevels2) {
            outofrange = true;
        } else  ++ans(v1[i]-1, v2[i]-1);
    }
    
    if (outofrange && !naRm) Rcpp::stop("there were NAs or values out of range");
    
    ans.attr("dimnames") = Rcpp::List::create(levels1, levels2);
    return ans;
}

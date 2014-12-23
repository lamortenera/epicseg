#include <Rcpp.h>
#include <iostream>
#include <fstream>


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
	std::vector<std::vector<int>> graph(n);
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


// [[Rcpp::export]]
Rcpp::IntegerMatrix bindCols(Rcpp::List vlist){
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
	for (int i = 0; i < ncol; ++i){
		Rcpp::IntegerVector v = vlist[i];
		memcpy(mat.begin() + nrow*i, v.begin(), sizeof(int)*nrow);
	}
	
	Rcpp::List dimnames(2);
	dimnames[1] = vlist.attr("names");
	mat.attr("dimnames") = dimnames;
	
	return mat;
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
#include <Rcpp.h>
#include <iostream>
#include <fstream>


//iterator through run-length-encoded values
class RleIter {
    Rcpp::IntegerVector rlens;
    Rcpp::IntegerVector values;
    std::vector<std::string> names;
    
    
    public:
    
        int run;
        int rlen;
        int rpos;
        bool valid;

        RleIter(Rcpp::RObject rle):
            rlens(Rcpp::as<Rcpp::IntegerVector>(rle.slot("lengths"))),
            values(Rcpp::as<Rcpp::IntegerVector>(rle.slot("values"))),
            names(Rcpp::as<std::vector<std::string> >(values.attr("levels"))),
            run(0), rpos(-1)
        {
            next();
        }
        
        bool next(){
            ++rpos;
            if (rpos == rlens[run]){ //end of the run, go to the next
                ++run; rpos = 0;
                if (run == rlens.length())
                    valid = false;
                    return valid;
            }
            valid = true;
            return valid;
        }
        
        std::string& getValue(){
            return names[values[run]-1];
        }
};

void printSegment(std::ofstream& out, std::string& chr, int start, int end, std::string& label, std::string& color){
    out << chr << "\t" << start << "\t" << end << "\t" << label << "\t0\t.\t" << start << "\t" << end << "\t" << color << "\n";
} 

//segment in bed format
struct Segment {
    std::string chr;
    int start; //0-indexed, inclusive
    int end; //0-indexed, exclusive
    int state; //from an R vector, so 1-indexed
    
    Segment(std::string& _chr, int _start, int _end, int _state) : chr(_chr) {
        start = _start;
        end = _end;
        state = _state;
    }
};

std::vector<Segment> getSegments(Rcpp::RObject gr, Rcpp::IntegerVector states){
    //argument checking
    if (not gr.inherits("GRanges")) Rcpp::stop("must provide a GRanges object");
    Rcpp::IntegerVector starts = Rcpp::as<Rcpp::IntegerVector>(Rcpp::as<Rcpp::RObject>(gr.slot("ranges")).slot("start"));
    Rcpp::IntegerVector lens =   Rcpp::as<Rcpp::IntegerVector>(Rcpp::as<Rcpp::RObject>(gr.slot("ranges")).slot("width"));
    int nbins = states.length();
    int nreg = lens.length();
    long totlen = 0; 
    for (int i = 0; i < nreg; ++i) totlen += lens[i];
    if (totlen % nbins != 0) Rcpp::stop("total genomic span is not a multiple of the number of states");
    int binsize = totlen / nbins;
    for (int i = 0; i < nreg; ++i) if (lens[i] % binsize != 0) Rcpp::stop("region lengths must be multiple of binsize");
    
    //iterates through the run-length-encoded chromosome names
    RleIter chrs(Rcpp::as<Rcpp::RObject>(gr.slot("seqnames")));
    
    //accumulate segments here
    std::vector<Segment> segments;
    int initSize = states.length()/4;
    segments.reserve(std::max(10, initSize));
    
    
    int stateidx = 0;
    //loop through each region
    for (int i = 0, e = lens.length(); i < e; ++i, chrs.next()){
        int rangestart = starts[i] - 1; //starts are one-indexed when coming from the GRanges format
        int rangebins = lens[i]/binsize;
        
        if (rangebins>0) {
            int segstart = rangestart; 
            int state = states[stateidx];
            //loop through the states in that region
            for (int j = 1; j < rangebins; ++j){
                int currstate = states[stateidx + j];
                if (currstate != state){
                    int segend = rangestart + j*binsize;
                    segments.push_back(Segment(chrs.getValue(), segstart, segend, state));
                    segstart = segend;
                    state = currstate;
                }
            }
            segments.push_back(Segment(chrs.getValue(), segstart, rangestart + rangebins*binsize, state));
            stateidx += rangebins;
        }
    }
    
    return segments;
}


std::vector<int> segmentNamesToInts(Rcpp::RObject segments, int maxInt){
    Rcpp::RObject ranges = Rcpp::as<Rcpp::RObject>(segments.slot("ranges"));
    Rcpp::CharacterVector names = Rcpp::as<Rcpp::CharacterVector>(ranges.slot("NAMES"));
    int len = names.length();
    std::vector<int> states(len);
    for (int i = 0; i < len; ++i){
        int state = atoi(names[i]);
        if (state <= 0 || state > maxInt) Rcpp::stop("names in the ranges object must be numbers from 1 to n (number of clusters)");
        states[i] = state;
    }
    return states;
}

// [[Rcpp::export]]
Rcpp::List statesToSegments_helper(Rcpp::RObject regions, Rcpp::IntegerVector states){
    std::vector<Segment> segments = getSegments(regions, states);
    int nsegm = segments.size();
    std::vector<std::string> chrs(nsegm);
    Rcpp::IntegerVector starts(nsegm);
    Rcpp::IntegerVector ends(nsegm);
    Rcpp::IntegerVector segstates(nsegm);
    for (int i = 0; i < nsegm; ++i){
        Segment& s = segments[i];
        chrs[i] = s.chr;
        starts[i] = s.start + 1;
        ends[i] = s.end;
        segstates[i] = s.state;
    }
    
    return Rcpp::List::create(Rcpp::Named("chrs")=chrs, Rcpp::Named("starts")=starts,
                    Rcpp::Named("ends")=ends, Rcpp::Named("states")=segstates);
}

// [[Rcpp::export]]
void segmentsToBed(Rcpp::RObject segments, std::vector<std::string> labels, std::vector<std::string> colors, std::string path){
    //argument checking
    if (colors.size() != labels.size()) Rcpp::stop("'labels' doens't match with 'colors'");
    int n = colors.size();
    if (not segments.inherits("GRanges")) Rcpp::stop("must provide a GRanges object");
    RleIter chrs(Rcpp::as<Rcpp::RObject>(segments.slot("seqnames")));
    Rcpp::RObject ranges = Rcpp::as<Rcpp::RObject>(segments.slot("ranges"));
    Rcpp::IntegerVector starts = Rcpp::as<Rcpp::IntegerVector>(ranges.slot("start"));
    Rcpp::IntegerVector lens =   Rcpp::as<Rcpp::IntegerVector>(ranges.slot("width"));
    std::vector<int> states = segmentNamesToInts(segments, n);
    int nsegm = starts.length();
    //open stream
    std::ofstream out; out.open(path);
    //print each segment
    for (int i = 0; i < nsegm; ++i, chrs.next()){
        int state = states[i]-1;
        int start = starts[i]-1;
        int end = start+lens[i];
        out << chrs.getValue() << "\t" << start << "\t" << end << "\t" << labels[state];
        out << "\t0\t+\t" << start << "\t" << end << "\t" << colors[state] << "\n";
    }
    //close stream
    out.close();
}

/*
// [[Rcpp::export]]
Rcpp::IntegerVector segmentsToStates(Rcpp::RObject segments, int nstates, int binsize){
    //argument checking
    if (not segments.inherits("GRanges")) Rcpp::stop("must provide a GRanges object");
    Rcpp::RObject ranges = Rcpp::as<Rcpp::RObject>(segments.slot("ranges"));
    Rcpp::IntegerVector lens =   Rcpp::as<Rcpp::IntegerVector>(ranges.slot("width"));
    std::vector<int> states = segmentNamesToInts(segments, nstates);
    
    std::vector
    
}
*/

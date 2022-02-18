/*
A KBase module: Snekmer
*/

module Snekmer {
    /*
    Run Model capabilities of Snekmer.

    workspace_name - the name of the workspace for input/output
    kmer - kmer length for features
    alphabet - mapping function for reduced amino acid sequences
    min_rep_thresh - min number of sequences to include feature for prefiltering
    processes - for parallelization

    */

    typedef structure {
        string workspace_name;
        int kmer;
        string alphabet;
        float min_rep_thresh;
        int processes;
    } SnekmerModelParams;

    typedef structure {
        string report_name;
        string report_ref;
    } SnekmerModelOutput;

    /*
        run_Snekmer_model accepts some of the model params for now, and returns results in a KBaseReport
    */
    funcdef run_Snekmer_model(SnekmerModelParams params) returns (SnekmerModelOutput output) authentication required;

};

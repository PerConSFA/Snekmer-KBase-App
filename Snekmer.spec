/*
A KBase module: Snekmer
*/

module Snekmer {
    /*
        Reference to a Genome object in the workspace
        @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;

    /*
    Input parameters for Snekmer Model.

    workspace_name - the name of the workspace for input/output
    k - kmer length for features
    alphabet - mapping function for reduced amino acid sequences
    min_rep_thresh - min number of sequences to include feature for prefiltering
    processes - for parallelization

    */

    typedef structure {
        string workspace_name;
        int k;
        string alphabet;
        float min_rep_thresh;
        int processes;
    } SnekmerModelParams;

    /*
    Output parameters for Snekmer Model.

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.

    */
    typedef structure {
        string report_name;
        string report_ref;
    } SnekmerModelOutput;

    /*
    Input parameters for Snekmer Search.

    workspace_name - the name of the workspace for input/output
    object_ref - GenomeSet
    k - kmer length for features
    alphabet - mapping function for reduced amino acid sequences
    output_genome_name - output object name
    */
    typedef string obj_ref;

    typedef structure {
        string workspace_name;
        string object_ref;
        int k;
        int alphabet;
        string output_genome_name;
    } SnekmerSearchParams;

    /*
    Output parameters for Snekmer Search.

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.

    */
    typedef structure {
        string report_name;
        string report_ref;
        genome_ref output_genome_ref;
    } SnekmerSearchOutput;


    /*
        run_Snekmer_model accepts some of the model params for now, and returns results in a KBaseReport
    */
    funcdef run_Snekmer_model(SnekmerModelParams params) returns (SnekmerModelOutput output) authentication required;

    /*
        run_Snekmer_search accepts some of the search params for now, and returns results in a KBaseReport
    */
    funcdef run_Snekmer_search(SnekmerSearchParams params) returns (SnekmerSearchOutput output) authentication required;
};

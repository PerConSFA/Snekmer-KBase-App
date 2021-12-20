/*
A KBase module: Snekmer
*/

module Snekmer {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_Snekmer(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};

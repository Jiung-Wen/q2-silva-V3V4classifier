#!/usr/bin/env python3

import qiime2
import argparse
from qiime2 import Artifact
from qiime2.plugins.feature_classifier.methods import classify_sklearn, \
                                                        extract_reads, \
                                                        fit_classifier_naive_bayes


def get_args():
    parser = argparse.ArgumentParser(description="Flexibly train qiime2 classifiers based on the 132 silva dataset.")
    parser.add_argument("--gene", help="Specify 16S or 18S, defaults to all (both genes).", default="all", choices=["16S", "18S", "all"], type=str, required=False)
    parser.add_argument("--perc_id", help="Specify a percentage id.", default="99", choices=["99", "97", "94", "90"], type=str, required=False)
    parser.add_argument("--f_primer", help="Specify a forward primer.", default="CCTACGGGNGGCWGCAG", type=str, required=False)
    parser.add_argument("--r_primer", help="Specify a reverse primer.", default="GACTACHVGGGTATCTAATCC", type=str,  required=False)
    args = parser.parse_args()
    return args

def get_file_loc(gene, perc_id):
    if gene == "16S" :
        seq_loc= "rep_set/rep_set_16S_only/" + perc_id + "/silva132_" + perc_id + "_16S.fna"
        tax_loc =  "taxonomy/16S_only/" + perc_id + "/majority_taxonomy_7_levels.txt"
        return seq_loc, tax_loc
    elif gene == "18S" :
        seq_loc= "rep_set/rep_set_18S_only/" + perc_id + "/silva132_" + perc_id + "_18S.fna"
        tax_loc =  "taxonomy/18S_only/" + perc_id + "/majority_taxonomy_7_levels.txt"
        return seq_loc, tax_loc
    elif gene == "all" :
        seq_loc= "rep_set/rep_set_all/" + perc_id + "/silva132_" + perc_id + ".fna"
        tax_loc =  "taxonomy/taxonomy_all/" + perc_id + "/majority_taxonomy_7_levels.txt"
        return seq_loc, tax_loc

def main():
    # Parse arguments
    options = get_args()

    seq_loc, tax_loc = get_file_loc(options.gene, options.perc_id)

    # Training feature classifiers with q2-feature-classifier
    # https://docs.qiime2.org/2019.1/tutorials/feature-classifier/
    silva_132 = Artifact.import_data('FeatureData[Sequence]', seq_loc )

    silva_132_taxonomy = Artifact.import_data('FeatureData[Taxonomy]', tax_loc,
                                           view_type = 'HeaderlessTSVTaxonomyFormat')
   
    # extract reference reads
    # V3-V4: 341f: CCTACGGGNGGCWGCAG; 806r: GACTACHVGGGTATCTAATCC
    ref_seqs_s = extract_reads(sequences = silva_132,
                           f_primer = options.f_primer,
                           r_primer = options.r_primer)

    # train the classifier
    silva_classifier = fit_classifier_naive_bayes(reference_reads = ref_seqs_s.reads,
                                              reference_taxonomy = silva_132_taxonomy)

    # save the classifier
    silva_classifier.classifier.save("silva132_" + options.perc_id + "_v3v4_" + options.gene)

if __name == "__main__":
    main()

# Bag-of-Segments
The regression-tree based segmentation for CNAs reported by LC-WGS, visualization and feature generation.

# input
The output file of the Wandy software available at: http://bioinformaticstools.mayo.edu/research/wandy

required columns:

is.reliable.bin: 0, 1

Chr: chr1, chr2,..., chrX

corrected.count: numeric

# output:
segment_info: a data frame for the detailed information of the segments

fullseq_plot: the segmentation full sequencing plot, when argument sequence_plot = TRUE

segment_2d_plot: 2d with-height plot

BOS_features, width_features, height_features: Bag-of-Segments features and reduced features, defined by the width and height quantiles (arguments width_thre, height_thre)


The R project contains the developed code for the differential analysis of some date collected from NGS regarding breast cancer. Before executing the R code, we ran the data through different programs such as:
- FastQC: Creates different graphs and did some brief data analysis;
- Trimmomatic: Cleans the primers of the existing reads from the results of FastQC;
- BWA aligner: Aligns the resulting reads with a reference genome;
- HTSeq (Python package): Aligns the best results of the different steps with the gtf file of the human being.

Due to the fact that these programs were executed in a SVN repository, the report of this phase contains the main commands used during this implementation.

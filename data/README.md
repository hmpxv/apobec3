## Data files for analysis

### Note:
XML files contain partitions for all taxa in analysis except EPI_ISL_13983888 because that must be sourced from GISAID directly and the partition then constructed manually. 
We have included a make_partitions.py script (in the scripts directory of this repository) that can be modified to construct a partition on a custom alignment. 
Please bear in mind that if the user wishes to slot the new partition into the rest of the xml the sequences must be aligned using squirrel as the alignment needs to have sequences all of the same length.
Alternatively the user may wish to comment out the EPI_ISL_13983888 taxon name and run the xml files as they are if they would rather work independent of GISAID.

# Workflow 
## File Read in
- read in `readcounts` file -> `ReadCounts` object to store relevant count data and relevant Transcript IDs
- Read in GTF - only safe transcripts that are relevant (skip others)
- create transcripts from transcript object exon regions and `fasta`, `fidx`
- iterate through all entries in `ReadCounts` and run parallel read creation


## Read Generation
- Draw random value with passed parameters until valid
- get random fragment from transcript
- get ends based on `length` parameter
- return data and write it collectively to output files

## Plotting
plots

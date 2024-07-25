## reads amplicon pairs file from circdna nextflow pipeline
## file : nextflowdirectory/results/ampliconsuite/ampliconclassifier/result/*_result_table.tsv
read_amp_class = function(file_path) {
    amp.dt = fread(file_path)
    names(amp.dt) = gsub(" ","_",names(amp.dt))
    amp.dt[, pair := Sample_name]
    return(amp.dt)
}

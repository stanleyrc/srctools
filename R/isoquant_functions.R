## title = isoquant_read_assignments
## read in read assignments - contains:
## assignment_events - events to classify the reads
## classification-novel_in_catalog,intergenic etc
## qname
## transcript_id (reference ID)
isoquant_read_assignments = function(read_assignments) {
  reads.dt = fread(read_assignments)
  reads.dt = reads.dt[2:nrow(reads.dt),]
  names(reads.dt) = c("qname", "chr", "strand", "transcript_id", "gene_id", "assignment_type", "assignment_events", "exons", "additional_info")
  reads.dt[, structural_category := tstrsplit(additional_info, "; Classification=", keep = 2)]
  reads.dt[!grepl("; Classification=",additional_info), structural_category := tstrsplit(additional_info, "Classification=", keep = 2)]
  reads.dt[, structural_category := gsub(";","", structural_category)]
  names(reads.dt)[names(reads.dt) == "strand"] = "strand_classification"
  names(reads.dt)[names(reads.dt) == "chr"] = "chr_classification"
  return(reads.dt)
}

## isoquant_transcript_counts
## returns feature_id (reference transcript) and count
## transcript_counts = pairs.dt[pair,transcript_counts]
isoquant_transcript_counts = function(transcript_counts) {
  counts.dt = fread(transcript_counts) %>% setnames(.,c("feature_id","count"))
  return(counts.dt)
}

## isoquant_transcript_tpm
## returns feature_id and tpm
## transcript_tpm = pairs.dt[pair,transcript_tpm]
isoquant_transcript_tpm = function(transcript_tpm) {
  tpm.dt = fread(transcript_tpm) %>% setnames(.,c("feature_id","tpm"))
  return(tpm.dt)
}

## isoquant_transcript_info
## uses isoquant_transcript_tpm and isoquant_transcript_counts to return feature_id count and tpm
isoquant_transcript_info = function(transcript_counts, transcript_tpm) {
  counts.dt = isoquant_transcript_counts(transcript_counts)
  tpm.dt = isoquant_transcript_tpm(transcript_tpm)
  tr.info.dt = merge.data.table(counts.dt,tpm.dt, by = "feature_id", all = TRUE)
  return(tr.info.dt)
}

## isoquant_bed
## outputs bed file with qname, start position of read and block sizes and starts
## corrected_reads_bed = pairs.dt[pair,corrected_reads_bed]
isoquant_bed = function(corrected_reads_bed) {
  bed.dt = fread(corrected_reads_bed)
  names(bed.dt) = gsub("#","",names(bed.dt))
  names(bed.dt)[1:3] = c("seqnames","start","end")
  return(bed.dt)
}

## isoquant_extended_annotation_gtf
## isoquant_extended_annotation_gtf = pairs.dt[pair,extended_annotation_gtf]
isoquant_extended_annotation_gtf = function(extended_annotation_gtf) {
  anno.dt = readGFF(extended_annotation_gtf) %>% as.data.table
  names(anno.dt)[1] = "seqnames"
  return(anno.dt)
}

## isoquant_transcript_models_gtf
## isoquant_transcript_models_gtf = pairs.dt[pair,transcript_models_gtf]
isoquant_transcript_models_gtf = function(transcript_models_gtf) {
  tr.dt = readGFF(transcript_models_gtf) %>% as.data.table
  names(tr.dt)[1] = "seqnames"
  return(tr.dt)
}

## isoquant_model_reads - transcript_id assignments for each read
## transcript_model_reads = pairs.dt[pair,transcript_model_reads]
isoquant_model_reads = function(transcript_model_reads) {
  tr.dt = fread(transcript_model_reads)
  names(tr.dt) = gsub("#","",names(tr.dt))
  return(tr.dt)
}

## isoquant_model_counts
## transcript_model_counts = pairs.dt[pair,transcript_model_counts]
isoquant_model_counts = function(transcript_model_counts) {
  counts.dt = fread(transcript_model_counts)
  names(counts.dt) = gsub("#","",names(counts.dt))
  return(counts.dt)
}

## isoquant_model_tpm
## transcript_model_tpm = pairs.dt[pair,transcript_model_tpm]
isoquant_model_tpm = function(transcript_model_tpm) {
  tpm.dt = fread(transcript_model_tpm)
  names(tpm.dt) = gsub("#","",names(tpm.dt))
  return(tpm.dt)
}

## isoquant_model_info
## transcript_model_counts = pairs.dt[pair,transcript_model_counts]
## transcript_model_tpm = pairs.dt[pair,transcript_model_tpm]
isoquant_model_info = function(transcript_model_counts, transcript_model_tpm) {
  counts.dt = isoquant_model_counts(transcript_model_counts)
  tpm.dt = isoquant_model_tpm(transcript_model_tpm)
  tr.info.dt = merge.data.table(counts.dt,tpm.dt, by = "feature_id", all = TRUE)
  return(tr.info.dt)
}

## isoquant_gene_counts
## gene_counts = pairs.dt[pair,gene_counts]
isoquant_gene_counts = function(gene_counts) {
  counts.dt = fread(gene_counts)
  names(counts.dt) = gsub("#","",names(counts.dt))
  return(counts.dt)
}

## isoquant_gene_tpm
## gene_tpm = pairs.dt[pair,gene_tpm]
isoquant_gene_tpm = function(gene_tpm) {
  tpm.dt = fread(gene_tpm)
  names(tpm.dt) = gsub("#","",names(tpm.dt))
  return(tpm.dt)
}

## isoquant_gene_info
##gene_counts = pairs.dt[pair,gene_counts]
## gene_tpm = pairs.dt[pair,gene_tpm]
isoquant_gene_info = function(gene_counts, gene_tpm) {
  counts.dt = isoquant_gene_counts(gene_counts)
  tpm.dt = isoquant_gene_tpm(gene_tpm)
  tr.info.dt = merge.data.table(counts.dt,tpm.dt, by = "feature_id", all = TRUE)
  return(tr.info.dt)
}

## gtf2dt for getting the annotation file used
## gtf = "/gpfs/data/imielinskilab/projects/melanoma_pdx/data/melanoma_cell_lines_isoseq_pbmm2/pacbio_ref_data/gencode.v39.annotation.sorted.gtf"
gtf2dt = function(gtf) {
  gtf.dt = readGFF(gtf) %>% as.data.table
  names(gtf.dt)[1] = "seqnames"
  return(gtf.dt)
}


## function to get the reads the map to a specific gene
## column first_read_assignment represents unique isoform reads to use when subsetting with get_iso_reads
## function to get the reads the mapp to a specific gene
## function to get the reads the mapp to a specific gene
reads_genes = function(ref_gtf.gr, sample_gtf.gr, reads.dt, genes, group_by = "exons", pad = 0) {
  if(inherits(ref_gtf.gr, "character")) {
    if(file.exists(ref_gtf.gr)) {
      ref_gtf.gr = gtf2dt(ref_gtf.gr) %>% GRanges(.,seqlengths = hg_seqlengths())
    }
  } else if (inherits(ref_gtf.gr, "GRanges")) {
    message("User supplied ref_gtf.gr as GRanges. Continuing...")
  } else if(inherits(ref_gtf.gr, "data.table")) {
    warning("User supplied ref_gtf.gr as data.table. Converting to a GRanges")
    ref_gtf.gr = GRanges(ref_gtf.gr, seqlengths = hg_seqlengths())
  }
  if(inherits(sample_gtf.gr, "character")) {
    if(file.exists(sample_gtf.gr)) {
      sample_gtf.gr = isoquant_extended_annotation_gtf(extended_annotation_gtf = sample_gtf.gr) %>% GRanges(.,seqlengths = hg_seqlengths())
    }
  } else if (inherits(sample_gtf.gr, "GRanges")) {
    message("User supplied sample_gtf.gr as GRanges. Continuing...")
  } else if(inherits(sample_gtf.gr, "data.table")) {
    message("User supplied sample_gtf.gr as data.table. Converting to a GRanges")
    sample_gtf.gr = GRanges(sample_gtf.gr, seqlengths = hg_seqlengths())
  } 

  if(inherits(reads.dt, "character")) {
    reads.dt = isoquant_read_assignments(reads.dt)
  } else if (inherits(reads.dt,"data.table")) {
    message("User supplied reads.dt as data.table")
  }
  sample_gtf_sub.gr = sample_gtf.gr %&% ((ref_gtf.gr %Q% (gene_name %in% genes)) + pad)
  ##get unique transcript ids contributing to this gene region
  uniq.gene.tr_ids.lst = sample_gtf_sub.gr$transcript_id %>% unique %>% na.omit %>% as.character
  ## need to add model_reads.dt to get all of the qnames here
  ## uniq.gene.tr_ids.lst
  ## get accurate counts and tpm by counting the number of reads with a specific transcript_id and assignment_events
  reads.dt[, count_assignment := .N, by = c("transcript_id","assignment_events")]
  reads.dt[, count_assignment_type := .N, by = c("transcript_id","assignment_events","assignment_type")]
  reads.dt[, count_assignment_type_exons := .N, by = c("transcript_id","assignment_events","assignment_type","exons")]
  gene.dt = reads.dt[transcript_id %in% uniq.gene.tr_ids.lst,]
  gene.dt = gene.dt[order(-count_assignment),]
  reads.dt[, grp_assignment := .GRP, by = c("transcript_id","assignment_events")]
  gene.dt = gene.dt[order(-count_assignment_type),]
  reads.dt[, grp_assignment_type := .GRP, by = c("transcript_id","assignment_events","assignment_type")]
  gene.dt = gene.dt[order(-count_assignment_type_exons),]
  reads.dt[, grp_assignment_type_exons := .GRP, by = c("transcript_id","assignment_events","assignment_type","exons")]

  ## 
  gene.dt = reads.dt[transcript_id %in% uniq.gene.tr_ids.lst,]
  gene.dt = gene.dt[order(-count_assignment_type_exons),]## [count_assignment_type_exons != 1,]
  gene.dt[, exons_list := list(strsplit(exons, ",", fixed = TRUE)), by = "qname"]
  gene.dt[, exons_minus_ends_list := lapply(exons_list, function(x) {
    return(x[2:(length(x)-1)])
  })]
  gene.dt[, exons_minus_ends := sapply(exons_minus_ends_list, function(x) {
    return(paste0(x, collapse = ","))
  })]
  gene.dt[, count_assignment_type_exons_minus_ends := .N, by = c("transcript_id","assignment_events","assignment_type","exons_minus_ends")]
  gene.dt = gene.dt[order(-count_assignment_type_exons_minus_ends),]
  gene.dt[, grp_assignment_type_exons_minus_ends := .GRP, by = c("transcript_id","assignment_events","assignment_type","exons_minus_ends")]
  gene.dt[, count_exons_only := .N, by = c("transcript_id","exons_minus_ends")]
  gene.dt = gene.dt[order(-count_exons_only),]
  gene.dt[, grp_exons_only := .GRP, by = c("transcript_id","exons_minus_ends")]
  if(group_by == "exons_assignment_events") {
  ## get all reads as a list object 
    gene.dt[, reads_grp_assignment_type_exons_minus_ends := list(list(sort(qname))), by = grp_assignment_type_exons_minus_ends]
    ## get first read for each so I can plot every GENE variant
    gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type_exons_minus_ends)[1], by = grp_assignment_type_exons]
    gene.dt[, final_grp_assignment := grp_assignment_type_exons_minus_ends]
    gene.dt[, final_count_assignment := count_assignment_type_exons_minus_ends]
  } else if (group_by == "assignment_events") {
    gene.dt[, reads_grp_assignment_type := list(list(sort(qname))), by = grp_assignment_type]
    gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type)[1], by = grp_assignment_type]
    gene.dt[, final_grp_assignment := grp_assignment_type]
    gene.dt[, final_count_assignment := count_assignment_type]
  } else if (group_by == "exons_only") {
    gene.dt[, reads_grp_assignment_type := list(list(sort(qname))), by = grp_exons_only]
    gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type)[1], by = grp_exons_only]
    gene.dt[, final_grp_assignment := grp_exons_only]
    gene.dt[, final_count_assignment := count_exons_only]
  }
  return(gene.dt)
}

## function to use reads_genes and get the gwalks for that
unique_reads_gw = function(bam, #aligned bam to get reads
                           ref_gtf.gr, # reference gtf-can be character to path, granges, or data.table that can be converted
                           sample_gtf.gr, # sample gtf output of isoquant-can be character to path, granges, or data.table that can be converted
                           reads.dt,      #reads.dt- output of isoquant_read_assignments after reading in read assignments
                           genes,         #what genes to get unique isoforms
                           cores = 2,     #cores for get_iso_reads and adding y
                           pad = 0,       #pad for genes around a specific region
                           group_by = "exons_assignment_events", #group by for grouping the reads-exons_assignment_events is the most strict-all exons except the first and last have to be correct. assignment_type is slightly less strict and uses the assignment_events from isoquant
                           min_reads_isoform = NULL,
                           return_type = "gw", #return as a gw or grl for the first element returned
                           get_reads = TRUE, #boolean of whether to just return the read assignments or reads assignments and reads as gw or grl
                           add_y = TRUE, # if true, forces return_type to gw but 
                           pad_get_iso = 10000) { #pad for reading in the specific qnames
  if(inherits(reads.dt, "character")) {
    reads.dt = isoquant_read_assignments(reads.dt)
  } else if (inherits(reads.dt,"data.table")) {
    message("User supplied reads.dt as data.table")
  }
  gene_reads.dt = reads_genes(ref_gtf.gr = ref_gtf.gr, sample_gtf.gr = sample_gtf.gr, reads.dt = reads.dt, genes = genes, pad = pad, group_by = group_by)
  if(!is.null(min_reads_isoform)) {
    ## if(group_by == "exons_assignment_events") {
    ##   gene_reads.dt = gene_reads.dt[count_assignment_type_exons_minus_ends > min_reads_isoform,]
    ## } else if (group_by == "assignment_type") {
    ##   gene_reads.dt = gene_reads.dt[count_assignment_type > min_reads_isoform,]
    ## }
    ## } else if (group_by == "exons_only") {
    ##   gene_reads.dt = gene_reads.dt[count_assignment_type > min_reads_isoform,]
    ## }
    gene_reads.dt = gene_reads.dt[final_count_assignment >= min_reads_isoform,]
  }
  qname_list = gene_reads.dt$first_read_assignment %>% unique
  gene.gr = (ref_gtf.gr %Q% (gene_name %in% genes)) %>% gr.reduce
  if(get_reads) {
    gene.gw = get_iso_reads(bam,
                            gr = (gene.gr + pad_get_iso),
                            gtf = ref_gtf.gr,
                            cores = cores,
                            reannotate = FALSE,
                            annotate_missing = TRUE,
                            read_assignments = reads.dt,
                            type = return_type,
                            pipeline = "isoquant",
                            subset_qname = qname_list)
    if(return_type == "gw") {
      ## add the transcript name and id that is plotted to gene_reads.dt
      gene.gw.dt = as.data.table(grl.unlist(gene.gw$grl))[,.(transcript_name,transcript_id,qname)] %>% setnames(., c("plot_transcript_name", "plot_transcript_id", "first_read_assignment")) %>% unique
    } else if(return_type =="grl") {
      ## add the transcript name and id that is plotted to gene_reads.dt
      gene.gw.dt = as.data.table(grl.unlist(gene.gw))[,.(transcript_name,transcript_id,qname)] %>% setnames(., c("plot_transcript_name", "plot_transcript_id", "first_read_assignment")) %>% unique
    }
    gene_reads.dt = merge.data.table(gene.gw.dt,gene_reads.dt, by = "first_read_assignment", all = TRUE)
    if(add_y) {
      if(inherits(gene.gw,"gWalk")) {
        gene.gw = gene.gw$grl
      }
      gene.gw.lst = mclapply(1:length(gene.gw), function(x) {
        gene.sub.gw = gene.gw[[x]]
        ## qname1 = names(gene.gw[x])
        ## gene_reads_sub.dt = gene_reads.dt[first_read_assignment == names(gene.gw[x]),]
        qname1 = gene.sub.gw$qname[1] %>% unique
        gene_reads_sub.dt = gene_reads.dt[first_read_assignment == qname1,]
        add_count = unique(gene_reads_sub.dt$final_count_assignment)
        ## sometimes the same read can be for multiple "isoforms" - if one read was assigned to multiple isoforms from isoqaunt
        ## I have only seen reads with one count like this
        if(length(add_count) > 1) {
          gene.sub.gw$count = paste0(add_count, collapse = ",")
        } else {
          gene.sub.gw$count = as.character(add_count)
        }
        add_grp = unique(gene_reads_sub.dt$final_grp_assignment)
        if(length(add_grp) > 1) {
          gene.sub.gw$grp = paste0(add_grp,collapse = ",")
        } else {
          gene.sub.gw$grp = as.character(add_grp)
        }
        ## add the grp assignment to the walks grl
        tr = gene.sub.gw$transcript_id[1] %>% unique
        grp_assignment = as.character(gene_reads_sub.dt[transcript_id == tr & first_read_assignment == qname1,]$final_grp_assignment %>% unique) %>% paste0("GRP.",.)
        qname_list = unlist(gene_reads_sub.dt[transcript_id == tr & first_read_assignment == qname1,]$reads_grp_assignment_type) %>% unique
        if(grp_assignment == "GRP." & is.null(qname_list)) {
          grp_assignment = as.character(gene_reads_sub.dt[first_read_assignment == qname1,]$final_grp_assignment %>% unique) %>% paste0("GRP.",.)
          qname_list = unlist(gene_reads_sub.dt[first_read_assignment == qname1,]$reads_grp_assignment_type) %>% unique
        }
        gene.sub.gw$grp_assignment = grp_assignment
        ## gene.sub.gw$qname_list = list(qname_list)
        return(as.data.table(gene.sub.gw))
      },mc.cores = cores)
      gene.gw.dt = rbindlist(gene.gw.lst, fill = TRUE)
      ## gene.gw.dt[, y := as.numeric(count)]
      gene.gw.gr = GRanges(gene.gw.dt, seqlengths = hg_seqlengths())
      gene.gw = rtracklayer::split(gene.gw.gr, f = mcols(gene.gw.gr)["qname"])
      if(return_type == "gw") {
        message("Converting to gW")
        gene.gw = gW(grl = gene.gw)
        ## add transcript names
        message("Adding transcript names as: transcript_name ; Count")
        name.gw(gene.gw,c("transcript_name","count"))
        ## transcript_names = sapply(gene.gw$grl, function(x) x$transcript_name[1])
        ## gene.gw$set(name = transcript_names)
        ## add y values in order of abundance
        message("Adding y values in order of abundance")
        abundance_order.dt = data.table(count = (sapply(gene.gw$grl, function(x) x$count[1]) %>% as.numeric))
        abundance_order.dt[,current_order := 1:.N]
        abundance_order.dt = abundance_order.dt[order(-count),]
        abundance_order.dt[, y_value := .N:1]
        abundance_order.dt = abundance_order.dt[order(current_order),]
        gene.gw$set(y = abundance_order.dt$y_value)
      }
    }
    return(list(gene.gw, gene_reads.dt))
  } else {
    return(gene_reads.dt)
  }
}


## function to group transcripts across samples
## use output of unique_reads_gw with get_reads FALSE for many samples as input this function
genes_reads_group_by_samples = function(genes_reads.dt, ref_gtf.dt, group_by = "exons_assignment_events") {
  genes_reads.dt[, count_assignment_all_samples := .N, by = c("transcript_id","assignment_events")]
  genes_reads.dt[, count_assignment_type_all_samples := .N, by = c("transcript_id","assignment_events","assignment_type")]
  genes_reads.dt[, count_assignment_type_exons_all_samples := .N, by = c("transcript_id","assignment_events","assignment_type","exons")]
  genes_reads.dt = genes_reads.dt[order(-count_assignment_type_exons),]## [count_assignment_type_exons != 1,]
  genes_reads.dt[, exons_list_all_samples := list(strsplit(exons, ",", fixed = TRUE)), by = "qname"]
  genes_reads.dt[, exons_minus_ends_list_all_samples := lapply(exons_list, function(x) {
    return(x[2:(length(x)-1)])
  })]
  genes_reads.dt[, exons_minus_ends_all_samples := sapply(exons_minus_ends_list, function(x) {
    return(paste0(x, collapse = ","))
  })]
  genes_reads.dt[, count_assignment_type_exons_minus_ends_all_samples := .N, by = c("transcript_id","assignment_events","assignment_type","exons_minus_ends")]
  genes_reads.dt = genes_reads.dt[order(-count_assignment_type_exons_minus_ends_all_samples),]
  genes_reads.dt[, grp_assignment_type_exons_minus_ends_all_samples := .GRP, by = c("transcript_id","assignment_events","assignment_type","exons_minus_ends")]
  genes_reads.dt[, count_exons_only_all_samples := .N, by = c("transcript_id","exons_minus_ends")]
  genes_reads.dt = genes_reads.dt[order(-count_exons_only_all_samples),]
  genes_reads.dt[, grp_exons_only_all_samples := .GRP, by = c("transcript_id","exons_minus_ends")]
  ## use group_by = "exons_assignment_events"
  ## genes_reads.dt[, reads_grp_assignment_type_exons_minus_ends := list(list(sort(qname))), by = grp_assignment_type_exons_minus_ends_all_samples]
  genes_reads.dt = genes_reads.dt[order(-count_assignment_type),]
  genes_reads.dt[, grp_assignment_type_all_samples := .GRP, by = c("transcript_id","assignment_events","assignment_type")]
  genes_reads.dt = genes_reads.dt[order(-count_assignment_type_exons_all_samples),]
  genes_reads.dt[, grp_assignment_type_exons := .GRP, by = c("transcript_id","assignment_events","assignment_type","exons")]

  ## get first read for each so I can plot every GENE variant
  if(group_by == "exons_assignment_events") {
    genes_reads.dt[, first_read_assignment := unlist(reads_grp_assignment_type_exons_minus_ends)[1], by = grp_assignment_type_exons_minus_ends_all_samples]
    genes_reads.dt[, final_grp_assignment := grp_assignment_type_exons_minus_ends_all_samples]
    genes_reads.dt[, final_count_assignment := count_assignment_type_exons_minus_ends_all_samples]
  } else if (group_by == "assignment_events") {
    genes_reads.dt[, first_read_assignment := unlist(reads_grp_assignment_type)[1], by = grp_assignment_type_all_samples]
    genes_reads.dt[, final_grp_assignment := grp_assignment_type_all_samples]
    genes_reads.dt[, final_count_assignment := count_assignment_type_all_samples]
  } else if (group_by == "exons_only") {
    stop("exons_only not implemented yet. Too strict")
  }
## ##########################
##   if(group_by == "exons_assignment_events") {
##   ## get all reads as a list object 
##     gene.dt[, reads_grp_assignment_type_exons_minus_ends := list(list(sort(qname))), by = grp_assignment_type_exons_minus_ends_all_samples]
##     ## get first read for each so I can plot every GENE variant
##     gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type_exons_minus_ends)[1], by = grp_assignment_type_exons_minus_ends_all_samples]
##     gene.dt[, final_grp_assignment := grp_assignment_type_exons_minus_ends]
##     gene.dt[, final_count_assignment := count_assignment_type_exons_minus_ends]
##   } else if (group_by == "assignment_events") {
##     gene.dt[, reads_grp_assignment_type := list(list(sort(qname))), by = grp_assignment_type]
##     gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type)[1], by = grp_assignment_type]
##     gene.dt[, final_grp_assignment := grp_assignment_type]
##     gene.dt[, final_count_assignment := count_assignment_type]
##   } else if (group_by == "exons_only") {
##     gene.dt[, reads_grp_assignment_type := list(list(sort(qname))), by = grp_exons_only]
##     gene.dt[, first_read_assignment := unlist(reads_grp_assignment_type)[1], by = grp_exons_only]
##     gene.dt[, final_grp_assignment := grp_exons_only]
##     gene.dt[, final_count_assignment := count_exons_only]
##   }

  ########################
  ## now add transcript name to the table as well
  gtf_sub.dt = ref_gtf.dt[,.(transcript_id,gene_name,transcript_name)] %>% unique
  gtf_sub.dt = gtf_sub.dt[!is.na(transcript_id),]
  genes_reads.dt = merge.data.table(genes_reads.dt, gtf_sub.dt, by = "transcript_id", all.x = TRUE)
  ## add what sample is for each read
  sample_qname.dt = unique(genes_reads.dt[,.(qname,sample)]) %>% setnames(.,c("qname","first_read_assignment_sample"))
  genes_reads.dt = merge.data.table(genes_reads.dt, sample_qname.dt, by.x = "first_read_assignment", by.y = "qname", all.x = TRUE)
  ## add counts per sample as separate columns for each patient
  if(group_by == "assignment_events_exons") {
    genes_reads.dt[, final_count_assignment := count_assignment_type_exons_minus_ends]
    genes_reads.dt[, final_count_assignment_all_samples := count_assignment_type_exons_minus_ends_all_samples]
    counts.dt = genes_reads.dt[, .N, by = c("transcript_id", "assignment_events", "assignment_type", "exons_minus_ends", "sample")]
    setnames(counts.dt, "N", "count_assignment_type_exons_minus_ends_all_samples")
    counts.dt[,sample := paste0("counts_",sample)]
    counts_wide.dt = dcast(counts.dt, transcript_id + assignment_events + assignment_type + exons_minus_ends ~ sample, value.var = "count_assignment_type_exons_minus_ends_all_samples", fill = 0)
    ## add columns with counts for each
    merged.dt = merge.data.table(genes_reads.dt, counts_wide.dt, by = c("transcript_id", "assignment_events", "assignment_type", "exons_minus_ends"), all.x = TRUE)
    merged.dt[, total_counts_by_gene_name := .N, by = "gene_name"]
  } else if (group_by == "assignment_events") {
    genes_reads.dt[, final_count_assignment := count_assignment_type]
    genes_reads.dt[, final_count_assignment_all_samples := count_assignment_type_all_samples]
    counts.dt = genes_reads.dt[, .N, by = c("transcript_id", "assignment_events", "assignment_type", "sample")]
    setnames(counts.dt, "N", "count_assignment_type_all_samples")
    counts.dt[,sample := paste0("counts_",sample)]
    counts_wide.dt = dcast.data.table(counts.dt, transcript_id + assignment_events + assignment_type ~ sample, value.var = "count_assignment_type_all_samples", fill = 0)
    ## add columns with counts for each
    merged.dt = merge.data.table(genes_reads.dt, counts_wide.dt, by = c("transcript_id", "assignment_events", "assignment_type"), all.x = TRUE)
    merged.dt[, total_counts_by_gene_name := .N, by = "gene_name"]
  }
  return(merged.dt)
}


## function to take genes_reads.dt and subset based on counts and return a gwalk with the isoforms
genes_reads2gw = function(genes_reads.dt, #output of genes_reads_group_by_samples
                          pairs.dt,       #pairs table to access bam
                          min_counts_isoform = NULL, #minimum counts for a given isoform assigned
                          min_perc_gene_reads_isoform = NULL, #minimum percent of reads for a gene is assigned to the isoform
                          get_reads = TRUE,                   #TRUE means return the gwalk as well, FALSE will just filter
                          cores = 4,                          #number of cores for getting the reads
                          ref_gtf.gr,                         #reference gtf as a granges
                          bam_column = "aligned_bam",         #column in pairs table for bam file
                          read_assignments_column = "read_assignments", #column in pairs table with read assignments
                          group_by = "assignment_events_exons",
                          return_type = "gw",                          #whether to return as a gw, else will return grl
                          pad_get_iso = 10000,                            #value to pad reading in the bam files
                          ignore_polyA_assignment = TRUE)                 #whether to ignore polyA=False and polyA=TRUE in additional_info
  {
  if(!is.null(min_counts_isoform)) {
    genes_reads.dt2 = genes_reads.dt[final_count_assignment_all_samples > min_counts_isoform,]
  }
  if(!is.null(min_perc_gene_reads_isoform)) {
    genes_reads.dt2 = genes_reads.dt[final_count_assignment_all_samples > (total_counts_by_gene_name *min_perc_gene_reads_isoform),]
  }
  if((!is.null(min_perc_gene_reads_isoform)) & (!is.null(min_counts_isoform))) {
    genes_reads.dt2 = genes_reads.dt
  }
  if(nrow(genes_reads.dt2) == 0) {
    return(NULL)
  }
  if(get_reads) {
    unique_samples_reads.dt = genes_reads.dt2[,.(first_read_assignment,gene_name,first_read_assignment_sample,final_grp_assignment,final_count_assignment_all_samples)] %>% unique
    uniq_samples = unique(unique_samples_reads.dt$first_read_assignment_sample)
    message("Getting one read per unique isoform")
    genes.gw.lst = mclapply(uniq_samples, function(sample) {
      unique.sub.dt = unique_samples_reads.dt[first_read_assignment_sample == sample,]
      sub.pairs.dt = pairs.dt[sample,]
      bam = sub.pairs.dt[[bam_column]]
      ## reads_file = sub.pairs.dt[[read_assignments_column]]
      ## reads.dt = isoquant_read_assignments(reads_file)
      qname_list = unique.sub.dt$first_read_assignment
      gene.gr = (ref_gtf.gr %Q% (gene_name %in% unique.sub.dt$gene_name & type == "transcript")) %>% gr.reduce
      gene.gr = (gene.gr + pad_get_iso) %>% gr.reduce
      cols_keep = c("qname", "chr_classification", "strand_classification", "transcript_id", "gene_id", "assignment_type", "assignment_events", "exons", "additional_info", "structural_category")
      reads_sub.dt = genes_reads.dt[,..cols_keep]
      gene.gw = get_iso_reads(bam,
                              gr = (gene.gr),
                              gtf = ref_gtf.gr,
                              cores = 1,
                              reannotate = FALSE,
                              annotate_missing = TRUE,
                              ## read_assignments = reads.dt,
                              read_assignments = reads_sub.dt,
                              type = "gw",
                              pipeline = "isoquant",
                              subset_qname = qname_list)
      gene.gw = gene.gw$grl
    ## if(return_type == "gw") {
      ## add the transcript name and id that is plotted to gene_reads.dt
      ## gene.gw.dt = as.data.table(grl.unlist(gene.gw$grl))[,.(transcript_name,transcript_id,qname)] %>% setnames(., c("plot_transcript_name", "plot_transcript_id", "first_read_assignment")) %>% unique
    ## } else if(return_type =="grl") {
    ##   ## add the transcript name and id that is plotted to gene_reads.dt
      gene.gw.dt = as.data.table(grl.unlist(gene.gw))[,.(transcript_name,transcript_id,qname)] %>% setnames(., c("plot_transcript_name", "plot_transcript_id", "first_read_assignment")) %>% unique
      ## }
      gene.gw.dt = merge.data.table(gene.gw.dt, unique.sub.dt[,.(first_read_assignment,final_grp_assignment)], by = "first_read_assignment", all = TRUE)
      genes_reads.dt3 = merge.data.table(gene.gw.dt,genes_reads.dt2, by = c("first_read_assignment", "final_grp_assignment"),all.x = TRUE)
      gene.gw.lst = mclapply(1:length(gene.gw), function(x) {
        gene.sub.gw = gene.gw[[x]]
        ## qname1 = names(gene.gw[x])
        ## genes_reads_sub.dt = genes_reads.dt3[first_read_assignment == names(gene.gw[x]),]
        qname1 = gene.sub.gw$qname[1] %>% unique
        genes_reads_sub.dt = genes_reads.dt3[first_read_assignment == qname1,]
        add_count = unique(genes_reads_sub.dt$final_count_assignment_all_samples)
        ## sometimes the same read can be for multiple "isoforms" - if one read was assigned to multiple isoforms from isoqaunt
        ## I have only seen reads with one count like this
        if(length(add_count) > 1) {
          gene.sub.gw$count = paste0(add_count, collapse = ",")
        } else {
          gene.sub.gw$count = as.character(add_count)
        }
        add_grp = unique(genes_reads_sub.dt$final_grp_assignment)
        if(length(add_grp) > 1) {
          gene.sub.gw$grp = paste0(add_grp,collapse = ",")
        } else {
          gene.sub.gw$grp = as.character(add_grp)
        }
        ## add the grp assignment to the walks grl
        tr = gene.sub.gw$transcript_id[1] %>% unique
        grp_assignment = as.character(genes_reads_sub.dt[transcript_id == tr & first_read_assignment == qname1,]$final_grp_assignment %>% unique) %>% paste0("GRP.",.)
        qname_list = unlist(genes_reads_sub.dt[transcript_id == tr & first_read_assignment == qname1,]$reads_grp_assignment_type) %>% unique
        if(grp_assignment == "GRP." & is.null(qname_list)) {
          grp_assignment = as.character(genes_reads_sub.dt[first_read_assignment == qname1,]$final_grp_assignment %>% unique) %>% paste0("GRP.",.)
          qname_list = unlist(genes_reads_sub.dt[first_read_assignment == qname1,]$reads_grp_assignment_type) %>% unique
        }
        if(length(add_grp) > 1) {
          gene.sub.gw$grp_assignment = paste0(grp_assignment, collapse = ",")
        } else {
          gene.sub.gw$grp_assignment = grp_assignment
        }
        ## gene.sub.gw$qname_list = list(qname_list)
        return(as.data.table(gene.sub.gw))
      },mc.cores = 1)
      gene.gw.dt = rbindlist(gene.gw.lst, fill = TRUE)
      return(gene.gw.dt)
    }, mc.cores = cores)
    
    message("Binding all unique walks per read")
    genes.gw.dt = rbindlist(genes.gw.lst,fill = TRUE)
    message("Converting reads to GRanges")    
    gene.gw.gr = GRanges(genes.gw.dt, seqlengths = hg_seqlengths())
    message("Converting GRanges to GRangesList")
    gene.gw = rtracklayer::split(gene.gw.gr, f = mcols(gene.gw.gr)["qname"])
    ## browser()
    if(return_type == "gw") {
      message("Converting to gW")
      gene.gw = gW(grl = gene.gw)
      ## add transcript names
      message("Adding transcript names as: transcript_name ; Count")
      if(length(gene.gw) > 1) {
        name.gw(gene.gw,c("transcript_name","count"))
      } else if (length(gene.gw) == 1){
        new_name = paste(gene.gw$grl[[1]]$transcript_name[1], gene.gw$grl[[1]]$count[1], sep = "; ")
        gene.gw$set(name = new_name)
      }
      ## transcript_names = sapply(gene.gw$grl, function(x) x$transcript_name[1])
      ## gene.gw$set(name = transcript_names)
      ## add y values in order of abundance
      ## browser()
      message("Adding y values in order of abundance. Adding count for all samples and the grp_assignment.")
      abundance_order.dt = data.table(count = (sapply(gene.gw$grl, function(x) x$count[[1]])), grp_assignment = (sapply(gene.gw$grl, function(x) x$grp_assignment[[1]])),transcript_name = (sapply(gene.gw$grl, function(x) x$transcript_name[[1]])), transcript_id = (sapply(gene.gw$grl, function(x) x$transcript_id[[1]])))
      abundance_order.dt[,current_order := 1:.N]
      abundance_order.dt[, max_count := sapply(count, function(x) {
        max(as.numeric(tstrsplit(x,",")))
      })]
      abundance_order.dt = abundance_order.dt[order(-max_count),]
      abundance_order.dt[, y_value := .N:1]
      abundance_order.dt = abundance_order.dt[order(current_order),]
      gene.gw$set(y = abundance_order.dt$y_value)
      gene.gw$set(count = abundance_order.dt$max_count)
      gene.gw$set(grp_assignment = abundance_order.dt$grp_assignment)
      gene.gw$set(transcript_name = abundance_order.dt$transcript_name)
      gene.gw$set(transcript_id = abundance_order.dt$transcript_id)
      message("Adding assignment_events, assignment_type, exons_minus_ends, structural_category, final_grp_assignment, additional_info, and counts for each sample.")
      gene.gw.dt = gene.gw$dt
      if(group_by == "assignment_events_exons") {
        cols_add = c("assignment_events", "assignment_type", "exons_minus_ends", "structural_category", "final_grp_assignment", "additional_info", "gene_name")
        cols_add = c(cols_add,grep("^counts_",names(genes_reads.dt2),value = TRUE))
      } else if(group_by == "assignment_events") {
        cols_add = c("assignment_events", "assignment_type","structural_category", "final_grp_assignment", "additional_info", "gene_name")
        cols_add = c(cols_add,grep("^counts_",names(genes_reads.dt2),value = TRUE))
      }
      genes_reads.dt2_merge = genes_reads.dt2[,..cols_add] %>% unique
      if(ignore_polyA_assignment) {
        genes_reads.dt2_merge[, additional_info := gsub("PolyA=False|PolyA=True", "PolyA=ignored",additional_info)]
        genes_reads.dt2_merge = unique(genes_reads.dt2_merge)
      }
      genes_reads.dt2_merge[,final_grp_assignment := paste0("GRP.",final_grp_assignment)]
      genes_reads.dt2_merge = genes_reads.dt2_merge[final_grp_assignment %in% gene.gw.dt$grp_assignment,]
      
      if(nrow(genes_reads.dt2_merge) != length(gene.gw)) {
        warning("Genes_reads.dt2 unique values are not the same length as the number of walks. Some transcripts have the same read assigned so these will not have counts") #
      }
      set.dt = merge.data.table(gene.gw.dt[,.(walk.id,grp_assignment)], genes_reads.dt2_merge, by.x = "grp_assignment", by.y = "final_grp_assignment", all = TRUE)
      set.dt = set.dt[order(walk.id),]
      cols_add2 = names(set.dt)
      cols_add2 = cols_add2[!(cols_add2 %in% c("walk.id","grp_assignment"))]
      for(x in cols_add2) {
        expr <- parse(text = paste0("gene.gw$set(", x, " = set.dt[['", x, "']])")) #set for all to add to the data.table
        eval(expr)
      }
      ## } else {
      ## warning("Genes_reads.dt2 unique values are not the same length as the number of walks. Extra annotations will not be added")
      ## }
      gene.gw = gene.gw[order(-count)]
    }
    return(list(gene.gw, genes_reads.dt2))
  } else {
    return(genes_reads.dt2)
  }
}


## function for plotting the gwalks that come out of genes_reads2gw
plot_isoforms_gw = function(gw,
                         gr,
                         samples,
                         gencode = NULL,
                         gw_plotname = "",
                         plot1 = "plot.png",
                         plot2 = "plot2.png",
                         plot3 = "plot3.png",
                         gw_height = 30,
                         gencode_height = 15,
                         plot_pad = 0,
                         subset_gencode = TRUE,
                         name = c("transcript_name","grp_assignment"),
                         sep_name = "; ",
                         counts_col_prefix = "counts_",
                         plot_title = "",
                         counts_theme = NULL, #theme for the counts plot
                         total_reads = NULL,  #if not null will plot a third plot with tpm-has to be the same length as samples
                         max_isoforms = NULL,
                         max_order_by = "Count", #relevant if max_isoforms is not NULL, takes the top N isoforms by either "Count" for ordering isoforms by single sample counts or "total_count" for ordering by the total counts across samples
                         gene = NULL, #if NULL will consider all genes within the specified granges, TRUE will filter for a specfied gene name
                         walk_height = 20000,
                         walk_width = 12000,
                         walk_res = 300,
                         counts_height = 5000,
                         counts_width = 8000,
                         counts_res = 300) {
  gw2 = gw %&% gr
  if(!is.null(gene)) {
    indices = gw2$dt[gene_name == gene]$walk.id
    gw2 = gw2[indices]
  }
  plot.dt = gw2$dt
  ## resetting the y based on the subset transcripts
  gw2$set(y = length(gw2):1)
  if(subset_gencode & !is.null(gencode)) {
    message("subsetting gencode to the transcripts in the gwalk")
    gencode.subset = subset_gencode(gencode = gencode, gw = gw2, height = gencode_height)
    gencode.subset$height = gencode_height
  } else if(!is.null(gencode)) {
    gencode.subset = gencode
    gencode.subset$height = gencode_height
  } else {
    gencode.subset = NULL
  }
  if(plot_pad != 0) {
    gr = gr + plot_pad
  }
  ## set the names specified
  name.dt = plot.dt[,..name]
  name.dt[, new_name := apply(.SD, 1, function(row) paste(row, collapse = sep_name))]
  gw2$set(name = name.dt$new_name)
  if(is.null(max_isoforms)) {
    ppng(plot(c(gw2$gtrack(name = gw_plotname,height = gw_height,y.field = "y",yaxis = FALSE), gencode.subset), gr), height = walk_height, width = walk_width, res = walk_res, file = plot1)
  }
  
  ## now make the barplot
  sample_cols = paste0(counts_col_prefix, samples)
  cols_keep = c("grp_assignment","transcript_name","gene_name",sample_cols)
  plot.dt2 = plot.dt[,..cols_keep]
  long.dt = melt.data.table(data = plot.dt2,id.vars = c("grp_assignment","transcript_name","gene_name"), value.name = "Count", variable.name = "Sample")
  long.dt[, total_count := sum(Count), by = grp_assignment]
  long.dt = long.dt[order(-total_count),]
  ## add the grp_assignment as a factor so the order of the plot is the same as the order of the gwalks
  long.dt[, grp_assignment := factor(grp_assignment, levels = unique(long.dt$grp_assignment))]
  sub.dt = long.dt[,..name] #subset data.table to add same name as walks
  sub.dt[, new_name := apply(.SD, 1, function(row) paste(row, collapse = sep_name))]
  long.dt[,new_name := sub.dt$new_name]
  long.dt[,Sample := gsub(counts_col_prefix,"",Sample)]
  long.dt = long.dt[order(-total_count),]
  long.dt[, new_name := factor(new_name, levels = unique(long.dt$new_name))]
  unique_grp_length = long.dt$grp_assignment %>% unique %>% length
  if(!is.null(max_isoforms) & (unique_grp_length > max_isoforms)) {
    counts.dt = long.dt[,.(grp_assignment,Count,total_count)]
    if(max_order_by == "Count") {
      counts.dt = counts.dt[order(-Count),]
    } else if (max_order_by == "total_count") {
      counts.dt = counts.dt[order(-total_count),]
    }
    grp_subset = counts.dt$grp_assignment %>% unique
    grp_subset = grp_subset[1:max_isoforms] %>% na.omit %>% as.character
    long.dt = long.dt[grp_assignment %in% grp_subset,]
  }
  if(!is.null(max_isoforms) & (unique_grp_length > max_isoforms)) {
    gw3 = gw2[grp_assignment %in% grp_subset]
    gw3$set(y = length(gw3):1)
    ppng(plot(c(gw3$gtrack(name = gw_plotname,height = gw_height,y.field = "y",yaxis = FALSE), gencode.subset), gr), height = walk_height, width = walk_width, res = walk_res, file = plot1)
  } else if(!is.null(max_isoforms)) {
    ppng(plot(c(gw2$gtrack(name = gw_plotname,height = gw_height,y.field = "y",yaxis = FALSE), gencode.subset), gr), height = walk_height, width = walk_width, res = walk_res, file = plot1)
  }
  pt = plot_isoform_counts(data = long.dt,
                           x = "new_name",
                           y = "Count",
                           plot_title = plot_title,
                           xlab = "",
                           ylab = "Count",
                           fill = "Sample",
                           filllab = "Sample",
                           type = "bar",
                           add_theme = counts_theme)
  ppng(print(pt),height = counts_height, width = counts_width, res = counts_res, file = plot2)
  if(!is.null(total_reads)) {
    sample_counts.dt = data.table(Sample = samples, total_read_counts = total_reads)
    long.dt = merge.data.table(long.dt, sample_counts.dt, by = "Sample", all.x = TRUE)
    long.dt[, TPM := Count / (total_read_counts / 1e6)]
    pt2 = plot_isoform_counts(data = long.dt,
                             x = "new_name",
                             y = "TPM",
                             plot_title = plot_title,
                             xlab = NULL,
                             ylab = "TPM",
                             fill = "Sample",
                             filllab = "Sample",
                             type = "bar",
                             add_theme = counts_theme)
    ppng(print(pt2),height = counts_height, width = counts_width, res = counts_res, file = plot3)
  }
  return(long.dt)
}


## function for plotting isoform counts
plot_isoform_counts = function(data, x, y, fill, plot_title, xlab = "", ylab = y, filllab = fill, add_theme = NULL, type = "bar") {
  if(type == "bar") {
    pt = ggplot(data = data, aes(x = get(x), y = get(y), fill = get(fill))) +
      geom_bar(stat = "identity", position = position_dodge()) +
      labs(title = plot_title,
           x = xlab,
           y = ylab,
           fill = filllab) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
            axis.text.y = element_text(size = 26),
            axis.title.y = element_text(size = 30),
            plot.title = element_text(size = 28, hjust = 0.5),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 16),
            legend.position="right", legend.direction="vertical",
            legend.key.size = grid::unit(0.5, "inch")) +
      facet_wrap(~new_name, scales = "free_x", nrow = 1) +
      theme(strip.text.x = element_text(size = 20))
  } else {
    stop("Only bar plots are currently implemented")
  }
  if(!is.null(add_theme)) {
    pt = pt + add_theme
  }
  return(pt)
}

  
  

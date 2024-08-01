## ## function to merge ggraph and fusions gwalk objects- could not get this to work well
## merge_junctions_gg_iso = function(fusions_gw, #fusions_gg from fusions task
##                                   iso_gw,  #gwalk after generating gwalks for all samples
##                                   pad = 1e3,
##                                   force_seqlengths = FALSE, #useful when gwalks do not have the same seqlengths
##                                   filter_type = "ALT") #can be a list of c("ALT","REF") to merge ref and alt junctions
## { 
##   if(inherits(fusions_gw,"character")) {
##     fusions_gw = readRDS(fusions_gw)
##   }
##   if(!inherits(fusions_gw, "gWalk")) {
##     stop("fusions_gw must be either a path to a gWalk or a gWalk object")
##   }
##   if(inherits(iso_gw,"character")) {
##     iso_gw = readRDS(iso_gw)
##   }
##   if(!inherits(iso_gw, "gWalk")) {
##     stop("iso_gw must be either a path to a gWalk or a gWalk object")
##   }
##   ## fusions junctions
##   ## fusions_gw$nodes$mark(fus_gw_id = paste0("ID",1:length(fusions_gw$nodes))) #add id, isoseq should already have a read
##   ## for each walk- need edges involved so I can mark the edges with walk ids
##   fusions_gw$set(fus_id = paste0("fus",1:length(fusions_gw)))
##   fusions.edges.dt = fusions_gw$dt[,.(fus_id, sedge.id)]
##   id_edges.dt = fusions.edges.dt[, .(sedge.id = unlist(sedge.id)), by = fus_id]
##   ## unique_edge.dt = id_edges.dt[, .(all_ids = list(unique(fus_id))), by = sedge.id]
##   unique_edge.dt = id_edges.dt[, .(all_ids = paste0((unique(fus_id)), collapse = ",")), by = sedge.id]
##   ## mark the edges and get the graph
##   mark.dt = fusions_gw$edgesdt[,.(sedge.id)]
##   mark.dt[, order := 1:.N]
##   mark.dt = merge.data.table(mark.dt,unique_edge.dt, by = "sedge.id", all.x = TRUE)
##   fusions_gw$edges$mark(all_ids = mark.dt$all_ids)
##   fusions.gg = fusions_gw$graph
##   ## fusions.jj = fusions.gg$junctions[type %in% filter_type,]
##   ## iso junctions
##   iso_gw$set(fus_id = paste0("Iso",1:length(iso_gw)))
##   iso.edges.dt = iso_gw$dt[,.(fus_id, sedge.id)]
##   id_edges.dt = iso.edges.dt[, .(sedge.id = unlist(sedge.id)), by = fus_id]
##   id_edges.dt[, all_ids := list(unique(fus_id)), by = "sedge.id"]
##   ## unique_edge.dt = id_edges.dt[, .(all_ids = list(unique(fus_id))), by = sedge.id]
##   unique_edge.dt = id_edges.dt[, .(all_ids = paste0((unique(fus_id)), collapse = ",")), by = sedge.id]
##   ## mark the edges and get the graph
##   mark.dt = iso_gw$edgesdt[,.(sedge.id)]
##   mark.dt[, order := 1:.N]
##   mark.dt = merge.data.table(mark.dt,unique_edge.dt, by = "sedge.id", all.x = TRUE)
##   iso_gw$edges$mark(all_ids = mark.dt$all_ids)
##   iso.gg = iso_gw$graph
##   ## iso.jj = iso.gg$junctions[type %in% filter_type,]
##   browser()
##   ## iso_gw$nodes$mark(iso_gw_id = paste0("ID",1:length(iso_gw$nodes)))
##   ## iso.gg = iso_gw$graph
##   ## iso.jj = iso.gg$junctions[type %in% filter_type,]
##   if(force_seqlengths) {
##     fusions.jj = jJ(gr.nochr(fusions.jj$grl))
##     iso.jj = jJ(gr.nochr(iso.jj$grl))
##   }
##   ## merge the junctions
##   merged.jj = gGnome::merge.Junction(fusions.jj, iso.jj, pad = pad,all = TRUE)
##   merged.jj$dt[,.(seen.by.ra1,seen.by.ra2)] %>% table
##   ## get the ids that are shared
##   merged.simple.dt = merged.jj$dt[,.(bp1.ra1,bp2.ra1,bp1.ra2,bp2.ra2, all_ids.ra1, all_ids.ra2, seen.by.ra1, seen.by.ra2)]
##   names(merged.simple.dt) = gsub(".ra1","_gg",names(merged.simple.dt))
##   names(merged.simple.dt) = gsub(".ra2","_iso",names(merged.simple.dt))
  
##   return()
  
## }



## get fusion data from isoseq task using fusion_breakpoints_groups
get_fus_data = function(fusion_breakpoints_groups, cores = 1) {
  fus.groups.dt = read_fusion_groups(fusion_breakpoints_groups)
  fus.groups.dt[, reads := tstrsplit(extra,"=",keep = 3, fixed = TRUE)]
  fus.groups.dt[, reads := tstrsplit(reads,";",keep = 1, fixed = TRUE)]
  fus.groups.dt[, reads2 := list(strsplit(reads,",",fixed=TRUE))]
  fus.groups.dt2 = copy(fus.groups.dt)
  fus.groups.dt = fus.groups.dt[, .(qname = unlist(reads2)), by = c("chr1", "start1", "end1", "chr2", "start2", "end2", "id", "score", "strand1", "strand2", "info", "row_number", "CB", "AC")]
  fus.sub.dt = fus.groups.dt2[,.(chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, row_number, CB, AC)]
  fus.sub.dt = unique(fus.sub.dt)
  fus.lst = mclapply(1:nrow(fus.sub.dt), function(x) {
    ## subset to get granges for each pair of breakpoints
    sub.dt = fus.sub.dt[x,]
    bp.id = sub.dt$id
    sub.dt[, reads := list(unique(fus.groups.dt[id == bp.id,]$qname))]
    sub.dt[, read_count := length(unique(fus.groups.dt[id == bp.id,]$qname))]
    sub.dt[, first_read := fus.groups.dt[id == bp.id,]$qname[1]]
    return(sub.dt)
  }, mc.cores = cores)
  fus.dt = rbindlist(fus.lst, fill = TRUE)
  fus.dt[, gene_names := sub(".*GN=([^;]+);.*", "\\1", info)]
  fus.dt[, gene_names := gsub("\\/","",gene_names)]
  fus.dt[, gene_names := gsub("NA,",",",gene_names)]
  fus.dt[, gene_names := gsub(",NA",",",gene_names)]
  max_genes = max(sapply(strsplit(fus.dt$gene_name, ","), length), na.rm = TRUE)
  fus.dt[, paste0("gene", 1:max_genes) := tstrsplit(gene_names, ",", fixed=TRUE, fill=NA)]
  fus.dt[gene1 == "", gene1 := NA]
  fus.dt[gene2 == "", gene2 := NA]
  if("gene2" %in% names(fus.dt)) {
    fus.dt[gene3 == "", gene3 := NA]
  }
  if("gene4" %in% names(fus.dt)) {
    fus.dt[gene4 == "", gene4 := NA]
  }
  fus.dt = fus.dt[order(-read_count),]
  fus.dt[, log_read_count := log(read_count)]
  fus.dt[, log10_read_count := log10(read_count)]
  return(fus.dt)
}

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
  ## browser()
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
    ## browser()
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
      ## gene.gw = get_iso_reads(bam,
      gene.gw = quick_ireads(bam,
                              gr = (gene.gr),
                              gtf = ref_gtf.gr,
                              cores = 1,
                              reannotate = FALSE,
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


## much cleaner version of get_iso_reads-will use already labeled transcripts-assigning a transcript to each read should be a separate function-it could then be implemented here much more easily
quick_ireads = function(bam,
                      gr,
                      gtf,
                      add_transcript_name = TRUE,
                      read_assignments = NULL, # data.table with qname and transcript_id, other meta data will be added to the gwalk, multiple transcript_ids per qname will only result in one walk by selecting the shorter reference transcript_id
                      type = "gw",
                      reannotate = FALSE,
                      remove_cigar_tags = c("S","D"),
                      unique_reads = TRUE,
                      subset_qname = NULL, #reads to subset on
                      filter_qname = NULL, # reads to filter out
                      pipeline = "isoquant",
                      cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  if(length(md.gr.bam) == 0) {
    warning("no reads in the specified region. Returning null")
    return(NULL)
  }
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  md.dt = md.dt[width != 0,]
  if(unique_reads) {
    md.dt = unique(md.dt)
  }
  if(!is.null(subset_qname)) {
    md.dt = md.dt[qname %in% subset_qname,]
  }
  if(!is.null(filter_qname)) {
    md.dt = md.dt[!(qname %in% filter_qname),]
  }
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
  message("Splicing cigar string")
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  message("Done Splicing cigar string")
  if(length(md.grl) == 0) {
    return(NULL)
  }
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  if ((!reannotate & pipeline == "isoquant")) {
    if(is.null(read_assignments)) {
      stop("Read_assignments now must be provided")
    }
    if(inherits(read_assignments,"character")) {
      reads.dt = isoquant_read_assignments(read_assignments)
    } else if (inherits(read_assignments,"data.table")) {
      message("User supplied read assignments as data.table. Continuing ...")
      reads.dt = read_assignments
    }
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    md.dt3 = merge.data.table(md.dt2, reads.dt, by = "qname", all = TRUE)[type != "N",]
    ## add check for multiple assigned transcripts of the same read-if so pick shortest
    check.dt = md.dt3[,.(qname,transcript_id)] %>% unique
    check.dt[, N_qname_tr_id := .N, by = "qname"]
    if(nrow(check.dt[N_qname_tr_id > 1,]) > 1) {
      check.sub.dt = check.dt[N_qname_tr_id > 1,]
      qnames_multiple_tr_ids = check.sub.dt$qname %>% unique
      warning(paste0("Some transcripts are assigned to multiple reads by isoquant. Will assign the shortest. Reads: ",qnames_multiple_tr_ids))
      md.lst = mclapply(unique(check.sub.dt$qname), function(x) {
        tr.test = check.sub.dt[qname == x,]$transcript_id
        tr.dt = as.data.table((gtf.gr %Q% (transcript_id %in% tr.test)) %Q% (type == "transcript"))
        tr.dt = tr.dt[order(width),]
        selected_transcript = tr.dt$transcript_id[1]
        sub.dt = md.dt3[qname == x & transcript_id == selected_transcript,]
        return(sub.dt)
      },mc.cores = cores)
      md.temp.dt = rbindlist(md.lst, fill = TRUE)
      md.dt3 = rbind(md.dt3[!(qname %in% md.temp.dt$qname),],md.temp.dt)
    }
  }
  ## md.dt3 - every read must now be assigned to a reference transcript
  if("associated_transcript" %in% names(md.dt3)) { ## for pipeline == "pacbio")
    md.full = md.dt3[structural_category == "full-splice_match" | structural_category == "full_splice_match",]
    md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
    md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
    if(nrow(md.novel) >0 ) {
      stop("Every read must now have a transcript assigned")
    }
  } else {
    md.full = md.dt3[structural_category == "full-splice_match" | structural_category == "full_splice_match",]
    md.sub = md.dt3[!(structural_category == "full-splice_match"| structural_category == "full_splice_match")]
    md.novel = md.dt3[structural_category != "full-splice_match" & is.null(transcript_id)]
    if(nrow(md.novel) >0 ) {
      stop("Every read must now have a transcript assigned")
    }
  }
  md.annotated = rbind(md.full,md.sub)
  
  if(class(gtf) == "character") {
    message(paste0("Reading in ",gtf)) 
    gtf = rtracklayer::readGFF(gtf) %>% as.data.table
    message(paste0("Done reading in ",gtf))
  }
  if(class(gtf)[1] != "data.table" & class(gtf)[1] != "GRanges") {
    stop("Gtf must be a character to a gtf file or a data.table")
  }
  if(any(class(gtf) != "GRanges")) {
    gtf.gr = GRanges(gtf, seqlengths = hg_seqlengths())
  } else {
    gtf.gr = gtf
  }
  ## md.annotated[seqnames != chr,] ## duplicate seqnames column
  md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
  md.annotated.gr = md.annotated.gr %Q% (!(type %in% remove_cigar_tags))
#############################NEWWWWWWWWWWWWWWWWWWWWWW
  ## now there are a few things to annotated for each read:
  ## (1) what type is every base pair (UTR, exon, intron, etc)
  ## (2) anything missing
  ## (3) add meta data
  ## get unique qnames and do it for every read
  unique_qnames = md.annotated.gr$qname %>% unique
  ## x = "m64466e_230825_205028/114753823/ccs"
  md.annotated.lst = mclapply(unique_qnames, function(x) {
    ## get the read alignments
    sub.gr = md.annotated.gr %Q% (qname == x)
    ## sub.gr = sub.gr %Q% (!(type %in% remove_cigar_tags))
    ## get the reference transcript
    ref.tr.id = unique(sub.gr$transcript_id)
    ref.tr.gr = ((gtf.gr %Q% (!is.na(transcript_id))) %Q% (transcript_id == ref.tr.id)) #remove NAs then subset granges
    ref.tr.gr = ref.tr.gr %Q% (type != "transcript") #this is the full reference coordiantes which I do not want for annotating
    ref.tr.gr$type = ref.tr.gr$type %>% as.character #this is sometimes a factor and messes up gr.val
    ## create tiled reference 
    tile.gr = gr.tile(ref.tr.gr, width = 1)
    tile.gr = gr.val(tile.gr,ref.tr.gr, "type")
    tile.dt = gr2dt(tile.gr)
    ## reduce the tile so no overlapping annotations
    tile.dt[, coding_type_vect := lapply(type, function(x) unlist(strsplit(x, ", ")))]
    tile.dt[, coding_type_simple := sapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                         ifelse("stop_codon" %in% x, "stop_codon",
                                                                         ifelse("UTR" %in% x, "UTR",
                                                                         ifelse("exon" %in% x, "exon",
                                                                         ifelse("CDS" %in% x, "exon","intron"))))))]
    tile.dt[, c("type","coding_type_vect") := NULL]
    tile.gr2 = GRanges(tile.dt, seqlengths = hg_seqlengths())
    tile.gr3 = gr.reduce(tile.gr2, by = "coding_type_simple", ignore.strand = FALSE)
    ##
    sub.gr2 = gr.breaks(query = sub.gr, bps = tile.gr3)
    ## add whether each exists
    ## sub.gr$mapped = TRUE
    sub.gr$mapped = "yes"
    sub.gr2 = gr.val(sub.gr2, sub.gr, "mapped") #riid should be in order but will have NAs so added a new order at the bottom here - need to switch order if transcript on negative strand
    ## new 7-31-24
    sub.gr2 = sub.gr2 %Q% (mapped != "")
    sub.gr2$mapped = NULL
    sub.gr3 = gr.val(sub.gr2, tile.gr3, "coding_type_simple")
    read.dt = gr2dt(sub.gr3)
    read.dt[is.na(coding_type_simple),coding_type_simple := "missing"]
    read.dt[coding_type_simple == "",coding_type_simple := "missing"]
    ## end new
    if(as.character(unique(ref.tr.gr@strand@values)) == "-") {
      read.dt = read.dt[, order := .N:1]
      read.dt = read.dt[order(order),]
    } else {
      read.dt[, order := 1:.N]
      read.dt = read.dt[order(order),]
    }
    return(read.dt)
  }, mc.cores = cores)
  reads.dt = rbindlist(md.annotated.lst, fill = TRUE)
  message("Done labeling each read to a transcript")
  reads.gr = GRanges(reads.dt, seqlengths = hg_seqlengths())
  reads.grl = rtracklayer::split(reads.gr, f = mcols(reads.gr)["qname"])
  message("Converted reads to grl. Adding meta data to grl")
  reads.dt
  meta.dt = copy(md.annotated)[, c("seqnames","start","end", "width","strand","type","rid","riid","fid") := NULL] %>% unique
  meta.dt[, qname := factor(qname, levels = names(reads.grl))]
  if(add_transcript_name) {
    meta.dt2 = merge.data.table(meta.dt, unique(gr2dt(gtf)[,.(transcript_id,transcript_name)]), by = "transcript_id")
    if(nrow(meta.dt2) > length(reads.grl)) {
      warning("Failed to add transcript_id. Likely more than one unique transcript_id and transcript_name combination in gtf")
      meta.dt2 = meta.dt
      add_transcript_name = FALSE
    }
  } else {
    meta.dt2 = meta.dt
  }
  meta.dt2 = meta.dt2[order(qname),]
  mcols(reads.grl) = meta.dt2
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    reads.gw = gW(grl = reads.grl)
    ## add coloring of nodes to match track.gencode
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "UTR_long","del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "red", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      reads.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
    if(add_transcript_name) {
      reads.gw$set(name = transcript_id)
    }
    message("done converting grl to gW object. returning gW")
    return(reads.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(reads.grl)
  }
}


## much cleaner version of get_ordered_fusions2 - way cooler too with labeling for each transcripts (exons, introns, UTR...)
quick_fusions = function(fusion_bam, #subset bam to only fusion reads
                         fusions_dt, #fusions data.table after running get_fus_data
                         first_read = TRUE,  #whether to just use the first read or all reads
                         pad = 1e4,                #pad for reading in bam
                         gtf,                      #gr gtf file with gene_id
                         remove_cigar_tags = c("S","D"),
                         sample_reads = FALSE,
                         sample_count = 10, #only used if sample reads is TRUE, requires read count column in fusions_dt
                         unique_reads = TRUE, #whether to grab unique after grabbing the reads
                         cores = 1)
{
  ## get the coordinates from the fusions_dt
  gr.dt1 = fusions_dt[,.(chr1,start1,end1)] %>% setnames(.,c("seqnames","start","end"))
  gr.dt2 = fusions_dt[,.(chr2,start2,end2)] %>% setnames(.,c("seqnames","start","end"))
  gr = rbind(gr.dt1,gr.dt2) %>% GRanges(.,seqlengths = hg_seqlengths())
  gr = gr + pad
  message(paste0("Reading in ",fusion_bam))
  fus.og.gr = bamUtils::read.bam(fusion_bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  if(unique_reads) {
    fus.og.gr = gr2dt(fus.og.gr) %>% unique %>% GRanges(., seqlengths = hg_seqlengths())
  }
  fus.gr = bamUtils::splice.cigar(fus.og.gr,get.seq = TRUE, rem.soft = FALSE, return.grl = FALSE)
  names(fus.gr) = NULL
  fus.gr = fus.gr %Q% (!(type %in% remove_cigar_tags))
  ## add back meta data from pre splicing cigar
  fus.meta.dt = as.data.table(mcols(fus.og.gr))
  fus.meta.dt[, new_row := 1:.N]
  names_add = unique(c("new_row",names(fus.meta.dt)[!(names(fus.meta.dt) %in% names(mcols(fus.gr)))]))
  fus.meta.dt = fus.meta.dt[,..names_add]
  add.meta.dt = data.table(new_row = fus.gr$rid)
  add.meta.dt = merge.data.table(add.meta.dt, fus.meta.dt, by = "new_row", all.x = TRUE, all.y = FALSE, allow.cartesian = FALSE)
  mcols(fus.gr) = cbind(mcols(fus.gr),add.meta.dt)
  if(length(fus.gr) == 0) {
    warning("no reads in the specified region. Returning null")
    return(NULL)
  }
  message("Done reading")
  if(first_read) { #if false subsetting to one read
    fus.gr = fus.gr %Q% (qname %in% fusions_dt$first_read)
  }
  if(sample_reads) {
    if(first_read) {
      stop("first_read must be FALSE if sample_reads is TRUE")
    }
    fusions_dt[read_count > sample_count,sample_reads := TRUE]
    fusions_dt[read_count <= sample_count,sample_reads := FALSE]
    ## sample if greater than 10
    fusions_dt[, sub_reads := lapply(.I, function(i) {
      if (sample_reads[i]) {
        list(sample(unlist(reads[[i]]), sample_count))
      } else {
        reads[i]
      }
    })]
    subset_reads_q = fusions_dt$sub_reads %>% unlist %>% unique
    fus.gr = fus.gr %Q% (qname %in% subset_reads_q)
  }
  ## process the multiple alignment differently so add a new id
  meta.dt = as.data.table(fus.gr[,c("qname","flag","mapq","cigar")])
  meta.dt[, id_alignment := .GRP, by = c("qname","flag","mapq","cigar")]
  meta.dt[, id_read := paste0(qname,"_",id_alignment)]
  fus.gr$id_read = meta.dt$id_read
  fus.gr$id_alignment = meta.dt$id_alignment
  fus.gr$og_order = 1:length(fus.gr)
  ## order these fusions to start- have to order the breakpoints - add which label
  split_cols = fusions_dt[, tstrsplit(AC, "/", fixed=TRUE)]
  col_names = paste0("bp", seq_along(split_cols))
  setnames(split_cols, col_names)
  split_cols[, bp1 := gsub("AC=","",bp1)]
  fusions_dt = cbind(fusions_dt,split_cols)
  ## fusions_dt[, AC := gsub("AC=", "", AC)]
  ## fusions_dt[, paste0("bp", 1:length(tstrsplit(AC[1], "/", fixed=TRUE))) := tstrsplit(AC, "/", fixed=TRUE)]
  ## get the reference gene identified
  fusions_dt[, GI := sub(".*GI=([^;]+);.*", "\\1", info)]
  gi_split = fusions_dt[, tstrsplit(GI, ",")]
  gi_col_names = paste0("gene_id", seq_along(gi_split))
  setnames(gi_split, gi_col_names)
  fusions_dt = cbind(fusions_dt,gi_split)
  message("Annotating each read with reference gene annotations")
  ## now annotate each alignment with each reference gene
  ## reads.lst = mclapply(unique(fus.gr$id_read)[c(1,2,4,5,6,8,9)], function(x) {
  reads.lst = mclapply(unique(fus.gr$id_read), function(x) {
    sub.gr = fus.gr %Q% (id_read == x)
    ## get which bp which also tells us which reference gene
    og_qname = sub.gr$qname %>% unique
    sub_fusions.dt = fusions_dt[sapply(reads, function(x) any(x %like% og_qname))]
    if(nrow(sub_fusions.dt) > 1) {
      warning("This read ",og_qname, "is found in multiple fusions inputted. Picking the first one. id = ", sub_fusions.dt[1,]$id, "; ", sub_fusions.dt[1,]$AC)
      sub_fusions.dt = sub_fusions.dt[1,]
    }
    ## which bp does this correspond to-not sure about more than two breakpoints- I guess I just add more pad here
    inter.gr1 = sub.gr %&% parse.gr(sub_fusions.dt$bp1)
    inter.gr2 = sub.gr %&% parse.gr(sub_fusions.dt$bp2)
    if(length(inter.gr1) == 0 & length(inter.gr2) == 0) {
      inter.gr1 = sub.gr %&% (parse.gr(sub_fusions.dt$bp1) + 500)
      inter.gr2 = sub.gr %&% (parse.gr(sub_fusions.dt$bp2) + 500)
      if(length(inter.gr1) == 0 & length(inter.gr2) == 0) {
        inter.gr1 = sub.gr %&% (parse.gr(sub_fusions.dt$bp1) + 1000)
        inter.gr2 = sub.gr %&% (parse.gr(sub_fusions.dt$bp2) + 1000)
      }
      if(length(inter.gr1) == 0 & length(inter.gr2) == 0) {
        sub.gr$bp_read = NA
        sub.gr$gene_id = NA
        sub.gr$breakpoint_number = NA
        warning("Read :", og_qname, "could not be split based on reference reads")
        sub.gr$coding_type_simple = "unassigned"
        ## have to subset the columns to what a working output will work with
        cols_keep = c("seqnames", "start", "end", "strand", "width", "type", "rid", "riid", "fid", "qname", "new_row", "flag", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual", "MD", "MQ", "id_read", "id_alignment", "og_order", "bp_read", "gene_id", "breakpoint_number", "qid", "breakpoint_id")
        sub.dt = gr2dt(sub.gr)
        cols_keep = cols_keep[cols_keep %in% names(sub.dt)]
        sub.dt = sub.dt[,..cols_keep]
        return(sub.dt)
      }
    }
    ## assign the breakpoint
    if(length(inter.gr1) != 0 & length(inter.gr2) != 0) { #both intersections non zero
      warning("This read ", og_qname, " with new id ", unique(sub.gr$id_read), "has multiple potential breakpoints. Selecting bp1,", sub_fusions.dt$bp1, "and reference gene ", sub_fusions.dt$gene_id1)
      bp_read = sub_fusions.dt$bp1
      gene_id_read = sub_fusions.dt$gene_id1
      gene_read = sub_fusions.dt$gene1
      breakpoint_number = "bp1"
    } else if (length(inter.gr1) != 0) {
      bp_read = sub_fusions.dt$bp1
      gene_id_read = sub_fusions.dt$gene_id1
      gene_read = sub_fusions.dt$gene1
      breakpoint_number = "bp1"
    } else if (length(inter.gr2) != 0) {
      bp_read = sub_fusions.dt$bp2
      gene_id_read = sub_fusions.dt$gene_id2
      gene_read = sub_fusions.dt$gene2
      breakpoint_number = "bp2"
    }
    sub.gr$bp_read = bp_read
    sub.gr$gene_id = gene_id_read
    sub.gr$breakpoint_number = breakpoint_number
    ## subset the gtf for this gene
    sub.gtf.gr = gtf %Q% (gene_name == gene_read)
    if(length(sub.gtf.gr) == 0) {
      warning("Read :", og_qname, "could not be split based on reference reads")
      sub.gr$coding_type_simple = "unassigned"
      ## have to subset the columns to what a working output will work with
      cols_keep = c("seqnames", "start", "end", "strand", "width", "type", "rid", "riid", "fid", "qname", "new_row", "flag", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize", "seq", "qual", "MD", "MQ", "id_read", "id_alignment", "og_order", "bp_read", "gene_id", "breakpoint_number", "qid", "breakpoint_id")
      sub.dt = gr2dt(sub.gr)
      cols_keep = cols_keep[cols_keep %in% names(sub.dt)]
      sub.dt = sub.dt[,..cols_keep]
      return(sub.dt)
    }
    ## construct non overlapping reference granges
    sub.gtf.gr = sub.gtf.gr %Q% (type != "transcript") #this is the full reference coordiantes which I do not want for annotating
    sub.gtf.gr = sub.gtf.gr %Q% (type != "intron")
    sub.gtf.gr$type = sub.gtf.gr$type %>% as.character #this is sometimes a factor and messes up gr.val
    ## create tiled reference 
    tile.gr = gr.tile(sub.gtf.gr, width = 1)
    tile.gr = gr.val(tile.gr,sub.gtf.gr, "type")
    tile.dt = gr2dt(tile.gr)
    tile.dt[, c("query.id","tile.id") := NULL]
    tile.dt = unique(tile.dt)
    ## reduce the tile so no overlapping annotations
    tile.dt[, coding_type_vect := lapply(type, function(x) unlist(strsplit(x, ", ")))]
    tile.dt[, coding_type_simple := sapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                         ifelse("stop_codon" %in% x, "stop_codon",
                                                                         ifelse("exon" %in% x, "exon",
                                                                         ifelse("UTR" %in% x, "UTR",
                                                                         ifelse("CDS" %in% x, "exon", "intron"))))))]
    ## had to change away from below which works well for isoform differences-using a different gtf here because it needs to be collaped 
    ## tile.dt[, coding_type_simple := sapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
    ##                                                                      ifelse("stop_codon" %in% x, "stop_codon",
    ##                                                                      ifelse("UTR" %in% x, "UTR",
    ##                                                                      ifelse("exon" %in% x, "exon",
    ##                                                                      ifelse("CDS" %in% x, "exon", "intron"))))))]
    tile.dt[, c("type","coding_type_vect") := NULL]
    tile.gr2 = GRanges(tile.dt, seqlengths = hg_seqlengths())
    tile.gr3 = gr.reduce(tile.gr2, by = "coding_type_simple", ignore.strand = FALSE)
    ## now break the read up
    sub.gr2 = gr.breaks(query = sub.gr, bps = tile.gr3)
    ## get only the ones that were actually mapped
    sub.gr$mapped = "yes"
    sub.gr2 = gr.val(sub.gr2, sub.gr, "mapped")
    sub.gr2 = sub.gr2 %Q% (mapped != "")
    sub.gr2$mapped = NULL
    ## 
    sub.gr3 = gr.val(sub.gr2, tile.gr3, "coding_type_simple")
    sub.dt = gr2dt(sub.gr3)
    sub.dt[, id_read := unique(sub.gr$id_read)]
    sub.dt[, bp_read := bp_read]
    sub.dt[, gene_id := gene_id_read]
    sub.dt[, breakpoint_id := breakpoint_number]
    sub.dt[is.na(coding_type_simple),coding_type_simple := "unassigned"]
    return(sub.dt)
  }, mc.cores = cores)
  reads.dt = rbindlist(reads.lst, fill = TRUE)
  ## ## add meta data back - don't need this, gr.breaks retains the meta data
  ## fus.dt = gr2dt(fus.gr)
  ## unique_meta.dt = fus.dt[,.(qname,id_alignment,cigar,mapq,seq, id_read)] %>% unique
  ## reads.dt2 = merge.data.table(reads.dt, unique_meta.dt, by = "id_read", all.x = TRUE)
  ## order the breakpoints
  reads.dt[, breakpoint_id := factor(breakpoint_id, levels = c("bp1","bp2","bp3","bp4","bp5","bp6"))]
  reads.dt[strand == "+", order_start := start]
  reads.dt[strand == "-", order_start := -start]
  reads.dt2 = reads.dt[order(breakpoint_id,order_start),]
  reads.dt2[is.na(coding_type_simple), coding_type_simple := "No_Value"]
  ## convert to granges and gwalk
  reads.gr = GRanges(reads.dt2, seqlengths = hg_seqlengths())
  reads.grl = rtracklayer::split(reads.gr, f = mcols(reads.gr)["qname"])
  message("Converted reads to grl. Adding meta data to grl")
  qnames_grl = data.table(qname = (reads.grl %>% names))
  mcols(reads.grl) = qnames_grl
  message("converting grl to gW object")
  reads.gw = gW(grl = reads.grl)
  ## add coloring of nodes to match track.gencode
  cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "UTR_long","del", "missing", "unassigned", "outside_ref"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "red", "orange", "orange", "blue", "blue"))
  for(x in 1:nrow(cmap.dt)) {
    cmap.sub = cmap.dt[x,]
    reads.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
  }
  return(reads.gw)
}

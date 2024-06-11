##qucikly rename a gwalk object to a specified first entry of a gwalk
## name.gw = function(gw, grl_col) {
##   new_names = sapply(gw$grl, function(x) mcols(x)[[grl_col]][1])
##   gw$set(name = new_names)
## }

## ## quickly rename a walk but also add ability to combine column from grl
## name.gw = function(gw, cols= NULL, div = "; ", col = NULL) {
##   ## if(!is.null(col)) {
##   ##   new_names = sapply(gw$grl, function(x) mcols(x)[[col]][1])
##   ##   gw$set(name = new_names)
##   ## } else if (!is.null(cols)) {
##   new_names_mat = sapply(cols, function(col) {
##     sapply(gw$grl, function(x) mcols(x)[[col]][1])
##   })
##   new_names = apply(new_names_mat, 1, function(row) {
##     paste(row, collapse = div)
##   })
##   names(new_names) = NULL
##   gw$set(name = new_names)
##   ## }
## }
## quickly rename a walk but also add ability to combine column from grl
name.gw = function(gw, cols= NULL, div = "; ", col = NULL) {
  new_names_mat = sapply(cols, function(col) {
    sapply(gw$grl, function(x) mcols(x)[[col]][1])
  })
  new_names = apply(new_names_mat, 1, function(row) {
    paste(row, collapse = div)
  })
  names(new_names) = NULL
  gw$set(name = new_names)
}

## function to reduce the transcripts plotted by combining the same transcript into one and returning the counts
collapse_same_tr = function(gw, keep_longest = TRUE, keep_incomplete_matches = TRUE) {
  if(!keep_incomplete_matches) {
    stop("this is not implemented are probably not going to be... Should keep all incomplete matches")
  }
  ##get the transcript assigned for each qname
  tr.qname.dt = lapply(gw$grl, function(x) as.data.table(unique(mcols(x)[c("qname","structural_category", "transcript_id")]))) %>% rbindlist %>% unique
  if(nrow(tr.qname.dt) != length(gw)) {
    stop("there is more than one transcript assigned for at least one qname or a nonunique entry based on qname and transcript_id")
  }
  tr.qname.dt[, ID := 1:.N]
  unique.trs.dt = tr.qname.dt[,.(structural_category,transcript_id)] %>% unique
  if(!(nrow(unique.trs.dt) > 0)) {
    warning("no full-splice_match to reduce")
    return(gw)
  }
  unique_tr.lst = unique.trs.dt[structural_category == "full-splice_match",]$transcript_id
  
  ids.keep.lst = lapply(unique_tr.lst, function(tr) {
    tr.collapse.dt = tr.qname.dt[transcript_id == tr & structural_category == "full-splice_match",]
    ##find which transcript has the largest width
    length.tr.dt = lapply(tr.collapse.dt$ID, function(x) {
      length.tr = gw$grl[[x]]@ranges@width %>% sum
      return(data.table(id = x, length_tr = length.tr))
    }) %>% rbindlist(.)
    keep.id = length.tr.dt[which.max(length_tr),]$id
    return(keep.id)
  }) %>% unlist
  ids.keep.lst = c(ids.keep.lst, tr.qname.dt[structural_category != "full-splice_match",]$ID) %>% sort %>% unique
  ##add counts of transcripts getting collapsed
  splice_match.dt = tr.qname.dt[structural_category == "full-splice_match",]
  splice_match.dt[, count := .N, by = c("structural_category", "transcript_id")]
  splice_match.dt2 = rbind(splice_match.dt,tr.qname.dt[structural_category != "full-splice_match",], fill = TRUE)
  splice_match.dt2[is.na(count), count := 1]  
  gw.grl = gw$grl
  sub.dt = lapply(ids.keep.lst, function(x) {
    sub.gr = gw.grl[[x]]
    sub.gr$count = splice_match.dt2[ID == x,]$count
    return(as.data.table(sub.gr))
  }) %>% rbindlist(.)
  sub.gr = GRanges(sub.dt, seqlengths = hg_seqlengths())
  sub.grl = rtracklayer::split(sub.gr, f = mcols(sub.gr)["qname"])
  sub.gw = gW(grl = sub.grl)
  ## setup gw same as in get_iso_reads
  ## add coloring of nodes to match track.gencode
  cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
  for(x in 1:nrow(cmap.dt)) {
    cmap.sub = cmap.dt[x,]
    sub.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
  }
  transcript_names = sapply(sub.gw$grl, function(x) x$transcript_name[1])
  sub.gw$set(name = transcript_names)
  return(sub.gw)
}


## generate walks from gff
isoseq_gff2gw = function(collapsed.gff, classification, chr = TRUE,subset_class = c("isoform", "structural_category", "associated_gene", "associated_transcript", "coding", "subcategory"), save_grl = NULL,save_gw = NULL, type = "walk", color = NULL, column = "structural_category", sorted = TRUE, filter_unique_isoforms = TRUE, filter_chromosomes = paste0("chr",c(1:22,"X","Y"))) {
  gff.dt = readGFF(collapsed.gff) %>% as.data.table
  names(gff.dt) = c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "gene_id", "isoform")
  class.dt = fread(classification)
  if(!is.null(subset_class)) {
    class.dt = class.dt[,..subset_class]
  }
  merged.dt = merge.data.table(gff.dt, class.dt, by = "isoform", all.x = TRUE)
  ##need to sort positive strand in order and negative strand in opposite of chromosome order
  merged.dt = merged.dt[order(get(column)),]
  pos.dt = merged.dt[strand == "+",][order(seqnames,start,end),]
  neg.dt = merged.dt[strand == "-",][order(seqnames,-start,-end),]
  merged.dt = rbind(pos.dt,neg.dt)
  merged.dt = merged.dt[type == "exon",]
  if(filter_unique_isoforms) {
    merged.dt[, length := .N, by = "isoform"]
    merged.dt = merged.dt[length != 1,]
  }
  if(!is.null(filter_chromosomes)) {
    merged.dt = merged.dt[seqnames %in% filter_chromosomes,]
  }
  merged.gr = GRanges(merged.dt,seqlengths = hg_seqlengths(genome = "BSgenome.Hsapiens.UCSC.hg38",chr = chr))
  message("Converting to grl")
  merged.grl = gr2grl2(merged.gr, "isoform") %>% GRangesList(.)
  message("Done converting")
  if(!is.null(save_gw) || type == "walk") {
    gff.gw = gW(grl = merged.grl)
    if(!is.null(color)) {
      for(x in 1:nrow(color)) {
        color1 = color[x,color]
        category = color[x,category]
        gff.gw$nodes[get(column) == category]$mark(col = color1)
      }
    }
  }
  if(!is.null(save_gw)) {
    message(paste0("Saving gwalks to ",save_gw))
    saveRDS(gff.gw, save_gw)
  }
  if(!is.null(save_grl)) {
    message(paste0("Saving grl to ",save_grl))
    saveRDS(merged.grl, save_grl)
  }
  if(type == "walk") {
    return(gff.gw)
  }
  if(type == "gr") {
    return(merged.gr)
  }
  if(type == "grl") {
    return(merged.grl)
  }
  if(type == "table") {
    return(as.data.table(merged.gr))
  }
}

## genearte grl from gr based on id column
gr2grl2 = function(gr, ID) {
  return(split(gr, f = mcols(gr)[[ID]]))
}


## parallel version of gr2grl2
gr2grl3 = function(gr, ID, cores = 1) {
  unique_values = unique(mcols(gr)[[ID]])
  gr_list = mclapply(unique_values, function(value) {
    subset_gr = gr[mcols(gr)[[ID]] == value]
    return(subset_gr)
  }, mc.cores = cores)
  
  return(gr_list)
}


## ## reorganize this function- should have same labeling for previously annotated and newly annotated; novel should just return the matching transcripts
get_iso_reads = function(bam,
                         gr,
                         gtf,
                         collapsed_group = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         collapsed_class = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for
                         read_assignments = NULL, # isoquant read assignements output, if reannotate = FALSE and this, will look for
                         type = "gw",
                         annotate_mismatch = TRUE,
                         annotate_mismatch_type = "smart",
                         reann_with = "both", #whether to match using both exons and UTRS or just exon ("exon")
                         reannotate = FALSE,
                         select_type = "bases", #only bases works for annotate_mismatch_type = "smart"
                         add_gencode_name = TRUE, #add the gencode name to the transcript for plotting
                         reverse_overlap = FALSE, #old method - need to remove
                         consider = "only_exon_alignment", #only applies to select_type = "fraction"
                         read_over_potential = 0.5, #only applies to annotate_mismatch_type = "percent"
                         potential_over_read = 0.8, #only applies to annotate_mismatch_type = "percent"
                         minimize = "tr_bases_missing", #only applies to annotate_mismatch_type = smart, typically works better then sum
                         annotate_missing = TRUE,
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
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
  if(!reannotate & pipeline == "pacbio") {
    ##look for collapsed_group and collapsed_class if missing
    if(is.null(collapsed_group) | is.null(collapsed_class)) {
        if(is.null(collapsed_group) && is.null(collapsed_class)) {
          message("looking for collapsed_group and collapsed class")
          folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
          collapsed_group = paste0(folder_find,"collapsed.group.txt")
          collapsed_class = paste0(folder_find,"collapsed_classification.txt")
          if(all(file.exists(c(collapsed_group,collapsed_class)))) {
            message("Found! collapsed.group.txt and collapsed_classification.txt")
          }
        }
    }
    group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
    class.dt = fread(collapsed_class)
    ## get labels for already identified transcripts
    group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
    group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
    group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
    expanded.dt = copy(group.sub.dt)
    expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
    group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
      md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
    } else {
      md.dt3 = group.sub.dt2
    }
    
  } else if ((!reannotate & pipeline == "isoquant")) {
    if(is.null(read_assignments)) {
      message("looking for read assignments")
      folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-2)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
      read_assignments = paste0(folder_find,"OUT.read_assignments.tsv.gz")
      if(file.exists(read_assignments)) {
        message("Found! OUT.read_assignments.tsv.gz")
      }
    }

    ## reads.dt = fread(read_assignments)
    ## ## very slow reads.dt = read.delim(read_assignments, header = TRUE, comment.char = "#")
    ## names(reads.dt) = gsub("#","",reads.dt[1,])
    ## reads.dt = reads.dt[2:nrow(reads.dt),]
    ## reads.dt = reads.dt[read_id %in% md.dt$qname,]
    ## names(reads.dt) = c("qname", "chr", "strand", "transcript_id", "gene_id", "assignment_type", "assignment_events", "exons", "additional_info")
    ## reads.dt[, structural_category := tstrsplit(additional_info, "; Classification=", keep = 2)]
    ## reads.dt[, structural_category := gsub(";","", structural_category)]
    ## names(reads.dt)[names(reads.dt) == "strand"] = "strand_classification"
    ## names(reads.dt)[names(reads.dt) == "chr"] = "chr_classification"
    ##
    if(inherits(read_assignments,"character")) {
      reads.dt = isoquant_read_assignments(read_assignments)
    } else if (inherits(read_assignments,"data.table")) {
      message("User supplied read assignments as data.table. Continuing ...")
      reads.dt = read_assignments
    }
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
      md.dt3 = merge.data.table(md.dt2, reads.dt, by = "qname", all = TRUE)[type != "N",]
    } else {
      md.dt3 = reads.dt
    }
  } else if(reannotate) {
    ## cols_keep = c("qname", "seqnames", "start", "end", "width", "strand", "type")
    ## cols_keep = cols_keep[cols_keep %in% names(md.dt3)]
    md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
  }
  if(annotate_mismatch) {
    if(!reannotate) {
      if("associated_transcript" %in% names(md.dt3)) { ## for pipeline == "pacbio")
        md.full = md.dt3[structural_category == "full-splice_match" | structural_category == "full_splice_match",]
        md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
        md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
      } else {
        md.full = md.dt3[structural_category == "full-splice_match" | structural_category == "full_splice_match",]
        md.sub = md.dt3[!(structural_category == "full-splice_match"| structural_category == "full_splice_match")]
        md.novel = md.dt3[structural_category != "full-splice_match" & is.null(transcript_id)]
      }
    } else if(reannotate) {
      md.novel = md.dt3
      md.full = NULL
    }
###################################
    ## first  label ones without a transcript annotation
    if(nrow(md.novel) == 0) {
      md.novel.gr = GRanges()
      } else {
        md.novel.gr = GRanges(md.novel)
      }
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
    ##add potential genes to md.novel
    if(length(md.novel.gr) > 0) {
      md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
      potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
      potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
      ## annotate to the longest transcript??
      if(annotate_mismatch_type == "longest") {
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "level") {
        ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "match") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        ##potential one
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
        ## get the transcript with the most N by transcript
        potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "percent") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              if(consider == "only_exon_alignment") {
                md.sub.gr = md.sub.gr %&% potential.sub.gr
              }
              md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
            } else if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                md.sub.gr$percent = x$width.ov
              } else {
                md.sub.gr$percent = 0
              }
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
              ## dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "fraction") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr2 := rank(pot_over_tr,-mean_overlap,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        if(select_type == "bases") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        ## need to rewrite to loop through the correct combinations
        if(reverse_overlap) { # this option should be removed- ultimately in the loop above
          message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
          rows.cycle = mean.overlap.dt[rank_tr <= 20,]
          mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
            rows.cycle.sub = rows.cycle[x,]
            tr = rows.cycle.sub$query
            pot.tr = rows.cycle.sub$subject
            md.sub.gr = md.novel.grl[[tr]]          
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
            } else if(select_type == "bases") {
              ## potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = potential.sub.gr
              y = md.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                potential.sub.gr$percent = x$width.ov
              } else {
                potential.sub.gr$percent = 0
              }

            }
            ##try returning the mean
            return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          if(select_type == "fraction") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
            mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
          if(select_type == "bases") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
            mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
        } else {
          mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
        if(is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
        }
        if(!is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
          add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
          md.novel.gr$gene_name = NULL
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
          if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
            warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
          }
        }
        ## now annotation the exons overlaps
        md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
        potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      } else if(annotate_mismatch_type == "smart") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        ## add section to only do overlap with overlapping transcripts - new 5-24-24
        ## only consider ones with any overlap for each transcript
        matched.potential.tr.dt = mclapply(names(md.novel.grl), function(tr) {
          tr.gr = md.novel.grl[[tr]]
          potent.tr = (potential_transcripts.gr2 %&% tr.gr)$transcript_id %>% unique
          dt1 = data.table(transcript_id = tr, potential_transcripts = list(potent.tr))
          return(dt1)
        },mc.cores = cores)
        matched.potential.tr.dt = rbindlist(matched.potential.tr.dt, fill = TRUE)
        ##
        ## get maximum overlap of each transcript with each potential transcript
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          potential_transcripts_tr = unlist(matched.potential.tr.dt[transcript_id == tr,]$potential_transcripts)
          ## mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
          mean.overlap.dt = lapply(potential_transcripts_tr, function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0) {
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                ##new part for annotate_mismatch_type == "smart"
                y_sum_bases = y@ranges@width %>% sum
                x_sum_bases = x@ranges@width %>% sum
                ov_sum_bases = x$width.ov %>% sum
                ## ov_sum_bases = ov$width %>% sum
                potential_leftout = y_sum_bases - ov_sum_bases
                tr_leftout = x_sum_bases - ov_sum_bases
                ##end new
              } else {
                tr_leftout = "ALL"
                potential_leftout = "ALL"
              }
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = tr_leftout, pot_over_tr = potential_leftout)
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "bases") {
          mean.overlap.dt[, mean_overlap := as.numeric(mean_overlap)]
          mean.overlap.dt[, pot_over_tr := as.numeric(pot_over_tr)]
          mean.overlap.dt[, sum_missing := (mean_overlap + pot_over_tr)]
          mean.overlap.dt[, rank_tr := rank(sum_missing,ties.method = "min"), by = "query"]
          mean.overlap.dt = mean.overlap.dt[order(mean_overlap, rank_tr),]
          if(minimize == "sum") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
          }
          if(minimize == "tr_bases_missing") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
          }
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
      if(is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
      }
      if(!is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
        add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
        md.novel.gr$gene_name = NULL
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
        if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
          warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
        }
      }
      ## now annotation the exons overlaps
      md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
      potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      }
    }
    if(exists("md.novel.dt2") & !reannotate) {
    gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
    md.annotated = rbind(md.full,md.sub)
    md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
    md.annotated[, transcript_id := associated_transcript]
    md.novel.dt2[, reannotated := TRUE]
    md.annotated[, reannotated := FALSE]
    md.annotated = rbind(md.novel.dt2,md.annotated)
    } else if (reannotate) {
      md.annotated = md.novel.dt2
      md.annotated[, reannotated := TRUE]
    } else {
      gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique
      md.annotated = rbind(md.full,md.sub)
      if(nrow(gtf.possible.tr.dt) > 0) {
        md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
        md.annotated[, transcript_id := associated_transcript]
        md.annotated[, reannotated := FALSE]
      } else {
      }
    }
    ## md.annotated[seqnames != chr,] ## duplicate seqnames column
    md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
    if(pipeline == "pacbio") {
      potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
      potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
      add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
    } else {
      potential.gtf.gr = gtf.gr %Q% (transcript_id %in% md.annotated.gr$transcript_id)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
    }
    md.annotated.dt2 = md.annotated
    md.annotated.dt2 = md.annotated.dt2[!(type %in% remove_cigar_tags),]
    md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
    ## now annotation the exons overlaps
    md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
    potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
    potential.transcripts = potential.gtf.labels$transcript_id %>% unique
    message("breaking reads to label each element based on reference transcripts")
    md.annotated.lst3 = mclapply(potential.transcripts, function(tr) {
      md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
      ## md_potential_tr.sub.gr = md_potential_tr.gr %Q% (qname == "m64466e_230629_174835/100076360/ccs")
      potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
      potential.gtf.labels.sub = potential.gtf.labels.sub %Q% (coding_type %in% c("CDS","UTR"))
      ##breaks isn't working properly-convert these to breakpoints
      potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end,strand, coding_type)]
      ##test end + 1 for pos
      bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt)))) %>% unique
##      bps.dt = copy(potential.gtf.labels.sub.dt)
      bps.dt[,start := end -1]
      bps.dt[start_end == "start", end := end - 1]
      bps.dt[start_end == "start", start := start - 1]
      bps.dt[start_end == "end", end := end]
      bps.dt[start_end == "end", start := start + 1]
      ## ## ##test labels
      ## ## ## bps.dt[,start := start + 1]
      ## ## ## bps.dt[,end := start]
      ## ## ##
      ## bps.dt[,start_end := NULL]
      ## bps.dt = unique(bps.dt)
      ## bps.dt[,strand := NULL]
      ## ## ## test no end
      ## bps.dt[, end := end]
      ## bps.dt[, start := end]
      ## bps.dt[,start_end := NULL]
      ## bps.dt = unique(bps.dt)
      bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths(chr = FALSE)) %>% gr.chr
      ## #################################################################################### 6-5-24
      ## ## new attempt at fixing breakpoints- don't have coding sequence overlap - exon should not overlap UTR
      ## ## gene transcript exon CDS start_codon stop_codon UTR


      ## potential.gtf.labels.sub$orig_order = 1:length(potential.gtf.labels.sub)
      ## ## granges for each type of the gene
      ## exon.gr = potential.gtf.labels.sub %Q% (type == "exon")
      ## transcript.gr = potential.gtf.labels.sub %Q% (type == "transcript")
      ## CDS.gr = potential.gtf.labels.sub %Q% (type == "CDS")
      ## start_codon.gr = potential.gtf.labels.sub %Q% (type == "start_codon")
      ## stop_codon.gr = potential.gtf.labels.sub %Q% (type == "stop_codon")

      ## potential.gtf.labels.sub
      ## tiled.gr = gr.tile(potential.gtf.labels.sub,1)
      ## tiled.gr = gr.val(tiled.gr,potential.gtf.labels.sub, "coding_type", mean = FALSE)
      ## tiled.dt = as.data.table(tiled.gr)
      ## tiled.dt[,coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
      ## tiled.dt[, coding_type_simple := sapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
      ##                                                                              ifelse("stop_codon" %in% x, "stop_codon",
      ##                                                                               ifelse("UTR" %in% x, "UTR",
      ##                                                                               ifelse("exon" %in% x, "exon","intron")))))]
      ## tiled.gr2 = GRanges(tiled.dt, seqlengths = hg_seqlengths())
      ## bps.gr = gr.reduce(tiled.gr2, by = c("coding_type_simple"))
      ## bps.gr$coding_type = bps.gr$coding_type_simple
      ## test.gr = potential.gtf.labels.sub %Q% (coding_type %in% c("CDS","UTR"))
      ## test.gt = gTrack(potential.gtf.labels.sub,gr.colorfield = "type")
      ## test.gt2 = gTrack(test.gr,gr.colorfield = "type")
      ## test.gt2 = gTrack(potential.gtf.labels.sub,gr.colorfield = "coding_type")
      ## test.gt3 = gTrack(tiled.gr3,gr.colorfield = "coding_type_simple")
      ## test.gt4 = gTrack(bps.gr2)
      ## ppng(plot(c(test.gt,test.gt2), gene.gr),height = 6000, width = 4000, res = 300)
      ## ppng(plot(c(test.gt,test.gt2,test.gt3), gene.gr),height = 6000, width = 4000, res = 300)
      ## ppng(plot(c(test.gt,test.gt2,test.gt3, test.gt4), gr.reduce(bps.gr2+5)),height = 6000, width = 4000, res = 300)




      ## bps.dt = as.data.table(bps.gr)[order(start),]
      ## bps.dt2 = data.table(seqnames = unique(bps.dt$seqnames), start = c((bps.dt$start[1]-1),bps.dt$end), end = c((bps.dt$start[1]),bps.dt$start[2:nrow(bps.dt)],(bps.dt$end[nrow(bps.dt)]+1)))

      ## bps.dt2 = data.table(seqnames = unique(bps.dt$seqnames), start = c((bps.dt$start -2)))
      ## bps.dt2[,end := start]



      ## bps.dt = as.data.table(sort(bps.gr))
      ## bps.dt2 = data.table(seqnames = unique(bps.dt$seqnames), start = c((bps.dt$start[1]-1),bps.dt$end), end = c((bps.dt$start[1]),bps.dt$start[2:nrow(bps.dt)],(bps.dt$end[nrow(bps.dt)]+1)))
      ## test.dt = bps.dt[coding_type_simple == "exon",]
      ## test.dt$start[1] = test.dt$end[1]
      ## test.dt$end[7] = test.dt$start[7]
      ## test.dt2 = bps.dt[coding_type_simple != "exon",]
      ## test.dt2 = data.table(seqnames = "chr12", start = c(57747727, 57748524,57748527,57751714,57751717,57751736,57752310))
      ## test.dt2[, end := start]
      ## test.gr = GRanges(rbind(test.dt,test.dt2,fill = TRUE),seqlengths = hg_seqlengths())

##      bps.dt = data.table(seqnames = c(bps.dt$seqnames, bps.dt$seqnames), end = c(bps.dt$start,bps.dt$end), strand = as.character(bps.dt$strand, bps.dt$strand), start_end = c(rep("start",nrow(bps.dt)),rep("end",nrow(bps.dt)))) %>% unique
      bps.dt = copy(potential.gtf.labels.sub.dt)
      bps.dt = data.table(seqnames = rep(bps.dt$seqnames,2), end = c(bps.dt$start,bps.dt$end), strand = as.character(bps.dt$strand, bps.dt$strand), start_end = c(rep("start",nrow(bps.dt)),rep("end",nrow(bps.dt)))) %>% unique
      bps.dt[, start := end]
      bps.gr2 = GRanges(bps.dt, seqlengths = hg_seqlengths())
      md_potential_tr.sub.gr = gr.chr(md_potential_tr.sub.gr)
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
      md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr2,query = md_potential_tr.sub.gr)
      ## ## think I need to do the breaks for each individual qname because gr.breaks is not working for all at once
      potential.gtf.labels.sub.gr = GRanges(potential.gtf.labels.sub.dt, seqlengths = hg_seqlengths(chr = TRUE))
      md.novel.lst3 = mclapply(unique(md_potential_tr.gr$qname), function(x) {
        ## tr.sub.gr = md_potential_tr.gr %Q% (qname == x)
        tr.sub.gr = md_potential_tr.sub.gr %Q% (qname == x)
        ##tr.sub.gr2 = gUtils::gr.breaks(bp = bps.gr2,tr.sub.gr)
        tr.sub.gr2 = gUtils::gr.breaks(bps = potential.gtf.labels.sub.gr,query = tr.sub.gr)
        tr.sub.gr2 = gr.val(tr.sub.gr2,potential.gtf.labels.sub,"coding_type")
        return(as.data.table(tr.sub.gr2))
      }, mc.cores = cores)
      ##md_potential_tr.sub.gr = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      ## md_potential_tr.sub.gr = gr.val(md_potential_tr.sub.gr,bps.gr,"coding_type")
      ##md.novel.dt3 = as.data.table(md_potential_tr.sub.gr)
      md.novel.dt3 = rbindlist(md.novel.lst3)
      return(md.novel.dt3)
####################################################################################################################################################################################################################################################################################
      ## fix bps.dt
      ## do not use first CDS end as breakpoint
##      if(potential.gtf.labels.sub.dt$strand[1] == "-") {
      ##        potential.gtf.labels.sub.dt = potential.gtf.labels.sub.dt[order(-start),]
      ## bps.dt = potential.gtf.labels.sub.dt[coding_type != "exon",]
      ## bps.dt[, coding_type_entry := 1:.N, by = "coding_type"]
      
      ## bps.dt = data.table(seqnames = c(bps.dt$seqnames, bps.dt$seqnames), start = c(bps.dt$start, bps.dt$end))
      ## bps.dt[, end := start,]
      ## bps.dt = unique(bps.dt)
      ## bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths(chr = FALSE)) %>% gr.chr
      ## ## starts = potential.gtf.labels.sub.dt$start
      ## ## end = potential.gtf.labels.sub.dt$start
      ## ## } else {
      ## ##   potential.gtf.labels.sub.dt = potential.gtf.labels.sub.dt[order(start),]
      ## ## }

      ## bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt))), type = )
####################################################################################################################################################################################################################################################################################
      ## md_potential_tr.sub.gr = gUtils::gr.breaks(bp = bps.gr2,md_potential_tr.sub.gr)
      ## new gr.val to match to specific transcripts
      ## potential.gtf.labels.sub = potential.gtf.labels.sub %+% 1 #Have to add one to the coordinates when doing gr.val here


      ## ############### newer commented out
      ## potential.gtf.labels.sub = sort(potential.gtf.labels.sub)
      ## ##sort((md_potential_tr.sub.gr %Q% (qname == "m64466e_230629_174835/109120158/ccs")) %&% (potential.gtf.labels.sub[9]))
      ## ##sort((md_potential_tr.sub.gr %Q% (qname == "m64466e_230629_174835/109120158/ccs")) %&% (bps.gr[15:16]))
      ## ## (md_potential_tr.sub.gr %Q% (qname == "m64466e_230629_174835/109120158/ccs")) %&% (potential.gtf.labels.sub[2:5])
      ## ## (potential.gtf.labels.sub[8:9])
      ## md_potential_tr.sub.gr = sort(md_potential_tr.sub.gr)
      ## md_potential_tr.sub.gr = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
      ## md_potential_tr.sub.gr = gUtils::gr.breaks(bp = (potential.gtf.labels.sub %Q% (coding_type == "UTR")),md_potential_tr.sub.gr)
      ## md_potential_tr.sub.gr = gUtils::gr.breaks(bp = (potential.gtf.labels.sub %Q% (coding_type == "CDS")),md_potential_tr.sub.gr)

      ## md_potential_tr.sub.gr = gUtils::gr.breaks(bp = (potential.gtf.labels.sub[8:9]),md_potential_tr.sub.gr)
      

      ## md_potential_tr.sub.gr = gr.val(md_potential_tr.sub.gr,potential.gtf.labels.sub,"coding_type")
      ## ## md_potential_tr.sub.gr = gr.val(md_potential_tr.sub.gr,bps.gr,"coding_type")
      ## md.novel.dt3 = as.data.table(md_potential_tr.sub.gr)

      ## end new commented out
######################################################################################################################################################################################################################################
      ## ## new attempt at these breaks-think the reason it does not break is that CDS overlaps stop_codon/ start_codon
      ## bps.no.st.gr = bps.gr %Q% (!type %in% c("start_codon","stop_codon"))
      ## bps.st.gr = bps.gr %Q% (type %in% c("start_codon","stop_codon"))
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.no.st.gr,md_potential_tr.sub.gr)
      ## potential.gtf.labels.sub$orig_order = 1:length(potential.gtf.labels.sub)
      ## potential.gtf.labels.sub.no.st = potential.gtf.labels.sub %Q% (!type %in% c("start_codon","stop_codon"))
      ## potential.gtf.labels.sub.st = potential.gtf.labels.sub %Q% (type %in% c("start_codon","stop_codon"))
      ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub.no.st,"coding_type")

      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.no.st.gr,md_potential_tr.gr2)
      ## potential.gtf.labels.sub.st$coding_type_st = potential.gtf.labels.sub.st$coding_type
      ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub.st,"coding_type_st")

      
      ## cds_utr.gr = potential.gtf.labels.sub %Q% (!type %in% c("start_codon","stop_codon"))
      ## ## st.gr = potential.gtf.labels.sub %Q% (type %in% c("start_codon","stop_codon"))

      
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
      ## ## new gr.val to match to specific transcripts
      ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      ## md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      ## potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)
      ## potential.gtf.labels.sub.dt[, og_order := 1:.N]
      ## potential.gtf.labels.sub = GRanges(potential.gtf.labels.sub.dt, seqlengths = hg_seqlengths())
      ## cds_utr.gr = potential.gtf.labels.sub %Q% (!type %in% c("start_codon","stop_codon"))
      ## st.gr = potential.gtf.labels.sub %Q% (type %in% c("start_codon","stop_codon"))
      ## gr.breaks(query = cds_utr.gr, bp = st.gr)
      ## ## cds_utr.gr2 = cds_utr.gr[cds_utr.gr %outside% st.gr]
      ## bp.gr = c(cds_utr.gr2,st.gr)
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub, md_potential_tr.sub.gr)
      ## ## new gr.val to match to specific transcripts
      ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      ## md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      ## return(md.novel.dt3)
    }, mc.cores = 1)
    md.annotated.dt3 = rbindlist(md.annotated.lst3)
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
################
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
    ## md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
    ##                                                                               ifelse("UTR" %in% x, "UTR",
    ##                                                                               ifelse("exon" %in% x, "exon","intron"))))]
    md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                                  ifelse("UTR" %in% x, "UTR",
                                                                                  ifelse("exon" %in% x, "exon",
                                                                                  ifelse("CDS" %in% x, "exon","intron")))))]
    ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3 = md.annotated.dt3[!(type %in% remove_cigar_tags),]
    ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3[type == "X",coding_type_simple := "del"]
    md.all.dt = md.annotated.dt3
    message("finished merging together")
    
    ## rbind all together
    ## sort based on strand to have walks align
    pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
    neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
    md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
    if(add_gencode_name) {
      gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
      gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
      md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
    }
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    ## convert to grl and gw
    message("converting to grl")
    md.gr2$qr_tr = paste0(gsub("transcript/","",md.gr2$qname),"; ", md.gr2$transcript_name)
    md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
    message("done coverting to grl")
    ## now annotating missing parts
    if(annotate_missing) {
      message("annotating missing regions from transcripts")
      md.all.dt2 = mclapply(1:length(md.grl2), function(x) {
        md.sub.tr = md.grl2[[x]]
        qr_tr1 = md.sub.tr$qr_tr[1]
        qname_tr1 = md.sub.tr$qname[1]
        uniq.tr = unique(md.sub.tr$transcript_id)
        if(length(uniq.tr) == 1) {
          tr.gr = potential.gtf.dt[transcript_id == uniq.tr,][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        } else {
          ## pick the longest if ambiguous
          uniq.lst = lapply(1:length(uniq.tr), function(x) {
            uniq.tr1 = uniq.tr[x]
            tr.gr = potential.gtf.dt[transcript_id == uniq.tr1,][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
            return(data.table(entry = x, tr = uniq.tr1, length = sum(tr.gr@ranges@width)))
          })
          uniq.dt = rbindlist(uniq.lst)[order(-length),]
          tr1 = uniq.dt$tr[1]
          tr.gr = potential.gtf.dt[transcript_id == tr1,][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        }
        ## tr.gr = potential.gtf.dt[transcript_id == unique(md.sub.tr$transcript_id),][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
                                        #tr.missing.gr = tr.gr[tr.gr %outside% md.sub.tr]
                                        #mcols(tr.missing.gr) = NULL
        ##attempt not using outside
        tiled.gr = gr.tile(tr.gr,width = 1)
        tiled.gr$overlap = tiled.gr %O% md.sub.tr
        ## tr.missing.gr = (tiled.gr %Q% (overlap < 1)) %>% gr.reduce
        tr.missing.gr = (tiled.gr %Q% (overlap < 1))
        ## add annotation for exon from tr.missing!
        tr.missing.gr = gr.val(tr.missing.gr,tr.gr, c("exon_id","exon_number"))
        tr.missing.gr = gr.reduce(tr.missing.gr,by = c("exon_id"))
        ## tr.gr$type = as.character(tr.gr$type)
        ##
        if(length(tr.missing.gr) == 0) {
          return(as.data.table(md.sub.tr))
        }
        tr.missing.gr = as.data.table(tr.missing.gr) %>% GRanges(., seqlengths = hg_seqlengths())
        ##end  attempt
        ##add meta data
        meta.add.dt = mcols(md.sub.tr)[1,] %>% as.data.table
        ##meta.add.dt[, c("qid", "type","rid","riid","fid","col",) := NULL]
        if("isoform" %in% names(meta.add.dt)) {
          meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
        } else if ("gene_name" %in% names(meta.add.dt)) {
          meta.add.dt = meta.add.dt[,.(transcript_id, qname, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
        } else {
          meta.add.dt = meta.add.dt[,.(transcript_id, qname, transcript_type, transcript_name, qr_tr)]
        }
        tr.missing.gr = unique(tr.missing.gr)
        ## save exon ids and number to add to meta data
        exon_ids = tr.missing.gr$exon_id
        exon_number = tr.missing.gr$exon_id
        mcols(tr.missing.gr) = meta.add.dt
        tr.missing.gr$exon_ids = exon_ids
        tr.missing.gr$exon_number = exon_number
        ##
        strand_gene = md.sub.tr@strand@values
        tr.missing.dt = as.data.table(tr.missing.gr)
        tr.missing.dt[, strand := strand_gene]
        tr.missing.gr = tr.missing.dt %>% GRanges(., seqlengths = hg_seqlengths())
        ##
        md.sub.tr2 = c(md.sub.tr,tr.missing.gr)
        md.sub.tr.dt = as.data.table(md.sub.tr2)
        ## md.sub.tr.dt[is.na(type), type := "missing"]
        md.sub.tr.dt[is.na(coding_type_simple), coding_type_simple := "missing"]
        pos.dt = md.sub.tr.dt[strand == "+" | strand == "*",][order(seqnames,start,end),]
        neg.dt = md.sub.tr.dt[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
        result.dt = rbind(pos.dt,neg.dt)
        return(result.dt)
      }, mc.cores = cores)
      md.all.dt2 = rbindlist(md.all.dt2, fill = TRUE)
      ## add intron labels
      message("adding labels for retained introns")
      md.all.dt2[, id_coord := 1:.N]
      sub.intron.dt = md.all.dt2[unlist(coding_type_simple) == "intron",]
      ## sub.intron.dt[qname %in% c("transcript/583727","transcript/598795"),]
      ## add labels for each intron by each transcript
      tr.dt = potential.gtf.dt[transcript_id %in% unique(sub.intron.dt$transcript_id),]
      tr.gr = tr.dt[type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      ## just make the data.table with transcript introns and exons before adding it to the reads
      tr.intron.lst = mclapply(unique(tr.dt$transcript_id), function(x) {
        full.tr.gr = tr.dt[transcript_id == x & type == "transcript",] %>% GRanges(., seqlengths = hg_seqlengths())
        parts.tr.gr = tr.dt[transcript_id == x & type != "transcript",] %>% GRanges(., seqlengths = hg_seqlengths())
        tile.gr = gr.tile(full.tr.gr, 1)
        tile.gr = gr.val(tile.gr,parts.tr.gr, c("exon_id","exon_number"))
        tile.gr = gr.reduce(tile.gr, by = "exon_number",ignore.strand = FALSE) %>% sort
        tile.dt = as.data.table(tile.gr)
        if(tile.dt$strand[1] == "-") {
          tile.dt = tile.dt[order(-start),]
          tile.dt[exon_number == "", intron_number := 1:.N]
        } else if (tile.dt$strand[1] == "+") {
          tile.dt = tile.dt[order(start),]
          tile.dt[exon_number == "", intron_number := 1:.N]
        }
        tile.dt[, transcript_id := x]
        return(tile.dt)
      }, mc.cores = cores)
      tr.intron.dt = rbindlist(tr.intron.lst)
      if(nrow(tr.intron.dt) > 0) {
        tr.intron.gr = GRanges(tr.intron.dt, seqlengths = hg_seqlengths())
      }
      ## add transcript labels to the reads
      ##sub.intron.lst = mclapply(unique(md.all.dt2$qname)[seq(4,400,16)], function(x) {
      message("adding transcript labels to reads")
      sub.intron.lst = mclapply(unique(md.all.dt2$qname), function(x) {
        sub.gr = md.all.dt2[qname == x,] %>% GRanges(., seqlengths = hg_seqlengths())
        if(nrow(tr.intron.dt[transcript_id == sub.gr$transcript_id[1],])>0) {
          sub.tr.gr = tr.intron.dt[transcript_id == sub.gr$transcript_id[1],] %>% GRanges(., seqlengths = hg_seqlengths())
          sub.gr = gr.val(sub.gr, sub.tr.gr, "intron_number")
          ## add all unique missing exon and retained intron numbers
          unique.missing.exon.numbers = unique(na.omit(sub.gr$exon_number))
          if(length(unique.missing.exon.numbers) == 0) {
            unique.missing.exon.numbers = NA
          }
          unique.retained.intron.numbers = unique(na.omit(sub.gr$intron_number))
          if(length(unique.retained.intron.numbers) == 0) {
            unique.retained.intron.numbers = NA
          }
          sub.gr$missing_exons = paste0(unique.missing.exon.numbers, collapse = ",")
          sub.gr$retained_introns = paste0(unique.retained.intron.numbers, collapse = ",")
          ## fix labels where extended UTR is called UTR
          whole_transcript.gr = sub.tr.gr %>% gr.reduce
          whole_transcript.gr$range = "in_gene"
          sub.gr = gr.val(sub.gr, whole_transcript.gr, "range")
          sub.dt = as.data.table(sub.gr)
        } else {
          sub.dt = as.data.table(sub.gr)
          sub.dt[, range := "in_gene"]
        }
        sub.dt[, coding_type_simple := as.character(coding_type_simple)]
        sub.dt[range == "" & coding_type_simple == "intron", coding_type_simple := "UTR_long"]
        sub.dt[, range := NULL]
        ## add footprint for whole read (not including "missing") and whole read without strictness of UTR
        whole.fp = GRanges(sub.dt[coding_type_simple != "missing",]) %>% gr.nochr %>% gr.string %>% paste0(.,collapse = ",")
        ## sub.fp = GRanges(sub.dt[coding_type_simple != "missing" & coding_type_simple != "UTR",]) %>% gr.nochr %>% gr.string %>% paste0(.,collapse = ",")
        sub.fp = GRanges(sub.dt[coding_type_simple != "missing" & coding_type_simple != "UTR" & coding_type_simple != "UTR_long",]) %>% gr.nochr %>% gr.string %>% paste0(.,collapse = ",")
        ## whole.fp = GRanges(sub.dt[coding_type_simple != "missing",]) %>% gr.nochr %>% gr.string
        ## sub.fp = GRanges(sub.dt[coding_type_simple != "missing" & coding_type_simple != "UTR",]) %>% gr.nochr %>% gr.string
        sub.dt[, footprint := list(whole.fp)]
        sub.dt[, footprint_minus_UTR := list(sub.fp)]
        return(sub.dt)
      }, mc.cores = cores)
      ## test.dt = rbindlist(sub.intron.lst)
      md.all.dt2 = rbindlist(sub.intron.lst, fill = TRUE)
      pos.dt = md.all.dt2[strand == "+" | strand == "*",][order(seqnames,start,end),]
      neg.dt = md.all.dt2[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
      md.all.dt2 = rbind(pos.dt, neg.dt)
      ##
    } else {
      md.all.dt2 = md.grl2 %>% unlist %>% as.data.table
    }
    message("done annotating missing regions")
    md.gr3 = GRanges(md.all.dt2, seqlengths = hg_seqlengths())
    md.grl3 = rtracklayer::split(md.gr3, f = mcols(md.gr3)["qname"])
  } else {
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    md.grl3 = split(md.gr2, f = mcols(md.gr2)["qname"])
    return(md.grl3)
  }
  ## sort positive strand forward and negative strand reverse
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    md.gw = gW(grl = md.grl3)
    ## add coloring of nodes to match track.gencode
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "UTR_long","del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "red", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
    transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
    md.gw$set(name = transcript_names)
    message("done converting grl to gW object. returning gW")
    return(md.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(md.grl3)
  }
}




get_iso_fusions = function(bam,
                         gr,
                         gtf,
                         collapsed_group = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         collapsed_class = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         type = "gw",
                         annotate_mismatch = TRUE,
                         annotate_mismatch_type = "smart",
                         reann_with = "both", #whether to match using both exons and UTRS or just exon ("exon")
                         reannotate = FALSE,
                         select_type = "bases", #only bases works for annotate_mismatch_type = "smart"
                         add_gencode_name = TRUE, #add the gencode name to the transcript for plotting
                         reverse_overlap = FALSE, #old method - need to remove
                         consider = "only_exon_alignment", #only applies to select_type = "fraction"
                         read_over_potential = 0.5, #only applies to annotate_mismatch_type = "percent"
                         potential_over_read = 0.8, #only applies to annotate_mismatch_type = "percent"
                         minimize = "tr_bases_missing", #only applies to annotate_mismatch_type = smart, typically works better then sum
                         annotate_missing = TRUE, ## match walk with transcript and annotate missing regions as missing
                         walk_per_alignment = FALSE, ## whether to make the qname unique for each alignment
                         cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  if(walk_per_alignment) {
    md.dt[, qname := paste0(qname,"_",1:.N), by = "qname"]
  }
  ########################
  ##attempt to order by the number of soft clipped bases at the beginning
  ##md.dt = md.dt[qname == "transcript/482749",]
  ## get number of soft clipped bases at beginning and order in that manner
  md.dt[, soft_clipped_bases := tstrsplit(cigar, "S", keep = 1, fixed = TRUE)]
  md.dt[, order_softclipped :=  as.integer(soft_clipped_bases)]## anything that is na is the first read
  md.dt[is.na(order_softclipped), order_softclipped := 0]
  md.dt = md.dt[order(order_softclipped),]
  md.dt[, n_order_softclipped := 1:.N, by = "qname"]
########################  
  ## md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##  md.gr = md.gr %Q% (!is.na(cigar))
  ##need to add n_order_softclipped after parsing cigar
  order.dt = as.data.table(md.gr)[,.(qname,n_order_softclipped,order_softclipped)]
  order.dt[, ID := 1:.N]
  order.dt[, qname := NULL]
  ##
  message("Splicing cigar string")
  ## md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = TRUE)
  message("Done Splicing cigar string")
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
  ## fix to add IDs for ordering
  md.grl.dt = grl.unlist(md.grl) %>% as.data.table
  names(md.grl.dt) = gsub("grl.iix", "order_id", names(md.grl.dt))
  md.grl.dt = merge.data.table(md.grl.dt, order.dt, by.x = "order_id", by.y = "ID", all.x = TRUE)
  md.grl.gr = GRanges(md.grl.dt, seqlengths = hg_seqlengths())
  ## md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["grl.ix"])
  md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["order_id"])
  
  if(!reannotate) {
    ##look for collapsed_group and collapsed_class if missing
    if(is.null(collapsed_group) | is.null(collapsed_class)) {
      if(is.null(collapsed_group) && is.null(collapsed_class)) {
        message("looking for collapsed_group and collapsed class")
        folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
        collapsed_group = paste0(folder_find,"collapsed.group.txt")
        collapsed_class = paste0(folder_find,"collapsed_classification.txt")
        if(all(file.exists(c(collapsed_group,collapsed_class)))) {
          message("Found! collapsed.group.txt and collapsed_classification.txt")
        }
      }
    }
    group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
    class.dt = fread(collapsed_class)
    ## get labels for already identified transcripts
    group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
    group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
    group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
    expanded.dt = copy(group.sub.dt)
    expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
    group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
      md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
    } else {
      md.dt3 = group.sub.dt2
    }
    
  } else if(reannotate) {
    if(length(md.grl) > 0 ) {##necessary for fusions to add this check
      md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type, grl.ix, order_id, order_softclipped,n_order_softclipped)][type != "N",]
    }
  }
  if(annotate_mismatch) {
    if(!reannotate) {
      md.full = md.dt3[structural_category == "full-splice_match",]
      md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
      md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
    } else if(reannotate) {
      md.novel = md.dt3
      md.full = NULL
    }
###################################
    ## first  label ones without a transcript annotation
    md.novel.gr = GRanges(md.novel)
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
    ##add potential genes to md.novel
    if(length(md.novel.gr) > 0) {
      md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
      potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
      potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
      ## annotate to the longest transcript??
      if(annotate_mismatch_type == "longest") {
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "level") {
        ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "match") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        ##potential one
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
        ## get the transcript with the most N by transcript
        potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "percent") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              if(consider == "only_exon_alignment") {
                md.sub.gr = md.sub.gr %&% potential.sub.gr
              }
              md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
            } else if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                md.sub.gr$percent = x$width.ov
              } else {
                md.sub.gr$percent = 0
              }
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
              ## dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "fraction") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr2 := rank(pot_over_tr,-mean_overlap,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        if(select_type == "bases") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        ## need to rewrite to loop through the correct combinations
        if(reverse_overlap) { # this option should be removed- ultimately in the loop above
          message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
          rows.cycle = mean.overlap.dt[rank_tr <= 20,]
          mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
            rows.cycle.sub = rows.cycle[x,]
            tr = rows.cycle.sub$query
            pot.tr = rows.cycle.sub$subject
            md.sub.gr = md.novel.grl[[tr]]          
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
            } else if(select_type == "bases") {
              ## potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = potential.sub.gr
              y = md.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                potential.sub.gr$percent = x$width.ov
              } else {
                potential.sub.gr$percent = 0
              }

            }
            ##try returning the mean
            return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          if(select_type == "fraction") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
            mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
          if(select_type == "bases") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
            mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
        } else {
          mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
        if(is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
        }
        if(!is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
          add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
          md.novel.gr$gene_name = NULL
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
          if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
            warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
          }
        }
        ## now annotation the exons overlaps
        md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
        potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      } else if(annotate_mismatch_type == "smart") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0) {
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                ##new part for annotate_mismatch_type == "smart"
                y_sum_bases = y@ranges@width %>% sum
                x_sum_bases = x@ranges@width %>% sum
                ov_sum_bases = x$width.ov %>% sum
                ## ov_sum_bases = ov$width %>% sum
                potential_leftout = y_sum_bases - ov_sum_bases
                tr_leftout = x_sum_bases - ov_sum_bases
                ##end new
              } else {
                tr_leftout = "ALL"
                potential_leftout = "ALL"
              }
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = tr_leftout, pot_over_tr = potential_leftout)
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "bases") {
          mean.overlap.dt[, mean_overlap := as.numeric(mean_overlap)]
          mean.overlap.dt[, pot_over_tr := as.numeric(pot_over_tr)]
          mean.overlap.dt[, sum_missing := (mean_overlap + pot_over_tr)]
          mean.overlap.dt[, rank_tr := rank(sum_missing,ties.method = "min"), by = "query"]
          mean.overlap.dt = mean.overlap.dt[order(mean_overlap, rank_tr),]
          if(minimize == "sum") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
          }
          if(minimize == "tr_bases_missing") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
          }
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
      if(is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
      }
      if(!is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
        add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
        md.novel.gr$gene_name = NULL
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
        if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
          warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
        }
      }
      ## now annotation the exons overlaps
      md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
      potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      }
    }
    if(exists("md.novel.dt2") & !reannotate) {
    gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
    md.annotated = rbind(md.full,md.sub)
    md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
    md.annotated[, transcript_id := associated_transcript]
    md.novel.dt2[, reannotated := TRUE]
    md.annotated[, reannotated := FALSE]
    md.annotated = rbind(md.novel.dt2,md.annotated)
    } else if (reannotate) {
      md.annotated = md.novel.dt2
      md.annotated[, reannotated := TRUE]
    } else {
      gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique
      md.annotated = rbind(md.full,md.sub)
      md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
      md.annotated[, transcript_id := associated_transcript]
      md.annotated[, reannotated := FALSE]
    }
    md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
    potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
    potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
    potential.gtf.dt = as.data.table(potential.gtf.gr)
    
    add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
    md.annotated.dt2 = md.annotated
    md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
    ## now annotation the exons overlaps
    md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
    potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
    potential.transcripts = potential.gtf.labels$transcript_id %>% unique
    md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
      md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
      potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
      ##breaks isn't working properly-convert these to breakpoints
      potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end,strand)]
      ##test end + 1 for pos
      bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt)))) %>% unique
      bps.dt[,start := end -1]
      bps.dt[start_end == "start", end := end - 1]
      bps.dt[start_end == "start", start := start - 1]
      bps.dt[start_end == "end", end := end + 1]
      bps.dt[start_end == "end", start := start + 1]
      ##test labels
      ## bps.dt[,start := start + 1]
      ## bps.dt[,end := start]
      ##
      bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths()) %>% gr.chr
      
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
      md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
      ## new gr.val to match to specific transcripts
      md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      return(md.novel.dt3)
    })
    md.annotated.dt3 = rbindlist(md.annotated.lst3)  
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
################
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
    md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                                  ifelse("UTR" %in% x, "UTR",
                                                                                  ifelse("exon" %in% x, "exon","intron"))))]
    md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3[type == "X",coding_type_simple := "del"]
    md.all.dt = md.annotated.dt3
    message("finished merging together")
    ## rbind all together
    ## sort based on strand to have walks align
    ## fix with ordering for fusions
    md.all.dt$strand = factor(md.all.dt$strand, levels = c("+","-"))
    ##pos.dt = md.all.dt[strand == "+",][order(grl.ix, strand, order_id),][qname == "transcript/611213",]
    md.all.dt[, order_start := ifelse(strand == "+", start, -start)]
    md.dt4 = md.all.dt[order(order_softclipped, strand, order_start),]
    
    ## md.dt4 = md.all.dt[order(strand, order_start, grl.ix,  order_id),] ##[qname == "transcript/611213",]

    ##### OLD
      ## md.all.dt[order(strand, grl.ix,  order_id),][qname == "transcript/611213",]
    ## pos.dt = md.all.dt[strand == "+",][order(order_id),]
    ## neg.dt = md.all.dt[strand == "-",][order(order_id),]
    ## pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
    ## neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
    ## md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
    ## END OLD
    if(add_gencode_name) {
      gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
      gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
      md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
    }
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    ## convert to grl and gw
    message("converting to grl")
    md.gr2$qr_tr = paste0(gsub("transcript/","",md.gr2$qname),"; ", md.gr2$transcript_name)
    md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
    message("done coverting to grl")
    ## now annotating missing parts
    if(annotate_missing) {
      message("annotating missing regions from transcripts")
      md.all.dt2 = mclapply(1:length(md.grl2), function(x) {
        md.sub.tr = md.grl2[[x]]
        qr_tr1 = md.sub.tr$qr_tr[1]
        qname_tr1 = md.sub.tr$qname[1]
        tr.gr = potential.gtf.dt[transcript_id == unique(md.sub.tr$transcript_id),][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
                                        #tr.missing.gr = tr.gr[tr.gr %outside% md.sub.tr]
                                        #mcols(tr.missing.gr) = NULL
        ##attempt not using outside
        tiled.gr = gr.tile(tr.gr,width = 1)
        tiled.gr$overlap = tiled.gr %O% md.sub.tr
        tr.missing.gr = (tiled.gr %Q% (overlap < 1)) %>% gr.reduce
        if(length(tr.missing.gr) == 0) {
          return(as.data.table(md.sub.tr))
        }
        tr.missing.gr = as.data.table(tr.missing.gr) %>% GRanges(., seqlengths = hg_seqlengths())
        ##end  attempt
        ##add meta data
        meta.add.dt = mcols(md.sub.tr)[1,] %>% as.data.table
        ##meta.add.dt[, c("qid", "type","rid","riid","fid","col",) := NULL]
        if("isoform" %in% names(meta.add.dt)) {
          cols_keep = c("transcript_id", "qname", "isoform", "structural_category", "subcategory", "associated_transcript", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        } else {
          cols_keep = c("transcript_id", "qname", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        }
        ## meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
        meta.add.dt = meta.add.dt[,..cols_keep]
        tr.missing.gr = unique(tr.missing.gr)
        mcols(tr.missing.gr) = meta.add.dt
        ##
        strand_gene = md.sub.tr@strand@values
        tr.missing.gr@strand@values = strand_gene
        ##
        md.sub.tr2 = c(md.sub.tr,tr.missing.gr)
        md.sub.tr.dt = as.data.table(md.sub.tr2)
        md.sub.tr.dt[is.na(type), type := "missing"]
        md.sub.tr.dt[is.na(coding_type_simple), coding_type_simple := "missing"]
        pos.dt = md.sub.tr.dt[strand == "+" | strand == "*",][order(seqnames,start,end),]
        neg.dt = md.sub.tr.dt[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
        return(rbind(pos.dt,neg.dt))
      }, mc.cores = cores)
      md.all.dt2 = rbindlist(md.all.dt2, fill = TRUE)
    } else {
      md.all.dt2 = md.grl2 %>% unlist %>% as.data.table
    }
    message("done annotating missing regions")
    md.gr3 = GRanges(md.all.dt2, seqlengths = hg_seqlengths())
    md.grl3 = rtracklayer::split(md.gr3, f = mcols(md.gr3)["qname"])
  } else {
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    md.grl3 = split(md.gr2, f = mcols(md.gr2)["qname"])
    return(md.grl3)
  }
  ## sort positive strand forward and negative strand reverse
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    md.gw = gW(grl = md.grl3)
    ## add coloring of nodes to match track.gencode
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
    transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
    md.gw$set(name = transcript_names)
    message("done converting grl to gW object. returning gW")
    return(md.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(md.grl3)
  }
}


##function for reading in pacbio's *.breakpoints.groups.bed
read_fusion_groups = function(file) {
  fus.group.dt = read.delim(file,header=F,comment.char='#') %>% as.data.table
  names(fus.group.dt) = c("chr1", "start1", "end1", "chr2", "start2", "end2", "id", "score", "strand1", "strand2", "info", "extra")
  fus.group.dt[,qname := tstrsplit(extra, ";", fixed = TRUE, keep = 1)]
  fus.group.dt[,qname := gsub("RN=","",qname)]
  fus.group.dt[, qname_list := strsplit(qname, ",")] #split transcripts into unique rows
  fus.group.dt[, row_number := 1:.N]
                                        #fus.group.dt[, id_number := .N, by = "row_number"]
  expanded.dt = copy(fus.group.dt)
  expanded.dt2 = expanded.dt[, .(qname = unlist(qname_list)), by = row_number]
  fus.group.dt2 = merge.data.table(expanded.dt[, !"qname"], expanded.dt2, by = "row_number")
  fus.group.dt2[,CB := tstrsplit(extra, ";", fixed = TRUE, keep = 3)]
  fus.group.dt2[, GN := tstrsplit(info, ";", fixed = TRUE, keep = 3)]
  fus.group.dt2[, GN := gsub("GN=","",GN)]
  fus.group.dt2[,AC := tstrsplit(info, ";", fixed = TRUE, keep = 13)]
  return(fus.group.dt2)
}

## function to reduce the walks so the different colors don't plot on different lines
reduce_gw = function(gw, return = "walk", type = "isoseq", cores = 1) {
  iso.grl = gw$grl
  reduced.grl.dt = mclapply(1:length(iso.grl), function(x) {
    sub.gr = iso.grl[[x]]
    sub.gr$original_node_order = as.numeric(1:length(sub.gr))
    ##attempt at padding by one to get the original order
    ## sub.pad.gr = sub.gr #
    ## sub.pad.gr$pad_order = 1:length(sub.pad.gr)
    ## sub.pad.gr = gr.reduce(sub.pad.gr, ignore.strand = FALSE, "pad_order")
    ##
    sub.gr2 = gr.reduce(sub.gr,ignore.strand = FALSE)
    sub.gr2$original_walk_id = x
    ## have to consider the original order
    sub.gr2$original_node_order = NULL
    sub.gr2 = gr.val(sub.gr2,sub.gr, "original_node_order", mean = FALSE, ignore.strand = FALSE, FUN = min, na.rm = TRUE) #
    sub.dt = as.data.table(sub.gr2)
    sub.dt[, original_node_order := as.integer(original_node_order)]
    sub.dt = sub.dt[order(original_node_order),]
    return(sub.dt)
  }, mc.cores = cores) %>% rbindlist
  final.gr = GRanges(reduced.grl.dt, seqlengths = hg_seqlengths())
  final.grl = rtracklayer::split(final.gr, f = mcols(final.gr)["original_walk_id"])
  final.gw = gW(grl = final.grl)
  ## add original names and color used in isoseq reading in
  original.names = gw$dt$name
  final.gw$set(name = original.names)
  if(type == "isoseq") {
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      final.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
  } else if (type == "fusions") {
  }
  return(final.gw)
}


get_ordered_fusions = function(bam,
                         gr,
                         gtf,
                         collapsed_group = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         collapsed_class = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         type = "gw",
                         annotate_mismatch = TRUE,
                         annotate_mismatch_type = "smart",
                         reann_with = "both", #whether to match using both exons and UTRS or just exon ("exon")
                         reannotate = FALSE,
                         select_type = "bases", #only bases works for annotate_mismatch_type = "smart"
                         add_gencode_name = TRUE, #add the gencode name to the transcript for plotting
                         reverse_overlap = FALSE, #old method - need to remove
                         consider = "only_exon_alignment", #only applies to select_type = "fraction"
                         read_over_potential = 0.5, #only applies to annotate_mismatch_type = "percent"
                         potential_over_read = 0.8, #only applies to annotate_mismatch_type = "percent"
                         minimize = "tr_bases_missing", #only applies to annotate_mismatch_type = smart, typically works better then sum
                         annotate_missing = TRUE, ## match walk with transcript and annotate missing regions as missing
                         walk_per_alignment = FALSE, ## whether to make the qname unique for each alignment
                         breakpoints = NULL,
                         pad = 10, #pad used for the breakpoints to identify the order of alignments
                         cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  if(length(md.gr.bam) == 0) {
    warning("no reads were read in this region")
    return(NULL)
  }
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  if(walk_per_alignment) {
    md.dt[, qname := paste0(qname,"_",1:.N), by = "qname"]
  }
########################
  md.dt[, soft_clipped_bases := tstrsplit(cigar, "S", keep = 1, fixed = TRUE)]
  md.dt[, order_softclipped :=  as.integer(soft_clipped_bases)]## anything that is na is the first read
  md.dt[is.na(order_softclipped), order_softclipped := 0]
  md.dt = md.dt[order(order_softclipped),]
  md.dt[, n_order_softclipped := 1:.N, by = "qname"]

  ## look for breakpoints if empty
  if(is.null(breakpoints)) {
    ##look for it
    folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-2)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
    if(file.exists(paste0(folder_find,"fusion.breakpoints.groups.bed"))) {
      fus.groups.dt = read_fusion_groups(paste0(folder_find,"fusion.breakpoints.groups.bed"))
      } else {
        warning("did not find fusion.breakpoints.groups.bed. Will order based on least soft clipped base pairs")
        ##attempt to order by the number of soft clipped bases at the beginning
        ##md.dt = md.dt[qname == "transcript/482749",]
        ## get number of soft clipped bases at beginning and order in that manner
      }
  } else {
    fus.groups.dt = read_fusion_groups(breakpoints)
  }
  
########################  
  ## md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##  md.gr = md.gr %Q% (!is.na(cigar))
  ##need to add n_order_softclipped after parsing cigar
  order.dt = as.data.table(md.gr)[,.(qname,n_order_softclipped,order_softclipped)]
  order.dt[, ID := 1:.N]
  order.dt[, qname := NULL]
  ##
  message("Splicing cigar string")
  ## md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = TRUE)
  message("Done Splicing cigar string")
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
  ## fix to add IDs for ordering
  md.grl.dt = grl.unlist(md.grl) %>% as.data.table
  names(md.grl.dt) = gsub("grl.iix", "order_id", names(md.grl.dt))
  md.grl.dt = merge.data.table(md.grl.dt, order.dt, by.x = "order_id", by.y = "ID", all.x = TRUE)
  md.grl.gr = GRanges(md.grl.dt, seqlengths = hg_seqlengths())
  ## md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["grl.ix"])
  md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["order_id"])
  
  if(!reannotate) {
    ##look for collapsed_group and collapsed_class if missing
    if(is.null(collapsed_group) | is.null(collapsed_class)) {
      if(is.null(collapsed_group) && is.null(collapsed_class)) {
        message("looking for collapsed_group and collapsed class")
        folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
        collapsed_group = paste0(folder_find,"collapsed.group.txt")
        collapsed_class = paste0(folder_find,"collapsed_classification.txt")
        if(all(file.exists(c(collapsed_group,collapsed_class)))) {
          message("Found! collapsed.group.txt and collapsed_classification.txt")
        }
      }
    }
    group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
    class.dt = fread(collapsed_class)
    ## get labels for already identified transcripts
    group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
    group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
    group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
    expanded.dt = copy(group.sub.dt)
    expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
    group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
      md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
    } else {
      md.dt3 = group.sub.dt2
    }
    
  } else if(reannotate) {
    if(length(md.grl) > 0 ) {##necessary for fusions to add this check
      md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type, grl.ix, order_id, order_softclipped,n_order_softclipped)][type != "N",]
    }
  }
  if(annotate_mismatch) {
    if(!reannotate) {
      md.full = md.dt3[structural_category == "full-splice_match",]
      md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
      md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
    } else if(reannotate) {
      md.novel = md.dt3
      md.full = NULL
    }
###################################
    ## first  label ones without a transcript annotation
    md.novel.gr = GRanges(md.novel)
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
    ##add potential genes to md.novel
    if(length(md.novel.gr) > 0) {
      md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
      potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
      potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
      ## annotate to the longest transcript??
      if(annotate_mismatch_type == "longest") {
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "level") {
        ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "match") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        ##potential one
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
        ## get the transcript with the most N by transcript
        potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "percent") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              if(consider == "only_exon_alignment") {
                md.sub.gr = md.sub.gr %&% potential.sub.gr
              }
              md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
            } else if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                md.sub.gr$percent = x$width.ov
              } else {
                md.sub.gr$percent = 0
              }
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
              ## dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "fraction") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr2 := rank(pot_over_tr,-mean_overlap,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        if(select_type == "bases") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        ## need to rewrite to loop through the correct combinations
        if(reverse_overlap) { # this option should be removed- ultimately in the loop above
          message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
          rows.cycle = mean.overlap.dt[rank_tr <= 20,]
          mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
            rows.cycle.sub = rows.cycle[x,]
            tr = rows.cycle.sub$query
            pot.tr = rows.cycle.sub$subject
            md.sub.gr = md.novel.grl[[tr]]          
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
            } else if(select_type == "bases") {
              ## potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = potential.sub.gr
              y = md.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                potential.sub.gr$percent = x$width.ov
              } else {
                potential.sub.gr$percent = 0
              }

            }
            ##try returning the mean
            return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          if(select_type == "fraction") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
            mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
          if(select_type == "bases") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
            mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
        } else {
          mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
        if(is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
        }
        if(!is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
          add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
          md.novel.gr$gene_name = NULL
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
          if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
            warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
          }
        }
        ## now annotation the exons overlaps
        md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
        potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      } else if(annotate_mismatch_type == "smart") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0) {
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                ##new part for annotate_mismatch_type == "smart"
                y_sum_bases = y@ranges@width %>% sum
                x_sum_bases = x@ranges@width %>% sum
                ov_sum_bases = x$width.ov %>% sum
                ## ov_sum_bases = ov$width %>% sum
                potential_leftout = y_sum_bases - ov_sum_bases
                tr_leftout = x_sum_bases - ov_sum_bases
                ##end new
              } else {
                tr_leftout = "ALL"
                potential_leftout = "ALL"
              }
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = tr_leftout, pot_over_tr = potential_leftout)
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "bases") {
          mean.overlap.dt[, mean_overlap := as.numeric(mean_overlap)]
          mean.overlap.dt[, pot_over_tr := as.numeric(pot_over_tr)]
          mean.overlap.dt[, sum_missing := (mean_overlap + pot_over_tr)]
          mean.overlap.dt[, rank_tr := rank(sum_missing,ties.method = "min"), by = "query"]
          mean.overlap.dt = mean.overlap.dt[order(mean_overlap, rank_tr),]
          if(minimize == "sum") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
          }
          if(minimize == "tr_bases_missing") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
          }
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
      if(is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
      }
      if(!is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
        add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
        md.novel.gr$gene_name = NULL
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
        if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
          warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
        }
      }
      ## now annotation the exons overlaps
      md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
      potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      }
    }
    if(exists("md.novel.dt2") & !reannotate) {
    gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
    md.annotated = rbind(md.full,md.sub)
    md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
    md.annotated[, transcript_id := associated_transcript]
    md.novel.dt2[, reannotated := TRUE]
    md.annotated[, reannotated := FALSE]
    md.annotated = rbind(md.novel.dt2,md.annotated)
    } else if (reannotate) {
      md.annotated = md.novel.dt2
      md.annotated[, reannotated := TRUE]
    } else {
      gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique
      md.annotated = rbind(md.full,md.sub)
      md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
      md.annotated[, transcript_id := associated_transcript]
      md.annotated[, reannotated := FALSE]
    }
    md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
    potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
    potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
    potential.gtf.dt = as.data.table(potential.gtf.gr)
    
    add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
    md.annotated.dt2 = md.annotated
    md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
    ## now annotation the exons overlaps
    md_potential_tr.gr = GRanges(md.annotated.dt2, seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    ## md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
    potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
    potential.transcripts = potential.gtf.labels$transcript_id %>% unique
    md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
      md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
      potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
      ##breaks isn't working properly-convert these to breakpoints
      potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end,strand)]
      ##test end + 1 for pos
      bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt)))) %>% unique
      bps.dt[,start := end -1]
      bps.dt[start_end == "start", end := end - 1]
      bps.dt[start_end == "start", start := start - 1]
      bps.dt[start_end == "end", end := end + 1]
      bps.dt[start_end == "end", start := start + 1]
      ##test labels
      ## bps.dt[,start := start + 1]
      ## bps.dt[,end := start]
      ##
      bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths()) %>% gr.chr
      
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
      md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
      ## new gr.val to match to specific transcripts
      md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      return(md.novel.dt3)
    })
    md.annotated.dt3 = rbindlist(md.annotated.lst3)
    ##add back ones that were not assigned a transcript
    md.annotated.dt3 = rbind(md.annotated.dt3, md.annotated.dt2[transcript_id == "multiple_potential_genes",], fill = TRUE)
    ##
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
################
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
    md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                                  ifelse("UTR" %in% x, "UTR",
                                                                                  ifelse("exon" %in% x, "exon","intron"))))]
    md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3[type == "X",coding_type_simple := "del"]
    md.all.dt = md.annotated.dt3
    md.all.dt[,qname_simple := tstrsplit(qname, "_", keep = 1)]
    message("finished merging together")
    ## rbind all together
    ## new ordering using the fusion groups
    if(exists("fus.groups.dt")) {
      ## fus.groups.dt[, extra := as.list(extra)]
      ## test.fus.dt = fus.groups.dt[qname == "transcript/482749",.(row_number, chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, qname_list, qname, CB, AC)]
      ## test.md.gr = GRanges(md.all.dt[grepl("transcript/482749", qname),], seqlengths = hg_seqlengths())
      
      ##using the inherent position 1 and 2
      ## bp1.gr = test.fus.dt[,.(chr1,start1,end1)] %>% setnames(.,c("seqnames","start","end")) %>% GRanges(., seqlengths = hg_seqlengths())
      ## bp2.gr = test.fus.dt[,.(chr2,start2,end2)] %>% setnames(.,c("seqnames","start","end")) %>% GRanges(., seqlengths = hg_seqlengths())
      ## parsing the AC with the more than two if there are
      qnames.lst = fus.groups.dt$qname %>% tstrsplit(., "_", keep = 1) %>% unique %>% unlist
      qnames.lst = qnames.lst[qnames.lst %in% md.all.dt$qname_simple]
      md.lst = mclapply(qnames.lst, function(x) {
        sub.fus.dt = fus.groups.dt[qname == x,.(row_number, chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, qname_list, qname, CB, AC)]
        sub.md.gr = GRanges(md.all.dt[qname_simple == x,], seqlengths = hg_seqlengths())

        bps.lst = sub.fus.dt$AC %>% gsub("AC=","",.) %>% tstrsplit(.,"/") %>% unlist
        qnames.dt = lapply(1:length(bps.lst), function(x) {
          pad1 = pad ## just returning this at this point for debugging purposes
          bp.gr = parse.gr(bps.lst[x]) + pad
          int.gr = sub.md.gr %&% (bp.gr)
          ## if no intersection- try adding a pad a few times
          if(length(int.gr) == 0) {
            pad2 = pad * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad2
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          if(length(int.gr) == 0) {
            pad3 = pad2 * 10
            bp.gr2 = parse.gr(bps.lst[x]) + pad3
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          if(length(int.gr) == 0) {
            pad4 = pad3 * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad4
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          qname1 = int.gr$qname %>% unique
          qnames.dt = data.table(qname = qname1, order = x)
          return(qnames.dt)
        }) %>% rbindlist(.)
        mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
        qnames.dt = qnames.dt[, .SD[which.min(order)], by = "qname"]
        if(nrow(qnames.dt) < length(bps.lst)) {
          ##get which qname is missing
          qnames.lst = sub.md.gr$qname %>% unique
          qname.missing.lst = qnames.lst[!(qnames.lst %in% qnames.dt$qname)]
          qname.missing.dt = data.table(qname = qname.missing.lst, order = 10)
          qnames.dt = rbind(qnames.dt, qname.missing.dt)
        }
        md.sub.dt = as.data.table(sub.md.gr)
        
        ## qname_1st = (test.md.gr %&% (bp1.gr+20))$qname %>% unique ## get the unique qname to paste in order
        ## qname_2nd = (test.md.gr %&% (bp2.gr+20))$qname %>% unique ## get the unique qname to paste in order
        
        md.dt4 = merge.data.table(md.sub.dt, qnames.dt, by = "qname", all.x = TRUE)
                                        #      md.dt4[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt4[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt4 = md.dt4[, c("start", "end") := .(ifelse(end < start, end, start), ifelse(end < start, start, end))]
        md.dt4 = md.dt4[order(order, order_start),]
        md.dt4[, qname := x]
        md.dt4[, final_order := 1:.N]
        return(md.dt4)
      }, mc.cores = cores)
      md.dt4 = rbindlist(md.lst)
      
      ## md.dt4 = md.dt4[order(order, order_id, order_start),]
      ## order grl.ix after ordering the qnames
    } else {
      ## sort based on strand to have walks align
      ## fix with ordering for fusions
      md.all.dt$strand = factor(md.all.dt$strand, levels = c("+","-"))
      ##pos.dt = md.all.dt[strand == "+",][order(grl.ix, strand, order_id),][qname == "transcript/611213",]
      md.all.dt[, order_start := ifelse(strand == "+", start, -start)]
      md.dt4 = md.all.dt[order(order_softclipped, strand, order_start),]
      ## md.dt4 = md.all.dt[order(strand, order_start, grl.ix,  order_id),] ##[qname == "transcript/611213",]

##### OLD
      ## md.all.dt[order(strand, grl.ix,  order_id),][qname == "transcript/611213",]
      ## pos.dt = md.all.dt[strand == "+",][order(order_id),]
      ## neg.dt = md.all.dt[strand == "-",][order(order_id),]
      ## pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
      ## neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
      ## md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
      ## END OLD
    }
    if(add_gencode_name) {
      gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
      gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
      md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
    }
    md.dt4 = md.dt4[order(qname,final_order),]
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    ## convert to grl and gw
    message("converting to grl")
    md.gr2$qr_tr = paste0(gsub("transcript/","",md.gr2$qname),"; ", md.gr2$transcript_name)
    md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
    message("done coverting to grl")
    ## now annotating missing parts
    if(annotate_missing) {
      message("annotating missing regions from transcripts")
      md.all.dt2 = mclapply(1:length(md.grl2), function(x) {
        md.sub.tr = md.grl2[[x]]
        qr_tr1 = md.sub.tr$qr_tr[1]
        qname_tr1 = md.sub.tr$qname[1]
        tr.gr = potential.gtf.dt[transcript_id == unique(md.sub.tr$transcript_id),][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
                                        #tr.missing.gr = tr.gr[tr.gr %outside% md.sub.tr]
                                        #mcols(tr.missing.gr) = NULL
        ##attempt not using outside
        tiled.gr = gr.tile(tr.gr,width = 1)
        tiled.gr$overlap = tiled.gr %O% md.sub.tr
        tr.missing.gr = (tiled.gr %Q% (overlap < 1)) %>% gr.reduce
        if(length(tr.missing.gr) == 0) {
          return(as.data.table(md.sub.tr))
        }
        tr.missing.gr = as.data.table(tr.missing.gr) %>% GRanges(., seqlengths = hg_seqlengths())
        ##end  attempt
        ##add meta data
        meta.add.dt = mcols(md.sub.tr)[1,] %>% as.data.table
        ##meta.add.dt[, c("qid", "type","rid","riid","fid","col",) := NULL]
        if("isoform" %in% names(meta.add.dt)) {
          cols_keep = c("transcript_id", "qname", "isoform", "structural_category", "subcategory", "associated_transcript", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        } else {
          cols_keep = c("transcript_id", "qname", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        }
        ## meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
        meta.add.dt = meta.add.dt[,..cols_keep]
        tr.missing.gr = unique(tr.missing.gr)
        mcols(tr.missing.gr) = meta.add.dt
        ##
        strand_gene = md.sub.tr@strand@values
        tr.missing.gr@strand@values = strand_gene
        ##
        md.sub.tr2 = c(md.sub.tr,tr.missing.gr)
        md.sub.tr.dt = as.data.table(md.sub.tr2)
        md.sub.tr.dt[is.na(type), type := "missing"]
        md.sub.tr.dt[is.na(coding_type_simple), coding_type_simple := "missing"]
        pos.dt = md.sub.tr.dt[strand == "+" | strand == "*",][order(seqnames,start,end),]
        neg.dt = md.sub.tr.dt[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
        return(rbind(pos.dt,neg.dt))
      }, mc.cores = cores)
      md.all.dt2 = rbindlist(md.all.dt2, fill = TRUE)
    } else {
      md.all.dt2 = md.grl2 %>% unlist %>% as.data.table
      md.all.dt2 = md.all.dt2[order(qname,final_order),]
    }
    message("done annotating missing regions")
    md.gr3 = GRanges(md.all.dt2, seqlengths = hg_seqlengths())
    md.grl3 = rtracklayer::split(md.gr3, f = mcols(md.gr3)["qname"])
  } else {
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    md.grl3 = split(md.gr2, f = mcols(md.gr2)["qname"])
    return(md.grl3)
  }
  ## sort positive strand forward and negative strand reverse
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    md.gw = gW(grl = md.grl3)
    ## add coloring of nodes to match track.gencode
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
    transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
    md.gw$set(name = transcript_names)
    message("done converting grl to gW object. returning gW")
    return(md.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(md.grl3)
  }
}


############################################################################################################################################################################################################################################################################
get_ordered_fusions2 = function(bam,
                         gr,
                         gtf,
                         collapsed_group = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         collapsed_class = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
                         type = "gw",
                         annotate_mismatch = FALSE,
                         annotate_mismatch_type = "none",
                         reann_with = "both", #whether to match using both exons and UTRS or just exon ("exon")
                         reannotate = FALSE,
                         select_type = "bases", #only bases works for annotate_mismatch_type = "smart"
                         add_gencode_name = TRUE, #add the gencode name to the transcript for plotting
                         reverse_overlap = FALSE, #old method - need to remove
                         consider = "only_exon_alignment", #only applies to select_type = "fraction"
                         read_over_potential = 0.5, #only applies to annotate_mismatch_type = "percent"
                         potential_over_read = 0.8, #only applies to annotate_mismatch_type = "percent"
                         minimize = "tr_bases_missing", #only applies to annotate_mismatch_type = smart, typically works better then sum
                         annotate_missing = FALSE, ## match walk with transcript and annotate missing regions as missing
                         walk_per_alignment = TRUE, ## whether to make the qname unique for each alignment
                         breakpoints = NULL,
                         pad = 10, #pad used for the breakpoints to identify the order of alignments
                         cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  if(length(md.gr.bam) == 0) {
    warning("no reads were read in this region")
    return(NULL)
  }
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  if(walk_per_alignment) {
    md.dt[, qname := paste0(qname,"_",1:.N), by = "qname"]
  }
########################
  md.dt[, soft_clipped_bases := tstrsplit(cigar, "S", keep = 1, fixed = TRUE)]
  md.dt[, order_softclipped :=  as.integer(soft_clipped_bases)]## anything that is na is the first read
  md.dt[is.na(order_softclipped), order_softclipped := 0]
  md.dt = md.dt[order(order_softclipped),]
  md.dt[, n_order_softclipped := 1:.N, by = "qname"]

  ## look for breakpoints if empty
  if(is.null(breakpoints)) {
    ##look for it
    folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-2)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
    if(file.exists(paste0(folder_find,"fusion.breakpoints.groups.bed"))) {
      fus.groups.dt = read_fusion_groups(paste0(folder_find,"fusion.breakpoints.groups.bed"))
      } else {
        warning("did not find fusion.breakpoints.groups.bed. Will order based on least soft clipped base pairs")
        ##attempt to order by the number of soft clipped bases at the beginning
        ##md.dt = md.dt[qname == "transcript/482749",]
        ## get number of soft clipped bases at beginning and order in that manner
      }
  } else {
    fus.groups.dt = read_fusion_groups(breakpoints)
  }
  
########################  
  ## md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##  md.gr = md.gr %Q% (!is.na(cigar))
  ##need to add n_order_softclipped after parsing cigar
  order.dt = as.data.table(md.gr)[,.(qname,n_order_softclipped,order_softclipped)]
  order.dt[, ID := 1:.N]
  order.dt[, qname := NULL]
  ##
  message("Splicing cigar string")
  ## md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = TRUE)
  message("Done Splicing cigar string")
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
  ## fix to add IDs for ordering
  md.grl.dt = grl.unlist(md.grl) %>% as.data.table
  names(md.grl.dt) = gsub("grl.iix", "order_id", names(md.grl.dt))
  md.grl.dt = merge.data.table(md.grl.dt, order.dt, by.x = "order_id", by.y = "ID", all.x = TRUE)
  md.grl.gr = GRanges(md.grl.dt, seqlengths = hg_seqlengths())
  ## md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["grl.ix"])
  md.grl = rtracklayer::split(md.grl.gr, f = mcols(md.grl.gr)["order_id"])
  
  if(!reannotate & annotate_mismatch_type != "none") {
    ##look for collapsed_group and collapsed_class if missing
    if(is.null(collapsed_group) | is.null(collapsed_class)) {
      if(is.null(collapsed_group) && is.null(collapsed_class)) {
        message("looking for collapsed_group and collapsed class")
        folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
        collapsed_group = paste0(folder_find,"collapsed.group.txt")
        collapsed_class = paste0(folder_find,"collapsed_classification.txt")
        if(all(file.exists(c(collapsed_group,collapsed_class)))) {
          message("Found! collapsed.group.txt and collapsed_classification.txt")
        }
      }
    }
    group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
    class.dt = fread(collapsed_class)
    ## get labels for already identified transcripts
    group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
    group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
    group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
    expanded.dt = copy(group.sub.dt)
    expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
    group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
    ##add transcript labels to bam reads
    md.dt2 = unlist(md.grl) %>% as.data.table()
    if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
      md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
    } else {
      md.dt3 = group.sub.dt2
    }
    
  } else if(reannotate & annotate_mismatch_type != "none") {
    if(length(md.grl) > 0 ) {##necessary for fusions to add this check
      md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type, grl.ix, order_id, order_softclipped,n_order_softclipped)][type != "N",]
    }
  } else if (annotate_mismatch_type == "none") {
    ## md.dt4 = copy(md.dt)
    ## md.dt4[, c("seq","qual") := NULL]
    md.dt4 = grl.unlist(md.grl) %>% as.data.table
  }
  if(annotate_mismatch) {
    if(!reannotate) {
      md.full = md.dt3[structural_category == "full-splice_match",]
      md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
      md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
    } else if(reannotate) {
      md.novel = md.dt3
      md.full = NULL
    }
###################################
    ## first  label ones without a transcript annotation
    md.novel.gr = GRanges(md.novel)
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
    ##add potential genes to md.novel
    if(length(md.novel.gr) > 0) {
      md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
      potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
      potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
      potential.gtf.dt = as.data.table(potential.gtf.gr)
      ## annotate to the longest transcript??
      if(annotate_mismatch_type == "longest") {
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "level") {
        ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
        potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "match") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        ##potential one
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
        ## get the transcript with the most N by transcript
        potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
      } else if(annotate_mismatch_type == "percent") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              if(consider == "only_exon_alignment") {
                md.sub.gr = md.sub.gr %&% potential.sub.gr
              }
              md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
            } else if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                md.sub.gr$percent = x$width.ov
              } else {
                md.sub.gr$percent = 0
              }
              ## test reverse
              md.sub.gr2 = md.novel.grl[[tr]]
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
              ##
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
              ## dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "fraction") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr2 := rank(pot_over_tr,-mean_overlap,ties.method = "min"), by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        if(select_type == "bases") {
          ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
          mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
          mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
          ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
          ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
        }
        ## need to rewrite to loop through the correct combinations
        if(reverse_overlap) { # this option should be removed- ultimately in the loop above
          message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
          rows.cycle = mean.overlap.dt[rank_tr <= 20,]
          mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
            rows.cycle.sub = rows.cycle[x,]
            tr = rows.cycle.sub$query
            pot.tr = rows.cycle.sub$subject
            md.sub.gr = md.novel.grl[[tr]]          
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "fraction") {
              potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
            } else if(select_type == "bases") {
              ## potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = potential.sub.gr
              y = md.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0){
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                potential.sub.gr$percent = x$width.ov
              } else {
                potential.sub.gr$percent = 0
              }

            }
            ##try returning the mean
            return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          if(select_type == "fraction") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
            mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
          if(select_type == "bases") {
            mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
            mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
            names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
            mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
            mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          }
        } else {
          mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
        if(is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
        }
        if(!is.null(potential_transcripts_merge)) {
          add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
          add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
          md.novel.gr$gene_name = NULL
          md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
          md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
          if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
            warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
          }
        }
        ## now annotation the exons overlaps
        md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
        potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
        potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      } else if(annotate_mismatch_type == "smart") {
        gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
        gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
        potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
        ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
        potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
        potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            ## annotate with utr and exons or just exons
            if(reann_with == "both") {
              potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            } else if (reann_with == "exon") {
              potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
            }
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            if(select_type == "bases") {
              ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
              ## %o% is not working hmmm
              x = md.sub.gr
              y = potential.sub.gr
              ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
              if (nrow(ov)>0) {
                ov = ov[ , sum(width), keyby = query.id]
                x$width.ov = 0
                x$width.ov[ov$query.id] = ov$V1
                ##new part for annotate_mismatch_type == "smart"
                y_sum_bases = y@ranges@width %>% sum
                x_sum_bases = x@ranges@width %>% sum
                ov_sum_bases = x$width.ov %>% sum
                ## ov_sum_bases = ov$width %>% sum
                potential_leftout = y_sum_bases - ov_sum_bases
                tr_leftout = x_sum_bases - ov_sum_bases
                ##end new
              } else {
                tr_leftout = "ALL"
                potential_leftout = "ALL"
              }
              dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = tr_leftout, pot_over_tr = potential_leftout)
            }
            ##try returning the mean
            return(dt_return)
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }, mc.cores = cores) %>% do.call(rbind,.)
        if(select_type == "bases") {
          mean.overlap.dt[, mean_overlap := as.numeric(mean_overlap)]
          mean.overlap.dt[, pot_over_tr := as.numeric(pot_over_tr)]
          mean.overlap.dt[, sum_missing := (mean_overlap + pot_over_tr)]
          mean.overlap.dt[, rank_tr := rank(sum_missing,ties.method = "min"), by = "query"]
          mean.overlap.dt = mean.overlap.dt[order(mean_overlap, rank_tr),]
          if(minimize == "sum") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
          }
          if(minimize == "tr_bases_missing") {
            mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
          }
        }
        mean.overlap.dt2[, N_gene := .N, by = "query"]
        if(any(mean.overlap.dt2$N_gene > 1)) {
          warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
          select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
          select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
          select.dt[, N_gene := .N, by = "gene_name"]
          if(any(select.dt$N_gene > 1)) {
            warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
            select.dt = select.dt[which(width == min(width)),]
            rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
          }
          potent_tr = select.dt$transcript_id
        }
        if(exists("potent_tr")) {
          if(length(potent_tr) > 0) {
            potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
          }
        } else {
          ## potential_transcripts = mean.overlap.dt2$query
          potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
        }
      if(is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
      }
      if(!is.null(potential_transcripts_merge)) {
        add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
        add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
        md.novel.gr$gene_name = NULL
        md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
        md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
        if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
          warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
        }
      }
      ## now annotation the exons overlaps
      md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
      potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      }
    } else if(annotate_mismatch_type == "none") {
      message("skipped assigning a transcript to fusions. Or unassigned transcripts")
    }
    
    if(exists("md.novel.dt2") & !reannotate) {
    gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
    md.annotated = rbind(md.full,md.sub)
    md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
    md.annotated[, transcript_id := associated_transcript]
    md.novel.dt2[, reannotated := TRUE]
    md.annotated[, reannotated := FALSE]
    md.annotated = rbind(md.novel.dt2,md.annotated)
    } else if (reannotate) {
      md.annotated = md.novel.dt2
      md.annotated[, reannotated := TRUE]
    } else {
      gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique
      md.annotated = rbind(md.full,md.sub)
      md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
      md.annotated[, transcript_id := associated_transcript]
      md.annotated[, reannotated := FALSE]
    }
    md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
    potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
    potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
    potential.gtf.dt = as.data.table(potential.gtf.gr)
    
    add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
    md.annotated.dt2 = md.annotated
    md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
    ## now annotation the exons overlaps
    md_potential_tr.gr = GRanges(md.annotated.dt2, seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    ## md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
    potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
    potential.transcripts = potential.gtf.labels$transcript_id %>% unique
    md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
      md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
      potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
      ##breaks isn't working properly-convert these to breakpoints
      potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end,strand)]
      ##test end + 1 for pos
      bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt)))) %>% unique
      bps.dt[,start := end -1]
      bps.dt[start_end == "start", end := end - 1]
      bps.dt[start_end == "start", start := start - 1]
      bps.dt[start_end == "end", end := end + 1]
      bps.dt[start_end == "end", start := start + 1]
      ##test labels
      ## bps.dt[,start := start + 1]
      ## bps.dt[,end := start]
      ##
      bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths()) %>% gr.chr
      
      ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
      md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
      ## new gr.val to match to specific transcripts
      md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
      md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      return(md.novel.dt3)
    })
    md.annotated.dt3 = rbindlist(md.annotated.lst3)
    ##add back ones that were not assigned a transcript
    md.annotated.dt3 = rbind(md.annotated.dt3, md.annotated.dt2[transcript_id == "multiple_potential_genes",], fill = TRUE)
    ##
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
################
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
    md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                                  ifelse("UTR" %in% x, "UTR",
                                                                                  ifelse("exon" %in% x, "exon","intron"))))]
    md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3[type == "X",coding_type_simple := "del"]
    md.all.dt = md.annotated.dt3
    md.all.dt[,qname_simple := tstrsplit(qname, "_", keep = 1)]
    message("finished merging together")
    ## rbind all together
    ## new ordering using the fusion groups
    if(exists("fus.groups.dt")) {
      ## fus.groups.dt[, extra := as.list(extra)]
      ## test.fus.dt = fus.groups.dt[qname == "transcript/482749",.(row_number, chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, qname_list, qname, CB, AC)]
      ## test.md.gr = GRanges(md.all.dt[grepl("transcript/482749", qname),], seqlengths = hg_seqlengths())
      
      ##using the inherent position 1 and 2
      ## bp1.gr = test.fus.dt[,.(chr1,start1,end1)] %>% setnames(.,c("seqnames","start","end")) %>% GRanges(., seqlengths = hg_seqlengths())
      ## bp2.gr = test.fus.dt[,.(chr2,start2,end2)] %>% setnames(.,c("seqnames","start","end")) %>% GRanges(., seqlengths = hg_seqlengths())
      ## parsing the AC with the more than two if there are
      qnames.lst = fus.groups.dt$qname %>% tstrsplit(., "_", keep = 1) %>% unique %>% unlist
      qnames.lst = qnames.lst[qnames.lst %in% md.all.dt$qname_simple]
      md.lst = mclapply(qnames.lst, function(x) {
        sub.fus.dt = fus.groups.dt[qname == x,.(row_number, chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, qname_list, qname, CB, AC)]
        sub.md.gr = GRanges(md.all.dt[qname_simple == x,], seqlengths = hg_seqlengths())

        bps.lst = sub.fus.dt$AC %>% gsub("AC=","",.) %>% tstrsplit(.,"/") %>% unlist
        qnames.dt = lapply(1:length(bps.lst), function(x) {
          pad1 = pad ## just returning this at this point for debugging purposes
          bp.gr = parse.gr(bps.lst[x]) + pad
          int.gr = sub.md.gr %&% (bp.gr)
          ## if no intersection- try adding a pad a few times
          if(length(int.gr) == 0) {
            pad2 = pad * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad2
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          if(length(int.gr) == 0) {
            pad3 = pad2 * 10
            bp.gr2 = parse.gr(bps.lst[x]) + pad3
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          if(length(int.gr) == 0) {
            pad4 = pad3 * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad4
            int.gr = sub.md.gr %&% (bp.gr2)
          }
          qname1 = int.gr$qname %>% unique
          qnames.dt = data.table(qname = qname1, order = x)
          return(qnames.dt)
        }) %>% rbindlist(.)
        mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
        qnames.dt = qnames.dt[, .SD[which.min(order)], by = "qname"]
        if(nrow(qnames.dt) < length(bps.lst)) {
          ##get which qname is missing
          qnames.lst = sub.md.gr$qname %>% unique
          qname.missing.lst = qnames.lst[!(qnames.lst %in% qnames.dt$qname)]
          qname.missing.dt = data.table(qname = qname.missing.lst, order = 10)
          qnames.dt = rbind(qnames.dt, qname.missing.dt)
        }
        md.sub.dt = as.data.table(sub.md.gr)
        
        ## qname_1st = (test.md.gr %&% (bp1.gr+20))$qname %>% unique ## get the unique qname to paste in order
        ## qname_2nd = (test.md.gr %&% (bp2.gr+20))$qname %>% unique ## get the unique qname to paste in order
        
        md.dt4 = merge.data.table(md.sub.dt, qnames.dt, by = "qname", all.x = TRUE)
                                        #      md.dt4[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt4[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt4 = md.dt4[, c("start", "end") := .(ifelse(end < start, end, start), ifelse(end < start, start, end))]
        md.dt4 = md.dt4[order(order, order_start),]
        md.dt4[, qname := x]
        md.dt4[, final_order := 1:.N]
        return(md.dt4)
      }, mc.cores = cores)
      md.dt4 = rbindlist(md.lst)
      
      ## md.dt4 = md.dt4[order(order, order_id, order_start),]
      ## order grl.ix after ordering the qnames
    } else {
      ## sort based on strand to have walks align
      ## fix with ordering for fusions
      md.all.dt$strand = factor(md.all.dt$strand, levels = c("+","-"))
      ##pos.dt = md.all.dt[strand == "+",][order(grl.ix, strand, order_id),][qname == "transcript/611213",]
      md.all.dt[, order_start := ifelse(strand == "+", start, -start)]
      md.dt4 = md.all.dt[order(order_softclipped, strand, order_start),]
      ## md.dt4 = md.all.dt[order(strand, order_start, grl.ix,  order_id),] ##[qname == "transcript/611213",]

##### OLD
      ## md.all.dt[order(strand, grl.ix,  order_id),][qname == "transcript/611213",]
      ## pos.dt = md.all.dt[strand == "+",][order(order_id),]
      ## neg.dt = md.all.dt[strand == "-",][order(order_id),]
      ## pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
      ## neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
      ## md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
      ## END OLD
    }
    if(add_gencode_name) {
      gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
      gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
      md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
    }
    md.dt4 = md.dt4[order(qname,final_order),]
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    ## convert to grl and gw
    message("converting to grl")
    md.gr2$qr_tr = paste0(gsub("transcript/","",md.gr2$qname),"; ", md.gr2$transcript_name)
    md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
    message("done coverting to grl")
    ## now annotating missing parts
    if(annotate_missing) {
      message("annotating missing regions from transcripts")
      md.all.dt2 = mclapply(1:length(md.grl2), function(x) {
        md.sub.tr = md.grl2[[x]]
        qr_tr1 = md.sub.tr$qr_tr[1]
        qname_tr1 = md.sub.tr$qname[1]
        tr.gr = potential.gtf.dt[transcript_id == unique(md.sub.tr$transcript_id),][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
                                        #tr.missing.gr = tr.gr[tr.gr %outside% md.sub.tr]
                                        #mcols(tr.missing.gr) = NULL
        ##attempt not using outside
        tiled.gr = gr.tile(tr.gr,width = 1)
        tiled.gr$overlap = tiled.gr %O% md.sub.tr
        tr.missing.gr = (tiled.gr %Q% (overlap < 1)) %>% gr.reduce
        if(length(tr.missing.gr) == 0) {
          return(as.data.table(md.sub.tr))
        }
        tr.missing.gr = as.data.table(tr.missing.gr) %>% GRanges(., seqlengths = hg_seqlengths())
        ##end  attempt
        ##add meta data
        meta.add.dt = mcols(md.sub.tr)[1,] %>% as.data.table
        ##meta.add.dt[, c("qid", "type","rid","riid","fid","col",) := NULL]
        if("isoform" %in% names(meta.add.dt)) {
          cols_keep = c("transcript_id", "qname", "isoform", "structural_category", "subcategory", "associated_transcript", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        } else {
          cols_keep = c("transcript_id", "qname", "gene_name", "reannotated", "transcript_type", "transcript_name", "qr_tr")
        }
        ## meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
        meta.add.dt = meta.add.dt[,..cols_keep]
        tr.missing.gr = unique(tr.missing.gr)
        mcols(tr.missing.gr) = meta.add.dt
        ##
        strand_gene = md.sub.tr@strand@values
        tr.missing.gr@strand@values = strand_gene
        ##
        md.sub.tr2 = c(md.sub.tr,tr.missing.gr)
        md.sub.tr.dt = as.data.table(md.sub.tr2)
        md.sub.tr.dt[is.na(type), type := "missing"]
        md.sub.tr.dt[is.na(coding_type_simple), coding_type_simple := "missing"]
        pos.dt = md.sub.tr.dt[strand == "+" | strand == "*",][order(seqnames,start,end),]
        neg.dt = md.sub.tr.dt[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
        return(rbind(pos.dt,neg.dt))
      }, mc.cores = cores)
      md.all.dt2 = rbindlist(md.all.dt2, fill = TRUE)
    } else {
      md.all.dt2 = md.grl2 %>% unlist %>% as.data.table
      md.all.dt2 = md.all.dt2[order(qname,final_order),]
    }
    message("done annotating missing regions")
    md.gr3 = GRanges(md.all.dt2, seqlengths = hg_seqlengths())
    md.grl3 = rtracklayer::split(md.gr3, f = mcols(md.gr3)["qname"])
  } else {
    md.dt5 = order_fusion_breakpoint(fusion_groups.dt = fus.groups.dt, md.all.dt = md.dt4, cores = cores, pad = pad)
    md.dt5 = md.dt5[!(type %in% c("N","X", "S")),]
    md.gr2 = GRanges(md.dt5, seqlengths = hg_seqlengths())
    md.grl3 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
    if(type == "grl") {
      return(md.grl3)
    } else if (type == "gw") {
      md.gw = gW(grl = md.grl3)
      transcript_names = sapply(md.gw$grl, function(x) x$qname_simple[1])
      md.gw$set(name = transcript_names)
      message("done converting grl to gW object. returning gW")
      return(md.gw)
    }
  }
  ## sort positive strand forward and negative strand reverse
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    md.gw = gW(grl = md.grl3)
    ## add coloring of nodes to match track.gencode
    cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
    for(x in 1:nrow(cmap.dt)) {
      cmap.sub = cmap.dt[x,]
      md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
    }
    transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
    md.gw$set(name = transcript_names)
    message("done converting grl to gW object. returning gW")
    return(md.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(md.grl3)
  }
}

##
order_fusion_breakpoint = function(fusion_groups.dt, md.all.dt, pad = 10, cores = 1) {
      qnames.lst = fusion_groups.dt$qname %>% tstrsplit(., "_", keep = 1) %>% unique %>% unlist
      md.all.dt[, qname_simple := tstrsplit(qname, "_", keep = 1)]
      qnames.lst = qnames.lst[qnames.lst %in% md.all.dt$qname_simple]
      ## qnames.lst = qnames.lst[qnames.lst %in% md.all.dt$qname]
      md.lst = mclapply(qnames.lst, function(x) {
        sub.fus.dt = fusion_groups.dt[qname == x,.(row_number, chr1, start1, end1, chr2, start2, end2, id, score, strand1, strand2, info, qname_list, qname, CB, AC)]
        sub.md.gr = GRanges(md.all.dt[qname_simple == x,], seqlengths = hg_seqlengths())
        sub.md.gr = sub.md.gr %Q% (!type %in% c("N","X", "S"))

        bps.lst = sub.fus.dt$AC %>% gsub("AC=","",.) %>% tstrsplit(.,"/") %>% unlist
        bps.lst = sub.fus.dt$AC %>% gsub("AC=","",.) %>% tstrsplit(.,"/|,") %>% unlist
        qnames.dt = lapply(1:length(bps.lst), function(x) {
          pad1 = pad ## just returning this at this point for debugging purposes
          bp.gr = parse.gr(bps.lst[x]) + pad
          int.gr = sub.md.gr %&% (bp.gr)
          if(length(int.gr) != 0) {
            int.gr$attempt = 1
          }
          ## if no intersection- try adding a pad a few times
          if(length(int.gr) == 0) {
            pad2 = pad * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad2
            int.gr = sub.md.gr %&% (bp.gr2)
            if(length(int.gr) != 0) {
              int.gr$attempt = 2
            }
          }
          if(length(int.gr) == 0) {
            pad3 = pad2 * 10
            bp.gr2 = parse.gr(bps.lst[x]) + pad3
            int.gr = sub.md.gr %&% (bp.gr2)
            if(length(int.gr) != 0) {
              int.gr$attempt = 2
            }
          }
          if(length(int.gr) == 0) {
            pad4 = pad3 * 100
            bp.gr2 = parse.gr(bps.lst[x]) + pad4
            int.gr = sub.md.gr %&% (bp.gr2)
            if(length(int.gr) != 0) {
              int.gr$attempt = 2
            }
          }
          if(length(int.gr) == 0) {
            return(NULL)
          }
          qname1 = int.gr$qname %>% unique
          qnames.dt = data.table(qname = qname1, order = x, attempt = unique(int.gr$attempt))
          return(qnames.dt)
        })
        qnames.dt = rbindlist(qnames.dt, fill = TRUE)
        ## mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
        qnames.dt = qnames.dt[, .SD[which.min(attempt)], by = "qname"]
        qnames.dt = qnames.dt[, .SD[which.min(order)], by = "qname"]
        if(nrow(qnames.dt) < length(bps.lst)) {
          ##get which qname is missing
          qnames.lst = sub.md.gr$qname %>% unique
          qname.missing.lst = qnames.lst[!(qnames.lst %in% qnames.dt$qname)]
          qname.missing.dt = data.table(qname = qname.missing.lst, order = 10)
          qnames.dt = rbind(qnames.dt, qname.missing.dt, fill = TRUE)
        }
        md.sub.dt = as.data.table(sub.md.gr)
        
        ## qname_1st = (test.md.gr %&% (bp1.gr+20))$qname %>% unique ## get the unique qname to paste in order
        ## qname_2nd = (test.md.gr %&% (bp2.gr+20))$qname %>% unique ## get the unique qname to paste in order
        
        md.dt3 = merge.data.table(md.sub.dt, qnames.dt, by = "qname", all.x = TRUE)
                                        #      md.dt3[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt3[, order_start := ifelse(strand == "+", order_id, -order_id)]
        md.dt3 = md.dt3[, c("start", "end") := .(ifelse(end < start, end, start), ifelse(end < start, start, end))]
        md.dt3 = md.dt3[order(order, order_start),]
        md.dt3[, qname := x]
        md.dt3[, final_order := 1:.N]
        return(md.dt3)
      }, mc.cores = cores)
      md.dt3 = rbindlist(md.lst)
      return(md.dt3)
      ## md.dt4 = md.dt4[order(order, order_id, order_start),]
      ## order grl.ix after ordering the qnames
}



collect_transcript_info = function(isoform_counts,
                                   collapsed_class,
                                   read_stat,
                                   class_subset = c("isoform", "structural_category", "subcategory")) {
  count.dt = fread(isoform_counts) %>% setnames(., c("isoform", "isoform_count"))
  stat.dt = fread(read_stat) %>% setnames(., c("id","length","isoform"))
  count.dt = merge.data.table(count.dt, stat.dt, by = "isoform")
  ## get the classification for each isoform
  if(!is.null(class_subset)) {
    class.dt = fread(collapsed_class)[,..class_subset]
  } else {
    class.dt = fread(collapsed_class)
  }
  ##add length as a list
  count.dt[, length := as.integer(length)]
  count.dt[, isoform := as.character(isoform)]
  ##
  count.dt[, length_per_read := toString(length), by = "isoform"]
  ##  count.dt[, length_per_read := list(list(length)), by = "isoform"]
  count.dt[, c("length","id") := NULL]
  count.dt = unique(count.dt)
  ##count.dt[, unique_id := .GRP, by = "isoform"]
  ##make unique isoform data.table and then merge back the length of the transcripts
  ## sub.dt = count.dt[,.(isoform,isoform_count,unique_id)] %>% unique
  ## uniq.length.dt = copy(count.dt)
  ## uniq.length.dt[,c("isoform_count","unique_id") := NULL]
  ## setkey(uniq.length.dt,isoform)
  ## uniq.length.dt[, length_per_read := sapply(length_per_read, toString), by = isoform]
  ## uniq.length.dt = unique(uniq.length.dt)
  ## merge.dat.table(sub.dt,  
  ## uniq.count.dt = count.dt[unique(unique_id),]
  ## uniq.count.dt[, unique_id := NULL]
  ## merge.dt = merge.data.table(class.dt, uniq.count.dt, by = "isoform", all = TRUE)
  merge.dt = merge.data.table(class.dt, count.dt, by = "isoform", all = TRUE)
  return(merge.dt)
}








################################
## old get_iso_reads- 5-24-24
## ## ## reorganize this function- should have same labeling for previously annotated and newly annotated; novel should just return the matching transcripts
## get_iso_reads = function(bam,
##                          gr,
##                          gtf,
##                          collapsed_group = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
##                          collapsed_class = NULL, #pacbio pigeon annotation output, if reannotate = FALSE and this, will look for 
##                          type = "gw",
##                          annotate_mismatch = TRUE,
##                          annotate_mismatch_type = "smart",
##                          reann_with = "both", #whether to match using both exons and UTRS or just exon ("exon")
##                          reannotate = FALSE,
##                          select_type = "bases", #only bases works for annotate_mismatch_type = "smart"
##                          add_gencode_name = TRUE, #add the gencode name to the transcript for plotting
##                          reverse_overlap = FALSE, #old method - need to remove
##                          consider = "only_exon_alignment", #only applies to select_type = "fraction"
##                          read_over_potential = 0.5, #only applies to annotate_mismatch_type = "percent"
##                          potential_over_read = 0.8, #only applies to annotate_mismatch_type = "percent"
##                          minimize = "tr_bases_missing", #only applies to annotate_mismatch_type = smart, typically works better then sum
##                          annotate_missing = TRUE,
##                          cores = 1) {
##   message(paste0("Reading in ",bam))
##   md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
##   if(length(md.gr.bam) == 0) {
##     warning("no reads in the specified region. Returning null")
##     return(NULL)
##   }
##   message("Done reading")
##   md.dt = as.data.table(md.gr.bam)
##   md.dt = md.dt[width != 0,]
##   md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##   message("Splicing cigar string")
##   md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
##   message("Done Splicing cigar string")
##   if(length(md.grl) == 0) {
##     return(NULL)
##   }
##   ## gencode.gr = gencode@data[[1]] %>% unlist
##   gencode.gr = gtf
##   gencode.gr$type2 = gencode.gr$type
##   ## if(exists("potential_transcript_merge")) {
##   potential_transcript_merge = NULL
##   ## }
##   if(!reannotate) {
##     ##look for collapsed_group and collapsed_class if missing
##     if(is.null(collapsed_group) | is.null(collapsed_class)) {
##       if(is.null(collapsed_group) && is.null(collapsed_class)) {
##         message("looking for collapsed_group and collapsed class")
##         folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
##         collapsed_group = paste0(folder_find,"collapsed.group.txt")
##         collapsed_class = paste0(folder_find,"collapsed_classification.txt")
##         if(all(file.exists(c(collapsed_group,collapsed_class)))) {
##           message("Found! collapsed.group.txt and collapsed_classification.txt")
##         }
##       }
##     }
##     group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
##     class.dt = fread(collapsed_class)
##     ## get labels for already identified transcripts
##     group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
##     group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
##     group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
##     expanded.dt = copy(group.sub.dt)
##     expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
##     group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
##     ##add transcript labels to bam reads
##     md.dt2 = unlist(md.grl) %>% as.data.table()
##     if(nrow(md.dt2) > 0) { ## necessary for fusions where md.dt2 is probably empty
##       md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
##     } else {
##       md.dt3 = group.sub.dt2
##     }
    
##   } else if(reannotate) {
##     md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
##   }
##   if(annotate_mismatch) {
##     if(!reannotate) {
##       md.full = md.dt3[structural_category == "full-splice_match",]
##       md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##       md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
##     } else if(reannotate) {
##       md.novel = md.dt3
##       md.full = NULL
##     }
## ###################################
##     ## first  label ones without a transcript annotation
##     md.novel.gr = GRanges(md.novel)
##     if(class(gtf) == "character") {
##       message(paste0("Reading in ",gtf)) 
##       gtf = rtracklayer::readGFF(gtf) %>% as.data.table
##       message(paste0("Done reading in ",gtf))
##     }
##     if(class(gtf)[1] != "data.table" & class(gtf)[1] != "GRanges") {
##       stop("Gtf must be a character to a gtf file or a data.table")
##     }
##     if(any(class(gtf) != "GRanges")) {
##       gtf.gr = GRanges(gtf, seqlengths = hg_seqlengths())
##     } else {
##       gtf.gr = gtf
##     }
##     ##add potential genes to md.novel
##     if(length(md.novel.gr) > 0) {
##       md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##       potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##       potential.gtf.dt = as.data.table(potential.gtf.gr)
##       ## annotate to the longest transcript??
##       if(annotate_mismatch_type == "longest") {
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "level") {
##         ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "match") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         ##potential one
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
##         ## get the transcript with the most N by transcript
##         potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "percent") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
##         potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
##         potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
##         ## get maximum overlap of each transcript with each potential transcript
##         md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        
##         message("getting overlap of transcripts with potential transcripts")
##         mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             ## annotate with utr and exons or just exons
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             if(select_type == "fraction") {
##               if(consider == "only_exon_alignment") {
##                 md.sub.gr = md.sub.gr %&% potential.sub.gr
##               }
##               md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##               ## test reverse
##               md.sub.gr2 = md.novel.grl[[tr]]
##               potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
##               ##
##               dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
##             } else if(select_type == "bases") {
##               ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
##               ## %o% is not working hmmm
##               x = md.sub.gr
##               y = potential.sub.gr
##               ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
##               if (nrow(ov)>0){
##                 ov = ov[ , sum(width), keyby = query.id]
##                 x$width.ov = 0
##                 x$width.ov[ov$query.id] = ov$V1
##                 md.sub.gr$percent = x$width.ov
##               } else {
##                 md.sub.gr$percent = 0
##               }
##               ## test reverse
##               md.sub.gr2 = md.novel.grl[[tr]]
##               potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr2
##               ##
##               dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE), pot_over_tr = mean(potential.sub.gr$percent, na.rm = TRUE))
##               ## dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
##             }
##             ##try returning the mean
##             return(dt_return)
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt)
##         }, mc.cores = cores) %>% do.call(rbind,.)
##         if(select_type == "fraction") {
##           ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
##           mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
##           mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
##           ## mean.overlap.dt[, rank_tr2 := rank(pot_over_tr,-mean_overlap,ties.method = "min"), by = "query"]
##           ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##           ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##         }
##         if(select_type == "bases") {
##           ## mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5 & pot_over_tr > 0.2,][order(-mean_overlap)]
##           mean.overlap.dt = mean.overlap.dt[mean_overlap > read_over_potential & pot_over_tr > potential_over_read,][order(-mean_overlap)]
##           mean.overlap.dt[, rank_tr := rank(-mean_overlap,-pot_over_tr,ties.method = "min"), by = "query"]
##           ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
##           ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##         }
##         ## need to rewrite to loop through the correct combinations
##         if(reverse_overlap) { # this option should be removed- ultimately in the loop above
##           message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
##           rows.cycle = mean.overlap.dt[rank_tr <= 20,]
##           mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
##             rows.cycle.sub = rows.cycle[x,]
##             tr = rows.cycle.sub$query
##             pot.tr = rows.cycle.sub$subject
##             md.sub.gr = md.novel.grl[[tr]]          
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             if(select_type == "fraction") {
##               potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
##             } else if(select_type == "bases") {
##               ## potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
##               ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
##               ## %o% is not working hmmm
##               x = potential.sub.gr
##               y = md.sub.gr
##               ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
##               if (nrow(ov)>0){
##                 ov = ov[ , sum(width), keyby = query.id]
##                 x$width.ov = 0
##                 x$width.ov[ov$query.id] = ov$V1
##                 potential.sub.gr$percent = x$width.ov
##               } else {
##                 potential.sub.gr$percent = 0
##               }

##             }
##             ##try returning the mean
##             return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           if(select_type == "fraction") {
##             mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
##             mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
##             names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
##             mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
##             mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
##           }
##           if(select_type == "bases") {
##             mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
##             mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
##             names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
##             mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
##             mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
##           }
##         } else {
##           mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
##         }
##         mean.overlap.dt2[, N_gene := .N, by = "query"]
##         if(any(mean.overlap.dt2$N_gene > 1)) {
##           warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
##           select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
##           select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
##           select.dt[, N_gene := .N, by = "gene_name"]
##           if(any(select.dt$N_gene > 1)) {
##             warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
##             select.dt = select.dt[which(width == min(width)),]
##             rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
##           }
##           potent_tr = select.dt$transcript_id
##         }
##         if(exists("potent_tr")) {
##           if(length(potent_tr) > 0) {
##             potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
##           }
##         } else {
##           ## potential_transcripts = mean.overlap.dt2$query
##           potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
##         }
##         if(is.null(potential_transcripts_merge)) {
##           add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
##           md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##           md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##         }
##         if(!is.null(potential_transcripts_merge)) {
##           add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
##           add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
##           md.novel.gr$gene_name = NULL
##           md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
##           md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
##           if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
##             warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
##           }
##         }
##         ## now annotation the exons overlaps
##         md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##         potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##         potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       } else if(annotate_mismatch_type == "smart") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
##         potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
##         potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
##         ## get maximum overlap of each transcript with each potential transcript
##         md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
##         message("getting overlap of transcripts with potential transcripts")
##         mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             ## annotate with utr and exons or just exons
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             if(select_type == "bases") {
##               ## md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
##               ## %o% is not working hmmm
##               x = md.sub.gr
##               y = potential.sub.gr
##               ov = gr2dt(gr.findoverlaps(x, reduce(gr.stripstrand(y))))
##               if (nrow(ov)>0) {
##                 ov = ov[ , sum(width), keyby = query.id]
##                 x$width.ov = 0
##                 x$width.ov[ov$query.id] = ov$V1
##                 ##new part for annotate_mismatch_type == "smart"
##                 y_sum_bases = y@ranges@width %>% sum
##                 x_sum_bases = x@ranges@width %>% sum
##                 ov_sum_bases = x$width.ov %>% sum
##                 ## ov_sum_bases = ov$width %>% sum
##                 potential_leftout = y_sum_bases - ov_sum_bases
##                 tr_leftout = x_sum_bases - ov_sum_bases
##                 ##end new
##               } else {
##                 tr_leftout = "ALL"
##                 potential_leftout = "ALL"
##               }
##               dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = tr_leftout, pot_over_tr = potential_leftout)
##             }
##             ##try returning the mean
##             return(dt_return)
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt)
##         }, mc.cores = cores) %>% do.call(rbind,.)
##         if(select_type == "bases") {
##           mean.overlap.dt[, mean_overlap := as.numeric(mean_overlap)]
##           mean.overlap.dt[, pot_over_tr := as.numeric(pot_over_tr)]
##           mean.overlap.dt[, sum_missing := (mean_overlap + pot_over_tr)]
##           mean.overlap.dt[, rank_tr := rank(sum_missing,ties.method = "min"), by = "query"]
##           mean.overlap.dt = mean.overlap.dt[order(mean_overlap, rank_tr),]
##           if(minimize == "sum") {
##             mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(rank_tr)], by = query]
##           }
##           if(minimize == "tr_bases_missing") {
##             mean.overlap.dt2 = mean.overlap.dt[, .SD[which.min(mean_overlap)], by = query]
##           }
##         }
##         mean.overlap.dt2[, N_gene := .N, by = "query"]
##         if(any(mean.overlap.dt2$N_gene > 1)) {
##           warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
##           select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
##           select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
##           select.dt[, N_gene := .N, by = "gene_name"]
##           if(any(select.dt$N_gene > 1)) {
##             warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
##             select.dt = select.dt[which(width == min(width)),]
##             rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
##           }
##           potent_tr = select.dt$transcript_id
##         }
##         if(exists("potent_tr")) {
##           if(length(potent_tr) > 0) {
##             potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
##           }
##         } else {
##           ## potential_transcripts = mean.overlap.dt2$query
##           potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
##         }
##       if(is.null(potential_transcripts_merge)) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##       }
##       if(!is.null(potential_transcripts_merge)) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
##         add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
##         md.novel.gr$gene_name = NULL
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
##         if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
##           warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
##         }
##       }
##       ## now annotation the exons overlaps
##       md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##       potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##       potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       }
##     }
##     if(exists("md.novel.dt2") & !reannotate) {
##     gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
##     md.annotated = rbind(md.full,md.sub)
##     md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
##     md.annotated[, transcript_id := associated_transcript]
##     md.novel.dt2[, reannotated := TRUE]
##     md.annotated[, reannotated := FALSE]
##     md.annotated = rbind(md.novel.dt2,md.annotated)
##     } else if (reannotate) {
##       md.annotated = md.novel.dt2
##       md.annotated[, reannotated := TRUE]
##     } else {
##       gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique
##       md.annotated = rbind(md.full,md.sub)
##       md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
##       md.annotated[, transcript_id := associated_transcript]
##       md.annotated[, reannotated := FALSE]
##     }
##     md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
##     potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##     potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##     potential.gtf.dt = as.data.table(potential.gtf.gr)
    
##     add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
##     md.annotated.dt2 = md.annotated
##     md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##     ## now annotation the exons overlaps
##     md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##     potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##     potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##     potential.transcripts = potential.gtf.labels$transcript_id %>% unique
##     md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
##       md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
##       potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
##       ##breaks isn't working properly-convert these to breakpoints
##       potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end,strand)]
##       ##test end + 1 for pos
##       bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end), strand = as.character(potential.gtf.labels.sub.dt$strand, potential.gtf.labels.sub.dt$strand), start_end = c(rep("start",nrow(potential.gtf.labels.sub.dt)),rep("end",nrow(potential.gtf.labels.sub.dt)))) %>% unique
##       bps.dt[,start := end -1]
##       bps.dt[start_end == "start", end := end - 1]
##       bps.dt[start_end == "start", start := start - 1]
##       bps.dt[start_end == "end", end := end + 1]
##       bps.dt[start_end == "end", start := start + 1]
##       ##test labels
##       ## bps.dt[,start := start + 1]
##       ## bps.dt[,end := start]
##       ##
##       bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths(chr = FALSE)) %>% gr.chr
      
##       ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
##       md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
##       ## new gr.val to match to specific transcripts
##       md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
##       md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       return(md.novel.dt3)
##     })
##     md.annotated.dt3 = rbindlist(md.annotated.lst3)  
##     md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
## ################
##     md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##     md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                                   ifelse("UTR" %in% x, "UTR",
##                                                                                   ifelse("exon" %in% x, "exon","intron"))))]
##     md.annotated.dt3 = md.annotated.dt3[type != "N",]
##     md.annotated.dt3[type == "X",coding_type_simple := "del"]
##     md.all.dt = md.annotated.dt3
##     message("finished merging together")
    
##     ## rbind all together
##     ## sort based on strand to have walks align
##     pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
##     neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
##     md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
##     if(add_gencode_name) {
##       gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
##       gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
##       md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
##     }
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     ## convert to grl and gw
##     message("converting to grl")
##     md.gr2$qr_tr = paste0(gsub("transcript/","",md.gr2$qname),"; ", md.gr2$transcript_name)
##     md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
##     message("done coverting to grl")
##     ## now annotating missing parts
##     if(annotate_missing) {
##       message("annotating missing regions from transcripts")
##       md.all.dt2 = mclapply(1:length(md.grl2), function(x) {
##         md.sub.tr = md.grl2[[x]]
##         qr_tr1 = md.sub.tr$qr_tr[1]
##         qname_tr1 = md.sub.tr$qname[1]
##         tr.gr = potential.gtf.dt[transcript_id == unique(md.sub.tr$transcript_id),][type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##                                         #tr.missing.gr = tr.gr[tr.gr %outside% md.sub.tr]
##                                         #mcols(tr.missing.gr) = NULL
##         ##attempt not using outside
##         tiled.gr = gr.tile(tr.gr,width = 1)
##         tiled.gr$overlap = tiled.gr %O% md.sub.tr
##         tr.missing.gr = (tiled.gr %Q% (overlap < 1)) %>% gr.reduce
##         if(length(tr.missing.gr) == 0) {
##           return(as.data.table(md.sub.tr))
##         }
##         tr.missing.gr = as.data.table(tr.missing.gr) %>% GRanges(., seqlengths = hg_seqlengths())
##         ##end  attempt
##         ##add meta data
##         meta.add.dt = mcols(md.sub.tr)[1,] %>% as.data.table
##         ##meta.add.dt[, c("qid", "type","rid","riid","fid","col",) := NULL]
##         meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
##         tr.missing.gr = unique(tr.missing.gr)
##         mcols(tr.missing.gr) = meta.add.dt
##         ##
##         strand_gene = md.sub.tr@strand@values
##         tr.missing.gr@strand@values = strand_gene
##         ##
##         md.sub.tr2 = c(md.sub.tr,tr.missing.gr)
##         md.sub.tr.dt = as.data.table(md.sub.tr2)
##         md.sub.tr.dt[is.na(type), type := "missing"]
##         md.sub.tr.dt[is.na(coding_type_simple), coding_type_simple := "missing"]
##         pos.dt = md.sub.tr.dt[strand == "+" | strand == "*",][order(seqnames,start,end),]
##         neg.dt = md.sub.tr.dt[strand == "-" | strand == "*",][order(seqnames,-start,-end),]
##         return(rbind(pos.dt,neg.dt))
##       }, mc.cores = cores)
##       md.all.dt2 = rbindlist(md.all.dt2, fill = TRUE)
##     } else {
##       md.all.dt2 = md.grl2 %>% unlist %>% as.data.table
##     }
##     message("done annotating missing regions")
##     md.gr3 = GRanges(md.all.dt2, seqlengths = hg_seqlengths())
##     md.grl3 = rtracklayer::split(md.gr3, f = mcols(md.gr3)["qname"])
##   } else {
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     md.grl3 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     return(md.grl3)
##   }
##   ## sort positive strand forward and negative strand reverse
##   if(type == "gw") {
##     message("type is 'gw' converting grl to gW object")
##     md.gw = gW(grl = md.grl3)
##     ## add coloring of nodes to match track.gencode
##     cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del", "missing"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange", "orange"))
##     for(x in 1:nrow(cmap.dt)) {
##       cmap.sub = cmap.dt[x,]
##       md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
##     }
##     transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
##     md.gw$set(name = transcript_names)
##     message("done converting grl to gW object. returning gW")
##     return(md.gw)
##   } else if(type == "grl") {
##     message("type is 'grl' returning grl")
##     return(md.grl3)
##   }
## }






## ## ## reorganize this function- should have same labeling for previously annotated and newly annotated; novel should just return the matching transcripts
## get_iso_reads = function(bam, gr, gtf, collapsed_group = NULL, collapsed_class = NULL, type = "gw", annotate_mismatch = TRUE, annotate_mismatch_type = "percent",reann_with = "exon", reannotate = FALSE, select_type = "fraction", add_gencode_name = TRUE, cores = 1) {
##   message(paste0("Reading in ",bam))
##   md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
##   message("Done reading")
##   md.dt = as.data.table(md.gr.bam)
##   ## pos.dt = md.dt[strand == "+",][order(seqnames,start,end),]
##   ## neg.dt = md.dt[strand == "-",][order(seqnames,-start,-end),]
##   ## md.dt = rbind(pos.dt,neg.dt)
##   md.dt = md.dt[width != 0,]
##   md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##   message("Splicing cigar string")
##   md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
##   message("Done Splicing cigar string")
##   ## gencode.gr = gencode@data[[1]] %>% unlist
##   gencode.gr = gtf
##   gencode.gr$type2 = gencode.gr$type
##   if(exists("potential_transcript_merge")) {
##     potential_transcript_merge = NULL
##   }
##   if(!reannotate) {
##     ##look for collapsed_group and collapsed_class if missing
##     if(is.null(collapsed_group) | is.null(collapsed_class)) {
##       if(is.null(collapsed_group) && is.null(collapsed_class)) {
##         message("looking for collapsed_group and collapsed class")
##         folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
##         collapsed_group = paste0(folder_find,"collapsed.group.txt")
##         collapsed_class = paste0(folder_find,"collapsed_classification.txt")
##         if(all(file.exists(c(collapsed_group,collapsed_class)))) {
##           message("Found! collapsed.group.txt and collapsed_classification.txt")
##         }
##       }
##     }
##     group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
##     class.dt = fread(collapsed_class)
##     ## get labels for already identified transcripts
##     group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
##     group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
##     group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
##     expanded.dt = copy(group.sub.dt)
##     expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
##     group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
##     ##add transcript labels to bam reads
##     md.dt2 = unlist(md.grl) %>% as.data.table()
##     md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")[type != "N",]
##   } else if(reannotate) {
##     ## md.dt3 = md.dt[,.(qname,seqnames,start,end,width,strand,type)]
##     md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
##   }
##   if(annotate_mismatch) {
##     if(!reannotate) {
##       md.full = md.dt3[structural_category == "full-splice_match",]
##       md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##       md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
##     } else if(reannotate) {
##       md.novel = md.dt3
##       md.full = NULL
##     }
## ###################################
##     ## first  label ones without a transcript annotation
##     md.novel.gr = GRanges(md.novel)
##     if(class(gtf) == "character") {
##       message(paste0("Reading in ",gtf)) 
##       gtf = rtracklayer::readGFF(gtf) %>% as.data.table
##       message(paste0("Done reading in ",gtf))
##     }
##     if(class(gtf)[1] != "data.table" & class(gtf)[1] != "GRanges") {
##       stop("Gtf must be a character to a gtf file or a data.table")
##     }
##     if(any(class(gtf) != "GRanges")) {
##       gtf.gr = GRanges(gtf, seqlengths = hg_seqlengths())
##     } else {
##       gtf.gr = gtf
##     }
##     ##add potential genes to md.novel
##     if(length(md.novel.gr) > 0) {
##       md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##       potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##       potential.gtf.dt = as.data.table(potential.gtf.gr)
##       ## annotate to the longest transcript??
##       if(annotate_mismatch_type == "longest") {
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "level") {
##         ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "match") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         ##potential one
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
##         ## get the transcript with the most N by transcript
##         potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "percent") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
##         potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
##         potential_transcripts.grl = rtracklayer::split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
##         ## get maximum overlap of each transcript with each potential transcript
##         md.novel.grl = rtracklayer::split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
##         message("getting overlap of transcripts with potential transcripts")
##         mean.overlap.dt = mclapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             ## annotate with utr and exons or just exons
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             if(select_type == "fraction") {
##               md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##               dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE))
##             } else if(select_type == "bases") {
##               md.sub.gr$percent = md.sub.gr %o% potential.sub.gr
##               dt_return = data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = sum(md.sub.gr$percent, na.rm = TRUE))
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             ##try returning the mean
##             return(dt_return)
##             ## return(data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt)
##         }, mc.cores = cores) %>% do.call(rbind,.)
##         ## mean.overlap.dt = mean.overlap.dt[, .SD[which.max(mean_overlap)], by = query]
##         if(select_type == "fraction") {
##           mean.overlap.dt = mean.overlap.dt[mean_overlap > 0.5,][order(-mean_overlap)]
##           mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##           #mean.overlap.dt = mean.overlap.dt[rank_tr <= 20,]
##         }
##         if(select_type == "bases") {
##           mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
##           mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##         }
##         ## mean.overlap.dt = mean.overlap.dt[, .SD[which.max(mean_overlap)], by = query]        
##         message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
##         ## mean.overlap.dt2 = mclapply(names(md.novel.grl), function(tr) {
##         ## need to rewrite to loop through the correct combinations
##         rows.cycle = mean.overlap.dt[rank_tr <= 20,]
##         ## mean.overlap.dt2 = mclapply(names(md.novel.grl), function(tr) {
##         mean.overlap.dt2 = mclapply(1:nrow(rows.cycle), function(x) {
##           rows.cycle.sub = rows.cycle[x,]
##           tr = rows.cycle.sub$query
##           pot.tr = rows.cycle.sub$subject
##           md.sub.gr = md.novel.grl[[tr]]          
##           potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##           if(reann_with == "both") {
##             potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##           } else if (reann_with == "exon") {
##             potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##           }
##           ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##           if(select_type == "fraction") {
##             potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
##           } else if(select_type == "bases") {
##             potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
##           }
##           ##try returning the mean
##           return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
##         }) %>% do.call(rbind,.)
##         ## mean.overlap.dt2 = mclapply(names(md.novel.grl), function(tr) {
##         ##   md.sub.gr = md.novel.grl[[tr]]
##         ##   mean.overlap.dt2 = lapply(mean.overlap.dt[rank_tr <= 20,]$subject, function(pot.tr) {
##         ##     potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##         ##     if(reann_with == "both") {
##         ##       potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##         ##     } else if (reann_with == "exon") {
##         ##       potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##         ##     }
##         ##     ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##         ##     if(select_type == "fraction") {
##         ##       potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
##         ##     } else if(select_type == "bases") {
##         ##       potential.sub.gr$percent = potential.sub.gr %o% md.sub.gr
##         ##     }
##         ##     ##try returning the mean
##         ##     return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
##         ##   }) %>% do.call(rbind,.)
##         ##   return(mean.overlap.dt2)
##         ## }, mc.cores = cores) %>% do.call(rbind,.)
##         ## mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = query] # used to be selecting for best transcript but I actually want to assign the best transcript to each read
##         if(select_type == "fraction") {
##           mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)] %>% unique
##           mean.overlap.dt2[order(-mean_overlap), rank_tr := 1:.N, by = "subject"]
##           ##mean.overlap.dt2 = mean.overlap.dt2[rank_tr <= 20,]
##           names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
##           ## mean.merged.dt = merge.data.table(mean.overlap.dt2, mean.overlap.dt, by = c("subject","query","gene_name"), suffixes = c("over_ref","over_tr"), all = TRUE)
##           mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
##           mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
##           ## mean.merged.dt[, avg_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = c("rank_trover_ref", "rank_trover_tr")]
##           ## mean.merged.dt[, .SD[which.max(rank_overlapover_ref)], by = query]
##           ## mean.merged.dt[, low_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = c("rank_trover_ref", "rank_trover_tr")]
##           ## mean.merged.dt = mean.merged.dt[, .SD[which.max(avg_rank)], by = subject]
##         }
##         if(select_type == "bases") {
##           mean.overlap.dt2 = mean.overlap.dt2[mean_overlap > 0.5,][order(-mean_overlap)]
##           mean.overlap.dt2[, rank_tr := 1:.N, by = "subject"]
##           names(mean.overlap.dt2) = c("subject", "query", "gene_name", "mean_overlap", "rank_tr")
##           mean.overlap.dt2 = rbind(mean.overlap.dt2, mean.overlap.dt)
##           mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.min(rank_tr)], by = query]
          
##           ## mean.overlap.dt = mean.overlap.dt[order(-mean_overlap)]
##           ## mean.overlap.dt[, rank_tr := 1:.N, by = "query"]
##           ## mean.overlap.dt = mean.overlap.dt[rank_tr <= 20,]
##         }
        
##         ## mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = subject]
##         ## mean.overlap.dt2[, N_gene := .N, by = "subject"]
##         mean.overlap.dt2[, N_gene := .N, by = "query"]
##         if(any(mean.overlap.dt2$N_gene > 1)) {
##           warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
##           select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
##           select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
##           select.dt[, N_gene := .N, by = "gene_name"]
##           if(any(select.dt$N_gene > 1)) {
##             warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
##             select.dt = select.dt[which(width == min(width)),]
##             rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
##           }
##           potent_tr = select.dt$transcript_id
##         }
##         if(exists("potent_tr")) {
##           if(length(potent_tr) > 0) {
##             potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
##           }
##         } else {
##           ## potential_transcripts = mean.overlap.dt2$query
##           potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
##         }
##           ## (md.sub.gr@ranges %>% as.data.table())[end < start,]
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = TRUE, exact = TRUE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
          
##           ## gr.val(md.sub.gr, unlist(potential_transcripts.grl), val = names(unlist(potential_transcripts.grl)))
##           ## (md.sub.gr %O% potential_transcripts.grl) %>% head
##           ## gr.val(using val = names(values(y)))
##       }
##       ##   else if(!annotate_mismatch_type %in% c("longest","level")) {
##       ##   stop("annotate_mismatch_type must be either longest or level")
##       ## }
##       if(is.null("potential_transcript_merge")) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##       }
##       if(!is.null("potential_transcript_merge")) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$subject,.(gene_name,transcript_id)] %>% unique
##         add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "subject", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
##         md.novel.gr$gene_name = NULL
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
##       }
##       if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
##         warning("Some transcript have multiple potential genes. They will not be returned. Use reannotate = TRUE to get better transcript predictions")
##       }
##       ## if(reannotate) {
##       ##   potential_genes = md.novel.dt2[is.na(transcript_id),]$gene_name %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       ##   potential_genes.gr = as.data.table(gtf.gr)[gene_name %in% potential_genes,]
##       ##   ## now get the transcript from potential transcripts
##       ##   potential.exons.gr = potential_genes.gr[transcript_id %in% potential_transcripts,][type %in% c("exon","CDS"),]
        
##       ## }
##       ## now annotation the exons overlaps
##       md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##       potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##       potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
##       ##start new
##       ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##       ## ## new gr.val to match to specific transcripts
##       ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type", by = "transcript_id")
##       ## md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       ## potential.transcripts = potential.gtf.labels$transcript_id %>% unique
##       ## ##match the transcript id with the gr.val output
##       ## md.novel.lst3 = lapply(potential.transcripts, function(tr) {
##       ##   md.sub = md.novel.dt3[transcript_id == tr,.(seqnames, start, end, width, strand, qname, type, transcript_id, gene_name, qid, get(paste0("coding_type.",tr)))]
##       ##   names(md.sub) = gsub("V11","coding_type", names(md.sub))
##       ##   return(md.sub)
##       ## })
##       ## md.novel.dt3 = rbindlist(md.novel.lst3)  
##       ## md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##       ##end new
##       ##start new new
##       ##match the transcript id with the gr.val output


## ################################################################################################################################################################################################################################################################
##       ##removing this section to do all together
##       ## potential.transcripts = potential.gtf.labels$transcript_id %>% unique
##       ## md.novel.lst3 = lapply(potential.transcripts, function(tr) {
##       ##   md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
##       ##   potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
##       ##   ##breaks isn't working properly-convert these to breakpoints
##       ##   ## potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end)]
##       ##   ## bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), start = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end))
##       ##   ## bps.dt[,end := start ]
##       ##   ## bps.gr = GRanges(bps.dt, seqlengths = hg_seqlengths())
##       ##   ##
##       ##   md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
##       ##   ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
##       ##   ## new gr.val to match to specific transcripts
##       ##   md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
##       ##   md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       ##   return(md.novel.dt3)
##       ## })
##       ## md.novel.dt3 = rbindlist(md.novel.lst3)  
##       ## md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
## ################################################################################################################################################################################################################################################################
##       ##end new new
##       ##end new gr.val
##       ##start old gr.val
##       ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##       ## md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       ## md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##       ##end old gr.val
##       ##md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% unlist(x),"start_codon",x))]
##       ##make start codon only start
##                                         #md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x,"start_codon",list(x)))]
##       ##make exon only exon
##       ##md.novel.dt3[,coding_type_simple := lapply(coding_type_simple, function(x) ifelse("exon" %in% x,"exon",list(x)))]
##       ## md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##       ##                                                                           ifelse("UTR" %in% x, "UTR",
##       ##                                                                           ifelse("exon" %in% x, "exon","intron"))))]
##       ## md.novel.dt3 = md.novel.dt3[type != "N",]
##       ## md.novel.dt3[type == "X",coding_type_simple := "del"]
##       ## message("done annotating novel transcripts")
##     }
##     ## ###########################################################################################################################################################################################################################################

##     ## done with novel, now do ones with matched transcripts
##     ## md.full = md.dt3[structural_category == "full-splice_match",]
##     ## md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##     ## add the columns in md.novel.dt2 not in md.full or md.sub
##     if(exists("md.novel.dt2") & !reannotate) {    
##     gtf.possible.tr.dt = as.data.table(gtf.gr)[transcript_id %in% md.full$associated_transcript | transcript_id %in% md.sub$associated_transcript, .(gene_name, transcript_id)] %>% unique      
##     md.annotated = rbind(md.full,md.sub)
##     md.annotated = merge.data.table(md.annotated, gtf.possible.tr.dt, by.x = "associated_transcript",by.y = "transcript_id",all.x = TRUE)
##     md.annotated[, transcript_id := associated_transcript]
##       md.annotated = rbind(md.novel.dt2,md.annotated)
##     } else if (reannotate) {
##       md.annotated = md.novel.dt2
##     }
##     ## if(!reannotate) {
##       ## md.annotated = rbind(md.full,md.sub)
##       ## md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##     md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths())
##     potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##     potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##     potential.gtf.dt = as.data.table(potential.gtf.gr)
    
##     add.potential.ts = potential.gtf.dt[transcript_id %in% md.annotated$transcript_id,.(gene_name,transcript_id)] %>% unique
##     ## add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts$transcript_id,.(gene_name,transcript_id)] %>% unique
##     ## md.annotated.dt2 = merge.data.table(as.data.table(md.annotated.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##     md.annotated.dt2 = md.annotated
##     md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##     ## now annotation the exons overlaps
##     md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##     potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##     potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##     ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
## #########
##     ## potential.transcripts = md.annotated$transcript_id %>% unique      
##     ## md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
##     ##   md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
##     ##   potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
##     ##   ##breaks isn't working properly-convert these to breakpoints
##     ##   potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end)]
##     ##   bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), start = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end))
##     ##   bps.dt[,end := start ]
##     ##   bps.gr = GRanges(bps.dt, seqlengths = hg_seqlengths())
    
##     ##   ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
##     ##   md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
##     ##   ## new gr.val to match to specific transcripts
##     ##   md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
##     ##   md.annotated.dt3 = as.data.table(md_potential_tr.gr2)
##     ##   return(md.annotated.dt3)
##     ## })
##     potential.transcripts = potential.gtf.labels$transcript_id %>% unique
##     md.annotated.lst3 = lapply(potential.transcripts, function(tr) {
##       md_potential_tr.sub.gr = md_potential_tr.gr %Q% (transcript_id == tr)
##       potential.gtf.labels.sub = potential.gtf.labels %Q% (transcript_id == tr)
##       ##breaks isn't working properly-convert these to breakpoints
##       potential.gtf.labels.sub.dt = as.data.table(potential.gtf.labels.sub)[,.(seqnames,start,end)]
##       bps.dt = data.table(seqnames = c(potential.gtf.labels.sub.dt$seqnames, potential.gtf.labels.sub.dt$seqnames), end = c(potential.gtf.labels.sub.dt$start,potential.gtf.labels.sub.dt$end)) %>% unique
##       bps.dt[,start := end -1]
##       bps.gr = GRanges(bps.dt[order(start)], seqlengths = hg_seqlengths()) %>% gr.chr
      
##       ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels.sub,md_potential_tr.sub.gr)
##       md_potential_tr.gr2 = gUtils::gr.breaks(bp = bps.gr,md_potential_tr.sub.gr)
##       ## new gr.val to match to specific transcripts
##       md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels.sub,"coding_type")
##       md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       return(md.novel.dt3)
##     })
##     md.annotated.dt3 = rbindlist(md.annotated.lst3)  
##     md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
## ################
##     ## md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##     ## md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##     ## md.annotated.dt3 = as.data.table(md_potential_tr.gr2)
##     md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##     md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                                   ifelse("UTR" %in% x, "UTR",
##                                                                                   ifelse("exon" %in% x, "exon","intron"))))]
##     ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
##     ## md.annotated.dt3[type == "X",coding_type_simple := "intron"]
##     md.annotated.dt3 = md.annotated.dt3[type != "N",]
##     md.annotated.dt3[type == "X",coding_type_simple := "del"]
##     md.all.dt = md.annotated.dt3
##     message("finished merging together")
##     ## rbind all together
##     ## md.all.dt = rbind(md.annotated.dt3,md.novel.dt3, fill = TRUE)
##                                         #    } #else if(reannotate) {
##                                         #md.all.dt = md.novel.dt3
##                                         #}
##     ## sort based on strand to have walks align
##     pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
##     neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
##     md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
##     if(add_gencode_name) {
##       gencode.dt = as.data.table(gencode.gr)[transcript_id %in% md.dt4$transcript_id,]
##       gencode.sub.dt = gencode.dt[,.(transcript_id,transcript_type,transcript_name)] %>% unique
##       md.dt4 = merge.data.table(md.dt4,gencode.sub.dt, by = "transcript_id", all.x = TRUE)
##     }
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     ## convert to grl and gw
##     message("converting to grl")
##     md.grl2 = rtracklayer::split(md.gr2, f = mcols(md.gr2)["qname"])
##     message("done coverting to grl")
##     ## md.novel[qname == "transcript/583547",]    
##   } else {
##     md.dt3
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     return(md.grl2)
##   }
##   ## sort positive strand forward and negative strand reverse
##   if(type == "gw") {
##     message("type is 'gw' converting grl to gW object")
##     md.gw = gW(grl = md.grl2)
##     ## add coloring of nodes to match track.gencode
##     cmap.dt = data.table(category = c("exon", "intron", "start_codon", "stop_codon", "UTR", "del"), color = c("#0000FF99", "#FF0000B3", "green", "red", "#A020F066", "orange"))
##     for(x in 1:nrow(cmap.dt)) {
##       cmap.sub = cmap.dt[x,]
##       md.gw$nodes[coding_type_simple == cmap.sub$category]$mark(col = cmap.sub$color)
##     }
##     transcript_names = sapply(md.gw$grl, function(x) x$transcript_name[1])
##     md.gw$set(name = transcript_names)
##     message("done converting grl to gW object. returning gW")
##     return(md.gw)
##   } else if(type == "grl") {
##     message("type is 'grl' returning grl")
##     return(md.grl2)
##   }
## }

## ## function to generate gwalks or grls from isoseq bams and a specified gr
## get_iso_reads = function(bam, gr, gtf, collapsed_group = NULL, collapsed_class = NULL, type = "gw", annotate_mismatch = TRUE, annotate_mismatch_type = "percent",reann_with = "both", reannotate = FALSE) {
##   message(paste0("Reading in ",bam))
##   md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
##   message("Done reading")
##   md.dt = as.data.table(md.gr.bam)
##   ## pos.dt = md.dt[strand == "+",][order(seqnames,start,end),]
##   ## neg.dt = md.dt[strand == "-",][order(seqnames,-start,-end),]
##   ## md.dt = rbind(pos.dt,neg.dt)
##   md.dt = md.dt[width != 0,]
##   md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##   message("Splicing cigar string")
##   md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
##   message("Done Splicing cigar string")
##   gencode.gr = gencode@data[[1]] %>% unlist
##   gencode.gr$type2 = gencode.gr$type
##   if(exists("potential_transcript_merge")) {
##     potential_transcript_merge = NULL
##   }
##   if(!reannotate) {
##     ##look for collapsed_group and collapsed_class if missing
##     if(is.null(collapsed_group) | is.null(collapsed_class)) {
##       if(is.null(collapsed_group) && is.null(collapsed_class)) {
##         message("looking for collapsed_group and collapsed class")
##         folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
##         collapsed_group = paste0(folder_find,"collapsed.group.txt")
##         collapsed_class = paste0(folder_find,"collapsed_classification.txt")
##         if(all(file.exists(c(collapsed_group,collapsed_class)))) {
##           message("Found! collapsed.group.txt and collapsed_classification.txt")
##         }
##       }
##     }
##     group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
##     class.dt = fread(collapsed_class)
##     ## get labels for already identified transcripts
##     group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
##     group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
##     group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
##     expanded.dt = copy(group.sub.dt)
##     expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
##     group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
##     ##add transcript labels to bam reads
##     md.dt2 = unlist(md.grl) %>% as.data.table()
##     md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")
##   } else if(reannotate) {
##     ## md.dt3 = md.dt[,.(qname,seqnames,start,end,width,strand,type)]
##     md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
##   }
##   if(annotate_mismatch) {
##     if(!reannotate) {
##       md.full = md.dt3[structural_category == "full-splice_match",]
##       md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##       md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
##     } else if(reannotate) {
##       md.novel = md.dt3
##       md.full = NULL
##     }
## ###################################
##     ## first  label ones without a transcript annotation
##     md.novel.gr = GRanges(md.novel)
##     if(class(gtf) == "character") {
##       message(paste0("Reading in ",gtf)) 
##       gtf = rtracklayer::readGFF(gtf) %>% as.data.table
##       message(paste0("Done reading in ",gtf))
##     }
##     if(class(gtf)[1] != "data.table" & class(gtf)[1] != "GRanges") {
##       stop("Gtf must be a character to a gtf file or a data.table")
##     }
##     if(any(class(gtf) != "GRanges")) {
##       gtf.gr = GRanges(gtf, seqlengths = hg_seqlengths())
##     } else {
##       gtf.gr = gtf
##     }
##     ##add potential genes to md.novel
##     if(length(md.novel.gr) > 0) {
##       md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##       potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##       potential.gtf.dt = as.data.table(potential.gtf.gr)
##       ## annotate to the longest transcript??
##       if(annotate_mismatch_type == "longest") {
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "level") {
##         ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "match") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         ##potential one
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
##         ## get the transcript with the most N by transcript
##         potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "percent") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
##         potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
##         potential_transcripts.grl = split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
##         ## get maximum overlap of each transcript with each potential transcript
##         md.novel.grl = split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
##         message("getting overlap of transcripts with potential transcripts")
##         mean.overlap.dt = lapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             ## annotate with utr and exons or just exons
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             ##try returning the mean
##             return(data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt)
##         }) %>% do.call(rbind,.)
##         mean.overlap.dt = mean.overlap.dt[, .SD[which.max(mean_overlap)], by = query]
##         message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
##         mean.overlap.dt2 = lapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt2 = lapply(mean.overlap.dt$subject, function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             if(reann_with == "both") {
##               potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             } else if (reann_with == "exon") {
##               potential.sub.gr = potential.sub.gr[!(potential.sub.gr$type %in% c("transcript","UTR"))]
##             }
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
##             ##try returning the mean
##             return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt2)
##         }) %>% do.call(rbind,.)
##         ## mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = query] # used to be selecting for best transcript but I actually want to assign the best transcript to each read
##         mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = subject]
##         mean.overlap.dt2[, N_gene := .N, by = "subject"]
##         if(any(mean.overlap.dt2$N_gene > 1)) {
##           warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
##           select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
##           select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
##           select.dt[, N_gene := .N, by = "gene_name"]
##           if(any(select.dt$N_gene > 1)) {
##             warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
##             select.dt = select.dt[which(width == min(width)),]
##             rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
##           }
##           potent_tr = select.dt$transcript_id
##         }
##         if(exists("potent_tr")) {
##           if(length(potent_tr) > 0) {
##             potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene == 1,]$query)
##           }
##         } else {
##           ## potential_transcripts = mean.overlap.dt2$query
##           potential_transcripts_merge = mean.overlap.dt2[,.(subject,query)] %>% unique()
##         }
##           ## (md.sub.gr@ranges %>% as.data.table())[end < start,]
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = TRUE, exact = TRUE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
          
##           ## gr.val(md.sub.gr, unlist(potential_transcripts.grl), val = names(unlist(potential_transcripts.grl)))
##           ## (md.sub.gr %O% potential_transcripts.grl) %>% head
##           ## gr.val(using val = names(values(y)))
##       }
##       ##   else if(!annotate_mismatch_type %in% c("longest","level")) {
##       ##   stop("annotate_mismatch_type must be either longest or level")
##       ## }
##       if(is.null("potential_transcript_merge")) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##       }
##       if(!is.null("potential_transcript_merge")) {
##         add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts_merge$query,.(gene_name,transcript_id)] %>% unique
##         add.potential.ts = merge.data.table(potential_transcripts_merge, add.potential.ts, by.x = "query", by.y = "transcript_id", all.x = TRUE) %>% setnames(., c("transcript_id","qname","gene_name"))
##         md.novel.gr$gene_name = NULL
##         md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "qname", all.x = TRUE) #add gene name & most likely transcript to novel
##         md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]        
##       }
##       if(nrow(md.novel.dt2[transcript_id == "multiple_potential_genes",]) > 0) {
##         warning("Some transcript have multiple potential genes. They will not be returned")
##         ## warning("Some transcript have multiple potential genes. If reannotate == TRUE, then we will select the gene with the most exon overlap. If false, they will not be returned")
##       }
##       ## if(reannotate) {
##       ##   potential_genes = md.novel.dt2[is.na(transcript_id),]$gene_name %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       ##   potential_genes.gr = as.data.table(gtf.gr)[gene_name %in% potential_genes,]
##       ##   ## now get the transcript from potential transcripts
##       ##   potential.exons.gr = potential_genes.gr[transcript_id %in% potential_transcripts,][type %in% c("exon","CDS"),]
        
##       ## }
##       ## now annotation the exons overlaps
##       md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##       potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##       potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
##       md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##       md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##       md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##       ##md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% unlist(x),"start_codon",x))]
##       ##make start codon only start
##                                         #md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x,"start_codon",list(x)))]
##       ##make exon only exon
##       ##md.novel.dt3[,coding_type_simple := lapply(coding_type_simple, function(x) ifelse("exon" %in% x,"exon",list(x)))]
##       md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                                 ifelse("UTR" %in% x, "UTR",
##                                                                                 ifelse("exon" %in% x, "exon","intron"))))]
##       md.novel.dt3 = md.novel.dt3[type != "N",]
##       md.novel.dt3[type == "X",coding_type_simple := "del"]
##       message("done annotating novel transcripts")
##     }
##     ## ###########################################################################################################################################################################################################################################

##     ## done with novel, now do ones with matched transcripts
##     ## md.full = md.dt3[structural_category == "full-splice_match",]
##     ## md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##     if(!reannotate) {
##       md.annotated = rbind(md.full,md.sub)
##       md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##       potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##       potential.gtf.dt = as.data.table(potential.gtf.gr)
      
##       add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
##       md.annotated.dt2 = merge.data.table(as.data.table(md.annotated.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##       md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##       ## now annotation the exons overlaps
##       md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##       potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##       potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
##       md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##       md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##       md.annotated.dt3 = as.data.table(md_potential_tr.gr2)
##       md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##       md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                                     ifelse("UTR" %in% x, "UTR",
##                                                                                     ifelse("exon" %in% x, "exon","intron"))))]
##       ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
##       ## md.annotated.dt3[type == "X",coding_type_simple := "intron"]
##       md.annotated.dt3 = md.annotated.dt3[type != "N",]
##       md.annotated.dt3[type == "X",coding_type_simple := "del"]
##       message("finished merging together")
##       ## rbind all together
##       md.all.dt = rbind(md.annotated.dt3,md.novel.dt3)
##     } else if(reannotate) {
##       md.all.dt = md.novel.dt3
##     }
##     ## sort based on strand to have walks align
##     pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
##     neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
##     md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     ## convert to grl and gw
##     message("converting to grl")
##     md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     message("done coverting to grl")
##     ## md.novel[qname == "transcript/583547",]    
##   } else {
##     md.dt3
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     return(md.grl2)
##   }
##   ## sort positive strand forward and negative strand reverse
##   if(type == "gw") {
##     message("type is 'gw' converting grl to gW object")
##     md.gw = gW(grl = md.grl2)
##     message("done converting grl to gW object. returning gW")
##     return(md.gw)
##   } else if(type == "grl") {
##     message("type is 'grl' returning grl")
##     return(md.grl2)
##   }
## }

## ## old version
## get_iso_reads = function(bam, gr, gtf, collapsed_group = NULL, collapsed_class = NULL, type = "gw", annotate_mismatch = TRUE, annotate_mismatch_type = "level") {
##   message(paste0("Reading in ",bam))
##   md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
##   message("Done reading")
##   md.dt = as.data.table(md.gr.bam)
##   ## pos.dt = md.dt[strand == "+",][order(seqnames,start,end),]
##   ## neg.dt = md.dt[strand == "-",][order(seqnames,-start,-end),]
##   ## md.dt = rbind(pos.dt,neg.dt)
##   md.dt = md.dt[width != 0,]
##   md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
##   message("Splicing cigar string")
##   md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
##   message("Done Splicing cigar string")
##   gencode.gr = gencode@data[[1]] %>% unlist
##   gencode.gr$type2 = gencode.gr$type
##   ##look for collapsed_group and collapsed_class if missing
##   if(is.null(collapsed_group) | is.null(collapsed_class)) {
##     if(is.null(collapsed_group) && is.null(collapsed_class)) {
##       message("looking for collapsed_group and collapsed class")
##       folder_find = tstrsplit(bam,"/") %>% .[1:(length(.)-1)] %>% paste0(., collapse = "/") %>% paste0(.,"/")
##       collapsed_group = paste0(folder_find,"collapsed.group.txt")
##       collapsed_class = paste0(folder_find,"collapsed_classification.txt")
##       if(all(file.exists(c(collapsed_group,collapsed_class)))) {
##         message("Found! collapsed.group.txt and collapsed_classification.txt")
##       }
##     }
##   }
##   group.dt = fread(collapsed_group, col.names = c("isoform","transcripts"))
##   class.dt = fread(collapsed_class)
##   ## get labels for already identified transcripts
##   group.dt = merge.data.table(group.dt,class.dt, by = "isoform", all.x = TRUE)
##   group.sub.dt = group.dt[,.(isoform,transcripts,structural_category,subcategory,associated_transcript)]
##   group.sub.dt[, transcripts := strsplit(transcripts, ",")] #split transcripts into unique rows
##   expanded.dt = copy(group.sub.dt)
##   expanded.dt2 = expanded.dt[, .(transcript = unlist(transcripts)), by = isoform]
##   group.sub.dt2 = merge.data.table(group.sub.dt[, !"transcripts"], expanded.dt2, by = "isoform")  
##   ##add transcript labels to bam reads
##   md.dt2 = unlist(md.grl) %>% as.data.table()
##   md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")
##   if(annotate_mismatch) {
##     md.full = md.dt3[structural_category == "full-splice_match",]
##     md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##     md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
## ###################################
##     ## first  label ones without a transcript annotation
##     md.novel.gr = GRanges(md.novel)
##     if(class(gtf) == "character") {
##       message(paste0("Reading in ",gtf)) 
##       gtf = rtracklayer::readGFF(gtf) %>% as.data.table
##       message(paste0("Done reading in ",gtf))
##     }
##     if(class(gtf)[1] != "data.table" & class(gtf)[1] != "GRanges") {
##       stop("Gtf must be a character to a gtf file or a data.table")
##     }
##     if(any(class(gtf) != "GRanges")) {
##       gtf.gr = GRanges(gtf, seqlengths = hg_seqlengths())
##     } else {
##       gtf.gr = gtf
##     }
##     ##add potential genes to md.novel
##     if(length(md.novel.gr) > 0) {
##       md.novel.gr = GRanges(md.novel,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##       potential_genes = md.novel.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##       potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##       potential.gtf.dt = as.data.table(potential.gtf.gr)
##       ## annotate to the longest transcript??
##       if(annotate_mismatch_type == "longest") {
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.max(width)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "level") {
##         ## annotate mismatch to the highest level transcript for that gene- if multiple pick the first
##         potential_transcripts = potential.gtf.dt[type == "transcript",][, .SD[which.min(level)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "match") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]        
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         ##potential one
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         potential_transcripts = as.data.table(potential_transcripts[gr.match(query = md.novel.gr, subject = potential_transcripts)])[,.(gene_name, transcript_id)] %>% table %>% t %>% as.data.table
##         ## get the transcript with the most N by transcript
##         potential_transcripts = potential_transcripts[, .SD[which.max(N)], by = gene_name]$transcript_id
##       } else if(annotate_mismatch_type == "percent") {
##         gtf.sub.gr = gtf.gr[gtf.gr$type == "transcript"]
##         gtf.sub.gr$perc_overlap = gtf.sub.gr %O% md.novel.gr
##         potential_transcripts = gtf.sub.gr %Q% (perc_overlap != 0)
##         ## get transcripts into grl rather than getting the max overlap with transcripts in annotate_mismatch_type = "match"
##         potential_transcripts.gr2 = gtf.gr[gtf.gr$transcript_id %in% unique(potential_transcripts$transcript_id) & type != "transcript"]
##         potential_transcripts.grl = split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
##         ## get maximum overlap of each transcript with each potential transcript
##         md.novel.grl = split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
##         message("getting overlap of transcripts with potential transcripts")
##         mean.overlap.dt = lapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             ##try returning the mean
##             return(data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt)
##         }) %>% do.call(rbind,.)
##         mean.overlap.dt = mean.overlap.dt[, .SD[which.max(mean_overlap)], by = query]
##         message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
##         mean.overlap.dt2 = lapply(names(md.novel.grl), function(tr) {
##           md.sub.gr = md.novel.grl[[tr]]
##           mean.overlap.dt2 = lapply(mean.overlap.dt$subject, function(pot.tr) {
##             potential.sub.gr = potential_transcripts.grl[[pot.tr]]
##             potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
##             ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
##             potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
##             ##try returning the mean
##             return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
##           }) %>% do.call(rbind,.)
##           return(mean.overlap.dt2)
##         }) %>% do.call(rbind,.)
##         mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = query]
##         mean.overlap.dt2[, N_gene := .N, by = "gene_name"]
##         if(any(mean.overlap.dt2$N_gene > 1)) {
##           warning(paste0("multiple potential transcripts found for "),unique(mean.overlap.dt2[N_gene > 1,]$gene_name), "will try to assign one based on the lowest support level in gtf (most support)")        
##           select.dt = potential.gtf.dt[type == "transcript",][transcript_id %in% (mean.overlap.dt2[N_gene > 1,]$query),]
##           select.dt = select.dt[which(transcript_support_level == min(transcript_support_level)),] #get the entries with the minimum
##           select.dt[, N_gene := .N, by = "gene_name"]
##           if(any(select.dt$N_gene > 1)) {
##             warning(paste0("Multiple transcripts still found, selecting the shorter one for ",unique(select.dt[N_gene > 1,]$gene_name)))
##             select.dt = select.dt[which(width == min(width)),]
##             rbind(select.dt,mean.overlap.dt2[N_gene  == 1,],fill = TRUE)
##           }
##           potent_tr = select.dt$transcript_id
##         }
##         if(length(potent_tr) > 0) {
##           potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene ==1,]$query)
##         } else {
##           potential_transcripts = mean.overlap.dt2$query
##         }
##           ## (md.sub.gr@ranges %>% as.data.table())[end < start,]
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = TRUE, exact = TRUE, maxgap = 20)
##           ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
          
##           ## gr.val(md.sub.gr, unlist(potential_transcripts.grl), val = names(unlist(potential_transcripts.grl)))
##           ## (md.sub.gr %O% potential_transcripts.grl) %>% head
##           ## gr.val(using val = names(values(y)))
##       }
      
##       ##   else if(!annotate_mismatch_type %in% c("longest","level")) {
##       ##   stop("annotate_mismatch_type must be either longest or level")
##       ## }
##       add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
##       md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##       md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##       ## now annotation the exons overlaps
##       md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##       potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##       potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##       ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
##       md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##       md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##       md.novel.dt3 = as.data.table(md_potential_tr.gr2)
##       md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##       ##md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% unlist(x),"start_codon",x))]
##       ##make start codon only start
##                                         #md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x,"start_codon",list(x)))]
##       ##make exon only exon
##       ##md.novel.dt3[,coding_type_simple := lapply(coding_type_simple, function(x) ifelse("exon" %in% x,"exon",list(x)))]
##       md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                                 ifelse("UTR" %in% x, "UTR",
##                                                                                 ifelse("exon" %in% x, "exon","intron"))))]
##       md.novel.dt3 = md.novel.dt3[type != "N",]
##       md.novel.dt3[type == "X",coding_type_simple := "del"]
##       message("done annotating novel transcripts")
##     }
##     ## ###########################################################################################################################################################################################################################################

##     ## done with novel, now do ones with matched transcripts
##     ## md.full = md.dt3[structural_category == "full-splice_match",]
##     ## md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
##     md.annotated = rbind(md.full,md.sub)
##     md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
##     potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
##     potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
##     potential.gtf.dt = as.data.table(potential.gtf.gr)
    
##     add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
##     md.annotated.dt2 = merge.data.table(as.data.table(md.annotated.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
##     md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
##     ## now annotation the exons overlaps
##     md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
##     potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
##     potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
##     ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
##     md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
##     md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
##     md.annotated.dt3 = as.data.table(md_potential_tr.gr2)
##     md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
##     md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
##                                                                               ifelse("UTR" %in% x, "UTR",
##                                                                               ifelse("exon" %in% x, "exon","intron"))))]
##     ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
##     ## md.annotated.dt3[type == "X",coding_type_simple := "intron"]
##     md.annotated.dt3 = md.annotated.dt3[type != "N",]
##     md.annotated.dt3[type == "X",coding_type_simple := "del"]
##     message("finished merging together")
##     ## rbind all together
##     md.all.dt = rbind(md.annotated.dt3,md.novel.dt3)
##     ## sort based on strand to have walks align
##     pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
##     neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
##     md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     ## convert to grl and gw
##     message("converting to grl")
##     md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     message("done coverting to grl")
##     ## md.novel[qname == "transcript/583547",]    
##   } else {
##     md.dt3
##     md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
##     md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
##     return(md.grl2)
##   }
##   ## sort positive strand forward and negative strand reverse
##   if(type == "gw") {
##     message("type is 'gw' converting grl to gW object")
##     md.gw = gW(grl = md.grl2)
##     message("done converting grl to gW object. returning gW")
##     return(md.gw)
##   } else if(type == "grl") {
##     message("type is 'grl' returning grl")
##     return(md.grl2)
##   }
## }


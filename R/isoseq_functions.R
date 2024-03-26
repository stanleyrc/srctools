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
                         cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
  message("Splicing cigar string")
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  message("Done Splicing cigar string")
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
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
    md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
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
        meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
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
                         annotate_missing = TRUE,
                         cores = 1) {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
  message("Splicing cigar string")
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  message("Done Splicing cigar string")
  ## gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr = gtf
  gencode.gr$type2 = gencode.gr$type
  ## if(exists("potential_transcript_merge")) {
  potential_transcript_merge = NULL
  ## }
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
      md.dt3 = as.data.table(unlist(md.grl))[,.(qname,seqnames,start,end,width,strand,type)][type != "N",]
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
        browser()
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
        meta.add.dt = meta.add.dt[,.(transcript_id, qname, isoform, structural_category, subcategory, associated_transcript, gene_name, reannotated, transcript_type, transcript_name, qr_tr)]
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



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

## function to generate gwalks or grls from isoseq bams and a specified gr

get_iso_reads = function(bam, gr, gtf, collapsed_group = NULL, collapsed_class = NULL, type = "gw", annotate_mismatch = TRUE, annotate_mismatch_type = "level") {
  message(paste0("Reading in ",bam))
  md.gr.bam = bamUtils::read.bam(bam, gr, stripstrand = FALSE, verbose = TRUE, pairs.grl.split = FALSE)
  message("Done reading")
  md.dt = as.data.table(md.gr.bam)
  ## pos.dt = md.dt[strand == "+",][order(seqnames,start,end),]
  ## neg.dt = md.dt[strand == "-",][order(seqnames,-start,-end),]
  ## md.dt = rbind(pos.dt,neg.dt)
  md.dt = md.dt[width != 0,]
  md.gr = GRanges(md.dt,seqlengths = hg_seqlengths())
  message("Splicing cigar string")
  md.grl = bamUtils::splice.cigar(md.gr,get.seq = TRUE, rem.soft = FALSE)
  message("Done Splicing cigar string")
  gencode.gr = gencode@data[[1]] %>% unlist
  gencode.gr$type2 = gencode.gr$type
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
  md.dt3 = merge.data.table(md.dt2, group.sub.dt2, by.x = "qname", by.y = "transcript")
  if(annotate_mismatch) {
    md.full = md.dt3[structural_category == "full-splice_match",]
    md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
    md.novel = md.dt3[structural_category != "full-splice_match" & associated_transcript == "novel"]
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
        potential_transcripts.grl = split(potential_transcripts.gr2, f = mcols(potential_transcripts.gr2)["transcript_id"])
        ## get maximum overlap of each transcript with each potential transcript
        md.novel.grl = split(md.novel.gr, f = mcols(md.novel.gr)["qname"])
        message("getting overlap of transcripts with potential transcripts")
        mean.overlap.dt = lapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt = lapply(names(potential_transcripts.grl), function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            ##try returning the mean
            return(data.table(query = tr, subject = pot.tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(md.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt)
        }) %>% do.call(rbind,.)
        mean.overlap.dt = mean.overlap.dt[, .SD[which.max(mean_overlap)], by = query]
        message("Now getting overlap of potential transcripts with the transcripts to select the best transcript")
        mean.overlap.dt2 = lapply(names(md.novel.grl), function(tr) {
          md.sub.gr = md.novel.grl[[tr]]
          mean.overlap.dt2 = lapply(mean.overlap.dt$subject, function(pot.tr) {
            potential.sub.gr = potential_transcripts.grl[[pot.tr]]
            potential.sub.gr = potential.sub.gr[potential.sub.gr$type != "transcript"]
            ## md.sub.gr$percent = md.sub.gr %O% potential.sub.gr
            potential.sub.gr$percent = potential.sub.gr %O% md.sub.gr
            ##try returning the mean
            return(data.table(query = pot.tr, subject = tr, gene_name = unique(potential.sub.gr$gene_name),mean_overlap = mean(potential.sub.gr$percent, na.rm = TRUE)))
          }) %>% do.call(rbind,.)
          return(mean.overlap.dt2)
        }) %>% do.call(rbind,.)
        mean.overlap.dt2 = mean.overlap.dt2[, .SD[which.max(mean_overlap)], by = query]
        mean.overlap.dt2[, N_gene := .N, by = "gene_name"]
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
        if(length(potent_tr) > 0) {
          potential_transcripts = c(potent_tr, mean.overlap.dt2[N_gene ==1,]$query)
        } else {
          potential_transcripts = mean.overlap.dt2$query
        }
          ## (md.sub.gr@ranges %>% as.data.table())[end < start,]
          ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
          ## grl.in(potential_transcripts.grl, md.sub.gr, logical = TRUE, exact = TRUE, maxgap = 20)
          ## grl.in(potential_transcripts.grl, md.sub.gr, logical = FALSE, maxgap = 20)
          
          ## gr.val(md.sub.gr, unlist(potential_transcripts.grl), val = names(unlist(potential_transcripts.grl)))
          ## (md.sub.gr %O% potential_transcripts.grl) %>% head
          ## gr.val(using val = names(values(y)))
      }
      
      ##   else if(!annotate_mismatch_type %in% c("longest","level")) {
      ##   stop("annotate_mismatch_type must be either longest or level")
      ## }
      add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
      md.novel.dt2 = merge.data.table(as.data.table(md.novel.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
      md.novel.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
      ## now annotation the exons overlaps
      md_potential_tr.gr = GRanges(md.novel.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
      potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
      potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
      ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
      md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
      md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
      md.novel.dt3 = as.data.table(md_potential_tr.gr2)
      md.novel.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
      ##md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% unlist(x),"start_codon",x))]
      ##make start codon only start
                                        #md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x,"start_codon",list(x)))]
      ##make exon only exon
      ##md.novel.dt3[,coding_type_simple := lapply(coding_type_simple, function(x) ifelse("exon" %in% x,"exon",list(x)))]
      md.novel.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                                ifelse("UTR" %in% x, "UTR",
                                                                                ifelse("exon" %in% x, "exon","intron"))))]
      md.novel.dt3 = md.novel.dt3[type != "N",]
      md.novel.dt3[type == "X",coding_type_simple := "del"]
      message("done annotating novel transcripts")
    }
    ## ###########################################################################################################################################################################################################################################

    ## done with novel, now do ones with matched transcripts
    ## md.full = md.dt3[structural_category == "full-splice_match",]
    ## md.sub = md.dt3[structural_category != "full-splice_match" & associated_transcript != "novel"]
    md.annotated = rbind(md.full,md.sub)
    md.annotated.gr = GRanges(md.annotated,seqlengths = hg_seqlengths()) %>% gUtils::gr.val(.,gtf.gr, "gene_name")
    potential_genes = md.annotated.gr$gene_name  %>% unique %>% strsplit(., ",") %>% unlist %>% gsub(" ","",.) %>% unique
    potential.gtf.gr = gtf.gr %Q% (gene_name %in% potential_genes)
    potential.gtf.dt = as.data.table(potential.gtf.gr)
    
    add.potential.ts = potential.gtf.dt[transcript_id %in% potential_transcripts,.(gene_name,transcript_id)] %>% unique
    md.annotated.dt2 = merge.data.table(as.data.table(md.annotated.gr),add.potential.ts, by = "gene_name", all.x = TRUE) #add gene name to novel
    md.annotated.dt2[is.na(transcript_id), transcript_id := "multiple_potential_genes"]
    ## now annotation the exons overlaps
    md_potential_tr.gr = GRanges(md.annotated.dt2[transcript_id != "multiple_potential_genes",], seqlengths = hg_seqlengths()) ## subset to ones attempting to annotate
    potential.gtf.labels = potential.gtf.dt[transcript_id %in% unique(md_potential_tr.gr$transcript_id) & type != "transcript",] %>% GRanges(.,seqlengths = hg_seqlengths())
    potential.gtf.labels$coding_type = as.character(potential.gtf.labels$type)
    ## md_potential_tr.gr = md_potential_tr.gr %Q% (qname == "transcript/583547")
    md_potential_tr.gr2 = gUtils::gr.breaks(bp = potential.gtf.labels,md_potential_tr.gr)
    md_potential_tr.gr2 = gr.val(md_potential_tr.gr2,potential.gtf.labels,"coding_type")
    md.annotated.dt3 = as.data.table(md_potential_tr.gr2)
    md.annotated.dt3[, coding_type_vect := lapply(coding_type, function(x) unlist(strsplit(x, ", ")))]
    md.annotated.dt3[, coding_type_simple := lapply(coding_type_vect, function(x) ifelse("start_codon" %in% x, "start_codon",
                                                                              ifelse("UTR" %in% x, "UTR",
                                                                              ifelse("exon" %in% x, "exon","intron"))))]
    ## md.annotated.dt3 = md.annotated.dt3[type != "N",]
    ## md.annotated.dt3[type == "X",coding_type_simple := "intron"]
    md.annotated.dt3 = md.annotated.dt3[type != "N",]
    md.annotated.dt3[type == "X",coding_type_simple := "del"]
    message("finished merging together")
    ## rbind all together
    md.all.dt = rbind(md.annotated.dt3,md.novel.dt3)
    ## sort based on strand to have walks align
    pos.dt = md.all.dt[strand == "+",][order(seqnames,start,end),]
    neg.dt = md.all.dt[strand == "-",][order(seqnames,-start,-end),]
    md.dt4 = rbind(pos.dt,neg.dt)#, test.tr, fill = TRUE)
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    ## convert to grl and gw
    message("converting to grl")
    md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
    message("done coverting to grl")
    ## md.novel[qname == "transcript/583547",]    
  } else {
    md.dt3
    md.gr2 = GRanges(md.dt4, seqlengths = hg_seqlengths())
    md.grl2 = split(md.gr2, f = mcols(md.gr2)["qname"])
    return(md.grl2)
  }
  ## sort positive strand forward and negative strand reverse
  if(type == "gw") {
    message("type is 'gw' converting grl to gW object")
    md.gw = gW(grl = md.grl2)
    message("done converting grl to gW object. returning gW")
    return(md.gw)
  } else if(type == "grl") {
    message("type is 'grl' returning grl")
    return(md.grl2)
  }
}


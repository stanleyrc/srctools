


getPGV = function(web = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/", pgv_sub_folder = "pgv/", json_file = paste0(web,pgv_sub_folder, "/public/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/public/data/"), settings = paste0(web,pgv_sub_folder, "/public/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}


getsharedPGV = function(web = "/gpfs/commons/home/sclarke/lab/pgv_content/", pgv_sub_folder = "", json_file = paste0(web,pgv_sub_folder, "/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/data/"), settings = paste0(web,pgv_sub_folder, "/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}

#create json from data.table
dt2json = function(dt,patient.id,ref,settings,file_name = paste(getwd(),"test.json",sep="/")) {
    #create vector of seqlengths
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[,.(chromosome,startPoint,endPoint)]
    colnames(chrom_lengths) = c("seqnames","start","end")
    chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M
#convert to ggraph and create json
    gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.chr()
    if(any(gr1@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(gr1))]$end)) {
        stop(paste("the seqlengths of your granges has ranges that are not contained in the seqlengths of",ref))
    }
    jab = gG(nodes=gr1)
#    jab = gGnome::refresh(jab)
    settings_y = list(y_axis = list(title = "copy number",
                                  visible = TRUE))
    node.json = gr2dt(jab$nodes$gr[, "snode.id"])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, y = 1,title=snode.id,type="interval",strand="*")]
    gg.js = list(intervals = node.json, connections = data.table())
    gg.js = c(list(settings = settings_y), gg.js)
    jsonlite::write_json(gg.js, file_name,
                         pretty=TRUE, auto_unbox=TRUE, digits=4)
    return(file_name)
}






#get hic resolutions printed in scientific
                                        #return resolutions of hic as scientific
hic_res = function(hic) {
                                        #return as scientific
    old_options = options("scipen", "digits")
    options(scipen = -999, digits = 7)
                                        #get resolutions
    reses = strawr::readHicBpResolutions(hichip.pairs[pair, contact.matrix]) %>% sort() %>% signif(., digits = 5)
    on.exit(options(old_options))
    return(print(reses))
}


hunt = function(huntdown = NULL) {
    ps = ps() %>% data.table()
    ps_all <- ps[,.(total_memory_GB = round(sum(rss / (1024 * 1024), na.rm = TRUE) / 1024, 4)), by = username][order(-total_memory_GB)]
    if(!is.null(huntdown)) {
        ps_user = ps[grep(huntdown, username),.(username, pid, ppid, name, status, resident_mem_gb = round(rss / (1024^3), 2), virtual_mem_gb = round(rss / (1024^3), 2), created)][order(-resident_mem_gb)]
    }
    if(is.null(huntdown)) {
        return(ps_all)
    } else {
        return(ps_user)
    }
}

get_lambda = function(pvals) {
    pval.dt = data.table(pval = pvals)
    names(pval.dt) = "pval"
    pval.dt = pval.dt[order(pval),]
    pval.dt[, rank := rank(log(pval)) / .N]
    pval.dt[, log.pval := -log10(pval)]
    pval.dt[, log.rank := -log10(rank)]
    lambda = lm(log.pval ~ log.rank - 1, pval.dt)$coefficients
}

hg19_seq = function() {
    hg_seqlengths("~/DB/UCSC/hg19.chrom.sizes")
}
hg38_seq = function() {
    hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")
}

hg19_gr = function(chr=TRUE,keep = c(1:22,"X","Y")) {
    gr = si2gr(hg_seqlengths("~/DB/UCSC/hg19.chrom.sizes"))  %>% gr.nochr() %>% sortSeqlevels() %>% sort()
    gr = gr %Q% (seqnames %in% keep)
    if(chr==TRUE) {
        return(gr.chr(gr))
    } else {
        return(gr)
    }
}

hg38_gr = function(chr=TRUE,keep = c(1:22,"X","Y")) {
    gr = si2gr(hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes"))  %>% gr.nochr() %>% sortSeqlevels() %>% sort()
    gr = gr %Q% (seqnames %in% keep)
    if(chr==TRUE) {
        return(gr.chr(gr))
    } else {
        return(gr)
    }
}


## hg19_gr = function(chr=FALSE,keep = c(1:22,"X","Y")) {
##     seq_lengths = hg_seqlengths("~/DB/UCSC/hg19.chrom.sizes")
##     gr = GRanges(seqnames = names(seq_lengths), ranges = IRanges(start = 1, end = seq_lengths)) %>% gr.nochr() %>% sortSeqlevels() %>% sort()
##     gr = gr %Q% (seqnames %in% keep)
##     if(chr==TRUE) {
##         return(gr.chr(gr))
##     } else {
##         return(gr)
##     }
## }

## hg38_gr = function(chr=TRUE,keep = c(1:22,"X","Y")) {
##     seq_lengths = hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")
##     gr = GRanges(seqnames = names(seq_lengths), ranges = IRanges(start = 1, end = seq_lengths)) %>% gr.nochr() %>% sortSeqlevels() %>% sort()
##     gr = gr %Q% (seqnames %in% keep)
##     if(chr==TRUE) {
##         return(gr.chr(gr))
##     } else {
##         return(gr)
##     }
## }


## getPGVold = function(json_file = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/datafiles.json", datadir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/data/",settings = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/settings.json") {
##     pgv = PGVdb$new(datafiles_json_path=json_file, datadir = datadir, settings = settings)
##     return(pgv)
## }

## getPGV = function(json_file = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/datafiles.json", datadir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/data/",settings = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/settings.json") {
##     pgv = PGVdb$new(datafiles_json_path=json_file, datadir = datadir, settings = settings)
##     return(pgv)
## }





## getPGVold = function(json_file = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/datafiles.json", datadir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/data/",publicdir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/public/",settings = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/settings.json") {
##     pgv = PGVdb$new(datafiles_json_path=json_file, datadir = datadir, publicdir = publicdir, settings = settings)
##     return(pgv)
## }




## getPGV = function(json_file = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/datafiles.json", datadir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/data/",publicdir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/public/",settings = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/PGV_NEW/pgv/build/settings.json") {
##     pgv = PGVdb$new(datafiles_json_path=json_file, datadir = datadir, publicdir = publicdir, settings = settings)
##     return(pgv)
## }

## getPGV = function() {

## json_file = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/datafiles.json"
## datadir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/data/"
## publicdir = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/public/"
## settings = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/pgv/build/settings.json" 
## pgv = PGVdb$new(datafiles_json_path=json_file, datadir = datadir, publicdir = publicdir, settings = settings)
##     return(pgv)
## }

    
#old method of pgv database
## #just making the function shorter since I'm always copying and running this
## getPGV = function() {
##     return_PGV_db(datafiles.json = "/gpfs/commons/groups/imielinski_lab/pgv_content/datafiles.json",
##                   data_folder = "/gpfs/commons/groups/imielinski_lab/pgv_content/data/",
##                   PGV_public_dir = "/gpfs/commons/home/cxanthopoulakis/pgv/public/")}

#convert table to ggraph json
## dt2json = function(dt,filename) {
##     gr1 = dt2gr(dt)
##     jab = gG(nodes=gr1)
##     jab = gGnome::refresh(jab)
##     settings = list(y_axis = list(title = "copy number",
##                                   visible = TRUE))
##     node.json = gr2dt(jab$nodes$gr[, "snode.id"])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, y = 1,title=snode.id,type="interval",strand="*")]
##     gg.js = list(intervals = node.json, connections = data.table())
##     gg.js = c(list(settings = settings), gg.js)
##     jsonlite::write_json(gg.js, filename,
##                          pretty=TRUE, auto_unbox=TRUE, digits=4)
## }

idgen = function(pair) {
    paste0(as.data.table(t(as.data.table(unlist(strsplit(pair,split = "-")))))[,c(1:3)],collapse = "-")
}

load_dev = function() {
    ## devtools::load_all("~/git/gGnome")
    ## devtools::load_all("~/git/pgvdb")
    library(skitools); devtools::load_all("~/git/GxG"); library(Flow) ; library(fishHook) ; library(rtracklayer) ; library(pbmcapply)
    setDTthreads(1)
}

## #Clean a model so I'm not saving all of the training data with it
## cleanModel1 = function(cm) {
##     cm$y = c()
##     cm$model = c()
##     cm$residuals = c()
##     cm$fitted.values = c()
##     cm$effects = c()
##     cm$qr$qr = c()
##     cm$linear.predictors = c()
##     cm$weights = c()
##     cm$prior.weights = c()
##     cm$data = c()
##     cm
## }



## #push PGVDb with ordering by the title
## pushDB = function(json_db, backup = TRUE) {
##     json_db = .valid_json_db(json_db)
##                                         # clean up and remove added columns
    
    
##                                         # plots
##     json_db$plots[,plot_id:=NULL]
##                                         # json_db$plots[,file.source:=NULL]
##                                         #
##     pt_ids = unique(json_db$references$patient.id)
##     json_format = lapply(pt_ids, function(x) {
##         description = json_db$descriptions[patient.id == x]$tags
##         plots_Df = json_db$plots[patient.id == x, !"patient.id"]
##         plots_Df = plots_Df %>% mutate(title=factor(title, levels = c(levels = c("Centromeres","ATAC_R2","ATAC_R1","Hits","Coverage","Cov_bw2","Cov_bw","peel","Genome2","Jameson_wg_binsets1","Jameson_wg_binsets2","Binsets_EP_Paper","Concatemers","Genome")))) %>% arrange(desc(title))
## #        plots_Df = plots_Df %>% mutate(title=factor(title, levels = c(levels = c("Genome","peel","Cov_bw","Coverage","Hits","ATAC_R1","ATAC_R2","Centromeres")))) %>% arrange(desc(title))
##         plots_Df$title = as.character(plots_Df$title)

##         return(list(description = description, 
##                     reference = json_db$references[patient.id == x]$reference, 
##                     plots = plots_Df))
##     })
##     names(json_format) = pt_ids
    
##                                         # here we save a backup version with date & time for now
##     file.rename(json_db$datafiles.json, 
##                 paste0(json_db$datafiles.json,
##                        format(Sys.time(), "%Y%m%d_%H%M%S")))
##                                         # re json
    
##     jsonlite::write_json(json_format, json_db$datafiles.json,
##                          pretty=TRUE, auto_unbox=TRUE, digits=4)
## }


## #functions for controlling higlass server :
##                                         #list tilesets
## list_tiles = function() {
##     cmd = ("sh ~/Projects/higlass_serv/UPLOAD_SCRIPTS/list_tilesets.sh")
##     tiles.dt = system(cmd,intern=T)
##     tiles.dt2 = rbindlist(lapply(4:length(tiles.dt),function(x) {
##         as.data.table(t(as.data.table(unlist(strsplit(gsub("]","",tiles.dt[x])," ")))))
##     }))
##     tiles.dt3 = tiles.dt2[,c(4,6,8)]
##     colnames(tiles.dt3) = c("tile","type","uuid")
##     return(tiles.dt3)
## }

##                                         #upload bigwigs
## upload_bigwigs = function(dt,bigwig.col,uuid.col,cores) {
##     dt[,run := paste0("sh ~/Projects/higlass_serv/UPLOAD_SCRIPTS/upload_bigwigs.sh ",dt[[bigwig.col]]," ",dt[[uuid.col]])]
##     upload_higlass_dt = mclapply(1:nrow(dt),function(pair) {
##         cmd=(dt$run[pair])
##         system(cmd,intern=T)
##     },mc.cores=cores)
##     return(upload_higlass_dt)
## }

##                                         #delete tiles
## delete_tiles = function(dt,uuid.col,cores) {
##     dt[,delete := paste0("sh ~/Projects/higlass_serv/UPLOAD_SCRIPTS/delete_tileset.sh ",dt[[uuid.col]])]
##     delete_higlass_dt = mclapply(1:nrow(dt),function(pair) {
##         cmd=(dt$delete[pair])
##         system(cmd,intern=T)
##     },mc.cores=cores)
##     return(delete_higlass_dt)
## }

##                                         #upload chromsizes file
## upload_chromsizes = function(chrom.file) {
##     cmd = paste0("sh ~/Projects/higlass_serv/UPLOAD_SCRIPTS/upload_chromsizes.sh"," ",chrom.file)
##     system(cmd,intern=T)
## }



## #json conversion for peel output with modifying stack gap
## json_gap = function (filename = ".", save = TRUE, verbose = FALSE, annotations = NULL, 
##     nfields = NULL, efields = NULL, stack.gap = 1e+05, include.graph = TRUE, 
##     settings = list(y_axis = list(title = "copy number", visible = TRUE)), 
##     cid.field = NULL, no.y = FALSE,file=NULL,stack.y=1,stack.x=0) {
##     self=file
##     if (length(self) == 0) {
##         warning("This is an empty gWalk so no JSON will be produced.")
##         return(NA)
##     }
##     if (length(self$edges) == 0) {
##         warning("There are no edges in this gWalk so no JSON will be produced.")
##         return(NA)
##     }
##     non.alt.exist = any(self$dt[, sapply(sedge.id, length) == 
##         0] | self$dt[, sapply(sedge.id, anyNA) == TRUE])
##     if (non.alt.exist) {
##         json_ret = refresh(self[self$dt[, sapply(sedge.id, length) >
##                                           0] & self$dt[, sapply(sedge.id, anyNA) == FALSE]])
##         return(json_gap(file= json_ret, filename = filename, 
##                                       save = save, verbose = verbose, annotations = annotations, 
##                                       nfields = nfields, efields = efields, stack.gap = stack.gap, 
##                                       include.graph = include.graph, settings = settings, 
##                                       no.y = no.y,stack.x=stack.x,stack.y=stack.y))
##     }
##     ## if (non.alt.exist) {
##     ##     return(refresh(self[self$dt[, sapply(sedge.id, length) > 
##     ##         0] & self$dt[, sapply(sedge.id, anyNA) == FALSE]])$json(filename = filename, 
##     ##         save = save, verbose = verbose, annotations = annotations, 
##     ##         nfields = nfields, efields = efields, stack.gap = stack.gap, 
##     ##         include.graph = include.graph, settings = settings, 
##     ##         no.y = no.y))
##     ## }
##     if (include.graph) {
##         graph.js = refresh(self$graph)$json(filename = NA, save = FALSE, 
##             verbose = verbose, annotations = annotations, nfields = nfields, 
##             efields = efields, settings = settings, no.y = no.y)
##     }
##     pids = split(self$dt[, .(pid = walk.id, strand = "+", type = ifelse(self$circular, 
##         "cycle", "path"))], 1:self$length)
##     efields = unique(c("type", efields))
##     protected_efields = c("cid", "source", "sink", "title", "weight")
##     rejected_efields = intersect(efields, protected_efields)
##     if (length(rejected_efields) > 0) {
##         warning(sprintf("The following fields were included in efields: \"%s\", but since these are conserved fields in the json walks output then they will be not be included in efields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.", 
##             paste(rejected_efields, collapse = "\" ,\"")))
##         efields = setdiff(efields, rejected_efields)
##     }
##     missing_efields = setdiff(efields, names(self$edges$dt))
##     if (length(missing_efields) > 0) {
##         warning(sprintf("Invalid efields value/s provided: \"%s\". These fields were not found in the gWalk and since will be ignored.", 
##             paste(missing_efields, collapse = "\" ,\"")))
##         efields = intersect(efields, names(self$edges$dt))
##     }
##     sedu = dunlist(self$sedge.id)
##     cids = lapply(unname(split(cbind(data.table(cid = sedu$V1, 
##         source = self$graph$edges[sedu$V1]$left$dt$snode.id, 
##         sink = -self$graph$edges[sedu$V1]$right$dt$snode.id, 
##         title = "", weight = 1), self$graph$edges[sedu$V1]$dt[, 
##         ..efields]), sedu$listid)), function(x) unname(split(x, 
##         1:nrow(x))))
##     snu = dunlist(self$snode.id)
##     snu$ys = gGnome:::draw.paths.y(self$grl,path.stack.x.gap = stack.x ,path.stack.y.gap = stack.y) %>% unlist
##     protected_nfields = c("chromosome", "startPoint", "endPoint", 
##         "y", "type", "strand", "title")
##     rejected_nfields = intersect(nfields, protected_nfields)
##     if (length(rejected_nfields) > 0) {
##         warning(sprintf("The following fields were included in nfields: \"%s\", but since these are conserved fields in the json walks output then they will be not be included in nfields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.", 
##             paste(rejected_nfields, collapse = "\" ,\"")))
##         nfields = setdiff(nfields, rejected_nfields)
##     }
##     missing_nfields = setdiff(nfields, names(self$nodes$dt))
##     if (length(missing_nfields) > 0) {
##         warning(sprintf("Invalid nfields value/s provided: \"%s\". These fields were not found in the gWalk and since will be ignored.", 
##             paste(missing_nfields, collapse = "\" ,\"")))
##         nfields = intersect(nfields, names(self$edges$dt))
##     }
##     iids = lapply(unname(split(cbind(data.table(iid = abs(snu$V1)), 
##         self$graph$nodes[snu$V1]$dt[, .(chromosome = seqnames, 
##             startPoint = start, endPoint = end, y = snu$ys, type = "interval", 
##             strand = ifelse(snu$V1 > 0, "+", "-"), title = abs(snu$V1))], 
##         self$graph$nodes[snu$V1]$dt[, ..nfields]), snu$listid)), 
##         function(x) unname(split(x, 1:nrow(x))))
##     walks.js = lapply(1:length(self), function(x) c(as.list(pids[[x]]), #had to change this
##         list(cids = rbindlist(cids[[x]])), list(iids = rbindlist(iids[[x]]))))
##     if (include.graph) {
##         out = c(graph.js, list(walks = walks.js))
##     }
##     else {
##         out = list(walks = walks.js)
##     }
##     if (save) {
##         if (verbose) {
##             message("Saving JSON to: ", filename)
##         }
##         jsonlite::write_json(out, filename, pretty = TRUE, auto_unbox = TRUE, 
##             digits = 4)
##         return(normalizePath(filename))
##     }
##     else {
##         return(out)
##     }
## }



## #Kat's version in git for json conversion
json_kat = function(filename = '.',
                                      save = TRUE,
                                      verbose = FALSE,
                                      annotations = NULL,
                                      nfields = NULL,
                                      efields = NULL,
                                      stack.gap = 20, 
                                      include.graph = TRUE,
                                      settings = list(y_axis = list(title = "copy number",
                                                                    visible = TRUE)),
                                      cid.field = NULL,
                                      no.y = FALSE,file=NULL) {
    self=file
    if (length(self) == 0) {
                            warning('This is an empty gWalk so no JSON will be produced.')
                            return(NA)
                        }
                        if (length(self$edges) == 0) {
                            warning('There are no edges in this gWalk so no JSON will be produced.')
                            return(NA)
                        }
                        # check if the gWalk includes walks with no ALT edges and if so then remove these and call json again
                        # non.alt.exist = any(self$dt[,sapply(sedge.id, length) == 0])
                        non.alt.exist = any(self$dt[,sapply(sedge.id, length) == 0] | self$dt[,sapply(sedge.id, anyNA) == TRUE])  ## kat: I had to add this bandaid because I was using subsetted gWalk objects, where walks of length one have NAs in place of NULLs for sedge.id... ultimately this should be patched under the hood in the gWalk subset function so that subsetted walk sedge.id output is EXACTLY the same style as original...

                        if (non.alt.exist) {
                          # we call json function again but only including the walks that include at least one ALT edge (anyNA condition added to fix subsetting issue as described in the comments just above)
                          return(refresh(self[self$dt[,sapply(sedge.id, length) > 0] & self$dt[,sapply(sedge.id, anyNA) == FALSE]])$json(filename = filename,
                                                                                         save = save,
                                                                                         verbose = verbose,
                                                                                         annotations = annotations,
                                                                                         nfields = nfields,
                                                                                         efields = efields,
                                                                                         stack.gap = stack.gap,
                                                                                         include.graph = include.graph,
                                                                                         settings = settings,
                                                                                         no.y = no.y))
                        }

                        if (include.graph) {
                            ## build graph level js 
                            graph.js = refresh(self$graph)$json(filename = NA, save = FALSE, verbose = verbose,
                                                 annotations = annotations, nfields = nfields, efields = efields,
                                                 settings = settings, no.y = no.y)
                        }


                        pids = split(self$dt[, .(pid = walk.id,
                                           strand = '+',
                                           type = ifelse(self$circular, 'cycle', 'path'))], 1:self$length)

                        efields = unique(c('type', efields))
                        protected_efields = c('cid', 'source', 'sink', 'title', 'weight')
                        rejected_efields = intersect(efields, protected_efields)
                        if (length(rejected_efields) > 0) {
                            warning(sprintf('The following fields were included in efields: "%s", but since these are conserved fields in the json walks output then they will be not be included in efields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.', paste(rejected_efields, collapse = '" ,"')))
                            efields = setdiff(efields, rejected_efields)
                        }
                        missing_efields = setdiff(efields, names(self$edges$dt))
                        if (length(missing_efields) > 0) {
                            warning(sprintf('Invalid efields value/s provided: "%s". These fields were not found in the gWalk and since will be ignored.', paste(missing_efields, collapse = '" ,"')))
                            efields = intersect(efields, names(self$edges$dt))
                        }

                        sedu = dunlist(self$sedge.id)
                        cids = lapply(unname(split(cbind(data.table(cid = sedu$V1,
                                                       source = self$graph$edges[sedu$V1]$left$dt$snode.id,
                                                       sink = -self$graph$edges[sedu$V1]$right$dt$snode.id, # notice that we need to add negative sign here to meet the gGnome.js expectations
                                                       title = "",
                                                       weight = 1),
                                                   self$graph$edges[sedu$V1]$dt[, ..efields]
                                                 ), sedu$listid)),
                                      function(x) unname(split(x, 1:nrow(x))))

                        snu = dunlist(self$snode.id)
                        snu$ys = gGnome:::draw.paths.y(self$grl, path.stack.y.gap=stack.gap) %>% unlist

                        protected_nfields = c('chromosome', 'startPoint', 'endPoint',
                                              'y', 'type', 'strand', 'title')
                        rejected_nfields = intersect(nfields, protected_nfields)
                        if (length(rejected_nfields) > 0) {
                            warning(sprintf('The following fields were included in nfields: "%s", but since these are conserved fields in the json walks output then they will be not be included in nfields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.', paste(rejected_nfields, collapse = '" ,"')))
                            nfields = setdiff(nfields, rejected_nfields)
                        }
                        missing_nfields = setdiff(nfields, names(self$nodes$dt))
                        if (length(missing_nfields) > 0) {
                            warning(sprintf('Invalid nfields value/s provided: "%s". These fields were not found in the gWalk and since will be ignored.', paste(missing_nfields, collapse = '" ,"')))
                            nfields = intersect(nfields, names(self$edges$dt))
                        }
                        iids = lapply(unname(split(cbind(
                          data.table(iid = abs(snu$V1)),
                          self$graph$nodes[snu$V1]$dt[,
                                                      .(chromosome = seqnames,
                                                        startPoint = start,
                                                        endPoint = end,
                                                        y = snu$ys,
                                                        type = "interval",
                                                        strand = ifelse(snu$V1 > 0, "+", "-"),
                                                        title = abs(snu$V1))],
                          self$graph$nodes[snu$V1]$dt[,..nfields]), snu$listid)),
                          function(x) unname(split(x, 1:nrow(x))))

                        walks.js = lapply(1:length(self), function(x)
                          c(as.list(pids[[x]]),
                            list(cids = rbindlist(cids[[x]])),
                            list(iids = rbindlist(iids[[x]]))))

                        if (include.graph) {
                            out = c(graph.js, list(walks = walks.js))
                        } else {
                            out = list(walks = walks.js)
                        }

                        if (save) {
                          if (verbose) {
                            message("Saving JSON to: ", filename)
                          }
                          jsonlite::write_json(out, filename,
                                               pretty=TRUE, auto_unbox=TRUE, digits=4)
                          return(normalizePath(filename))
                        } else {
                          return(out)
                        }
                      }

## fish_event_covar = function(pairs,hits_gr,event_types,mcores=1,idcap=1,eligible,tiles) {
##     mclapply(setNames(nm = pairs), function(pair) {
##         jab = readRDS(atac.pairs3[pair,complex])
##         hits_sub_gr = hits_gr %Q% (sample == pair)
##         if((length(hits_sub_gr) > 0) & (nrow(jab$meta$events) >0)) {
##             if(length(hits_sub_gr) > 100) {mcores_fish=10} else {mcores_fish=2}
##             events.dt = jab$meta$events[type %in% event_types,]
##             covariates_fish = lapply(unique(events.dt$type), function(svtype) {
##                 grl1 = parse.gr(events.dt[type==svtype,]$footprint,seqlengths = hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")) %>% sortSeqlevels()        
##                 grl1$type = svtype
##                 grl1$grl.ix = NULL
##                 grl1$grl.iix = NULL
##                 if(svtype == "chromoplexy" | svtype == "tic") {
##                     grl1 = (grl1 + 1e4) %>% trim()
##                     }
##                 fish_svtype = fishHook::Cov(data = grl1,name = svtype,pad = 1e4)
##                 return(fish_svtype)
##             })
##             fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=eligible, idcol="sample", mc.cores=mcores_fish,idcap=idcap)
##                                         #    fish_per_pat$score()
##             if(length(covariates_fish)==0) {
##                 return(NULL)
##             }
##             if(length(covariates_fish)>0) {
##                 for (x in 1:length(covariates_fish)) { # have to use a for loop here since lapply keeps it out of memory
##                     fish_per_pat$merge(covariates_fish[[x]])
##                 }
##                 fish_per_pat$score()
##                 return(fish_per_pat)
##             }} else {return(NULL)}
##     },mc.cores=mcores)
## }

## fish_nocovar = function(pairs,hits_gr,event_types,mcores=1,idcap=1,eligible,tiles) {
##     mclapply(setNames(nm = pairs), function(pair) {
##         hits_sub_gr = hits_gr %Q% (sample == pair)
##         if(length(hits_sub_gr) > 100) {mcores_fish=10} else {mcores_fish=2}
##         fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=eligible, idcol="sample", mc.cores=mcores_fish,idcap=idcap)
##         fish_per_pat$score()
##         return(fish_per_pat)
##     },mc.cores=mcores)
## }

#added ability to rerun and rescore only using significant covariates and to start from a saved data table and set threshold for pvalue column
## fish_event_covar = function(pairs,hits_gr,event_types,mcores=1,idcap=1,eligible,tiles,rescore=FALSE) {
##     mclapply(setNames(nm = pairs), function(pair) {
##         jab = readRDS(atac.pairs3[pair,complex])
##         hits_sub_gr = hits_gr %Q% (sample == pair)
##         if((length(hits_sub_gr) > 0) & (nrow(jab$meta$events) >0)) {
##             if(length(hits_sub_gr) > 100) {mcores_fish=10} else {mcores_fish=2}
##             events.dt = jab$meta$events[type %in% event_types,]
##             covariates_fish = lapply(setNames(nm = unique(events.dt$type)), function(svtype) {
##                 grl1 = parse.gr(events.dt[type==svtype,]$footprint,seqlengths = hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")) %>% sortSeqlevels()        
##                 grl1$type = svtype
##                 grl1$grl.ix = NULL
##                 grl1$grl.iix = NULL
##                 if(svtype == "chromoplexy" | svtype == "tic") {
##                     grl1 = (grl1 + 1e4) %>% trim()
##                     }
##                 fish_svtype = fishHook::Cov(data = grl1,name = svtype,pad = 1e4)
##                 return(fish_svtype)
##             })
##             fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=gr.chr(eligible), idcol="sample", mc.cores=mcores_fish,idcap=idcap)
##                                         #    fish_per_pat$score()
##             if(length(covariates_fish)==0) {
##                 return(NULL)
##             }
##             if(length(covariates_fish)>0) {
##                 for (x in 1:length(covariates_fish)) { # have to use a for loop here since lapply keeps it out of memory
##                     fish_per_pat$merge(covariates_fish[[x]])
##                 }
##                 fish_per_pat$score()
##                 if(!rescore) {
##                     return(fish_per_pat)
##                 } else {
##                     summary_model = summary(fish_per_pat$model)
##                     coefficient_summary.dt = as.data.table(summary_model$coefficients,keep.rownames="variable")
##                     names(coefficient_summary.dt) = c("variable","Estimate","Std_Error","t_value","Pr_greater_t")
##                     all_variables = coefficient_summary.dt[!(grepl("Intercept",variable)),]$variable
##                     variables.keep = coefficient_summary.dt[!(grepl("Intercept",variable)) & (Pr_greater_t < 0.05),]$variable
##                     if(length(all_variables) > length(variables.keep)) {
##                         fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=gr.chr(eligible), idcol="sample", mc.cores=mcores_fish,idcap=idcap)
##                         for (x in variables.keep) {
##                             fish_per_pat$merge(covariates_fish[[x]])
##                         }
##                         fish_per_pat$score()
##                     }
##                     return(fish_per_pat)}
##                 }
##             }} else {return(NULL)}
##     },mc.cores=mcores)
## }

fish_event_covar = function(pairs_table,event_types,mcores=1,idcap=1,eligible,tiles,rescore=FALSE,hits_column,pvalue_column,pvalue_threshold,rescore_threshold=0.05) {
    mclapply(setNames(nm = pairs_table$pair), function(pair) {
        jab = readRDS(pairs_table[pair,complex])
        hits_sub_gr = readRDS(pairs_table[pair,get(hits_column)])[get(pvalue_column) <= pvalue_threshold,] %>% GRanges()
        if((length(hits_sub_gr) > 0) & (nrow(jab$meta$events) >0)) {
            if(length(hits_sub_gr) > 100) {mcores_fish=4} else {mcores_fish=1}
            events.dt = jab$meta$events[type %in% event_types,]
            covariates_fish = lapply(setNames(nm = unique(events.dt$type)), function(svtype) {
                grl1 = parse.gr(events.dt[type==svtype,]$footprint,seqlengths = hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")) %>% sortSeqlevels()        
                grl1$type = svtype
                grl1$grl.ix = NULL
                grl1$grl.iix = NULL
                if(svtype == "chromoplexy" | svtype == "tic") {
                    grl1 = (grl1 + 1e4) %>% trim()
                    }
                fish_svtype = fishHook::Cov(data = grl1,name = svtype,pad = 1e4)
                return(fish_svtype)
            })
            fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=gr.chr(eligible), idcol="sample", mc.cores=mcores_fish,idcap=idcap)
                                        #    fish_per_pat$score()
            if(length(covariates_fish)==0) {
                return(NULL)
            }
            if(length(covariates_fish)>0) {
                for (x in 1:length(covariates_fish)) { # have to use a for loop here since lapply keeps it out of memory
                    fish_per_pat$merge(covariates_fish[[x]])
                }
                fish_per_pat$score()
                if(!rescore) {
                    return(fish_per_pat)
                } else {
                    summary_model = summary(fish_per_pat$model)
                    coefficient_summary.dt = as.data.table(summary_model$coefficients,keep.rownames="variable")
                    names(coefficient_summary.dt) = c("variable","Estimate","Std_Error","t_value","Pr_greater_t")
                    all_variables = coefficient_summary.dt[!(grepl("Intercept",variable)),]$variable
                    variables.keep = coefficient_summary.dt[!(grepl("Intercept",variable)) & (Pr_greater_t < rescore_threshold),]$variable
                    if(length(all_variables) > length(variables.keep)) {
                        fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=gr.chr(eligible), idcol="sample", mc.cores=mcores_fish,idcap=idcap)
                        for (x in variables.keep) {
                            fish_per_pat$merge(covariates_fish[[x]])
                        }
                        fish_per_pat$score()
                    }
                    return(fish_per_pat)}
                }
            } else {return(NULL)}
    },mc.cores=mcores)
}

#gethits
fish_hits = function(pairs,fish_covariate_obj,mcores=1) {
    mclapply(setNames(nm=pairs), function(pair) {
        fish_obj = fish_covariate_obj[[pair]]
        if(!is.null(fish_obj)) {
            fish.dt = gr2dt(fish_obj$res)[order(p), ][fdr<0.25, ]
            fish.dt[,sample := pair]            
            fish.dt[,covariates_used := paste0(fish_obj$coefficients$covariate,collapse=",")]
            return(fish.dt)
            } else {return(NULL)}
    },mc.cores=mcores)
}


##get coefficeints from fishhook
get_fish_coef = function(pairs,fish_covariate_obj,mcores=1) {
    coefficients_lst = mclapply(setNames(nm=pairs), function(pair) {
        fish_obj = fish_covariate_obj[[pair]]
        if(!is.null(fish_obj)) {
            coef.dt = fish_obj$coefficients
            coef.dt[, sample := pair]
            return(coef.dt)
        } else { return(NULL)}
    },mc.cores=20)
    coef.dt = rbindlist(coefficients_lst)
    return(coef.dt)
    }

######plot qq plots for the model
fish_qq = function(pairs,fish_covariate_obj,save_file="plot.png",label_count=20,subsample=1e5,max_x = 10,max_y=NULL,mcores=1) {
    fish_res.lst = mclapply(setNames(nm=pairs), function(pair) {
        fish_obj = fish_covariate_obj[[pair]]
        if(!is.null(fish_obj)) {
            dt1 = as.data.table(fish_obj$res)
            return(dt1)
            } else {return(NULL)}
    },mc.cores=mcores)
    fish_res.dt = rbindlist(fish_res.lst,fill=TRUE)
    fish_res.dt = fish_res.dt[order(p),]
    fish_res.dt[1:label_count, label_qq := nearest.gene]
    ppng(qq_pval(fish_res.dt$p,subsample=subsample,max.x=max_x,max.y = max_y,label=fish_res.dt$label_qq,repel=TRUE),file=save_file)
}

##match hits up with whether they are in specific complex events
match_hits_events = function(pairs,fish_hits_gr,mcores=1) {
    fish_hits_events_lst = mclapply(setNames(nm=pairs), function(pair) {
        jab = readRDS(atac.pairs3[pair,complex])
        fish_hits_sub_gr = fish_hits_gr %Q% (sample == pair)
        if(nrow(jab$meta$events) >0) {
            events.dt = jab$meta$events[type %in% event_types,]
            if(length(fish_hits_sub_gr)==0) {
                return(NULL)
            } else {
                if(nrow(events.dt)==0) {
                    return(as.data.table(fish_hits_sub_gr))
                } else {
                    events_gr = lapply(unique(events.dt$type), function(svtype) {
                        grl1 = parse.gr(events.dt[type==svtype,]$footprint,seqlengths = hg_seqlengths("~/DB/UCSC/hg38.chrom.sizes")) %>% sortSeqlevels()        
                        grl1$type = svtype
                        grl1$grl.ix = NULL
                        grl1$grl.iix = NULL
                        if(svtype == "chromoplexy" | svtype == "tic") {
                            grl1 = (grl1 + 1e4) %>% trim()
                        }
                        return(grl1)
                    }) %>% do.call(c,.)
                    fish_hits_sub_dt = as.data.table(gr.val(fish_hits_sub_gr, events_gr, "type"))
                    fish_hits_sub_dt[type == "", type := NA]
                    fish_hits_sub_dt[covariates_used == "", covariates_used := NA]
                    return(fish_hits_sub_dt)
                }}} else {
                      if(length(fish_hits_sub_gr)==0) {
                          return(NULL)
                      } else {return(fish_hits_sub_gr)}
                      }
    },mc.cores=mcores)
    return(rbindlist(fish_hits_events_lst,fill=TRUE))
}

#####plot a qq_plot for each patient
qq_per_pair = function(pairs,fish_covariate_obj,folder,suffix="_qq.png",mcores=1) {
    mclapply(setNames(nm=pairs), function(pair) {
        fish_obj = fish_covariate_obj[[pair]]
        if(!is.null(fish_obj)) {
            ppng(fish_obj$qqp(plotly=FALSE),file=paste0(folder,pair,suffix))
        }
    },mc.cores=mcores)
}

fish_nocovar = function(pairs_table,mcores=1,idcap=1,eligible,tiles,hits_column,pvalue_column,pvalue_threshold) {
    mclapply(setNames(nm = pairs_table$pair), function(pair) {
        hits_sub_gr = readRDS(pairs_table[pair,get(hits_column)])[get(pvalue_column) <= pvalue_threshold,] %>% GRanges()
        if(length(hits_sub_gr) > 100) {mcores_fish=10} else {mcores_fish=2}
        fish_per_pat = Fish(hypotheses=tiles, events=hits_sub_gr, eligible=eligible, idcol="sample", mc.cores=mcores_fish,idcap=idcap)
        fish_per_pat$score()
        return(fish_per_pat)
    },mc.cores=mcores)
}

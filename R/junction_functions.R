
#create function that swaps peaks or swaps junctions for permutting peaks (normalized intensity y)
alift_permute_peaks2 = function(random_pairs,bin_size,counts_file_col,rand_counts_file_col,junctions_gr,last_bin = -5e7,pval_column,padj_column,FDR_thresh = 0.25,swap = "peaks",mcores=1) {
    mclapply(setNames(nm = random_pairs$pair),function(pair) {
        dt.sub = readRDS(random_pairs[pair,get(counts_file_col)])
        dt.sub[get(padj_column) <= FDR_thresh,type_peak := "Significant Peak"]
        dt.sub[get(padj_column) > FDR_thresh,type_peak := "Not Significant Peak"]
        #dt.sub = counts.dt[pair,]
        rand.pair = random_pairs[pair,random_pair]
                                        #junctions from pair
        junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
        counts.gr = dt2gr(dt.sub)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr,
                                          window = 5e7)
        alift.pks.dt = as.data.table(alift.pks.gr)
        alift.pks.dt[,permute := "not_permute"]
                                        #junctions from random pair within project
        ## junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
                                        #swapping peaks not junction
        if(swap == "peaks") {
            dt.sub2 = readRDS(random_pairs[random_pair==rand.pair,get(rand_counts_file_col)])
            counts.gr2 = dt2gr(dt.sub2)
            alift.pks.gr = gUtils::anchorlift(counts.gr2,
                                              junctions.sub.gr,
                                              window = -last_bin)
        }
        if(swap == "junctions") {
            junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
            alift.pks.gr = gUtils::anchorlift(counts.gr,
                                              junctions.sub.gr2,
                                              window = -last_bin)
        }
        alift.pks.dt2 = as.data.table(alift.pks.gr)
        alift.pks.dt2[,permute := "permute"]
        alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
        alift.pks.dt[,uniq := "not_unique"]
        if(length(names(alift.pks.dt)) > 4) {
            setDT(alift.pks.dt)
            alift.pks.dt[, dis := ((start + end)/2)]
            alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
            alift.pks.dt = alift.pks.dt[,.(sample,type,permute,uniq,get(pval_column),get(padj_column),bin,type_peak)]
            names(alift.pks.dt) = c("sample","type","permute","uniq",pval_column,padj_column,"bin","type_peak")
            return(alift.pks.dt)
        } else {return(NULL)}
    }, mc.cores=mcores)
}

#create function that swaps peaks or swaps junctions for permutting lambda-does not require input of padj_column which is required to get peak enrichment (defining peaks requires the adjusted p value where as calculating lambda does not)
alift_permute_lambda2 = function(random_pairs,bin_size,counts_file_col,rand_counts_file_col,junctions_gr,last_bin = -5e7,pval_column,swap = "peaks",mcores=1) {
    mclapply(setNames(nm = random_pairs$pair),function(pair) {
        dt.sub = readRDS(random_pairs[pair,get(counts_file_col)])
        rand.pair = random_pairs[pair,random_pair]
                                        #junctions from pair
        junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
        counts.gr = dt2gr(dt.sub)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr,
                                          window = 5e7)
        alift.pks.dt = as.data.table(alift.pks.gr)
        alift.pks.dt[,permute := "not_permute"]
                                        #junctions from random pair within project
        ## junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
                                        #counts from random pair within project
        if(swap == "peaks") {
            dt.sub2 = readRDS(random_pairs[random_pair==rand.pair,get(rand_counts_file_col)])
            counts.gr2 = dt2gr(dt.sub2)
            alift.pks.gr = gUtils::anchorlift(counts.gr2,
                                              junctions.sub.gr,
                                              window = -last_bin)
        }
        if(swap == "junctions") {
            junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
            alift.pks.gr = gUtils::anchorlift(counts.gr,
                                              junctions.sub.gr2,
                                              window = -last_bin)
        }
        alift.pks.dt2 = as.data.table(alift.pks.gr)
        alift.pks.dt2[,permute := "permute"]
        alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
        alift.pks.dt[,uniq := "not_unique"]
        if(length(names(alift.pks.dt)) > 4) {
            setDT(alift.pks.dt)
            alift.pks.dt[, dis := ((start + end)/2)]
            alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
            alift.pks.dt = alift.pks.dt[,.(sample,permute,uniq,get(pval_column),bin)]
            names(alift.pks.dt) = c("sample","permute","uniq",pval_column,"bin")
            return(alift.pks.dt)
        } else {return(NULL)}
    }, mc.cores=mcores)
}

#get junctions granges and add data from the dt such as SPAN and class
get_breakends = function(complex, return.type = "data.table", verbose = FALSE) {
                                        #gg = gG(jabba = jabba_rds)
    gg = readRDS(complex)
    loose = gg$loose %Q% (!terminal)
    if (length(loose)) {
        loose = loose[, c("node.id", "node.cn", "cn")]
    }
    if (length(gg$junctions) > 0 && length(gg$junctions[cn > 0 & class != "REF"]) > 0) {
        junctions = stack(gg$junctions[cn > 0 & class != "REF"]$grl)
        cols.keep = c('node.id', 'edge.id', 'linkedsv','bfb','chromoplexy','chromothripsis','del','dm','dup','pyrgo','rigma','simple','tic','tyfonas','coord','mcoord','ecluster','cpxdm','svtype','sedge.id','cluster')
        cols = intersect(cols.keep, names(values(junctions)))
#        cols = intersect(c("node.id", "edge.id", "linkedsv"), names(values(junctions)))
        junctions = junctions[, cols]
        values(junctions)[, "node.cn"] = gg$nodes$dt[values(junctions)[, "node.id"], cn]
        values(junctions)[, "cn"] = gg$edges$dt[values(junctions)[, "edge.id"], cn]
        values(junctions)[, "jstring"] = rep(grl.string(gg$junctions[cn > 0 & class != "REF"]$grl),
                                             each = 2)
                                        #add to breakend code to get class and span so I can filter on that
        bp.dt = gg$junctions[type=="ALT" & cn > 0]$dt
        values(junctions)[,"SPAN"] = rep(bp.dt[rep(from %in% values(junctions)[, "node.id"] | to  %in% values(junctions)[, "node.id"],length=.N),SPAN],each=2) # adds the span for each junction to both of the +/- node edges involved
        values(junctions)[, "class"] = gg$edges$dt[values(junctions)[, "edge.id"], class]
                                        #add the distance to the nearest junction that's not the same edge id
        junctions = lapply(unique(junctions$edge.id), function(edge_id) {
            junctions_sub = junctions %Q% (edge.id %in% edge_id)
            junctions_not_same_edge = junctions %Q% (edge.id != edge_id)
            distance1 = distanceToNearest(junctions_sub,junctions_not_same_edge,ignore.strand=TRUE)@elementMetadata$distance
            if(length(distance1)==2) {
                junctions_sub$dist_nearest_j = distance1
            } else { junctions_sub$dist_nearest_j = NA}
            return(junctions_sub)
        }) %>% do.call(c,.)
                                        #add loose ends and distance to nearest junction
        merge.dt = data.table(value =1:length(junctions))
        distance.dt = distanceToNearest(junctions,loose) %>% as.data.table()
        distance.dt = merge.data.table(merge.dt,distance.dt,all=TRUE,by.x = "value", by.y = "queryHits")
        junctions$dist_to_loose = distance.dt$distance
        out = grbind(loose, junctions)
        values(out)[, "loose.or.junction"] = c(rep("loose", times = length(loose)),
                                               rep("junction", times = length(junctions)))
    }     else {
              out = loose
              values(out)[, "loose.or.junction"] = rep("loose", times = length(loose))
          }
if (return.type == "data.table") {
    return(as.data.table(out))
}
return(out)
}


## #need to modify such that it's not swapping junctions but instead swapping the counts
## alift_permute_lambda2 = function(random_pairs,bin_size,counts_file_col,rand_counts_file_col,junctions_gr,last_bin = -5e7,pval_column,mcores=1) {
##     mclapply(setNames(nm = random_pairs$pair),function(pair) {
##         dt.sub = readRDS(random_pairs[pair,get(counts_file_col)])
##         rand.pair = random_pairs[pair,random_pair]
##                                         #junctions from pair
##         junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
##         counts.gr = dt2gr(dt.sub)
##         alift.pks.gr = gUtils::anchorlift(counts.gr,
##                                           junctions.sub.gr,
##                                           window = 5e7)
##         alift.pks.dt = as.data.table(alift.pks.gr)
##         alift.pks.dt[,permute := "not_permute"]
##                                         #junctions from random pair within project
##         ## junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
##                                         #counts from random pair within project
##         dt.sub2 = readRDS(random_pairs[random_pair==rand.pair,get(rand_counts_file_col)])
##         counts.gr2 = dt2gr(dt.sub2)
##         alift.pks.gr = gUtils::anchorlift(counts.gr2,
##                                           junctions.sub.gr,
##                                           window = 5e7)
##         alift.pks.dt2 = as.data.table(alift.pks.gr)
##         alift.pks.dt2[,permute := "permute"]
##         alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
##         alift.pks.dt[,uniq := "not_unique"]
##         if(length(names(alift.pks.dt)) > 4) {
##             setDT(alift.pks.dt)
##             alift.pks.dt[, dis := ((start + end)/2)]
##             alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
##             alift.pks.dt = alift.pks.dt[,.(sample,permute,uniq,get(pval_column),bin)]
##             names(alift.pks.dt) = c("sample","permute","uniq",pval_column,"bin")
##             return(alift.pks.dt)
##         } else {return(NULL)}
##     }, mc.cores=mcores)
## }


#I don't have to load each patient in at once but I need to rewrite the alift_permute_peaks and alift_permute_lambda functions
alift_permute_lambda = function(random_pairs,bin_size,counts_file_col,junctions_gr,last_bin = -5e7,pval_column,mcores=1) {
    mclapply(setNames(nm = random_pairs$pair),function(pair) {
        dt.sub = readRDS(random_pairs[pair,get(counts_file_col)])
        rand.pair = random_pairs[pair,random_pair]
                                        #junctions from pair
        junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
        counts.gr = dt2gr(dt.sub)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr,
                                          window = 5e7)
        alift.pks.dt = as.data.table(alift.pks.gr)
        alift.pks.dt[,permute := "not_permute"]
                                        #junctions from random pair within project
        junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr2,
                                          window = 5e7)
        alift.pks.dt2 = as.data.table(alift.pks.gr)
        alift.pks.dt2[,permute := "permute"]
        alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
        alift.pks.dt[,uniq := "not_unique"]
        if(length(names(alift.pks.dt)) > 4) {
            setDT(alift.pks.dt)
            alift.pks.dt[, dis := ((start + end)/2)]
            alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
            alift.pks.dt = alift.pks.dt[,.(sample,type,permute,uniq,get(pval_column),bin)]
            names(alift.pks.dt) = c("sample","type","permute","uniq",pval_column,"bin")
            return(alift.pks.dt)
        } else {return(NULL)}
    }, mc.cores=mcores)
}

#modify the alift_permute lambda function to work with significant and not significant peaks
alift_permute_peaks = function(random_pairs,bin_size,counts_file_col,junctions_gr,last_bin = -5e7,pval_column,mcores=1) {
    mclapply(setNames(nm = random_pairs$pair),function(pair) {
        dt.sub = readRDS(random_pairs[pair,get(counts_file_col)])
        #dt.sub = counts.dt[pair,]
        rand.pair = random_pairs[pair,random_pair]
                                        #junctions from pair
        junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
        counts.gr = dt2gr(dt.sub)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr,
                                          window = 5e7)
        alift.pks.dt = as.data.table(alift.pks.gr)
        alift.pks.dt[,permute := "not_permute"]
                                        #junctions from random pair within project
        junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
        alift.pks.gr = gUtils::anchorlift(counts.gr,
                                          junctions.sub.gr2,
                                          window = 5e7)
        alift.pks.dt2 = as.data.table(alift.pks.gr)
        alift.pks.dt2[,permute := "permute"]
        alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
        alift.pks.dt[,uniq := "not_unique"]
        if(length(names(alift.pks.dt)) > 4) {
            setDT(alift.pks.dt)
            alift.pks.dt[, dis := ((start + end)/2)]
            alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
            alift.pks.dt = alift.pks.dt[,.(sample,type,permute,uniq,get(pval_column),bin,type_peak)]
            names(alift.pks.dt) = c("sample","type","permute","uniq",pval_column,"bin","type_peak")
            return(alift.pks.dt)
        } else {return(NULL)}
    }, mc.cores=mcores)
}
## alift_permute_lambda = function(random_pairs,bin_size,counts_dt,junctions_gr,last_bin = -5e7,pval_column,mcores=1) {
##     mclapply(setNames(nm = random_pairs$pair),function(pair) {
##         dt.sub = counts.dt[pair,]
##         rand.pair = random_pairs[pair,random_pair]
##                                         #junctions from pair
##         junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
##         counts.gr = dt2gr(dt.sub)
##         alift.pks.gr = gUtils::anchorlift(counts.gr,
##                                           junctions.sub.gr,
##                                           window = 5e7)
##         alift.pks.dt = as.data.table(alift.pks.gr)
##         alift.pks.dt[,permute := "not_permute"]
##                                         #junctions from random pair within project
##         junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
##         alift.pks.gr = gUtils::anchorlift(counts.gr,
##                                           junctions.sub.gr2,
##                                           window = 5e7)
##         alift.pks.dt2 = as.data.table(alift.pks.gr)
##         alift.pks.dt2[,permute := "permute"]
##         alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
##         alift.pks.dt[,uniq := "not_unique"]
##         if(length(names(alift.pks.dt)) > 4) {
##             setDT(alift.pks.dt)
##             alift.pks.dt[, dis := ((start + end)/2)]
##             alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
##             alift.pks.dt = alift.pks.dt[,.(sample,type,permute,uniq,get(pval_column),bin)]
##             names(alift.pks.dt) = c("sample","type","permute","uniq",pval_column,"bin")
##             return(alift.pks.dt)
##         } else {return(NULL)}
##     }, mc.cores=mcores)
##     }

#calculate dispersion function
get_dispersion = function(bins_list, anchor_lift_dt,pval_column, mcores=1) {
    mclapply(bins_list, function(bin1) {
        lambda_not_permute1 = lambda_pval(with(anchor_lift_dt[bin == bin1 & permute == "not_permute",],get(pval_column)),plotly=TRUE)
        lambda_permute1 = lambda_pval(with(anchor_lift_dt[bin == bin1 & permute == "permute",],get(pval_column)),plotly=TRUE)
        dt.long = data.table(lambda = c(lambda_not_permute1,lambda_permute1),bin = bin1, permute = c("not_permuted","permuted"))
        return(dt.long)
    },mc.cores=mcores)
    }


get_aggregate_permute = function(bins_list,anchor_lift_dt,normalized_distance=2e7) {
                                        #anchor_lift_dt[, type_peak_junction := paste(type,type.peak)]
    pt.pks.dt = anchor_lift_dt[, .(N = .N), by = .(bin, type_peak,permute)]
    pt.pks.dt[,type.peak.permute := (paste0(type_peak,"_",permute))]
    pt.pks.dt$type.peak.permute = gsub(" ","_",pt.pks.dt$type.peak.permute)
#    pt.pks.dt[, norm.factor := mean(.SD[bin > 2.5e7 | bin < -2.5e7,N]), by = .(type.peak.permute)]
    pt.pks.dt[, norm.factor := mean(.SD[bin > normalized_distance | bin < -normalized_distance,N]), by = .(type.peak.permute)]
    pt.pks.dt[, y_value := N/norm.factor,by=.(type.peak.permute)]
    return(pt.pks.dt)
    }

## get_aggregate_permute = function(bins_list,anchor_lift_dt) {
##                                         #anchor_lift_dt[, type_peak_junction := paste(type,type.peak)]
##     pt.pks.dt = anchor_lift_dt[, .(N = .N), by = .(bin, type_peak,permute)]
##     pt.pks.dt[,type.peak.permute := (paste0(type_peak,"_",permute))]
##     pt.pks.dt$type.peak.permute = gsub(" ","_",pt.pks.dt$type.peak.permute)
##     pt.pks.dt[, norm.factor := mean(.SD[bin > 2.5e7 | bin < -2.5e7,N]), by = .(type.peak.permute)]
##     pt.pks.dt[, y_value := N/norm.factor,by=.(type.peak.permute)]
##     return(pt.pks.dt)
## }


                                        #add better ability to control y axis limits
plot_junctions = function(plot_dt,x = "bin",y, color_by,save_file="plot.png",xlab="distance from junction",ylab = "normalized intensity",plot_title=NULL,y0 = NULL, y1 = NULL,res1=300,height1=3000,width1=4000) {
    if(is.null(y1)) {
        y1 = as.integer(max(with(plot_dt,get(y)))+1)
        if(y1 < 2) {
            y1 =2
        }
    }
    if(is.null(y0)) {
        y0 = 0
    }
    pt = ggplot(plot_dt, aes(x = get(x), y = get(y), color = get(color_by))) +
        geom_line() +
        xlim(min(plot_dt$bin), max(plot_dt$bin)) +
        ylim(y0, y1) +
        ggpubr::theme_pubr() +
        labs(x = xlab, y = ylab) +
        theme(legend.position = "right") +
        geom_vline(xintercept=0, linetype="dotted") +
        scale_color_discrete(name = color_by) +
        ggtitle(plot_title)
    ppng(print(pt),file=save_file,res=res1,height=height1,width=width1)
}


## plot_junctions = function(plot_dt,x = "bin",y, color_by,save_file="plot.png",xlab="distance from junction",ylab = "normalized intensity",plot_title=NULL,res1=300,height1=3000,width1=4000) {
##     max.value = as.integer(max(with(plot_dt,get(y)))+1)
##     if(max.value < 2) {max.value =2}
##     pt = ggplot(plot_dt, aes(x = get(x), y = get(y), color = get(color_by))) +
##         geom_line() +
##         xlim(min(plot_dt$bin), max(plot_dt$bin)) +
##         ylim(0, max.value) +
##         ggpubr::theme_pubr() +
##         labs(x = xlab, y = ylab) +
##         theme(legend.position = "right") +
##         geom_vline(xintercept=0, linetype="dotted") +
##         scale_color_discrete(name = color_by) +
##         ggtitle(plot_title)
##     ppng(print(pt),file=save_file,res=res1,height=height1,width=width1)
## }

## plot_junctions = function(plot_dt,x = "bin",y, color_by,save_file="plot.png",xlab="distance from junction",ylab = "normalized intensity",plot_title=NULL,res1=300,height1=3000,width1=4000) {
##     max.value = as.integer(max(with(plot_dt,get(y)))+1)
##     if(max.value < 2) {max.value =2}
##     pt = ggplot(plot_dt, aes(x = get(x), y = get(y), color = get(color_by))) +
##         geom_line() +
##         xlim(-5e7, 5e7) +
##         ylim(0, max.value) +
##         ggpubr::theme_pubr() +
##         labs(x = xlab, y = ylab) +
##         theme(legend.position = "right") +
##         geom_vline(xintercept=0, linetype="dotted") +
##         scale_color_discrete(name = color_by) +
##         ggtitle(plot_title)
##     ppng(print(pt),file=save_file,res=res1,height=height1,width=width1)
## }

## #modify the alift_permute lambda function to work with significant and not significant peaks
## alift_permute_peaks = function(random_pairs,bin_size,counts_dt,junctions_gr,last_bin = -5e7,pval_column,mcores=1) {
##     mclapply(setNames(nm = random_pairs$pair),function(pair) {
##         dt.sub = counts.dt[pair,]
##         rand.pair = random_pairs[pair,random_pair]
##                                         #junctions from pair
##         junctions.sub.gr = junctions_gr %Q% (sample %in% pair)
##         counts.gr = dt2gr(dt.sub)
##         alift.pks.gr = gUtils::anchorlift(counts.gr,
##                                           junctions.sub.gr,
##                                           window = 5e7)
##         alift.pks.dt = as.data.table(alift.pks.gr)
##         alift.pks.dt[,permute := "not_permute"]
##                                         #junctions from random pair within project
##         junctions.sub.gr2 = junctions_gr %Q% (sample==rand.pair)
##         alift.pks.gr = gUtils::anchorlift(counts.gr,
##                                           junctions.sub.gr2,
##                                           window = 5e7)
##         alift.pks.dt2 = as.data.table(alift.pks.gr)
##         alift.pks.dt2[,permute := "permute"]
##         alift.pks.dt = rbind(alift.pks.dt,alift.pks.dt2,fill=TRUE)
##         alift.pks.dt[,uniq := "not_unique"]
##         if(length(names(alift.pks.dt)) > 4) {
##             setDT(alift.pks.dt)
##             alift.pks.dt[, dis := ((start + end)/2)]
##             alift.pks.dt[, bin := zitools:::easy.cut((start + end)/2, last_bin, step = bin_size)]
##             alift.pks.dt = alift.pks.dt[,.(sample,type,permute,uniq,get(pval_column),bin,type_peak)]
##             names(alift.pks.dt) = c("sample","type","permute","uniq",pval_column,"bin","type_peak")
##             return(alift.pks.dt)
##         } else {return(NULL)}
##     }, mc.cores=mcores)
## }




lambda_pval = function(obs, highlight = c(), exp = NULL, lwd = 1, bestfit=T, col = NULL, col.bg='black', pch=18, cex=1, conf.lines=FALSE, max=NULL, max.x = NULL, max.y = NULL, qvalues=NULL, label = NULL, repel = FALSE, plotly = FALSE, annotations = list(), gradient = list(), titleText = "", subsample = NA, ...)
{
  if(!(plotly))
  {
    is.exp.null = is.null(exp)

    if (is.null(col))
      col = rep('black', length(obs))

    ix1 = !is.na(obs)
    if (!is.null(exp))
      if (length(exp) != length(obs))
        stop('length of exp must be = length(obs)')
      else
        ix1 = ix1 & !is.na(exp)

    if (is.null(highlight))
      highlight = rep(FALSE, length(obs))

    if (is.null(label))
      label = rep('', length(label))

    else if (is.logical(highlight))
    {
      if (length(highlight) != length(obs))
        stop('highlight must be either logical vector of same length as obs or a vector of indices')
    }
    else
      highlight = 1:length(obs) %in% highlight

    obs = -log10(obs[ix1])
    col = col[ix1]
    highlight = highlight[ix1]
    label = label[ix1]
    
    if (!is.null(exp))
      exp = -log10(exp[ix1])

    ix2 = !is.infinite(obs)
    if (!is.null(exp))
      ix2 = ix2 &  !is.infinite(exp)

    obs = obs[ix2]
    col = col[ix2]
    highlight = highlight[ix2]
    label = label[ix2]

    if (!is.null(exp))
      exp = exp[ix2]

    N <- length(obs)
    ## create the null distribution
    ## (-log10 of the uniform)

    if (is.null(exp))
      exp <- -log(1:N/N,10)
    else
      exp = sort(exp)

    if (is.null(max))
      max = max(obs,exp) + 0.5

    if (!is.null(max) & is.null(max.x))
      max.x = max

    if (!is.null(max) & is.null(max.y))
      max.y  = max

    if (is.null(max.x))
      max.x <- max(obs,exp) + 0.5

    if (is.null(max.y))
      max.y <- max(obs,exp) + 0.5

    if (is.exp.null)
    {
      tmp.exp = rev(seq(0, 7, 0.01))
      ix = 10^(-tmp.exp)*N
      c95 <-  qbeta(0.975,ix,N-ix+1)
      c05 <-  qbeta(0.025,ix,N-ix+1)

      if (conf.lines){
        ## plot the two confidence lines
        plot(tmp.exp, -log(c95,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)
        plot(tmp.exp, -log(c05,10), ylim=c(0,max.y), xlim=c(0,max.x), type="l", axes=FALSE, xlab="", ylab="")
        par(new=T)

        p1 <- rep(tmp.exp[1], 2)
        p2 <- c(-log(c95,10)[1], -log(c05,10)[1])

        lines(x=p1, y=p2)
        x.coords <- c(tmp.exp,rev(tmp.exp))
        y.coords <- c(-log(c95,10),rev(-log(c05,10)))
        polygon(x.coords, y.coords, col='light gray', border=NA)
        par(new=T)
      }
    }

    ord = order(obs)

                                        #colors = vector(mode = "character", length = length(obs)); colors[] = "black";

    colors = col
    colors[highlight] = "red";

    dat = data.table(x = sort(exp), y = obs[ord], colors = colors[ord], label = label[ord], pch = pch, cex = cex)

    if (!is.null(names(obs)))
    {
      dat$names = names(obs[ord])
      setkey(dat, names)
    }

    if (nrow(dat)>1e5) ## rough guide to subsampling the lower p value part of the plot
      subsample = 5e4/nrow(dat)

    lambda = lm(y ~ x-1, dat)$coefficients;

    if (is.na(subsample[1]))
      dat[, plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
    else
    {
      subsample = pmin(pmax(0, subsample[1]), 1)
      dat[ifelse(x<=2, ifelse(runif(length(x))<subsample, TRUE, FALSE), TRUE), plot(x, y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, max.x), col = colors, ylim = c(0, max.y), pch=pch, cex=cex, bg=col.bg, ...)]
    }
    
    if (!is.null(dat$label) && any(nchar(dat$label)>0, na.rm = TRUE))
    {
      dat[nchar(label)>0, text(x, y, labels=label, pos=3)];
    }
    
    lines(x=c(0, max(max.y, max.x)), y = c(0, max(max.x, max.y)), col = "black", lwd = lwd)
    
    if (!is.na(subsample))
      dat = dat[sample(nrow(dat), subsample*nrow(dat)), ]


    lines(x=c(0, max.x), y = c(0, lambda*max.x), col = "red", lty = 2, lwd = lwd);
    legend('bottomright',sprintf('lambda=\n %.2f', lambda), text.col='red', bty='n')
  }
  else{

                                        #browser()

        if(length(annotations) < 1){
            hover <- do.call(cbind.data.frame, list(p = obs))
        }
        else{
            hover <- do.call(cbind.data.frame, list(annotations, p = obs))
        }
        hover <- as.data.table(hover)


        is.exp.null = is.null(exp)
        if (is.null(col))
            col = "black"
        ix1 = !is.na(hover$p)
        if (!is.null(exp))
            if (length(exp) != length(hover$p))
                stop("length of exp must be = length(hover$obs)")
            else ix1 = ix1 & !is.na(exp)
        if (is.null(highlight))
            highlight = rep(FALSE, length(hover$p))
        else if (is.logical(highlight)) {
            if (length(highlight) != length(hover$p))
                stop("highlight must be either logical vector of same length as obs or a vector of indices")
        }
        else highlight = 1:length(hover$p) %in% highlight
        hover$obs = -log10(hover$p[ix1])
        hover = hover[ix1]
        highlight = highlight[ix1]
        if (!is.null(exp))
            exp = -log10(exp[ix1])
        ix2 = !is.infinite(hover$obs)
        if (!is.null(exp))
            ix2 = ix2 & !is.infinite(exp)
        hover = hover[ix2]
        highlight = highlight[ix2]
        if (!is.null(exp))
            exp = exp[ix2]
        N <- length(hover$obs)
        if (is.null(exp))
            exp <- -log(1:N/N, 10)
        else exp = sort(exp)
        if (is.null(max))
            max <- max(hover$obs, exp) + 0.5
        else max <- max
        if (is.exp.null) {
            tmp.exp = rev(seq(0, 7, 0.01))
            ix = 10^(-tmp.exp) * N
            c95 <- qbeta(0.975, ix, N - ix + 1)
            c05 <- qbeta(0.025, ix, N - ix + 1)
            if (FALSE) {   ##Don't need if not using conf.line (might put this in the future)
                plot(tmp.exp, -log(c95, 10), ylim = c(0, max), xlim = c(0, max),
                     type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                plot(tmp.exp, -log(c05, 10), ylim = c(0, max), xlim = c(0, max),
                     type = "l", axes = FALSE, xlab = "", ylab = "")

                par(new = T)
                p1 <- rep(tmp.exp[1], 2)
                p2 <- c(-log(c95, 10)[1], -log(c05, 10)[1])
                lines(x = p1, y = p2)
                x.coords <- c(tmp.exp, rev(tmp.exp))
                y.coords <- c(-log(c95, 10), rev(-log(c05, 10)))
                polygon(x.coords, y.coords, col = "light gray", border = NA)
                par(new = T)
            }
        }

                                        #creating the ploting data.table (dat) and organizing the annotations to create hover text
        ord = order(hover$obs)
        hover = hover[ord]
        dat = hover
        hover$obs = NULL

                                        #Creating the hover text
        if(length(colnames(hover)) > 1){
            annotation_names  = sapply(colnames(hover), paste0, " : ")
            annotation_names_wLineBreak  = paste("<br>", annotation_names[2:length(annotation_names)],
                                                 sep = "")
            annotation_names = c(annotation_names[1], annotation_names_wLineBreak)
        }
        else{
            annotation_names  = sapply(colnames(hover), paste0, " : ")
        }

                                        #Checking if there is a gradient and if so adding it to the plotting data.table (dat)
        gradient_control = FALSE
        if(length(gradient )!= 0){
            dat$grad = gradient[[1]][ord]
            gradient_control = TRUE
        }
        else {
            dat$grad = c()
        }


        dat$x = sort(exp)
        dat$y = dat$obs

                                        #declare so we can use in If statement
        p <- NULL

                                        #hacky subsampling but works really well, just maxing out the number of points at 8k
                                        #and removing the extra from the non-sig
                                        #(looks to be -logp of 2.6 here can make this more dynamic later )

        if (nrow(dat) <=  8000){

            dat4 = dat
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }
                                        #            browser()
            if(gradient_control){
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text, color = grad,
                                   colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                   mode = "markers",type = 'scatter')
                         %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                    yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat[, plot_ly(data = dat, x=x, y=y, hoverinfo = "text",text = hover_text,
                                   mode = "markers",type = 'scatter')
                         %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                    yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
        }


        else {

            dat$ID = c(1:nrow(dat))
            dat2 = dat[ y < 2.6,]
            dat3 = as.data.frame(dat2)
            dat3 = as.data.table(dat3[ sample(nrow(dat3), min(4000,nrow(dat3))), ])
            dat2 = rbind(dat3,dat[!(ID%in%dat2$ID),])
            dat2$ID = NULL

            dat4 = dat2
            dat4$obs = NULL
            dat4$x = NULL
            dat4$y = NULL
            dat4$grad = NULL

            trans = t(dat4)
            hover_text = c()
            for (i in 1:dim(trans)[2]){
                outstr = paste(c(rbind(annotation_names, trans[,i])), sep = "", collapse = "")
                hover_text = c(hover_text,outstr)
            }

            if(gradient_control){
                p <- dat2[, plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text, color = grad,
                                    colors = c("blue2","gold"),marker = list(colorbar = list(title = names(gradient[1]))),
                                    mode = "markers",type = 'scatter')
                          %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                     yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }
            else{
                p <- dat2[,  plot_ly(data = dat2, x=x, y=y,hoverinfo = "text", text = hover_text,
                                     mode = "markers",type = 'scatter')
                          %>% layout(xaxis = list(title = "<i>Expected -log<sub>10</sub>(P)</i>"),
                                     yaxis = list(title = "<i>Observed -log<sub>10</sub>(P)</i>")) ]
            }

        }

                                        #       browser()

                                        #Calculating lambda, Note that this is using the whole data set not the subsampled one
        lambda = lm(y ~ x - 1, dat)$coefficients
      lambda_max = max*as.numeric(lambda)
      return(lambda)
    }
}

## function to merge ggraphs
pan_gg = function(ggs, #list of ggraphs to merge- will return the same length vector. Can be gGraphs or paths to ggraphs
                          fix = FALSE, #whether to fix the seqlengths of the ggraph so that there are not nodes on non standard chromosomes in one graphs
                          return_disjoin = FALSE,
                  seqlengths = NULL, #if seqlengths are provided and fix == TRUE, the these seqlengths will be used instead of hg_seqlengths
                  fixloose = TRUE, #Primarily for debugging purposes
                  fixmeta = TRUE,
                    cores = 1         #number of cores for generating the individual graphs
                  ) {
    ## add ids to the graphs before concatenating and fix if specified
    ggs = mclapply(1:length(ggs), function(x) {
        if(inherits(ggs[[x]],"character")) {
            new_gg = readRDS(ggs[[x]])
        } else if (inherits(ggs[[x]], "gGraph")) {
            new_gg = ggs[[x]]$copy
        }
        ##marking allows the gg_new_nodes function to stand on it's own and not require a value for which graph it is (graph.id in gg_new_nodes) - maybe that should be a requirement?
        new_gg$nodes$mark(ggraph_id = as.character(x))
        new_gg$edges$mark(ggraph_id = as.character(x))
        if(fix & is.null(seqlengths)) {
            new_gg$fix(seqlengths = gUtils::hg_seqlengths(chr = FALSE))
        } else if (fix & !is.null(seqlengths)) {
            new_gg$fix(seqlengths = gUtils::hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
        }
        return(new_gg)
    }, mc.cores = cores)
    ## concatenate and mark
    merged.gg = do.call(c,ggs)
    merged.gg$nodes$mark(node_disjoin_id = as.character(1:length(merged.gg$nodes)))
    merged.gg$edges$mark(edge_disjoin_id = as.character(1:length(merged.gg$edges)))
    ## browser()
    disjoin.gg = merged.gg$copy$disjoin()
    ## now get the individual graphs from the disjoin
    ggs.new.lst = mclapply(1:length(ggs), function(x) {
        ##run function to lift over edges to new node ids
        gg2 = gg_new_nodes(gg = ggs[[x]], disjoin.gg = disjoin.gg, merged.gg = merged.gg, fixloose = fixloose, fixmeta = fixmeta)
        return(gg2)
    }, mc.cores = cores)
    if(!return_disjoin) {
        return(ggs.new.lst)
    } else {
        return(list(ggs.new.lst,disjoin.gg))
    }
}


gg_new_nodes = function(gg,
                        disjoin.gg, 
                        merged.gg,
                        fixloose = TRUE,
                        fixmeta = TRUE) {
    disjoin.gg2 = disjoin.gg$copy
    ## add which parent graph based on the ggraph_id
    ggraph.id = gg$dt$ggraph_id %>% unique
    ## set node and ref count to zero
    disjoin.gg2$edges$mark(cn = NA)
    disjoin.gg2$nodes$mark(cn = NA)
    #######################################################
    ## get the copy number of the nodes from the parent ids
    merged.gg2 = merged.gg$copy
    cn.dt = merged.gg2$nodes$dt[ggraph_id == ggraph.id, .(node_disjoin_id,cn)]
    ## browser()
    ## get the ids and split into more rows to get the node.id for the node_disjoin_id in this sample
    disjoin.sub.dt = disjoin.gg2$nodes$dt[,.(node_disjoin_id,node.id)]
    disjoin.sub.dt = disjoin.sub.dt[, .(node_disjoin_id = unlist(strsplit(node_disjoin_id, ","))), by = node.id]
    disjoin.cn.dt = merge.data.table(cn.dt, disjoin.sub.dt, by = "node_disjoin_id", all.x = TRUE)
    ## mark the node cn
    node_order = disjoin.gg2$nodes[node.id %in% disjoin.cn.dt$node.id]$dt$node.id
    disjoin.cn.dt[, node_factor := factor(node.id, levels = node_order)] #have to get nodes in the same order
    disjoin.cn.dt = disjoin.cn.dt[order(node_factor),]
    disjoin.gg2$nodes[node.id %in% disjoin.cn.dt$node.id,]$mark(cn = disjoin.cn.dt$cn)
    #######################################################
    ## get the copy number of the edges from the parent ids
    cn.edge.dt = merged.gg2$edges$dt[ggraph_id == ggraph.id, .(edge_disjoin_id,cn)]
    ## get the ids for the edges in the new graph
    disjoin.sub.dt = disjoin.gg2$edgesdt[,.(edge_disjoin_id,edge.id)]
    disjoin.sub.dt = disjoin.sub.dt[, .(edge_disjoin_id = unlist(strsplit(edge_disjoin_id, ","))), by = edge.id] #put into long format by splitting the ids
    disjoin.edge.cn.dt = merge.data.table(cn.edge.dt, disjoin.sub.dt, by = "edge_disjoin_id", all.x = TRUE)
    ## mark the edge cn
    edge_order = disjoin.gg2$edges[edge.id %in% disjoin.edge.cn.dt$edge.id]$dt$edge.id
    disjoin.edge.cn.dt[, edge.id := factor(edge.id, levels = edge_order)]
    disjoin.edge.cn.dt = disjoin.edge.cn.dt[order(edge.id),]
    disjoin.gg2$edges[edge.id %in% disjoin.edge.cn.dt$edge.id]$mark(cn = disjoin.edge.cn.dt$cn)
    disjoin.gg2$edges[edge.id %in% disjoin.edge.cn.dt$edge.id]$mark(cn = disjoin.edge.cn.dt$cn)
    disjoin.gg2$edges[!(edge.id %in% disjoin.edge.cn.dt$edge.id)]$mark(cn = NA) #not sure why I have to add this-the disjoin.gg2$edges$mark(cn=0) at the top should have done this but there are edges that were not added that have values
    ## fix reference copy number before running loosefix
    ## edges.dt = disjoin.gg2$edgesdt[,.(type,cn,n1,n2,n1.side,n2.side)]
    ## add reference edges not in the original parent graph
    edge.ref.dt = disjoin.gg2$edges$dt[ggraph_id != ggraph.id, .(edge.id, n1,n2,n1.side,n2.side,cn, type)][type == "REF" & is.na(cn),]
    ##add copy number of n1 and n2 here
    edge.ref.dt = merge.data.table(edge.ref.dt,disjoin.gg2$nodes$dt[,.(n1 = node.id,n1.cn = cn)], by = "n1")
    edge.ref.dt = merge.data.table(edge.ref.dt,disjoin.gg2$nodes$dt[,.(n2 = node.id,n2.cn = cn)], by = "n2")
    edge.ref.dt[, cn := ifelse(n1.cn == n2.cn, n1.cn, NA)]
    ## browser()
    ## have to get it in the same order
    edge_order = disjoin.gg2$edges[edge.id %in% edge.ref.dt$edge.id]$dt$edge.id
    edge.ref.dt[, edge_factor := factor(edge.id, levels = edge_order)]
    edge.ref.dt = edge.ref.dt[order(edge_factor),]
    ## mark these reference edges
    disjoin.gg2$edges[edge.id %in% edge.ref.dt$edge.id]$mark(cn = edge.ref.dt$cn)
    ## browser()
    ## add the metadata back for nodes
    if(fixmeta) {
        og.nodes.gr = gg$nodes$gr
        new.nodes.gr = disjoin.gg2$nodes$gr
        names_meta_og.nodes = names(mcols(og.nodes.gr))
        names_meta_new.nodes = names(mcols(new.nodes.gr))
        meta_keep = c("loose.left", "loose.right", "cn", "start.ix", "end.ix", "eslack.in", "eslack.out", "loose", "edges.in", "edges.out", "tile.id", "cl", "loose.cn.left", "loose.cn.right", "node.id", "snode.id","index", "cluster","pcluster","ncluster")#,"col")
        meta_strip = names_meta_new.nodes[names_meta_new.nodes %in% names_meta_og.nodes]
        meta_strip = meta_strip[!(meta_strip %in% meta_keep)]
        ## add meta data to new nodes
        meta_keep_merged = c("node.id","snode.id","index", meta_strip)
        new.meta.nodes.dt = gr.val(new.nodes.gr[,meta_keep], og.nodes.gr[,meta_strip], meta_strip) %>% gr2dt(.) %>% .[,..meta_keep_merged]
        ## make sure this dt is in the same order as the nodes$node.id
        node_meta_order = disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id]$dt$node.id
        new.meta.nodes.dt[, node_factor := factor(node.id, levels = node_order)] #have to get nodes in the same order
        new.meta.nodes.dt = new.meta.nodes.dt[order(node_factor),]
        ## add missing meta so mark works with all when hard coded below
        ## strip meta in new graph and mark from old graph
        ## disjoin.gg$nodes$mark(
        ##                      for(x in meta_strip) {
        ##                          disjoin.gg2$annotate(colName = x, data = new.meta.nodes.dt[[x]], id = "node.id", class = "node")
        ##                          ## disjoin.gg2$nodes$mark(x = new.meta.nodes.dt[[x]])
        ##                          disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(x = new.meta.nodes.dt[[x]])
        ##                          ##
        ##                          do.call(disjoin.gg2$nodes$mark, list(x = NA))
        ##                      }
        ## ## alternative - swap breaks here and annotation does not work
        ## new.nodes.meta.gr = gr.val(new.nodes.gr[,meta_keep], og.nodes.gr[,meta_strip],  meta_strip)##  %>% gr2dt(.) %>% .[,..meta_keep_merged]
        ## disjoin.gg2$nodes$swap(new.nodes.meta.gr)
        ## ##
        ## mark.dt = new.meta.nodes.dt[,.(snode.id,index,epgap,simple,cpxdm,tyfonas,dm,bfb,chromothripsis,rigma,del,pyrgo,dup,tic,ecluster,qrpmix,qrpmin,qrppos,ggraph_id,node_factor)]
        meta_mark = c("node.id", "snode.id", "epgap", "simple", "cpxdm", "tyfonas", "dm", "bfb", "chromothripsis", "rigma", "pyrgo", "del", "dup", "tic", "ecluster", "qrpmix", "qrpmic", "qrpos")
        names_add_dt = meta_mark[!(meta_mark %in% names(new.meta.nodes.dt))]
        new.meta.nodes.dt[, (names_add_dt) := NA]
        ## hard code for now
        ## remove
        disjoin.gg2$nodes$mark(epgap = NA)
        disjoin.gg2$nodes$mark(simple = NA)
        disjoin.gg2$nodes$mark(cpxdm = NA)
        disjoin.gg2$nodes$mark(tyfonas = NA)
        disjoin.gg2$nodes$mark(dm = NA)
        disjoin.gg2$nodes$mark(bfb = NA)
        disjoin.gg2$nodes$mark(chromothripsis = NA)
        disjoin.gg2$nodes$mark(rigma = NA)
        disjoin.gg2$nodes$mark(pyrgo = NA)
        disjoin.gg2$nodes$mark(del = NA)
        disjoin.gg2$nodes$mark(dup = NA)
        disjoin.gg2$nodes$mark(tic = NA)
        disjoin.gg2$nodes$mark(qrpmix = NA)
        disjoin.gg2$nodes$mark(qrpos = NA)
        ##
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(epgap = new.meta.nodes.dt[["epgap"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(simple = new.meta.nodes.dt[["simple"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(cpxdm = new.meta.nodes.dt[["cpxdm"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(tyfonas = new.meta.nodes.dt[["tyfonas"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(dm = new.meta.nodes.dt[["dm"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(bfb = new.meta.nodes.dt[["bfb"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(chromothripsis = new.meta.nodes.dt[["chromothripsis"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(rigma = new.meta.nodes.dt[["rigma"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(pyrgo = new.meta.nodes.dt[["pyrgo"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(del = new.meta.nodes.dt[["del"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(dup = new.meta.nodes.dt[["dup"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(tic = new.meta.nodes.dt[["tic"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(qrpmix = new.meta.nodes.dt[["qrpmix"]])
        disjoin.gg2$nodes[node.id %in% new.meta.nodes.dt$node.id,]$mark(qrpos = new.meta.nodes.dt[["qrpos"]])
    }
    if(fixloose) {
        disjoin.gg3 = loosefix(disjoin.gg2$copy)
    } else {
        disjoin.gg3 = disjoin.gg2$copy
    }
    return(disjoin.gg3)
}

## function to compare cn of two ggs outputted by pan_gG
markdiff_pan_gg = function(ggs,
                          type = "pairwise", #can be pairwise if just comparing two samples, multi if more but have to have a reference sample
                          reference_sample = NULL  #if type = multi, then the reference sample is the sample that is subtracted from
                       ) {
    if(type == "pairwise") {
        ## nodes
        nodes.gr1 = ggs[[1]]$copy$nodes$gr
        nodes.gr2 = ggs[[2]]$copy$nodes$gr
        nodes.dt1 = gr2dt(nodes.gr1)
        nodes.dt2 = gr2dt(nodes.gr2)
        merge.nodes.dt = merge.data.table(nodes.dt1, nodes.dt2, by = "node.id", suffixes = c("_gg1","_gg2"), all = TRUE)
        ## Reduce(function(x,y) merge(x, y, by = "node.id", suffixes = c("_gg1","_gg2"), all = TRUE), list(nodes.dt1,nodes.dt2))
        diff_nodes.dt = merge.nodes.dt[,.(node.id, cn_gg1, cn_gg2)][cn_gg1 != cn_gg2,]
        ## edges
        edges.dt1 = ggs[[1]]$copy$edgesdt ## %>% gr ##.unlist
        edges.dt2 = ggs[[2]]$copy$edgesdt ## %>% gr ##.unlist
        merge.edges.dt = merge.data.table(edges.dt1, edges.dt2, by = "edge.id", suffixes = c("_gg1","_gg2"), all = TRUE)
        diff_edges.dt = merge.edges.dt[,.(edge.id, cn_gg1, cn_gg2)][cn_gg1 != cn_gg2,]

        ## nodes mark group_dif and color for nodes that are different cn
        diff_nodes.dt[, diff_cn := abs(cn_gg1 - cn_gg2)]
        diff_nodes.dt[, `:=`(group_diff = fifelse(diff_cn == 1, "one",
                                                  fifelse(diff_cn > 1 & diff_cn <= 5, "one_five",
                                                          fifelse(diff_cn > 5 & diff_cn <= 10, "five_ten", "ten_greater"))),
                             color = fifelse(diff_cn == 1, rgb(t(col2rgb("blue")), maxColorValue = 255),
                                             fifelse(diff_cn > 1 & diff_cn <= 5, rgb(t(col2rgb("orange")), maxColorValue = 255),
                                                     fifelse(diff_cn > 5 & diff_cn <= 10, rgb(t(col2rgb("red")), maxColorValue = 255), rgb(t(col2rgb("pink")), maxColorValue = 255)))))]
        ## edges mark group_dif and color for edges that are different cn
        diff_edges.dt[, diff_cn := abs(cn_gg1 - cn_gg2)]
        diff_edges.dt[, `:=`(group_diff = fifelse(diff_cn == 1, "one",
                                                  fifelse(diff_cn > 1 & diff_cn <= 5, "one_five",
                                                          fifelse(diff_cn > 5 & diff_cn <= 10, "five_ten", "ten_greater"))),
                             color = fifelse(diff_cn == 1, rgb(t(col2rgb("blue")), maxColorValue = 255),
                                             fifelse(diff_cn > 1 & diff_cn <= 5, rgb(t(col2rgb("orange")), maxColorValue = 255),
                                                     fifelse(diff_cn > 5 & diff_cn <= 10, rgb(t(col2rgb("red")), maxColorValue = 255), rgb(t(col2rgb("pink")), maxColorValue = 255)))))]
        ## mark these nodes with different colors
        gg1 = ggs[[1]]$copy
        gg2 = ggs[[2]]$copy
        ## nodes
        gg1$nodes$mark(col = "grey")
        gg1$nodes[node.id %in% diff_nodes.dt$node.id]$mark(col = diff_nodes.dt$color)
        gg2$nodes$mark(col = "grey")
        gg2$nodes[node.id %in% diff_nodes.dt$node.id]$mark(col = diff_nodes.dt$color)
        ## edges
        gg1$edges$mark(col = "grey")
        gg1$edges[edge.id %in% diff_edges.dt$edge.id]$mark(col = diff_edges.dt$color)
        gg2$edges$mark(col = "grey")
        gg2$edges[edge.id %in% diff_edges.dt$edge.id]$mark(col = diff_edges.dt$color)
        return(list(gg1,gg2))
    } else {
        stop("Not implemented")
    }
}


## not sure if I'll use these functions again
## function to subtract gg2 cn for edges and nodes from gg1
subtract_cn_gg = function(gg1, gg2, add_log_cn = TRUE) {
    gg1 = gg1$copy
    gg2 = gg2$copy
    ## nodes
    gg1.nodes.dt = gg1$nodes$dt[,.(node.id, cn)]
    gg2.nodes.dt = gg2$nodes$dt[,.(node.id, cn)]
    ## merge nodes
    merge.nodes.dt = merge.data.table(gg1.nodes.dt, gg2.nodes.dt, by = "node.id", suffixes = c("_gg1","_gg2"))
    merge.nodes.dt[, new_cn := (cn_gg1 - cn_gg2)]
    ## edges
    gg1.edges.dt = gg1$edges$dt[,.(edge.id, cn)]
    gg2.edges.dt = gg2$edges$dt[,.(edge.id, cn)]
    ## merge edges
    merge.edges.dt = merge.data.table(gg1.edges.dt, gg2.edges.dt, by = "edge.id", suffixes = c("_gg1","_gg2"))
    merge.edges.dt[, new_cn := (cn_gg1 - cn_gg2)]
    ## now mark
    gg1$nodes$mark(cn = merge.nodes.dt$new_cn)
    gg1$edges$mark(cn = merge.edges.dt$new_cn)
    if(add_log_cn) {
        merge.nodes.dt[new_cn != 0, log_new_cn := log10(abs(new_cn))]
        merge.nodes.dt[new_cn < 0, log_new_cn := -log_new_cn]
        merge.nodes.dt[new_cn == 0, log_new_cn := 0]
        gg1$nodes$mark(log_new_cn = merge.nodes.dt$log_new_cn)
    }
    return(gg1)
}

## function to return a gg after detecting a minimum cn within a grange
min_gg_gt = function(gg, gr, pad = 1) {
    gg = gg$copy
    nodes.gr = gg$nodes$gr %&% gr
    if(length(nodes.gr) == 0)
        stop("No nodes overlap the specific gr. (chr) must be the same between both")
    yf = gg$meta$y.field
    miny = mcols(nodes.gr)[,yf] %>% min
    maxy = mcols(nodes.gr)[,yf] %>% max
    gg$set(y0 = (miny-pad))
    gg$set(y1 = (maxy-pad))
    return(gg)
}

## calculates the euclidean distance and sum square (euclidean squared) in cn for events-will likely add other things later
distance_cn = function(patient,
                       pair_names = as.character(1:length(complexes)),
                       complexes,  #ggraphs to get events from
                       merged_ggs, #output of pan_gg
                       tile_bins = TRUE, #TRUE will tile so the bins are all the same size for comparing cn, if false, will just use the nodes from the merged_ggs
                       cores = 2,
                       event_pairs = NULL, #only get events from these pair_names
                       events_keep = c("tyfonas", "bfb", "dm", "cpxdm")) { #can be euclidean
    lengths_check = c(length(pair_names),length(complexes), length(merged_ggs)) %>% unique
    if(length(lengths_check) != 1)
        stop("pair_names, complexes, and merged_ggs have to have the same length")
    message("Starting ",patient)
    ## get the events
    ## make a pairs table
    pairs.dt = data.table(pair = pair_names, complex = complexes, merged_gg = merged_ggs)
    setkey(pairs.dt, pair)
    ## events.dt = mclapply(setNames(nm = pair_names), function(x) {
    events.dt = mclapply(setNames(nm = pairs.dt$pair), function(pr) {
        gg = readRDS(pairs.dt[pr,complex])
        if(nrow(gg$meta$events) == 0) {
            return(NULL)
        }
        events.dt = gg$meta$events[type %in% events_keep,]
        events.dt[, sample := pr]
        return(events.dt)
    }, mc.cores = cores)
    events.dt = rbindlist(events.dt, fill = TRUE)
    if(!is.null(event_pairs) & (nrow(events.dt) > 0)) {
        events.dt = events.dt[sample %in% event_pairs,]
    }
    if(nrow(events.dt) == 0) {
        if(!is.null(event_pairs)) {
            warning("No events in the complexes for event_pairs. Returning NULL.")
        } else {
            warning("No events in any of the complexes. Returning NULL.")
        }
        return(NULL)
    }
    ## footprints-will duplicate if multiple sames have the same events
    ## if exactly the same I do not want to duplicate
    ## unique_fp.dt = events.dt[,.(type,footprint)] %>% unique
    unique_fp.dt = events.dt[,.(type,footprint, sample)]
    unique_fp.dt[, event_numb := paste0(type,"_",1:.N), by = c("type","sample")]
    unique_fp.dt[, event_numb_pat := paste0(type,"_",1:.N), by = c("type","sample")]
    unique_fp.dt[, event_row := 1:.N]
    ## now use the merged ggraph as ggs
    ggs = mclapply(setNames(nm = pairs.dt$pair), function(pr) {
        gg = readRDS(pairs.dt[pr, merged_gg])
        return(gg)
    }, mc.cores = cores)
    all_amps_fp.lst = mclapply(1:nrow(unique_fp.dt), function(x) {
        amp_fp.gr = parse.gr(unique_fp.dt$footprint[x])
        amp_fp.dt = as.data.table(amp_fp.gr)
        amp_fp.dt[, event_row := x,]
        amp_fp.dt = merge.data.table(amp_fp.dt, unique_fp.dt[event_row == x,.(sample, event_row,event_numb,footprint)], by = "event_row")
        return(amp_fp.dt)
    }, mc.cores = cores)
    all_amps_fp.dt = rbindlist(all_amps_fp.lst, fill = TRUE)
    all_amps_fp.dt[, sample_event := sample]
    all_amps_fp.gr = GRanges(all_amps_fp.dt, seqlengths = hg_seqlengths())
    ## now get the cn of this region-using the merged ggraphs so they have the same nodes
    cn.lst = mclapply(setNames(nm = names(ggs)), function(x) {
        ## nodes.dt = (ggs[[x]]$nodes$gr %&% all_amps_fp.gr) %>% gr2dt
        nodes.gr = (ggs[[x]]$nodes$gr %&% all_amps_fp.gr)
        ## have to do gr.val for each event - I know this is lazy with all the loopsssss
        unique_events.dt = unique(all_amps_fp.dt[,.(sample,event_numb)])
        nodes.lst = lapply(1:nrow(unique_events.dt), function(y) {
            sub.event.dt = unique_events.dt[y,]
            sub.amp.gr = all_amps_fp.gr %Q% (event_numb == sub.event.dt$event_numb & sample == sub.event.dt$sample)
            nodes.dt = gr.val(nodes.gr, sub.amp.gr, "event_numb") %>% gr2dt
            cols_keep = c("seqnames", "start", "end", "node.id", "cn", "event_numb", events_keep)
            cols_keep = cols_keep[cols_keep %in% names(nodes.dt)]
            nodes.dt = nodes.dt[, ..cols_keep]
            nodes.dt = nodes.dt[event_numb == sub.event.dt$event_numb,] #subset to nodes with this event
            nodes.dt[,sample_event := sub.event.dt$sample]
            if(tile_bins) {
                nodes.gr = GRanges(nodes.dt)
                nodes.gr2 = gr.tile(nodes.gr, width = 1e4)
                nodes.gr2 = gr.val(nodes.gr2,nodes.gr, cols_keep[!(cols_keep %in% c("seqnames","start","end"))])
                nodes.dt = gr2dt(nodes.gr2)
                nodes.dt[, c("strand","width","query.id") := NULL]
                nodes.dt[,sample_event := sub.event.dt$sample]
                return(nodes.dt)
            } else {
                return(nodes.dt)
            }
        })
        nodes.dt = rbindlist(nodes.lst, fill = TRUE)
        return(nodes.dt)
    }, mc.cores = cores)
    ## merge the cn
    add_suffix = paste0("_",names(cn.lst))
    cols_keep = c("seqnames", "start", "end", "node.id", "event_numb", "sample_event")#, events_keep)
    cn.dt = cn.lst[[1]]
    ## cn.dt = cn.dt[,..cols_keep]
    ## add suffix to this first sample and every sample
    sample1_suffix = names(cn.dt)[!(names(cn.dt) %in% cols_keep)]
    names(cn.dt)[names(cn.dt) %in% sample1_suffix] = paste0(sample1_suffix,add_suffix[1])
    names_keep = c(cols_keep, grep("cn",names(cn.dt),value = TRUE))
    cn.dt = cn.dt[,..names_keep]
    ## merge with loop and add specific suffix
    for (i in 2:length(cn.lst)) {
        cn.new.dt = cn.lst[[i]]
        suffix = add_suffix[i]
        new_sample_suffix = names(cn.new.dt)[!(cn.new.dt %in% cols_keep)]
        names(cn.new.dt)[names(cn.new.dt) %in% sample1_suffix] = paste0(sample1_suffix,add_suffix[i])
        names_keep = c(cols_keep, grep("cn", names(cn.new.dt), value = TRUE))
        cn.new.dt = cn.new.dt[,..names_keep]
        cn.dt = merge.data.table(cn.dt, cn.new.dt, by = cols_keep, all = TRUE, suffixes = c("", "")) #will give warning if something wrong here
    }
    ## calculate euclidean distance for the nodes in each event
    cn.dt[, sample := patient]
    names(cn.dt) = gsub("-","_",names(cn.dt))
    unique_events.dt = unique(cn.dt[,.(event_numb, sample_event)])
    ## final.lst = mclapply(unique(cn.dt$event_numb), function(x) {
    final.lst = mclapply(1:nrow(unique_events.dt), function(x) {
        sub.event.dt = unique_events.dt[x,]
        sub.cn.dt = cn.dt[event_numb == sub.event.dt$event_numb & sample_event == sub.event.dt$sample_event,]
        cn_cols = grep("cn.", names(sub.cn.dt),value = TRUE)
        dist.dt = sub.cn.dt[,..cn_cols]
        euclidean_distances.mat = as.matrix(dist(t(dist.dt), method = "euclidean"))
        diag(euclidean_distances.mat) = NA #diagonal is all 0 but we want max and min when more than two samples
        final.dt = sub.cn.dt[,.(sample,event_numb,sample_event)] %>% unique
        final.dt[, min_euc_dist := min(euclidean_distances.mat,na.rm = TRUE)]
        final.dt[, max_euc_dist := max(euclidean_distances.mat,na.rm = TRUE)]
        final.dt[, mean_euc_dist := mean(euclidean_distances.mat,na.rm = TRUE)]
        final.dt = merge.data.table(final.dt,unique(all_amps_fp.dt[,.(event_numb,footprint, sample_event)]), by = c("event_numb","sample_event"), all.x = TRUE)
        final.dt[, min_square_dist := min_euc_dist^2]
        final.dt[, max_square_dist := max_euc_dist^2]
        final.dt[, mean_square_dist := mean_euc_dist^2]
        ## add width of the footprint to filter on
        final.dt[, width := sum(parse.gr(final.dt$footprint)@ranges@width)]
    }, mc.cores = cores)
    final.dt = rbindlist(final.lst, fill = TRUE)
    final.dt[, N_samples := length(ggs)]
    return(final.dt)
}


## ## old function to merge ggraphs
## merge_gg = function(gg1,
##                     gg2,
##                     fix = FALSE, #whether to fix the seqlengths of the ggraph so that there are not nodes on non standard chromosomes in one graphs
##                     seqlengths = NULL #if seqlengths are provided and fix == TRUE, the these seqlengths will be used instead of hg_seqlengths
##                     ) {
##     ## make sure not to modify the original object
##     gg1 = gg1$copy
##     gg2 = gg2$copy
##     ## fix the seqlengths if specified so there are not nodes on extra chromosomes that will only be added for one sample
##     if(fix & is.null(seqlengths)) {
##         gg1$fix(seqlengths = hg_seqlengths(chr = FALSE))
##         gg2$fix(seqlengths = hg_seqlengths(chr = FALSE))
##     } else if (fix & !is.null(seqlengths)) {
##         gg1$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
##         gg2$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
##     }
##     ## mark the nodes and edges so cn is carried to the new object and can be set
##     gg1$nodes$mark(cn_gg1 = gg1$nodes$dt$cn)
##     gg1$edges$mark(cn_gg1 = gg1$edges$dt$cn)
##     gg2$nodes$mark(cn_gg2 = gg2$nodes$dt$cn)
##     gg2$edges$mark(cn_gg2 = gg2$edges$dt$cn)
##     ## merge the ggraphs
##     merged.gg = c(gg1,gg2)
##     ## disjoin to get merged nodes
##     merged.gg2 = merged.gg$copy$disjoin()
##     ## select only parts of the graph for each gg
##     browser()
##     gg1.only.gg = merged.gg2[grepl("gGraph2",parent.graph)]
##     gg2.only.gg = merged.gg2[grepl("gGraph1",parent.graph)]
##     ## had different columns for y.field
##     ## gg1.only.gg$set(y.field = "cn_gg1")
##     ## gg2.only.gg$set(y.field = "cn_gg2")
##     return(list(gg1.only.gg, gg2.only.gg))
## }


#modified the amp function to return all of the stats used to separate the junctions
amp_stats = function(gg, jcn.thresh = 8, cn.thresh = 2, fbi.cn.thresh = 0.5,  
               n.jun.high.bfb.thresh = 26, n.jun.high.dm.thresh = 31, width.thresh = 1e5, 
               fbi.width.thresh = 1e5, mc.cores = 1, mark = TRUE, mark.col = 'purple', 
               mark.nos = FALSE, min.nodes = 3, min.jun = 2) {
    if (mark.nos) {
        gg$nodes$mark(nos = as.integer(NA))
        gg$edges$mark(nos = as.integer(NA))
    }

    gg$nodes$mark(cpxdm = as.integer(NA))
    gg$edges$mark(cpxdm = as.integer(NA))
    gg$nodes$mark(tyfonas = as.integer(NA))
    gg$edges$mark(tyfonas = as.integer(NA))
    gg$nodes$mark(dm = as.integer(NA))
    gg$edges$mark(dm = as.integer(NA))
    gg$nodes$mark(bfb = as.integer(NA))
    gg$edges$mark(bfb = as.integer(NA))
    gg$edges$mark(fbi = gg$edges$class == 'INV-like' & gg$edges$span < fbi.width.thresh)
    gg$set(amp = data.table())
    ploidy = gg$nodes$dt[!is.na(cn), sum(cn*as.numeric(width))/sum(as.numeric(width))]
    keep = (gg$nodes$dt$cn/ploidy) > cn.thresh
    gg$clusters(keep)
    if (!any(!is.na(gg$nodes$dt$cluster)))
        return(gg)

  tiny = gg$edges$mark(tiny = gg$edges$dt$class %in% c('DEL-like', 'DUP-like') & gg$edges$span <1e4)
  ucl = gg$nodes$dt[!is.na(cluster), .(wid = sum(width)), by = cluster][wid > width.thresh, cluster] %>% sort

  amps = mclapply(ucl, function(cl, ploidy) {
    cl.nodes = gg$nodes[cluster == cl]
    cl.edges = cl.nodes$edges[type == "ALT" & tiny == FALSE]
    if (!length(cl.edges)) 
      return(NULL)
    if (!length(cl.edges)) 
      return(NULL)
    data.table(cluster = cl,
               nodes = paste(cl.nodes$dt$node.id,
                             collapse = ","),
               edges = paste(cl.edges$dt$edge.id, 
                             collapse = ","),
               fbi.cn = 2 * sum(cl.edges$dt[fbi == TRUE, sum(cn)]),
               n.jun = length(cl.edges),
               n.jun.high = sum(cl.edges$dt[, sum(cn > 3)]), 
               max.jcn = max(c(0, cl.edges$dt$cn)),
               max.cn = max(cl.nodes$dt$cn),
               footprint = paste(gr.string(cl.nodes$footprint),
                                 collapse = ","))
  }, ploidy, mc.cores = mc.cores) %>% rbindlist

  if (nrow(amps)) {
      if (!mark.nos) {
          amps = amps[max.jcn >= jcn.thresh,]
          ##amps[max.jcn >= jcn.thresh, ]
      } else {
          ## keep only clusters with a sufficient number of nodes but don't filter by jcn
          amps = amps[max.jcn >= jcn.thresh | 
                      (n.jun >= min.jun &
                       (strsplit(nodes, ",") %>% lapply(length) %>% unlist) >= min.nodes),]
      }
  }
  ## implementing decision tree in https://tinyurl.com/srlbkh2
  if (nrow(amps)) {
      ## order / rename and mark
      gg$set(amp = amps)

      ## call and mark event types
      amps[, type := ifelse(
                 max.jcn < jcn.thresh,
                 "nos",
                     ifelse(
                         ## few high copy junctions, high fbi cn -> BFB, otherwise tyfonas
                         fbi.cn / max.cn >= fbi.cn.thresh,
                     ifelse(n.jun.high < n.jun.high.bfb.thresh, 
                            'bfb',
                            'tyfonas'),
                     ## few high copy junctions, low fbi cn -> DM, otherwise CPXDM
                     ifelse(n.jun.high >= n.jun.high.dm.thresh, 
                            'cpxdm',     
                            'dm')
                     )
             )]
      
      amps[, ev.id := 1:.N, by = type]

      ## unlist node and edge ids and map back to type and ev label
      nodelist = strsplit(amps$nodes, ',') %>% lapply(as.integer) %>% dunlist
      edgelist = strsplit(amps$edges, ',') %>% lapply(as.integer) %>% dunlist
      nodelist = cbind(nodelist, amps[nodelist$listid, .(type, ev.id)])
      edgelist = cbind(edgelist, amps[edgelist$listid, .(type, ev.id)])
      return(amps)
  }
}

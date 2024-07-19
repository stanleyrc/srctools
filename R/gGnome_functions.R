## function to merge ggraphs
pan_gg = function(ggs, #list of ggraphs to merge- will return the same length vector. Can be gGraphs or paths to ggraphs
                          fix = FALSE, #whether to fix the seqlengths of the ggraph so that there are not nodes on non standard chromosomes in one graphs
                          return_disjoin = FALSE,
                  seqlengths = NULL, #if seqlengths are provided and fix == TRUE, the these seqlengths will be used instead of hg_seqlengths
                  fixloose = TRUE, #Primarily for debugging purposes
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
            new_gg$fix(seqlengths = hg_seqlengths(chr = FALSE))
        } else if (fix & !is.null(seqlengths)) {
            new_gg$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
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
        gg2 = gg_new_nodes(gg = ggs[[x]], disjoin.gg = disjoin.gg, merged.gg = merged.gg, fixloose = fixloose)
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
                        fixloose = TRUE) {
    disjoin.gg2 = disjoin.gg$copy
    ## add which parent graph based on the ggraph_id
    ggraph.id = gg$dt$ggraph_id %>% unique
    ## browser()
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
    ## have to get it in the same order
    edge_order = disjoin.gg2$edges[edge.id %in% edge.ref.dt$edge.id]$dt$edge.id
    edge.ref.dt[, edge_factor := factor(edge.id, levels = edge_order)]
    edge.ref.dt = edge.ref.dt[order(edge_factor),]
    ## mark these reference edges
    disjoin.gg2$edges[edge.id %in% edge.ref.dt$edge.id]$mark(cn = edge.ref.dt$cn)
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


## ## function to merge ggraphs
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

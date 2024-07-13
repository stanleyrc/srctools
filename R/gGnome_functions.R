## function to merge ggraphs
merge_gg_nodes = function(ggs, #list of ggraphs to merge- will return the same length vector
                    fix = FALSE, #whether to fix the seqlengths of the ggraph so that there are not nodes on non standard chromosomes in one graphs
                    seqlengths = NULL, #if seqlengths are provided and fix == TRUE, the these seqlengths will be used instead of hg_seqlengths
                    cores = 1         #number of cores for generating the individual graphs
                    ) {
    ## fix the seqlengths if specified so there are not nodes on extra chromosomes that will only be added for one sample
    if(fix & is.null(seqlengths)) {
        gg1$fix(seqlengths = hg_seqlengths(chr = FALSE))
        gg2$fix(seqlengths = hg_seqlengths(chr = FALSE))
    } else if (fix & !is.null(seqlengths)) {
        gg1$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
        gg2$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
    }
    number_ggs = length(ggs)
    ## add ids to the graphs before concatenating
    ggs = mclapply(1:length(ggs), function(x) {
        new_gg = ggs[[x]]$copy
        new_gg$nodes$mark(ggraph_id = as.character(x))
        new_gg$edges$mark(ggraph_id = as.character(x))
        return(new_gg)
    }, mc.cores = cores)
    ## concatenate and mark
    merged.gg = do.call(c,ggs)
    merged.gg$nodes$mark(node_disjoin_id = as.character(1:length(merged.gg$nodes)))
    merged.gg$edges$mark(edge_disjoin_id = as.character(1:length(merged.gg$edges)))
    disjoin.gg = merged.gg$copy$disjoin()

    ## now get the individual graphs from the disjoin
    ggs.new.lst = mclapply(1:number_ggs, function(x) {
        ##run function to lift over edges to new node ids
        gg2 = gg_new_nodes(gg = ggs[[x]], disjoin.gg = disjoin.gg, merged.gg = merged.gg)
        return(gg2)
    }, mc.cores = cores)
    return(ggs.new.lst)
}


gg_new_nodes = function(gg,
                        disjoin.gg, 
                        merged.gg) {
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
    ## get the ids and split into more rows to get the node.id for the node_disjoin_id in this sample
    disjoin.sub.dt = disjoin.gg2$nodes$dt[,.(node_disjoin_id,node.id)]
    disjoin.sub.dt = disjoin.sub.dt[, .(node_disjoin_id = unlist(strsplit(node_disjoin_id, ","))), by = node.id]
    disjoin.cn.dt = merge.data.table(cn.dt, disjoin.sub.dt, by = "node_disjoin_id", all.x = TRUE)
    ## mark the cn
    disjoin.cn.dt = disjoin.cn.dt[order(node.id)]
    disjoin.gg2$nodes[node.id %in% disjoin.cn.dt$node.id,]$mark(cn = disjoin.cn.dt$cn)
    #######################################################
    ## get the copy number of the edges from the parent ids
    cn.edge.dt = merged.gg2$edges$dt[ggraph_id == ggraph.id, .(edge_disjoin_id,cn)]
    ## get the ids for the edges in the new graph
    disjoin.sub.dt = disjoin.gg2$edgesdt[,.(edge_disjoin_id,edge.id)]
    disjoin.sub.dt = disjoin.sub.dt[, .(edge_disjoin_id = unlist(strsplit(edge_disjoin_id, ","))), by = edge.id] #put into long format by splitting the ids
    disjoin.edge.cn.dt = merge.data.table(cn.edge.dt, disjoin.sub.dt, by = "edge_disjoin_id", all.x = TRUE)
    ## mark the edge cn
    disjoin.edge.cn.dt = disjoin.edge.cn.dt[order(edge.id),]
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
    ## mark these reference edges
    disjoin.gg2$edges[edge.id %in% edge.ref.dt$edge.id]$mark(cn = edge.ref.dt$cn)
    disjoin.gg3 = loosefix(disjoin.gg2$copy)
    return(disjoin.gg3)
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

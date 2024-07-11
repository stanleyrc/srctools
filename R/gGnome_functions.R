
## function to merge ggraphs
merge_gg = function(gg1,
                    gg2,
                    fix = FALSE, #whether to fix the seqlengths of the ggraph so that there are not nodes on non standard chromosomes in one graphs
                    seqlengths = NULL #if seqlengths are provided and fix == TRUE, the these seqlengths will be used instead of hg_seqlengths
                    ) {
    ## make sure not to modify the original object
    gg1 = gg1$copy
    gg2 = gg2$copy
    ## fix the seqlengths if specified so there are not nodes on extra chromosomes that will only be added for one sample
    if(fix & is.null(seqlengths)) {
        gg1$fix(seqlengths = hg_seqlengths(chr = FALSE))
        gg2$fix(seqlengths = hg_seqlengths(chr = FALSE))
    } else if (fix & !is.null(seqlengths)) {
        gg1$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
        gg2$fix(seqlengths = hg_seqlengths(seqlengths = seqlengths, chr = FALSE))
    }
    ## mark the nodes and edges so cn is carried to the new object and can be set
    gg1$nodes$mark(cn_gg1 = gg1$nodes$dt$cn)
    gg1$edges$mark(cn_gg1 = gg1$edges$dt$cn)
    gg2$nodes$mark(cn_gg2 = gg2$nodes$dt$cn)
    gg2$edges$mark(cn_gg2 = gg2$edges$dt$cn)
    ## merge the ggraphs
    merged.gg = c(gg1,gg2)
    ## disjoin to get merged nodes
    merged.gg2 = merged.gg$copy$disjoin()
    ## select only parts of the graph for each gg
    gg1.only.gg = merged.gg2[grepl("gGraph2",parent.graph)]
    gg1.only.gg$set(y.field = "cn_gg1")
    gg2.only.gg = merged.gg2[grepl("gGraph1",parent.graph)]
    gg2.only.gg$set(y.field = "cn_gg2")
    return(list(gg1.only.gg, gg2.only.gg))
}

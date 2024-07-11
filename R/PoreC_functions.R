## function for getting concatemers that interescet another gr with multiple ranges- returned gr has all concatemers interseting all of the ranges in intersect.gr
intersect_parq = function(parquet.gr,
                      intersect.gr,
                      idx_col = "read_idx",
                      return_type = "gr", #can be ids or gr
                      cores = 1) {
    parquet.sub.gr = parquet.gr %&% intersect.gr
    if(class(mcols(parquet.sub.gr)[[idx_col]]) != "character") {
        warning("idx_col is not a character. Will convert to a character")
        orig_class = class(mcols(parquet.sub.gr)[[idx_col]])
        mcols(parquet.sub.gr)[[idx_col]] = as.character(mcols(parquet.sub.gr)[[idx_col]])
    } else {
        orig_class = "character"
    }
    ## get ids intersecting each range
    intersect.gr = gr.val(intersect.gr, parquet.sub.gr, idx_col, mc.cores = cores)
    read_idx.lst = mclapply(1:length(intersect.gr), function(x) {
        read_idx = tstrsplit(mcols(intersect.gr[x])[[idx_col]],", ") %>% unlist
        dt1 = data.table(range = x, ids = read_idx)
        ## return(list(read_idx))
        return(dt1)
    }, mc.cores = cores)
    ids.dt = rbindlist(read_idx.lst) %>% unique
    ids.dt[, count_id := .N, by = "ids"]
    ## get the ids that are in all of the granges
    ids.common.dt = ids.dt[count_id >= length(intersect.gr),]
    common_ids = ids.common.dt$ids %>% unique
    ## now return the ids or the granges
    if(return_type == "ids") {
        if(orig_class == "integer") {
            common_ids = as.integer(common_ids)
        }
        return(common_ids)
    } else if(return_type == "gr") {
        if(orig_class == "integer") {
            common_ids = as.integer(common_ids)
        }
        final.gr = parquet.gr[mcols(parquet.gr)[[idx_col]] %in% common_ids]
        return(final.gr)
    }
}


##jamesons cocount function with a unique
jameson_cocount = function(events, bins = disjoin(events), by = names(values(events))[1], weight = NULL, frac = FALSE, full = FALSE, fill = 0, na.rm = TRUE)
{
  if (length(events)==0)
    return(gM(gr = bins, full = full, fill = fill, agg.fun = agg.fun, na.rm = na.rm))
  
  if (is.na(by))
    stop('by must be specified and a metadata column of events GRanges')
           
  if (!(by %in% names(values(events))))
    stop('by must be a metadata column of events GRanges')

  events$group = values(events)[[by]]
  if (!is.null(weight))
    tmp = gr2dt(events[, c("group", weight)] %*% bins)
  else
    tmp = gr2dt(events[, c("group")] %*% bins)

  if (nrow(tmp)>0)
    tmp = tmp[!is.na(group), ]

  if (!nrow(tmp))
    {
      if (length(bins)>0)
        return(gM(bins))
      else
        return(gM())
    }

  if (!is.null(weight))
    {
      tmp$weight = tmp[[weight]]
    }
  else if (frac)
    tmp[, weight := width/sum(width), by = group]
  else
    tmp[, weight := 1, by = group]
  
  tmp = tmp[, .(bid = subject.id, group = as.integer(factor(group)), weight)]

####what if we dedupe here
  ##browser()
  tmp = unique(tmp, by=c('bid','group'))
  
  ## sum weights inside bin pairs that share a group
  dat = merge(tmp, tmp, by = c('group'), allow.cartesian = TRUE)[, .(value = sum(weight.x * weight.y)), by = .(i = bid.x, j = bid.y)]
###get rid of self edges how about that
##  dat = dat[i != j]
  
  return(gM(bins, dat, full = full, fill = fill, na.rm = na.rm, agg.fun = sum))
}

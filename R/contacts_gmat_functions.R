#function for getting the contacts of an event with other regions
contacts_events = function(pair, res = 1e6, gg, hic, hic_type = ".hic",gr_seqlengths = hg_seqlengths(),straw.cores=1, event.cores = 4,event_types = c("bfb","chromoplexy","chromothripsis","del","cpxdm","dm","dup","pyrgo","rigma","simple","tic","tyfonas"),add_event_event_contact_labels=FALSE, annotate_sum = TRUE, annotate_cn = TRUE, fill_missing_bins = TRUE) {
    ##maybe add ability to pad to events later??,pad = NULL
    if(class(gg) == "character") {
        gg = readRDS(gg)
    }
    ## if("gGraph" %in% class(gg)) {
    ## }
    ##read in hic
    if(hic_type == ".hic") {
        hic = GxG::straw(hic = hic,
                         res = res,
                         gr = gr_seqlengths,
                         mc.cores = straw.cores)
        message("Done reading in .hic")
    } else if(hic_type == ".mcool") {
        hic = GxG::cooler(file = hic,
                    res = res,
                    gr = gr_seqlengths,
                    mc.cores = straw.cores)
        message("Done reading in .mcool")
    } else if (hic_type == "gm") {
    }
    sum_contacts = hic$dat$value %>% sum
    if(annotate_sum) {
        i.contacts = hic$dat[,.(i,value)]
        j.contacts = hic$dat[,.(j,value)] %>% setnames(.,c("i","value"))
        ## i.contacts = events_contacts.dt[,.(i,value)]
        ## j.contacts = events_contacts.dt[,.(j,value)] %>% setnames(.,c("i","value"))
        contacts_i_j = rbind(i.contacts,j.contacts)
        contacts_i_j[, sum_bin := sum(value), by = "i"]
        contacts_i_j2 = contacts_i_j[,.(i,sum_bin)] %>% unique
    }
    hic.gr = hic$gr %>% gr.nochr()
    hic.gr$i = 1:length(hic.gr)
    gg = gGnome::refresh(gg)
    events.dt = gg$meta$events
    events.dt = events.dt[type %in% event_types,][,.(type,footprint)] #complex_events
    ##get granges with all events labeled
    ## events_unique.dt =  events.dt[,.(type,footprint)] %>% unique
    events_unique.dt =  unique(events.dt)
    if(nrow(events_unique.dt) != 0) {
        events_unique.dt[,event.number := paste0(type,"_",1:.N), by = "type"] #add index to events to separate if multiple of same event
                                        #add events to hic gr
        events_contacts.lst = mclapply(1:nrow(events_unique.dt), function(x) {
            gr.event1 = parse.gr(events_unique.dt[x,]$footprint)
            dt.event1 = as.data.table(gr.event1)
            dt.event1[,seqnames := gsub("chr","",seqnames)]
            gr.event1 = GRanges(dt.event1,seqlengths = seqlengths(gr.nochr(gr_seqlengths)))  %>% sortSeqlevels()  %>% sort()
                                        #gr.event1 = gr.fix(gr.event1,hic.gr)
            gr.event1$event_type = events_unique.dt[x,]$event.number
            ## gr.event1$footprint = events_unique.dt[x,]$footprint
            ## gr.event1$type = events_unique.dt[x,]$type
            ## hic.gr = gr.val(hic.gr,gr.event1,c("event_type","footprint","type"))
            hic.gr = gr.val(hic.gr,gr.event1,"event_type")
            event.filter.gr = hic.gr %Q% (!is.na(event_type))
            if(length(event.filter.gr) == 0) {
                return(NULL)
            }
            event.filter.gr = hic.gr %Q% (event_type != "")
                                        #    hic.gr = hic.gr %Q% (event_type != "")
            ## hic_events.dt = as.data.table(hic.gr)[,.(seqnames,start, end, i,type, event_type,footprint)]
            hic_events.dt = as.data.table(hic.gr)[,.(seqnames,start, end, i,event_type)]
            dat1 = hic$dat[i %in% event.filter.gr$i | j %in% event.filter.gr$i,] #subset data to only ones in i or j in the event
            ## add missing bins with zeroes
            if(fill_missing_bins) {
                ##all combos of i and j
                all_combos.dt = expand.grid(i = 1:max(hic.gr$i),j = event.filter.gr$i) %>% as.data.table
                sub.dat = dat1[,.(i,j)]
                sub.dat[, i_j := paste0(i,"_",j)]
                all_combos.dt[, i_j := paste0(i,"_",j)]
                missing_combos.dt = all_combos.dt[!(i_j %in% sub.dat$i_j),]
                missing_combos.dt[, i_j := NULL]
                missing_combos.dt[, value := 0]
                dat1 = rbind(dat1,missing_combos.dt,fill = TRUE)
            }
            dat1[i %in% event.filter.gr$i,event_type.i := gr.event1$event_type[1]]
            dat1[j %in% event.filter.gr$i,event_type.j := gr.event1$event_type[1]]
            dat1[!is.na(event_type.i) & is.na(event_type.j), c("i","j","event_type.i","event_type.j") := .(j,i,event_type.j,event_type.i)] #move all of the contacts with the event to i and move the event to j
            ##add coordinates for both i and j
            dat1 = merge.data.table(dat1, hic_events.dt[,.(seqnames,start,end,i)], by = "i")
            names(dat1) = gsub("seqnames","seqnames.i",names(dat1)) %>% gsub("start","start.i",.) %>% gsub("end","end.i",.)
            dat1 = merge.data.table(dat1, hic_events.dt[,.(seqnames,start,end,i)], by.x = "j", by.y = "i")
            names(dat1) = gsub("seqnames$","seqnames.j",names(dat1)) %>% gsub("start$","start.j",.) %>% gsub("end$","end.j",.)
            ##add whether the contact is intra chromosomal
            dat1[, intra_chrom := ifelse(seqnames.i == seqnames.j,TRUE,FALSE)]
            ##add other meta data
            dat1[, event.numb := events_unique.dt[x,]$event.numb]
            dat1[, footprint := events_unique.dt[x,]$footprint]
            dat1[, type := events_unique.dt[x,]$type]
            ##add whether i and j are both on a chromosome with the event
            seqnames_fp = parse.gr(unique(dat1$footprint))@seqnames@values %>% gsub("chr","",.)
            dat1[, chrom_in_event := ifelse((seqnames.i %in% seqnames_fp),TRUE,FALSE)]

            return(dat1)
        },mc.cores = event.cores)#  %>% do.call(c,.)  %>% sortSeqlevels() %>% sort()
        events_contacts.dt = rbindlist(events_contacts.lst)
        if(add_event_event_contact_labels) {
            ##add all event types for i and j, j already added but might want to look at contacts between different types of events later - i event needs to added where event_type.i is not na
            events.indices.i = events_contacts.dt[!is.na(event_type.i),.(i,event_type.i)] %>% unique %>% setnames(.,c("i","event_type"))
            events.indices.j = events_contacts.dt[!is.na(event_type.j),.(j,event_type.j)] %>% unique %>% setnames(.,c("i","event_type"))
            events.indices = rbind(events.indices.i,events.indices.j) %>% unique
            events.indices = events.indices[order(i),]
            events.indices = events.indices[, .(event.numb_all = paste0(event_type, collapse = ",")), by = i]
            #merge to i
            events_contacts.dt = merge.data.table(events_contacts.dt,events.indices, by = "i", all.x = TRUE)
            events_contacts.dt[is.na(event_type.i) & !is.na(event.numb_all), event_type.i := event.numb_all]
            events_contacts.dt[!is.na(event.numb_all) & (event.numb_all != event_type.j), event.numb_all := paste0(event.numb,",",event.numb_all)]
            events_contacts.dt[is.na(event.numb_all), event.numb_all := event_type.j]
            events_contacts.dt[, chrom_in_event := as.character(chrom_in_event)]
            events_contacts.dt[event_type.i == event_type.j, chrom_in_event := "contact_with_self"]
            events_contacts.dt[, event_type_all := gsub(paste0("_",seq(1,100,by = 1),collapse = "|"),"",event.numb_all)]
        }
        ##add the total sum contacts for this sample
        events_contacts.dt[, sum_contacts := sum_contacts]
        if(annotate_sum) {
            names(contacts_i_j2) = c("i"," sum_i_contacts")
            events_contacts.dt = merge.data.table(events_contacts.dt, contacts_i_j2,by = "i")
            names(contacts_i_j2) = c("j","sum_j_contacts")
            events_contacts.dt = merge.data.table(events_contacts.dt, contacts_i_j2,by = "j")
        }
        if(annotate_cn) {
            nodes.dt = as.data.table(gr.chr(gg$nodes$gr))
            nodes.gr = nodes.dt %>% .[seqnames %in% paste0("chr",c(1:22,"X","Y")),] %>% GRanges(.,seqlengths = hg38_seq()) %>% gr.nochr
            events_contacts.dt = annotate_gmat_dt(events_contacts.dt, nodes.gr, "cn")
        }
            
        
########should I add contacts that were not in the original gxg object as value of 1? or 0?
                                        #Maybe also add option to return all of the data from the .hic file labeled or not labels
        return(events_contacts.dt)
    } else {
        return(NA)
    }
}


#function for getting the contacts with the event- input is the output of contacts_events and jabba for annotating the copy number of the target loci which is i
contacts_events2gr = function(contacts_event.dt,jabba_gg,seqlength, jabba_chr = TRUE) {
                                        #make seqnames i,start.i, end.i as just seqnames start and end
                                        #replace seqnames, start, end with seqnames.i.. and seqnames.j with seqnames
    names(contacts_event.dt) = names(contacts_event.dt) %>% sub("^seqnames.i$", "seqnames", .) %>% sub("^start.i$", "start", .) %>% sub("^end.i$", "end", .)
    contacts_event.gr = GRanges(contacts_event.dt, seqlengths = seqlength)
    jab = readRDS(jabba_gg)
    nodes.dt = jab$nodes$dt
    ## subset
    if(jabba_chr) {
        nodes.dt = nodes.dt[seqnames %in% paste0("chr",c(1:22,"X","Y")),]
        nodes.dt[,seqnames := gsub("chr","",seqnames)]
    } else {
        nodes.dt = nodes.dt[seqnames %in% c(1:22,"X","Y"),]
    }
    nodes.gr = GRanges(nodes.dt,seqlengths = seqlength) %>% gr.nochr
    contacts_event.gr = gr.val(contacts_event.gr,nodes.gr, "cn")
    return(contacts_event.gr)
}


#function for getting the contacts of the SV with the potential integration location
integration_site_contacts = function(contacts_event.dt, integration_gr, seqlength, jabba_gg, event_number) {
    seed.dt = contacts_event.dt
                                        #get footprint for event
    jab = readRDS(jabba_gg)
    events.dt = jab$meta$events
    events.dt[,event.numb := paste0(type,"_",1:.N), by = "type"]
    events.sub.dt = events.dt[event.numb == event_number,]
    fp = events.sub.dt$footprint
    names(seed.dt) = gsub("seqnames.i","seqnames",names(seed.dt)) %>% gsub("start.i","start",.) %>% gsub("end.i","end",.)
    seed.dt[, seqnames := paste0("chr",seqnames)]
    seed.gr = GRanges(seed.dt, seqlengths = seqlength)
                                        #adds copy number of 
    seed.gr = gr.val(gr.chr(seed.gr),jab$nodes$gr,"cn")
    fp.gr = parse.gr(events.sub.dt$footprint, seqlengths = hg38_seq())
                                        #subset i to only the node- i is other j are nodes in the event
    seed.sub.gr = seed.gr %&% integration_gr
                                        #now swap i and j to some on the events-j
    seed.sub.dt = as.data.table(seed.sub.gr)
    #replace seqnames, start, end with seqnames.i.. and seqnames.j with seqnames
    names(seed.sub.dt) = names(seed.sub.dt) %>% sub("^seqnames$", "seqnames.i", .) %>% sub("^start$", "start.i", .) %>% sub("^end$", "end.i", .) %>% sub("^seqnames.j$", "seqnames", .) %>% sub("^start.j$", "start", .) %>% sub("^end.j$", "end", .)
    seed.sub.gr2 = GRanges(seed.sub.dt, seqlengths = seqlength)
                                        #now normalize this sum to the copy number in the SV
                                        #seed.sum.gr = seed.sum.gr %Q% (!grepl("chr",seqnames ))
    nodes.dt = jab$nodes$gr %>% as.data.table()
    nodes.dt = nodes.dt[seqnames %in% paste0("chr",c(1:22,"X","Y")),]
    nodes.dt[,seqnames := gsub("chr","",seqnames)]
    nodes.gr = GRanges(nodes.dt,seqlengths = seqlength) %>% gr.nochr
    seed.sub.gr2 = gr.val(seed.sub.gr2,nodes.gr,"cn")
#add the footprint so it's easy to access
    seed.sub.gr2$footprint = paste0(gr.string(fp.gr), collapse = ",") %>% gsub("chr","",.)
                                        #return without normalizing, can do that after
    return(seed.sub.gr2)
}

##function for transforming gmat into data table with seqnames start and end for i and 
gmat_alldt = function(gmat) {
    dat = gmat$dat
    hic.gr = gmat$gr
    gr.dt = as.data.table(hic.gr)[,.(seqnames,start,end)]
    gr.dt[,index := 1:.N]
    gr.dt.i = gr.dt
    names(gr.dt.i) = c("seqnames.i","start.i","end.i","index")
    gr.dt.j = gr.dt
    names(gr.dt.j) = c("seqnames.j","start.j","end.j","index")
    dat2 = merge.data.table(dat,gr.dt.i, by.x = "i", by.y = "index")
    dat3 = merge.data.table(dat2,gr.dt.j, by.x = "j", by.y = "index")[,.(i,j,value,seqnames.i,start.i,end.i,seqnames.j,start.j,end.j)][order(i,j),]
    return(dat3)
}

##function for annotating i and j with another granges
annotate_gmat_dt = function(gmat_dt, annotate_gr, column, annotate.i = TRUE, annotate.j = TRUE) {
    if(annotate.i) {
        ## annotate i
        gmat_dt[, c("seqnames","start","end") := .(seqnames.i, start.i, end.i)]
        i.gr = GRanges(gmat_dt)
        i.gr = gr.val(i.gr, annotate_gr, column)
        colnames(mcols(i.gr))[colnames(mcols(i.gr)) == column] = paste0(column,".i") #rename column to column.i
        if(!annotate.j) {
            i.dt = as.data.table(i.gr)
            i.dt[, c("seqnames","start","end","strand","width") := NULL]
            return(i.dt)
        }
    }
    if(annotate.i && annotate.j) {
        j.dt = as.data.table(i.gr)[, c("seqnames","start","end","strand","width") := NULL]
        j.dt[, c("seqnames","start","end") := .(seqnames.j, start.j, end.j)]
        j.gr = GRanges(j.dt)
        j.gr = gr.val(j.gr, annotate_gr, column)
        colnames(mcols(j.gr))[colnames(mcols(j.gr)) == column] = paste0(column,".j") #rename column to column.j
        final.dt = as.data.table(j.gr)[, c("seqnames","start","end","strand","width") := NULL]
        return(final.dt)
    }
    if(!annotate.i && annotate.j) {
        gmat_dt[, c("seqnames","start","end") := .(seqnames.j, start.j, end.j)]
        j.gr = GRanges(gmat_dt)
        j.gr = gr.val(j.gr, annotate_gr, column)
        colnames(mcols(j.gr))[colnames(mcols(j.gr)) == column] = paste0(column,".j") #rename column to column.j
        final.dt = as.data.table(j.gr)[, c("seqnames","start","end","strand","width") := NULL]
        return(final.dt)
    }
}


## convert from contacts_events to oned tracks for pgv or gtrack
contacts2_1D_track = function(contacts.dt,
                              pairs, #with column complex and key set to get purity and ploidy
                              cores = 4,
                              ref = "hg38", #reference if return_type == "pgv"
                              return_type = "pgv", #can be data.table or pgv which will return data.tables for uploading to pgv
                              scale_factor = 1e5,
                              remove_self_amplicon_contacts = FALSE
                              ) {
    ## get the unique events present
    unique_events.dt = contacts.dt[,.(event.numb,sample)] %>% unique
    ## now sum 
    contact.1d_lst = mclapply(1:nrow(unique_events.dt), function(x) {
        event.dt = unique_events.dt[x,]
        sub.dt = contacts.dt[sample == event.dt$sample & event.numb == event.dt$event.numb,]
        sub.dt[, sum_by_i := sum(value), by = c("seqnames.i","start.i","end.i")]
        sub.dt[, c("seqnames.j","start.j","end.j","cn.j","sum_j_contacts","j","value","id") := NULL]
        ## dim(sub.dt)
        sub.dt = unique(sub.dt)[order(i),]
        return(sub.dt)
    },mc.cores = cores)
    oned.dt = rbindlist(contact.1d_lst)
    ##
    ## sum_unique.dt = oned.dt[,.(sample,sum_contacts)] %>% unique
    ## sum_unique.dt[, scaled_sum2 := scale(sum_contacts, center = 1)]
    ## sum unique events
    ## sum_events.dt = oned.dt[,.(sample,event.numb,sum_all_j_contacts)] %>% unique
    ## sum_events.dt[, sum_all_j_scaled := scale(sum_all_j_contacts, center = 1)]
    ## sum_events.dt = merge.data.table(sum_events.dt,sum_unique.dt, by = "sample")
    ## ##
    ## oned.dt = merge.data.table(oned.dt, sum_events.dt[,.(sample, event.numb, scaled_sum2, sum_all_j_scaled)], by = c("sample","event.numb"))

    ## ## add purity and ploidy
    ## oned.dt2[, event_contact_fraction := (sum_all_j_contacts / sum_contacts)]
    ##get purity and ploidy for all to then calculate relative copynumber
    if(remove_self_amplicon_contacts) {
        message("remove_self_amplicon_contacts is TRUE, removing self amplicon contacts")
        oned.dt2 = oned.dt[event_type.i != event_type.j | is.na(event_type.i),] 
    } else {
        message("remove_self_amplicon_contacts is FALSE, not removing self amplicon contacts")
        oned.dt2 = oned.dt
    }
    ## 
    samples.lst = oned.dt2$sample %>% unique
    purs.lst = mclapply(samples.lst, function(pair) {
        gg = readRDS(pairs[pair,complex])
        return(data.table(sample = pair, purity = gg$meta$purity, ploidy = gg$meta$ploidy))
    },mc.cores = cores)
    purs.dt = rbindlist(purs.lst)
    oned.dt2 = merge.data.table(oned.dt2,purs.dt, by = "sample", all.x = TRUE)
    ##
    oned.dt2[, cn.rel.i := (cn.i * purity + 2 * (1 - purity)) / (2 * (1-purity) + purity * ploidy)]
    ## now create a pgv upload temp file for each event
    sub.lst = mclapply(1:nrow(unique_events.dt), function(x) {
        event.dt = unique_events.dt[x,]
        sub.dt = oned.dt2[sample == event.dt$sample & event.numb == event.dt$event.numb,]
        sub.dt2 = sub.dt[,.(seqnames.i,start.i,end.i,event.numb,sum_by_i,sample, cn.i,sum_all_j_contacts)] %>% unique
        sub.dt2[,sum_by_i_cn_norm := (sum_by_i / (cn.i + 1))]
        sub.dt2[,sum_by_i_cn_norm_sum_norm := scale_factor * ((sum_by_i / sum_all_j_contacts) / (cn.i + 1))]
        sub.dt2[,sum_by_i_cn_norm_sum_norm_log := log(sum_by_i_cn_norm_sum_norm + 1)]
        names(sub.dt2) = c("seqnames", "start", "end", "event.numb", "sum_by_i", "sample", "cn.i", "sum_all_j_contacts","sum_by_i_cn_norm","sum_by_i_cn_norm_sum_norm","sum_by_i_cn_norm_sum_norm_log")
        return(sub.dt2)
    }, mc.cores = cores)
    sub.dt = rbindlist(sub.lst)
    if(return_type == "data.table") {
        return(sub.dt)
    } else if (return_type == "pgv") {
        arrow.add.lst = mclapply(1:nrow(unique_events.dt), function(x) {
            event.dt = unique_events.dt[x,]
            sub.dt2 = sub.dt[sample == event.dt$sample & event.numb == event.dt$event.numb,]
            sub.dt2[, seqnames := as.character(seqnames)]
            sub.gr = GRanges(sub.dt2, seqlengths = hg_seqlengths(chr = FALSE))
            arrow.add = Skilift::arrow_temp(patient_id = event.dt$sample,
                                            order = 30,
                                            x = list(sub.gr),
                                            ref = ref,
                                            field = "sum_by_i_cn_norm",
                                            title = paste0(event.dt$event.numb, " Sum by i relative CN Normalized"),
                                            overwrite = FALSE)
            return(arrow.add)
        },mc.cores = cores)
        arrow.add.dt = rbindlist(arrow.add.lst)
        return(arrow.add.dt)
    }
}

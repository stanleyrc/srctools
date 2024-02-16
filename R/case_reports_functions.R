
#create json from data.table for mutations
dt2json_mut = function(dt,patient.id,ref,settings,file_name = paste(getwd(),"test.json",sep="/"), meta_data = NULL, y_col = NULL) {
    #create vector of seqlengths
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(settings_data$coordinates$sets[[ref]])[,.(chromosome,startPoint,endPoint)]
    colnames(chrom_lengths) = c("seqnames","start","end")

    if(nrow(chrom_lengths[grepl("chr",seqnames),]) > 0) {
        chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M
    }
                                        #add y value specified
    if(is.null(y_col)) {
        dt$y_value = 1
    } else {
        dt[,y_value := get(y_col)]
        #dt[,y_value := y_value]
    }
#convert to ggraph and create json
    ## gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.chr()
    if(nrow(chrom_lengths[grepl("chr",seqnames),]) > 0) {
        gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.chr()
    } else {
        gr1 = dt2gr(dt[order(seqnames,start),]) %>% sortSeqlevels() %>% gr.nochr()
    }
    if(any(gr1@seqinfo@seqlengths > chrom_lengths[seqnames %in% names(seqlengths(gr1))]$end)) {
        stop(paste("the seqlengths of your granges has ranges that are not contained in the seqlengths of",ref))
    }
    jab = gG(nodes=gr1)
#    jab = gGnome::refresh(jab)
    settings_y = list(y_axis = list(title = "copy number",
                                    visible = TRUE))
    node.json = gr2dt(jab$nodes$gr[, c("snode.id","y_value",meta_data)])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id,title=snode.id,type="interval", y = y_value, .SD,.SDcols = meta_data)]
    names(node.json) = gsub(".SD.","",names(node.json))
    node.json[, start := NULL]
    node.json[, end := NULL]
    gg.js = list(intervals = node.json, connections = data.table())
    gg.js = c(list(settings = settings_y), gg.js)
    jsonlite::write_json(gg.js, file_name,
                         pretty=TRUE, auto_unbox=TRUE, digits=4)
    return(file_name)
}

create_somatic_json = function(somatic_snv_cn, out_file, pair, pgv_settings, return_table_pgv = FALSE, meta_keep = NULL) {
    som.dt = readRDS(somatic_snv_cn)
    som.dt = som.dt[!is.na(est_cn_llrm),]
    som.dt[start == end, end := end +1]
    som.dt[, strand := NULL]
    som.dt[variant.p != "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; Protein_variant: ", variant.p, "; VAF: ",vaf)]
    som.dt[variant.p == "",annotation := paste0("Type: ", annotation, "; Gene: ", gene, "; Variant: ",variant.c, "; VAF: ",vaf)]
    dt2json_mut(dt = som.dt, ref = "hg19",settings = pgv_settings, meta_data = meta_keep, y_col = "est_cn",
                file_name = out_file)
    if(return_table_pgv) {
        dt.add = data.table(patient.id = pair, type = "genome",visible = TRUE, title = "Copy Number Mutations", source = "mutations.json")
    }
}



#create the filtered events json for case reports
filtered_events_json = function(pair, oncotable, jabba_gg, out_file, cgc_file = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv", return_table = FALSE) {
    ##Driver CNA windows
    ##Load details from oncotable
    ot = readRDS(oncotable)
    snvs = ot[grepl('frameshift|missense|stop|disruptive', annotation)]
    snvs = snvs[!duplicated(variant.p)]
    ##Note here probably have to crossreference these missense muts with hetdels
    #hetdel_snvs = snvs[gene %in% ot[type == 'hetdel',gene]]
                                        #possible_drivers = rbind(hetdel_snvs,homdels)
    homdels = ot[type == 'homdel']
    amps = ot[type == 'amp']
    jab = readRDS(jabba_gg)
    possible_drivers = rbind(snvs,homdels,amps)
    cgc = fread(cgc_file)
    names(cgc) = gsub(' ','.', names(cgc))
    cgc$gene = cgc$Gene.Symbol
    ## longlist = merge.data.table(possible_drivers, cgc, by = 'gene', all.x = TRUE)
    longlist = merge.data.table(possible_drivers, cgc, by = 'gene')
    res = longlist[ ,.(gene, id, type, variant.p, Name, Genome.Location, Tier, Role.in.Cancer)]
    names(res) = c("gene", "id", "type", "Variant", "Name", "Genome_Location", "Tier", "Role_in_Cancer")
                                        #add copy number to homdels
    res = res %>% unique
    if(nrow(res) > 0) {
        res[,seqnames := tstrsplit(Genome_Location,":",fixed=TRUE,keep=1)]
        res[,start := tstrsplit(Genome_Location,"-",fixed=TRUE,keep=1)]
        res[,start := tstrsplit(start,":",fixed=TRUE,keep=2)]
        res[,end := tstrsplit(Genome_Location,"-",fixed=TRUE,keep=2)]
        res.mut = res[!is.na(Variant),]
        if(nrow(res.mut) > 0) {
            res.mut[,Variant := gsub("p.","",Variant)]
        }
        res.cn = res[is.na(Variant),]
        if(nrow(res.cn) >0) {
            res.cn.gr = GRanges(res.cn)
            res.cn.gr = gr.val(res.cn.gr,jab$nodes$gr,c("cn","cn.low","cn.high"))
            res.cn.dt = as.data.table(res.cn.gr)
            res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3),"; CN Low:",round(cn.low,digits = 3),"; CN High:",round(cn.high,digits = 3))]
            res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:",round(cn,digits = 3))]
            res.cn.dt[,cn := NULL]
            res.cn.dt[,cn.high := NULL]
            res.cn.dt[,cn.low := NULL]
            res.cn.dt[,width := NULL]
            res.cn.dt[,strand := NULL]
            res.final = rbind(res.mut,res.cn.dt)
        } else {
            res.final = res.mut
            res.final[,seqnames := NULL]
            res.final[,start := NULL]
            res.final[,end := NULL]
        }
        write_json(res.final, out_file, pretty=TRUE)
        res.final[,sample := pair]
        if(return_table) {
            return(res.final)
        } else {
            return(NULL)
        }
    }
}

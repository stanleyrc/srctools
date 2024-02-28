getPGV = function(web = paste0("/gpfs/commons/groups/imielinski_lab/mskiweb/",Sys.getenv("USER"),"/"), pgv_sub_folder = "pgv/", json_file = paste0(web,pgv_sub_folder, "/public/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/public/data/"), settings = paste0(web,pgv_sub_folder, "/public/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}



getsharedPGV = function(web = "/gpfs/commons/home/sclarke/lab/pgv_content/", pgv_sub_folder = "", json_file = paste0(web,pgv_sub_folder, "/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/data/"), settings = paste0(web,pgv_sub_folder, "/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}


## ## function to rebin after mask - or not mask
## cov2arrow_mask_rebin = function(pair, dryclean_cov, mask = NULL, jabba_gg, purity = NULL, ploidy = NULL) {
##     if(!is.null(jabba_gg)) {
##         cov_gr = cov2abs(dryclean_cov, jabba_gg)
##     }
##     if(!is.null(purity) && !is.null(ploidy)) {
##         cov_gr = cov2abs(dryclean_cov, purity = purity, ploidy = ploidy)
##     }
##     if(!is.null(mask)) {
##         cov_gr = gr.val(cov_gr,mask, "mask")
##         cov_gr = cov_gr %Q% (mask != "mask")
##         cov_gr$mask = NULL
##     }
##     cov_gr2 = rebin(cov_gr, 1e4, field = "foregroundabs")
##     cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0] = 0
##     add.dt = arrow_temp(patient_id = pair, ref = "hg19", field = "foregroundabs", x = list(cov_gr2), title = "mask then rebin")
##     return(add.dt)
## }


## ## function to not rebin using higlass but mask
## cov2bw = function(pair, dryclean_cov, jabba_gg, mask = NULL, seq.fix, purity = NULL, ploidy = NULL) {
##     cov_gr = cov2abs(dryclean_cov, jabba_gg)
##     cov_gr = cov2abs(dryclean_cov, jabba_gg)
##     if(!is.null(mask)) {
##         cov_gr = gr.val(cov_gr,mask, "mask")
##         cov_gr = cov_gr %Q% (mask != "mask")
##         cov_gr$mask = NULL
##     }
##     cov_gr$foregroundabs[cov_gr$foregroundabs < 0] = 0
##     ## fix seqlengths to hg19
##     cov_gr = GRanges(as.data.table(cov_gr),seqlengths = seq.fix) %>% trim()
##     add.dt = bw_temp(patient_id = pair, ref = "hg19", field = "foregroundabs", x = list(cov_gr), title = "mask no rebin")
##     return(add.dt)
## }

## #####add coverages -- let's compare binning strategies with higlass
## ## function to convert to absolute coverage ( same units as jabba graph). Using either a ggraph to get purity/plody or specified purity and ploidy
## cov2abs = function(dryclean_cov, jabba_gg = NULL, purity = NULL, ploidy = NULL, field = "foreground", new_col = "foregroundabs") {
##     cov_gr = readRDS(dryclean_cov)
##     if(!is.null(jabba_gg)) {
##         gg = readRDS(jabba_gg)
##         purity = gg$meta$purity
##         ploidy = gg$meta$ploidy
##     }
##     if(!is.null(purity) && !is.null(ploidy)) {
##         purity = purity
##         ploidy = ploidy
##     }
##     mcols(cov_gr)[new_col] = rel2abs(gr = cov_gr,
##                                      purity = purity,
##                                      ploidy = ploidy,
##                                      field = field
##                                      )
##     return(cov_gr)
## }




## ## function to not rebin using higlass, not mask and not make 0
## cov2bw = function(pair, dryclean_cov, jabba_gg, seq.fix, purity = NULL, ploidy = NULL) {
##     cov_gr = cov2abs(dryclean_cov, jabba_gg)
##     cov_gr = GRanges(as.data.table(cov_gr),seqlengths = seq.fix) %>% trim()
##     add.dt = bw_temp(patient_id = pair, ref = "hg19", field = "foregroundabs", x = list(cov_gr), title = "no mask no rebin")
##     return(add.dt)
## }

## function to rebin before mask- should not do this
## cov2arrow_rebin_mask = function(pair, dryclean_cov, mask, jabba_gg, purity = NULL, ploidy = NULL) {
##     cov_gr = cov2abs(dryclean_cov, jabba_gg)
##     cov_gr2 = rebin(cov_gr, 1e4, field = "foregroundabs")
##     cov_gr2 = gr.val(cov_gr2,mask, "mask")
##     cov_gr2 = cov_gr2 %Q% (mask != "mask")
##     cov_gr2$mask = NULL
##     cov_gr2$foregroundabs[cov_gr2$foregroundabs < 0] = 0
##     add.dt = arrow_temp(patient_id = pair, ref = "hg19", field = "foregroundabs", x = list(cov_gr2), title = "rebin then mask")
##     return(add.dt)
## }


##                                         #get templates for different plots so I don't have to remember everything that the plots need
## bw_temp = function(patient_id = NA,order = NA, x = list(NA), ref = NA, chart_type = "area", visible = TRUE, title = NA, type = "bigwig", field = "foreground", overwrite = FALSE) {
##     dt1 = data.table(patient.id = patient_id,
##                      visible = visible,
##                      x = x,
##                      type = type,
##                      field = field,
##                      ref = ref,
##                      title = title,
##                      order = order,
##                      defaultChartType = chart_type,
##                      overwrite = overwrite
##                      )
##     return(dt1)
## }

## arrow_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, chart_type = "scatterplot", visible = TRUE, title = NA, type = "scatterplot", field = "foreground", overwrite = FALSE) {
##     dt1 = data.table(patient.id = patient_id,
##                      visible = visible,
##                      x = x,
##                      type = type,
##                      field = field,
##                      ref = ref,
##                      title = title,
##                      order = order,
##                      defaultChartType = chart_type,
##                      overwrite = overwrite
##                      )
##     return(dt1)
## }

## genome_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE, title = NA, annotation = list(c('bfb','chromoplexy','chromothripsis','del','dm','dup','pyrgo','rigma','simple','tic','tyfonas')), overwrite = FALSE) {
##                                         #use type = allelic to make a color a genome graph
##     dt1 = data.table(patient.id = patient_id,
##                      type = type,
##                      visible = visible,
##                      title = title,
##                      x = x,
##                      ref = ref,
##                      order = order,
##                      annotation = annotation,
##                      overwrite = overwrite
##                      )

##     return(dt1)
## }

## walks_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "walks", visible = TRUE, title = NA, overwrite = FALSE) {
##     dt1 = data.table(patient.id = patient_id, 
##                      visible = visible,
##                      x = x,
##                      order = order,
##                      ref = ref,
##                      title = title,
##                      overwrite = overwrite
##                      )
##     return(dt1)
## }


pgv_reorder = function(pgv, update = FALSE) {
    if(!update) {
        pgv2 = PGVdb$new(datafiles_json_path = gsub("settings.json","datafiles.json", pgv$settings),datadir = pgv$datadir,settings = pgv$settings,higlass_metadata = pgv$higlass_metadata)
    } else {
        pgv2 = pgv
    }
    pgv2$plots = pgv2$plots[order(order),]
    pgv2$update_datafiles_json()
    return(pgv2)
}


## getPGV = function(web = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/", pgv_sub_folder = "pgv/", json_file = paste0(web,pgv_sub_folder, "/public/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/public/data/"), settings = paste0(web,pgv_sub_folder, "/public/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
##     ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
##     pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
##     return(pgv)
## }

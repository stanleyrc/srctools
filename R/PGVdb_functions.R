getPGV = function(web = paste0("/gpfs/commons/groups/imielinski_lab/mskiweb/",Sys.getenv("USER"),"/"), pgv_sub_folder = "pgv/", json_file = paste0(web,pgv_sub_folder, "/public/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/public/data/"), settings = paste0(web,pgv_sub_folder, "/public/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}



getsharedPGV = function(web = "/gpfs/commons/home/sclarke/lab/pgv_content/", pgv_sub_folder = "", json_file = paste0(web,pgv_sub_folder, "/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/data/"), settings = paste0(web,pgv_sub_folder, "/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
    ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}


                                        #get templates for different plots so I don't have to remember everything that the plots need
bw_temp = function(patient_id = NA,order = NA, x = list(NA), ref = NA, chart_type = "area", visible = TRUE, title = NA, type = "bigwig", field = "foreground", overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     title = title,
                     order = order,
                     defaultChartType = chart_type,
                     overwrite = overwrite
                     )
    return(dt1)
}

arrow_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, chart_type = "scatterplot", visible = TRUE, title = NA, type = "scatterplot", field = "foreground", overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     title = title,
                     order = order,
                     defaultChartType = chart_type,
                     overwrite = overwrite
                     )
    return(dt1)
}

genome_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE, title = NA, annotation = list(c('bfb','chromoplexy','chromothripsis','del','dm','dup','pyrgo','rigma','simple','tic','tyfonas')), overwrite = FALSE) {
                                        #use type = allelic to make a color a genome graph
    dt1 = data.table(patient.id = patient_id,
                     type = type,
                     visible = visible,
                     title = title,
                     x = x,
                     ref = ref,
                     order = order,
                     annotation = annotation,
                     overwrite = overwrite
                     )

    return(dt1)
}

walks_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE, title = NA, overwrite = FALSE) {
    dt1 = data.table(patient.id = patient_id, 
                     visible = visible,
                     x = x,
                     order = order,
                     ref = ref,
                     title = title,
                     overwrite = overwrite
                     )
    return(dt1)
}


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

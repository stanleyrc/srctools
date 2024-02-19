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
bw_temp = function(patient_id = NA,order = NA, x = list(NA), ref = NA, chart_type = "area", visible = TRUE, type = "bigwig", field = "foreground") {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     order = order,
                     defaultChartType = chart_type
                     )
    return(dt1)
}

arrow_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, chart_type = "scatterplot", visible = TRUE, type = "scatterplot", field = "foreground") {
    dt1 = data.table(patient.id = patient_id,
                     visible = visible,
                     x = x,
                     type = type,
                     field = field,
                     ref = ref,
                     order = order,
                     defaultChartType = chart_type
                     )
    return(dt1)
}

genome_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE) {
                                        #use type = allelic to make a color a genome graph
    dt1 = data.table(patient.id = NA,
                     type = type,
                     visible = TRUE,
                     x = NA,
                     ref=ref,
                     order = order,
                     )

    return(dt1)
}

walks_temp = function(patient_id = NA, order = NA, x = list(NA), ref = NA, type = "genome", visible = TRUE) {
    dt1 = data.table(patient.id = patient_id, 
                     visible = visible,
                     x = list(NA),
                     order = order,
                     ref = ref
                     )
    return(dt1)
}




## getPGV = function(web = "/gpfs/commons/groups/imielinski_lab/mskiweb/sclarke/", pgv_sub_folder = "pgv/", json_file = paste0(web,pgv_sub_folder, "/public/datafiles.json"), datadir = paste0(web,pgv_sub_folder, "/public/data/"), settings = paste0(web,pgv_sub_folder, "/public/settings.json"), higlass.list = list(endpoint = "https://higlass01.nygenome.org/")) {
##     ## build.dir = paste0(web,pgv_sub_folder, "/public/settings.json")
##     pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
##     return(pgv)
## }

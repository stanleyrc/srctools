
## example
## gm2hicool(gm = gm, hic_out = "~/Projects/temp_files/test.hic", juicer_jar = "~/modules/Juicer/juicer/scripts/juicer_tools_linux_0.8.jar", mcool = TRUE, res = 50000, hic_res = c(50000,100000,500000))

## wrapper function to convert from gmatrix to .hic (and cooler)
gm2hicool = function(gm, #gmatrix to convert to .hic
                     hic_out, #hic file out
                     mcool = TRUE, #mcool file out, will not be made if null, if true, will just 
                     temp_file = NULL, # by default will just make a temp file with the name of hic_out plus _temp.tsv
                     keep_temp = FALSE, #if TRUE will keep the temporary valid pairs that is inputted to juicer
                     chrom_sizes = NULL, # if null will create one in the same folder from the seqlengths of the gr in the gm
                     res = NULL,        #if NULL will try to get the smallest res in the .hic file
                     hic_res = NULL,           #list of resolutions for the .hic file, mcool will be the same if not specified, will be converted to an integer
                     mem = 8, #memory for java for validPairs2hic
                     juicer_jar = "~/modules/Juicer/juicer/scripts/juicer_tools_linux_0.8.jar",        #path to jar for juicer
                     python_path = "/gpfs/commons/home/sclarke/miniconda3/envs/hic2cool/bin/python3.12" #python path with hic2cool installed
                     ) {
    ## get res for gm2ValidPairs and validPairs2hic
    if(is.null(res)) {
        ## find res
        gr.dt = as.data.table(gm$gr)
        res = gr.dt$width %>% unique
        if(length(res) != 1) {
            stop("res was not provided and there is not a consistent width of the granges in the gM")
        } else {
            message(paste0("res of your gM was found to be: ",res))
        }
    }
    out.dt = gm2ValidPairs(gm, res = res)
    if(is.null(temp_file)) {
        temp_file = paste0(hic_out,"_temp.tsv")
    }
    message(paste0("Writing temporary valid pairs for juicer to: ",temp_file))
    fwrite(out.dt, temp_file, col.names = FALSE, sep = "\t")
    ## make chromosome sizes tsv if null with validPairs2hic
    if(is.null(chrom_sizes)) {
        chrom_sizes = seqlengths(gm$gr) %>% si2gr %>% as.data.table
        chrom_sizes = chrom_sizes[,.(seqnames,end)]
        chrom_sizes[, seqnames := gsub("chr","",seqnames)]
        chrom_file = paste0(hic_out,"_chromosome.sizes_temp.tsv")
        message(paste0("Writing temporary chromosome lengths to: ", chrom_file))
        fwrite(chrom_sizes, chrom_file, col.names = FALSE, sep = "\t")
    } else { #assume the chromosome sizes is a tsv of chromosome sizes
        chrom_file = chrom_sizes
    }
    validPairs2hic(validPairs_path = temp_file, hic_out = hic_out, hic_res = hic_res, chrom_sizes = chrom_file, juicer_jar = juicer_jar, mem = mem)
    if(!file.exists(hic_out)) {
        stop("hic_out does not exist")
    }
    if(is.null(mcool)) {
        message("mcool is NULL, will not convert to mcool")
    } else {
        if(isTRUE(mcool)) {
            ##name mcool if just TRUE
            mcool = unlist(strsplit(x = hic_out,".", fixed = TRUE))
            mcool = paste0(mcool[1:(length(mcool)-1)], collapse = ".")
            mcool = paste0(mcool,".mcool")
        }        ## if mcool is just a string, it will write there
        hic2mcool(hic = hic_out, mcool_out = mcool, python_path = python_path)
    }
    if(!keep_temp) {
        ## will remove temp_file (temporary valid pairs file for hic) and chrom_file
        chrom_file = normalizePath(chrom_file)
        temp_file = normalizePath(temp_file)
        message(paste0("Removing temporary files : ", chrom_file, " and ", temp_file))
        system(paste0("rm ",chrom_file), intern = FALSE)
        system(paste0("rm ",temp_file), intern = FALSE)
    }
    if(file.exists(mcool) & file.exists(hic_out)) {
        message(paste0(".mcool and .hic are located: ", mcool, " and ", hic_out))
    } else if (file.exists(hic_out)) {
        message(paste0(".hic file is located: ", hic_out))
    }
}


## returns a validPairs data.table that can be saved after
gm2ValidPairs = function(gm, #gmatrix object
                         res = NULL #resolution of the gmatrix object
                         ) {
    dat1 = gm$dat
    gr.dt = as.data.table(gm$gr)
    if(is.null(res)) {
        ## find res
        res = gr.dt$width %>% unique
        if(length(res) != 1) {
            stop("res was not provided and there is not a consistent width of the granges in the gM")
        } else {
            message(paste0("res of your gM was found to be: ",res))
        }
    }
    ## always remove chr here
    gr.dt[, seqnames := gsub("chr","",seqnames)]
    gr.dt[, id := 1:.N]
    gr.dt[, c("width","strand") := NULL]
    gr.dt = gr.dt[,.(seqnames,start,end,id)]
    names(gr.dt) = c("seqnames.i","start.i", "end.i","id")
    gr.dt2 = copy(gr.dt)
    names(gr.dt2) = gsub(".i",".j",names(gr.dt2))
    ## merge coordinates onto data.table
    dt1 = merge.data.table(dat1, gr.dt, by.x = "i", by.y = "id")
    dt2 = merge.data.table(dt1, gr.dt2, by.x = "j", by.y = "id")
    ## add dummy values- necessary for .hic
    dt2[, str1 := 1]
    dt2[, str2 := 1]
    dt2[, frag1 := 0]
    dt2[, frag2 := 1]
    dt3 = dt2[,.(str1, seqnames.i, start.i, frag1, str2, seqnames.j, start.j, frag2, value)]
    ## dt3[, c("end.i","end.j") := NULL]
    dt3[, seqnames.i := factor(seqnames.i,levels = c(1:22,"X","Y"))]
    dt3[, seqnames.j := factor(seqnames.j,levels = c(1:22,"X","Y"))]
    dt3 = dt3[order(seqnames.i,seqnames.j,start.i,start.j),]
    return(dt3)
}

## converts a written validPairs that is a tsv to .hic
validPairs2hic = function(validPairs_path, #path to valid pairs tsv
                          hic_out,         #path to create the .hic file at
                          res,             #resolution of the .hic file-necessary
                          hic_res = NULL,  #resolutions to output
                          chrom_sizes = NULL, #chromosome sizes tsv, format: chromosome length
                          juicer_jar,         #jar for juicer
                          mem = 8 #memory for java
                          ) {
    ## make seqlengths file if NULL
    if(is.null(chrom_sizes)) {
        chrom_sizes = chrom_sizes[,.(seqnames,end)]
        chrom_file = paste0(hic_out,"_chromosome.sizes_temp.tsv")
        message(paste0("Writing temporary chromosome lengths to: ", chrom_file))
        fwrite(chrom_sizes, chrom_file, col.names = FALSE, sep = "\t")
    } else {
        chrom_file = chrom_sizes
    }
    if(is.null(hic_res)) {
        hic_res = res * (10^(0:6))
        hic_res = hic_res[hic_res < 1e8]
        hic_res = as.integer(hic_res)
    } else {
        hic_res = as.integer(hic_res)
    }
    if(length(hic_res) > 1) {
        hic_res = paste0(hic_res, collapse = ",")
    }
    message("Converting ", validPairs_path, " to .hic at ", hic_out)
    ## now convert to .hic
    ## cmd = paste0("java -Xmx2g -jar ",juicer_jar, " pre -r ",hic_res, " ",validPairs_path," ", hic_out, " ",chrom_file," -n")
    cmd = paste0("java -Xmx",mem,"g -jar ",juicer_jar, " pre -r ",hic_res, " ",validPairs_path," ", hic_out, " ",chrom_file," -n")
    message(paste0("Running ",cmd))
    system(command=cmd,intern=TRUE)
    message("Done converting .hic")
}


hic2mcool = function(hic, #input .hic file
                     mcool_out, #.mcool path to write to
                     python_path = "/gpfs/commons/home/sclarke/miniconda3/envs/hic2cool/bin/python3.12" #python path with hic2cool
                     ) {
    hic = normalizePath(hic, mustWork = TRUE) #has to exist
    mcool_out = normalizePath(mcool_out,mustWork = FALSE) #does not have to exist
    resolution = 0 #will convert all resolutions present in the .hic to multi resolution mcool file
    python_command = sprintf(
        "from hic2cool import hic2cool_convert; hic2cool_convert('%s', '%s', resolution=%d)",
        hic, mcool_out, resolution
    )
    ## command for converting from .hic to mcool
    cmd2 = sprintf('%s -c "%s"', python_path, python_command)
    ## Print the command to check it
    message("Running ", cmd2, " to convert")
    ## Run the command using system()
    output = system(cmd2, intern = TRUE)
}

## cannot upload coolers through the api to higlass at the moment but this works!
upload_cooler = function(path, name, uuid = NULL) {
    if(is.null(uuid)) {
        uuid = paste0(name, "_",sample(1:1e7,size = 1))
    }
    cmd <- paste0('ssh -t higlass01.nygenome.org  "module purge ; module load higlass-server/1.7.2 ; cd /opt/software/higlass-server/higlass-server-1.7.2/higlass-server ; \n /opt/software/higlass-server/higlass-server-1.7.2/bin/python /opt/software/higlass-server/higlass-server-1.7.2/higlass-server/manage.py ingest_tileset --filename ', path, ' --filetype cooler --datatype matrix --name ', name, ' --uid ',uuid,'"')
    system(cmd,intern=T)
}
## upload_cooler(mcool, name = "june_upload_cooler_test")






bind_years <- function() {}





clean_flows <- function(raw, clean) {
    files <- list.files(raw)
    ## omit years before 2000
    for (dir in files[7:length(files)]) {
        message(dir)
        fs <- list.files(file.path(raw, dir), pattern = ".xls", full.names = TRUE)
        kws <- lapply(fs, function(f) read_and_clean(f))
        kws <- rbindlist(kws)
        year <- gsub("kws_", "", dir)
        kws[, "year" := as.integer(year)]
        fwrite(kws, file.path(clean, paste0(dir, ".csv")), quote = TRUE)
        message(sprintf("file %s written to disk", dir))
    }
    return(NULL)
}

read_and_clean <- function(file) {
    ## thanks chatGPT!
    file_short <- rev(strsplit(rev(file), "/", fixed = TRUE)[[1]])[1]
    message(sprintf("File %s is being processed", file_short))
    sheets <- readxl::excel_sheets(file)
    sheets <- sheets[!grepl("Außen", sheets)]
    message(sheets)
    kws <- list()
    for (s in sheets) {
        message(s)
        kws[[s]] <- readxl::read_excel(file, col_names = FALSE, sheet = s)
    }
    #    kws <- lapply(sheets, function(s) (readxl::read_excel(file, col_names = FALSE, sheet = s)))
    kws <- lapply(kws, function(x) clean_districts(data.table::setDT(x)))
    kws <- data.table::rbindlist(kws)
    return(kws)
}

clean_districts <- function(kws, cols) {
    kws <- drop_cols_only_na(kws)
    kws <- renew_colnames(kws)
    keep <- c("ags_o", "ags_d", "age_group", "wins_g_a", "losses_g_a")
    kws <- na.omit(kws[-c(1:4), ..keep])
    old_cols <- c("ags_o", "ags_d", "wins_g_a", "losses_g_a")
    new_cols <- c("origin", "destination", "wins", "losses")
    setnames(kws, old_cols, new_cols)
    kws[, "origin" := as.integer(origin)]
    kws[, "destination" := as.integer(destination)]
    kws[, "wins" := as.integer(wins)]
    kws[, "losses" := as.integer(losses)]
    ## Add checks if wins and losses by age_group are the same
    kws <-  kws[, .(origin, destination, age_group, flow = wins)]
    kws <- kws[order(origin, destination)]
    return(kws)
}

drop_cols_only_na <- function(kws) {
    cols <- colnames(kws)
    idx <- kws[, lapply(.SD, function(x) nrow(kws) == sum(is.na(x))), 
               .SDcols = cols]
    idx <- as.logical(idx)
    if (length(cols[idx]) == 0) {
        message("No columns dropped")
    }
    if (length(cols[idx]) > 0) {
        set(x = kws, j = cols[idx], value = NULL)
        message(sprintf("dropped %s only na columns.", length(cols) - ncol(kws)))
    }
    return(kws)
}

renew_colnames <- function(kws) {
    n_cols_kws <- ncol(kws)
    cols <- c("bl_o", "ags_o", "h-ld", "h-rb", "h-kr", "name_o", "bl_d",
              "ags_d", "z-ld", "z-rb", "z-kr", "name_d", "reg_type",
              "age_group", "wins_a_a", "wins_a_m", "wins_a_f", "wins_g_a",
              "wins_g_m", "wins_g_f", "wins_i_a", "wins_i_m", "wins_i_f",
              "losses_a_a", "losses_a_m", "losses_a_f", "losses_g_a",
              "losses_g_m", "losses_g_f", "losses_i_a", "losses_i_m",
              "losses_i_f")
    n_cols_new <- length(cols)
    if (n_cols_new == n_cols_kws) {
        colnames(kws) <- cols
        message("set column names")
    }
    if (n_cols_new - 1 == n_cols_kws) {
        ## then there is column "reg_type" missing
        idx <- grep("reg_type", cols)
        colnames(kws) <- cols[- idx]
        message("set column names without reg_type")
    }
    if (n_cols_new != n_cols_kws &
        n_cols_new - 1 != n_cols_kws) {
        stop("Unknown number of columns in data!")
    }
    return(kws)
}


## The following is out commented because I did this once and now, if
## I want to start a new project, I just copy the unzipped files with
## the unified file names

### This I only do once if the files are newly unzipped
## the unzipping I did by hand. Not sure if I can enter a passphrase
## in R
## p <- file.path(paths$raw_flows_di, "unzipped/flows_districts")
## files <- list.files(p, pattern = "KWS")
## ### I move files to the unzipped folder
## lapply(files, function(f) file.rename(from = file.path(p, f),
##                                       to = file.path(p_uz, f)))
## ## no I unify file names. Probably more clear in loop. I think I can
## ## do both steps in one step
## lapply(files, function(f) file.rename(from = file.path(p_uz, f),
##                                       to = file.path(p_uz,
##                                                      gsub("[^0-9]+", "kws_", f))))



### clean_districts only works from 2000 onwards, before different format

## for (dir in files[7:length(files)]) {
##     message(dir)
##     fs <- list.files(file.path(p_uz, dir), pattern = ".xls", full.names = TRUE)
##     print(lapply(fs, excel_sheets))
## }

## fs <- list.files(file.path(p_uz, files[9]), pattern = ".xls", full.names = TRUE)
## sheets <- lapply(fs, excel_sheets)


### es sieht aus, als fehlt bei einigen Dateien die Spalte "regionstyp"
## test <- setDT(read_excel(file.path(p_uz, "kws_2000", 'Saarland 2000.xls'), col_names = FALSE))


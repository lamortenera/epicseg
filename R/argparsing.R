parseArgs <- function(optList, args, launcherName="<script>"){
    opts <- makeOptions(optList)
    #parse the arguments
    tryCatch(values <- matchOptsWithArgs(opts, args),
        error=function(e){
            cat("@@@ ERROR parsing command line arguments @@@", fill=TRUE)
            cat(e$message, fill=TRUE)
            cat("@@@", fill=TRUE)
            printUsage(launcherName, opts)
            stop()
        })
    
    values
}


argName <- function(arg){
    if (is.null(arg) || !grepl("^--", arg) || nchar(arg) < 3) return(NA)
    sub("^--", "", arg)
}

checkOption <- function(opt){
    name <- argName(opt$arg)
    if (is.na(name)) stop("invalid 'arg' argument (no 'arg' argument? forgot the '--'?)")
    if (is.null(opt$type)) opt$type <- "character"
    if (is.null(opt$flag)) opt$flag <- FALSE
    if (is.null(opt$help)) opt$help <- "undocumented option"
    if (is.null(opt$required)) opt$required <- FALSE
    if (is.null(opt$vectorial)) opt$vectorial <- FALSE
    if (is.null(opt$parser)) opt$parser <- identity
    if (is.null(opt$meta)) opt$meta <- name
    if (is.null(opt$short)) opt$short <- paste0("-", substr(name, 1,1))
    opt$meta <- toupper(opt$meta)
    if (opt$vectorial && opt$flag) stop("a flag cannot be vectorial")
    #check if the default makes sense
    if (!is.null(opt$default)){
        tryCatch({
        opt$default <- convertArg(opt$default, opt$type)
        opt$parser(opt$default)}
        ,error = function(e){ stop(paste0("the default argument for option ", 
            name,  " causes an error"))})
        #if there is a default the option cannot be required
        opt$required <- FALSE
    }
    
    opt
}

makeOptions <- function(optList){
    optNames <- lapply(optList, function(opt) argName(opt$arg))
    if (any(is.na(optNames))) stop("provide a 'arg' argument starting with '--' for each option")
    if (anyDuplicated(optNames)) stop("dupicated option names are invalid")
    opts <- lapply(optList, checkOption)
    names(opts) <- optNames
    
    #make sure that the short flags are unique
    shorts <- lapply(opts, function(opt) opt$short)
    shorts[duplicated(shorts)] <- NA
    for (i in seq_along(opts)) opts[[i]]$short <- shorts[[i]]
    
    opts
}

convertArg <- function(arg, type){
    res <- as(arg, type)
    if (is.na(res)) stop(paste0("failed to convert argument ", arg, " to type ", type))
    res
}

matchOptsWithArgs <- function(opts, args){
    values <- list()
    i <- 1
    len <- length(args)
    
    #make aliases based on short flags
    aliases <- list()
    for (name in names(opts)){
        shname <- opts[[name]]$short
        if (!is.na(shname)) aliases[[shname]] <- name
    }
    
    allowedNames <- names(opts)
    allowedAliases <- names(aliases)
    
    while (i <= len){
        arg <- args[i]
        #get option name taking care of the aliases
        name <- argName(arg)
        if (is.na(name)){
            if (!arg %in% allowedAliases) {
                stop(paste0("'", arg, "' is neither an option name nor an alias"))
            } else {
                name <- aliases[[arg]]
            }
        } 
        if (!name %in% allowedNames) stop(paste0("invalid option name: ", arg))
        
        #parse the argument value
        if (opts[[name]]$flag){
            #the option is a flag
            if (!is.null(values[[name]])) stop("a flag can be specified only once")
            values[[name]] <- TRUE
        } else {
            #the option is not a flag
            i <- i+1
            if (i > len) stop("missing value for the last option")
            value <- convertArg(args[i], opts[[name]]$type)
            values[[name]] <- c(values[[name]], value)
        }
        i <- i+1
    }
    
    #add default arguments
    for (name in setdiff(allowedNames, names(values))){
        def <- opts[[name]]$default
        if (!is.null(def)){
            values[[name]] <- def
        }
    }
    
    #check required args
    requiredNames <- allowedNames[sapply(opts, function(opt) opt$required)]
    missingArgs <- setdiff(requiredNames, names(values))
    if (length(missingArgs)>0) stop(paste("Missing required arguments:", paste(collapse=", ", missingArgs)))
    
    #apply functions to each argument
    for (name in names(values)){
        value <- values[[name]]
        opt <- opts[[name]]
        #check if the vectorial constraint is satisfied
        if (!opt$vectorial && length(value)>1) stop(paste0("the field '", name, "' can be specified only once"))
        values[[name]] <- opt$parser(value)
    }
    
    values
}

toString_safe <- function(obj, maxChars=20){
    s <- tryCatch(
        toString(obj), error=function(e) paste0(class(obj), " object")
    )
    #remove multiple lines, if any
    s <- gsub("\n.*$", "", s)
    #bound the number of characters
    if (nchar(s) > maxChars) s <- paste0(substr(s, 1, maxChars), "...")
    s
}

printUsage <- function(launcherName, opts, progDesc=NULL){
    if (!is.null(progDesc)) cat(progDesc, fill=TRUE)
    #divide arguments into required and not required
    reqArgs <- list()
    optArgs <- list()
    for (name in names(opts)){
        if (opts[[name]]$required){
            reqArgs[[name]] <- opts[[name]]
        } else {
            optArgs[[name]] <- opts[[name]]
        }
    }
    
    cat(sub("%prog", launcherName, "usage: %prog [arguments]"), fill=TRUE)
    cat("\n")
    
    printArgument <- function(name, opt){
        cat(sep="", "  --", name)
        if (!opt$flag) cat(sep="", " ", opt$meta)
        if (!is.na(opt$short)) {
            cat(sep="", ", ", opt$short)
            if (!opt$flag) cat(sep="", " ", opt$meta)
        }
        if (opt$vectorial) cat(sep="", " (repeatable)")
        if (!is.null(opt$default)) {
            cat(sep="", " (default: ", toString_safe(opt$default), ")")
        }
        cat("\n    ")
        cat(opt$help, fill=TRUE)
        cat("\n")
    }
    
    #print required argument first
    if (length(reqArgs) > 0){
        cat("REQUIRED ARGUMENTS:\n")
        for (name in names(reqArgs)) printArgument(name, reqArgs[[name]])
    }
    #print optional arguments
    if (length(optArgs) > 0){
        cat("OPTIONAL ARGUMENTS:\n")
        for (name in names(optArgs)) printArgument(name, optArgs[[name]])
    }
    
}

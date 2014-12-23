parseArgs <- function(optList, args, launcherName="<script>"){
	opts <- makeOptions(optList)
	#parse the arguments
	tryCatch(values <- matchOptsWithArgs(opts, args),
		error=function(e){
			cat("@@@ ERROR parsing command line arguments @@@", fill=TRUE)
			cat(e$message, fill=TRUE)
			cat("@@@", fill=TRUE)
			printUsage(launcherName, opts)
			quit(status=1)
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

matchOptsWithArgs <- function(opts, args){
	values <- list()
	i <- 1
	len <- length(args)
	
	#make aliases based on short flags
	aliases <- list()
	for (name in names(opts)){
		opt <- opts[[name]]
		if (!is.na(opt$short)) aliases[[opt$short]] <- name
	}
	
	allowedNames <- names(opts)
	allowedAliases <- names(aliases)
	
	while (i < len){
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
		
		if (opts[[name]]$flag){
			#the option is a flag
			if (!is.null(values[[name]])) stop("a flag can be specified only once")
			values[[name]] <- TRUE
		} else {
			#the option is not a flag
			i <- i+1
			if (i > len) stop("missing value for the last option")
			value <- as(args[i], opts[[name]]$type)
			values[[name]] <- c(values[[name]], value)
		}
		i <- i+1
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

printUsage <- function(launcherName, opts, progDesc=NULL){
	if (!is.null(progDesc)) cat(progDesc, fill=TRUE)
	cat(sub("%prog", launcherName, "usage: %prog [options]"), fill=TRUE)
	cat("\n")
	cat("Options:\n")
	for (name in names(opts)){
		opt <- opts[[name]]
		cat(sep="", "  --", name)
		if (!opt$flag) cat(sep="", " ", opt$meta)
		if (!is.na(opt$short)) {
			cat(sep="", ", ", opt$short)
			if (!opt$flag) cat(sep="", " ", opt$meta)
		}
		if (opt$required) cat(sep="", " (required)")
		if (opt$vectorial) cat(sep="", " (repeatable)")
		cat("\n\t")
		cat(opt$help, fill=TRUE)
		cat("\n")
	}
}

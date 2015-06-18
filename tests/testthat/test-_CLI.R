context("CLI")

#this sets up the CLI interface and the shortcuts ${indir} and ${outdir}
source("utils.R")
configureSys()

#this is to fix some obscure behaviour of R CMD check
Sys.setenv("R_TESTS" = "")

#run a command without displaying absolutely anything
runQuiet <- function(cmd){
    code <- paste('-c', shQuote(cmd))
    #errors can happen at the R level, or the command can cause the error
    #the first case is handled by this try catch,
    #the second case is handled by checking if out has attribute 'status'
    #in both cases out should contain the error message and a status flag
    out <- suppressWarnings(
        tryCatch(
        system2("bash", code, stdout=FALSE, stderr=TRUE),
        error=function(e){
            msg <- e$msg
            attr(msg, "status") <- 1
            msg
    }))
    
    if (!is.null(attr(out, "status")))  stop(paste(out, collapse = '\n'))
}

#test that a command runs
truns <- function(cmd) expect_runs(runQuiet(cmd), label=cmd)
#test that a command fails
tfails <- function(cmd) expect_error(runQuiet(cmd), label=cmd)

test_that("Command line interface",{
    #GETCOUNTS
    #this creates $outdir/countmat.txt, 
    #which we will base other tests on
    truns(
    "epicseg.R getcounts --target ${outdir}/countmat.txt \\
    --regions ${indir}/contigs.bed \\
    --mark H3K27ac:${indir}/H3K4me3.bam \\
    --mark H3K27me3:${indir}/H3K36me3.bam \\
    --mark H3K36me3:${indir}/H3K9me3.bam \\
    --mark H3K4me1:${indir}/H3K4me3.bam \\
    --mark H3K4me3:${indir}/H3K36me3.bam \\
    --mark H3K9me3:${indir}/H3K9me3.bam \\
    --mark Input:${indir}/H3K4me3.bam")
    
    
    #check if the automatic region trimming works
    truns(
    "epicseg.R getcounts -r ${indir}/contigs.bed \\
    --mark H3K4me3:${indir}/H3K4me3.bam --mark H3K36me3:${indir}/H3K36me3.bam \\
     -b 157 -t ${outdir}/counts_prova.rda")
    
    #let's try to give three bam files and two pairedends options (it should throw an error)
    tfails(
    "epicseg.R getcounts -r ${indir}/contigs.bed \\
    -m H3K4me3:${indir}/H3K4me3.bam \\
    -m H3K36me3:${indir}/H3K36me3.bam \\
    -m H3K27me3:${indir}/H3K9me3.bam \\
    -b 157 -t ${outdir}/counts_prova.txt -p T -p F")
    
    #let's try to give the -p option without a value (it should throw an error)
    tfails(
    "epicseg.R getcounts -r ${indir}/contigs.bed \\
    -m H3K4me3:${indir}/H3K4me3.bam \\
    -b 157 -t ${outdir}/counts_prova.txt -p")
    
    #check if saving to txt works
    truns(
    "epicseg.R getcounts -r ${indir}/contigs.bed \\
    --mark H3K4me3:${indir}/H3K4me3.bam \\
    --mark H3K36me3:${indir}/H3K36me3.bam \\
    -b 157 -t ${outdir}/counts_prova.txt")
    
    #check if it runs with replicate experiments
    truns(
    "epicseg.R getcounts -r ${indir}/contigs.bed \\
    --mark H3K4me3:${indir}/H3K4me3.bam \\
    --mark H3K36me3:${indir}/H3K36me3.bam \\
    --mark H3K36me3:${indir}/H3K9me3.bam \\
     -b 157 -t ${outdir}/counts_prova.txt")
    
    
    #SEGMENT
    #check if reading counts from a text file works
    #we will use the output file ${outdir}/model.txt for later tests
    truns(
    "epicseg.R segment -c ${outdir}/countmat.txt -r \\
    ${indir}/contigs.bed -n 10 --nthreads 4 \\
    -a genes:${indir}/genes.bed --maxiter 20  --outdir ${outdir}")
    
    #segment with --model option
    truns(
    "epicseg.R segment -n 10 -m ${outdir}/model.txt \\
    -c ${outdir}/countmat.txt -r ${indir}/contigs.bed --outdir ${outdir}/")
    
    #let's try with an incomplete model (we delete initP and transP from model.txt)
    system("head ${outdir}/model.txt -n 15 > ${outdir}/incomplete_model.txt")
    
    truns(
    "epicseg.R segment -n 10 -m ${outdir}/incomplete_model.txt \\
    -c ${outdir}/countmat.txt -r ${indir}/contigs.bed --outdir ${outdir}")
    
    #let's try out the predict, collapseInitP and save_rdata flags
    truns(
    "epicseg.R segment --collapseInitP T --save_rdata \\
    --notrain -m ${outdir}/model.txt -c ${outdir}/countmat.txt \\
    -r ${indir}/contigs.bed -n 10 --outdir ${outdir}")
    
    #check whether the rdata was created
    truns("ls ${outdir}/rdata.Rdata")
    
    #REPORT
    truns(
    "epicseg.R report -m ${outdir}/model.txt -s ${outdir}/segmentation.bed \\
    --outdir ${outdir}")
    
    #let's try the colors option
    truns(
    "epicseg.R report --colors ${outdir}/colors.txt \\
    -m ${outdir}/model.txt -s ${outdir}/segmentation.bed --outdir ${outdir}")
    #let's try the labels option
    
    writeLines("ale\nmarco\nandrea\ngiacomo\nanna\njuliane\nbarbara\nisabella\nfrancesco\nchiara", 
    file.path(Sys.getenv("outdir"), "labels.txt"))
    truns(
    "epicseg.R report --labels ${outdir}/labels.txt \\
    -m ${outdir}/model.txt -s ${outdir}/segmentation.bed --outdir ${outdir}")
    truns("ls ${outdir}/segmentation_labelled.bed")
        
    #let's try multiple annotations
    truns(
    "epicseg.R report \\
    -a genes:${indir}/genes.bed --annot genes2:${indir}/genes.bed \\
    -m ${outdir}/model.txt -s ${outdir}/segmentation.bed --outdir ${outdir}")
    
    truns("ls ${outdir}/annot_genes.txt")
    truns("ls ${outdir}/annot_genes2.txt")
    
    #now try multiple segmentations
    truns(
    "epicseg.R report --labels ${outdir}/labels.txt  -o ${outdir} -m ${outdir}/model.txt \\
    -s seg1:${outdir}/segmentation.bed -s seg2:${outdir}/segmentation.bed \\
    -a genes:${indir}/genes.bed")
    
})


test_that("multiple datasets",{
    #make 3 count matrices
    os <- list(c(1,2,3), c(2,1,3), c(3,2,1))
    os <- lapply(os, function(o) c("H3K4me3.bam", "H3K36me3.bam", "H3K9me3.bam")[o])
    dsets <- paste0("cmat", 1:3)
    targets <- file.path("${outdir}", paste0(dsets, ".txt"))
    for (i in seq_along(os)){
        o <- os[[i]]
        mline <- paste0(collapse=" ", "-m ", 
        c("H3K4me3", "H3K36me", "H3K9me3"), ":${indir}/", o)
        truns(paste("epicseg.R getcounts -r ${indir}/contigs.bed -t ", targets[i], mline))
    }
    
    #run normalize counts
    cline <- paste0(collapse=" ", "-c ", targets)
    truns(paste0("epicseg.R normalizecounts ", cline))
    
    #check that the new matrices have been created according to the suffix
    newtargets <- gsub(".txt$", paste0(defSuffix, ".txt"), targets)
    for (t in newtargets) truns(paste0("ls ", t))
    
    #check that using the 'triggerOverwrite' suffix no new files are created 
    #(in the temporary directory)
    outdir <- Sys.getenv("outdir")
    lfOld <- list.files(outdir)
    
    cline <- paste0(collapse=" ", "-c ", newtargets)
    truns(paste0("epicseg.R normalizecounts ", cline, " -s ", triggerOverwrite))
    expect_true(setequal(list.files(outdir), lfOld))
    
    #check that the segmentation runs with multiple datasets
    goodcline <- paste0(collapse=" ", "-c ", dsets, ":", newtargets)
    truns(paste0("epicseg.R segment ", goodcline, " -n 5 -r ${indir}/contigs.bed -o ${outdir}"))
    
    #check that without labels it fails!
    badcline <- paste0(collapse=" ", "-c ", newtargets)
    tfails(paste0("epicseg.R segment ", badcline, " -n 5 -r ${indir}/contigs.bed -o ${outdir}"))
})

test_that("usage examples", {
    #get location of the Rmarkdown file
    rmd <- system.file("extdata", "cliExamples.Rmd", package="epicseg")
    library(knitr)
    #write markdown output in a temporary file
    tmp <- tempfile("cliexamples", tmpdir=Sys.getenv("outdir"), fileext=".md")
    expect_runs(knit(rmd, tmp))
})

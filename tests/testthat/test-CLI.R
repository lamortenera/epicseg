context("CLI")
source("utils.R")

cmd_sc <- function(cmd){
    l <- list(cmd)
    for (i in names(shortcuts)){
        l[[1]] <- gsub(i, shortcuts[[i]], l[[1]])
    }
    l[[1]]
}

simuCLI_shortcuts <- function(cmd){
    simulateCLI(cmd_sc(cmd))
}
#tries to remove any kind of output
no <- function(expr){
    suppressMessages(capture.output(expr))
}
#test that a command with shortcuts runs
truns <- function(cmd) expect_runs(no(simuCLI_shortcuts(cmd)), label=cmd)
#test that a command with shortcuts fails
tfails <- function(cmd) expect_error(no(simuCLI_shortcuts(cmd)), label=cmd)

test_that("Command line interface",{
    #GETCOUNTS
    truns(
    "init getcounts --target tmpdir/countmat.txt 
    --regions regs.bed \
    --mark H3K27ac:bam1 \
    --mark H3K27me3:bam2 \
    --mark H3K36me3:bam3 \
    --mark H3K4me1:bam1 \
    --mark H3K4me3:bam2 \
    --mark H3K9me3:bam3 \
    --mark Input:bam1")

    #check if the automatic region trimming works
    truns(
    "init getcounts -r regs.bed \
    --mark H3K4me3:bam1 --mark H3K36me3:bam2  -b 157 -t tmpdir/counts_prova.rda")
    
    #let's try to give three bam files and two pairedends options (it should throw an error)
    tfails(
    "init getcounts -r tests/data/regions.bed \
    -m H3K4me3:bam1 \
    -m H3K36me3:bam2 \
    -m H3K27me3:bam3 \
    -b 157 -t tmpdir/counts_prova.txt -p T -p F")
    
    #let's try to give the -p option without a value (it should throw an error)
    tfails(
    "init getcounts -r tests/data/regions.bed \
    -m H3K4me3:bam1 \
    -b 157 -t tmpdir/counts_prova.txt -p")
    
    #check if saving to txt works
    truns(
    "init getcounts -r tests/data/regions.bed \
    --mark H3K4me3:bam1 \
    --mark H3K36me3:bam2 \
    -b 157 -t tmpdir/counts_prova.txt")
    
    #check if it runs with replicate experiments
    truns(
    "init getcounts -r tests/data/regions.bed \
    --mark H3K4me3:bam1 \
    --mark H3K36me3:bam2 \
    --mark H3K36me3:bam3 \
     -b 157 -t tmpdir/counts_prova.txt")
    
    #SEGMENT
    #check if reading counts from a text file works
    truns(
    "init segment -c tmpdir/counts_prova.txt -r \
    tmpdir/counts_prova_refined_regions.bed -n 10 --nthreads 4 
    -a genes:bed1 --maxiter 20  --outdir tmpdir")
    
    truns(
    "init segment -c tests/data/randomCounts.txt -r \
    tests/data/regions.bed -n 10 --nthreads 4 \
    -a genes:bed1  --maxiter 20 --outdir tmpdir")
    
    #segment with --model option (the likelihood should start from a lower value than before)
    truns(
    "init segment -n 10 -m tmpdir/model.txt \
    -c tests/data/randomCounts.txt -r tests/data/regions.bed --outdir tmpdir/")
    
    #let's try with an incomplete model (we delete initP and transP from model.txt)
    system(gsub("tmpdir", tmpdir, "head tmpdir/model.txt -n 15 > tmpdir/incomplete_model.txt"))
    truns(
    "init segment -n 10 -m tmpdir/incomplete_model.txt \
    -c tests/data/randomCounts.txt -r tests/data/regions.bed --outdir tmpdir")
    
    #let's try out the predict, collapseInitP and save_rdata flags
    truns(
    "init segment --collapseInitP T --save_rdata \
    --notrain -m tmpdir/model.txt -c tests/data/randomCounts.txt \
    -r tests/data/regions.bed -n 10 --outdir tmpdir")
    
    #check whether the rdata was created
    expect_true(file.exists(gsub("tmpdir", tmpdir, "tmpdir/rdata.Rdata")))
    
    #REPORT
    truns(
    "init report -m tmpdir/model.txt 
    -c tests/data/randomCounts.txt -s tmpdir/segmentation.bed 
    --outdir tmpdir")
    
    #let's try the colors option
    truns(
    "init report --colors tmpdir/colors.txt 
    -m tmpdir/model.txt -c tests/data/randomCounts.txt 
    -s tmpdir/segmentation.bed --outdir tmpdir")
    #let's try the labels option
    
    writeLines("ale\nmarco\nandrea\ngiacomo\nanna\njuliane\nbarbara\nisabella\nfrancesco\nchiara", 
    file.path(tmpdir, "labels.txt"))
    truns(
    "init report --labels tmpdir/labels.txt 
    -m tmpdir/model.txt -c tests/data/randomCounts.txt 
    -s tmpdir/segmentation.bed --outdir tmpdir")
    
    #let's try multiple annotations
    truns(
    "init report 
    -a genes:bed1 --annot mansegm:bed2
    -m tmpdir/model.txt -c tests/data/randomCounts.txt 
    -s tmpdir/segmentation.bed --outdir tmpdir")
    
    expect_true(file.exists(cmd_sc("tmpdir/annot_genes.txt")))
    expect_true(file.exists(cmd_sc("tmpdir/annot_mansegm.txt")))
    
})


test_that("normalizecounts",{
    #make 3 count matrices
    os <- list(c(1,2,3), c(2,1,3), c(3,2,1))
    targets <- paste0("tmpdir/cmat", 1:3, ".txt")
    for (i in seq_along(os)){
        o <- os[[i]]
        mline <- paste0(collapse=" ", "-m ", 
        c("H3K4me3", "H3K36me", "H3K9me3"), ":bam", o)
        truns(paste("init getcounts -r regs.bed -t ", targets[i], mline))
    }
    
    #run normalize counts
    cline <- paste0(collapse=" ", "-c ", targets)
    truns(paste0("init normalizecounts ", cline))
    
    #check that the new matrices have been created according to the suffix
    newtargets <- gsub(".txt$", paste0(defSuffix, ".txt"), targets)
    expect_true(all(file.exists(cmd_sc(newtargets))))
    
    #check that using the 'triggerOverwrite' suffix no new files are created 
    #(in the temporary directory)
    lfOld <- list.files(tmpdir)
    
    cline <- paste0(collapse=" ", "-c ", newtargets)
    truns(paste0("init normalizecounts ", cline, " -s ", triggerOverwrite))
    expect_true(setequal(list.files(tmpdir), lfOld))
})

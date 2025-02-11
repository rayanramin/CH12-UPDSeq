#!/usr/bin/env Rscript
# Retrieve arguments passed to the script
.libPaths(.libPaths()[grep(".conda",.libPaths())])
args <- commandArgs(trailingOnly = TRUE)
# Initialize an empty list to store argument values
arg_list <- list()

error_message  <- "use -h or --help to see the usage"
help_message <- "Usage: Rscript script.R --bams=[] --regions=[] --ref=[] --dbsnp=[] --contexts=[] --outdir=[] --prefix=[]\n
bams = string path/to/bam/file or list.txt \n
regions = string chr:from-to or Regions.bed \n
ref = string path/to/reference file\n
dbsnp = optional; string path/to/dbsnp file\n
contexts = optional; string list of contexts to analyze; default is 'NC,WRCY,nWRCY'\n
outdir = optional; string path/to/output/directory; default is PWD\n
prefix = optional; string prefix for output files; default is 'output'\n"



# Loop over arguments and parse out the name and value
for (arg in args) {
  if(arg == "-h" | arg == "--help") { cat(help_message); quit() }
  # Split argument into name and value
  splitted <- strsplit(sub("^--", "", arg), "=")[[1]]
  if (length(splitted) == 2) {
    arg_list[[splitted[1]]] <- splitted[2]
  }
}

if( is.null(arg_list[["bams"]])) {
    stop(paste("bams file not provided","\n", error_message))
} else { 
    if( grepl(".txt$", arg_list[["bams"]]) & file.exists(arg_list[["bams"]])) {
        bams <- readLines(arg_list[["bams"]])
        for(i in 1:length(bams)) {
            if( !file.exists(bams[i])) {
                stop(paste("Error parsing bam file.", bams[i], " does not exist", "\n", error_message))
            }
        }
    } else if( grepl(".bam$", arg_list[["bams"]])) {
        bams <- arg_list[["bams"]]
        if(! file.exists(bams)) {
            stop(paste("Error parsing bam file.", bams, " does not exist", "\n", error_message))
        }
    } else {
        stop(paste("Error parsing bam file.","\n", error_message))
    }
}
cat(length(bams), " bams loaded\n")
print(bams)

if( is.null(arg_list[["regions"]])) {
    stop(paste("regions file not provided","\n", error_message))
} else { 
    if( grepl(".bed$", arg_list[["regions"]]) & file.exists(arg_list[["regions"]])) {
        regions <- read.table(arg_list[["regions"]], header = FALSE, sep = "\t")
        if(ncol(regions) < 3) {
            stop(paste("Error parsing regions","\n", error_message))
        }
        names(regions)[1:3] <- c("chr", "start", "end")
        if(ncol(regions) >=4) {names(regions)[4] <- "regions"} else {regions$regions <- paste0("region_", 1:nrow(regions))}
    } else if( grepl("^.*:[0-9]+-[0-9]+$", arg_list[["regions"]])) {
        tmp <- arg_list[["regions"]]
        regions <- data.frame(chr= substr(tmp, 1, regexpr(":", tmp)-1), 
                   start= as.numeric(sub("^.*:", "", sub("-.*$", "", tmp))), 
                   end= as.numeric(sub("^.*-", "", tmp)))
    } else {
        stop(paste("Error parsing regions","\n", error_message))
    }
}
cat(nrow(regions), " regions loaded\n") 

if( is.null(arg_list[["ref"]])) {
    stop(paste("ref file not provided","\n", error_message))
} else { 
    if( file.exists(arg_list[["ref"]])) {
        ref <- arg_list[["ref"]]
    } else {
        stop(paste("Error accessing ref file","\n", error_message))
    }
}
cat("ref: ", ref, "\n")

if(!is.null(arg_list[["dbsnp"]])) {
    if( file.exists(arg_list[["dbsnp"]])) {
        dbsnp <- arg_list[["dbsnp"]]
    } else {
        stop(paste("Error accessing dbsnp file","\n", error_message))
    }
} else {
    dbsnp <- NULL
}

if(!is.null(arg_list[["contexts"]])) {
    CONTEXTS <- strsplit(arg_list[["contexts"]], ",")[[1]]
} else {
    CONTEXTS <- c("NC","WRCY","nWRCY")
}

cat("CONTEXTS: ", CONTEXTS, "\n")

if(!is.null(arg_list[["outdir"]])) {
    if( dir.exists(arg_list[["outdir"]])) {
        outdir <- arg_list[["outdir"]]
    } else {
        stop(paste("Error parsing outdir","\n", error_message))
    }
} else {
    outdir <- getwd()
}
outdir <- normalizePath(outdir)

if(!is.null(arg_list[["prefix"]])) {
    prefix <- arg_list[["prefix"]]
} else {
    prefix <- "output"
}
############################################
suppressMessages({
    library(dplyr, quietly=TRUE) ;
    library(stringr, quietly=TRUE)
})
############################################
robust.system <- function (cmd) {
    stdoutFile = tempfile(pattern="R_robust.system_stdout", fileext=as.character(Sys.getpid()))
    retval = list()
    retval$exitStatus = system(paste0(cmd, " > ", shQuote(stdoutFile)))
    retval$stdout = readLines(stdoutFile)
    unlink(stdoutFile)
    return(retval)
}
read.readcount <- function(X){
    read.table(X, sep="\t", col.names=c("Chr.","POS","REF","COV","A","C","G","T"))
}
Detect_Context <- function(ref_fasta,chr,pos,context){
  if (missing(ref_fasta)) {
    stop("ERROR: input ref_fasta is required")}
  if (missing(chr)) {
    stop("ERROR: input chr is required")}
  if (missing(pos)) {
    stop("ERROR: input pos is required")}
  if (missing(context)) {
    stop("ERROR: input context is required")}
  if(stringr::str_count(context,"C")>1 & stringr::str_count(context,"X")==0){
    stop('ERROR: wrong context, if more than 1 "C" in the context, you should specify the target "C" by an "X"')}
  if(stringr::str_count(context,"X")>1){
    stop('ERROR: wrong context, you can only use 1 "X" to denote the target "C"')}
  if(stringr::str_count(context,"X")==1){
  context <- toupper(context)
    l <- nchar(context)
    xp <- stringr::str_locate(context,"X")[1]
    ba = l-xp; bb = xp-1;
    if(bb < ba){ context <- paste0(paste0(rep("N",ba-bb),collapse = ""),context)}
    if(bb > ba){ context <- paste0(context, paste0(rep("N",bb-ba),collapse = ""))}
    context <- stringr::str_replace(context,"X","C")
  } else if(stringr::str_count(context,"X")==0){
    l <- nchar(context)
    xp <- stringr::str_locate(context,"C")[1]
    ba = l-xp; bb = xp-1;
    if(bb < ba){ context <- paste0(paste0(rep("N",ba-bb),collapse = ""),context)}
    if(bb > ba){ context <- paste0(context, paste0(rep("N",bb-ba),collapse = ""))}
  }
  flank <-  max(ba,bb)
  context <- bioseq::dna(context)
  revcom_context <-  bioseq::seq_complement(bioseq::seq_reverse(context))
  GenomicRanges::GRanges(paste0(chr,":",pos-flank,"-",pos+flank)) -> g
  Rsamtools::scanFa(ref_fasta, param=g)  %>%
  as.character %>% unname %>% bioseq::dna() %>%
    {case_when(bioseq::seq_detect_pattern(.,pattern = context) ~ 1,
               bioseq::seq_detect_pattern(.,pattern =  revcom_context) ~ -1,
               TRUE ~ 0)}
}
UX_v4.5 <- function (X, test_vector) {
    cntxt <- test_vector
	#####
	nC2T <- ifelse(cntxt == 1, as.integer(X$T), 0)
    t <- nC2T/X$COV
    nc <- sum(cntxt == 1)
    AllC <- sum(t)/nc * 1000
	#####
	nG2A <- ifelse(cntxt == -1, as.integer(X$A), 0)
    a <- nG2A/X$COV
    ng <- sum(cntxt == -1)
    AllG <- sum(a)/ng * 1000
	#####
    u <- (a + t)
    z <- sum(cntxt %in% c(-1, 1))
    AllU <- sum(u)/z * 1000
    eps <- 0.001
    bias <- log((AllC + eps)/(AllG + eps), 2)
	#####
    return(list(
		Number_of_available_context_top_strand = nc, 
        Number_of_available_context_bottom_strand = ng,
		Number_of_available_context_both_strands = z,
		################
		Number_of_CtoT_mismatches = sum(nC2T),
		Number_of_GtoA_mismatches = sum(nG2A),
		################
		Nubmer_of_sites_with_any_CtoT_mismatches = sum(nC2T>0),
		Number_of_sites_with_any_GtoA_mismatches = sum(nG2A>0),
		################
		CtoT_UI = AllC, 
		GtoA_UI = AllG, 
		UI = AllU, 
        strand_bias = bias
		))
}
############################################
# Check if bam-readcount is available on the system
if (suppressWarnings({! grepl("bam-readcount version", robust.system("bam-readcount --version")$stdout)})) {
    stop("bam-readcount not available on the system")	
}

write.table(regions, file=paste0(outdir,"/TMP_regions.bed"), sep="\t", quote=F, row.names=F, col.names=F)
sink(paste0(outdir,"/TMP_brc.sh"))
brc_files <- c()
for(i in 1:length(bams)){
    bai <- c(str_replace(normalizePath(bams[i]), ".bam$",".bai"), paste0(normalizePath(bams[i]), ".bai"))
    if(!file.exists(bai[1]) & !file.exists(bai[2])){
        cat(paste("index file not available for ", bams[i], "--> skipping. \n" ),file = stderr())
        next
    }
    brc_outfile <- paste0(outdir,"/",prefix,"_", str_replace(basename(bams[i]), ".bam$",".brc"))
    brc_files <- c(brc_files, brc_outfile)
    if(file.exists(brc_outfile)){
        cat(paste("brc file already exists for ", bams[i], "--> skipping. \n" ),file = stderr())
        next
    }
	cat(paste0("bam-readcount -w0 -l ", outdir, "/TMP_regions.bed -f ",normalizePath(ref)," ", normalizePath(bams[i]), " | awk -F \":|\\t|=\" 'BEGIN {OFS = \"\\t\"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > ", brc_files[i]))
    cat("\n \n")
}
sink()
cat("running bam-readcount\n")
suppressWarnings({brc_job=robust.system(paste0("sh ", outdir, "/TMP_brc.sh"))})
if(brc_job$exitStatus != 0){
    stop("bam-readcount failed")
}
system(paste0("rm ", outdir, "/TMP_brc.sh"))
system(paste0("rm ", outdir, "/TMP_regions.bed"))
cat("bam-readcount done\n")
############################################

Sample_Names <- gsub(".brc","",basename(brc_files))

Xdbsnp <- NULL
if(!is.null(dbsnp)){ 
    cat("loading dbsnp\n")
    Xdbsnp <- read.table(dbsnp, header = FALSE, sep = "\t", comment.char = "#")
}
## import igh brc data

cat("processing data\n")

Data <- list()
for(i in 1:length(brc_files)){
	read.readcount(brc_files[i] ) %>%
    mutate( sample = Sample_Names[i] , REF=toupper(REF)) %>%
    filter( REF %in% c("C","G")) -> tmp
    if(!is.null(Xdbsnp)){
        tmp[which(!(paste(tmp$Chr., tmp$POS) %in% paste(Xdbsnp$V1, Xdbsnp$V2))),] -> tmp
    }
    tmp ->  Data[[Sample_Names[i]]]
}

Data.summary <- data.frame()

for(i in seq_along(Data)){
    for(j in seq_along(regions)){
        Data[[Sample_Names[i]]] %>% filter(POS >= regions$start[j] & POS <= regions$end[j]) -> tmp
        if(nrow(tmp)>0){
            for(context in CONTEXTS ){
                if(str_detect(context,"^n")){ 
                tmp %>% {ifelse(Detect_Context(ref_fasta=ref, chr=.$Chr., pos=.$POS, context=str_remove(context,"^n"))!=0,0,Detect_Context(ref_fasta=ref, chr=.$Chr., pos=.$POS, context="C"))} -> condet
                } else {
                tmp %>% {Detect_Context(ref_fasta=ref, chr=.$Chr., pos=.$POS, context=context)} -> condet
                }
                tmp %>% UX_v4.5(.,condet) %>% 
                as.data.frame %>% mutate(Sample = Sample_Names[i], regions = regions$region[j] , context=context) -> tmp2
                if(nrow(Data.summary)==0){
                    Data.summary <- tmp2
                } else {
                    Data.summary <- rbind(Data.summary, tmp2)
                }
            }
        }
    }
} 

cat("writing output to ", paste0(outdir,"/",prefix,".txt"), "\n")
left_join(regions,Data.summary, by="regions") %>% 
write.table(paste0(outdir,"/",prefix,".txt"),sep="\t", col.names = T, row.names = F, quote = F )

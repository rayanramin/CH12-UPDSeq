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
contexts = optional; string list of contexts to check;'\n
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
        dbsnp <- normalizePath(dbsnp)
    } else {
        stop(paste("Error accessing dbsnp file","\n", error_message))
    }
} else {
    dbsnp <- NULL
}


if(!is.null(arg_list[["contexts"]])) {
    CONTEXTS <- strsplit(arg_list[["contexts"]], ",")[[1]]
    cat("CONTEXTS: ", CONTEXTS, "\n")
} else {
    CONTEXTS <- NULL
}


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
    library(stringr, quietly=TRUE);
    library(GenomicRanges, quietly=TRUE);
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
UX <- function(Y) {
    X <- Y
    # defining conditions    
    cond1 <- logical(nrow(X)) ; cond1 <- X$REF == "C"
    cond2 <- logical(nrow(X)) ; cond2 <- X$REF == "G"
    # C to T changes only
    t <- ifelse(cond1, as.integer(X$T), 0) / X$COV
    # G to A changes only
    a <- ifelse(cond2, as.integer(X$A), 0) / X$COV
    # All changes
    c(a+t) 
}
get_mm_data <- function(X){
	for(i in 1:length(X[[1]]$seq)){
	seq <- X[[1]]$seq[i]
	cigar <- X[[1]]$cigar[i]
	mdz <- X[[1]]$tag$MD[i]
	pos <- X[[1]]$pos[i]

	x <- NA
	gsub("([\\^]*[ACGT]+)[0]*", " \\1 ", mdz) %>%
	lapply(., function(i) strsplit(i, "[ ]+")) %>% unlist() %>%
	lapply(., function(y) { 
			if (suppressWarnings(!is.na(as.numeric(y)))) {
				o <- rep("M", as.numeric(y))
			} else if( length(grep("\\^", y)) > 0) {
				o <- rep("D", nchar(y) - 1)
			} else if (nchar(y) == 1 &  (y %in% c("A","T","C","G"))) {
				o <- rep(y, 1)
			} else if (nchar(y) == 1 & !(y %in% c("A","T","C","G"))) {
				o <- rep("N", 1)
			}  
		}) %>% do.call(c, .) -> x

	matches <- str_match_all(cigar, "(\\d+)([MIDNSHP=X])")[[1]]
	cigar_data <- list(lengths = as.integer(matches[, 2]), types = matches[, 3])

	### sanity check 1
	if (sum(cigar_data$lengths[cigar_data$type %in% c("S","I","M","=","X")]) != nchar(seq)){stop("Error in parsing CIGAR string")}
	### sanity check 2
	if (sum(cigar_data$lengths[cigar_data$type %in% c("M","D")]) != length(x)){stop("Error in parsing MD tag")}

	gpos <- pos # genome position
	rpos <- 1 # read position # can be different that 1-length(mdz) if there are soft clips
	xpos <- 1 # query position

	mismatches <- data.frame( chr = character(),gpos = numeric(), rpos = numeric(), base = character(), alt = character()) # Initialize mismatches as an empty data frame
	# Loop through each CIGAR operation
	for (k in seq_along(cigar_data$types)) {
		op <- cigar_data$types[k]
		len <- cigar_data$lengths[k]

		if(op == "S") { #soft clip
			rpos <- rpos + len # move read position
		} else if(op %in% c("M","=","X")) { #match
			for (j in seq_len(len)) {
				if (x[xpos] %in% c("A","C","G","T")) {
					mismatch <- data.frame(chr = X[[1]]$rname[i], gpos = gpos, rpos = rpos, base = x[xpos], alt = substr(seq, rpos, rpos))
					rbind(mismatches, mismatch) -> mismatches
				} 
				gpos <- gpos + 1
				rpos <- rpos + 1
				xpos <- xpos + 1
			}
		} else if ( op == "I" ){ #insertion in the read
			rpos <- rpos + len
			# xpos <- xpos + len ##(‘H’, ‘S’, ‘P’, ‘N’, and ‘I’ CIGAR operations) are not represented in the MD string.
		} else if ( op == "D" ){ #deletion in the read
			gpos <- gpos + len
			xpos <- xpos + len
		}
	}
	## sanity check 3
	# if (nrow(mismatches) > 0){
	# mismatches$REF <- NA
	# for (m in seq_len(nrow(mismatches))) {
	# mismatches$REF[m] <- robust.system(paste0("module load samtools; samtools faidx ", ref_fasta, " ", X[[1]]$rname[i], ":", mismatches$gpos[m], "-", mismatches$gpos[m]))$stdout[2] %>% toupper()
	# }
	# if(all(mismatches$REF == mismatches$base)){ cat("check pass \t",i,"\n")}else{stop("check fail \t")}
	# ## in test runs it returns no error. comment out to speed up
	# }
	##########

	if(nrow(mismatches)>0){ mismatches$readnum <- i}
	if (i == 1) {out <- mismatches } else {out <- rbind(out,mismatches)}
	}
	return(out)
}
get_GC2_mm_distance <- function(X){
	X$C2T <- ifelse(X$base == "C" & X$alt == "T", T, F)
	X$G2A <- ifelse(X$base == "G" & X$alt == "A", T, F)
	X$CGmmDist <- NA
	X$GCmmStrand <- NA
	X$mid <- 1:nrow(X)

	lapply(1:length(unique(X$readnum)), function(i){
		subset(X, readnum == i & (C2T | G2A)) -> tmp
		if (nrow(tmp) > 1){ 
			tmp$CGmmDist  <- c(NA, diff(tmp$rpos))
			for(j in 2:nrow(tmp)){
				if ((tmp$C2T[j-1] & tmp$C2T[j]) | (tmp$G2A[j-1] & tmp$G2A[j])){
					tmp$GCmmStrand[j] <- "same"
				} else if(((tmp$C2T[j-1] & tmp$G2A[j]) | (tmp$G2A[j-1] & tmp$C2T[j]))){
					tmp$GCmmStrand[j] <- "opposite"
				}
			}}
			tmp
		}) %>% do.call(rbind, .) -> tmp
	X$CGmmDist <-  mapacross(X$mid, tmp$mid, tmp$CGmmDist, NA)
	X$GCmmStrand <- mapacross(X$mid, tmp$mid, tmp$GCmmStrand, NA)
	return(X)
}
mapacross <- function(x, A, B, m=NA) {
  ifelse(x %in% A, B[match(x, A)], m)
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

Sample_Names <- gsub(paste0(prefix,"_"),"",gsub(".brc","",basename(brc_files)))

lapply(brc_files, function(x) read.readcount(x) ) -> brc_data

snp_min_cov <- 5
snp_min_AF <- 0.5
lapply(brc_data, function(x) {x$U <- UX(x);  filter(x,COV > snp_min_cov & U >=snp_min_AF)}) %>% 
    do.call(rbind, .) %>% 
    filter(!duplicated(paste(.$Chr,.$POS))) -> snps

Xdbsnp <- NULL
if(!is.null(dbsnp)){ 
    cat("loading dbsnp\n")
    Xdbsnp <- read.table(dbsnp, header = FALSE, sep = "\t", comment.char = "#")
}
## import igh brc data

cat("processing data ...\n")

Data <- list()
for (i in seq_along(brc_files)){
    brc_data[[i]] %>%
    mutate(REF=toupper(REF)) -> tmp
    # filter(tmp, REF %in% c("C","G")) -> tmp
    if(nrow(snps) > 0){
    tmp[which(!(paste(tmp$Chr., tmp$POS) %in% paste(snps$Chr., snps$POS))),] -> tmp
    }
    if(!is.null(Xdbsnp)){
        tmp[which(!(paste(tmp$Chr., tmp$POS) %in% paste(Xdbsnp$V1, Xdbsnp$V2))),] -> tmp
    }
    tmp$U <- UX(tmp)
    tmp ->  Data[[Sample_Names[i]]]
}

X <- list()
for (i in seq_along(brc_files)){
    Sample_Names[i] -> sname
	print(sname)
	Data[[sname]] -> tmpbrc
	X[[ sname ]] <- list()
	for(j in 1:nrow(regions)){
		print(regions$region[j])
		param <- Rsamtools::ScanBamParam(what = c("seq","rname","pos","cigar","qwidth","mapq"),tag = c("NM","MD"),
								 which=GRanges("chr12", IRanges(start=regions$start[j],end=regions$end[j])))
		Rsamtools::scanBam(bams[grep(sname,bams)] , param = param ) -> tmp
		if(length(tmp[[1]]$seq) == 0){next}
		get_mm_data(tmp) -> tmpMM0
        if(!is.null(Xdbsnp)){
            tmpMM0[which(!(paste(tmpMM0$chr, tmpMM0$gpos) %in% paste(Xdbsnp$V1, Xdbsnp$V2))),] -> tmpMM0
        }
		filter(tmpMM0 , gpos > regions$start[j] & gpos < regions$end[j] ) %>% # remove mm outside the region
        # filter(base %in% c("C","G")) %>% 
		get_GC2_mm_distance() %>% 
		mutate( region= regions$region[j], sample = sname, nreads_in_region = length(tmp[[1]]$seq) ) -> tmpMM
		### getting depth of coverage for these positions
        inner_join(
            dplyr::rename(tmpMM, "CHR"="chr", "POS"="gpos"),
            dplyr::rename(tmpbrc, "CHR"="Chr."),
            by=c("CHR","POS")) -> MM
        ### check for contexts
        if(!is.null(CONTEXTS)){
            for(context in CONTEXTS ){
                if(str_detect(context,"^n")){ 
                MM %>% {ifelse(Detect_Context(ref_fasta=ref, chr=.$CHR, pos=.$POS, context=str_remove(context,"^n"))!=0,0,Detect_Context(ref_fasta=ref, chr=.$CHR, pos=.$POS, context="C"))} -> condet
                } else {
                MM %>% {Detect_Context(ref_fasta=ref, chr=.$CHR, pos=.$POS, context=context)} -> condet
                }
                MM[[paste0("context_is_",context)]] <- NA
                MM[[paste0("context_is_",context)]] <- condet
            }
        }
		## sanity check
		if( ! all(MM$REF == MM$base)){stop("Error in joining brc data")}
		select(MM, sample, region,read_id=readnum, chr=CHR, pos=POS, ref=base, alt,rpos, C2T, G2A, CGmmDist, GCmmStrand,mm_id=mid,nreads_in_region, COV,A,C,G,T,U, starts_with("context_is_")) ->  X[[ sname ]][[regions$region[j]]]
	}
}
cat("done\n")


cat("writing output to ", paste0(outdir,"/",prefix,".txt"), "\n")
do.call(rbind, lapply(X, function(a) do.call(rbind, a)))  %>% 
write.table(paste0(outdir,"/",prefix,".txt"),sep="\t", col.names = T, row.names = F, quote = F )

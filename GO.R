#!/usr/bin/env Rscript

options(stringsAsFactors=F)
pseudocount = 1e-04

#---------------------------------------------------------
# OPTIONS
#---------------------------------------------------------

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-u", "--universe"), help="a list of gene identifiers (ensEMBL ids), NO header"),
make_option(c("-G", "--genes"), default="stdin",
	help="a list of gene identifiers for the foreground (ensEMBL ids), WITHOUT header [default=%default]"),
make_option(c("-c", "--categ"), help="choose the GO category < BP | MF | CC > [default=%default]", default="BP"),
make_option(c("-s", "--species"), help="choose the species < human | mouse | dmel | danio > [default=%default]", default="human"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=%default]", default="out"),
make_option(c("-d", "--output_dir"), default="./", help="directory for the output [default=%default]"),
make_option(c("-v", "--verbose"), default=FALSE, action="store_true", help="verbose")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}



#---------------------------------------------------------
# LIBRARIES
#--------------------------------------------------------- 

if (opt$verbose) {cat("Loading libraries... ")}


suppressPackageStartupMessages(library("GO.db"))
if (opt$species == "human") {suppressPackageStartupMessages(library("org.Hs.eg.db"))}
if (opt$species == "mouse") {suppressPackageStartupMessages(library("org.Mm.eg.db"))}
if (opt$species == "dmel") {suppressPackageStartupMessages(library("org.Dm.eg.db"))}
if (opt$species == "danio") {suppressPackageStartupMessages(library("org.Dr.eg.db"))}
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))

if (opt$verbose) {cat("DONE\n\n")}

#---------------------------------------------------------
# THE CODE
#---------------------------------------------------------

U = read.table(opt$universe, h=F, col.names='hs')

if (opt$genes == "stdin") {
	G = read.table(file("stdin"), h=F, col.names='hs')	
} else {
	G = read.table(opt$genes, h=F, col.names='hs')
}

if (opt$species == "human") {
universe = unlist(mget(U$hs, org.Hs.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "mouse") {
universe = unlist(mget(U$hs, org.Mm.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "dmel") {
universe = unlist(mget(U$hs, org.Dm.egENSEMBL2EG, ifnotfound=NA))}

if (opt$species == "danio") {
universe = unlist(mget(U$hs, org.Dr.egENSEMBL2EG, ifnotfound=NA))}

if (opt$verbose) {sprintf("%s background genes; %s with a corresponding entrez id", nrow(U), length(unique(universe)))}

createParams = function(x, species="human") {
	if (species == "human") {
		ann = "org.Hs.eg.db"
		geneset = unlist(mget(x, org.Hs.egENSEMBL2EG, ifnotfound=NA))
	}
	if (species == "mouse") {
		ann = "org.Mm.eg.db"
		geneset = unlist(mget(x, org.Mm.egENSEMBL2EG, ifnotfound=NA))
	}
	if (species == "dmel") {
		ann = "org.Dm.eg.db"
		geneset = unlist(mget(x, org.Dm.egENSEMBL2EG, ifnotfound=NA))
	}
        if (species == "danio") {
		ann = "org.Dr.eg.db"
		geneset = unlist(mget(x, org.Dr.egENSEMBL2EG, ifnotfound=NA))
	}
	sprintf("%s foreground genes; %s with a corresponding entrez id", length(x), length(unique(geneset)))
	pv = 1-(1-0.05)**(1/length(x))
	params = new("GOHyperGParams",
		geneIds = geneset,
		universeGeneIds = universe,
		annotation = ann,
		ontology = opt$categ,
		pvalueCutoff = pv,
		conditional = TRUE,
		testDirection='over')
	return(params)}

res = hyperGTest(createParams(G$hs, opt$species))

df = summary(res)
df$Pvalue = format(df$Pvalue, digits=1) 
df$OddsRatio <- round(df$OddsRatio, 2)
df$ExpCount <- round(df$ExpCount, 2)

#---------------------------------------------------------
# OUTPUT
#---------------------------------------------------------
output = sprintf("%s/%s.%s", opt$output_dir, opt$output, opt$categ)
write.table(df, file=sprintf("%s.tsv", output), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.html", output))

q(save='no')



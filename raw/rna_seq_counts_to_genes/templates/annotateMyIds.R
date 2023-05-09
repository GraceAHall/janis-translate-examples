
args <- commandArgs(trailingOnly=TRUE)

options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

id_type <- "ENTREZID"
organism <- "Mm"
output_cols <- "ENTREZID,SYMBOL,GENENAME"
file_has_header <- TRUE
remove_dups <- FALSE

ids <- as.character(read.table(args[1], header=file_has_header)[,1], sep="\t", quote="")

if(organism == "Hs"){
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db <- org.Hs.eg.db
} else if (organism == "Mm"){
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    db <- org.Mm.eg.db
} else if (organism == "Dm"){
    suppressPackageStartupMessages(library(org.Dm.eg.db))
    db <- org.Dm.eg.db
} else if (organism == "Dr"){
    suppressPackageStartupMessages(library(org.Dr.eg.db))
    db <- org.Dr.eg.db
} else {
    cat(paste("Organism type not supported", organism))
}

cols <- unlist(strsplit(output_cols, ","))
result <- select(db, keys=ids, keytype=id_type, columns=cols)

if(remove_dups) {
    result <- result[!duplicated(result$args[1]),]
}

write.table(result, file=args[2], sep="\t", row.names=FALSE, quote=FALSE)
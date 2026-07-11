context("src_organism-select")

suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(txdbmaker)
})

.loadLightTxDb <- function(dbpath) {
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), dbpath)
    tx_df <- RSQLite::dbReadTable(conn, "ranges_tx")
    exon_df <- RSQLite::dbReadTable(conn, "ranges_exon")
    cds_df <- RSQLite::dbReadTable(conn, "ranges_cds")
    seq_df <- RSQLite::dbReadTable(conn, "seqinfo")
    RSQLite::dbDisconnect(conn)
    transcripts <- unique(tx_df[, c("tx_id", "tx_name", "tx_chrom", "tx_strand", "tx_start", "tx_end")])
    splicings <- merge(exon_df, cds_df, by=c("tx_id", "exon_rank"), all.x=TRUE)
    splicings <- splicings[, c("tx_id", "exon_rank", "exon_id", "exon_start", "exon_end", "cds_id", "cds_start", "cds_end")]
    chrominfo <- seq_df[, c("seqnames", "seqlengths", "isCircular")]
    colnames(chrominfo) <- c("chrom", "length", "is_circular")
    chrominfo$is_circular <- as.logical(chrominfo$is_circular)
    genes <- unique(tx_df[!is.na(tx_df$entrez), c("tx_id", "entrez")])
    colnames(genes) <- c("tx_id", "gene_id")
    genes$gene_id <- as.character(genes$gene_id)
    txdb <- suppressWarnings(makeTxDb(transcripts, splicings, genes=genes, chrominfo=chrominfo))
    suppressWarnings(GenomeInfoDb::genome(txdb) <- "hg38")
    txdb
}

hg38light <- hg38light()
txdb <- .loadLightTxDb(hg38light)
src <- src_organism(dbpath=hg38light)

test_that("keytypes", {
    expect_equal(class(keytypes(src)), "character")
    expect_true(length(keytypes(src)) > 0)
})

test_that("columns", {
    expect_equal(class(columns(src)), "character")
    expect_true(length(columns(src)) > 0)
})

test_that("keys", {
    expect_error(keys(src, "foo"))
    expect_equal(class(keys(src)), "character")
    expect_true(length(keys(src)) > 0)
    expect_true(all(keys(src, "tx_id") %in% keys(txdb, "TXID")))
})

test_that("select", {
    columns_src <- c("entrez", "tx_id", "tx_name","exon_id")
    keytype_src <- "tx_name"
    columns_txdb <- c("GENEID", "TXID", "TXNAME","EXONID")
    keytype_txdb <- "TXNAME"

    keys <- head(keys(src, keytype_src))
    rs_src <- select(src, keys, columns_src, keytype_src) %>% collect()
    rs_txdb <- select(txdb, keys, columns_txdb, keytype_txdb)

    expect_equal(dim(rs_src), dim(rs_txdb))
    expect_equal(
        rs_src[order(rs_src[,keytype_src]),]$tx_id,
        rs_txdb[order(rs_txdb[,keytype_txdb]),]$TXID
    )
})

test_that("mapIds", {
    keys <- head(keys(src, "tx_name"))

    rs_src <- mapIds(src, keys, "tx_id", "tx_name")
    rs_txdb <- mapIds(txdb, keys, "TXID", "TXNAME")

    expect_equal(rs_src, rs_txdb)
})

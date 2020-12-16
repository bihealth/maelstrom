require(logger)
require(cn.mops)
require(genomation)

log_threshold(DEBUG)
log_info('cnmops_call.R script starting up')

# load input ------------------------------------------------------------------

log_info('reading TSV file from {snakemake@input[["tsv"]]}...')
input = readGeneric(
    snakemake@input[["tsv"]],
    header=TRUE,
    keep.all.metadata=TRUE
)
log_info('... successfully loaded {length(input)} records')

# perform CNV calling ---------------------------------------------------------

log_info('running cn.mops() with {snakemake@threads} threads...')
result = cn.mops(
    input=input,
    parallel=snakemake@threads
)
log_info('... done running cn.mops()')

log_info('computing integer copy numbers...')
integer_result = calcIntegerCopyNumbers(result)
log_info('... done computing integer copy numbers')

# write output ----------------------------------------------------------------

log_info('writing output files...')
segm = as.data.frame(segmentation(integer_result))
cnvs = as.data.frame(cnvs(integer_result))
cnv_regions = as.data.frame(cnvr(integer_result))

log_info(' => to {snakemake@output[["segmentation"]]}')
write.table(
    segm,
    file=snakemake@output[["segmentation"]],
    sep="\t"
)
log_info(' => to {snakemake@output[["cnvs"]]}')
write.table(
    cnvs,
    file=snakemake@output[["cnvs"]],
    sep="\t"
)
log_info(' => to {snakemake@output[["cnv_regions"]]}')
write.table(
    cnv_regions,
    file=snakemake@output[["cnv_regions"]],
    sep="\t"
)
log_info('done writing output files')

# compute checksums -----------------------------------------------------------

log_info('computing md5 sums...')
system(
    sprintf(
        "cd %s; md5sum %s >%s.md5; cd -",
        dirname(snakemake@output[["segmentation"]]),
        basename(snakemake@output[["segmentation"]]),
        basename(snakemake@output[["segmentation"]])
    )
)
system(
    sprintf(
        "cd %s; md5sum %s >%s.md5; cd -",
        dirname(snakemake@output[["cnvs"]]),
        basename(snakemake@output[["cnvs"]]),
        basename(snakemake@output[["cnvs"]])
    )
)
system(
    sprintf(
        "cd %s; md5sum %s >%s.md5; cd -",
        dirname(snakemake@output[["cnv_regions"]]),
        basename(snakemake@output[["cnv_regions"]]),
        basename(snakemake@output[["cnv_regions"]])
    )
)
log_info('... done computing md5 sums')

log_info('All done. Have a nice day!')

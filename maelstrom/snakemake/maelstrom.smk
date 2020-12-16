import textwrap

# def maelstrom_bam_scan_input(wildcards):
#     return BY_NAME[wildcards.sample].bam


# rule maelstrom_bam_scan:
#     input:
#         maelstrom_bam_scan_input,
#     output:
#         multiext(
#             "work/maelstrom_bam_scan/{sample}/{sample}",
#             ".bam",
#             ".bam.md5",
#             ".bam.bai",
#             ".bam.bai.md5",
#         ),
#     shell:
#         """
#         set -x

#         export TMPDIR=$(mktemp -d)
#         # trap "rm -rf $TMP DIR" ERR EXIT

#         ( \\
#             echo "blocked_regions_bed = \\"{BLOCKLIST}\\""
#             echo "htslib_io_threads = 4"; \\
#         ) > $TMPDIR/maelstrom.toml
#         cat $TMPDIR/maelstrom.toml

#         output=$(echo {output} | tr ' ' '\\n' | grep 'bam$')

#         rm -rf $(dirname $output)
#         mkdir -p $(dirname $output)

#         maelstrom-bam-scan -c $TMPDIR/maelstrom.toml {input} - \\
#         | samtools sort -T $output.by_name.tmp -n -O BAM -o - \\
#         | maelstrom-bam-unique -- - - \\
#         | samtools sort -T $output.by_coord.tmp -O BAM -o $output
#         samtools index $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).bai >$(basename $output).bai.md5
#         """


# rule maelstrom_bam_collect_pesr:
#     input:
#         "work/maelstrom_bam_scan/{sample}/{sample}.bam",
#     output:
#         multiext(
#             "work/maelstrom_bam_collect_pesr/{sample}/{sample}",
#             ".tsv.gz",
#             ".tsv.gz.md5",
#             ".tsv.gz.tbi",
#             ".tsv.gz.tbi.md5",
#         ),
#     shell:
#         """
#         set -x

#         output=$(echo {output} | tr ' ' '\\n' | grep 'tsv.gz$')

#         maelstrom-bam-collect-pesr {input} - \\
#         | maelstrom-bed-sort \\
#         | bgzip -c >$output
#         tabix -p bed -f $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).tbi >$(basename $output).tbi.md5
#         """


# def maelstrom_vcf_standardize_input(wildcards):
#     return BY_NAME[wildcards.sample].calls[wildcards.tool]


# rule maelstrom_vcf_standardize:
#     input:
#         maelstrom_vcf_standardize_input,
#     output:
#         multiext(
#             "work/maelstrom_vcf_standardize/{tool}/{sample}/{sample}",
#             ".vcf.gz",
#             ".vcf.gz.md5",
#             ".vcf.gz.tbi",
#             ".vcf.gz.tbi.md5",
#         ),
#     shell:
#         """
#         set -x

#         output=$(echo {output} | tr ' ' '\\n' | grep 'vcf.gz$')

#         export TMPDIR=$(mktemp -d)
#         # trap "rm -rf $TMP DIR" ERR EXIT

#         ( \\
#             echo "stdvcf_apply_filters = $(if [[ "{tool}" == delly* ]]; then echo true; else echo false; fi) \\
#         ) > $TMPDIR/maelstrom.toml
#         cat $TMPDIR/maelstrom.toml

#         maelstrom-vcf-standardize -c $TMPDIR/maelstrom.toml --tool={wildcards.tool} {input} - \\
#         | bcftools sort -O z -o $output
#         tabix -f $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).tbi >$(basename $output).tbi.md5
#         """


# def maelstrom_vcf_cluster_for_tool_input(wildcards):
#     return [
#         "work/maelstrom_vcf_standardize/{tool}/{sample}/{sample}.vcf.gz".format(
#             tool=wildcards.tool, sample=sample.name
#         )
#         for sample in SAMPLE_INFOS
#     ]


# rule maelstrom_vcf_cluster_for_tool:
#     input:
#         maelstrom_vcf_cluster_for_tool_input,
#     output:
#         multiext(
#             "work/maelstrom_vcf_cluster_for_tool/{tool}/clustered",
#             ".vcf.gz",
#             ".vcf.gz.md5",
#             ".vcf.gz.tbi",
#             ".vcf.gz.tbi.md5",
#         ),
#     shell:
#         """
#         set -x

#         output=$(echo {output} | tr ' ' '\\n' | grep 'vcf.gz$')

#         maelstrom-vcf-cluster --overwrite {input} $output
#         tabix -f $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).tbi >$(basename $output).tbi.md5
#         """


# rule maelstrom_vcf_extract:
#     input:
#         vcf="work/maelstrom_vcf_cluster_for_tool/{tool}/clustered.vcf.gz",
#     output:
#         vcf="work/maelstrom_vcf_extract/{tool}/{sample}/{tool}.{sample}.vcf.gz",
#         vcf_md5="work/maelstrom_vcf_extract/{tool}/{sample}/{tool}.{sample}.vcf.gz.md5",
#         tbi="work/maelstrom_vcf_extract/{tool}/{sample}/{tool}.{sample}.vcf.gz.tbi",
#         tbi_md5="work/maelstrom_vcf_extract/{tool}/{sample}/{tool}.{sample}.vcf.gz.tbi.md5",
#     shell:
#         """
#         set -x

#         output=$(echo {output.vcf} | tr ' ' '\\n' | grep 'vcf.gz$')
#         mkdir -p $(dirname $output)

#         bcftools view -s "{wildcards.sample}-N1-DNA1-WGS1" -O z -o {output.vcf} {input.vcf}

#         tabix -f {output.vcf}
#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).tbi >$(basename $output).tbi.md5
#         """


def maelstrom_vcf_annotate_mem_mb(_wildcards, _input=None, _threads=None, attempt=None):
    return 32000  # 2000 * (2 ** (2 * (attempt - 1)))


rule maelstrom_vcf_annotate:
    input:
        pesr_tsv="work/maelstrom_bam_collect_pesr/{sample}/{sample}.tsv.gz",
        doc_bcf="work/maelstrom_bam_collect_doc/100/{sample}/{sample}.bcf",
        snv_vcf="work/gatk_hc_gvcf_sample_filtered/{sample}/{sample}_filtered.vcf.gz",
        vcf="work/maelstrom_vcf_extract/{tool}/{sample}/{tool}.{sample}.vcf.gz",
    output:
        bcf="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.bcf",
        bcf_md5="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.bcf.md5",
        csi="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.bcf.csi",
        csi_md5="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.bcf.csi.md5",
        doc_tsv="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.doc.tsv",
        doc_tsv_md5="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.doc.tsv.md5",
        baf_tsv="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.baf.tsv",
        baf_tsv_md5="work/maelstrom_vcf_annotate/{tool}/{sample}/{sample}.baf.tsv.md5",
    resources:
        mem_mb=maelstrom_vcf_annotate_mem_mb,
    shell:
        """
        set -x

        maelstrom-vcf-annotate \
            --path-pesr-evidence {input.pesr_tsv} \
            --path-doc-evidence {input.doc_bcf} \
            --path-snv-vcf {input.snv_vcf} \
            --overwrite \
            --sample "{wildcards.sample}" \
            {input.vcf} \
            --path-out-doc-summary {output.doc_tsv} \
            --path-out-baf-snvs {output.baf_tsv} \
            {output.bcf}
        tabix -f {output.bcf}

        pushd $(dirname {output.bcf})
        md5sum $(basename {output.bcf}) >$(basename {output.bcf}).md5
        md5sum $(basename {output.csi}) >$(basename {output.csi}).md5
        md5sum $(basename {output.doc_tsv}) >$(basename {output.doc_tsv}).md5
        md5sum $(basename {output.baf_tsv}) >$(basename {output.baf_tsv}).md5
        """


# rule maelstrom_bam_collect_doc:
#     input:
#         maelstrom_bam_scan_input,
#     output:
#         multiext(
#             "work/maelstrom_bam_collect_doc/{window_length}/{sample}/{sample}",
#             ".bcf",
#             ".bcf.md5",
#             ".bcf.csi",
#             ".bcf.csi.md5",
#         ),
#     shell:
#         """
#         set -x

#         export TMPDIR=$(mktemp -d)
#         trap "rm -rf $TMPDIR" ERR EXIT

#         ( \\
#             echo "blocked_regions_bed = \\"{BLOCKLIST}\\""; \\
#             echo "htslib_io_threads = 4"; \\
#             echo "path_reference_fasta = \\"{REFERENCE}\\""; \\
#             echo ""; \\
#             echo "[collect_doc_config]"; \\
#             echo "min_mapq = 0"; \\
#             echo "min_unclipped = 0.6"; \\
#             echo "count_kind = \\"coverage\\"";\\
#             echo "window_length = {wildcards.window_length}"; \\
#         ) > $TMPDIR/maelstrom.toml
#         cat $TMPDIR/maelstrom.toml

#         output=$(echo {output} | tr ' ' '\\n' | grep 'bcf$')
#         mkdir -p $(dirname $output)

#         maelstrom-bam-collect-doc -c $TMPDIR/maelstrom.toml {input} $output;
#         tabix -f $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).csi >$(basename $output).csi.md5
#         """


# rule maelstrom_doc_summary:
#     input:
#         expand(
#             "work/maelstrom_vcf_annotate/{{tool}}/{sample}/{sample}.doc.tsv",
#             sample=SAMPLES,
#         )
#     output:
#         multiext(
#             "work/maelstrom_doc_summary/{tool}/all",
#             ".tsv",
#             ".tsv.md5",
#         )
#     shell:
#         """
#         set -x

#         output=$(echo {output} | tr ' ' '\\n' | grep 'tsv$')

#         paste {input} \
#         | awk -F $'\\t' 'BEGIN {{ OFS=FS; }} {{
#                 printf("%s", $1);
#                 for (i = 2; i <= NF; i += 2) {{
#                     printf("\\t%s", $i);
#                 }}
#                 printf("\\n");
#             }}' \
#         > $output

#         pushd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         """


# TODO: import sorted by SVID instead of by sample...
rule maelstrom_baf_gather:
    input:
        expand(
            "work/maelstrom_vcf_annotate/{{tool}}/{sample}/{sample}.baf.tsv", sample=SAMPLES,
        ),
    output:
        multiext(
            "work/maelstrom_baf_gather/{tool}/all", ".db", ".db.md5",
        ),
    shell:
        textwrap.dedent(
            """
        set -x

        export TMPDIR=$(mktemp -d)

        output=$(echo {output} | tr ' ' '\\n' | grep 'db$')
        input_1=$(echo {input} | tr ' ' '\\n' | head -n 1 || true)

        trap "rm -rf $output $TMPDIR" ERR
        trap "rm -rf $TMPDIR" EXIT

        rm -f $output*

        >&2 echo "Preparing database..."
        sqlite3 -batch $output <<EOF
        CREATE TABLE baf_snvs (
            svid CHARACTER(10),
            sample VARCHAR(100),
            bafs TEXT
        );
        EOF

        >&2 echo "Writing command scripts..."
        cat >$TMPDIR/cmds.txt <<EOF
        .mode tabs
        .import /dev/stdin baf_snvs
        EOF

        >&2 echo "Importing all data..."
        for i in {input}; do
            tail -n +2 $i
        done \
        | sort --compress-program=gzip -s -S 30G --parallel=8 -k 1,1 \
        | time sqlite3 --init $TMPDIR/cmds.txt $output

        >&2 echo "Creating indices..."
        time sqlite3 -batch $output <<EOF
        CREATE INDEX idx_svid ON baf_snvs (svid);
        EOF

        pushd $(dirname $output)
        md5sum $(basename $output) >$(basename $output).md5
        """
        )

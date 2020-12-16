SPLIT_COUNT = 100
SPLIT_PADDING = 10_000

# rule gatk_hc_gvcf_split_genome:
#     output:
#         expand("work/gatk_hc_gvcf_intervals/out/{i:04d}-scattered.interval_list", i=range(SPLIT_COUNT))
#     shell:
#         """
#         set -x

#         export TMPDIR=$(mktemp -d)
#         trap "rm -rf $TMPDIR" EXIT

#         awk -F $'\t' \\
#             'BEGIN {{ OFS=FS }} ($1 != "hs37d5" && $1 != "NC_007605") {{ print $1, 0, $2 }}' \\
#             {REFERENCE}.fai \\
#         > $TMPDIR/genome.bed

#         output_first=$(echo {output} | awk '{{print $1}}')

#         gatk SplitIntervals \\
#             -R {REFERENCE} \\
#             -L $TMPDIR/genome.bed \\
#             --scatter-count {SPLIT_COUNT} \\
#             -O $(dirname $output_first) \\
#             --interval-padding {SPLIT_PADDING}
#         """


# def gatk_hc_gvcf_single_on_interval_input(wildcards):
#     itv_tpl = "work/gatk_hc_gvcf_intervals/out/{:04d}-scattered.interval_list"
#     return {
#         "bam": config["bams"][wildcards.sample],
#         "interval": itv_tpl.format(int(wildcards.interval)),
#     }

# rule gatk_hc_gvcf_single_on_interval:
#     input: unpack(gatk_hc_gvcf_single_on_interval_input)
#     output:
#         multiext("work/gatk_hc_gvcf_single_on_interval/{sample}.{interval}/{sample}.{interval}", ".g.vcf.gz", ".g.vcf.gz.tbi", ".g.vcf.gz.md5",    ".g.vcf.gz.tbi.md5")
#     shell:
#         """
#         set -x

#         output=$(echo {output} | tr ' ' '\\n' | grep 'vcf.gz$')

#         gatk --java-options "-Xmx4g" HaplotypeCaller \\
#             --native-pair-hmm-threads 4 \\
#             -L {input.interval} \\
#             -R {REFERENCE} \\
#             -I {input.bam} \\
#             -O $output \\
#             -ERC GVCF

#         tabix -f $output
#         cd $(dirname $output)
#         md5sum $(basename $output) >$(basename $output).md5
#         md5sum $(basename $output).tbi >$(basename $output).tbi.md5
#         """


# rule gatk_hc_gvcf_all_on_interval:
#     input:
#         gvcfs=expand(
#             "work/gatk_hc_gvcf_single_on_interval/{sample}.{{interval}}/{sample}.{{interval}}.g.vcf.gz",
#             sample=SAMPLES,
#         ),
#     output:
#         gvcf="work/gatk_hc_gvcf_all_on_interval/{interval}/all.{interval}.g.vcf.gz",
#     shell:
#         """
#         set -x

#         gatk --java-options "-Xmx4g" CombineGVCFs \\
#             $(for f in {input.gvcfs}; do echo -V $f; done) \\
#             -R {REFERENCE} \\
#             -O {output.gvcf}

#         tabix -f {output.gvcf}
#         cd $(dirname {output.gvcf})
#         md5sum $(basename {output.gvcf}) >$(basename {output.gvcf}).md5
#         md5sum $(basename {output.gvcf}).tbi >$(basename {output.gvcf}).tbi.md5
#         """

# rule gatk_hc_gvcf_genotype_on_interval:
#     input:
#         gvcf="work/gatk_hc_gvcf_all_on_interval/{interval}/all.{interval}.g.vcf.gz",
#     output:
#         vcf="work/gatk_hc_gvcf_vcf/{interval}/all.{interval}.vcf.gz",
#     shell:
#         """
#         set -x

#         gatk --java-options "-Xmx4g" GenotypeGVCFs \\
#             -V {input.gvcf} \\
#             -R {REFERENCE} \\
#             -O {output.vcf}

#         tabix -f {output.vcf}
#         cd $(dirname {output.vcf})
#         md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
#         md5sum $(basename {output.vcf}).tbi >$(basename {output.vcf}).tbi.md5
#         """

# rule gatk_hc_gvcf_vcf:
#     input:
#         expand(
#             "work/gatk_hc_gvcf_vcf/{i:04d}/all.{i:04d}.vcf.gz",
#             i=range(SPLIT_COUNT)
#         )
#     output:
#         vcf="work/gatk_hc_gvcf_all/vcf/all.vcf.gz",
#     shell:
#         """
#         set -x

#         bcftools concat -O z -o {output.vcf} \\
#             {input}

#         tabix -f {output.vcf}
#         cd $(dirname {output.vcf})
#         md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
#         md5sum $(basename {output.vcf}).tbi >$(basename {output.vcf}).tbi.md5
#         """


# rule gatk_hc_gvcf_vcf_filter_all:
#     input:
#         vcf="work/gatk_hc_gvcf_all/vcf/all.vcf.gz",
#     output:
#         vcf="work/gatk_hc_gvcf_all_filtered/vcf/all_filtered.vcf.gz",
#         bcf="work/gatk_hc_gvcf_all_filtered/vcf/all_filtered.bcf",
#     params:
#         filter_snv='(type == "snp" && QD >= 2.0 && FS <= 6.0 && MQ >= 40.0 && MQRankSum >= -12.5)',
#         filter_indel='(type != "snp" && QD >= 2.0 && FS <= 200.0 && ReadPosRankSum >= -20.0)',
#     shell:
#         """
#         set -x

#         bcftools view -O b -o {output.bcf} {input.vcf} \
#             -i '{params.filter_snv} || {params.filter_indel}'

#         bcftools view -O z -o {output.vcf} {output.bcf}

#         tabix -f {output.bcf}
#         tabix -f {output.vcf}

#         cd $(dirname {output.vcf})
#         md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
#         md5sum $(basename {output.vcf}).tbi >$(basename {output.vcf}).tbi.md5
#         md5sum $(basename {output.bcf}) >$(basename {output.vcf}).csi
#         md5sum $(basename {output.bcf}).csi >$(basename {output.vcf}).csi.md5
#         """


# rule gatk_hc_gvcf_sample:
#     input:
#         vcf="work/gatk_hc_gvcf_all_filtered/vcf/all_filtered.vcf.gz",
#     output:
#         bcf="work/gatk_hc_gvcf_sample_filtered/{sample}/{sample}_filtered.bcf",
#         vcf="work/gatk_hc_gvcf_sample_filtered/{sample}/{sample}_filtered.vcf.gz",
#     shell:
#         """
#         set -x

#         suffix=-N1-DNA1-WGS1

#         bcftools view -s {wildcards.sample}$suffix -O u {input.vcf} \\
#         | bcftools view -i 'GT="alt"' -O b -o {output.bcf}

#         bcftools view -O z -o {output.vcf} {output.bcf}

#         tabix -f {output.bcf}
#         tabix -f {output.vcf}

#         cd $(dirname {output.vcf})
#         md5sum $(basename {output.vcf}) >$(basename {output.vcf}).md5
#         md5sum $(basename {output.vcf}).tbi >$(basename {output.vcf}).tbi.md5
#         md5sum $(basename {output.bcf}) >$(basename {output.bcf}).md5
#         md5sum $(basename {output.bcf}).csi >$(basename {output.bcf}).csi.md5
#         """


# rule gatk_hc_gvcf_bcf:
#     input:
#         vcf="work/gatk_hc_gvcf_all/vcf/all.vcf.gz",
#     output:
#         bcf="work/gatk_hc_gvcf_all/bcf/all.bcf",
#     shell:
#         """
#         set -x

#         bcftools view -O b -o {output.bcf} {input.vcf}
#         tabix -f {output.bcf}
#         cd $(dirname {output.bcf})
#         md5sum $(basename {output.bcf}) >$(basename {output.bcf}).md5
#         md5sum $(basename {output.bcf}).csi >$(basename {output.bcf}).csi.md5
#         """

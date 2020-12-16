@attr.s(frozen=True, auto_attribs=True)
class SampleBatch:
    #: The short name (hash) of the batch.
    name: str
    #: The full name (hash) of the batch.
    full_name: str
    #: The sample names.
    samples: typing.Tuple[str]

    def summary(self) -> str:
        """Return short string summary"""
        return "SampleBatch(name=%s, num_samples=%d)" % (self.name, len(self.samples))

    @staticmethod
    def from_samples(samples: typing.Iterable[str], short_len: int = 16):
        samples = tuple(sorted(samples))
        name_hash = sha256(";".join(samples).encode()).hexdigest()
        return SampleBatch(name=name_hash[:short_len], full_name=name_hash, samples=samples)


def compute_batches(samples: typing.List[str], batch_size: int) -> typing.Dict[str, SampleBatch]:
    samples = list(sorted(samples))
    num_batches = (len(samples) + batch_size - 1) // batch_size
    sample_to_batch = {s: i % num_batches for i, s in enumerate(samples)}
    batched_samples = {i: set() for i in range(num_batches)}
    for s, batch_no in sample_to_batch.items():
        batched_samples[batch_no].add(s)
    batches = [SampleBatch.from_samples(samples) for samples in batched_samples.values()]
    return {b.name: b for b in batches}


#: Batch mapping from samples
BATCHES_CNMOPS = {
    sex: compute_batches(
        samples=[s.name for s in samples], batch_size=config["settings"]["cn.mops"]["batch_size"]
    )
    for sex, samples in BY_SEX.items()
}

# Print out batches.
print("BATCHES\n=======\n", file=sys.stderr)
for sex, batches in BATCHES_CNMOPS.items():
    print("%s\n%s\n" % (sex, "-" * len(sex)), file=sys.stderr)
    print("  count: %d\n" % len(batches), file=sys.stderr)
    for batch in batches.values():
        print("- %s" % batch.summary(), file=sys.stderr)
        for sample in batch.samples:
            print("    - %s" % sample, file=sys.stderr)
    print(file=sys.stderr)


def cnmops_prepare_batch_input(wildcards):
    return [
        "work/maelstrom_bam_collect_doc/{window_length}/{sample}/{sample}.bcf".format(
            window_length=wildcards.window_length, sample=sample
        )
        for sample in BATCHES_CNMOPS[wildcards.sex][wildcards.batch].samples
    ]


rule cnmops_prepare_batch_bcf:
    input:
        cnmops_prepare_batch_input,
    output:
        multiext(
            "work/cnmops_prepare_batch_bcf/{window_length}/{sex}.{batch}/{sex}.{batch}",
            ".bcf",
            ".bcf.md5",
            ".bcf.csi",
            ".bcf.csi.md5",
        ),
    shell:
        """
        set -x

        output_bcf=$(echo {output} | tr ' ' '\\n' | grep 'bcf$')

        bcftools merge -m id -O b -o $output_bcf {input}
        tabix -f $output_bcf

        pushd $(dirname $output_bcf)
        md5sum $(basename $output_bcf) >$(basename $output_bcf).md5
        md5sum $(basename $output_bcf).csi >$(basename $output_bcf).csi.md5
        """


rule cnmops_prepare_batch_tsv:
    input:
        bcf="work/cnmops_prepare_batch_bcf/{window_length}/{sex}.{batch}/{sex}.{batch}.bcf",
    output:
        multiext(
            "work/cnmops_prepare_batch_tsv/{window_length}/{sex}.{batch}/{sex}.{batch}",
            ".tsv",
            ".tsv.md5",
        ),
    shell:
        """
        set -x

        output_tsv=$(echo {output} | tr ' ' '\\n' | grep 'tsv$')

        echo -en "CHROM\\tSTART\\tEND\\t" >$output_tsv
        bcftools view -h {input.bcf} | grep '^#CHROM' | cut -f 10- >>$output_tsv
        bcftools query -f "%CHROM\\t%POS\\t%INFO/END[\\t%RCV]\\n" {input.bcf} \\
        | egrep "$(if [[ "{wildcards.sex}" == "male" ]]; then echo '^CHROM|^[1-9M]|^X|^Y'; else echo '^CHROM|^[1-9M]|^X'; fi)" \\
        | sed -e 's/\\([0-9]\\+\\)\\.[0-9]*/\\1/g' \\
        >>$output_tsv

        pushd $(dirname $output_tsv)
        md5sum $(basename $output_tsv) >$(basename $output_tsv).md5
        """


rule cnmops_call:
    input:
        tsv="work/cnmops_prepare_batch_tsv/{window_length}/{sex}.{batch}/{sex}.{batch}.tsv",
    output:
        segmentation=(
            "work/cnmops_call/{window_length}/{sex}.{batch}/{sex}.{batch}.segmentation.tsv"
        ),
        cnvs="work/cnmops_call/{window_length}/{sex}.{batch}/{sex}.{batch}.cnvs.tsv",
        cnv_regions="work/cnmops_call/{window_length}/{sex}.{batch}/{sex}.{batch}.cnv_regions.tsv",
    threads: config.get("settings", {}).get("cn.mops", {}).get("threads", 1)
    script:
        "scripts/cnmops_call.R"


def cnmops_to_vcf_input(wildcards):
    for sex, batches in BATCHES_CNMOPS.items():
        for batch in batches.values():
            if wildcards.sample in batch.samples:
                tpl = "work/cnmops_call/{window_length}/{sex}.{batch}/{sex}.{batch}.cnvs.tsv"
                return tpl.format(batch=batch.name, window_length=wildcards.window_length, sex=sex)


rule cnmops_to_vcf:
    input:
        cnmops_to_vcf_input,
    output:
        vcf="work/cnmops_to_vcf/{window_length}/{sample}/{sample}.vcf.gz",
        vcf_md5="work/cnmops_to_vcf/{window_length}/{sample}/{sample}.vcf.gz.md5",
        tbi="work/cnmops_to_vcf/{window_length}/{sample}/{sample}.vcf.gz.tbi",
        tbi_md5="work/cnmops_to_vcf/{window_length}/{sample}/{sample}.vcf.gz.tbi.md5",
    script:
        "scripts/cnmops_cnvs_to_vcf.py"

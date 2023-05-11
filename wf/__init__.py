"""
Minimal template workflow to show the structure of a Latch workflow

For a more comprehensive example, see the assemble_and_sort workflow
For examples on how to use the Latch SDK, see https://docs.latch.bio/examples/workflows_examples.html
"""
import glob
import subprocess
from pathlib import Path

from latch import small_task, workflow
from latch.resources.launch_plan import LaunchPlan
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchFile,
    LatchMetadata,
    LatchRule,
    LatchOutputDir,
    LatchParameter,
)
from latch.verified import trim_galore
from latch.verified.trim_galore import AdapterSequence, BaseQualityEncoding

"""Minimal metadata object - fill in fields with your own values"""
metadata = LatchMetadata(
    display_name="Example WGS Workflow",
    documentation="CHANGE ME",
    author=LatchAuthor(
        name="CHANGE ME",
        email="CHANGE ME",
        github="CHANGE ME",
    ),
    repository="CHANGE ME",
    license="CHANGE ME",
    parameters={
        "input_forward": LatchParameter(
            display_name="Read 1",
            description="The paired-end sequence forward read (often has _R1 in name)",
            rules=[
                LatchRule(
                    regex="(.fastq|.fastq.gz|.fq|.fq.gz)$",
                    message="Only fastq, fastq.gz, fq, fq.gz extensions are valid",
                )
            ],
            batch_table_column=True,  # Show this parameter in batched mode.
        ),
        "input_reverse": LatchParameter(
            display_name="Read 2",
            description="The paired-end sequence reverse read (often has _R2 in name)",
            batch_table_column=True,  # Show this parameter in batched mode.
            rules=[
                LatchRule(
                    regex="(.fastq|.fastq.gz|.fq|.fq.gz)$",
                    message="Only fastq, fastq.gz, fq, fq.gz extensions are valid",
                )
            ],
        ),
        "ref_genome": LatchParameter(
            display_name="Reference Genome File",
            description="The reference genome file for reads to be aligned to",
            batch_table_column=True,
            rules=[
                LatchRule(
                    regex="(.fasta)$",
                    message="Only .fasta extension is valid",
                )
            ],
        ),
        "output_dir": LatchParameter(
            display_name="Output Directory",
            batch_table_column=True
        ),
    },
    tags=[],
)


@small_task
def build_index(
    output_dir: LatchOutputDir,
    ref_genome: LatchFile = LatchFile("latch:///wgs/ref_genome/ecoli_rel606.fasta"),
) -> LatchDir:
    _bwa_cmd = ["bwa", "index", ref_genome.local_path]

    subprocess.run(_bwa_cmd, check=True)

    output = Path(ref_genome.local_path).parent.resolve()

    return LatchDir(str(output), f"{output_dir.remote_path}/ref_genome")


@small_task
def align_reads(
    reads: LatchDir,
    output_dir: LatchOutputDir,
    ref_genome_dir: LatchDir = LatchDir("latch:///wgs/ref_genome"),
) -> LatchFile:

    local_reads_directory = reads.local_path

    # Search for files ending with "_val_2.fq" or "_val_1.fq" in the directory
    file_pattern = local_reads_directory + "/*_val_[1-2].fq"
    matching_files = glob.glob(file_pattern)

    matching_files.sort()

    # Make sure exactly two matching files are found
    if len(matching_files) != 2:
        print("Error: Expected two matching files, but found", len(matching_files))
    else:
        # Assign the absolute filepaths to variables read1 and read2
        read1 = matching_files[0]
        read2 = matching_files[1]

        print("Filepath of read1:", read1)
        print("Filepath of read2:", read2)

    local_ref_dir = ref_genome_dir.local_path
    fastas = glob.glob(f"{local_ref_dir}/*.fasta")
    sam_file = Path("aligned.sam").resolve()

    print(fastas, flush=True)

    if len(fastas) > 0:
        ref_genome = str(Path(fastas[0]).resolve())
        cmd = [
            "bwa",
            "mem",
            ref_genome,
            read1,
            read2,
            "-o",
            str(sam_file),
        ]
        subprocess.run(cmd, check=True)
    else:
        print(
            "No *.fasta file found in directory. Please ensure that your reference genome directory contains a *.fasta file."
        )

    return LatchFile(str(sam_file), f"{output_dir.remote_path}/aligned.sam")


@small_task
def convert_to_bam(
    output_dir: LatchOutputDir,
    sam: LatchFile = LatchFile("latch:///wgs/results/aligned.sam"),
) -> LatchFile:
    bam_file = Path("aligned.bam").resolve()
    _samtools_cmd = [
        "samtools",
        "view",
        "-S",
        "-b",
        sam.local_path,
        "-o",
        str(bam_file),
    ]

    subprocess.run(_samtools_cmd, check=True)

    return LatchFile(str(bam_file), f"{output_dir.remote_path}/aligned.bam")


@small_task
def sort_bam(output_dir: LatchOutputDir, bam: LatchFile) -> LatchFile:
    sorted_bam = Path("aligned.sorted.bam").resolve()
    _sort_cmd = ["samtools", "sort", "-o", str(sorted_bam), bam.local_path]

    subprocess.run(_sort_cmd, check=True)

    return LatchFile(str(sorted_bam), f"{output_dir.remote_path}/aligned.sorted.bam")


@small_task
def variant_calling(
    ref_genome: LatchFile,
    sorted_bam: LatchFile,
    output_dir: LatchOutputDir,
) -> LatchFile:

    # Calculate read coverage
    bcf = Path("raw.bcf").resolve()
    _read_coverage = [
        "bcftools",
        "mpileup",
        "-O",
        "b",
        "-o",
        str(bcf),
        "-f",
        ref_genome.local_path,
        sorted_bam.local_path,
    ]
    subprocess.run(_read_coverage, check=True)

    # Detect SNVs
    vcf = Path("variants.vcf").resolve()
    _snv_detection = [
        "bcftools",
        "call",
        "--ploidy",
        "1",
        "-m",
        "-v",
        "-o",
        str(vcf),
        str(bcf),
    ]
    subprocess.run(_snv_detection, check=True)

    # Filter and report the SNV variants in variant calling format (VCF)
    final_vcf = Path("final_variants.vcf").resolve()
    f = open(str(final_vcf), "w")
    _vcfutils_cmd = [
        "vcfutils.pl",
        "varFilter",
        str(vcf),
    ]
    subprocess.run(_vcfutils_cmd, stdout=f, check=True)

    return LatchFile(str(final_vcf), f"{output_dir.remote_path}/final_variants.vcf")


# change the name of this function to something more descriptive
@workflow(metadata)
def wgs_wf(
    input_forward: LatchFile,
    input_reverse: LatchFile,
    output_dir: LatchOutputDir,
    ref_genome: LatchFile,
) -> LatchFile:
    """
    An example WGS workflow
    """
    ref_genome_dir = build_index(ref_genome=ref_genome, output_dir=output_dir)
    trim_galore_outputs = trim_galore(
        input_forward=input_forward,
        input_reverse=input_reverse,
        output_directory=output_dir,
        base_out=None,
        fastqc_args=None,
        adapter=None,
        adapter2=None,
        consider_already_trimmed=None,
        max_length=None,
        max_n=None,
        clip_R1=None,
        clip_R2=None,
        three_prime_clip_R1=None,
        three_prime_clip_R2=None,
        hardtrim5=None,
        hardtrim3=None,
        quality=20,
        base_quality_encoding=BaseQualityEncoding.phred33,
        fastqc=True,
        adapter_sequence=AdapterSequence.auto,
        stringency=1,
        error_rate=0.01,
        gzip_output_files=False,
        length=20,
        trim_n=False,
        report_file=True,
        polyA=False,
        implicon=False,
        retain_unpaired=True,
        length_1=35,
        length_2=35,
    )
    sam = align_reads(
        ref_genome_dir=ref_genome_dir, reads=trim_galore_outputs, output_dir=output_dir
    )
    bam = convert_to_bam(sam=sam, output_dir=output_dir)
    sorted_bam = sort_bam(bam=bam, output_dir=output_dir)
    return variant_calling(
        ref_genome=ref_genome, sorted_bam=sorted_bam, output_dir=output_dir
    )


"""
Add test data with a LaunchPlan. Provide default values in a dictionary with
the parameter names as the keys. These default values will be available under
the 'Test Data' dropdown at console.latch.bio.
"""
LaunchPlan(
    wgs_wf,
    "Test Data",
    {
        "input_forward": LatchFile(
            "latch:///Bn510_local/raw/Bn0510A-01-x-Bg121-b01-t08-sWD1422_S1_L001_R1_001.fastq.gz"
        ),
        "input_reverse": LatchFile(
            "latch:///Bn510_local/raw/Bn0510A-01-x-Bg121-b01-t08-sWD1422_S1_L001_R2_001.fastq.gz"
        ),
        "ref_genome": LatchFile(
            "latch:///Bn510_local/references/fasta+gtf/sWD1422_20240414-EY_BVU_2094.fasta"
        ),
        "output_dir": LatchDir("latch:///Bn510_local/test_wgs_output"),
    },
)

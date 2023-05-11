from wf import (
    wgs_wf,
    build_index,
    align_reads,
    convert_to_bam,
    sort_bam,
    variant_calling,
)
from latch.types import LatchFile, LatchDir

# build_index(
#     ref_genome=LatchFile(
#         "latch:///Bn510_local/references/fasta+gtf/sWD1422_20240414-EY_BVU_2094.fasta"
#     ),
#     output_dir=LatchDir("latch:///Test WGS Outputs"),
# )

# align_reads(
#     reads=LatchDir("latch:///Trim Galore Output/Bn0510A-01-x-Bg121-b01-t08-sWD1422"),
#     output_dir=LatchDir("latch:///Test WGS Outputs"),
#     ref_genome_dir=LatchDir("latch:///Test WGS Outputs/ref_genome"),
# )


# convert_to_bam(
#     sam=LatchFile("latch:///Test WGS Outputs/aligned.sam"),
#     output_dir=LatchDir("latch:///Test WGS Outputs"),
# )

# sort_bam(
#     bam=LatchFile("latch:///Test WGS Outputs/aligned.bam"),
#     output_dir=LatchDir("latch:///Test WGS Outputs")
# )

# variant_calling(
#     ref_genome=LatchFile("latch:///Bn510_local/references/fasta+gtf/sWD1422_20240414-EY_BVU_2094.fasta"),
#     sorted_bam=LatchFile("latch:///Test WGS Outputs/aligned.sorted.bam"),
#     output_dir=LatchDir("latch:///Test WGS Outputs")
# )

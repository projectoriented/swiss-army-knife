import pandas as pd
from itertools import product
import re


# --- Config constants --- #
configfile: "config.yaml"


manifest = config.get("manifest", "manifest.tab")
ref = config.get("ref")
MIN_BASE_Q = config.get("min_base_q")


manifest_df = pd.read_csv(manifest, sep="\t", index_col="sample")

containers = {
    "gatk4": config["img"]["gatk4"],
}


# Regular functions
def amnt_flank(len_diff: int) -> int:
    if len_diff < 100:
        return 20
    elif 100 <= len_diff < 200:
        return 25
    elif 200 <= len_diff < 500:
        return 40
    elif 500 <= len_diff < 1000:
        return 100
    elif 1000 <= len_diff < 2000:
        return 100
    elif 2000 <= len_diff < 4000:
        return 250
    elif len_diff >= 4000:
        return 300


# Input functions
def get_aln(wildcards):
    return manifest_df.at[wildcards.sample, "aln"]


# def get_interval(wildcards):
#     # right now we assume just one region
#     # def inner(wildcards):
#     #     region_str = str(manifest_df.at[wildcards.sample, 'region'])
#     #     region_df = pd.MultiIndex.from_tuples(product([wildcards.sample],region_str.split(',')), names=('sample','region'))
#     region_str = str(manifest_df.at[wildcards.sample, "region"])
#
#     def parse_region(string):
#         pattern = re.compile(
#             r"(?P<contig>.+)-(?P<spos>\d+)-(?P<svtype>\S+)-(?P<svlen>\d+)"
#         )
#         match = pattern.match(string).groupdict()
#         end_pos = match["spos"] + match["svlen"]
#         return match["contig"] + ":" + match["spos"] + "-" + end_pos
#
#     return region_str if region_str.find(":") else parse_region(region_str)


def final_output(wildcards):
    pass


wildcard_constraints:
    sample="|".join(manifest_df.index),


rule gatk4_depth:
    input:
        fa=ref,
        bam=get_aln,
    output:
        sample_depth=temp(
            "{sample}/gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}.csv"
        ),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    container:
        containers["gatk4"]
    params:
        interval=lambda wildcards: f"{wildcards.contig}:{wildcards.pos}-{wildcards.end}",
        flank=lambda wildcards: amnt_flank(int(wildcards.flank)),
    shell:
        """
        gatk DepthOfCoverage \
            --omit-locus-table --omit-interval-statistics --omit-per-sample-statistics \
            --print-base-counts \
            --min-base-quality {MIN_BASE_Q} \
            --interval-padding {params.flank} -L {params.interval} \
            -R {input.fa} -I {input.bam} -O {output.sample_depth}
        """


rule transform_output:
    input:
        sample_depth=rules.gatk4_depth.output.sample_depth,
    output:
        region_stats=temp(
            "{sample}/gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}.tsv"
        ),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    params:
        svlen="",
    run:
        df = pd.read_csv(input.sample_depth, header=0)

        flank = int(wildcards.flank)
        total = int(df.shape[0])
        end_pos = int(wildcards.end)
        target_column = "Total_Depth"
        interval = (
            f"{wildcards.contig}:{wildcards.pos}-{end_pos}"
            if not params.svlen
            else params.svlen
        )

        # Get the depth (mean + min) for both left and right flanking regions
        left_side = (
            df.loc[0:flank, target_column]
            .agg({"L_MEAN": "mean", "L_MIN": "min"})
            .to_dict()
        )
        right_side = (
            df.loc[flank:total, target_column]
            .agg({"R_MEAN": "mean", "R_MIN": "min"})
            .to_dict()
        )

        # Get the mean depth for target region
        target_region = (
            df.loc[flank:end_pos, target_column]
            .agg({"MEAN": "mean", "MIN": "min"})
            .to_dict()
        )

        # Clean up
        del df

        # Combine the dictionaries
        final_dict = {
                "SAMPLE": wildcards.sample,
                "REGION": interval,
                **left_side,
                **target_region,
                **right_side,
            }

        # Write out
        pd.DataFrame.from_dict([final_dict]).to_csv(
            output.region_stats, index=False, sep="\t"
        )


rule merge:
    input:
        region_stats=expand(rules.transform_output.output.region_stats),
    output:
        merged_stats="depth_stats.tsv.gz",
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    run:
        df = pd.concat(
            [pd.read_csv(item, header=0, sep="\t") for item in input.region_stats]
        )

        # Implement parent stuff here in the future
        # Write out
        df.to_csv(output.merged_stats, index=False, header=True, sep="\t")

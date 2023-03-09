import pandas as pd
from itertools import product, chain
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


# --- Regular Functions --- #
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


def construct_file_names(l: list) -> list:
    def parse_region(string):
        pattern_one = re.compile(r"(?P<contig>.+):(?P<pos>\d+)-(?P<end>\d+)")
        pattern_two = re.compile(r"(?P<contig>.+)-(?P<pos>\d+)-(?P<svtype>\S+)-(?P<svlen>\d+)")
        if string.find(":") >= 0:
            match = pattern_one.match(string).groupdict()
        else:
            match = pattern_two.match(string).groupdict()
            match["end"] = int(match["pos"]) + int(match["svlen"])
        return match

    file_names = []
    for entry in l:
        d = parse_region(entry)

        len_diff = int(d["end"]) - int(d["pos"])
        d['flank'] = amnt_flank(len_diff=len_diff)

        formatted_fn = "stats_{contig}_{pos}-{end}_F-{flank}".format(**d)
        file_names.append(formatted_fn)
    return file_names


# --- Input functions --- #
def get_aln(wildcards):
    return manifest_df.at[wildcards.sample, "aln"]

def collect_all_files(_) -> list:
    sample_list = manifest_df.index
    def get_formatted(s):
        region_list = manifest_df.at[s, 'region'].split(',')
        nested_list = [['gatk4_depth'],[s], construct_file_names(region_list), ['tsv']]
        return ['{}/{}_{}.{}'.format(*x) for x in product(*nested_list)]
    return list(chain.from_iterable(map(get_formatted, sample_list)))


wildcard_constraints:
    sample="|".join(manifest_df.index),


rule gatk4_depth:
    input:
        fa=ref,
        bam=get_aln,
    output:
        sample_depth=temp(
            "gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}.csv"
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
            "gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}.tsv"
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
        region_stats=collect_all_files,
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

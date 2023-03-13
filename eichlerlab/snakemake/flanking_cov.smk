import pandas as pd
from itertools import product, chain
import re


# --- Config constants --- #
configfile: "config.yaml"


manifest = config.get("manifest", "manifest.tab")
ref = config.get("ref")
MIN_BASE_Q = config.get("min_base_q", 40)


manifest_df = pd.read_csv(manifest, sep="\t", index_col="sample")

containers = {
    "gatk4": config["img"]["gatk4"],
}


wildcard_constraints:
    sample="|".join(manifest_df.index),


# --- Regular Functions --- #
def amnt_flank(len_diff: int) -> int:
    """
    This function is tuned for structural variants and taken from ebert 2021 science paper.
    """
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
    def parse_region(string: str) -> dict:
        pattern_one = re.compile(r"(?P<contig>.+):(?P<pos>\d+)-(?P<end>\d+)")
        pattern_two = re.compile(r"(?P<contig>.+)-(?P<pos>\d+)-(?P<svtype>\S+)-(?P<svlen>\d+)")
        if string.find(":") >= 0:
            match = pattern_one.match(string).groupdict()
            match['alt_id'] = 'NA'
        else:
            match = pattern_two.match(string).groupdict()
            match["end"] = int(match["pos"]) + int(match["svlen"])
            match['alt_id'] = string
        return match
    file_names = []
    for entry in l:
        d = parse_region(entry)

        len_diff = int(d["end"]) - int(d["pos"])
        d['flank'] = amnt_flank(len_diff=len_diff)

        formatted_fn = "stats_{contig}_{pos}-{end}_F-{flank}_{alt_id}".format(**d)
        file_names.append(formatted_fn)
    return file_names


# --- Input functions --- #
def get_aln(wildcards):
    return manifest_df.at[wildcards.sample, "aln"]


def collect_all_files(_) -> list:
    sample_list = manifest_df.index
    def get_formatted(s):
        try:
            region_list = manifest_df.at[s, 'region'].split(',')
        except AttributeError:
            region_list = manifest_df.at[s, 'region']
        nested_list = [['gatk4_depth'],[s], construct_file_names(region_list), ['tsv']]
        return ['{}/{}_{}.{}'.format(*x) for x in product(*nested_list)]
    return list(chain.from_iterable(map(get_formatted, sample_list)))


def get_interval(wildcards):
    return f"{wildcards.contig}:{wildcards.pos}-{wildcards.end}"


# --- Target rule --- #
rule all:
    input:
        'depth_stats.tsv.gz'


# --- Action rules --- #
rule gatk4_depth:
    input:
        fa=ref,
        bam=get_aln,
    output:
        sample_depth=temp(
            "gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}_{alt_id}.csv"
        ),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    container:
        containers["gatk4"]
    params:
        interval=get_interval,
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
            "gatk4_depth/{sample}_stats_{contig}_{pos}-{end}_F-{flank}_{alt_id}.tsv"
        ),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: attempt * 2,
        hrs=72,
    params:
        native_id=lambda wildcards: get_interval if wildcards.alt_id == 'NA' else wildcards.alt_id
    run:
        df = pd.read_csv(input.sample_depth, header=0)

        # Definitions
        flank = int(wildcards.flank)
        total = int(df.shape[0])
        len_diff = int(wildcards.end) - int(wildcards.pos)

        target_column = "Total_Depth"
        id = params.native_id
        order_cols = ':'.join(['L_MEAN', 'L_MIN', 'SD', 'MEAN', 'R_MIN', 'R_MEAN'])

        # Get the depth (mean + min) for both left and right flanking regions
        left_side = (
            df.loc[0:(flank - 1), target_column]
            .agg(['mean', 'min'])
            .astype(int)
            .to_list()
        )
        right_side = (
            df.loc[(flank + 1):total, target_column]
            .agg(['min', 'mean'])
            .astype(int)
            .to_list()
        )
        # Get the depth (sd + mean) for target region
        target_region = (
            df.loc[flank:(len_diff + flank), target_column]
            .agg(['std', 'mean'])
            .astype(int)
            .to_list()
        )

        # Clean up
        del df

        # Combine
        final_list = left_side + target_region + right_side
        final_dict = {
            "REGION": id,
            "INFO": order_cols,
            wildcards.sample: ':'.join(map(str,final_list)),
            "FLANK (bp)": wildcards.flank
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
            [pd.read_csv(item, header=0, sep="\t", index_col=['REGION', 'INFO', 'FLANK (bp)']) for item in input.region_stats]
            ,axis=1
        )

        # Merge duplicated columns into one
        df = df.groupby(df.columns,axis=1).first()

        # Implement parent stuff here in the future
        # Write out
        df.to_csv(output.merged_stats, header=True, sep="\t")

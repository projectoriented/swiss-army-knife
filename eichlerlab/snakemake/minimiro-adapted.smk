import pandas as pd

configfile: 'config.yaml'

manifest = config.get("manifest", "manifest.tab") # manifest: sample,aln,target-regions.bed3,
manifest_df = pd.read_csv(manifest, sep='\t', index_col='sample')

# steps
# 1. get the bam/cram/sam file and convert to paf
# 2. subset paf to target regions using rustybam
# 3. break subset.paf at the size of the SVLEN to create miropeat intervals
# 4. repeatmasker for ref and query
# 5. dupmasker for ref and query
# 6. gene track for ref and query

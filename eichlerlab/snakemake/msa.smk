import pandas as pd


configfile: 'config.yaml'
manifest = config.get('manifest', 'manifest.tab')
ref=config.get('ref')
bed=config.get('regions', '')

manifest_df = pd.read_csv(manifest, sep='\t', index_col='sample')

# I want to make this pipeline to use http://yulab-smu.top/ggmsa/articles/ggmsa.html

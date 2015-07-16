import pandas as pd

chr1 = pd.read_table('dbscSNV1.1.chr1', sep = '\t', na_values = '.')


#rf_score histogram by region
chr1.rf_score.hist(by=chr1.Ensembl_region, bins=20)

#ada_score histogram by region
chr1.ada_score.hist(by=chr1.Ensembl_region, bins=20)

#Number and type of vars by region
chr1['VAR'] = chr1['ref'] + chr1['alt']
chr1.VAR.value_counts().plot(by=chr1.Ensembl_region)

#Number and type of vars by region
region = chr1.groupby(chr1.Ensembl_region)
region.VAR.value_counts().plot(kind='bar')

#number of SNVs by region
chr1.RefSeq_region.value_counts().plot(kind='barh')

#correlation between ada_score and rf_score
chr1.plot(x='ada_score',y='rf_score', kind='scatter')

#ada_score, rf_score scatter for the splicing region
splicing = region.get_group('splicing')
splicing.plot(x='ada_score',y='rf_score',kind='scatter',title='splicing')





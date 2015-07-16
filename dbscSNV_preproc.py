import pandas as pd
import numpy as np


def clean_import_scsnv():
    """read each chr 1-22 and X/Y from the dbscSNV download into a dict of dataframes 
    for further processing"""
    
    chrom_db = {}
    cols = [0,1,2,3,8,16,17]
    col_names = ['chr', 'hg19_pos', 'ref', 'alt', 'RefSeq_region', 'ada_score', 'rf_score']

    for i in range(23):
        if i > 0:
            chrom_db[str(i)]  = pd.read_table('dbscSNV1.1.chr'+str(i), sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=500)
                
    chrom_db.setdefault('X', pd.read_table('dbscSNV1.1.chrX', sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=500))
                
    chrom_db.setdefault('Y', pd.read_table('dbscSNV1.1.chrY', sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=500))
    
    return chrom_db  
    
    
def import_exome_capture_bed():
    """import exome regions into a dataframe. modify -/+ 200 bps around region start
    and end.  export the dataframe"""
    
    col_names = ['chr','hg19_start','hg19_end']
    cols = [0,1,2]
    
    exomes = pd.read_table('../Agilent_SureSelect_ExomeV5_Covered.bed.txt', names=col_names,
                usecols = cols, nrows=100)
    
    #pad 200 bp from start and end of regions to account for PCR variability            
    exomes['hg19_start'] = exomes['hg19_start'] - 200
    exomes['hg19_end'] = exomes['hg19_end'] + 200
    
    exomes.chr = exomes.chr.str.strip('chr')   #change to chrom numbers only (as strings) for downstream compatibility
    return exomes

    
    
def comp_ada_rf_scores(chrom_db):
    """compare the ada_score and rf_score to a threshold. Use a logical OR to indicate
    if value is above threshold for one or both scores. Add new column to database"""
    
    ada_score_true = chrom_db['ada_score'] >= 0.8     #these values are chosen somewhat arbitrarily and 
    rf_score_true = chrom_db['rf_score'] >= 0.65      #can change to reflect our experience with the variant calling workflow
    
    splice_alt = []
    
    for i in range(len(rf_score_true)):
        if rf_score_true[i] or ada_score_true[i]:
            splice_alt.append(1)
        else: splice_alt.append(0)
        
    chrom_db['splice_alt'] = splice_alt
    return chrom_db
    

def filter_scsnv_by_exomes(chrom_db, exome_db):
    """filter each chr dataframe by the exome capture regions in the bed file. return
    a dict of the new dataframes with only snvs within exome regions"""
    
    exome_regions = []
    
    for i in range(len(exome_db.hg19_start)):  
        for j in range(exome_db.hg19_start[i], exome_db.hg19_end[i]):
            exome_regions.append(j)

    chrom_db[chrom_db.hg19_pos.isin(exome_regions)]
    
    return chrom_db
    
    
def write_flat_file_dbscSNV():
    """write out a complete flat file in VCF format of the filtered and processed 
    dbscSNV database at exome capture regions only"""
    pass
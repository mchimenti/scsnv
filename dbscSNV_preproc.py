"""Import the dbscSNV database of known SNVs with their associated ensemble splice-altering scores.
Clean the database, through away regions not in exome workflow. Calculate yes/no score from both 
ensemble scores.  Write out flat file at the end."""

import pandas as pd


def clean_import_scsnv():
    """read each chr 1-22 and X/Y from the dbscSNV download into a dict of dataframes 
    for further processing"""
    
    chrom_dict = {}
    cols = [0,1,2,3,8,16,17]
    col_names = ['chr', 'hg19_pos', 'ref', 'alt', 'RefSeq_region', 'ada_score', 'rf_score']

    #for i in range(2):
    #    if i > 0:
    #        chrom_dict[str(i)]  = pd.read_table('dbscSNV1.1.chr'+str(i), sep = '\t', 
    #            na_values = '.', usecols=cols, names=col_names, header=0)
    #            
    #chrom_dict.setdefault('X', pd.read_table('dbscSNV1.1.chrX', sep = '\t', 
    #            na_values = '.', usecols=cols, names=col_names, header=0))
    #            
    #chrom_dict.setdefault('Y', pd.read_table('dbscSNV1.1.chrY', sep = '\t', 
    #            na_values = '.', usecols=cols, names=col_names, header=0))
                
    for i in range(23):
        if i > 0:
            chrom_dict[str(i)]  = pd.read_table('dbscSNV1.1.chr'+str(i), sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=10000)
                
    chrom_dict.setdefault('X', pd.read_table('dbscSNV1.1.chrX', sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=10000))
                
    chrom_dict.setdefault('Y', pd.read_table('dbscSNV1.1.chrY', sep = '\t', 
                na_values = '.', usecols=cols, names=col_names, header=0, nrows=10000))
    
    return chrom_dict  
    
    
def import_exome_capture_bed():
    """import exome regions into a dataframe. modify -/+ 200 bps around region start
    and end.  export the dataframe"""
    
    col_names = ['chr','hg19_start','hg19_end']
    cols = [0,1,2]
    exomes_dict = {}
    
    exomes = pd.read_table('../Agilent_SureSelect_ExomeV5_Covered.bed.txt', names=col_names,
                usecols = cols)
    
    #pad 200 bp from start and end of regions to account for PCR variability            
    exomes['hg19_start'] = exomes['hg19_start'] - 200
    exomes['hg19_end'] = exomes['hg19_end'] + 200
    
    exomes.chr = exomes.chr.str.strip('chr')   #change to chrom numbers only (as strings) for downstream compatibility
    exomes_grouped = exomes.groupby('chr')
    
    for i,j in exomes_grouped:
        exomes_dict.setdefault(i,j)    #populate the dict with key:value -> chrom:dataframe
    
    for i in exomes_dict:
        exomes_dict[i] = exomes_dict[str(i)].reset_index(drop=True)   #need to reset the index or create_exome_list() won't work
    
    return exomes_dict

        
def create_exome_list(chrom):
    """enumerate the regions contained within the exome capture kit in order
    to take advantage of pandas 'isin' method on dataframes"""
    
    exome_regions = []
    
    for i in range(len(chrom.hg19_start)):  
        for j in range(chrom.hg19_start[i], chrom.hg19_end[i]):
            exome_regions.append(j)
            
    return exome_regions


def filter_scsnv_by_exome_list(chrom_db, bed_chrom_dict):
    """filter each chr dataframe by the exome capture regions in the bed file. return
    the new dataframes with only snvs within exome regions"""

    chrom_db_filtered = chrom_db[chrom_db.hg19_pos.isin(create_exome_list(bed_chrom_dict))]
    chrom_db_filtered = chrom_db_filtered.reset_index(drop=True)
    
    return chrom_db_filtered
    

def comp_ada_rf_scores(chrom_db):
    """compare the ada_score and rf_score to a threshold. Use a logical OR to indicate
    if value is above threshold for one or both scores. Add new column to database"""
    
    ada_score_true = chrom_db['ada_score'] >= 0.6     #these values are chosen somewhat arbitrarily and 
    rf_score_true = chrom_db['rf_score'] >= 0.5      #can change to reflect our experience with the variant calling workflow
    
    splice_alt = []
    
    for i in range(len(rf_score_true)):
        if rf_score_true[i] or ada_score_true[i]:
            splice_alt.append(1)
        else: splice_alt.append(0)
        
    chrom_db['splice_alt'] = splice_alt
    
    chrom_db = chrom_db[chrom_db.splice_alt == 1]  #drop rows that calculate to zero
    chrom_db = chrom_db.reset_index(drop=True)     #reset index again after dropping rows
    
    return chrom_db
    
def format_flat_file_dbscSNV(chrom_preformat):
    """formatting manipulations of the dataframe to get into form for merge into VCF file
    final format: 'chr#  startpos  endpos  info_field' """
    
    chrom_final = comp_ada_rf_scores(chrom_preformat)  #add the 'splice_alt' summary score column
    
    chrom_final["hg19_pos_end"] = chrom_final["hg19_pos"]
    
    chrom_final["info"] = chrom_final.apply(lambda r: (r["ref"] + "->"+ r["alt"]), axis=1) #create info field with ref>alt

    chrom_final["chr"] = "chr" + str(chrom_final.chr[1])
    
    chrom_final.drop(["ref","alt","RefSeq_region","ada_score","rf_score","splice_alt"], axis=1, inplace=True)
    
    return chrom_final
    
def write_flat_file_dbscSNV(chrom_final_format, chrom_num):
    """write out a complete flat file in VCF format of the filtered and processed 
    dbscSNV database at exome capture regions only"""
    print chrom_final_format.head(n=25)
    #chrom_final_format.to_csv("dbscSNV_proc_"+str(chrom_num)+".txt", sep='\t', index=False, header=False)
    

def main():
    snv_dict = clean_import_scsnv()    #dict with key:chromosome, value:dataframe of SNVs
    exomes_dict = import_exome_capture_bed()  #dict with key:chromosome, value:dataframe of BED regions for EXome capture
    assert snv_dict.keys() == exomes_dict.keys()  #keys should match: chrom 1-22 and X and Y
    
    snv_list = []   #list to collect the filtered chromosome SNV dataframes
    
    for i in snv_dict:
        snv_list.append(filter_scsnv_by_exome_list(snv_dict[str(i)], exomes_dict[str(i)]))
        
    for i in range(len(snv_list)):
        write_flat_file_dbscSNV(format_flat_file_dbscSNV(snv_list[i]), snv_list[i].chr[1])

if __name__ == '__main__':
    
    main()
    #pass    
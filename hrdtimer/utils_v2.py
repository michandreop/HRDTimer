# -------------------- Imports -------------------------
import re
import pandas as pd
import vcfpy
import logging

# --------- VCF parsing/handling -----------
def read_vcf(record):
    return {
        "Chromosome": record.CHROM,
        "Position": record.POS,
        "ID": record.ID,
        "Reference Allele": record.REF,
        "Alternate Alleles": record.ALT,
        "Quality Score": record.QUAL,
        "Info": record.INFO
    }

def dict_to_vcf_INFO(d):
    parts = []
    for k, v in d.items():
        if isinstance(v, list):
            v_str = ','.join(str(i) for i in v)
        else:
            v_str = str(v)
        parts.append(f"{k}:{v_str}")
    return ';'.join(parts)

def INFO_field_to_dict(s):
    d = {}
    for item in s.split(';'):
        if not item:
            continue
        key, val = item.split(':', 1)
        if ',' in val:
            val_list = val.split(',')
            # Try to infer types (int, float)
            try:
                val_list = [x for x in val_list]
            except:
                pass
            d[key] = val_list
        else:
            # Try to infer scalar types
            if val == 'TRUE':
                d[key] = True
            elif val == 'FALSE':
                d[key] = False
            else:
                try:
                    d[key] = int(val)
                except ValueError:
                    try:
                        d[key] = float(val)
                    except ValueError:
                        d[key] = val
    return d

def vcf_to_dataframe(vcf_file):
    data = []
    with open(vcf_file, 'r') as f:
        reader = vcfpy.Reader.from_stream(f)
        for record in reader:
            data.append(read_vcf(record))
    return pd.DataFrame(data)

def dataframe_to_vcf(df, output_vcf):
    cols = df.columns[:df.columns.get_loc('INFO') + 1]
    df = df[cols].copy()
    df['INFO'] = df['INFO'].apply(dict_to_vcf_INFO)
    with open(output_vcf, 'w') as vcf:
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tINFO\n")
    df.to_csv(output_vcf, sep="\t", mode='a', index=False, header=False)

def extract_clonality(value):
    v = str(value)
    
    if v.startswith("clonal"):
        # The regex looks for content inside []
        match = re.search(r"\[(.+?)\]", v)
        return match.group(1) if match else None
        
    if v == "subclonal":
        return "subclonal"
        
    return None

def process_vcf_dataframe(df, time_analysis=False):
    df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'INFO']
    df['ID'] = df['ID'].apply(lambda x: x[0] if x else '.')
    df['ALT'] = df['ALT'].apply(lambda x: x[0].value if x else '.')

    get_info = lambda info, key: info.get(key, '.')
    fields = ['MajCN', 'MinCN', 'MutCN', 'context', 'CLS', 'inWGDregion', 'powr']
    for f in fields:
        df[f] = df['INFO'].apply(lambda info: get_info(info, f))

    df['CLS'] = df['CLS'].apply(extract_clonality)

    # Keep only rows with valid context and clonality
    df = df[(df['context'] != '.') & df['CLS'].notna()]

    if time_analysis:
        df = df[(df['MajCN'] == 2) & (df['inWGDregion'] == 'TRUE')]

    df.rename(columns={'powr': 'pow'}, inplace=True)
    return df
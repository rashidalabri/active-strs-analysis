def clean_chroms(df, chrs):
    return df[df['chr'].isin(chrs)]

def clean_motif_len(df, min, max):
    return df[(df['motif'].str.len() >= min) & (df['motif'].str.len() <= max)]

def clean_ref_len(df, min, max):
    ref_len = df['stop'] - df['start']
    return df[(ref_len >= min) & (ref_len <= max)]

def apply_filters(df, config):
    df = clean_chroms(df, config['filter_contigs'])

    filter_motif_len = config['filter_motif_len']
    df = clean_motif_len(df, filter_motif_len['min'], filter_motif_len['max'])

    filter_ref_len = config['filter_ref_len']
    df = clean_ref_len(df, filter_ref_len['min'], filter_ref_len['max'])

    df['motif'] = df['motif'].str.upper()

    return df

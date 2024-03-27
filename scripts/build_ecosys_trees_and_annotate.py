#!/usr/bin/env python3
'''Split S3 clusters by ecosystem, build trees, and generate iTOL annotations'''

import os, sys, csv, json, argparse, subprocess
import pandas as pd
from Bio import SeqIO

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--s3_clusters', type=str, help='S3 clusters file', default='/groups/banfield/scratch/projects/environmental/spot/int/2023/assembly.d/S3_diversity/results/S3c/all_S3c_pfam_clusters.txt')
    parser.add_argument('--s3_fasta', type=str, help='S3 fasta file', default='/groups/banfield/scratch/projects/environmental/spot/int/2023/assembly.d/S3_diversity/results/S3c/all_S3c_pfam_centroids.faa')
    parser.add_argument('--s3_cluster_NR', type=str, help='NR annotations for S3 clusters', default='/groups/banfield/scratch/projects/environmental/spot/int/2023/assembly.d/S3_diversity/results/S3c/all_S3c_centroids_NR_taxa.tsv')
    parser.add_argument('--gg_phylum_colors', type=str, help='GGkBase colors for phylum', default='./S3c_full_WY/data/ggkbase_color_scheme_phylum.csv')
    parser.add_argument('--phylogeny_level', type=str, help='Phylogeny level to use for iTOL annotation and alignment split if specified.', default='phylum')
    parser.add_argument('--split_by_phylogeny_level', action='store_true', help='Construct alignments based on phylogeny level instead of ecosystem', default=True)
    parser.add_argument('--ref_aln', type=str, help='Reference alignment file', default='./S3c_full_WY/data/LHUG_S3_B+A.faa')
    parser.add_argument('--eaf_loc_file', type=str, help='File listing ecosystem to EAF file locations.', default='./S3c_full_WY/data/eco_to_eaf.csv')
    parser.add_argument('--run_fasttree', action='store_true', help='Run FastTree on the merged files.', default=False)
    parser.add_argument('--filter_by_eaf', action='store_true', help='Remove clusters without calculated EAF.', default=False)
    parser.add_argument('--min_eaf', type=float, help='Minimum EAF value to include in the phylogeny', default=0.0)
    parser.add_argument('--output_dir', type=str, help='Output directory', default='./S3c_full_WY/')
    return parser.parse_args(args)

TIME_DICT = {'April': 4, 'May': 5, 'June': 6, 'July': 7,'September': 9}
COLORS = ["#294A3A", "#294A3A", "#294A3A", "#294A3A", "#294A3A",
          "#6B3C22", "#6B3C22", "#6B3C22", "#6B3C22", "#6B3C22",
          "#96711A", "#96711A", "#96711A", "#96711A", "#96711A"]

COLORS = {"0-10cm": "#294A3A", "20-30cm": "#6B3C22", "50-80cm": "#96711A"}

# Load phylum dict from ./S3c_full_WY/data/name_map.json
PHYLUM_DICT = json.load(open('./S3c_full_WY/data/name_map.json'))

def read_S3_clusters(path_to_clusters):
    '''Read the S3 clusters file into a dictionary of representative genes to list of genes within that cluster.'''
    gene_dict = dict()
    with open(path_to_clusters, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == 'C':
                continue
            if row[0] == 'S':
                gene_dict[row[8]] = [row[8]]
            if row[0] == 'H':
                gene_dict[row[9]].append(row[8])
    return gene_dict

def read_NR_annotations(path_to_NR, phylogeny_levels=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
    '''Read the NR annotations for S3 clusters into a dictionary of representative genes to NR annotations.'''
    temp_dict, gene_dict = dict(), dict()
    with open(path_to_NR, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            temp_dict[row[0]] = (row[3], row[4])
    #For each entry in gene_dict, split the lineage and lineagerank by ; then zip them together 
    #and add columns named after the lineage_rank splits with corresponding lineage values
    for gene, (lineage, lineagerank) in temp_dict.items():
        lineage_split = lineage.split(';')
        lineagerank_split = lineagerank.split(';')
        for rank, line in zip(lineagerank_split, lineage_split):
            if gene not in gene_dict:
                gene_dict[gene] = {}
            gene_dict[gene][rank] = line
    # Make a dataframe from the dictionary
    gene_df = pd.DataFrame.from_dict(gene_dict, orient='index')
    gene_df = modify_na_values(gene_df, phylogeny_levels)
    return gene_df

def modify_na_values(df, column_names):
    for i in range(1, len(column_names)):
        # Replace the NA values with Unclassified_ + the previous column value
        df[column_names[i]] = df[column_names[i]].fillna('Unclassified_' + df[column_names[i-1]])

    for column in df.columns:
        df[column] = df[column].str.replace(r'(Unclassified_)+', 'Unclassified ', regex=True)

    return df[column_names]

def split_clusters(s3_clusters_dict, do_not_split=False):
    '''Split the clusters by ecosystem and return a dictionary of ecosystem to clusters.'''
    eco_dict = dict()
    if do_not_split:
        eco_dict['ALL'] = list(s3_clusters_dict.keys())
        return eco_dict
    for cl, genes in s3_clusters_dict.items():
        gene_ecos = set([x[10] for x in genes])
        if len(gene_ecos) > 1:
            for eco in gene_ecos:
                if eco not in eco_dict:
                    eco_dict[eco] = []
                eco_dict[eco].append(cl)
            continue
        ecosys = list(gene_ecos)[0]
        if ecosys not in eco_dict:
            eco_dict[ecosys] = []
        eco_dict[ecosys].append(cl)
    return eco_dict

def read_eaf_files(eaf_loc_file):
    '''Read the eaf files into named dictionary'''
    eaf_dict = dict()
    with open(eaf_loc_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            eaf_list = [*csv.DictReader(open(row[1]))]
            eaf_dict[row[0]] = [dict(x, ecosys=row[0]) for x in eaf_list]
    return eaf_dict

def parse_eaf_per_ecosystem(eaf_dict):
    '''Parse the eaf dict into a dictionary of gene to eaf'''
    eco_gene_eaf, eco_cores = dict(), dict()
    for eco, eaf_list in eaf_dict.items():
        gene_to_eaf, all_cores = dict(), list()
        for entry in eaf_list:
            if entry['cluster_gene'] not in gene_to_eaf:
                gene_to_eaf[entry['cluster_gene']] = {}
            gene_to_eaf[entry['cluster_gene']][entry['core']] = entry['mean_resampled_EAF']
            all_cores.append(entry['core'])
        eco_gene_eaf[eco] = gene_to_eaf
        eco_cores[eco] = list(set(all_cores))
    return eco_gene_eaf, eco_cores

def parse_eaf_by_taxonomy(eaf_list, lineage='phylum'):
    '''Parse the eaf dict into a dictionary of phylum to gene to eaf'''
    taxon_gene_eaf, eco_cores, all_cores = dict(), dict(), list()
    for entry in eaf_list:
        if entry[lineage] not in taxon_gene_eaf:
            taxon_gene_eaf[entry[lineage]] = {}
        if entry['cluster_gene'] not in taxon_gene_eaf[entry[lineage]]:
            taxon_gene_eaf[entry[lineage]][entry['cluster_gene']] = {}
        taxon_gene_eaf[entry[lineage]][entry['cluster_gene']][entry['core']] = entry['mean_resampled_EAF']
        all_cores.append(entry['core'])
    eco_cores['ALL'] = list(set(all_cores))
    return taxon_gene_eaf, eco_cores

def generate_eco_to_seqio_objects(s3_fasta, eco_to_genes, eaf_dict):
    '''Generate a dictionary of ecosystem to SeqIO objects for each gene in the ecosystem.
       Adds the phylogeny in the description field.'''
    s3_seqs = SeqIO.to_dict(SeqIO.parse(s3_fasta, "fasta"))
    eco_to_seqdict, eco_to_seqdict_filt = dict(), dict()
    for eco, genes in eco_to_genes.items():
        for gene in genes:
            if gene not in s3_seqs:
                continue
            if eco not in eco_to_seqdict:
                eco_to_seqdict[eco], eco_to_seqdict_filt[eco] = list(), list()
            # Find the gene in eaf_dict
            cluster_seq_object = s3_seqs[gene]
            for eaf in eaf_dict[eco]:
                if gene == eaf['cluster_gene']:
                    phylogeny = f"{eaf['superkingdom']}_{eaf['phylum']}_{eaf['class']}_{eaf['order']}"
                    cluster_seq_object.description = phylogeny
                    eco_to_seqdict_filt[eco].append(cluster_seq_object)
                    break
            eco_to_seqdict[eco].append(cluster_seq_object)
    return eco_to_seqdict, eco_to_seqdict_filt

def generate_tax_to_seqio_objects(s3_fasta, genes_list, eaf_list, lineage='phylum'):
    '''Generate a dictionary of taxonomy to SeqIO objects for each gene in the ecosystem.
       Adds the phylogeny in the description field.'''
    s3_seqs = SeqIO.to_dict(SeqIO.parse(s3_fasta, "fasta"))
    tax_to_seqdict, tax_to_seqdict_filt = dict(), dict()
    for gene in genes_list:
        if gene not in s3_seqs:
            continue
        eaf = [x for x in eaf_list if x['cluster_gene'] == gene]
        if len(eaf) == 0:
            continue
        tax_rank = eaf[0][lineage].replace(' ','_').replace('Candidatus_','')
        if tax_rank not in tax_to_seqdict:
            tax_to_seqdict[tax_rank], tax_to_seqdict_filt[tax_rank] = list(), list()
        cluster_seq_object = s3_seqs[gene]
        phylogeny = f"{eaf[0]['superkingdom']}_{eaf[0]['phylum']}_{eaf[0]['class']}_{eaf[0]['order']}"
        cluster_seq_object.description = phylogeny
        tax_to_seqdict[tax_rank].append(cluster_seq_object)
        tax_to_seqdict_filt[tax_rank].append(cluster_seq_object)
    return tax_to_seqdict, tax_to_seqdict_filt 

def save_split_seq_files(eco_to_seqs, output_dir, run=False):
    files = list()
    for eco, seqs in eco_to_seqs.items():
        output_tagged_dir = f"{output_dir}/{eco}"
        if run:
            # Check if dir exist, make it if it doesn't
            if not os.path.exists(output_tagged_dir):
                os.makedirs(output_tagged_dir)
            with open(f'{output_tagged_dir}/{eco}.faa', 'w') as output_file:
                SeqIO.write(seqs, output_file, "fasta")
        files.append(f'{output_tagged_dir}/{eco}.faa')
    return files
def split_ref_aln(ref_aln, list_taxa, output_dir):
    """Split the reference alignment in separate alignments for each taxa in list_taxa"""
    ref_aln_dict = SeqIO.to_dict(SeqIO.parse(ref_aln, "fasta"))
    tax_to_ref = dict()
    for taxon in list_taxa:
        aln_to_save = [ref_aln_dict[x] for x in ref_aln_dict if taxon in x]
        if len(aln_to_save) == 0 and taxon in PHYLUM_DICT:
            aln_to_save = [ref_aln_dict[x] for x in ref_aln_dict if PHYLUM_DICT[taxon] in x]
        if len(aln_to_save) == 0:
            aln_to_save = [ref_aln_dict[x] for x in ref_aln_dict if 'Escherichia' in x]
        tax_to_ref[taxon] = f"{output_dir}{taxon}/ref.faa"
        with open(f"{output_dir}{taxon}/ref.faa", "w") as output_file:
            SeqIO.write(aln_to_save, output_file, "fasta")
    return tax_to_ref
    
def mafft_addfull(ecosys, eco_ref_aln, output_dir, run=False, split_ref=False):
    '''Run mafft addfull to add ecosys specific sequences to the reference alignment. 
       Yields the names of merged files.'''
    for eco in ecosys:
        # Get reference filename
        ref_name = eco_ref_aln[eco].split('/')[-1].split('.')[0]
        output_tagged_dir = f"{output_dir}{eco}/"
        # Using mafft add ecosys files to ref_aln
        if run:
            merge_file = open(f"{output_tagged_dir}{ref_name}+{eco}.faa", "w")
            subprocess.run(['mafft', '--addfull', f'{output_tagged_dir}{eco}.faa', eco_ref_aln[eco]], check=True, stdout=merge_file, text=True)
        yield f"{output_tagged_dir}{ref_name}+{eco}.faa"

def cut_gaps_trimal(merged_files, run=False):
    '''Run trimal to remove gappy columns from the merged files. Yields the names of the cut files.'''
    for file in merged_files:
        if run:
            subprocess.run(['trimal', '-in', file, '-out', file.replace('.faa', '_cg.faa'), '-gappyout'], check=True, text=True)
        yield file.replace('.faa', '_cg.faa')

def run_fasttree(merged_cut_files, run=False):
    '''Run FastTree on the cut files.'''
    for file in merged_cut_files:
        if run:
            subprocess.run(['FastTree', '-gamma', '-log', f'{file.replace(".faa",".log")}', '-out', f'{file.replace(".faa",".nwk")}', file], check=True, text=True)
    return True

def create_EAF_itol_annotation(eco_gene_eaf, eco_cores, output_dir, split_by_tax=False):
    for eco, gene_eaf in eco_gene_eaf.items():
        #Substitute time for month and sort by time then depth
        if split_by_tax:
            # Sort by last element
            eco_core = eco_cores['ALL']
        else:
            eco_core = eco_cores[eco]
        eco_cores_time = [x for x in sorted(eco_core, key=lambda x: TIME_DICT[x.split('_')[0]])]
        eco_cores_sorted = [x for x in sorted(eco_cores_time, key=lambda x: x.split('_')[1])]
        if split_by_tax:
            eco_cores_sorted = [x for x in sorted(eco_cores_sorted, key=lambda x: x.split('_')[-1])]
        depths = [x.split('_')[1] for x in eco_cores_sorted]
        times = list(set([x.split('_')[0] for x in eco_cores_sorted]))
        with open(f"{output_dir}{eco}_EAF_anno.txt", "w") as output_file:
            output_file.write("DATASET_EXTERNALSHAPE\n")
            output_file.write("SEPARATOR COMMA\n")
            output_file.write(f"DATASET_LABEL,EAF {eco} by time and depth\n")
            output_file.write("COLOR,#efa000\n")
            output_file.write(f"FIELD_COLORS,{','.join([COLORS[x] for x in depths])}\n")
            output_file.write(f"FIELD_LABELS,{','.join(eco_cores_sorted)}\n")
            output_file.write("SHAPE_SPACING,-2\n")
            output_file.write("SIZE_FACTOR,1.2\n")
            output_file.write("DASHED_LINES,1\n")
            output_file.write("DATA\n")
            for gene, eaf_dict in gene_eaf.items():
                data_string = f"{gene},"
                for core in eco_cores_sorted:
                    if core not in eaf_dict:
                        eaf_value = 0
                    else:
                        eaf_value = (abs(float(eaf_dict[core]))+float(eaf_dict[core]))/2
                    data_string += f"{eaf_value},"
                output_file.write(f"{data_string}\n")
    return True

def create_phylogeny_itol_annotation(eco_to_seqs, s3_NR_df, output_dir, gg_ph_colors, phylogeny_level):
    # Read in the ggkbase color scheme for phylum as a dictionary
    ggkbase_phylum_colors = pd.read_csv(gg_ph_colors)
    ggkbase_phylum_colors = pd.concat([ggkbase_phylum_colors, pd.DataFrame({'Phylum': 'Unclassified', 'Color': '#757575'}, index=[0])])
    for eco, seqs in eco_to_seqs.items():
        with open(f"{output_dir}{eco}_phylogeny_anno.txt", "w") as output_file:
            output_file.write("DATASET_COLORSTRIP\n")
            output_file.write("SEPARATOR COMMA\n")
            output_file.write(f"DATASET_LABEL,Phylogeny {eco}\n")
            output_file.write("COLOR,#4b5fde\n")
            output_file.write(f"COLOR_BRANCHES,0\n")
            output_file.write(f"DATA\n")
            for seq in seqs:
                gene = seq.id
                lineage = s3_NR_df.loc[gene, phylogeny_level]
                color = ggkbase_phylum_colors[ggkbase_phylum_colors['Phylum'] == lineage]['Color'].tolist()
                if len(color) == 0:
                    color = ['#757575']
                output_file.write(f"{gene},{color[0]}\n")
    return True
def debugger_is_active() -> bool:
    """Return if the debugger is currently active"""
    return hasattr(sys, 'gettrace') and sys.gettrace() is not None

def main(args):
    cl_args = parse_args(args)
    # Save arguments including defaults as text file in the output dir if running from command line
    if not debugger_is_active():
        import json
        with open(f"{cl_args.output_dir}args.json", 'w') as f: 
            json.dump(vars(cl_args), f)

    s3_clusters_dict = read_S3_clusters(cl_args.s3_clusters)
    s3_NR_df = read_NR_annotations(cl_args.s3_cluster_NR)
    eco_to_genes = split_clusters(s3_clusters_dict, do_not_split=cl_args.split_by_phylogeny_level)
    eaf_dict = read_eaf_files(cl_args.eaf_loc_file)
    if cl_args.split_by_phylogeny_level:
        new_eaf_dict = {'ALL': []}
        for eco, eaf_list in eaf_dict.items():
            update_core = [dict(x, core=f"{x['core']}_{x['ecosys']}") for x in eaf_list]
            new_eaf_dict['ALL'].extend(update_core)
        eaf_dict = new_eaf_dict
    if cl_args.filter_by_eaf:
        for eco, eaf_list in eaf_dict.items():
            eaf_dict[eco] = [x for x in eaf_list if float(x['mean_resampled_EAF']) >= cl_args.min_eaf]
    if cl_args.split_by_phylogeny_level:
        eco_gene_eaf, eco_cores = parse_eaf_by_taxonomy(eaf_dict['ALL'], lineage=cl_args.phylogeny_level)
        eco_to_seqs, eco_to_seqdict_filt = generate_tax_to_seqio_objects(cl_args.s3_fasta, eco_to_genes['ALL'], eaf_dict['ALL'], lineage=cl_args.phylogeny_level)
    else:
        eco_gene_eaf, eco_cores = parse_eaf_per_ecosystem(eaf_dict)
        eco_to_seqs, eco_to_seqdict_filt = generate_eco_to_seqio_objects(cl_args.s3_fasta, eco_to_genes, eaf_dict)
    # Output the sequences to two files
    if cl_args.filter_by_eaf:
        eco_to_seqs = eco_to_seqdict_filt
    split_seqs = save_split_seq_files(eco_to_seqs, cl_args.output_dir, run=cl_args.run_fasttree)
    # Run mafft addfull
    tax_to_ref = {x: cl_args.ref_aln for x in eco_to_seqs.keys()}
    if cl_args.split_by_phylogeny_level:
        tax_to_ref = split_ref_aln(cl_args.ref_aln, eco_to_seqs.keys(), cl_args.output_dir)
    merged_files = mafft_addfull(eco_to_seqs.keys(), tax_to_ref, cl_args.output_dir, run=cl_args.run_fasttree)
    merged_cut_files = cut_gaps_trimal(merged_files, run=cl_args.run_fasttree)
    run_fasttree(merged_cut_files, run=cl_args.run_fasttree)
    create_EAF_itol_annotation(eco_gene_eaf, eco_cores, cl_args.output_dir, cl_args.split_by_phylogeny_level)
    create_phylogeny_itol_annotation(eco_to_seqs, s3_NR_df, cl_args.output_dir, cl_args.gg_phylum_colors, cl_args.phylogeny_level)
    print(eco_gene_eaf)

if __name__ == '__main__':
    main(sys.argv[1:])
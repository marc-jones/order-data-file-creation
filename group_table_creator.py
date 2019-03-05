import re

base_filepath = ('/run/media/marcjones/Seagate Expansion Drive/' +
    'morris_cluster/')
annotation_filepath = (base_filepath + '2016_08_02_final_sequencing_results/' +
    'out_2016_08_15_19_05/merged_genes_annotation.tsv')
gtf_filepath = (base_filepath + '2016_08_02_final_sequencing_results/' +
    'out_2016_08_15_19_05/merged_asm/merged.gtf')

record_details_output_filepath = ('/home/marcjones/Documents/' +
    'order-data-files/record_details.tsv')
groups_output_filepath = ('/home/marcjones/Documents/' +
    'order-data-files/groups.tsv')

record_details_dict = {}
groups_dict = {}

with open(gtf_filepath) as f:
    for line in f:
        (chromosome, cufflinks, feature_type, start, end, dot1, strand,
            dot2, additional_info) = line.strip().split('\t')
        gene = additional_info.split('; ')[0]
        gene = gene.replace('"', '').replace('gene_id ', '')
        genoscope_link_string = ('<a href="http://www.genoscope.cns.fr/' +
            'brassicanapus/cgi-bin/gbrowse/colza/?name=' + chromosome +
            '%3A' + start + '..' + end + '" target="_blank">' + gene + '</a>')
        record_details_dict.setdefault(gene, {
            'name': gene,
            'label_tooltip': chromosome,
            'Brassica napus gene': genoscope_link_string,
            'Chromosome': chromosome,
            'Start (bp)': int(start),
            'End (bp)': int(end)
        })
        if int(start) < record_details_dict[gene]['Start (bp)']:
            record_details_dict[gene]['Start (bp)'] = int(start)
        if int(end) > record_details_dict[gene]['End (bp)']:
            record_details_dict[gene]['End (bp)'] = int(end)

with open(annotation_filepath) as f:
    header = f.readline()
    for line in f:
        if len(line.strip().split('\t')) == 10:
            (gene, chromosome, start, end, igv_format, identity, hsp_bit_score,
                length_of_hsp, agi_code, symbols) = line.strip().split('\t')
        elif len(line.strip().split('\t')) == 9:
            (gene, chromosome, start, end, igv_format, identity, hsp_bit_score,
                length_of_hsp, agi_code) = line.strip().split('\t')
            symbols = ''
        tair_link_string = ('<a href="https://www.arabidopsis.org/servlets/' +
            'TairObject?name=' + re.sub('\.[0-9]', '', agi_code) +
            '&amp;type=locus" target="_blank">' + agi_code +'</a>')
        # Check that it's consistent with the record details dictionary
        assert(
            record_details_dict[gene]['name'] == gene and
            record_details_dict[gene]['label_tooltip'] == chromosome and
            record_details_dict[gene]['Chromosome'] == chromosome and
            record_details_dict[gene]['Start (bp)'] == int(start) and
            record_details_dict[gene]['End (bp)'] == int(end)
        )
        groups_dict.setdefault(gene, [])
        groups_dict[gene].append({
            'name': gene,
            'group': agi_code,
            'nicknames': ','.join(symbols.split(', ')),
            'label_tooltip': chromosome,
            'label_colour': 'default',
            'Arabidopsis gene': tair_link_string,
            'Abbreviation': symbols,
            'BLAST Identity': float(identity),
            'BLAST HSP Bit Score': float(hsp_bit_score),
            'BLAST HSP Length': int(float(length_of_hsp))
        })

for gene in groups_dict.keys():
    bit_scores = [row_dict['BLAST HSP Bit Score']
        for row_dict in groups_dict[gene]]
    agi_codes = [row_dict['group']
        for row_dict in groups_dict[gene]]
    top_gene = sorted(zip(bit_scores, agi_codes), reverse=True)[0][1]
    for row_idx in range(len(groups_dict[gene])):
        if (re.sub('\.[0-9]', '', top_gene) in
            groups_dict[gene][row_idx]['group']):
            groups_dict[gene][row_idx]['label_colour'] = 'warning'
        if top_gene == groups_dict[gene][row_idx]['group']:
            groups_dict[gene][row_idx]['label_colour'] = 'success'

record_details_headers = ['name', 'label_tooltip', 'Brassica napus gene',
    'Chromosome', 'Start (bp)', 'End (bp)']

with open(record_details_output_filepath, 'w') as f_out:
    f_out.write('\t'.join(record_details_headers) + '\n')
    for gene in record_details_dict.keys():
        f_out.write('\t'.join([str(record_details_dict[gene][key])
            for key in record_details_headers]) + '\n')

groups_headers = ['name', 'group', 'nicknames', 'label_tooltip', 'label_colour',
    'Arabidopsis gene', 'Abbreviation', 'BLAST Identity',
    'BLAST HSP Bit Score', 'BLAST HSP Length']

with open(groups_output_filepath, 'w') as f_out:
    f_out.write('\t'.join(groups_headers) + '\n')
    for gene in groups_dict.keys():
        for row_dict in groups_dict[gene]:
            f_out.write('\t'.join([str(row_dict[key])
                for key in groups_headers]) + '\n')

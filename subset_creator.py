import csv
import os
import re
import shutil

flor_id_filepath = 'flor_id_combined_2015_01_27.csv'

order_data_filepath = '~/Documents/order-data-files'
order_data_subset_filepath = '~/Documents/order-data-files-subset'

agi_set = set([])

csv_reader = csv.reader(open(flor_id_filepath))
# Skip header
csv_reader.next()
for row in csv_reader:
    agi_set.add(row[6])

agi_to_xloc = {}
xloc_to_agi = {}
with open(order_data_filepath + '/groups.tsv') as groups:
    # Skip header
    groups.readline()
    for line in groups:
        xloc = line.strip().split('\t')[0]
        agi = re.sub('\.[0-9]', '', line.strip().split('\t')[1])
        agi_to_xloc.setdefault(agi, set([]))
        agi_to_xloc[agi].add(xloc)
        xloc_to_agi.setdefault(xloc, set([]))
        xloc_to_agi[xloc].add(agi)

xloc_set = set([])
xloc_set_length = -1

new_agi_set = agi_set.copy()

print('Total XLOC number: ' + str(len(xloc_to_agi.keys())))

while xloc_set_length < len(xloc_set):
    xloc_set_length = len(xloc_set)
    print(str(xloc_set_length))
    new_xloc_set = set([])
    for agi in new_agi_set:
        if agi in agi_to_xloc.keys():
            new_xloc_set = new_xloc_set.union(
                agi_to_xloc[agi].difference(xloc_set))
            xloc_set = xloc_set.union(agi_to_xloc[agi])
    new_agi_set = set([])
    for xloc in new_xloc_set:
        if xloc in xloc_to_agi.keys():
            new_agi_set = new_agi_set.union(
                xloc_to_agi[xloc].difference(agi_set))
            agi_set = agi_set.union(xloc_to_agi[xloc])

if os.path.isdir(order_data_subset_filepath):
    shutil.rmtree(order_data_subset_filepath)

user_content_name = 'user_content'
shutil.copytree(order_data_filepath + '/' + user_content_name,
    order_data_subset_filepath + '/' + user_content_name)

for filename in ['README.md', 'website_information.yaml', 'plot_regions.tsv']:
    shutil.copy2(order_data_filepath + '/' + filename,
        order_data_subset_filepath + '/')

fasta_name = 'genes.fasta'
with open(order_data_subset_filepath + '/' + fasta_name, 'w') as out_fasta:
    with open(order_data_filepath + '/' + fasta_name) as in_fasta:
        record = False
        for line in in_fasta:
            if line.startswith('>'):
                if line.strip().replace('>', '') in xloc_set:
                    record = True
                else:
                    record = False
            if record:
                out_fasta.write(line)

def filter_tsv_table(filename, column_idx):
    with open(order_data_subset_filepath + '/' + filename, 'w') as out_fasta:
        with open(order_data_filepath + '/' + filename) as in_fasta:
            out_fasta.write(in_fasta.readline())
            for line in in_fasta:
                row = line.strip().split('\t')
                if row[column_idx] in xloc_set:
                    out_fasta.write(line)

filter_tsv_table('groups.tsv', 0)
filter_tsv_table('record_details.tsv', 0)
filter_tsv_table('time_series_data.tsv', 0)

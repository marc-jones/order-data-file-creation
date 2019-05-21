import urllib2
import json

domain_name = 'http://order.jic.ac.uk'
domain_name = 'http://192.168.99.10'

xloc_set = set([])

with open('/home/marcjones/Documents/order-data-files/genes.fasta') as f:
    for line in f:
        if line.startswith('>'):
            xloc_set.add(line.strip().replace('>', ''))

with open('/home/marcjones/Downloads/verify_' + domain_name.replace('/', '') + '.tsv', 'w') as out:
    for xloc in sorted(list(xloc_set)):
        print(xloc)
        try:
            order_request = urllib2.urlopen(
                domain_name + '/postcheckboxchange?names=' + xloc + ';')
        except:
            out.write(xloc + '\tfailed_to_connect\n')
            continue
        if order_request.getcode() == 200:
            json_data = json.loads(order_request.read())
            if not len(json_data[0]['measurements']) == 27:
                out.write(xloc + '\tnot_all_measurements\n')
        else:
            out.write(xloc + '\tfailed_with_error_' + str(order_request.getcode()) + '\n')

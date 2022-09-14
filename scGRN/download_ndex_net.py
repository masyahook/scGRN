import argparse
import sys
import os

import json
import ndex2
import networkx as nx

if __name__ == '__main__':
    
    # Setting working directory as home
    home_dir = os.path.expanduser('~')
    os.chdir(os.path.expanduser('~'))
    
    parser = argparse.ArgumentParser(description='Download NDEx network using UUID.')
    parser.add_argument('-i', '--id', type=str, help='The UUID of the NDExx network', required=True)
    args = parser.parse_args()

    # Create NDEx2 python client
    client = ndex2.client.Ndex2()

    # Download BioGRID: Protein-Protein Interactions (SARS-CoV) from NDEx
    # http://ndexbio.org/viewer/networks/669f30a3-cee6-11ea-aaef-0ac135e8bacf
    client_resp = client.get_network_as_cx_stream(args.id)

    # Convert downloaded network to NiceCXNetwork object
    net_cx = ndex2.create_nice_cx_from_raw_cx(json.loads(client_resp.content))

    # Create Networkx network
    g = net_cx.to_networkx()
    print(g.edges())
    print(g.nodes())

    print('Downloaded network with name: ' + str(g))
    print('Number of nodes: ' + str(g.number_of_nodes()))
    print('Number of edges: ' + str(g.number_of_edges()))
    print('Network annotations: ' + str(g.graph))
    
    nx.write_gpickle(g, f'res/ndex/{str(g)}')
    print()
    print(f'Saved to ~/res/ndex/{str(g)}')
    print(f'To load run: G = nx.read_gpickle("res/ndex/{str(g)}")')

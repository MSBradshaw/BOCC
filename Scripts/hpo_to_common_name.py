import os
import sys
import BOCC
import networkx
import obonet
import requests
import json


"""
STDIN: list of HPO terms
STDOUT: list of HPO terms and their common names
"""

url = 'http://purl.obolibrary.org/obo/hp.obo'
G = obonet.read_obo(url)
G_dict = G.nodes(data=True)

for item in sys.stdin:
    item = item.strip()
    print(item, G.nodes(data=True)[item]['name'])



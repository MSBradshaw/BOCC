import networkx as nx
import matplotlib.pyplot as plt
import obonet
import BOCC
import pandas as pd
import numpy as np
import math
from webweb import Web

G = nx.read_edgelist('Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt')

nx.draw(G, node_size=1)
plt.show()

x = ["a", "b""c", "d"]

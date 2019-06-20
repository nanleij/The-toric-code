# The-toric-code
The Toric Code (surface code for quantum error correction) together with the MWPM decoder (need a package) and the union-find decoder programmed in python 3

Files:
Toric_code.py       : contains all basic functions to make the grid, simulate errors and check if an correction is correct
MWPM_decoder.py     : contains functions to simulate the MWPM decoder on the toric code (package networkx needed)
Peeling_decoder.py  : contains functions needed to simulate the peeling decoder on the toric code, including functions to generate erasure errors
UF_decoder.py       : contains functions needed to simulate the Union-Find decoder on the toric code
simulate.py         : can be used to simulate one of the decoders on the toric code N times, shows progress while running, saves data and makes a plot of the results   

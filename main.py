from Number import Number
from GpnN import GpnN

import matplotlib.pyplot as plt
import time

start_time = time.time()
G2_22 = GpnN(2,0,5) #initialization
G2_22.generate_numbers()
G2_22.ODESols()
#G2_22.matrix()
#G2_22.representation_tree()
#G2_22.console_printing()

#G2_22.representation_tree()

#G2_22.matrix()

#G2_22.ODESols()
G2_22.export_gif()

print("--- %s seconds ---" % (time.time() - start_time))
  

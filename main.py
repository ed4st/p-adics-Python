from Number import Number
from GpnN import GpnN

import matplotlib.pyplot as plt



G2_22 = GpnN(2,-2,0) #initialization
G2_22.generate_numbers()
#G2_22.representation_tree()
#G2_22.console_printing()
'''
i_0 = 0
for i in G2_22.numbers:
    j_0=0
    for j in G2_22.numbers:
        
        
        print(i_0,"-",j_0,i.digits,"-",j.digits," = ", G2_22.p_sub(i,j).digits,"; order: ", G2_22.p_sub(i,j).order(),"; Norm: ", G2_22.p_sub(i,j).norm())
        j_0 += 1
    print("-------------------------------")
    i_0 += 1'''
#G2_22.representation_tree()

G2_22.matrix()

#G2_22.ODESols()
#G2_22.export_gif()

  

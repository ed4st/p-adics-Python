#NOTE: We suppose that all Numbers in GmMp have the same Number of digits
#therefore, following algorithms follow this fact 
from Number import Number
import random
from random import randint
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import networkx as nx
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import gif
from matplotlib import cm
from numpy import linspace

class GpnN:
  p = 0 #prime p
  __m = 0  #minimum positive index where the sum can start 
  __M = 0  #maximum positive index where the sum finish 
  n = 0 #negative number such that the norm of the number is greater that p^n
  N = 0 #positive number such that the norm of the number is lower that p^N
  numbers = []


  def __init__(self,p ,n ,N):
    self.p = p 
    self.n = n 
    self.N = N
    self.__m = self.N
    self.__M = -self.n 


#---------------Getters and setters-------------------------

  def getm(self):
    return self.__m
  def setm(self, m):
    self.__m = m

  def getM(self):
    return self.__M
  def setM(self, M):
    self.__M = M 
#-------------------------numbers_initialization-------------------------- 
  def generate_numbers(self):
    to_permute  =  [i for i in range(self.p)]
    perms  =  product(to_permute,repeat = self.__m + self.__M+1)
    i  =  0
    numbers_aux  =  []
    for p in perms:
      arr  =  []
      arr.extend(p) 
          
      numbers_aux.append(Number(self.p,self.n,self.N,arr)) 
          
      self.numbers  =  numbers_aux
      i  =  i + 1
#-------------------------GmMp_export-------------------------- 
  def GmMp_export(self, name):
    f  =  open(name + ".txt","w+")
    for num in self.numbers:
      f.write(" ".join(str(num.digits)) + "   "+str(num.norm()) +"\n")
    
    f.close()
#-------------------------console_printing--------------------------
  def console_printing(self):
    for i in self.numbers:
      i.show()
#-------------------------sum---------------------------
  #following function computes the sum of two p-adic 
  #numbers and truncates the result to the first M-m digits
  def p_sum(self, n1: Number, n2 : Number):
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]

    psum = []
    carry  =  0
    for i in range (len(n1.digits)):
      sum = n1.digits[i]+n2.digits[i]+ carry
      r = sum %n1.p
      psum.append(r)
      carry  =  int((sum)/n1.p)
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]
    return Number(n1.p, n1.n, n1.N, psum[::-1] )
#-------------------------substraction-----------------------
  #following function computes the substraction of two
  #p-adic numbers and truncates the result to the first M-m digits
  def p_sub(self,n1: Number, n2 : Number):
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    result = []#[0 for i in range(n1.digits)][1,0]
    carry  =  0
    for i in range (len(n1.digits)):
      sub  =  n1.digits[i]-n2.digits[i]+carry
      #sub = 0 - 1 + 0 = -1 
       
      r  =  sub%n1.p
      #r = -1%2 = 1 
      
      result.append(r)
       
      carry  =  int((sub-r)/n1.p)
      # int(-1/2)=0

    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    return Number(n1.p, n1.n, n1.N, result[::-1] )
#-------------------------product-------------------------
  #following function computes the multiplication of two
  #p-adic numbers and truncates the result to the first M-m digits
  def p_mul(self, n1: Number, n2 : Number):

    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  

    resultados  =  []
    for i in range(len(n2.digits)):
      resultado  =  []
      carry  =  0
      for j in range(len(n2.digits)):
        prod  =  n2.digits[j] * n1.digits[i] + carry
        resultado.append(prod % n2.p)
        carry  =  int(prod / n2.p)
      resultados.append(resultado)
      
    aux  =  Number(n1.p,n1.n,n1.N,resultados[0])
    for i in range(1,len(resultados),1):
      for j in range(i):
        resultados[i].pop(len(resultados[i])-1)
        resultados[i].insert(0,0)
      aux  = self.p_sum(aux,Number(n1.p,n1.n,n1.N,resultados[i]))
    aux.digits = aux.digits[::-1]
    n1.digits = n1.digits[::-1]  
    n2.digits = n2.digits[::-1]  
    return(aux)
#-------------------------division---------------------------
  #following function returns the solution of congruence equation ax = b(mod p)
  #where a,b are the first nonzero digit of divisor and dividend
  #respectively. This solution can be given because p is prime

  def first_dividend_digit_anihilator(self, dividend : Number,divisor : Number):
    dr  =  divisor.digits[0]
    dn  =  dividend.digits[0]
    if(dr == dn and(dr != 0 and dn != 0)):
      return 1
    elif(dn == 0 ):
      return 0
    else:          
      for i in range(dividend.p):
        if((dr*i-dn)and((dr*i-dn)%dividend.p)  ==  0):
          return i
        if(dr*i-dn == 0):
          return i
  #following funtion returns p-adic division when p is prime. 
  # In other case we cannot guarantee the correct result of 
  # first_dividend_digit_anihilator function
  def p_div(self,dividend: Number, divisor : Number):
    dividend_copy = dividend
    divisor_copy = divisor
    dividend_copy .digits = dividend_copy .digits[::-1]  
    divisor_copy .digits = divisor_copy .digits[::-1]  
    i = 0
    while(divisor_copy .digits[i] == 0):
      divisor_copy .digits.pop(0)
      divisor_copy .digits.append(0)
    resul = []
    for i in range(len(dividend_copy .digits)):
      r_temp  =  self.first_dividend_digit_anihilator(dividend_copy ,divisor_copy )
      c_tem_digits  =  [0 for i in range(len(dividend_copy .digits))]
      c_tem_digits[len(dividend_copy .digits)-1]  =  r_temp
      
      
      c_temp_number  =  Number(dividend_copy .p, dividend_copy .n, dividend_copy .N,c_tem_digits)
      divisor_copy_reverse = Number(divisor_copy .p,divisor_copy .n,divisor_copy .N,divisor_copy .digits[::-1])

      
      dividend_copy .digits = dividend_copy .digits[::-1] 
      product_aux = self.p_mul(c_temp_number,divisor_copy_reverse)
      substraction_aux = self.p_sub(dividend_copy ,product_aux)
      substraction_aux.digits.pop(len(substraction_aux.digits)-1)
      substraction_aux.digits.insert(0,0)

      dividend_copy  = substraction_aux
      dividend_copy .digits = dividend_copy .digits[::-1]
      resul.append(r_temp)
    
    return Number(dividend_copy .p,dividend_copy .n,dividend_copy .N,resul[::-1])     
#-------------------------inverse------------------------
  def p_inverse(self,n : Number):
    identity_digits  =  [0 for i in range(len(n.digits))]
    identity_digits[self.__M] = 1
    identity  =  Number(n.p,n.n,n.N,identity_digits)
    return self.p_div(identity,n)
#-------------------------GmMp_representation_tree---------------------------
  def representation_tree(self):
    try:
        import pygraphviz
        from networkx.drawing.nx_agraph import graphviz_layout
    except ImportError:
        try:
            import pydot
            from networkx.drawing.nx_pydot import graphviz_layout
        except ImportError:
            raise ImportError("This example needs Graphviz and either "
                              "PyGraphviz or pydot")
    
    G  =  nx.balanced_tree(self.p,self.__m + self.__M + 1)
    for u,v,d in G.edges(data = True):
      d['weight']  =  (v%self.p)/self.p #coloring edges               
    leaf_nodes  =  []
    for node in G.nodes():
        if nx.degree(G,node)  ==  1:
            leaf_nodes.append(node)


    
    norm  =  0.0
    nx.set_node_attributes(G,norm,'norm')
    for i in range(len(G.nodes)):
          G.nodes[i]['norm']  =  0.0
    
    for i in range(len(leaf_nodes)):
      G.nodes[leaf_nodes[i]]['norm']  =  self.numbers[i].norm()
      
    norms  =  []
    edges,weights  =  zip(*nx.get_edge_attributes(G,'weight').items())
    nodes,norms  =  zip(*nx.get_node_attributes(G,'norm').items())
    
    
    #Creating a list of p colors
    basic_colors = ['k','r','b','gray','green','c','y','m']
    if(self.p>len(basic_colors)-1):
      
      cm_linspace = linspace(0.0, 1.0, self.p-len(basic_colors))
      basic_colors = basic_colors + [cm.brg(x) for x in cm_linspace] 
    color_of_edges = []

    for i in range(0,len(edges),3):
      for j in range(self.p):
        color_of_edges.append(basic_colors[j])

    pos  =  graphviz_layout(G, prog = 'twopi', args = '')
    plt.figure(figsize = (8,8))

    name = 'G' + str(self.p) + '_' + str(abs(self.n)) + str(self.N)
    plt.title(name)
    
    norms_set = set(norms)
    unique_norms = []
    for i in norms_set: 
      unique_norms.append(i)
    #creating a list of colors based 
    # on how many diferents norms are
    start = 0.0
    stop = 1.0
    cm_subsection = linspace(start, stop, len(norms_set)) 
    norms_color_set = [ cm.jet(x) for x in cm_subsection ]
    
    unique_color = []  
    for i in norms_color_set:
      unique_color.append(i)
    
    unique_norms.sort()
    color_norms = []
    for i in norms:
      if(i == 0.0):
        color_norms.append('white')
      else:
        index = unique_norms.index(i)
        color_norms.append(unique_color[index])
    
    cm_aux = ListedColormap(unique_color)#creating a color map
                        
    nx.draw(G, pos, node_size = 100, alpha = 0.7, node_color  =  color_norms, edge_color  = color_of_edges,with_labels = False)
    plt.axis('equal')

    #vertical colorbar
    sm = plt.cm.ScalarMappable(cmap = cm_aux)
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.ax.set_yticklabels(unique_norms)
    cbar.set_label('Norm', rotation=270)
    plt.show()
    plt.savefig(name + '.png')    
    
#--------------------------Parisi_Matrix--------------------------
  #following function returns the (i,j)-th value of Parisi Matrix
  def fij(self, alpha, i: Number, j: Number):
    sub = self.p_sub(i,j)
    norm = sub.norm()
    order = sub.order()
    return order#1/(self.p_sub(i,j).norm()**alpha + 1)

  def matrix(self):
    W = []
    for i in self.numbers:
      aux = []
      for j in self.numbers:
        aux.append(self.fij(2,i,j))
      W.append(aux)
    
    for fila in W:
      print(fila)
    
    sum_row1 = sum(W[0]) 
    for i in range(len(W)):
      W[i] = [j*(sum_row1**(-1)) for j in (W[i])]
    #upgrading the diagonal of W
    row_copy = []
    row_copy = W[0].copy()
    row_copy.pop(0)
    w0 =  - sum(row_copy)
    for i in range(len(W)):  
      W[i][i] = w0 
      
    #matrix visualization
    fig, ax = plt.subplots()
    ax.matshow(W, cmap = plt.cm.get_cmap("jet"))
    for i in range(len(self.numbers)):
      for j in range(len(self.numbers)):
        c = W[i][j]
        #ax.text(i, j, f"{c:.2f}", va='center', ha='center')#text(i, j, f"{c:.2f}", va='center', ha='center')
    plt.show()
    return W
#-------------------------Solving the Master Equation System---------------------------
  def __model(self, u, t):
    W = self.matrix()
    dudts = []
    for i in range(len(W)):
      auxeq = 0
      for j in range(len(W)):
        auxeq += W[i][j]*u[j]
      dudts.append(auxeq)
    return dudts
  
  def ODESols(self):
    #initial condition
    
    #normal
    #u0 = np.random.normal(len(self.numbers)/2,1,len(self.numbers))
    
    #ones
    u0 = [1 for i in range(len(self.numbers))]
    
    #median = 0, variance = 1
    #u0 = np.random.normal(0,1,len(self.numbers))
    
    #random list from 0 to 1
    #u0 = random.sample(range(0,1),len(self.numbers))
    
    #time points
    t = np.linspace(0,1,5)

    u = odeint(self.__model,u0,t)
    return [u,t]
  #--------------------------Creating image transitions---------------
  @gif.frame
  def animate(self, boundary, time):
    try:
        import pygraphviz
        from networkx.drawing.nx_agraph import graphviz_layout
    except ImportError:
        try:
            import pydot
            from networkx.drawing.nx_pydot import graphviz_layout
        except ImportError:
            raise ImportError("This example needs Graphviz and either "
                              "PyGraphviz or pydot")
    
    G  =  nx.balanced_tree(self.p,self.__m + self.__M + 1)
    for u,v,d in G.edges(data = True):
      d['weight']  =  (v%self.p)/self.p #coloring edges               
    leaf_nodes  =  []
    for node in G.nodes():
        if nx.degree(G,node)  ==  1:
            leaf_nodes.append(node)
    
    norm  =  0.0
    nx.set_node_attributes(G,norm,'norm')
    for i in range(len(G.nodes)):
          G.nodes[i]['norm']  =  0.0
    
    for i in range(len(leaf_nodes)):
      G.nodes[leaf_nodes[i]]['norm']  =  boundary[i]
      
    norms  =  []
    edges,weights  =  zip(*nx.get_edge_attributes(G,'weight').items())
    nodes,norms  =  zip(*nx.get_node_attributes(G,'norm').items())
    
    
    #Creating a list of p colors
    basic_colors = ['k','r','b','gray','green','c','y','m']
    if(self.p>len(basic_colors)-1):
      
      cm_linspace = linspace(0.0, 1.0, self.p-len(basic_colors))
      basic_colors = basic_colors + [cm.brg(x) for x in cm_linspace] 
    color_of_edges = []

    for i in range(0,len(edges),3):
      for j in range(self.p):
        color_of_edges.append(basic_colors[j])

    pos  =  graphviz_layout(G, prog = 'twopi', args = '')
    plt.figure(figsize = (8,8))

    name = 'G' + str(self.p) + '_' + str(abs(self.n)) + str(self.N)
    plt.title(name)
    plt.text(2, 6, '$t =' + f"{time:.2f}" + '$', fontsize=15)

    norms_set = set(norms)
    unique_norms = []
    for i in norms_set: 
      unique_norms.append(i)
    #creating a list of colors based 
    # on how many diferents norms are
    start = 0.0
    stop = 1.0
    cm_subsection = linspace(start, stop, len(norms_set)) 
    norms_color_set = [ cm.jet(x) for x in cm_subsection ]
    
    unique_color = []  
    for i in norms_color_set:
      unique_color.append(i)
    
    unique_norms.sort()
    color_norms = []
    for i in norms:
      if(i == 0.0):
        color_norms.append('white')
      else:
        index = unique_norms.index(i)
        color_norms.append(unique_color[index])
    
    cm_aux = ListedColormap(unique_color)#creating a color map
                        
    nx.draw(G, pos, node_size = 100, alpha = 0.7, node_color  =  color_norms, edge_color  = color_of_edges,with_labels = False)
    plt.axis('equal')

    #vertical colorbar
    sm = plt.cm.ScalarMappable(cmap = cm_aux)
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label('', rotation=270)
    cbar.set_clim(-2.0, 2.0)


    #plt.show()

  def export_gif(self):
    frames = []
    u = self.ODESols()
    for u_i, t_i in zip(u[0],u[1]):
      print(u_i)
      print("------------------------")
      frame = self.animate(u_i,t_i)
      frames.append(frame)
    gif.save(frames,"ojala.gif",duration=3000)
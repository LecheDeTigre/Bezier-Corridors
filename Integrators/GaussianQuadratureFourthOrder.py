import numpy as np

def integrate(differential, limits):
    sum = 0
    # import pdb; pdb.set_trace()
    
    a = limits[0]
    b = limits[1]
    
    weights = [(18 + np.sqrt(30))/36, (18 + np.sqrt(30))/36, (18 - np.sqrt(30))/36, (18 - np.sqrt(30))/36]
    nodes = [np.sqrt(3/7 - 2/7*np.sqrt(6/5)), -np.sqrt(3/7 - 2/7*np.sqrt(6/5)), np.sqrt(3/7 + 2/7*np.sqrt(6/5)), -np.sqrt(3/7 + 2/7*np.sqrt(6/5))]
    
    for (weight, node) in zip(weights, nodes):
        sum += weight * differential(0.5*(b-a)*node + 0.5*(b+a))
        
    sum = 0.5*sum*(b-a)
    
    return sum
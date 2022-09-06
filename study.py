import sys
import pandas as pd
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

def plot_fit(x,y,a, b,x_label='x',y_label='y'):
    fig = plt.figure()
    ax = plt.gca()
        
    # Plot data points
    ax.scatter(x,y, color="red", marker="o", label="Original data")
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    # Plot linear regression line.
    y_pred = b*x**a
    ax.semilogy(x,y_pred, color="green", label="Fitted line "+"{:1.1e}".format(b)+"x^"+"{:1.1e}".format(a))

    # Set labels
    plt.legend(loc='best')
    plt.xlabel(x_label) 
    plt.ylabel(y_label) 

    plt.show()



csv_file=sys.argv[1]
df=pd.read_csv(csv_file)
df.columns = df.columns.str.replace(' ', '')
#print('inner A /outer=')
#print(df['inner']/df['outer'])
print('total outer = ',df['outer'][0:5].sum())
print('total inner A/total outer =',df[0:5]['inner'].sum()/df['outer'][0:5].sum())
#print('inner S /outer=')
#print(df['inner2']/df['outer'])
print('total inner S/total outer =',df['inner2'][0:5].sum()/df['outer'][0:5].sum())
CPU_linalg=df['cpulinsys'][0:5].sum()+df['cpuprec'][0:5].sum()
CPU_iterative=df['cpulinsys'][0:5].sum()
CPU_preprocess=df['cpuprec'][0:5].sum()
print("CPU ="+"{:1.2e}".format(CPU_linalg)+
      " Iterative = "+"{:.2f}".format(100*CPU_iterative/CPU_linalg)+
      "% -Preprocess ="+"{:.2f}".format(100*CPU_preprocess/CPU_linalg))

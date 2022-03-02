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
df['h']=1.0/df['Np']
df['dt']=1.0/df['N']

print('fixed H')
# fix h
for Np in [1,2,4]:
    df_dt=df[df['Np']==Np]
    x=df_dt['dt'].to_numpy()
    ymin=df_dt['min'].to_numpy()
    ymax=df_dt['max'].to_numpy()

    # fit y=b*x^a
    # slope=a
    # intercept=log(b)
    slope, intercept, r, p, se = linregress(np.log10(x), np.log10(ymin))
    alpha=slope
    beta=10**(intercept)
    #print('ymin~',"{:1.1e}".format(beta),'dt^',"{:1.1e}".format(alpha))
    #plot_fit(x,ymin,alpha,beta,'dt','ymin')

    slope, intercept, r, p, se = linregress(np.log10(x), np.log10(ymax))
    alpha=slope
    beta=10**(intercept)
    print('ymax~',"{:1.1e}".format(beta),'dt^',"{:1.1e}".format(alpha))
    #plot_fit(x,ymax,alpha,beta,'dt','ymax')

print('fixed dt')    

# fix dt
for N in [1,2,4]:
    df_h=df[df['N']==N]
    x=df_h['h'].to_numpy()
    ymin=df_h['min'].to_numpy()

    slope, intercept, r, p, se = linregress(np.log10(x), np.log10(ymin))
    alpha=slope
    beta=10**(intercept)
    print(x,ymin)
    print('ymin~',"{:1.1e}".format(beta),'h^',"{:1.1e}".format(alpha))
    #plot_fit(x,ymin,alpha,beta,'h','ymin')

for N in [1,2,4]:
    df_h=df[df['N']==N]
    x=df_h['h'].to_numpy()
    ymax=df_h['max'].to_numpy()

    
    slope, intercept, r, p, se = linregress(np.log10(x), np.log10(ymax))
    alpha=slope
    beta=10**(intercept)
    print(x,ymax)
    print('ymax~',"{:1.1e}".format(beta),'h^',"{:1.1e}".format(alpha))
    #plot_fit(x,ymax,alpha,beta,'h','ymax')

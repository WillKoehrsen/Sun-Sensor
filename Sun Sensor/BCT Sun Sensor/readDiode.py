import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv('diode_readings.csv')
df.set_index('time', inplace=True)

df.plot()
plt.show() 

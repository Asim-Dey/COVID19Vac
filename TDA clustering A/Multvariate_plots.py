%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_context('notebook')

import arviz

import seaborn as sns
############################# Implementing the Node2vec algorithm #################################

data = pd.read_csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv') # Load data
print(data)


data.columns = ["V"+str(i) for i in range(1, len(data.columns)+1)]  # rename column names to be similar to R naming convention
data.V1 = data.V1.astype(str)

print(data)

data.loc[:, "V9":"V14"]

pd.plotting.scatter_matrix(data.loc[:, "V14":"V16"], diagonal="kde")
plt.tight_layout()
plt.show()

  

#####################################################################################################

ax = data[["V10","V12","V14","V16","V19", "V21","V22","V23","V24","V25"]]
print(ax)

#Series_Complete_65PlusPop_Pct -- V16
#Percent.of.adults.with.a.bachelor.s.degree.or.higher.2015.19 --- V21
#covid_cases_per_100k--- V19


sns.pairplot(ax.loc[:, "V4":"V5"])

with sns.color_palette("Set3"):
sns.jointplot(data=ax, x="V5", y="V3", hue="V1")

with sns.color_palette("Set3"):
sns.jointplot(data=ax, x="V6", y="V3", hue="V1")


with sns.color_palette("Set3"):
sns.jointplot(data=ax, x="V7", y="V3", hue="V1")



sns.relplot(data=ax, x="V5", y="V3", hue="V1")


################ hue --clusters N 18 ###############

data_18A = pd.DataFrame({
    'Vaccinces Complete 65 Plus Pop, %': ax['V16'],
    'Vaccinces Complete 18 Plus Pop, %': ax['V14'],
    
    'Adults with a bachelors degree or higher, %': ax['V21'],
     'Covid cases (per 100k)': ax['V19'],
      'Median Household Income': ax['V23'],
      
     'Cluster18': ax['V24'],
      'Cluster65': ax['V25']
      
})
 



################ hue --clusters N 18 #################################################
sns.jointplot(data=data_18A, x='Median Household Income', 
              y='Vaccinces Complete 18 Plus Pop, %', hue="Cluster18",
               xlim=(0,125000),
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------
sns.jointplot(data=data_18A, x='Covid cases (per 100k)', 
              y='Vaccinces Complete 18 Plus Pop, %', hue="Cluster18",
               xlim=(0,600),
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------

sns.jointplot(data=data_18A, x='Adults with a bachelors degree or higher, %', 
              y='Vaccinces Complete 18 Plus Pop, %', hue="Cluster18",
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')


################ hue --clusters N 65 ##########################################################

sns.jointplot(data=data_18A, x='Median Household Income', 
              y='Vaccinces Complete 65 Plus Pop, %', hue="Cluster65",
              xlim=(0,125000),
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------
sns.jointplot(data=data_18A, x='Covid cases (per 100k)', 
              y='Vaccinces Complete 65 Plus Pop, %', hue="Cluster65",
               xlim=(0,600),
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------

sns.jointplot(data=data_18A, x='Covid cases (per 100k)', 
              y='Vaccinces Complete 65 Plus Pop, %', hue="Cluster65",
               xlim=(0,600),
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------

sns.jointplot(data=data_18A, x='Adults with a bachelors degree or higher, %', 
              y='Vaccinces Complete 65 Plus Pop, %', hue="Cluster65",
              palette=sns.color_palette("bright",4)) # bright

plt.legend(loc='lower right', title='Cluster')

#----------------------------------------------------------------------

sns.jointplot(data=data_18A, x="V6", y="V3", hue="V8",palette=sns.color_palette("bright",4))
plt.legend(loc='lower right', title='Cluster')
 
sns.jointplot(data=data_18A, x="V7", y="V3", hue="V8",palette=sns.color_palette("bright",4))
plt.legend(loc='lower right', title='Cluster')



################ hue --clusters N 65 ###############

sns.jointplot(data=ax, x="V5", y="V3", hue="V9",palette=sns.color_palette("bright",4))
plt.legend(loc='lower right', title='Cluster')
 
sns.jointplot(data=ax, x="V6", y="V3", hue="V9",palette=sns.color_palette("bright",4))
plt.legend(loc='lower right', title='Cluster')
 
sns.jointplot(data=ax, x="V7", y="V3", hue="V9",palette=sns.color_palette("bright",4))
plt.legend(loc='lower right', title='Cluster')

















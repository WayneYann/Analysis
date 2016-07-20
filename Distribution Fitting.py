# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 08:55:22 2016

@author: lansford
"""

import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib
import matplotlib.pyplot as plt
import os
vibration_file = os.path.expanduser('~/Box Sync/Synced_Files/Coding/Research/Analysis/Diatomic hydrogen and oxygen metals.csv')
Metal_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,0]
Data_Labels = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[0,1:]
Data = np.genfromtxt(vibration_file, delimiter=',')[1:,1:]
Metal_Info = np.genfromtxt(vibration_file, delimiter=',', dtype=str)[1:,0]
DataLabels = []
for i in range(0,len(Data_Labels)):
    DataLabels.append(Data_Labels[i])
    DataLabels.append(i)
matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
matplotlib.style.use('ggplot')

"""Getting Frequencies and residuals"""
vMOres =Data[:,2]  - Data[:,3]
vMHres = Data[:,7] - Data[:,8]
idxO = np.isfinite(vMOres)
vMOres = vMOres[idxO]
vMO = Data[:,3][idxO]
idxH = np.isfinite(vMHres)
vMHres = vMHres[idxH]
vMH = Data[:,8][idxH]

"""Getting Lengths and Residuals"""
rMOres = Data[:,4] - Data[:,5]
rMHres = Data[:,9] - Data[:,10]  
idxO = np.isfinite(rMOres)
rMOres = rMOres[idxO]
rMO = Data[:,5][idxO]
idxH = np.isfinite(rMHres)
rMHres = rMHres[idxH]
rMH = Data[:,10][idxH]

"""Getting percentage errors"""
pvMOres = vMOres/vMO
pvMHres = vMHres/vMH
prMOres = rMOres/rMO
prMHres = rMHres/rMH
prMO_Hres = np.concatenate((prMOres,prMHres))
pvMO_Hres = np.concatenate((pvMOres,pvMHres))

data = prMO_Hres*100

binnum = 9

# Create models from data
def best_fit_distribution(data, bins, ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    y, x = np.histogram(data, bins=binnum, normed=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Distributions to check
    DISTRIBUTIONS = [        
        st.alpha,st.anglit,st.arcsine,st.beta,st.betaprime,st.bradford,st.burr,st.cauchy,st.chi,st.chi2,st.cosine,
        st.dgamma,st.dweibull,st.erlang,st.expon,st.exponnorm,st.exponweib,st.exponpow,st.f,st.fatiguelife,st.fisk,
        st.foldcauchy,st.foldnorm,st.frechet_r,st.frechet_l,st.genlogistic,st.genpareto,st.gennorm,st.genexpon,
        st.genextreme,st.gausshyper,st.gamma,st.gengamma,st.genhalflogistic,st.gilbrat,st.gompertz,st.gumbel_r,
        st.gumbel_l,st.halfcauchy,st.halflogistic,st.halfnorm,st.halfgennorm,st.hypsecant,st.invgamma,st.invgauss,
        st.invweibull,st.johnsonsb,st.johnsonsu,st.ksone,st.kstwobign,st.laplace,st.levy,st.levy_l,st.levy_stable,
        st.logistic,st.loggamma,st.loglaplace,st.lognorm,st.lomax,st.maxwell,st.mielke,st.nakagami,st.ncx2,st.ncf,
        st.nct,st.norm,st.pareto,st.pearson3,st.powerlaw,st.powerlognorm,st.powernorm,st.rdist,st.reciprocal,
        st.rayleigh,st.rice,st.recipinvgauss,st.semicircular,st.t,st.triang,st.truncexpon,st.truncnorm,st.tukeylambda,
        st.uniform,st.vonmises,st.vonmises_line,st.wald,st.weibull_min,st.weibull_max,st.wrapcauchy
    ]

    # Best holders
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf

    # Estimate distribution parameters from data
    for distribution in DISTRIBUTIONS:

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # if axis pass in add to plot
                try:
                    if ax:
                        pd.Series(pdf, x).plot(ax=ax)
                    
                except Exception:
                    pass

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution
                    best_params = params
                    best_sse = sse

        except Exception:
            pass

    return (best_distribution.name, best_params)

def make_pdf(dist, params, size=1000):
    """Generate distributions's Propbability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf

# Load data from statsmodels datasets
#data = pd.Series(sm.datasets.elnino.load_pandas().data.set_index('YEAR').values.ravel())


data = pd.Series(data)

# Plot for comparison
plt.figure(figsize=(12,8))
ax = data.plot(kind='hist', bins=binnum, normed=False, alpha=0.5, color=plt.rcParams['axes.color_cycle'][1])
# Save plot limits
dataYLim = ax.get_ylim()

# Find best fit distribution
best_fit_name, best_fit_paramms = best_fit_distribution(data, binnum, ax)
best_dist = getattr(st, best_fit_name)
hist = np.histogram(data,bins=binnum)
# Update plots
ax.set_ylim(dataYLim)
ax.set_title(u'El Niño sea temp.\n All Fitted Distributions')
ax.set_xlabel(u'Temp (°C)')
ax.set_ylabel('Occurences')

# Make PDF
pdf = make_pdf(best_dist, best_fit_paramms)
pdf = pdf*max(hist[0])/max(pdf)
# Display
plt.figure(figsize=(14,10))
ax = pdf.plot(lw=2, label='PDF', legend=True)
data.plot(kind='hist', bins=binnum, normed=False, alpha=0.5, label='Data', legend=True, ax=ax)

param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
param_str = ', '.join(['{}={:0.3f}'.format(k,v) for k,v in zip(param_names, best_fit_paramms)])
dist_str = '{}({})'.format(best_fit_name, param_str)

ax.set_title(u'Histogram and Theoretical Distribution for Errors in Bond Length of Metal-Oxygen/Hydrogen Diatomics\n' + dist_str)
ax.set_xlabel('Error in Bond Length (%)')
ax.set_ylabel('Occurences')

"""Goodness of Fit test for distributions"""
print('Kolmogorov-Smirnov Goodness of fit test')
print('Oxygen length')
print('Hydrogen length')
print(st.kstest(prMHres*100,best_fit_name,best_fit_paramms))
print(st.kstest(prMOres*100,best_fit_name,best_fit_paramms))
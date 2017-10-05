import os
import sys
import pdb
import matplotlib
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
#import seaborn as sns
#sns.set_style("white")


def odds_ratios(file_name,output_file):
    plt.clf()
    fig = plt.figure()
    f = open(file_name)
    oddsratios = []
    ymin = []
    ymax = []
    count = 0
    tissues = []
    for line in f:
        line = line.rstrip()
        data = line.split()
        tissues.append(data[0])
        a = float(data[1])
        b = float(data[2])
        c = float(data[3])
        d = float(data[4])
        orat = (a/b)/(c/d)
        oddsratios.append(orat)
        log_or = np.log(orat)
        bounds = 1.96*np.sqrt((1.0/a) + (1.0/b) + (1.0/c) + (1.0/d))
        upper_bound = np.exp(log_or + bounds)
        lower_bound = np.exp(log_or - bounds)
        ymin.append(lower_bound)
        ymax.append(upper_bound)
    x = np.arange((len(oddsratios)))
    y = np.asarray(oddsratios)
    ytop = np.asarray(ymax)-y
    ybot = y - np.asarray(ymin)

    n = len(tissues)
    for i in range(n):
        plt.errorbar(x[i], y[i], yerr=([ybot[i]], [ytop[i]]), fmt = 'o')

    plt.xticks(range(len(oddsratios)),tissues,rotation='vertical',size = 14)
    plt.yticks(size=14)
    plt.plot([-1,len(oddsratios) ],[1,1],color = 'k')
    plt.subplots_adjust(bottom=.43)
    plt.xlabel('Pvalue Threshold',size=19)
    plt.ylabel('Enrichment',size=19)
    plt.title('DM Outlier Calling',size=19)

    fig.savefig(output_file)




input_file=sys.argv[1]
output_file = sys.argv[2]
odds_ratios(input_file,output_file)

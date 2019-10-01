#Name: RateSummary.py
#Date: 10 May 2018
#Creator: Joseph Hicks
#Python Version: 3.6
#Purpose: Sort and remove zeroes from actualRates (from a BSSVS analysis).
#         Summarize the rates with mean, median and 95%HPD.
#Input file: Remove all but "actual rates" from the rates log file. Then
#            rename the column names as "source%sink" but using find and
#            replace in excel to replace "." with "%".
########################################################################################################################
import pandas as pd
import numpy as np
import os

"""
This code was taken form the PyMC library https://github.com/pymc-devs/pymc
"""


def calc_min_interval(x, alpha):
    """Internal method to determine the minimum interval of a given width
    Assumes that x is sorted numpy array.
    """

    n = len(x)
    cred_mass = 1.0-alpha

    interval_idx_inc = int(np.floor(cred_mass*n))
    n_intervals = n - interval_idx_inc
    interval_width = x[interval_idx_inc:] - x[:n_intervals]

    if len(interval_width) == 0:
        raise ValueError('Too few elements for interval calculation')

    min_idx = np.argmin(interval_width)
    hdi_min = x[min_idx]
    hdi_max = x[min_idx+interval_idx_inc]
    return hdi_min, hdi_max


def hpd(x, alpha=0.05):
    """Calculate highest posterior density (HPD) of array for given alpha.
    The HPD is the minimum width Bayesian credible interval (BCI).
    :Arguments:
        x : Numpy array
        An array containing MCMC samples
        alpha : float
        Desired probability of type I error (defaults to 0.05)
    """

    # Make a copy of trace
    x = x.copy()
    # For multivariate node
    if x.ndim > 1:
        # Transpose first, then sort
        tx = np.transpose(x, list(range(x.ndim))[1:]+[0])
        dims = np.shape(tx)
        # Container list for intervals
        intervals = np.resize(0.0, dims[:-1]+(2,))

        for index in make_indices(dims[:-1]):
            try:
                index = tuple(index)
            except TypeError:
                pass

            # Sort trace
            sx = np.sort(tx[index])
            # Append to list
            intervals[index] = calc_min_interval(sx, alpha)
        # Transpose back before returning
        return np.array(intervals)
    else:
        # Sort univariate node
        sx = np.sort(x)
        return np.array(calc_min_interval(sx, alpha))

def BayesFactor(indicators, k):
    from math import exp, log
    prior = (log(2) + k - 1)/(k*(k-2)/2)
    posterior = indicators.mean()
    if posterior == 1:
        posterior = (len(indicators) - 1) / len(indicators)
    BF = (posterior/(1 - posterior))/(prior/(1 - prior))
    return BF


def TraitMatrix(list):
    matkey = {}
    length = len(list)
    size = int(length*(length-1)/2)
#    print(length)
#    print(size)
    index = 1
    for row in range(0, length-1):
        for col in range(row+1, length):
            matkey[index] = [list[row],list[col]]
#            print(index, matkey[index], row, col, sep='\t')
            index += 1
    for i in range(1, size+1):
        pair = matkey[i]
        matkey[i+size] = [pair[1],pair[0]]
#        print(i+size, matkey[i+size], sep='\t')
#    print(matkey)
#    x = matkey.get(1)
#    print(x[0])
    return(matkey)
#TraitMatrix(["a","b","c","d"])


def main():
    # Identify user selected directory path; only ".log" files from Beast v1.10 should be within the directory.
    dirname = input("Enter directory path: ")
    os.chdir(dirname)
    print(os.getcwd())
    # Create empty dataframe with desired columns


    # Loop through files in the directory
    for filename in os.listdir('./Logs/'):
        # Ignore the hidden '.DS_Store' file used in Mac OS and the output file
        if filename != '.DS_Store':
            output = pd.DataFrame(
                columns=['Epoch', 'Model', 'RateID', 'Source', 'Sink', 'Rate', 'CoefLowerHPD', 'CoefUpperHPD', 'pp', 'BF'])
            segmod = filename.split(sep='.')[0]
            splitfile = segmod.split(sep='_')
            epoch = splitfile[0]
            model = 'CountyCat'
            print(filename, epoch, model)
            for logs in os.listdir('./Logs/'):
                if epoch+'_' in logs:
                    logfile = logs
            print(epoch, model, logfile, sep='\t')

            Ratedata = pd.read_table('./Logs/'+logfile)

            with open('./Traits/APMV_CountyCat.txt') as f:
                traits = f.read().splitlines()
            K = len(traits)
            rateNum = K*(K-1)
            traitDict = TraitMatrix(traits)
 #           print(traitDict)
            # print(varNum)
            for i in range(rateNum):
                rateID = str(i + 1)
                rateIDint = i + 1
                indcol = model+'.indicators.'+epoch+rateID
                ratecol = model+'.actualRates.'+epoch+rateID
                indicators = Ratedata[indcol]

                if max(indicators) == 0:
                    interval = [0, 0]
                    median = 0
                else:
                    nozero = Ratedata[ratecol].replace(0, np.NaN)
                    nonan = nozero.dropna()
                    interval = hpd(nonan)
                    median = nonan.median()

                pp = indicators.mean()
                bf = BayesFactor(indicators, K)
                sousink = traitDict.get(rateIDint)
#                print(rateID)
#                print(sousink)
                output = output.append(pd.Series([epoch, model, rateID, sousink[0], sousink[1], median, interval[0], interval[1], pp, bf], index=output.columns), ignore_index=True)

            output.to_csv(dirname+'/RateSummary_'+epoch+'_'+model+'.txt',sep='\t')
    print("Done!")




main()

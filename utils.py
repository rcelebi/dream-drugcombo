import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
import math
from sklearn import metrics
import random
import numbers
import math

def primaryMetric(pred, obs):
    pred = pred.merge(obs, on=['COMBINATION_ID','CELL_LINE'])
    R = []
    w = []
    for g,df in pred.groupby('COMBINATION_ID'):
        x= df['PREDICTION']
        y= df['SYNERGY_SCORE']
        n = len(x)
        if  len(x) > 1:
            if np.sum(x==0) == n:
                r=0
            else:
                p = pearsonr(x=x,y=y)[0]
                r = math.sqrt(n-1)*p
            n = math.sqrt(n-1)
            R.append(r)
            w.append(n)
    return (sum(R)/sum(w))

def primaryMetric(pred):
    #pred = pred.merge(obs, on=['COMBINATION_ID','CELL_LINE'])
    R = []
    w = []
    for g,df in pred.groupby('COMBINATION_ID'):
        x= df['PREDICTION']
        y= df['SYNERGY_SCORE']
        n = len(x)
        if  len(x) > 1:
            if np.sum(x==0) == n:
                r=0
            else:
                p = pearsonr(x=x,y=y)[0]
                r = math.sqrt(n-1)*p
            n = math.sqrt(n-1)
            R.append(r)
            w.append(n)
    return (sum(R)/sum(w))


def multimetric_score(estimator, X_test, y_test, scorers):
    """Return a dict of score for multimetric scoring"""
    scores = {}
    for name, scorer in scorers.items():
        if y_test is None:
            score = scorer(estimator, X_test)
        else:
            score = scorer(estimator, X_test, y_test)

        if hasattr(score, 'item'):
            try:
                # e.g. unwrap memmapped scalars
                score = score.item()
            except ValueError:
                # non-scalar?
                pass
        scores[name] = score

        if not isinstance(score, numbers.Number):
            raise ValueError("scoring must return a number, got %s (%s) "
                             "instead. (scorer=%s)"
                             % (str(score), type(score), name))
    return scores

def get_scores(clf, X_new, y_new):

    scoring = ['r2','neg_mean_squared_error']
    scorers, multimetric = metrics.scorer._check_multimetric_scoring(clf, scoring=scoring)
    #print(scorers)
    scores = multimetric_score(clf, X_new, y_new, scorers)
    return scores

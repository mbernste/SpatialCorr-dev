import os
from os.path import join
import json
import numpy as np

SCHOT_DIR = './results_timing_rand_pairs/scHOT'
SPATIALDC_NO_MC_DIR = './results_timing_rand_pairs/SpatialDC_no_mc'
SPATIALDC_DIR = './results_timing_rand_pairs/SpatialDC'
SPATIALCORR_BR_NO_MC_DIR = './results_timing_rand_pairs/SpatialCorr_br_no_mc'
SPATIALCORR_BR_DIR = './results_timing_rand_pairs/SpatialCorr_br'

all_mins = []
for elem in os.listdir(SCHOT_DIR):
    fname = join(SCHOT_DIR, elem)
    with open(fname, 'r') as f:
        res = json.load(f)
        elapsed_time = res['elapsed_time']
        mins = elapsed_time/60
        all_mins.append(mins)
assert len(all_mins) == 20
print('scHOT average elapsed time: ', np.mean(all_mins))

all_mins = []
for elem in os.listdir(SPATIALDC_NO_MC_DIR):
    fname = join(SPATIALDC_NO_MC_DIR, elem)
    with open(fname, 'r') as f:
        res = json.load(f)
        elapsed_time = res['elapsed_time']
        mins = elapsed_time/60
        all_mins.append(mins)
assert len(all_mins) == 20
print()
print('SpatialCorr (no SMC) total elapsed time: ', np.sum(all_mins))
print('SpatialCorr (no SMC) average elapsed time: ', np.mean(all_mins))

all_mins = []
for elem in os.listdir(SPATIALDC_DIR):
    fname = join(SPATIALDC_DIR, elem)
    with open(fname, 'r') as f:
        res = json.load(f)
        elapsed_time = res['elapsed_time']
        mins = elapsed_time/60
        all_mins.append(mins)
assert len(all_mins) == 20
print()
print('SpatialCorr total elapsed time: ', np.sum(all_mins))
print('SpatialCorr average elapsed time: ', np.mean(all_mins))

all_mins = []
for elem in os.listdir(SPATIALCORR_BR_NO_MC_DIR):
    fname = join(SPATIALCORR_BR_NO_MC_DIR, elem)
    with open(fname, 'r') as f:
        res = json.load(f)
        elapsed_time = res['elapsed_time']
        mins = elapsed_time/60
        all_mins.append(mins)
assert len(all_mins) == 20
print()
print('SpatialCorr (BR-test; no SMC) total elapsed time: ', np.sum(all_mins))
print('SpatialCorr (BR-test; no SMC) average elapsed time: ', np.mean(all_mins))

all_mins = []
for elem in os.listdir(SPATIALCORR_BR_DIR):
    fname = join(SPATIALCORR_BR_DIR, elem)
    with open(fname, 'r') as f:
        res = json.load(f)
        elapsed_time = res['elapsed_time']
        mins = elapsed_time/60
        all_mins.append(mins)
assert len(all_mins) == 20
print()
print('SpatialCorr (BR-test) total elapsed time: ', np.sum(all_mins))
print('SpatialCorr (BR-test) average elapsed time: ', np.mean(all_mins))


###################
# Noah Burget
# Last updated: 12/7/24
# Driver script for the 1-Dimensional Portion of the loop extrusion simulations
# Originally written as a jupyter notebook, but converted to .py script for simplicity
###################

### 1-D loop extrusion simulation
from extruder import Extruder
import numpy as np
from pathlib import Path
import seaborn as sns
import h5py
import sys

def main():
    ###############################################################
    ### Translocating and writing trajectories using custom extruder class
    # Basic parameters for simulation run
    ###############################################################
    RUN_NAME = "MYC_Granta519_KO_EBF1Loading" #### Define a name for this simulation run, sowe can save parameters
    N1_pol = 900 # Size of 1 system
    front_buffer = 10 # 'buffer' zone for poltting purposes
    end_buffer = 10
    N1 = N1_pol + front_buffer + end_buffer
    M = 1 # No. of systems
    N = N1*M # Total size of system, in momomers
    occupied = np.zeros(N) # List to tell if current monomer is occupied by an extruder
    occupied[0] = 1 
    occupied[-1] = 1
    steps = 50000 # Timesteps for 1D sim.
    num_chunks = 50 # No. of chunks to write trajectories in
    LIFETIME = 800 # Cohesin lifetime
    LIFETIME_STALLED = LIFETIME // 10 # Cohesin lifetime when stalled

    ### Blockers (i.e. CTCF) - {pos. : prob.}
    left_blockers_capture = {}
    right_blockers_capture = {}
    left_blockers_release = {}
    right_blockers_release = {}

    TIER1_BLOCK = [0.99, 0.001]
    TIER2_BLOCK = [0.75, 0.001]
    TIER3_BLOCK = [0.5, 0.015] ### Increasing probability of release
    TIERB3_BLOCK = [0.45, 0.03]
    TIER4_BLOCK = [0.1, 0.1]

    #####################################################################
    ### Defining 'blocking regions' - consecutive monomers that can block
    # Each call below defines the blocking region as (low, high), EXCLUSIVE on the upper bound
    #####################################################################
    E1_blockingRegion_1 = np.arange(181+front_buffer, 185+front_buffer) # Boundary, not EBF1 involved
    E1_blockingRegion_2 = np.arange(210+front_buffer, 220+front_buffer) # E1, EBF1 involved (50, STRONG)
    E2_blockingRegion = np.arange(304+front_buffer, 307+front_buffer) # E2, EBF1 involved (15, MED)
    B3_blockingRegion_forward = np.arange(568+front_buffer, 570+front_buffer) # B3+, EBF1 involved (8, WEAK)
    B3_blockingRegion_reverse = np.arange(576+front_buffer, 578+front_buffer) # B3-, not EBF1 involved
    MYC_blockingRegion = np.arange(737+front_buffer, 742+front_buffer) # MYC, not EBF1 involved

    ###########################################################################
    ### Any EBF1-associated blockers below should have '_EBF1' as the suffix in the keys of the dictionary below
    ###########################################################################
    blockingRegions = {'E1_1':E1_blockingRegion_1, 'E1_2':E1_blockingRegion_2, 'E2':E2_blockingRegion, 'B3+': B3_blockingRegion_forward, 'B3-':B3_blockingRegion_reverse, 'MYC':MYC_blockingRegion}

    ### Assign each blocking region its capture and release probabilities (0<=x<=1)
    for br in blockingRegions:
        # Boundary blocker
        if br == 'E1_1': # Will never be EBF1 involved
            cap = 0.99 # probability of a monomer capturing a LEF
            rel = 0.001 # probability of a monomer releasing a LEF, typically much lower than capture prob. 
        # E1 blocker
        if 'E1_2' in br: # EBF1 involved
            cap = 0.5
            rel = 0.03
        # E2 blocker
        elif 'E2' in br: # EBF1 involved
            cap = 0.35
            rel = 0.03
        # B3 blocker
        elif 'B3-' in br: # EBF1 involved
            cap = 0.35
            rel = 0.03
        elif 'B3+' in br:
            cap = 0.2
            rel = 0.03
        # MYC promoter blocker
        elif br == 'MYC': # Unchanged
            cap = 0.75
            rel = 0.001
        # Now, also according to the dictionary key, we assign the directional blockers with the previously defined probabilities
        if 'EBF1' in br: # IF EBF1 blocker, bidirectional (in the EBF1 loading system, and the perturbed EBF1 blocking system, this will never happen)
            for loc in blockingRegions[br]:
                left_blockers_capture[loc] = cap
                left_blockers_release[loc] = rel
                right_blockers_capture[loc] = cap
                right_blockers_release[loc] = rel
        # The below two conditions will only be evaluated if the blocker is NOT EBF1-associated
        elif br == 'MYC' or br == 'B3-': # The two negatively-oriented blockers
            for loc in blockingRegions[br]:
                right_blockers_capture[loc] = cap
                right_blockers_release[loc] = rel
        else: # Just for E1_1
            for loc in blockingRegions[br]:
                left_blockers_capture[loc] = cap
                left_blockers_release[loc] = rel

    # Manually assigning blockers in dict
    # These are the old monomer blockers.
    #left_blockers_capture[180] = TIER1_BLOCK[0] # BOUNDARY -- SHOULD NOT CHANGE **
    #left_blockers_release[180] = TIER1_BLOCK[1]
    #left_blockers_capture[210] = TIER1_BLOCK[0] # E1.2 -- BEGINNING OF PROBE 8 (BiD) -- DECREASES
    #left_blockers_release[210] = TIER1_BLOCK[1]
    #left_blockers_capture[315] = TIER3_BLOCK[0] # E2 -- EBF1 -- CENTER OF PROBE 11 (BiD) - DRAMATICALLY DECREASES
    #left_blockers_release[315] = TIER3_BLOCK[1]
    left_blockers_capture[405+front_buffer] = TIER4_BLOCK[0] # B1 -- CENTER OF PROBE 14 -- SHOULD NOT CHANGE **
    left_blockers_release[405+front_buffer] = TIER4_BLOCK[1]
    #left_blockers_capture[555] = TIERB3_BLOCK[0] # B3 -- EBF1 -- CENTER OF PROBE 19 (BiD) -- DECREASES
    #left_blockers_release[555] = TIERB3_BLOCK[1]


    #right_blockers_capture[210] = TIER1_BLOCK[0] # E1.2 -- BEGINNING OF PROBE 8 (BiD) -- DECREASES
    #right_blockers_release[210] = TIER1_BLOCK[1]
    #right_blockers_capture[315] = TIER3_BLOCK[0] # E2 -- EBF1 -- CENTER OF PROBE 11 (BiD) -- DRAMATICALLY DECREASES
    #right_blockers_release[315] = TIER3_BLOCK[1]
    right_blockers_capture[495+front_buffer] = TIER4_BLOCK[0] # B2 -- CENTER OF PROBE 17 -- SHOULD NOT CHANGE **
    right_blockers_release[495+front_buffer] = TIER4_BLOCK[1]
    #right_blockers_capture[555] = TIERB3_BLOCK[0] # B3 -- EBF1 -- CENTER OF PROBE 19 (BiD) -- DECREASES
    #right_blockers_release[555] = TIERB3_BLOCK[1]
    #right_blockers_capture[735] = TIER1_BLOCK[0] # MYC promoter -- CENTER OF PROBE 25 -- SHOULD NOT CHANGE **
    #right_blockers_release[735] = TIER1_BLOCK[1]

    ########################################
    ### Now we define the cohesion LOADING regions
    ########################################
    E1_loading_1 = [186+front_buffer, 210+end_buffer] # Region 1
    E1_loading_2 = [221+front_buffer, 241+end_buffer] # Region 2
    E2_loading = [310+front_buffer, 318+end_buffer] # Region 3
    B3_loading_1 = [553+front_buffer, 570+end_buffer] # Region 4
    MYC_loading = [697+front_buffer, 737+end_buffer] # Region 5
    wholePolymer_loading = [0+front_buffer, N1-(end_buffer)-1] # Entire polymer

    ### List of loading regions for EBF1 loading
    loading_regions = [E1_loading_1, E1_loading_2, E2_loading, B3_loading_1, MYC_loading, wholePolymer_loading]

    ### List of loading regions for EBF1 blocking (just the entire polymer)
    #loading_regions = [wholePolymer_loading]

    ### Number of cohesins per loading region (each index is associated with the corresponding index of the loading_regions list)
    loading_region_freqs = [1,1,1,1,1,4]
    LEFNum = np.sum(loading_region_freqs)

    ### Beyond this point, you need not change any variable values.
    ###############################################################

    # Generate (initial) loading region spots for the loading region LEFs and list of loading site for each LEF
    LOADING_SPOTS = []
    REGIONS_INDEX = []
    for i, region in enumerate(loading_regions):
        print(i, region)
        for j in range(loading_region_freqs[i]):
            while True:
                spot = np.random.randint(low=region[0], high=region[1])
                if spot not in LOADING_SPOTS and spot+1 not in LOADING_SPOTS and spot-1 not in LOADING_SPOTS:
                    print(spot)
                    break
            LOADING_SPOTS.append(spot)
            REGIONS_INDEX.append(loading_regions[i])
    print(REGIONS_INDEX)

    EXTRUDERS = []
    targeted_idx = 0
    for i in range(len(LOADING_SPOTS)):
        #leg = np.random.randint(N)
        #print('Loading extruder at {}'.format(LOADING_SPOTS[i]))
        leg = LOADING_SPOTS[i]
        print('Loading extruder at {}, its loading region is {}'.format(leg, REGIONS_INDEX[i]))
        EXTRUDERS.append(Extruder(
                extruder_index = i,
                leg1 = leg,
                leg2 = leg+1,
                left_blockers_capture = left_blockers_capture,
                right_blockers_capture = right_blockers_capture,
                left_blockers_release = left_blockers_release,
                right_blockers_release = right_blockers_release,
                extrusion_occupancy = occupied,
                loading_region = REGIONS_INDEX[i],
                lifetime = LIFETIME,
                lifetime_stalled = LIFETIME_STALLED)
        )

    ### Write parameters to text file
    with open('{}_params.txt'.format(RUN_NAME),'w+') as pf:
        pf.write("N: {}\n1D Steps: {}\nTotal LEF: {}\nLoading Regions: {} LEFs per loading region: {}\n Lifetime: {} Lifetime stalled: {}\nBlocking Regions: {}\nLeft cap: {}, Left rel: {}\nRight cap: {}, Right rel: {}".format(
                                                                                                N1,steps,LEFNum,
                                                                                                REGIONS_INDEX,loading_region_freqs,
                                                                                                LIFETIME,LIFETIME_STALLED,blockingRegions,left_blockers_capture,left_blockers_release,
                                                                                                right_blockers_capture,right_blockers_release))
    outf = "trajectory/LEFPositions.h5"
    p = Path(outf)
    if p.exists():
        p.unlink()
    with h5py.File(outf, mode='w') as f:
        dset = f.create_dataset("positions", 
                shape=(steps, LEFNum, 2), 
                dtype=np.int32, 
                compression="gzip")
        bins = np.linspace(0, steps, num_chunks, dtype=int)
        for st,end in zip(bins[:-1], bins[1:]): # Loop through bins
            cur = []
            for i in range(st,end): # For bin in bins 
                positions = [(extruder.leg1.pos, extruder.leg2.pos) for extruder in EXTRUDERS] # Get both leg positions for all extruders
                cur.append(positions)
                for extruder in EXTRUDERS:
                    occupied = extruder.translocate(occupied) # Translocate extruder
            cur = np.array(cur)
            dset[st:end] = np.array(cur)
        f.attrs["N"] = N
        f.attrs["LEFNum"] = LEFNum
    del EXTRUDERS

if __name__ == '__main__':
    main()


#          _
#      .__(.)< (MEOW)
#       \___)

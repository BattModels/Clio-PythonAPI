from ElectrolyteComposition import ElectrolyteComposition
import pandas as pd
import numpy as np
import timeit
import itertools

##Utility function definitions

def fixed_length_partitions(n,L):
    """
    Integer partitions of n into L parts, in colex order.
    The algorithm follows Knuth v4 fasc3 p38 in rough outline;
    Knuth credits it to Hindenburg, 1779.
    """
    
    # guard against special cases
    if L == 0:
        if n == 0:
            yield []
        return
    if L == 1:
        if n > 0:
            yield [n]
        return
    if n < L:
        return

    partition = [n - L + 1] + (L-1)*[1]
    while True:
        yield partition
        if partition[0] - 1 > partition[1]:
            partition[0] -= 1
            partition[1] += 1
            continue
        j = 2
        s = partition[0] + partition[1] - 1
        while j < L and partition[j] >= partition[0] - 1:
            s += partition[j]
            j += 1
        if j >= L:
            return
        partition[j] = x = partition[j] + 1
        j -= 1
        while j > 0:
            partition[j] = x
            s -= x
            j -= 1
        partition[0] = s

def make_parts(dims=4,total=19):
    parts=[]
    for m in range(1,dims+1):
        for p in fixed_length_partitions(total,m):
            p += [0] * (dims-len(p))
            p=[int(x) for x in p]
            for pi in set(itertools.permutations(p)):
                parts.append(pi)
    return np.array(parts)


## Settings for solvents, salts, feeder solutions, volume step size, and total mix volume 
#for creating
all_solvents = ["EC","DMC"]
salts = {"LiTFSI":287.1} #(key,value) : (salt specie, molar mass (g/mol) )
all_salts=list(salts.keys())
densities_all={
	"DMC|100|LiTFSI|2":1.321,
	"DMC_EC|50_50|LiTFSI|2":1.436,
	"DMC|100":1.071,
	"DMC_EC|50_50":1.192,
}
volume_step=200
total=2200
fn="dmc_ec_tfsi_grid.h5"

#STEP 1: make volume grid
#Generate grid of all combinations of components in solvent_disc steps up to mix_volume
components=list(densities_all.keys())
df=pd.DataFrame()
start_time = timeit.default_timer()
grid=make_parts(dims=len(components),total=total/volume_step)
for i in range(grid.shape[1]):
  df[components[i]]=volume_step*grid[:,i]
elapsed = timeit.default_timer() - start_time
print("Step 1 {:.3f}".format(elapsed))

#STEP 2: calculate resulting electrolytes
#make source vectors
start_time = timeit.default_timer()
source_electrolytes = [ElectrolyteComposition.by_CompositionID(components[i])for i in range(len(components))]
source_solvent_mfs = np.array([[el.solvents[solvent]/el.solvent_precision  if solvent in el.solvents else 0 for solvent in all_solvents] for el in source_electrolytes])
source_salt_molarities = np.array([[el.salts[salt]  if salt in el.salts else 0 for salt in all_salts] for el in source_electrolytes])
salt_mms = np.array(list(salts.values()))
source_salt_mass = source_salt_molarities@salt_mms #assume 1kg solvent
source_mf_salt= source_salt_mass/(source_salt_mass+1000)
source_mf_solvent = np.ones(len(source_mf_salt))-source_mf_salt
#Calculate resulting electrolytes (vectorized)
volumes = df.to_numpy()
densities = np.array(list(densities_all.values()))
total_mass_per_source = 1/1000*np.multiply(volumes,densities) #total mass per source (grams)
total_salt_mass_per_source = np.multiply(total_mass_per_source,source_mf_salt)
total_solvent_mass_per_source=np.multiply(total_mass_per_source,source_mf_solvent)
total_salt_masses = np.sum(total_salt_mass_per_source,axis=1)
total_solvent_masses  = np.sum(total_solvent_mass_per_source,axis=1)
salt_molarities = np.divide(total_salt_masses/salt_mms,total_solvent_masses/1000)
solvent_mass_fractions=np.divide(total_solvent_mass_per_source@source_solvent_mfs,np.vstack([total_solvent_masses]*len(all_solvents)).T)
elapsed = timeit.default_timer() - start_time
print("Step 2 {:.3f}".format(elapsed))

#Write to dataframe
df["LiTFSI_m"] = salt_molarities
df["EC_mf"] = solvent_mass_fractions[:,0]
df["DMC_ratio"] = solvent_mass_fractions[:,1]/(np.ones(solvent_mass_fractions[:,0].shape) - solvent_mass_fractions[:,0])
df.to_hdf(fn,key='df')


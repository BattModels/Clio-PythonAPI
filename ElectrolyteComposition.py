import subprocess as sp
import pandas as pd
import re
import collections
from collections import OrderedDict
import pdb
import datetime
import json
import os
# from ElectrolyteComposition import data_dir
data_dir = "/Users/adarsh.r.dave/repos/natcomm_repo/data"

delim1="|"
delim2="_"
default_salt_decimals=2
default_solvent_precision=100
vals2str = lambda ls : [str(x) for x in ls]

class ElectrolyteComposition:
    # * TODOs: 
    #     * init from JSON specifier
    def __init__(self,solvents=None,salts=None,specified_from=None,solvent_DB=None,salt_DB=None,CompositionID=None,
                solvent_precision=default_solvent_precision, salt_decimals=default_salt_decimals):
        self.solvents=dict() if solvents==None else solvents
        self.salts=dict() if salts==None else salts
        self.specified_from="" if specified_from==None else specified_from #will contain info from classmethod used to create
        self.solvent_DB=self.load_solvent_DB() if solvent_DB==None else solvent_DB
        self.salt_DB=self.load_salt_DB() if salt_DB==None else salt_DB
        self.CompositionID="" if CompositionID==None else CompositionID
        self.solvent_precision=solvent_precision
        self.salt_decimals=salt_decimals
        #Get date of composition made
        self.date=datetime.datetime.now().strftime("%m/%d/%Y")
        if len(salts)>1:
            raise NotImplementedError("Binary salts not implemented, yet - Ady")
    def __repr__(self):
        return self.CompositionID
    def dump_info(self):
        solvent_DB=self.solvent_DB
        salt_DB=self.salt_DB
        filt_solvent_info=solvent_DB[solvent_DB.name.isin(self.solvents.keys())].to_json(orient="records")
        filt_salt_info=salt_DB[salt_DB.name.isin(self.salts.keys())].to_json(orient="records")
        #r_d={"chemicals":{"solvents":filt_solvent_info,"salts":filt_salt_info}}
        return {"solvents":filt_solvent_info,"salts":filt_salt_info}
    def name_composition(self):
        rep=self.CompositionID.replace("_","")
        return rep.replace("|","")
    def to_solution_volume(self):
        raise NotImplementedError()
        return
    def to_mole_fraction(self):
        missing_solvents=[solvent for solvent in self.solvents.keys() if solvent not in self.solvent_DB.name.to_list()]
        assert (len(missing_solvents)==0), f"Solvents {missing_solvents} not present in {data_dir}/solvent_DB.csv"
        missing_salts=[salt for salt in self.salts.keys() if salt not in self.salt_DB.name.to_list()]
        assert (len(missing_salts)==0), f"Salts {missing_salts} not present in {data_dir}/salt_DB.csv"
        mol_dict={}
        #calculate total mols across solvents and salts
        for solvent,mf_dec in self.solvents.items():
            mass_fraction=mf_dec/self.solvent_precision #convert to true mass fraction (in kg) from imprecise representation
            #kg/mol = g/mol * 1 kg / 1000 g
            mol_dict[solvent] = mass_fraction/(self.solvent_DB.loc[self.solvent_DB.name==solvent]["molar_mass"].values[0]/1000)
        for salt,molality in self.salts.items():
            mol_dict[salt] = molality
        return mol_dict        

    @staticmethod
    def cid_to_parsable(cid):
        rep=cid.replace("_","")
        return rep.replace("|","")
    @staticmethod
    def normalize_solvent_dictionary(solvents,solvent_precision):
        #pdb.set_trace()
        total=float(sum(solvents.values()))
        _solvents={solvent:int(round(solvents[solvent]/total*solvent_precision)) for i,solvent in enumerate(solvents.keys())}#round
        _solvents_nonzero={solvent:_solvents[solvent] for solvent in _solvents.keys() if _solvents[solvent]!=0}#filter
        ordered_solvents=OrderedDict(sorted(_solvents_nonzero.items(), key=lambda tup: tup[0]))
        return ordered_solvents
    @staticmethod
    def normalize_salt_dictionary(salts,salt_decimal):
        _salts={salt:round(salts[salt],salt_decimal) for salt in salts.keys()} #round
        _salts_nonzero={salt:_salts[salt] for salt in _salts.keys() if _salts[salt]!=0.0} #filter
        ordered_salts=OrderedDict(sorted(_salts.items(), key=lambda tup: tup[0]))
        return ordered_salts
    @staticmethod
    def load_solvent_DB(filename="solventDB.csv"):
        path=os.path.join(data_dir,filename)
        return pd.read_csv(path)
    @staticmethod
    def load_salt_DB(filename="saltDB.csv"):
        path=os.path.join(data_dir,filename)
        return pd.read_csv(path)
    @staticmethod
    def dicts_to_CompositionID(solvents={},salts={},solvent_precision=default_solvent_precision,salt_decimals=default_salt_decimals):
        #Filter
        solvents_normalized=ElectrolyteComposition.normalize_solvent_dictionary(solvents,solvent_precision)
        if len(salts)!=0:
            salts_normalized=ElectrolyteComposition.normalize_salt_dictionary(salts,salt_decimals)
            return delim1.join([delim2.join(x) for x in [solvents_normalized.keys(),vals2str(solvents_normalized.values()),salts_normalized.keys(),vals2str(salts_normalized.values())]])
        else:
            return delim1.join([delim2.join(x) for x in [solvents_normalized.keys(),vals2str(solvents_normalized.values())]])
    @staticmethod
    def CompositionID_to_dicts(CompositionID):
        ls=CompositionID.split(delim1)
        solvent_names=ls[0].split(delim2)
        solvent_mfs=[i for i in ls[1].split(delim2)] #normalize?
        assert len(solvent_names)==len(solvent_mfs), "CompositionID is invalid, different lengths for solvent_names vs solvent_mfs"
        assert 0 not in solvent_mfs, "Zeros not allowed in defining composition" #still caught by normalizing
        solvent_mfs_precisions=list(set([len(i) if len(i)>1 else 2 for i in solvent_mfs]))
        assert len(solvent_mfs_precisions)==1, "Length (precision) of solvent mass fractions must be identical: {}".format(solvent_mfs_precisions)
        #pdb.set_trace()
        solvent_precision=int(10**int(solvent_mfs_precisions[0]))
        solvent_mfs=[float(i) for i in solvent_mfs]
        solvents=ElectrolyteComposition.normalize_solvent_dictionary({solvent_names[i]:solvent_mfs[i] for i in range(len(solvent_names))},solvent_precision)
        salt_decimals=default_salt_decimals #temporary!
        if len(ls)>2:
            assert len(ls)==4, "If salts are added, must define molality"
            salt_names=ls[2].split(delim2)
            molality=[float(i) for i in ls[3].split(delim2)]
            assert len(salt_names)==len(molality), "CompositionID is invalid, different lengths for salt_names vs molality"
            salts=ElectrolyteComposition.normalize_salt_dictionary({salt_names[i]:molality[i] for i in range(len(salt_names))},salt_decimals)
        else:
            salts={}
        return {"solvents":solvents,"salts":salts,"solvent_precision":solvent_precision,"salt_decimals":salt_decimals}
    @classmethod
    def by_CompositionID(cls, CompositionID):
        dicts=cls.CompositionID_to_dicts(CompositionID)
        return cls(**dicts,CompositionID=CompositionID,specified_from=json.dumps({"CompositionID":CompositionID}))
    @classmethod
    def by_mass(cls,solvents={},salts={}): #solvent mass, salt mass
        raise NotImplementedError
    @classmethod
    def by_mass_fraction_and_molality(cls,solvents={},salts={},solvent_precision=default_solvent_precision,salt_decimals=default_salt_decimals,
        specified_from=None): 
        #solvent mass fraction, salt molality
        solvents_orig=solvents.copy()
        salts_orig=salts.copy()
        solvents_normalized=cls.normalize_solvent_dictionary(solvents,solvent_precision)
        if len(salts)!=0:
            salts_normalized=cls.normalize_salt_dictionary(salts,salt_decimals)
        else:
            salts_normalized=salts
        cid=cls.dicts_to_CompositionID(solvents=solvents_normalized,salts=salts_normalized,solvent_precision=solvent_precision,salt_decimals=salt_decimals)
        d={"solvents":solvents_normalized,"salts":salts_normalized,"CompositionID":cid,"solvent_precision":solvent_precision,"salt_decimals":salt_decimals}
        if specified_from==None:
            specified_from=json.dumps({"by_mass_fraction_and_molality":{"solvents":solvents_orig,"salts":salts_orig}})
        return cls(**d,specified_from=specified_from)
    @classmethod
    def by_mole_fraction(cls,solvents={},salts={},solvent_precision=default_solvent_precision,salt_decimals=default_salt_decimals):
        #Ensure data is present in DB
        solvent_DB=cls.load_solvent_DB()
        salt_DB=cls.load_salt_DB()
        missing_solvents=[solvent for solvent in solvents.keys() if solvent not in solvent_DB.name.to_list()]
        assert (len(missing_solvents)==0), f"Solvents {missing_solvents} not present solvent_DB.csv"
        missing_salts=[salt for salt in salts.keys() if salt not in salt_DB.name.to_list()]
        assert (len(missing_salts)==0), f"Salts {missing_salts} not present in salt_DB.csv"

        #Ensure mole fractions appropriately normalized such that sum is 1
        total_moles=sum(solvents.values())+sum(salts.values())
        solvents={solvent:solvents[solvent]/total_moles for solvent in solvents.keys()}
        salts={salt:salts[salt]/total_moles for salt in salts.keys()}

        #calculate solvent masses, all in kg
        solvent_masses= {solvent:solvents[solvent]*(solvent_DB.loc[solvent_DB.name==solvent]["molar_mass"].values[0]/1000.) for solvent in solvents.keys()}
        total_solvent_mass=sum(solvent_masses.values())

        salt_free_solvent_mass_fractions = {solvent:solvent_masses[solvent]/total_solvent_mass for solvent in solvent_masses.keys()}
        salt_molalities={salt:salts[salt]/total_solvent_mass if salts[salt]!=0 else 0 for salt in salts.keys()}

        return cls.by_mass_fraction_and_molality(solvents=salt_free_solvent_mass_fractions,salts=salt_molalities,
            specified_from=json.dumps({"by_mole_fraction":{"solvents":solvents,"salts":salts}}))

    @classmethod
    def by_solution_volume(cls,volumes={},densities={},solvent_precision=default_solvent_precision,salt_decimals=default_salt_decimals):
        solvent_DB=cls.load_solvent_DB()
        salt_DB=cls.load_salt_DB()
        volumes={k:int(v) for k,v in volumes.items()}
        densities={k:float(v) for k,v in densities.items()}
        specified_from=json.dumps({"by_solution_volume":{"volumes":volumes.copy(),"densities":densities.copy()}})
        solvents={} #mass fraction
        salts={} #molality
        solvents_mass={}
        salts_moles={}
        assert set(volumes.keys())==set(densities.keys()), "Same keys must be in each of volumes and densities"
        total_dose_masses={solution:volumes[solution]/1000*densities[solution] for solution in volumes.keys()} #this is in grams
        for solution in total_dose_masses.keys():
            solution_comp=cls.CompositionID_to_dicts(solution) #solution_comp["solvents"] is m.f.; '' ["salts"] is molal
            source_solvent_precision=int(solution_comp["solvent_precision"])

            #BOTTLE LEVEL
            solution_total_salt_mass=0
            if len(solution_comp["salts"])!=0:
                for salt in solution_comp["salts"].keys():
                    assert salt in salt_DB.name.values, "Salt proposed that is not in salt_DB, please check! - {}".format(salt)
                    mm=float(salt_DB[salt_DB.name==salt]["molar mass"].iloc[0])
                    m=solution_comp["salts"][salt]
                    solution_total_salt_mass += mm*m #g/mol * molality of single salt = mass of this salt in bottle
            solution_salt_mass_fraction = solution_total_salt_mass / (solution_total_salt_mass+1000)
            solution_solvent_mass_fraction = 1-solution_salt_mass_fraction

            #DOSE LEVEL
            dose_total_solvent_mass=solution_solvent_mass_fraction*total_dose_masses[solution]
            for solvent in solution_comp["solvents"].keys():
                if solvent not in solvents_mass:
                    solvents_mass[solvent]=dose_total_solvent_mass*solution_comp["solvents"][solvent]/source_solvent_precision
                else:
                    solvents_mass[solvent]+=dose_total_solvent_mass*solution_comp["solvents"][solvent]/source_solvent_precision
            if len(solution_comp["salts"])!=0:
                for salt in solution_comp["salts"].keys():
                    m=solution_comp["salts"][salt]
                    if salt not in salts_moles:
                        salts_moles[salt]=m*dose_total_solvent_mass/1000
                    else:
                        salts_moles[salt]+=m*dose_total_solvent_mass/1000
        
        #EVERYTHING HAS BEEN TOTALED
        solvents=cls.normalize_solvent_dictionary(solvents_mass,solvent_precision)
        salts=cls.normalize_salt_dictionary({salt:salts_moles[salt]/(sum(list(solvents_mass.values())))*1000 for salt in salts_moles},salt_decimals)
        cid=cls.dicts_to_CompositionID(solvents=solvents,salts=salts,solvent_precision=solvent_precision,salt_decimals=salt_decimals)
        #total salts into moles each, total solvents into mass each, give molality.
        d={"solvents":solvents,"salts":salts,"CompositionID":cid,"solvent_precision":solvent_precision,"salt_decimals":salt_decimals}
        return cls(**d,specified_from=specified_from)
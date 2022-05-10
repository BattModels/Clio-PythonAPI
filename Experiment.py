import requests
import pandas as pd
import json
import string
import random
import datetime
from decimal import Decimal
import pdb

class InventorySpeciesError(Exception):
    pass

class InventoryVolumeError(Exception):
    pass

class Clio:
    def __init__(self,run_type,person,mix_volume=1900,inventory=None,retain_inventory=None,verbose=True,
                wash_solvent="MECN",open_valve_name="OPEN2",waste_valve_name="WASTE",wash_volume=6000): #vols now in uL
        #Core attributes
        self.verbose=verbose #printing Otto's initialization message with inventory
        self.floor_volume=2000 #least amount of vol in a bottle before calling it empty
        self.preprime_floor=1000 #a way to use up finished solution
        self.precision=0 #How many decimal points for rounding out floating point errors?
        self.priming_volume=800 #uL volume for priming step <<< INPUT HERE
        self.timeout=3600 #seconds

        #URLs
        self.base_url = ""
        self.experiment_url='RUN'
        self.experiment_url_fast="RUNfast_derated2"
        self.retain_url="RUN_retain"
        self.wash_url='RINSE'

        #Volumes
        self.mix_volume=mix_volume #total volume

        #Wash params
        self.wash_solvent=wash_solvent
        self.open_valve_name=open_valve_name
        self.waste_valve_name=waste_valve_name

        self.wash_volume=wash_volume

        #Matr.io required items
        self.ingest_URL=''
        self.api_key=''
        #self.headers={'x-api-key': self.api_key}
        self.matrioproject="MUSE_CLIO"

        #various utility strings for processing
        ##NOTHING FOR NOW

        #Test connection
        ##NOTHING FOR NOW
        
        #Generate RunID and RunType
        length=6
        date=datetime.datetime.now().strftime("%m%d%Y")
        time=datetime.datetime.now().strftime("%H:%M:%S")
        self.start_date=date
        self.start_time=time
        s=person+date
        for i in range(length):
            s+=random.choice(string.ascii_lowercase)
        self.run_id=s
        self.run_type=run_type

        #Load in inventory - can pass string filepath, or (default) just read from labview
        if inventory:
            #Read inventory from local file
            self.inventory=pd.read_csv(inventory)#.set_index(['chemical'])
        else:
            #Get from Otto
            self.inventory = pd.read_csv(self.base_url+"inventory.csv")

        #Load in retain_inventory - can pass string filepath, or (default) just read from labview
        if retain_inventory:
            #Read inventory from local file
            self.retain_inventory=pd.read_csv(retain_inventory)#.set_index(['chemical'])
        else:
            #Get from Otto
            self.retain_inventory = pd.read_csv(self.base_url+"retain_inventory.csv")

        fmt="{} uL of {} at valve {}"
        for i,row in self.inventory.iterrows():
            if self.verbose:
                print(fmt.format(row.volume,row.serial,row.valve))

        #Run initialization message
        if self.verbose:
            print("Starting on {} at {}".format(date,time))
            print("RunType={}".format(run_type))
            print("RunID={}".format(s))

    def actual_volume(self,volume,priming=True):#calculator for actual fluid cost from requested volume
        if volume != 0:
            return volume+priming*self.priming_volume
        else:
            return 0

    def determine_row(self,name,volume):
        df=self.inventory
        possibles=df.loc[(df['chemical']==name) & (df['volume']>volume) & (df['volume']>self.floor_volume)]
        if len(possibles)==0:
            raise InventoryVolumeError("Cannot find requisite volume {} ml of chemical {} in inventory".format(volume,name))
        else:
            return possibles.iloc[0]

    def test_inventory(self,name,volume):
        df=self.inventory
        possibles=df.loc[(df['chemical']==name) & (df['volume']>(volume)) & (df['volume']>self.floor_volume)]
        return possibles

    def determine_preprime_row(self,name):
        df=self.inventory
        possibles=df.loc[(df['chemical']==name) & (df['volume']>self.preprime_floor)]
        if len(possibles)==0:
            raise InventoryVolumeError("Cannot find requisite volume {} ml of chemical {} in inventory".format(volume,name))
        else:
            return possibles.iloc[0]

    def check_payload(self,payload,correct_total_volume=1000):
        # round, abs to deal with any negative zeros.
        new={}
        for key in payload:
            new[key]=abs(round(payload[key],self.precision))
        if len(new)==0:
            raise ValueError("Please supply dictionary specifying mixture to measure")
        elif False in [key in list(self.inventory.chemical) for key in new]:
            #print([key in list(self.inventory.chemical) for key in new])
            raise InventorySpeciesError("Chemical supplied in payload that is NOT in Clio's inventory")
        # elif sum(payload.values()) < 5:
        #   return {self.conductivity_string:0}
        total=round(sum(new.values()),self.precision)

        if total != correct_total_volume and total != 0:
            print("Requested volumes do not sum to desired 'correct_total_volume={}mL'; using as proportions instead".format(correct_total_volume))
            for key in new:
                new[key]=new[key]*correct_total_volume/total
        return new

    def run_experiment_fast(self,payload={},priming={},run_order=1,RPM_derate=True,reverse=False,preprime="",mix_time=5000,retain_sample=False,transfer_viscosity=10): #correct_total_volume now in uL
        #HANDLE MALFORMED PAYLOADS
        payload=self.check_payload(payload,correct_total_volume=self.mix_volume)
        volumes_only=payload.copy()

        if len(priming)==0:
            priming = {k:True for k in payload.keys()}
        #ASSUMES NOW-WELL-FORMED PAYLOAD (DECIMAL numbers)
        deductions={}
        concentrations={}
        densities={}
        serials={}
        #payload={'water':3.5,'NaNO3':3.5}
        #water=3.5&NaNO3=3.5
        #REFORM PAYLOAD TO HAVE VALVE
        for key in payload:
            full_volume=self.actual_volume(payload[key],priming=priming[key])
            bottle=self.determine_row(key,full_volume)
            #"{:.2f}".format(n)
            #payload[key]='{}-{}'.format(payload[key],bottle.valve) #FORMAT THE POST DATA
            format_string="{:."+str(self.precision)+"f}"
            if RPM_derate==True:
                rpm=bottle.rpm_derate
            else:
                rpm=1

            payload[key]='{}-{}-{}'.format(format_string.format(payload[key]),bottle.valve,rpm)
            deductions[bottle.name]=full_volume #save index of bottle for deduction from main inventory
            concentrations[bottle.chemical]=bottle.concentration
            densities[bottle.chemical]=bottle["density (g/mL)"]
            serials[bottle.chemical]=bottle.serial
        #Convert dictionary to list of tuples sorted by solution density
        density= lambda k:float(densities[k])
        payload=[(key,payload[key]) for key in sorted(payload.keys(),key=density,reverse=reverse)] #reverse here

        #Get open valve
        df=self.inventory
        open_valve=df.loc[(df['chemical']==self.open_valve_name)].iloc[0].valve

        #Take deductions before posting.
        for i in deductions:
            self.inventory.loc[i,'volume'] = self.inventory.iloc[i].volume - float(deductions[i])

        #Get preprime bottle
        if len(preprime)!=0:
            #Find ingredient valve
            preprime_bottle=self.determine_preprime_row(preprime)
            preprime_valve=preprime_bottle.valve
        else:
            preprime_bottle={}
            preprime_valve=0

        #Determine sample retain
        if retain_sample==False:
            retain_row=None
            retain_valve=0
        else:
            df=self.retain_inventory
            possibles=df.loc[df['volume']==0]
            if len(possibles)==0:
                raise InventoryVolumeError("Cannot find empty retain vial in inventory")
            else:
                retain_row=possibles.iloc[0]
                retain_valve=retain_row.valve

        #ADD INFO KEY TO PAYLOAD
        payload.append(("info","{}-{}-{}-{}-{}-{}-{}-{}".format(self.run_type,self.run_id,run_order,open_valve,preprime_valve,mix_time,max(1,int(transfer_viscosity)),retain_valve)))

        #return payload
        #ASSUME FULLY-FORMED PAYlOAD
        r=requests.post(self.base_url+self.experiment_url_fast,data=payload,timeout=self.timeout)
        r.raise_for_status()

        result=r.json(strict=False) #Returns payload with all results accessible via strings stored in session attributes
        #reverse=False,preprime="",mix_time=5000
        result["reverse"]=reverse
        result["preprime"]=preprime
        result["mix_time"]=mix_time
        result['input_concentrations']=concentrations
        result['input_densities']=densities
        result['input_serials']=serials
        result['postdata']=payload
        result['endpoint']=self.experiment_url_fast
        #with retain_key, decorate the row in retain_inventory with new information from result
        # valve  volume contents  density (g/mL)  conductivity  viscosity bottle_concentrations bottle_densities  qr
        if retain_sample==True:
            self.retain_inventory.loc[retain_row.name,"volume"] += self.mix_volume
            self.retain_inventory.loc[retain_row.name,"contents"] = volumes_only
            self.retain_inventory.loc[retain_row.name,"density"] = result['density (g/mL)']
            self.retain_inventory.loc[retain_row.name,"conductivity"] = result['Cond (mS)']
            self.retain_inventory.loc[retain_row.name,"viscosity"] = result['cP_mean']
            self.retain_inventory.loc[retain_row.name,"bottle_concentrations"] = result['Bottle concentrations']
            self.retain_inventory.loc[retain_row.name,"bottle_densities"] = result['Bottle densities']
        return result

    def sync(self):
        #Check if inventory has changed
        server_inventory = pd.read_csv(self.base_url+"inventory.csv")
        if not server_inventory.equals(self.inventory):
            csv_name='session_inventory.csv'
            #Write current inventory to csv
            self.inventory.set_index('chemical').to_csv(csv_name)
            files = {'file': open(csv_name, 'rb')}
            r=requests.post(self.base_url+'sync',files=files,timeout=self.timeout)
            r.raise_for_status()
            if self.verbose:
                print("Server inventory updated!")
            return
        else:
            print("Server inventory appears to be in sync, no update performed!")
            return

    def sync_retain(self):
        #Check if inventory has changed
        server_inventory = pd.read_csv(self.base_url+"retain_inventory.csv")
        if not server_inventory.equals(self.retain_inventory):
            csv_name='retain_inventory.csv'
            #Write current inventory to csv
            self.retain_inventory.set_index('tray_position').to_csv(csv_name)
            files = {'file': open(csv_name, 'rb')}
            r=requests.post(self.base_url+'sync_retain',files=files,timeout=self.timeout)
            r.raise_for_status()
            if self.verbose:
                print("Server retain_inventory updated!")
            return
        else:
            print("Server retain_inventory appears to be in sync, no update performed!")
            return

    def run_wash(self,run_order=1,retain=False):
        payload={self.wash_solvent:self.wash_volume}
        #HANDLE MALFORMED PAYLOADS
        #ASSUMES NOW-WELL-FORMED PAYLOAD (DECIMAL numbers)
        deductions={}
        concentrations={}
        densities={}
        #REFORM PAYLOAD TO HAVE VALVE
        for key in payload:
            volume=payload[key]
            bottle=self.determine_row(key,payload[key])
            #"{:.2f}".format(n)
            #payload[key]='{}-{}'.format(payload[key],bottle.valve) #FORMAT THE POST DATA
            format_string="{:."+str(self.precision)+"f}"
            payload[key]='{}-{}'.format(format_string.format(payload[key]),bottle.valve)
            deductions[bottle.name]=volume #save index of bottle for deduction from main inventory
            concentrations[bottle.chemical]=bottle.concentration
            densities[bottle.chemical]=bottle["density (g/mL)"]

        #Convert dictionary to list of tuples sorted by solution density
        density= lambda k:float(densities[k])
        payload=[(key,payload[key]) for key in sorted(payload.keys(),key=density)]

        #Get open and waste valve 
        df=self.inventory
        open_valve=df.loc[(df['chemical']==self.open_valve_name)].iloc[0].valve
        waste_valve=df.loc[(df['chemical']==self.waste_valve_name)].iloc[0].valve

        #ADD INFO KEY TO PAYLOAD
        payload.append(("info","{}-{}-{}-{}-{}".format(self.run_type,self.run_id,run_order,open_valve,waste_valve)))
        #ASSUME FULLY-FORMED PAYlOAD
        if retain:
            endpoint_url="RINSE_retain"
        else:
            endpoint_url=self.wash_url
        r=requests.post(self.base_url+endpoint_url,data=payload,timeout=self.timeout)
        r.raise_for_status()

        for i in deductions:
            #self.inventory.loc[i,'volume'] = self.inventory.iloc[i].volume - float(self.actual_volume(deductions[i]))
            self.inventory.loc[i,'volume'] = self.inventory.iloc[i].volume - float(deductions[i]) #don't need to correct for prime volume for washs
        result=r.json(strict=False) #Returns payload with all results accessible via strings stored in session attributes
        result['Bottle concentrations']=concentrations
        result['Bottle densities']=densities
        result['postdata']=payload
        result['endpoint']=self.wash_url
        #result['method']=method #CURRENTLY UNUSED IN CLIO WITHOUT VOLTAGE
        return result


    def branin(self,x1=0,x2=8):
        if (x1<-5) or (x1>10):
            raise ValueError("check x1 domain")
        if (x2<0) or (x2>15):
            raise ValueError("check x2 domain")
        r=requests.post(self.base_url+"branin",data={"x1":x1,"x2":x2})
        r.raise_for_status()
        return r.json(strict=False)

    def run_retain(self,run_volume=1800,payload={},run_order=1,RPM_derate=True,reverse=False,transfer_viscosity=10,output_cid="",prime_line=False): #correct_total_volume now in uL
        #HANDLE MALFORMED PAYLOADS
        payload=self.check_payload(payload,correct_total_volume=run_volume)
        volumes_only=payload.copy()

        #ASSUMES NOW-WELL-FORMED PAYLOAD (DECIMAL numbers)
        deductions={}
        concentrations={}
        densities={}
        serials={}
        #payload={'water':3.5,'NaNO3':3.5}
        #water=3.5&NaNO3=3.5
        #REFORM PAYLOAD TO HAVE VALVE
        for key in payload:
            full_volume=self.actual_volume(payload[key],priming=True)
            bottle=self.determine_row(key,full_volume)
            #"{:.2f}".format(n)
            #payload[key]='{}-{}'.format(payload[key],bottle.valve) #FORMAT THE POST DATA
            format_string="{:."+str(self.precision)+"f}"
            if RPM_derate==True:
                rpm=bottle.rpm_derate
            else:
                rpm=1

            payload[key]='{}-{}-{}'.format(format_string.format(payload[key]),bottle.valve,rpm)
            deductions[bottle.name]=full_volume #save index of bottle for deduction from main inventory
            concentrations[bottle.chemical]=bottle.concentration
            densities[bottle.chemical]=bottle["density (g/mL)"]
            serials[bottle.chemical]=bottle.serial
        #Convert dictionary to list of tuples sorted by solution density
        density= lambda k:float(densities[k])
        payload=[(key,payload[key]) for key in sorted(payload.keys(),key=density,reverse=reverse)] #reverse here

        #Get open valve
        df=self.inventory
        open_valve=df.loc[(df['chemical']==self.open_valve_name)].iloc[0].valve

        #Take deductions before posting.
        for i in deductions:
            self.inventory.loc[i,'volume'] = self.inventory.iloc[i].volume - float(deductions[i])

        df=self.retain_inventory
        if prime_line:
            retain_valve=13
        else:
            possibles=df.loc[df['volume']==0]
            if len(possibles)==0:
                raise InventoryVolumeError("Cannot find empty retain vial in inventory")
            else:
                retain_row=possibles.iloc[0]
                retain_valve=retain_row.valve

        #ADD INFO KEY TO PAYLOAD
        payload.append(("info","{}-{}-{}-{}-{}-{}".format(self.run_type,self.run_id,run_order,open_valve,retain_valve,max(1,int(transfer_viscosity)))))

        #return payload
        #ASSUME FULLY-FORMED PAYlOAD
        #pdb.set_trace()

        r=requests.post(self.base_url+self.retain_url,data=payload,timeout=self.timeout)
        r.raise_for_status()

        result=r.json(strict=False) #Returns payload with all results accessible via strings stored in session attributes
        #reverse=False,preprime="",mix_time=5000
        result["reverse"]=reverse
        result['input_concentrations']=concentrations
        result['input_densities']=densities
        result['input_serials']=serials
        result['postdata']=payload
        result['endpoint']=self.retain_url
        result['retain_valve']=retain_valve
        result["output_cid"]=output_cid
        #with retain_key, decorate the row in retain_inventory with new information from result
        # valve  volume contents  density (g/mL)  conductivity  viscosity bottle_concentrations bottle_densities  qr
        if not prime_line:
            new_row={}
            new_row["volume"]=run_volume
            new_row["contents"]=str(volumes_only)
            new_row["conductivity"]=result["Cond (mS) 2"]
            new_row["input_concentrations"]=result["input_concentrations"]
            new_row["input_densities"]=result["input_densities"]
            new_row["output_cid"]=output_cid
            new_row["input_serials"]=result["input_serials"]
            self.retain_inventory.loc[retain_row.name,new_row.keys()] = new_row.values()
        return result

#!/usr/bin/env python
import numpy as np
import gsd
import gsd.hoomd
import hoomd, hoomd.md as md
from hoomd import azplugins
import argparse


parser = argparse.ArgumentParser(description='Performs a test NVT run at correct temerature')
parser.add_argument('-P',dest='P',action='store',required=True,help='Pressure')
args = parser.parse_args()
Temp =250
Press = float(args.P)

def get_param_dict():
    aaparams = {}
    with open('stats_module.dat','r') as f:
        for line in f:
            if line[0]!='#':
                tmp = line.split()
                aaparams[tmp[0]]=np.array(tmp[1:],dtype='float')
    return aaparams

def paramparse():
    aamass=[]
    aacharge=[]
    aaradius=[]
    aahps=[]
    aaparams = get_param_dict();
    aakeys = list(aaparams.keys())
    for i in aakeys:
        aamass.append(aaparams[str(i)][0])
        aacharge.append(aaparams[str(i)][1])
        aaradius.append(aaparams[str(i)][2])
        aahps.append(aaparams[str(i)][3])
    return aamass,aacharge,aaradius,aahps,aakeys,aaparams


def chain_parse(seqfile):
    dictseq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP','Z':'STR'}

    with open(seqfile,'r') as f:
        data = f.readlines()
    init_seq1 = data[0].strip()
    seq1 = []
    for i in init_seq1:
        seq1.append(dictseq[i])
    aamass,aacharge,aaradius,aahps,aakeys,aaparams = paramparse()
    ##Translate sequeunce ##
    chain_id=[]
    chain_mass=[]
    chain_charge=[]
    for i in seq1:
        index = aakeys.index(i)
        chain_id.append(index)
        chain_mass.append(aamass[index])
        chain_charge.append(aacharge[index])
    return chain_id,chain_mass,chain_charge,aakeys,aaparams


#Production run parameters
production_dt=0.01 # Time step for production run in picoseconds
production_steps=100000000 # Total number of steps
production_T=Temp # Temperature for production run in Kelvin

seqfile = 'sv28.dat'
chain_id, chain_mass, chain_charge,aakeys,aaparams = chain_parse(seqfile)
##Create the snapshot ##
bond_length=0.38
chain_length=len(chain_id)
box_length=bond_length*chain_length+10

#################################################################################################
# ### Minimized slab formed and saved in minimize.gsd
#---------------------------------------------------------------------------------------------
# ## 4.0. Run a production slab simulation using minimize.gsd from previous step
################################################################################################
hoomd.context.initialize("--notice-level=2")
sim = hoomd.context.SimulationContext()
system = hoomd.init.read_gsd('restart_tmp1_%3i.gsd'%(Temp),frame=-1)

n_steps = production_steps # 1 microseconds

fileroot = 'Production'
nl = hoomd.md.nlist.cell()

## Bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond', k=8360, r0=0.381)
harmonic.bond_coeff.set('STR_bond',k=8368,r0=1.0)
## Nonbonded
nl.reset_exclusions(exclusions=['1-2', 'body'])
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aaparams[i][3]+aaparams[j][3])/2.,epsilon=0.8368, sigma=(aaparams[i][2]+aaparams[j][2])/10./2.,r_cut=2.0)

## Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
for i,atom1 in enumerate(aakeys):
    for j,atom2 in enumerate(aakeys):
        yukawa.pair_coeff.set(atom1,atom2,epsilon=aaparams[atom1][1]*aaparams[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)

## Group Particles
all = hoomd.group.all()

## Set up integrator
hoomd.md.integrate.mode_standard(dt=production_dt) # Time units in ps
temp = production_T*0.00831446
integrator = hoomd.md.integrate.npt(group=all, kT=temp, tau = 0.5, tauP = 0.5, P = Press, couple='none', x=False, y=False, z=True) # Temp is kT/0.00831446
## Outputs
hoomd.analyze.log(filename=fileroot+'P%.2f_R%.3f.log'%(Press,Temp), quantities=['potential_energy', 'temperature','lx','ly','lz'], period=10000, overwrite=False, header_prefix='#')
hoomd.dump.gsd('restartP_tmp2_P%.2f_R%.3f.gsd'%(Press,Temp), period=100000, group=all)
hoomd.analyze.log(filename='stress_P%.2f_%.3f.log'%(Press,Temp), quantities=['pressure_xx', 'pressure_yy', 'pressure_zz','lx','ly','lz'], period=50, overwrite=False, header_prefix='#') # Output stress tensor
## Run simulation
hoomd.run(production_steps)
########################################################################################################

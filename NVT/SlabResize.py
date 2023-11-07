#!/usr/bin/env python
import numpy as np
import gsd
import gsd.pygsd
import gsd.hoomd
import argparse
import itertools as it
import argparse
import hoomd
from hoomd import azplugins
parser = argparse.ArgumentParser(description='Performs a test NVT run at correct density')
parser.add_argument('-T',dest='T',action='store',required=True,help='Temp')
args = parser.parse_args()

## Parameters ##
Temp = float(args.T)
nstar = 147
Narm = 3
Nchain2 = 576

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


#Resizing run parameters
resize_dt=0.01 # Time step for production run in picoseconds
resize_steps=500000 # Total number of steps
resize_T=Temp # Temperature for production run in Kelvin

seqfile1 = 'sv1.dat'
chain_id1, chain_mass1, chain_charge1,aakeys,aaparams = chain_parse(seqfile1)

seqfile2 = 'sv28.dat'
chain_id2, chain_mass2, chain_charge2,aakeys,aaparams = chain_parse(seqfile2)

bond_length=0.38
chain_length1=len(chain_id1)
chain_length2=len(chain_id2)

# #################################################################################################
# # start_sv1.gsd contains both sv1x3 and sv28 in a cubic box. This simulation compresses the box
# # and elongates it to perform a direct coexistence simulation for equilibration of the system.
# ################################################################################################

hoomd.context.initialize("--notice-level=2")
sim = hoomd.context.SimulationContext()
system = hoomd.init.read_gsd('start_sv1.gsd')

boxsize = 20 # Size of the compressed simulation box

harmonic=hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('AA_bond',k=8368,r0=bond_length)
harmonic.bond_coeff.set('STR_bond',k=8368,r0=1.0)
# #### Neighborlist and exclusions
nl = hoomd.md.nlist.cell()
nl.reset_exclusions(exclusions=['1-2', 'body'])

#### Pairwise interactions
nb = azplugins.pair.ashbaugh(r_cut=0, nlist=nl)
for i in aakeys:
    for j in aakeys:
        nb.pair_coeff.set(i,j,lam=(aaparams[i][3]+aaparams[j][3])/2.,epsilon=0.8368, sigma=(aaparams[i][2]+aaparams[j][2])/10./2.,r_cut=2.0)

# #### Electrostatics
yukawa = hoomd.md.pair.yukawa(r_cut=0.0, nlist=nl)
for i,atom1 in enumerate(aakeys):
    for j,atom2 in enumerate(aakeys):
        yukawa.pair_coeff.set(atom1,atom2,epsilon=aaparams[atom1][1]*aaparams[atom2][1]*1.73136, kappa=1.0, r_cut=3.5)
# #### Group Particles
all = hoomd.group.all()
# #### Set up integrator
hoomd.md.integrate.mode_standard(dt=resize_dt) # Time units in ps
kTinput=resize_T * 8.3144598/1000.
integrator = hoomd.md.integrate.langevin(group=all, kT=kTinput, seed=63535)
# #### Resize the box after replication to 15x15x15nm
hoomd.update.box_resize(L=hoomd.variant.linear_interp([(0,system.box.Lx),(resize_steps-500,boxsize)]),scale_particles=True)
for cnt,i in enumerate(aakeys):
    integrator.set_gamma(i,gamma=aaparams[aakeys[cnt]][0]/1000.0)
# #### Output log file with box dimensions and restart file after box resizing
hoomd.analyze.log(filename='resize.log', quantities=['potential_energy','kinetic_energy','temperature','pressure_xx','pressure_yy','pressure_zz','lx','ly','lz'], period=10000, overwrite=True, header_prefix='#')
hoomd.dump.gsd('resize_%3i.gsd'%(Temp), period=10000, group=all)
# #### Run resizing simulation
hoomd.run(tsteps=resize_steps)


slab_z_length = 120
###############################################################################################
# # resize.gsd contains the resized cubic box, which is then extended in z to slab_z_length and
# # periodic boundary conditions are corrected by unwrapping system.
###############################################################################################

def extend(s):
    nres = chain_length1
    boxdim = s.configuration.box[:3]
    zmin,zmax,dz = -boxdim[2]/2., boxdim[2]/2., boxdim[2]
    print(dz)
    pos1 =  s.particles.position
    pos = pos1.copy()
    pos2 = pos1.copy()[:,2]
    skip=0
    for i in range(nstar):
        centerid = (nres*Narm+1)*i #Central particle id
        for j in range(Narm):
            mol_id = np.arange(centerid+nres*j+1,centerid+nres*j+1+nres)
            mol_id = np.append(centerid,mol_id)
            mol_pos = pos2[mol_id]
            for k in range(1,nres+1):
                dist2 = (mol_pos[k] - mol_pos[k-1])**2
                if dist2 > 8:
                    excess = np.sign(mol_pos[k] - mol_pos[k-1])*dz
                    mol_pos[k] = mol_pos[k] - excess
                pos[centerid+nres*j+1:centerid+nres*j+1+nres,2] = mol_pos[1:]
        skip += Narm*nres+1
    nchain = Nchain2
    nres = chain_length2
    for i in range(nchain):
        mol_coord = pos[i*nres+skip:(i+1)*nres+skip,2]
        for j in range(1,nres):
            dist2 = (mol_coord[j] - mol_coord[j-1])**2
            if dist2 > 8:
                excess = np.sign(mol_coord[j] - mol_coord[j-1])*dz
                mol_coord[j] = mol_coord[j] - excess
            com = np.mean(mol_coord)
            if com < zmin:
                mol_coord += dz
            elif com > zmax:
                mol_coord -= dz
        pos[i*nres+skip:(i+1)*nres+skip,2] = mol_coord
    return pos

f = gsd.pygsd.GSDFile(open('resize_%3i.gsd'%(Temp),'rb'))
t = gsd.hoomd.HOOMDTrajectory(f)
s1 = t[-1]
s = gsd.hoomd.Snapshot()
s.particles.N = s1.particles.N
s.particles.types = s1.particles.types
s.particles.typeid = s1.particles.typeid
s.particles.mass = s1.particles.mass
s.particles.charge = s1.particles.charge
s.particles.position = extend(s1)
s.bonds.N = s1.bonds.N
s.bonds.types = s1.bonds.types
s.bonds.typeid = s1.bonds.typeid
s.bonds.group = s1.bonds.group
s.configuration.box = s1.configuration.box
s.configuration.dimensions=3
s.configuration.box = [s1.configuration.box[0],s1.configuration.box[1],slab_z_length,0,0,0]
s.configuration.step = 0
outfile = gsd.hoomd.open('box2slab_extend_%3i.gsd'%(Temp),'wb')
outfile.append(s)
outfile.close()

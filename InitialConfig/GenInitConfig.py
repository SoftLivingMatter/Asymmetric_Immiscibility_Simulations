#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
import numpy as np
import gsd
import gsd.pygsd
import gsd.hoomd
import argparse
import itertools as it


# In[2]:


seqname = 'ke1'
Narm = 3


# In[3]:


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

def chain_parse(seqfile,valence):
    dictseq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP','Z':'STR'}

    with open(seqfile,'r') as f:
        data = f.readlines()
    init_seq1 = data[0].strip()*valence
    if valence != 1:
        seq1 = ['STR']
    else:
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


# In[4]:

seqfile1 = seqname + '.dat'
chain_id1, chain_mass1, chain_charge1,aakeys,aaparams = chain_parse(seqfile1,Narm)

seqfile2 = 'ke7.dat'
chain_id2, chain_mass2, chain_charge2,aakeys,aaparams = chain_parse(seqfile2,1)


# In[5]:


## Reading initial coordinates of folded ke1 star polymer generated from preliminary simulations##
f = gsd.hoomd.open(name='start_ke1_3arm.gsd', mode='rb')
s1 = f[0]
pos = s1.particles.position
Nparticles_old = s1.particles.N
config = s1.configuration.box[0]
Nstar = 200
Nstarres = len(pos)
print(len(pos))


# In[6]:


##Positioning of Stars in upper half of cube##
m = int(np.floor((2*Nstar)**(1/3)))
NewCubeL = m*config
diff = (NewCubeL - config)/2
print(m,NewCubeL)
newpos = []
ctr = 0
for i in range(m):
    for j in range(m):
        for k in range(int(m/2)):
            ctr += 1
            newpos.append(pos - [-diff+i*config,-diff + j*config,-diff+k*config])
Nstar = ctr
newpos = np.array(newpos).reshape((Nstar*Nstarres,3))
print(Nstar)

## Getting type, mass and charge data for star polymers ##
arrtypeids = np.tile(chain_id1,Nstar)
arrmass = np.tile(chain_mass1,Nstar)
arrcharge = np.tile(chain_charge1,Nstar)
N1 = len(arrtypeids)


# In[7]:


## No of sv28 particles ##
chain_length2 = len(chain_id2)
Nchain2 = 460
Nparticles2 = chain_length2*Nchain2
bond_length = 0.38


# In[8]:


##Positioning of Linear Poly in lower half of cube##
Nparticles2 = Nchain2*chain_length2
layerspacing2 = bond_length*chain_length2
Nlayers2 = np.ceil((NewCubeL/2-bond_length*chain_length2)/layerspacing2)
Npartlayer2 = Nchain2/Nlayers2
planespacing2 = (NewCubeL)/np.sqrt(Npartlayer2)

xposarr2 = np.arange(-NewCubeL/2+1,NewCubeL/2-1,planespacing2)
yposarr2 = np.arange(-NewCubeL/2+1,NewCubeL/2-1,planespacing2)
zposarr2 = np.arange(-NewCubeL/2+1,-bond_length*chain_length2,layerspacing2)

j = 0
for x in xposarr2:
    for y in yposarr2:
        for z in zposarr2:
            for k in range(chain_length2):
                if j >= Nparticles2:
                    break
                j = j+1
                newpos = np.vstack((newpos,[x,y,z + k*bond_length]))

print(Nchain2,j/50,len(newpos))
Nchain2 = int(j/50)
print(Nchain2)
## Getting type, mass and charge data for linear polymers ##
for Nctr in range(Nchain2):
        arrtypeids = np.hstack((arrtypeids,chain_id2))
        arrmass = np.hstack((arrmass,chain_mass2))
        arrcharge = np.hstack((arrcharge,chain_charge2))
N2 = len(arrtypeids)
print(N1,N2,N1/N2) #Atoms of comp 1, atoms of comp2 and mass fraction


# In[9]:


##Bonding of Stars in upper half of cube##
Nbonds = s1.bonds.N
Nbonds_type = s1.bonds.types
Nbonds_typeid = s1.bonds.typeid
Nbonds_group_old = s1.bonds.group
Nbonds_group = np.copy(Nbonds_group_old)
for i in range(1,Nstar):
    grp_append = Nbonds_group_old+i*np.array([151,151])
    Nbonds_group = np.vstack((Nbonds_group,grp_append))
Nbonds_typeid1 = np.tile(Nbonds_typeid,Nstar)


# In[11]:


##Bonding of linear polymers in lower half of cube##
Nparticles1 = Nstar*Nparticles_old
Nparticles = Nparticles1 + Nchain2*chain_length2
Nids2 = np.arange(Nparticles1,Nparticles)
init_blist2 = zip(Nids2,Nids2[1:])
bond_pairs2 = list(it.filterfalse(lambda i: (np.sum(i)+1-2*Nparticles1)%(2*chain_length2)==0, init_blist2))
Nbonds_typeid2 = len(bond_pairs2)*[0]

bond_pairs = np.vstack((Nbonds_group,bond_pairs2))
Nbonds_new = len(bond_pairs)
Nbonds_typeid = np.hstack((Nbonds_typeid1,Nbonds_typeid2))
print(len(Nbonds_typeid),Nbonds_new)


# In[12]:


s = gsd.hoomd.Snapshot()
s.particles.N = Nstar*Nparticles_old + Nchain2*chain_length2
s.particles.types = s1.particles.types
s.particles.typeid = arrtypeids
s.particles.mass = arrmass
s.particles.charge = arrcharge
s.particles.position = newpos
s.bonds.N = Nbonds_new
s.bonds.types = s1.bonds.types
s.bonds.typeid = Nbonds_typeid
s.bonds.group = bond_pairs
s.configuration.box = s1.configuration.box
s.configuration.dimensions=3
s.configuration.box = [NewCubeL,NewCubeL,NewCubeL,0,0,0]
s.configuration.step = 0
outfile = gsd.hoomd.open('start_'+seqname+'.gsd','wb')
outfile.append(s)
outfile.close()

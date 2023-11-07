import numpy as np
import gsd
import gsd.hoomd
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
from matplotlib.image import imread
from tempfile import NamedTemporaryFile
import glob

seqname = 'sv1'

def getmean(z_pos, L, masswt):
    Box_min = -L/2
    Box_max = L/2
    Box_width = Box_max - Box_min
    Points_map = (z_pos - Box_min)/Box_width*2*np.pi-np.pi
    A = np.dot(np.sin(Points_map),masswt)/np.sum(masswt)
    B = np.dot(np.cos(Points_map),masswt)/np.sum(masswt)
    Pmean_map = np.arctan2(A,B)
    Pmean = (Pmean_map+np.pi)/(2*np.pi) * Box_width + Box_min
    return Pmean

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

Ratio = 0.5
Temp = 250
traj = gsd.hoomd.open('restart_tmp1_%3i.gsd'%(Temp))
s=traj[0]
Nparticles = s.particles.N
Narm =3

seqfile1 = seqname +'.dat'
chain_id1, chain_mass1, chain_charge1,aakeys,aaparams = chain_parse(seqfile1,Narm)

seqfile2 = 'sv28.dat'
chain_id2, chain_mass2, chain_charge2,aakeys,aaparams = chain_parse(seqfile2,1)

Nchains1 = 147
Nchains2 = 576
chain_length1=len(chain_id1)
chain_length2=len(chain_id2)

part_types = s.particles.typeid
body_arr = np.ones(len(part_types))
body_arr[Nchains1*chain_length1:] = 2

ff_para = 'stats_module.dat'
aalist={}
with open(ff_para,'r') as fid:
    for i in fid:
        if i[0]!='#':
            tmp=i.rsplit()
            aalist[tmp[0]]=np.loadtxt(tmp[1:],dtype=float)
aakeys=list(aalist.keys())
# This translates each amino acid type into a number, which will be used in HOOMD
# For example, GLY is with an ID of 10
aamass=[]
for i in aakeys:
    aamass.append(aalist[i][0])

skip = 500
slices = 5
traj = gsd.hoomd.open('restart_tmp1_%3i.gsd'%(Temp))[skip:]
length = len(traj)
s1=traj[0]
L = s1.configuration.box[2]
Lx = s1.configuration.box[0]
Box_min = -L/2
Box_max = L/2
Nparticles = s1.particles.N
zbins = np.arange(-L/2, L/2, 1)
cumulative_hist_arr = np.zeros((slices,len(zbins)-1))
cumulative_hist_arr1 = np.zeros((slices,len(zbins)-1))

for sl in range(slices):
    print(sl*length/slices,(sl+1)*length/slices)
    traj = gsd.hoomd.open('restart_tmp1_%3i.gsd'%(Temp))[skip:][int(sl*length/slices):int((sl+1)*length/slices)]
    cumulative_hist = np.zeros(len(zbins)-1)
    for frame in traj:
        z_pos = frame.particles.position[:, 2]
        z_pos_ids = frame.particles.typeid
        masswt = [aamass[x] for x in z_pos_ids]
        z_recenter = z_pos - getmean(z_pos,L,masswt)
        z_recenter = np.array([x - L if x>Box_max else x+L if x<Box_min else x for x in z_recenter])
        for i in range(21):
            mask = z_pos_ids == i
            if not any(mask):
                continue
            mask2 = body_arr[mask] == 1
            z_recenter_mask = z_recenter[mask][mask2]
            y, binedges = np.histogram(z_recenter_mask, bins=zbins)
            cumulative_hist += 1.66054*aamass[i]*(y/(Lx**2))
    bincenters = 0.5*(zbins[1:]+zbins[:-1])
    cumulative_hist = cumulative_hist/len(traj)
    cumulative_hist_arr[sl] = cumulative_hist

for sl in range(slices):
    print(sl*length/slices,(sl+1)*length/slices)
    traj = gsd.hoomd.open('restart_tmp1_%3i.gsd'%(Temp))[skip:][int(sl*length/slices):int((sl+1)*length/slices)]
    cumulative_hist1 = np.zeros(len(zbins)-1)
    for frame in traj:
        z_pos = frame.particles.position[:, 2]
        z_pos_ids = frame.particles.typeid
        masswt = [aamass[x] for x in z_pos_ids]
        z_recenter = z_pos - getmean(z_pos,L,masswt)
        z_recenter = np.array([x - L if x>Box_max else x+L if x<Box_min else x for x in z_recenter])
        for i in range(21):
            mask = z_pos_ids == i
            if not any(mask):
                continue
            mask2 = body_arr[mask] == 2
            z_recenter_mask = z_recenter[mask][mask2]
            y, binedges = np.histogram(z_recenter_mask, bins=zbins)
            cumulative_hist1 += 1.66054*aamass[i]*(y/(Lx**2))
    bincenters = 0.5*(zbins[1:]+zbins[:-1])
    cumulative_hist1 = cumulative_hist1/len(traj)
    cumulative_hist_arr1[sl] = cumulative_hist1

fname = 'z_dens_'+seqname
np.save(fname,[cumulative_hist_arr,cumulative_hist_arr1])
np.save('bincenters',bincenters)

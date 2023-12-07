import numpy as np


def get_param_dict(paramfile='stats_module.dat'):
    aaparams = {}
    for line in paramfile:
        if line.startswith('#'):
            continue
        aa, *params = line.split()
        aaparams[aa]=np.array(params, dtype='float')
    return aaparams


def paramparse(paramfile):
    aas = []
    masses = []
    charges = []
    radii = []
    hpss = []
    aaparams = get_param_dict(paramfile);
    for aa, params in aaparams.items():
        mass, charge, radius, hps = params
        aas.append(aa)
        masses.append(mass)
        charges.append(charge)
        radii.append(radius)
        hpss.append(hps)
    return masses, charges, radii, hpss, aas, aaparams


def chain_parse(seqfile, valence, paramfile):
    dictseq={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
     'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
     'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
     'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
     'W':'TRP','Z':'STR'}

    data = seqfile.readlines()
    print(data)
    init_seq1 = data[0].strip()*valence
    if valence != 1:
        seq1 = ['STR']
    else:
        seq1 = []
    for i in init_seq1:
        seq1.append(dictseq[i])
    aamass,aacharge,aaradius,aahps,aakeys,aaparams = paramparse(paramfile)
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


def chain_parse(seqfile, valence):
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

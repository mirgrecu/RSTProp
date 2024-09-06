from netCDF4 import Dataset
class scattTables:
    fh=Dataset("lookupTables/scatteringTablesGPM.nc")
    fhBB=Dataset("lookupTables/bbScatteringProp.nc")
    fhGMI=Dataset("lookupTables/scatteringTablesGMI.nc")
    #print(fhGMI)
    #freqs={10.97,19,22,37,85,166,187} [3] is Ka
    zKuR=fh["zKuR"][:289]
    zKaR=fh["zKaR"][:289]
    dmr=fh["dmr"][:289]
    rainRate=fh["rainRate"][:289]
    attKuR=fh["attKuR"][:289]
    attKaR=fh["attKaR"][:289]
    dmr=fh["dmr"][:289]
    kextR=fhGMI["kextR"][:289,:7]
    salbR=fhGMI["salbR"][:289,:7]
    asymR=fhGMI["asymR"][:289,:7]
    rwc=fh["rwc"][:289]
    #--------------------#
    zKuS=fh["zKuS"][:253]
    zKaS=fh["zKaS"][:253]
    dms=fh["dms"][:253]
    snowRate=fh["snowRate"][:253]
    swc=fh["swc"][:253]
    attKuS=fh["attKuS"][:253]
    attKaS=fh["attKaS"][:253]
    kextS=fhGMI["kextS"][:253,:7]
    salbS=fhGMI["salbS"][:253,:7]
    asymS=fhGMI["asymS"][:253,:7]
    #---------------------#
    zKuBB=fhBB["zKuBB"][:289]
    zKaBB=fhBB["zKaBB"][:289]
    dmBB=fhBB["dmBB"][:289]
    precRateBB=fhBB["pRateBB"][:289]
    attKuBB=fhBB["attKuBB"][:289]
    attKaBB=fhBB["attKaBB"][:289]
    kextBB=fhBB["kextBB"][:289,:7]
    salbBB=fhBB["salbBB"][:289,:7]
    asymBB=fhBB["asymBB"][:289,:7]
    
    #------------------#
    dmg=fh["dmg"][:272]
    zKuG=fh["zKuG"][:272]
    zKaG=fh["zKaG"][:272]
    dmg=fh["dmg"][:272]
    graupRate=fh["graupRate"][:272]
    gwc=fh["gwc"][:272]
    attKuG=fh["attKuG"][:272]
    attKaG=fh["attKaG"][:272]
    dmg=fh["dmg"][:272]
    kextG=fh["kextG"][:272,:7]
    salbG=fh["salbG"][:272,:7]
    asymG=fh["asymG"][:272,:7]
    kextG=fhGMI["kextG"][:272,:7]
    salbG=fhGMI["salbG"][:272,:7]
    asymG=fhGMI["asymG"][:272,:7]
    fhGMI.close()

    
def getGraupProp(zc,dnw,lkT):
    if zc>12:
        ibin=int((zc-10*dnw+12)/0.25)
        if ibin<=0:
            ibin=0
            dnw=(zc+12)/10.
        if ibin>=271:
            ibin=271
            dnw=(zc-lkT.zKaG[271])/10.
        zka=lkT.zKaG[ibin]+10*dnw
        attKa=lkT.attKaG[ibin]*10**dnw
        pRate=lkT.graupRate[ibin]*10**dnw
    else:
        zka=-99
        attKa=0
        pRate=0
    return zka,attKa,pRate

def getSnowProp(zc,dnw,lkT):
    if zc>12:
        ibin=int((zc-10*dnw+12)/0.25)
        if ibin<=0:
            ibin=0
            dnw=(zc+12)/10.
        if ibin>=252:
            ibin=252
            dnw=(zc-lkT.zKaS[252])/10.
        zka=lkT.zKaS[ibin]+10*dnw
        attKa=lkT.attKaS[ibin]*10**dnw
        pRate=lkT.snowRate[ibin]*10**dnw
        kext_Ka=lkT.kextS[ibin,3]*10**dnw
        salb_Ka=lkT.salbS[ibin,3]
        asym_Ka=lkT.asymS[ibin,3]
    else:
        zka=-99
        attKa=0
        pRate=0
        kext_Ka,salb_Ka,asym_Ka=0,0,0
    return zka,attKa,pRate,kext_Ka,salb_Ka,asym_Ka


def getBBProp(zc,dnw,lkT):
    if zc>12:
        ibin=int((zc-10*dnw+12)/0.25)
        if ibin<=0:
            ibin=0
            dnw=(zc+12)/10.
        if ibin>=271:
            ibin=0
            dnw=(zc-lkT.zKuBB[271])/10.
        zka=lkT.zKaBB[ibin]+10*dnw
        attKa=lkT.attKaBB[ibin]*10**dnw
        pRate=lkT.precRateBB[ibin]*10**dnw
        kext_Ka=lkT.kextBB[ibin,3]*10**dnw
        salb_Ka=lkT.salbBB[ibin,3]
        asym_Ka=lkT.asymBB[ibin,3]
    else:
        zka=-99
        attKa=0
        pRate=0
        kext_Ka,salb_Ka,asym_Ka=0,0,0
    return zka,attKa,pRate,kext_Ka,salb_Ka,asym_Ka

def getRainProp(zc,dnw,lkT):
    if zc>12:
        ibin=int((zc-10*dnw+12)/0.25)
        if ibin<=0:
            ibin=0
            dnw=(zc+12)/10.
        if ibin>=288:
            ibin=0
            dnw=(zc-lkT.zKaR[288])/10.
        zka=lkT.zKaR[ibin]+10*dnw
        attKa=lkT.attKaR[ibin]*10**dnw
        pRate=lkT.rainRate[ibin]*10**dnw
        kext_Ka=lkT.kextR[ibin,3]*10**dnw
        salb_Ka=lkT.salbR[ibin,3]
        asym_Ka=lkT.asymR[ibin,3]
    else:
        zka=-99
        attKa=0
        pRate=0
        kext_Ka,salb_Ka,asym_Ka=0,0,0
    return zka,attKa,pRate,kext_Ka,salb_Ka,asym_Ka
        

import numpy as np

def HI_mass(S,dl,Abeam,dv):
    Sv = np.sum(np.sum(S, axis=1))
    MHI = 2.35*10**5*(dl)**2*(Sv)*(1/Abeam)
    logMHI = np.log10(MHI)
    N = 0
    for i in range(0,len(S)):
        for j in range(0,len(S[0])):
                N += 1
    errsv = np.sqrt(N/Abeam)*dv*np.sqrt(np.mean(S**2))
    errMHI = (logMHI/np.log(10))*(errsv/Sv)
    return str(np.round(logMHI, 2)) + ' $\pm$ ' + str(np.round(errMHI, 2)) + '  not log = ' + str(np.round(MHI,2)) + ', ' + str(np.round(errsv,2))

def flux(S,Abeam,dv):
    Sv = np.sum(np.sum(S, axis=1))
    F = Sv*Abeam/(dv)
    errsv = np.sqrt((len(S)*len(S[0]))/Abeam)*dv*np.sqrt(np.mean(S**2))
    logF = np.log10(F)
    errF = np.abs((logF/np.log(10))*(errsv/Sv))
    return  '$ ' + str(np.round(logF, 5)) + ' \pm ' + str(np.round(errF, 2)) + ' $  not log = ' + str(np.round(F,2))

def galNHI(NHI):
    mNHI = NHI[NHI > 0]
    totNHI = np.mean(mNHI)
    return(str(np.round(totNHI, 2)))
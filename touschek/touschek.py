from scipy import interpolate
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from config import get_from_config
import fur.path_assistant as path_assistant
from lattice import lattice
from scipy import integrate
from scipy.special import iv


def ExI0(a, b):
    if b < 40:
        return np.exp(-a)*iv(0, b)
    else:
        return np.exp(-a+b)/np.sqrt(2*np.pi*b)\
            * (1+1/8/b+9/128/b**2+9*25/6/8**3/b**3)


def Ku(u, um, B1, B2):
    return np.sqrt(u/(1+u))*ExI0(B1*u, B2*u)\
        * ((2+1/u)**2*(u/um/(1+u)-1)+1-np.sqrt(um*(1+u)/u)
           - 1/2/u*(4+1/u)*np.log(u/um/(1+u)))


def Itsk(um, B1, B2):
    return integrate.quad(lambda u: Ku(u, um, B1, B2), um, np.inf)


gamma0 = get_from_config("gamma")
e = get_from_config("e")
c = get_from_config("c_m/s")*100
me = get_from_config("me_MeV")*1e6
alpha = get_from_config("ring_alpha")
q = get_from_config("RF_q")
VSR = get_from_config("VSR")
f0 = 1/get_from_config("IOTA_revolution_period")
fRF = q*f0
lamRF = c/fRF


re = 2.817941e-13  # cm


touschek_func_df = pd.read_csv(os.path.join(
    os.path.dirname(__file__), "touschek_func.csv"), header=None)
CI = interpolate.interp1d(touschek_func_df.iloc[:, 0],
                          touschek_func_df.iloc[:, 1])


def get_LamTska(lattice_df, Vrf, sp, ex, sz, Ibeam,
                aperture_factor=1.0,
                use_transverse_acceptance_octupole=False,
                use_constant_transverse_acceptance=False,
                constant_acceptance_cm=5.0,
                use_transverse_acceptance_by_emittance=False,
                emittance_acceptance_um=22.1,
                use_variable_aperture_df=False,
                variable_aperture_df=None,
                gamma=gamma0):
    """The result assumes ey=1. So it has to be devided by sqrt(ey) in units of um"""
    V0 = Vrf
    etas = alpha-1/gamma**2
    Ks = (alpha*gamma**2-1)/(gamma**2-1)
    phiacc = np.arcsin(VSR/V0)
    nus0 = np.sqrt(q*V0*np.abs(Ks)/2/np.pi/me/gamma)
    nus = nus0*np.sqrt(np.cos(phiacc))
    fs = f0*nus
    dP_Psep = aperture_factor*2*nus0/q/np.abs(etas)*np.sqrt(np.cos(phiacc)
                                            - (np.pi/2-phiacc)*np.sin(phiacc))
    ldf = lattice_df
    if use_transverse_acceptance_octupole:
        idx_ac = 793
        a_ac_cm = 0.38
        beta_ac_cm = lattice_df.loc[idx_ac, 'Beta_cm_X']
        dispersion_ac_cm = lattice_df.loc[idx_ac, 'Dispersion_cm_X']
        daccL = (a_ac_cm\
            / (np.sqrt(lattice_df['H']*beta_ac_cm)\
            + np.abs(dispersion_ac_cm))).values
        dP_Psep = np.where(daccL < dP_Psep, daccL, dP_Psep)
    if use_constant_transverse_acceptance:
        daccL = \
        [(constant_acceptance_cm\
            / (np.sqrt(ldf.loc[i, 'H']*ldf['Beta_cm_X'])
               + np.abs(ldf['Dispersion_cm_X']))).min() for i in ldf.index]
        daccL = np.array(daccL)
        dP_Psep = np.where(daccL < dP_Psep, daccL, dP_Psep)
    if use_transverse_acceptance_by_emittance:
        ap_vs_S = (1e-4*lattice.get_sigma_um(ldf['Beta_cm_X'],   
                emittance_acceptance_um, ldf['Dispersion_cm_X'], sp))
        daccL = np.asarray([
            np.abs(ap_vs_S/(np.sqrt(ldf['Beta_cm_X']*ldf.loc[i, 'H'])+ldf['Dispersion_cm_X'].abs())).min()
            for i in ldf.index])
        dP_Psep = np.where(daccL < dP_Psep, daccL, dP_Psep)
    if use_variable_aperture_df:
        ap_vs_S = np.interp(ldf['S_cm'], variable_aperture_df['S_cm'],
                    variable_aperture_df['Aperture_cm_X'],
                    left=variable_aperture_df['Aperture_cm_X'].values[0],
                    right=variable_aperture_df['Aperture_cm_X'].values[-1])
        daccL = np.asarray([
            np.abs(ap_vs_S/(np.sqrt(ldf['Beta_cm_X']*ldf.loc[i, 'H'])\
                + ldf['Dispersion_cm_X'].abs())).min() for i in ldf.index])
        dP_Psep = np.where(daccL < dP_Psep, daccL, dP_Psep)

    Gs = lamRF/2/np.pi*q/nus*np.abs(etas)
    ldf['Sigma_um_X'] = lattice.get_sigma_um(ldf['Beta_cm_X'], ex,
                                             ldf['Dispersion_cm_X'], sp)
    em = dP_Psep**2*ldf['Beta_cm_X']/gamma**2/ex/1e-4\
        / (1+(sp*ldf['Beta_cm_X']*ldf["Phi_X"]/ldf["Sigma_um_X"]/1e-4)**2)
    Ne = Ibeam*1e-3/f0/e
    ldf['Sigma_um_Y'] = lattice.get_sigma_um(ldf['Beta_cm_Y'], 1,
                                            0, 0)
    CIem=0
    try:
        CIem = CI(em)
    except:
        print(em)
    sum_components = (ldf['Beta_cm_X']/ex/1e-4)**1.5*CIem\
        / (1+(sp*ldf['Beta_cm_X']*ldf['Phi_X']/ldf['Sigma_um_X']/1e-4)**2)**1.5\
        / (ldf['Sigma_um_X']*1e-4*1e-4*ldf['Sigma_um_Y']*em)*ldf['dS']/(1/f0)
    LamTska = Ne*re**2/(8*np.pi*gamma**5*sz)*sum_components.sum()
    return LamTska


def get_Touchek_Lifetime(lattice_df, Vrf, sp, ex, ey, sz, Ibeam,
                         aperture_factor=1.0,
                         gamma=gamma0):
    """The result assumes ey=1. So it has to be devided by sqrt(ey) in units of um"""
    V0 = Vrf
    ldf = lattice_df
    etas = alpha-1/gamma**2
    Ks = (alpha*gamma**2-1)/(gamma**2-1)
    phiacc = np.arcsin(VSR/V0)
    nus0 = np.sqrt(q*V0*np.abs(Ks)/2/np.pi/me/gamma)
    nus = nus0*np.sqrt(np.cos(phiacc))
    fs = f0*nus
    dP_Psep = aperture_factor*2*nus0/q/np.abs(etas)*np.sqrt(np.cos(phiacc)
                                            - (np.pi/2-phiacc)*np.sin(phiacc))
    Gs = lamRF/2/np.pi*q/nus*np.abs(etas)


    ldf['Sigma_um_X'] = lattice.get_sigma_um(ldf['Beta_cm_X'], ex,
                                            ldf['Dispersion_cm_X'], sp)


    aux1 = ldf['Beta_cm_X']/ex/1e-4/(
            1+(sp*ldf['Beta_cm_X']*ldf["Phi_X"]/ldf["Sigma_um_X"]/1e-4)**2)

    aux2 = ldf['Beta_cm_Y']/ey/1e-4

    B1 = 1/2/gamma**2*np.abs(aux1 + aux2)
    B2 = 1/2/gamma**2*np.abs(aux1 - aux2)

    S33 = 1/sp**2 + ldf['H']/ex/1e-4

    Ne = Ibeam*1e-3/f0/e

    um = dP_Psep**2
    Itski = np.array([Itsk(um, b1, b2)[0] for b1, b2 in zip(B1, B2)])

    sum_components = Itski/np.sqrt(S33)*ldf['dS']/(1/f0)
    LamTska = Ne*re**2/(8*np.sqrt(np.pi)*gamma**4*ex*1e-4*ey*1e-4*sz*sp)\
        * sum_components.sum()

    return 1/LamTska

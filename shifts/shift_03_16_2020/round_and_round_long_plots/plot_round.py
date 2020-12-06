#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import os
from config import get_from_config
import fur.path_assistant as path_assistant
shift = path_assistant.PathAssistant('shift_03_16_2020')
from wiggler_radiation.number_of_coherent_modes.coherent_modes import get_M_interpolator_at_fixed_energy
cur_to_sum_channel = get_from_config("Beam_current_to_Sum_channel_ampl_V/mA")
sum_channel_to_photoelectrons = get_from_config('sum_channel_to_photoelectrons')
N_to_I = 1/sum_channel_to_photoelectrons/cur_to_sum_channel
meas_ROUND = pd.read_csv(shift.get_results_dir().fi('meas_ROUND_03_16_2020.csv'), index_col=0)
meas_ROUND = meas_ROUND.sort_values(by='N',ignore_index=True)




# In[3]:


def f(x, alpha):
    return x+alpha*x**2
import scipy.optimize
alpha = scipy.optimize.curve_fit(f, meas_ROUND['N'], meas_ROUND['varN'])[0][0]
def NfromVarN(vn):
    return (-1+np.sqrt(1+4*alpha*vn))/2/alpha
NfromVarN(1.7e8)


# In[4]:


df = pd.read_csv("M_on_grid_precalc.csv", index_col=0)
df_EB = pd.read_csv("M_on_grid_precalc_EB.csv", index_col=0)
df_ET = pd.read_csv("M_on_grid_precalc_ET.csv", index_col=0)


# In[5]:


df.head()


# In[6]:


es = df.columns.values.astype(np.float64)


# In[7]:


es


# In[8]:


def reconstruct_e(avN, varN, en='0'):
    Mexp = avN**2/(varN-avN)
    df0 = {'0':df, 'B':df_EB, 'T':df_ET}[en]
    return np.interp(Mexp,df0.loc[avN,:],es)


# In[9]:


reconstruct_e(meas_ROUND['N'][0], meas_ROUND['varN'][0]+0.01e8)


# In[10]:


reconstruct_e(meas_ROUND['N'][0], meas_ROUND['varN'][0])


# In[11]:


meas_ROUND['e_rec'] = meas_ROUND.apply(
    lambda row: reconstruct_e(*row[['N', 'varN']]), axis=1)
meas_ROUND['e_rec_EB'] = meas_ROUND.apply(
    lambda row: reconstruct_e(*row[['N', 'varN']], en='B'), axis=1)
meas_ROUND['e_rec_ET'] = meas_ROUND.apply(
    lambda row: reconstruct_e(*row[['N', 'varN']], en='T'), axis=1)
meas_ROUND['varNtop'] = meas_ROUND['varN']+meas_ROUND['errorbar']
meas_ROUND['varNbottom'] = meas_ROUND['varN']-meas_ROUND['errorbar']
meas_ROUND['e_error_top'] = meas_ROUND.apply(
    lambda row: reconstruct_e(
        *row[['N', 'varNbottom']]), axis=1)\
        - meas_ROUND['e_rec']
meas_ROUND['e_error_bottom'] = -meas_ROUND.apply(
    lambda row: reconstruct_e(
        *row[['N', 'varNtop']]), axis=1)\
        + meas_ROUND['e_rec']


# # Systematic error of reconstructed e due to beam energy uncertainty

# In[12]:


# 1e3*meas_ROUND.loc[:,['e_rec_EB', 'e_rec_ET']].apply(lambda col:col-meas_ROUND.loc[:,'e_rec'],axis=0)


# In[13]:


meas_ROUND['I_mA'] = meas_ROUND['N']*N_to_I


# In[14]:


# calculation of reconstructed sy
import lattice.lattice as lattice
#from lattice.summary_in_undulator import CalcTransverseBeamParams
lattice_df =     lattice.read_lattice_file(shift.get_6dsim_dir()    .fi("IOTA_1NL_100MeV_v8.6.1.4.6ds_data.txt"))
from scipy.interpolate import interp1d
bp_df = pd.read_csv(shift.get_results_dir().fi("beam_params_vs_current_round.csv"))
ex_func = interp1d(bp_df["N"], bp_df['ex_um'],
                   bounds_error=False, fill_value="extrapolate")
dpp_func = interp1d(bp_df["N"], bp_df['dp/p'],
                    bounds_error=False, fill_value="extrapolate")
# def get_sy(row):
#     Sx, Sy, dx, dy, sxp, syp = CalcTransverseBeamParams(
#         lattice_df, ex_func(row['N']), row['ey_rec'], dpp_func(row['N']))
#     return np.sqrt(Sy**2+syp**2*dy**2)
# meas_FLAT["sy_rec"] = meas_FLAT.apply(get_sy, axis=1)


# In[ ]:





# In[15]:


# # save new meas_FLAT with ey and I_mA
# meas_FLAT.to_csv(shift.get_results_dir().fi('meas_FLAT_03_16_2020.csv'))


# In[16]:


round_df = pd.read_csv(shift.get_results_dir().fi("full_beam_params_vs_current_round.csv"))


# In[17]:


meas_FLAT = pd.read_csv(shift.get_results_dir().fi('meas_FLAT_03_16_2020.csv'), index_col=0)


# In[18]:


flat_df = pd.read_csv(shift.get_results_dir().fi("full_beam_params_vs_current_flat.csv"),
                     index_col=0)


# In[19]:


from scipy.interpolate import interp1d
bp_df = pd.read_csv(
    shift.get_results_dir().fi("beam_params_vs_current_round.csv"))
ex_r_func = interp1d(bp_df["N"], bp_df['ex_um'],
                   bounds_error=False, fill_value="extrapolate")
ey_r_func = interp1d(bp_df["N"], bp_df['ey_um'],
                   bounds_error=False, fill_value="extrapolate")


bp_df = pd.read_csv(
    shift.get_results_dir().fi("beam_params_vs_current_flat.csv"))
sx_f_func = interp1d(bp_df["N"], bp_df['Sigma_um_X'],
                   bounds_error=False, fill_value="extrapolate")
sy_f_func = interp1d(bp_df["N"], bp_df['Sigma_um_Y'],
                   bounds_error=False, fill_value="extrapolate")
sz_f_func = interp1d(bp_df["N"], bp_df['sz_um'],
                   bounds_error=False, fill_value="extrapolate")
ex_f_func = interp1d(bp_df["N"], bp_df['ex_um'],
                   bounds_error=False, fill_value="extrapolate")
ey_f_func = interp1d(bp_df["N"], bp_df['ey_um'],
                   bounds_error=False, fill_value="extrapolate")
dpp_f_func = interp1d(bp_df["N"], bp_df['dp/p'],
                    bounds_error=False, fill_value="extrapolate")
Vrf_f_func = interp1d(bp_df["N"], bp_df['N:IRFEPA'],
                    bounds_error=False, fill_value="extrapolate")


# In[20]:


def remove_outliers_omce(df, c1, l1):
    di = df[c1].diff().abs()
    return df[di < l1]
    
def remove_outliers(df, c1, l1, niter):
    for _ in range(niter):
        df = remove_outliers_omce(df, c1, l1)
    return df
round_df = remove_outliers(round_df,'N:IWCMI_recalibrated_to_IWCMI_absolute', 0.1, 10)
round_df = remove_outliers(round_df,'ex_um', 0.01, 10)
round_df = remove_outliers(round_df,'ey_um', 0.01, 10)

flat_df = remove_outliers(flat_df,'N:IWCMI_recalibrated_to_IWCMI_absolute', 0.1, 10)


# In[21]:


lt_f_df = pd.read_csv(shift.get_results_dir().fi("life_time_flat_03_16_2020.csv"),
            index_col=0)
lt_f_df


# In[22]:


et_df = pd.read_csv("e_from_touschek_round.csv", index_col=0)


# In[23]:



# prop = fm.FontProperties(fname='/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf')


# In[24]:


from fur.path_assistant import get_plot_style_sheet
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use("pgf")
plt.style.use(get_plot_style_sheet("prl"))

# In[100]:

elinewidth=2
markeredgewidth=2
bap=0.1
htp = 0.3
hl = 1
powerfs = 30
powerpos = 0.93, 1.38
subtitlepos = 0.55, 0.05
subfs = 28
avNlabel = r"Photoelectron count mean $\langle\mathcal{N}\rangle$"
beamCurLabel = r"Beam current (mA)"
hybrid_balance_error = get_from_config("HybridBalanceError")
fit_errorbar = get_from_config("varN_error_fit")
df = round_df

plt.rcParams.update({'figure.figsize':(30,8)})
fs = 32


i_to_photoel = cur_to_sum_channel*sum_channel_to_photoelectrons
df = pd.concat([
    pd.DataFrame({'N:IWCMI_recalibrated_to_IWCMI_absolute':[-meas_ROUND['N'].max()/i_to_photoel],
                 'ex_um': ex_r_func(meas_ROUND['N'].max()),
                 'ey_um': ey_r_func(meas_ROUND['N'].max())
                 }),
    df], ignore_index=True)


fig, ax = plt.subplots(2,2, gridspec_kw={'height_ratios': [2, 3]})


photoel = i_to_photoel*(-df['N:IWCMI_recalibrated_to_IWCMI_absolute'])
axRE = ax[1][0]
axRE.text(*subtitlepos, r"Round beam ($\epsilon_1=\epsilon_2=\epsilon$)", transform=axRE.transAxes,
      fontsize=subfs, va='bottom', ha='right')
axRE.plot(photoel/i_to_photoel, 0.5*1000*(df['ex_um']+df['ey_um']),'-', linewidth=5 ,
        label=r'$\epsilon$ via SLMs')
axRE.plot(et_df['I (mA)'], 1000*et_df['e'],'.-', color='tab:orange',
                 marker='o', linewidth=3, label=r'$\epsilon$ via Touschek lifetime',
                zorder=99)

axRE.set_xlabel(beamCurLabel, fontsize=fs)
axRE.set_ylabel(r"Emittance (nm)", fontsize=fs)

color='tab:red'
yerr = 1e3*np.array([meas_ROUND['e_error_bottom'],meas_ROUND['e_error_top']])
#[meas_FLAT['Sigma_um_Y_Meas_Bottom'], meas_FLAT['Sigma_um_Y_Meas_Top']]
axRE.errorbar(meas_ROUND['N']/i_to_photoel, 1e3*meas_ROUND['e_rec'],
             marker='o', linestyle='None',
             elinewidth=elinewidth,
              markeredgewidth=markeredgewidth,
             color=color, yerr=yerr, zorder=100,
             label=r"$\epsilon$ via fluctuations")
axRE.set_ylim(0, 1.05*axRE.get_ylim()[1])
handles, labels = axRE.get_legend_handles_labels()
order = [0,2,1]
# order = [0,1]
axRE.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right',
#            bbox_to_anchor=(1.01,-0.045),
            borderaxespad=bap,
            labelspacing=0.0, handlelength=hl,
           handletextpad=htp)




axRF = ax[0][0]
pw=1e7
hybrid_errorbar = 2*hybrid_balance_error    *np.absolute(meas_ROUND['varN'])
error_bar = np.sqrt(0*hybrid_errorbar**2+fit_errorbar**2)
axRF.errorbar(meas_ROUND['N']/pw, meas_ROUND['varN'],
            marker='.',linestyle='None',yerr=error_bar,color='b',
            label = r'Round beam fluctuations')
axRF.plot(meas_ROUND['N']/pw, meas_ROUND['N'],color='blue',linestyle='--',
        label = r"$\mathrm{var}\left(\mathcal{N}\right)=\langle\mathcal{N}\rangle$")
axRF.set_ylabel(r"$\mathrm{var}\left(\mathcal{N}\right)$", fontsize=fs)
axRF.set_xlabel(avNlabel, fontsize=fs)
handles, labels = axRF.get_legend_handles_labels()
order = [1,0]
axRF.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
            loc='upper left',borderaxespad=bap,
            labelspacing=0.0, handlelength=hl,
           handletextpad=htp)
axRF.set_xlim(np.asarray(axRE.get_xlim())*i_to_photoel/pw)
axRF.xaxis.tick_top()
axRF.xaxis.set_label_position('top') 
axRF.set_xticks(ticks=axRF.get_xticks()[2:-1])
axRF.text(*powerpos, r"$\times\SI{e7}{}$", transform=axRF.transAxes,
      fontsize=powerfs, fontweight='bold', va='top', ha='left')
axRF.get_yaxis().get_offset_text().set_x(-0.05)
axRF.xaxis.labelpad = 0




axFE = ax[1][1]
axFE.set_ylabel(r"$\epsilon_x$ (nm)", fontsize=fs)

df = flat_df


# extrapolation:
df = pd.concat([
    pd.DataFrame({'N:IWCMI_recalibrated_to_IWCMI_absolute':[-meas_FLAT['N'].max()/i_to_photoel],
                 'ex_um': ex_f_func(meas_FLAT['N'].max())}),
    df], ignore_index=True)
photoel = i_to_photoel*(-df['N:IWCMI_recalibrated_to_IWCMI_absolute'])
ln1 = axFE.plot(photoel/i_to_photoel, 1000*df['ex_um'],
                    '-', linewidth=5, label=r'$\epsilon_x$ via SLMs', zorder=99)
axFE.set_ylim(0, 1.05*axFE.get_ylim()[1])
axFE.text(*subtitlepos, r"Flat beam ($\epsilon_x\gg\epsilon_y$)", transform=axFE.transAxes,
      fontsize=subfs, va='bottom', ha='right')
axFE.set_xlabel(beamCurLabel, fontsize=fs)

axFEy = axFE.twinx()
color='tab:red'
yerr = 1e3*np.array([meas_FLAT['ey_error_bottom'],meas_FLAT['ey_error_top']])
ln2 = axFEy.errorbar(meas_FLAT['N']/i_to_photoel, 1e3*meas_FLAT['ey_rec'],
             marker='o', linestyle='None',
             color=color, yerr=yerr, label=r'$\epsilon_y$ via fluctuations',
                    zorder=100)

yrange = axFEy.get_ylim()

ln3 = axFEy.plot(lt_f_df['I (mA)'], lt_f_df['ey_rec_no_transv_acc'],
                '.-', color='tab:orange',
                 marker='o', linewidth=3, label=r'$\epsilon_y$ via Touschek lifetime',
                zorder=1)

axFEy.set_ylim(yrange)

axFEy.set_ylabel(r"$\epsilon_y$ (nm)", fontsize=fs)

lns = [ln1[0], ln2, ln3[0]]
labs = [l.get_label() for l in lns]
axFE.legend(lns, labs, loc='lower right',
#            bbox_to_anchor=(1.01,-0.045),
           borderaxespad=bap,
            labelspacing=0.0, handlelength=hl,
           handletextpad=htp)




axFF = ax[0][1]
pw=1e7
meas_FLAT =     pd.read_csv(shift.get_results_dir().fi('meas_FLAT_03_16_2020.csv'),
                index_col=0)
hybrid_errorbar = 2*hybrid_balance_error    *np.absolute(meas_FLAT['varN'])
error_bar = np.sqrt(0*hybrid_errorbar**2+fit_errorbar**2)
axFF.errorbar(meas_FLAT['N']/pw,meas_FLAT['varN'],
            marker='.',linestyle='None',yerr=error_bar,color='b',
            label = r'Flat beam fluctuations')
axFF.plot(meas_FLAT['N']/pw, meas_FLAT['N'],color='blue',linestyle='--',
        label = r"$\mathrm{var}\left(\mathcal{N}\right)=\langle\mathcal{N}\rangle$")
axFF.set_ylabel(r"$\mathrm{var}\left(\mathcal{N}\right)$", fontsize=fs)
axFF.set_xlabel(avNlabel, fontsize=fs)
handles, labels = axFF.get_legend_handles_labels()
order = [1,0]
axFF.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
            loc='upper left',borderaxespad=bap,
            labelspacing=0.0, handlelength=hl,
           handletextpad=htp)
axFF.set_xlim(np.asarray(axFE.get_xlim())*i_to_photoel/pw)
axFF.xaxis.tick_top()
axFF.xaxis.set_label_position('top')
axFF.set_xticks(ticks=np.arange(0.8,2.41,0.2))
axFF.text(*powerpos, r"$\times\SI{e7}{}$", transform=axFF.transAxes,
      fontsize=powerfs, fontweight='bold', va='top', ha='left')
axFF.get_yaxis().get_offset_text().set_x(-0.05)
axFF.xaxis.labelpad = 0


fig.subplots_adjust(hspace=0.15)
for i, label in enumerate(('(a)', '(b)', '(c)', '(d)')):
    x = i // 2
    y = i - 2*x
    axx = ax[x][y] 
    axx.text(-0.12, 1.15, label, transform=axx.transAxes,
      fontsize=36, va='top', ha='left')
print("Remember that a small protion of SLM data at large Ibeam is extrapolated")
plt.savefig(path_assistant.get_PRL_images_dir().fi("iota_measurements.png"),
            dpi=300, bbox_inches='tight')

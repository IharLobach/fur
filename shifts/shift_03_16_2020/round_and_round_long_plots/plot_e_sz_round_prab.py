#!/usr/bin/env python
# coding: utf-8

# In[1]:


from fur.path_assistant import get_plot_style_sheet
import scipy.optimize
import matplotlib.font_manager as fm
import matplotlib as mpl
from scipy.interpolate import interp1d
import lattice.lattice as lattice
from wiggler_radiation.number_of_coherent_modes.coherent_modes import get_M_interpolator_at_fixed_energy
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from config import get_from_config
import fur.path_assistant as path_assistant
shift = path_assistant.PathAssistant('shift_03_16_2020')
cur_to_sum_channel = get_from_config("Beam_current_to_Sum_channel_ampl_V/mA")
sum_channel_to_photoelectrons = get_from_config(
    'sum_channel_to_photoelectrons')
N_to_I = 1/sum_channel_to_photoelectrons/cur_to_sum_channel
meas_ROUND = pd.read_csv(shift.get_results_dir().fi(
    'meas_ROUND_03_16_2020.csv'), index_col=0)
meas_ROUND = meas_ROUND.sort_values(by='N', ignore_index=True)


# In[2]:


colors = {"FLAT": 'blue'}

fig, ax = plt.subplots(figsize=(15, 7.5))
ax.errorbar(meas_ROUND['N'], meas_ROUND['varN'], marker='o', linestyle='None',
            yerr=meas_ROUND['errorbar'], color='b',
            label=r'Measurements in IOTA')
x_aux = np.linspace(0, ax.get_xlim()[1], 100)
ax.plot(x_aux, x_aux, color='green', linestyle='--',
        label=r"var$\left(\mathcal{N}\right)=\langle\mathcal{N}\rangle$")
ax.set_ylabel(r"Photoelectron count variance var$\left(\mathcal{N}\right)$")
ax.set_xlabel(r"Photoelectron count mean $\langle\mathcal{N}\rangle$")
handles, labels = ax.get_legend_handles_labels()
order = [1, 0]
ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order])
ax.set_xlim(0, 1.05*meas_ROUND['N'].max())
ax.set_ylim(0, 1.1*meas_ROUND['varN'].max())
ax1 = ax.twiny()
ax1.set_xlabel('Beam current, \SI{}{mA}')
ax1.set_xlim(N_to_I*np.asarray(ax.get_xlim()))
ax1.set_xticks(ticks=ax1.get_xticks()[1:-1])
plt.show()


# In[3]:


def f(x, alpha):
    return x+alpha*x**2


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
    df0 = {'0': df, 'B': df_EB, 'T': df_ET}[en]
    return np.interp(Mexp, df0.loc[avN, :], es)


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
#from lattice.summary_in_undulator import CalcTransverseBeamParams
lattice_df = lattice.read_lattice_file(
    shift.get_6dsim_dir()    .fi("IOTA_1NL_100MeV_v8.6.1.4.6ds_data.txt"))
bp_df = pd.read_csv(shift.get_results_dir().fi(
    "beam_params_vs_current_round.csv"))
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


round_df = pd.read_csv(shift.get_results_dir().fi(
    "full_beam_params_vs_current_round.csv"))


# In[17]:


meas_FLAT = pd.read_csv(shift.get_results_dir().fi(
    'meas_FLAT_03_16_2020.csv'), index_col=0)


# In[18]:


flat_df = pd.read_csv(shift.get_results_dir().fi("full_beam_params_vs_current_flat.csv"),
                      index_col=0)


# In[19]:


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


round_df = remove_outliers(
    round_df, 'N:IWCMI_recalibrated_to_IWCMI_absolute', 0.1, 10)
round_df = remove_outliers(round_df, 'ex_um', 0.01, 10)
round_df = remove_outliers(round_df, 'ey_um', 0.01, 10)

flat_df = remove_outliers(
    flat_df, 'N:IWCMI_recalibrated_to_IWCMI_absolute', 0.1, 10)


# In[21]:


lt_f_df = pd.read_csv(shift.get_results_dir().fi("life_time_flat_03_16_2020.csv"),
                      index_col=0)
lt_f_df


# In[22]:


et_df = pd.read_csv("e_from_touschek_round.csv", index_col=0)


mpl.use("pgf")
plt.style.use(get_plot_style_sheet("prab"))


i_to_photoel = cur_to_sum_channel*sum_channel_to_photoelectrons



# In[73]:


fig, axs = plt.subplots(2, figsize=(15, 8))
axRE, ax2 = axs
df = round_df
df = pd.concat([
    pd.DataFrame({'N:IWCMI_recalibrated_to_IWCMI_absolute': [-meas_ROUND['N'].max()/i_to_photoel],
                  'ex_um': ex_r_func(meas_ROUND['N'].max()),
                  'ey_um': ey_r_func(meas_ROUND['N'].max())
                  }),
    df], ignore_index=True)
photoel = i_to_photoel*(-df['N:IWCMI_recalibrated_to_IWCMI_absolute'])
axRE.text(0.5, 0.4, r"Round beam (by design $\epsilon_1=\epsilon_2=\epsilon$)", transform=axRE.transAxes, va='bottom', ha='center')
axRE.plot(photoel/i_to_photoel, 0.5*1000*(df['ex_um']+df['ey_um']), '-', linewidth=5,
          label=r'$\epsilon$ via SLMs')
# axRE.plot(photoel/i_to_photoel, 1000*df['ey_um'],'-', linewidth=5 ,
#         label=r'$\epsilon_y$ via SLMs')


axRE.set_ylabel(r"Emittance $(\SI{}{nm})$")

axRE.set_ylim(0, 1.05*axRE.get_ylim()[1])
axRE.legend(loc='lower right')
axRE.get_xaxis().set_visible(False)

ax2.plot(photoel/i_to_photoel, df['N:IWCMBR'], '-', linewidth=5,
         label=r'RMS bunch length $\sigma_z$', color='tab:red')
ax2.plot(photoel/i_to_photoel, df['N:IWCMBE'], '-', linewidth=5,
         label=r'Effective bunch length $\sigma_z^{\mathrm{eff}}$',
         color='purple')
ax2.legend(loc='lower right')
ax2.set_ylabel("Bunch length $(\SI{}{cm})$")
ax2.set_yticks(ticks=np.arange(24, 31, 2))


ax2.set_xlabel(r"Beam current $(\SI{}{mA})$")
# ax2.set_ylim(0, 1.05*ax2.get_ylim()[1])
fig.subplots_adjust(hspace=0.1)
plt.savefig(path_assistant.get_PRL_images_dir().fi("ex_ey_round_beam.png"),
            dpi=300, bbox_inches='tight')
plt.show()


# In[ ]:

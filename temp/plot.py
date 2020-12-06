from fur.path_assistant import get_plot_style_sheet
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
fm._rebuild()
# -*- coding:utf-8 -*-

mpl.use("pgf")
plt.style.use(get_plot_style_sheet("prl"))


plt.plot([1, 2, 3])
plt.xlabel(
    r'Photoelectron count variance $\mathrm{var}(\mathcal{N})$, $\epsilon_x$, (\SI{}{nm})', fontsize=14)
plt.ylabel(r'Normal text ...  $\epsilon_x$', fontsize=28)
plt.savefig("plot.png",
            dpi=300, bbox_inches='tight')

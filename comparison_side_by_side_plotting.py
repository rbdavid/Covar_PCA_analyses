#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# USAGE:

# PREAMBLE:

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.ticker import NullFormatter, MultipleLocator, AutoMinorLocator, FixedLocator
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.gridspec as gridspec

nRes = int(sys.argv[1])
res_offset = int(sys.argv[2])	# zero indexed

first_residue = res_offset+1
last_residue = res_offset+nRes+1
print first_residue, last_residue

sys_list = []
sys_list.append(['Apo',21,150])
sys_list.append(['ATP',21,150])
sys_list.append(['ssRNA',21,200])
sys_list.append(['ssRNA+ATP',21,150])

res_highlight_list = [[(first_residue,last_residue),(first_residue,last_residue)],
		[(first_residue,315+res_offset+1),(first_residue,315+res_offset+1)],
		[(187,192),(187,192)],
		[(285,291),(285,291)],
		[(410,416),(410,416)],
		[(285,291),(410,416)]]	# start and stop residue numbers (including the  offset); just subtract the res_offset from these numbers to get the correct index 

my_cmap = 'bwr'
minorLocator = MultipleLocator(10)

# -----------------
# FUNCTIONS:
def draw_text(ax,string):
	from matplotlib.offsetbox import AnchoredText
	at = AnchoredText(string,loc=4, prop=dict(size=15), frameon=True)
	at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
	ax.add_artist(at)

# -----------------
# MAIN:
residue_range = np.array(range(first_residue,last_residue))		# resid 187 is residue 1 in the analysis and has indexing of zero in the data file...
print residue_range[0],residue_range[-1], len(residue_range)

fig_list = []
for i in res_highlight_list:
	temp = plt.figure(figsize=(12,9))
	temp_grid = ImageGrid(temp,111,nrows_ncols=(2,2),axes_pad=0.15,share_all=True,label_mode='L',cbar_location='right',cbar_mode='single',cbar_size='7%',cbar_pad=0.15)
	fig_list.append([temp,temp_grid])

for i in range(len(sys_list)):
	distance_correlation_data = np.loadtxt('%03d.%03d.%s.average_distance_correlation.dat'%(sys_list[i][1],sys_list[i][2],sys_list[i][0]),dtype=np.float64)
	for j in range(len(res_highlight_list)):
		temp = fig_list[j][1][i].pcolormesh(residue_range,residue_range,distance_correlation_data,cmap=my_cmap,vmin=-1.0,vmax=1.0)
		plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
		fig_list[j][1][i].set_xlim((res_highlight_list[j][0][0],res_highlight_list[j][0][1]))
		fig_list[j][1][i].set_ylim((res_highlight_list[j][1][0],res_highlight_list[j][1][1]))

		if j>0 and i == 0:
			xlabels = [str(int(x)) for x in temp.axes.get_xticks()[:-1]]
			ylabels = [str(int(y)) for y in temp.axes.get_yticks()[:-1]]
			
			temp.axes.set_xticks(temp.axes.get_xticks()[:-1]+0.5)
			temp.axes.set_yticks(temp.axes.get_yticks()[:-1]+0.5)
			temp.axes.set_xticklabels(xlabels)
			temp.axes.set_yticklabels(ylabels)
			temp.axes.set_aspect('equal')

		draw_text(fig_list[j][1][i],sys_list[i][0])
		if i == 0:
			fig_list[j].append(temp)

for i in range(len(res_highlight_list)):
	fig_list[i][1][-1].cax.colorbar(fig_list[i][2])
	fig_list[i][0].text(0.5,0.05,'Residue Number',ha='center',va='center',fontsize=20)
	fig_list[i][0].text(0.135,0.5,'Residue Number',ha='center',va='center',rotation='vertical',fontsize=20)
	fig_list[i][0].text(0.89,0.5,'Correlation',ha='center',va='center',rotation='vertical',fontsize=20)
	fig_list[i][0].savefig('Comparison.%03d.%03d.%03d.%03d.heatmap.png'%(res_highlight_list[i][0][0],res_highlight_list[i][0][1],res_highlight_list[i][1][0],res_highlight_list[i][1][1]),dpi=300,transparent=True)

plt.close()
sys.exit()


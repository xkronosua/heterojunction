#!/usr/bin/python
# _*_ coding: utf-8 _*_
#Wed, 02 Jan 2013 21:38:16 +0200 

from pylab import *
from scipy import interpolate

# Fermi-Dirac half
FD = loadtxt('fermi.txt')
def fermi_i(Y):
	x, y = FD[:,1], FD[:,0]
	return interpolate.interp1d(x, y)(Y)

################################################################################
fname = './data/my_tmp/PGS'
x1 = 0.35
x2 = 0.35
Na_Nd = 1.2 * 10**24	# m^-3
n = 2. * 10**24			# m^-3
p = n + Na_Nd
N1 = 1. * 10**24	# 0.1
P1 = 1. * 10**24
N3, P3 = N1, P1
T = 297.				# K
da1 = 0.15 * 10**-6		# m

Width = da1 * 3
k = 1.38 * 10**-23		# Дж/K
e_ = 1.6 * 10**-19		# Кл
kT = k * T / e_			# eV

# Ширини заборонених зон (eB)
Eg1 = 1.424 + 1.247 * x1 
Eg2 = 1.424
Eg3 = 1.424 + 1.247 * x2

# маси електронів та дірок
m0 = 9.11 * 10**-31		# kg
mp1 = (0.48 + 0.31*x1) * m0
mp2 =  0.48 * m0
mp3 = (0.48 + 0.31*x2) * m0

mnG1 = (0.067 + 0.083 * x1) * m0
mnG2 =  0.067 * m0
mnG3 = (0.067 + 0.083 * x2) * m0

# Ефективні густини станів
Nc1 = 4.83 * 10**21 * (mnG1 * T / m0)**(3./2.)
Nc2 = 4.83 * 10**21 * (mnG2 * T / m0)**(3./2.)
Nc3 = 4.83 * 10**21 * (mnG3 * T / m0)**(3./2.)

Nv1 = 4.83 * 10**21 * (mp1 * T / m0)**(3./2.)
Nv2 = 4.83 * 10**21 * (mp2 * T / m0)**(3./2.)
Nv3 = 4.83 * 10**21 * (mp3 * T / m0)**(3./2.)

# Відносні рівні Фермі
Fn_Ec1 = fermi_i(N1 / Nc1) * kT
Fn_Ec2 = fermi_i( n / Nc2) * kT
Fn_Ec3 = fermi_i(N3 / Nc3) * kT

Ev_Fp1 = -fermi_i(P1 / Nv1) * kT
Ev_Fp2 = -fermi_i( p / Nv2) * kT
Ev_Fp3 = -fermi_i(P3 / Nv3) * kT

# Спорідненості (eB)
hi1 = 4.07 - 1.1 * x1
hi2 = 4.07
hi3 = 4.07 - 1.1 * x2

# Контактні різниці потенціалів
V_D12 = hi1 + Eg1 + Ev_Fp1 - (hi2 + Eg2 + Ev_Fp2) 
V_D23 = hi2 - Fn_Ec2 - (hi3 - Fn_Ec3) 

# Діелектрична проникність
e1 = 13.18 - 3.12 * x1
e2 = 13.18
e3 = 13.18 - 3.12 * x2
# Діел. прон. вакууму
e0 = 8.85 * 10**-12		# (Кл/м)^2 * Н


N_N1 = 2.3 * 10**24
N_N2 = 1.1 * p
N_N3 = 3. * 10**24

# Контактні різниці потенціалів, для кожної з сторін переходу
fi12 = V_D12 / (1 + e1 * N_N1 / e2 / N_N2 )
fi21 = V_D12 / (1 + e2 * N_N2 / e1 / N_N1 )

fi23 = V_D23 / (1 + e2 * N_N2 / e3 / N_N3 )
fi32 = V_D23 / (1 + e3 * N_N3 / e2 / N_N2 )

# Ширини ОПЗ
d12 = sqrt(2 * e0 * e1 / e_ / N_N1 * fi12)
d21 = sqrt(2 * e0 * e2 / e_ / N_N2	* fi21)
d23 = sqrt(2 * e0 * e2 / e_ / N_N2	* fi23)
d32 = sqrt(2 * e0 * e3 / e_ / N_N3 * fi32)

# Розрахунок форми викривлень
d_12 = linspace(0, d12, 10)
fi_12 = e_ * N_N1 / (2. * e0 * e1) * (d_12)**2
d_21 = linspace(0, d21, 10)
fi_21 = e_ * N_N2 / (2. * e0 * e2) * (d_21 - d21)**2
d_23 = linspace(0, d23, 10)
fi_23 = e_ * N_N2 / (2. * e0 * e2) * (d_23)**2
d_32 = linspace(0, d32, 10)
fi_32 = e_ * N_N3 / (2. * e0 * e3) * (d_32 - d32)**2


# Значення відносно рівня вакууму
V_D23_v = V_D12 + V_D23

hi1_v = hi1
hi2_v = hi2 + V_D12
hi3_v = hi3 + V_D23_v
Ec1_v, Ec2_v, Ec3_v = hi1_v, hi2_v, hi3_v

Ev1_v = Eg1 + hi1_v
Ev2_v = Eg2 + hi2_v
Ev3_v = Eg3 + hi3_v

Fn_v = Ec2_v - Fn_Ec2
Fp_v = Ev2_v + Ev_Fp2

# розриви
dEc1 = abs(hi2-hi1)
dEc2 = abs(hi3-hi2)

dEv1 = abs(Eg2 - Eg1) - dEc1
dEv2 = abs(Eg3 - Eg2) - dEc2

# повні ширини ОПЗ
d_OPZ1 = d12 + d21
d_OPZ2 = d23 + d32
################################################################################


d_2 = da1 / 2.
W_2 = Width / 2.
Height = Ev3_v
def Plot_hetero():
	plot(
		[-d_2]*2,[0.1, -Height * 1.2],'-.k',
		[d_2]*2,[0.1, -Height * 1.2],'-.k')
	def plot_base(Add=[], style=''):
		plot(
			[-W_2, -(d_2 +d12)], [-Add[0]] * 2, style,
			-(d_2 + d12 - d_12), -(Add[0] + fi_12),style,
			-(d_2 - d_21), -(-fi_21 + V_D12 + Add[1]), style,
			[-d_2] * 2, [-(Add[0] + fi_12).max(),
				-(-fi_21 + V_D12 + Add[1]).min()], style,
			[-(d_2 - d21), (d_2 - d23)], [-V_D12 - Add[1]]*2,style,
			[d_2] * 2, [-(V_D12 + Add[1] + fi_23).max(),
				-(-fi_32 + V_D12 + V_D23 + Add[2]).min()], style,
			
			(d_2 - d23 + d_23), -(V_D12 + fi_23 + Add[1]),style,
			(d_2 + d_32), -(V_D23 + V_D12- fi_32 + Add[2]),style,
			[(d_2 + d32), W_2], [-V_D23 - V_D12 - Add[2]]*2,style,
			)

	# V_D
	plot([-W_2, W_2], [0]*2, 'k', alpha=0.4)
	plot_base(Add=[0,0,0], style='c')

	# Ec, Ev
	plot_base(Add=[hi1, hi2, hi3], style='b')
	plot_base(Add=[hi1 + Eg1, hi2 + Eg2, hi3 + Eg3], style='r')

	# Рівні Фермі

	f_x = linspace(d_2, W_2, 8)
	y_x = ((f_x - d_2) / W_2)**2. 
	plot(
		-f_x, -(Fn_v + y_x),'--g',
		[-d_2, d_2], [-(Fn_v)] * 2,'--g',
		[d_2, W_2], [-(Fn_v)] * 2,'--g')

	plot(
		[-d_2, d_2], [-(Fp_v)] * 2,'--g',
		[-d_2, -W_2], [-( Fp_v)] * 2,'--g',
		f_x, -(Fp_v - y_x),'--g')
		
################################################################################
def Arrows(X=[], Y=[], xy=[], xytext=[], color='k', Text='',
			val=0, Round=3, offset=[0, 0], style="<|-|>"):
	if len(xy) == 4:
		X, Y = xy[:2], xy[2:]
	annotate("",
			xy=(X[0], Y[0]), xycoords='data',
			xytext=(X[1], Y[1]), textcoords='data',
			arrowprops=dict(arrowstyle=style,
							connectionstyle="arc3"),
			)
	# these are matplotlib.patch.Patch properies
	props = dict(boxstyle='square,pad=0.', facecolor='w', alpha=0.5, )
	if Text:
		if len(xytext) != 2:
			xytext = [(X[1] + X[0])/2., (Y[1] + Y[0])/2.]
		if val:
			Text += " = %.3f"%(val) 
		text(xytext[0] + offset[0], xytext[1] + offset[1], Text, bbox=props)
			
import matplotlib.gridspec as gridspec
fig = figure(num=0, dpi=80, figsize=(15,8))
conv=1/2.58
fig.set_size_inches(28 * conv,20 * conv)
gs = gridspec.GridSpec(4, 4)

params = {
		'axes.labelsize': 15,
		'text.fontsize': 25,
		'xtick.labelsize': 16,
		'ytick.labelsize': 16,}
rcParams.update(params)	
rc("text", fontsize=15)
rc("lines", linewidth=2)

# 0
ax0 = plt.subplot(gs[:, :3])
Plot_hetero()
xlabel(r'$distance, m$', ha='left')
ylabel(r'$E, eV$', va='bottom')
xl, yl = ax0.get_xlim(), ax0.get_ylim()
Ofs = abs(array([xl[1]-xl[0], yl[1]-yl[0]]))/100.
text(d_2*4./3, 0, r'Vacuum level', bbox=dict(boxstyle='round',
 facecolor='w', alpha=0.9))
xlim(-W_2, W_2)
ylim( -(Height + Ofs[1]*5), Ofs[1]*3)

al = 7 # arrow length in points
arrowprops=dict(clip_on=False, # plotting outside axes on purpose
		frac=1., # make end arrowhead the whole size of arrow
		headwidth=al, # in points
		facecolor='k')
kwargs = dict(	
					xycoords='axes fraction',
					textcoords='offset points',
					arrowprops= arrowprops,
				 )

annotate("",(1,0),xytext=(-al,0), **kwargs) # bottom spine arrow
annotate("",(0,1),xytext=(0,-al), **kwargs) # left spin arrow
ax0.yaxis.tick_left()
ax0.xaxis.tick_bottom()
"""""""""""""""""""""" Підписи """""""""""""""""""""""""""
# hi
Arrows(X=[-W_2 + Ofs[0] * 3.]*2, Y=[0, -hi1_v], Text=r"$\chi_1$", val=hi1,
	offset=[Ofs[0],0])
Arrows(X=[-d_2 + Ofs[0] * 3.]*2, Y=[-V_D12, -hi2_v], Text=r"$\chi_2$", val=hi2,
	offset=[Ofs[0],0])
Arrows(X=[d_2 + Ofs[0] * 4.]*2, Y=[-V_D23_v, -hi3_v], Text=r"$\chi_3$", val=hi3,
	offset=[Ofs[0],0])
		
# Eg
		
Arrows(X=[-W_2 + Ofs[0] * 3.]*2, Y=[-Ec1_v, -Ev1_v], Text=r"$E_{g1}$", val=Eg1,
	offset=[Ofs[0],0])
Arrows(X=[-d_2 + Ofs[0] * 3.]*2, Y=[-Ec2_v, -Ev2_v], Text=r"$E_{g2}$", val=Eg2,
	offset=[Ofs[0],0])
Arrows(X=[d_2 + Ofs[0] * 4.]*2, Y=[-Ec3_v, -Ev3_v], Text=r"$E_{g3}$", val=Eg3,
	offset=[Ofs[0],0])
		
# da
Arrows(X=[-d_2, d_2], Y=[-(Ev3_v + Ofs[1])]*2,
		xytext=[-d_2 + Ofs[0]*3, -(Ev3_v + Ofs[1]*4)],
		Text=r"$d_{a} = %.2e$"%(da1))


# V_D
props = dict(boxstyle='square,pad=0.', facecolor='w', alpha=0.5, )
annotate(r"$V_{D12} = %.3f$"%V_D12,
			(d_2 - Ofs[0]*3, -V_D12/2.),
			(-d_2 + Ofs[0]*5, -V_D12/2. - Ofs[1]*5),
			arrowprops=dict(arrowstyle="-|>",
							connectionstyle="angle3, angleA=90, angleB=0"),
			bbox=props)
annotate(r"$V_{D23} + V_{D12} = %.3f$"%V_D23_v,
			(d_2 + Ofs[0]*5, -V_D23_v/2.),
			(d_2 + Ofs[0]*10, -V_D23_v/2. - Ofs[1]*5),
			arrowprops=dict(arrowstyle="-[",
							connectionstyle="angle3, angleA=90, angleB=0"),
			bbox=props)

# Рівні Фермі
annotate(r"${F_{n2}}^\ast-Ec_2 = %.3f$"%Fn_Ec2,
			(- Ofs[0]*3, -Fn_v - Fn_Ec2/2.),
			( -Ofs[0]*4, -Fn_v - Fn_Ec2/2. + Ofs[1]*10	),
			arrowprops=dict(arrowstyle="-|>",
							connectionstyle="angle3, angleA=90, angleB=0"),
			bbox=props)
annotate(r"$Ev_2-{F_{p2}}^\ast = %.3f$"%Ev_Fp2,
			(- Ofs[0]*3, -Fp_v + Ev_Fp2/2.),
			( -Ofs[0]*4, -Fp_v + Ev_Fp2/2. + Ofs[1]*4	),
			arrowprops=dict(arrowstyle="-|>",
							connectionstyle="angle3, angleA=90, angleB=0"),
			bbox=props)

# Підписи зон
props_b = dict(boxstyle='square,pad=0.', facecolor='b', alpha=0.2, )
props_r = dict(boxstyle='square,pad=0.', facecolor='r', alpha=0.2, )
text(-d_2-Ofs[0]*10, -Ec1_v+Ofs[1]*1.5,r'$E_{c1}$', bbox=props_b)
text(-d_2-Ofs[0]*10, -Ev1_v+Ofs[1]*1.5,r'$E_{v1}$', bbox=props_r)

text(Ofs[0]*5, -Ec2_v+Ofs[1]*1.5,r'$E_{c2}$', bbox=props_b)
text(Ofs[0]*5, -Ev2_v-Ofs[1]*2,r'$E_{v2}$', bbox=props_r)

text(d_2+Ofs[0]*10, -Ec3_v+Ofs[1]*1.5,r'$E_{c3}$', bbox=props_b)
text(d_2+Ofs[0]*10, -Ev3_v+Ofs[1]*1.5,r'$E_{v3}$', bbox=props_r)

# P-p-N
props = dict(boxstyle='round', facecolor='r', alpha=0.5, )
text(-(W_2 + d_2)/2, -hi1_v/3, r'$P$', bbox=props)
text(0, -hi1_v/3, r'$p$', bbox=props)
text((W_2 + d_2)/2, -hi1_v/3, r'$N$', bbox=props)

# c1 ===========================================================================
ax1 = plt.subplot(gs[0, 3])
title(r"$View_{c1}$")
def hide_axis(ax,x,y):
	ax.get_xaxis().set_visible(x)
	ax.get_yaxis().set_visible(y)
hide_axis(ax1, 0, 0)
xlim((-(d_2 + d12)*1.05, -(d_2 - d21)*0.9))
ylim((-(Ec2_v + Ofs[1]/2.), -(Ec1_v - Ofs[1]/6.)))
Plot_hetero()
yl=ax1.get_ylim()
plot([-d_2-d12]*2, yl, 'm', alpha=0.2)
plot([-d_2+d21]*2, yl, 'm', alpha=0.2)
Arrows(X=[-d_2-d12, -d_2], Y=[(yl[1] + yl[0])/2]*2, Text=r'$d_{12}$',
		offset=(0,Ofs[1]/3))
Arrows(X=[-d_2, -d_2+d21], Y=[(yl[1] + yl[0])/2]*2, Text=r'$d_{21}$',
		offset=(0,Ofs[1]/3))
delta = -Ec1_v - fi_12.max() - ( -Ec2_v + fi_21.max())
Arrows(X=[-d_2+d21 + Ofs[0]/3]*2, Y=[-Ec1_v - fi_12.max(), -Ec2_v + fi_21.max()],
	Text=r'$\Delta E_{c1}$', style='|-|', val=abs(delta))

# v1 ===========================================================================
ax2 = plt.subplot(gs[1, 3])
title(r"$View_{v1}$")
hide_axis(ax2, 0, 0)
xlim((-(d_2 + d12)*1.05, -(d_2 - d21)*0.9))
ylim((-(Ev2_v + Ofs[1]/2.), -(Ev1_v - Ofs[1]/6.)))
Plot_hetero()
xl, yl = ax2.get_xlim(), ax2.get_ylim()

plot([-d_2-d12]*2, yl, 'm', alpha=0.2)
plot([-d_2+d21]*2, yl, 'm', alpha=0.2)
delta = -Ev1_v - fi_12.max() - ( -Ev2_v + fi_21.max())
Arrows(X=[-d_2+d21/2]*2, Y=[-Ev1_v - fi_12.max(), -Ev2_v + fi_21.max()],
	Text=r'$\Delta E_{v1}$', style='|-|', val=abs(delta), offset=[0,-Ofs[1]/10])

# c2 ===========================================================================
ax3 = plt.subplot(gs[2, 3])
title(r"$View_{c2}$")
hide_axis(ax3, 0, 0)
xlim(((d_2 - d23)*0.9, (d_2 + d32)*1.1))
ylim((-(Ec2_v + Ofs[1]*2), -(Ec3_v - Ofs[1]*3)))
Plot_hetero()
yl=ax3.get_ylim()
plot([d_2-d23]*2, yl, 'm', alpha=0.2)
plot([d_2+d32]*2, yl, 'm', alpha=0.2)

delta = -Ec2_v - fi_23.max() - ( -Ec3_v + fi_32.max())
Arrows(X=[d_2+d32/3]*2, Y=[-Ec2_v - fi_23.max(), -Ec3_v + fi_32.max()],
	Text=r'$\Delta E_{c2}$', style='|-|', val=abs(delta), offset=[0,-Ofs[1]/10])
		
# v2 ===========================================================================
ax4 = plt.subplot(gs[3, 3])
title(r"$View_{v2}$")
hide_axis(ax4, 0, 0)
xlim(((d_2 - d23)*0.9, (d_2 + d32)*1.1))
ylim(	-(Ev3_v + Ofs[1]), -(Ev2_v - Ofs[1]*2))
Plot_hetero()
yl=ax4.get_ylim()
plot([d_2-d23]*2, yl, 'm', alpha=0.2)
plot([d_2+d32]*2, yl, 'm', alpha=0.2)
Arrows(X=[d_2-d23, d_2], Y=[(yl[1] + yl[0])/2 - Ofs[1]*2]*2, Text=r'$d_{23}$',
		offset=(0,Ofs[1]/3))
Arrows(X=[d_2, d_2+d32], Y=[(yl[1] + yl[0])/2 - Ofs[1]*2]*2, Text=r'$d_{32}$',
		offset=(0,Ofs[1]/3))		


delta = -Ev2_v - fi_23.max() - ( -Ev3_v + fi_32.max())
#text(X=[d_2+d32/3]*2, Y=[-Ev2_v - fi_23.max(), -Ev3_v + fi_32.max()],
#	Text=r'$\Delta E_{v2}$', style='|-|', val=abs(delta), offset=[0,-Ofs[1]/10])
props = dict(boxstyle='square,pad=0.', facecolor='w', alpha=0.5, )
annotate(r"$\Delta E_{v2} = %.3f$"%delta,
			(d_2, -Ev2_v - fi_23.max() - delta/2.),
			( d_2 + d32/3., -Fp_v),
			arrowprops=dict(arrowstyle="-|>",
							connectionstyle="angle3, angleA=90, angleB=0"),
			bbox=props)

savefig(fname + '.pdf', dpi=300, facecolor='w', edgecolor='w',
	orientation='landscape', papertype='a4', format='pdf',
	transparent=False, bbox_inches=None, pad_inches=0.1)

savefig(fname + '.png', dpi=300, facecolor='w', edgecolor='w',
	orientation='landscape', papertype='a4', format='png',
	transparent=False, bbox_inches=None, pad_inches=0.1)	
		
show()


output = [
	['Ec1_v' ,  Ec1_v],
	['Ec2_v' ,  Ec2_v],
	['Ec3_v' ,  Ec3_v ],
	['Eg1' ,  Eg1 ],
	['Eg2' ,  Eg2 ],
	['Eg3' ,  Eg3 ],
	['Ev1_v' ,  Ev1_v],
	['Ev2_v' ,  Ev2_v],
	['Ev3_v' ,  Ev3_v ],
	['Na_Nd' ,  Na_Nd],
	['Nc1' ,  Nc1 ],
	['Nc2' ,  Nc2 ],
	['Nc3' ,  Nc3 ],
	['Nv1' ,  Nv1 ],
	['Nv2' ,  Nv2 ],
	['Nv3' ,  Nv3 ],
	['T' ,  T ],
	['Width' ,  Width ],
	['conv ' ,  conv  ],
	['da1' ,  da1 ],
	['e0' ,  e0 ],
	['e1' ,  e1 ],
	['e2' ,  e2 ],
	['e3' ,  e3 ],
	['e_' ,  e_],
	['hi1' ,  hi1 ],
	['hi2' ,  hi2 ],
	['hi3' ,  hi3 ],
	['k' ,  k ],
	['kT' ,  kT ],
	['m0' ,  m0 ],
	['mnG1' ,  mnG1 ],
	['mnG2' ,  mnG2],
	['mnG3' ,  mnG3 ],
	['mp1' ,  mp1 ],
	['mp2' ,  mp2 ],
	['mp3' ,  mp3 ],
	['n' ,  n ],
	['p' ,  p  ],
	['x1' ,  x1 ],
	['x2' ,  x2],
	['Fn_Ec1' ,  Fn_Ec1],
	['Fn_Ec2' ,  Fn_Ec2],
	['Fn_Ec3' ,  Fn_Ec3],
	['Ev_Fp1' ,  Ev_Fp1],
	['Ev_Fp2' ,  Ev_Fp2],
	['Ev_Fp3' ,  Ev_Fp3],
	['V_D12' ,  V_D12],
	['V_D23' ,  V_D23],
	['fi12' ,  fi12],
	['fi21' ,  fi21],
	['fi23' ,  fi23],
	['fi32' ,  fi32],
	['d12' ,  d12],
	['d21' ,  d21],
	['d23' ,  d23],
	['d32' ,  d32],
	['Fn_v' ,  Fn_v],
	['Fp_v' ,  Fp_v],
	['dEc1' ,  dEc1],
	['dEc2' ,  dEc2],
	['dEv1' ,  dEv1],
	['dEv2' ,  dEv2],
	['d_OPZ1' ,  d_OPZ1],
	['d_OPZ2' ,  d_OPZ2]
]

import csv
with open(fname + '.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in output:
    	print '\t', i[0],'\t\t=\t', i[1]
    	spamwriter.writerow(i)

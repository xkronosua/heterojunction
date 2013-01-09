#!/usr/bin/python
# _*_ coding: utf-8 _*_
# чт, 03 січ 2013 16:16:36 +0200 

from pylab import *
from scipy import interpolate

# Fermi-Dirac half
FD = loadtxt('fermi.txt')
def fermi_i(Y):
	x, y = FD[:,1], FD[:,0]
	return interpolate.interp1d(x, y)(Y)

################################################################################
V = 45. * 10**-3		# eV
da = 10. * 10**-9		# m

fname = './data/my_tmp/RO_PGS'
x1 = 0.35
x2 = 0.35

#-------------------------------------------------------------------------------
x_Al_Ga_As = 2. * V / 1.247


#-------------------------------------------------------------------------------

Na_Nd = 1.2 * 10**24	# m^-3
n = 2. * 10**24			# m^-3
p = n + Na_Nd
N1 = 1. * 10**24
P1 = 1. * 10**24
N3, P3 = N1, P1
T = 297.				# K
da1 = 0.15 * 10**-6		# m

Width = da1 * 3.
k = 1.38 * 10**-23		# Дж/K
e_ = 1.6 * 10**-19		# Кл
kT = k * T / e_			# eV

# Ширини заборонених зон (eB)
Eg1 = 1.424 + 1.247 * x1
Eg2 = 1.424 + 1.247 * x_Al_Ga_As
Eg3 = 1.424 + 1.247 * x2

# маси електронів та дірок
m0 = 9.11 * 10**-31		# kg
mp1 = (0.48 + 0.31*x1) * m0
mp2 = (0.48 + 0.31*x_Al_Ga_As) * m0
mp3 = (0.48 + 0.31*x2) * m0

mnG1 = (0.067 + 0.083 * x1) * m0
mnG2 = (0.067 + 0.083 * x_Al_Ga_As) * m0
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
hi2 = 4.07 - 1.1 * x_Al_Ga_As
hi3 = 4.07 - 1.1 * x2

# Контактні різниці потенціалів
V_D12 = hi1 + Eg1 + Ev_Fp1 - (hi2 + Eg2 + Ev_Fp2) 
V_D23 = hi2 - Fn_Ec2 - (hi3 - Fn_Ec3) 

# Діелектрична проникність
e1 = 13.18 - 3.12 * x1
e2 = 13.18 - 3.12 * x_Al_Ga_As
e3 = 13.18 - 3.12 * x2
# Діел. прон. вакууму
e0 = 8.85 * 10**-12		# (Кл/м)^2 * Н

N_N1 = 2.3 * 10**24
N_N2 = 2.3 * 10**24
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
Ec0_v = Ec2_v + V

Ev1_v = Eg1 + hi1_v
Ev2_v = Eg2 + hi2_v
Ev3_v = Eg3 + hi3_v
Ev0_v = Ev2_v - V

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

"""++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"""
#h = 6.63 * 10**-34	# Дж*с
h_ = 6.58 * 10**-16		# eV*c
m0_eV = 0.51 * 10**6	# eV*C^2
c = 3 * 10**8			# m/c
################################################################################
################################################################################
################################################################################
# Спектр
def En(n, m, V0=V, step=10**-6, plot=False):
	m *= 5.6095 * 10**35 / c**2		# eV
	out = []
	for E in arange(0, V0, step):
		E_i = (h_ * pi)**2 / 2. / m / da**2 *\
			(n - 2. / pi * arctan(sqrt(E / (V0 - E))))**2
		if plot: out.append([E, E_i])
		if (E - E_i)**2 <= step: break
	if plot: return E, array(out)
	else: return E
	
# Хвильові ф-ції
def psi(n, m, V0=V, da=da, step=10**-4, N = 150.):
	E = En(n, m, V0=V0,step=step)
	m *= 5.6095 * 10**35 / c**2		# eV
	K = sqrt(2 * m * (V0 - E) / h_**2)
	a = sqrt(2 * m * E / h_**2)
	A2 = 1. / sqrt(da/2. + 1. / K/2.)/2
	def f(n):
		if n%2 == 0:
			return sin, -1.
		else:
			return cos, 1.
	x = [ linspace(-da, -da/2., N), linspace(-da/2., da/2, N), linspace(da/2., da, N)]
	F, t = f(n)
	y = [
		t * 2 * A2 * exp(K*(da/2.+x[0])) * F(a*da/2.),
		2 * A2 * F(a * x[1]),
		2 * A2 * exp(K*(da/2.-x[2])) * F(a*da/2.) ]
	x = array(x[0].tolist() + x[1].tolist() + x[2].tolist())
	y = array(y[0].tolist() + y[1].tolist() + y[2].tolist())
	return x, y, E, K, a, A2
################################################################################


d_2 = da1 / 2.
W_2 = Width / 2.
Height = Ev3_v

p_well_e_ = []
p_well_h = []

step = 10**-5
for i in range(1,50):
	tmp = psi(i, 0.067 * m0, step=step)
	if (V - tmp[2])**2 <= step:
		break
	p_well_e_.append(tmp)

for i in range(1,50):
	tmp = psi(i, 0.48 * m0,step=step)
	if (V - tmp[2])**2 <= step:
		break
	p_well_h.append(tmp)


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
	
	# Ями
	plot(
		[-da/2., da/2.], [-Ec2_v]*2,'w',
		[-da/2., da/2.], [-Ec2_v - V]*2,'b',
		[-da/2.]*2, [-Ec2_v - V, -Ec2_v],'b',
		[da/2.]*2, [-Ec2_v - V, -Ec2_v],'b',
		
		[-da/2., da/2.], [-Ev2_v]*2,'w',
		[-da/2., da/2.], [-Ev2_v + V]*2,'r',
		[-da/2.]*2, [-Ev2_v, -Ev2_v + V],'r',
		[da/2.]*2, [-Ev2_v, -Ev2_v + V],'r')
	plot([-da/2]*2,[0, -Height],'m', alpha=0.4)
	plot([da/2]*2,[0, -Height],'m', alpha=0.4)
	
	
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
fig.set_size_inches(29 * conv,21 * conv)
gs = gridspec.GridSpec(4, 4)

params = {
		'axes.labelsize': 15,
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

"""""""""""""""""""""" Підписи """""""""""""""""""""""""""""""""""""""""""""""
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

Arrows(X=[da/2*1.5 ]*2, Y=[-Ec0_v, -Ev0_v], Text=r"$E_{g0}$",
	val=(Eg2 - 2. * V),
	offset=[Ofs[0],(Eg2 - 2. * V)/3.])


# da1
Arrows(X=[-d_2, d_2], Y=[-(Ev3_v + Ofs[1])]*2,
		xytext=[-d_2 + Ofs[0]*3, -(Ev3_v + Ofs[1]*4)],
		Text=r"$d_{a1} = %.2e$"%(da1))
# da
'''
Arrows(X=[-da/2, da/2], Y=[-(hi2_v*2/3 + Ofs[1])]*2,
		Text=r"$d_{a} = %.2e$"%(da),
		offset=[da])
'''
props = dict(boxstyle='square,pad=0.', facecolor='w', alpha=0.5, )
annotate(r"$d_{a} = %.3e$"%da,
			(0, -hi2_v*2/3),
			(da/2 + Ofs[0], -hi2_v/2.7),
			arrowprops=dict(arrowstyle="-[",
							connectionstyle="angle3, angleA=10, angleB=90"),
			bbox=props)

# V_D
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
props_R = dict(boxstyle='round', facecolor='r', alpha=0.5, )
text(-(W_2 + d_2)/2, -hi1_v/3, r'$P$', bbox=props_R)
text(0, -hi1_v/3, r'$p$', bbox=props_R)
text((W_2 + d_2)/2, -hi1_v/3, r'$N$', bbox=props_R)

def hide_axis(ax,x,y):
	ax.get_xaxis().set_visible(x)
	ax.get_yaxis().set_visible(y)
	
# Ями ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def plot_well(p, xl):
	# Ями
	if not p:
		plot(
			[xl[0], -da/2.], [-Ec2_v]*2,'b',
			[da/2., xl[1]], [-Ec2_v]*2,'b',
			[-da/2., da/2.], [-Ec2_v - V]*2,'b',
			[-da/2.]*2, [-Ec2_v - V, -Ec2_v],'b',
			[da/2.]*2, [-Ec2_v - V, -Ec2_v],'b'
		)
	else:
		plot(
			[xl[0], -da/2.], [-Ev2_v]*2,'r',
			[da/2., xl[1]], [-Ev2_v]*2,'r',
			[-da/2., da/2.], [-Ev2_v + V]*2,'r',
			[-da/2.]*2, [-Ev2_v, -Ev2_v + V],'r',
			[da/2.]*2, [-Ev2_v, -Ev2_v + V],'r'
		)
# рівні Фермі в ямах		
n_2D = n * da
p_2D = p * da
En1_Fn =  kT * log(2. * exp(n_2D * h_**2 * pi / (m0_eV * 0.067 / c**2) / kT) - 1.)
Fp_Ep1 = -kT * log(2. * exp(p_2D * h_**2 * pi / (m0_eV * 0.48  / c**2) / kT) - 1.)

# 2De_ =========================================================================
ax1 = plt.subplot(gs[:2, 3])
title(r"$2D_{e^-}$")
ax1.yaxis.tick_right()
hide_axis(ax1, 0, 1)
xlim((-(da/2. + Ofs[0]), (da/2. + Ofs[0])))
ylim((-(Ec2_v + V + Ofs[1]/10),
	-min((Ec2_v + V - p_well_e_[0][2] - En1_Fn - Ofs[1]/10),
		Ec2_v - Ofs[1]/10)	))
xl, yl = ax1.get_xlim(), ax1.get_ylim()
plot_well(0, xl)
plot(xl, [-(Ec2_v + V - p_well_e_[0][2])]*2,'-.m')
plot(xl, [-(Ec2_v + V - p_well_e_[0][2] - En1_Fn)]*2,'--g')

Arrows(X=[ -da/3]*2, Y=[-(Ec2_v + V), -(Ec2_v + V - p_well_e_[0][2])],
	Text=r"$E_{n1}$", val=p_well_e_[0][2], offset=[Ofs[0]/5, 0])

Arrows(X=[ xl[0] + Ofs[0]/3]*2,
	Y=[-(Ec2_v + V - p_well_e_[0][2]), -(Ec2_v + V - p_well_e_[0][2] - En1_Fn)],
	Text=r"${F_n}^\ast - E_{n1}$", val=abs(En1_Fn))


# 2Dh ==========================================================================
ax2 = plt.subplot(gs[2:, 3])
title(r"$2D_{h}$")
ax2.yaxis.tick_right()
hide_axis(ax2, 0, 1)
xlim((-(da/2. + Ofs[0]), (da/2. + Ofs[0])))
ylim((-(Ev2_v + Ofs[1]/10), -(Ev2_v - V - Ofs[1]/10)))
xl, yl = ax2.get_xlim(), ax2.get_ylim()
plot_well(1, xl)
plot(xl, [-(Ev2_v - V + p_well_h[0][2])]*2,'-.m')
plot(xl, [-(Ev2_v - V + p_well_h[0][2] - Fp_Ep1)]*2,'--g')

Arrows(X=[ -da/3]*2, Y=[ -(Ev2_v - V + p_well_h[0][2]), -(Ev2_v - V)],
	Text=r"$E_{p1}$", val=p_well_h[0][2], offset=[Ofs[0]/5, 0])

Arrows(X=[ xl[0] + Ofs[0]/3]*2,
	Y=[-(Ev2_v - V + p_well_h[0][2] - Fp_Ep1 ), -(Ev2_v - V + p_well_h[0][2])],
	Text=r"$ E_{p1} - {F_p}^\ast$", val=abs(Fp_Ep1))

savefig(fname + '.pdf', dpi=300, facecolor='w', edgecolor='w',
	orientation='landscape', papertype='a4', format='pdf',
	transparent=False, bbox_inches=None, pad_inches=0.1)		

savefig(fname + '.png', dpi=300, facecolor='w', edgecolor='w',
	orientation='landscape', papertype='a4', format='png',
	transparent=False, bbox_inches=None, pad_inches=0.1)

fig1 = figure(num=1, dpi=80, figsize=(15,8))
conv=1/2.58
fig.set_size_inches(28 * conv,20 * conv)

# 2De_
ax0 = subplot(121)
annotate("",(1,0),xytext=(-al,0), **kwargs) # bottom spine arrow
annotate("",(0,1),xytext=(0,-al), **kwargs) # left spin arrow
plot(
	[-da, -da/2.], [-Ec2_v]*2,'b',
	[da/2., da], [-Ec2_v]*2,'b',
	[-da/2., da/2.], [-Ec2_v - V]*2,'b',
	[-da/2.]*2, [-Ec2_v - V, -Ec2_v],'b',
	[da/2.]*2, [-Ec2_v - V, -Ec2_v],'b',
	[da, da/2], [-Ec2_v - V]*2, '--k'
	)
xlim((-da, da))
yl = ax0.get_ylim()
t_scale = (yl[1] - yl[0])/100.
hold(True)
for j,i in enumerate(p_well_e_):
	plot(
	[-da, da], [-(Ec2_v + V - i[2])]*2, '--m',
	i[0], -(Ec2_v + V - i[2] - i[1]/abs(i[1]).max()*10.*t_scale), 'g')
	text(-da/2 + Ofs[0]/5, -(Ec2_v + V - i[2] + Ofs[1]/50),
		r"$E_{n%d} - V_0$ = %.4f" % (j, i[2]), bbox=props)
text(da , -Ec2_v - V , r"$V_0$", bbox=props_R)
grid(True)
xlabel(r'$Distance, m$')
ylabel(r'$E, eV$')
yl = ax0.get_ylim()

text(-da/2, yl[1] + Ofs[1]/45, r'$-\frac{d}{2}$', fontsize=21, bbox=props)
text( da/2, yl[1] + Ofs[1]/45, r'$\frac{d}{2}$', fontsize=21, bbox=props)

# 2Dh
ax1 = subplot(122)
annotate("",(1,0),xytext=(-al,0), **kwargs) # bottom spine arrow
annotate("",(0,1),xytext=(0,-al), **kwargs) # left spin arrow
plot(
	[-da, -da/2.], [-Ev2_v]*2,'r',
	[da/2., da], [-Ev2_v]*2,'r',
	[-da/2., da/2.], [-Ev2_v + V]*2,'r',
	[-da/2.]*2, [-Ev2_v + V, -Ev2_v],'r',
	[da/2.]*2, [-Ev2_v + V, -Ev2_v],'r',
	[-da/2, -da], [-Ev2_v + V]*2, '--k')
xlim((-da, da))
ax1.yaxis.tick_right()
for j,i in enumerate(p_well_h):
	plot(
	[-da, da], [-(Ev2_v - V + i[2])]*2, '--m',
	i[0], -(Ev2_v - V + i[2] - i[1]/abs(i[1]).max()*10.*t_scale), 'g')
	text(-da/2 + Ofs[0]/5, -(Ev2_v - V + i[2] + Ofs[1]/50),
		r"$V_0-E_{p%d}$ = %.4f" % (j, i[2]), bbox=props)
text(-da , -Ev2_v + V, r"$V_0$", bbox=props_R)

grid(True)
xlabel(r'$Distance, m$')
ylabel(r'$E, eV$')
yl = ax1.get_ylim()

text(-da/2, yl[1] + Ofs[1]/45, r'$-\frac{d}{2}$', fontsize=21, bbox=props)
text( da/2, yl[1] + Ofs[1]/45, r'$\frac{d}{2}$', fontsize=21, bbox=props)

savefig(fname + '_p_well.pdf', dpi=300, facecolor='w', edgecolor='w',
	orientation='landscape', papertype='a4', format='pdf',
	transparent=False, bbox_inches=None, pad_inches=0.1)

savefig(fname + '_p_well.png', dpi=200, facecolor='w', edgecolor='w',
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
	['x_Al_Ga_As' ,  x_Al_Ga_As],
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
	['d_OPZ2' ,  d_OPZ2],
	['V' ,  V],
	['da' ,  da],
	['En1_Fn' ,  En1_Fn],
	['Fp_Ep1' ,  Fp_Ep1],
	['n_2D' , n_2D],
	['p_2D' ,  p_2D]
]

for j,i in enumerate(p_well_e_):
	for name, val in (('E', i[2]), ("K", i[3]), ('a', i[4]), ('A2', i[5])):
		output.append([ name + "_e_" + str(j+1), val])

for j,i in enumerate(p_well_h):
	for name, val in (('E', i[2]), ("K", i[3]), ('a', i[4]), ('A2', i[5])):
		output.append([ name + "_h" + str(j+1), val])


	
import csv
with open(fname + '.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=' ',
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
    for i in output:
    	print '\t', i[0],'\t\t=\t', i[1]
    	spamwriter.writerow(i)


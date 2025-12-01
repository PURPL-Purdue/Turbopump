import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import yaml
import ezdxf

from itertools import chain

from matplotlib.patches import Arc
from bisect import bisect_left

"""
 Implemented from the following technical notes 
 The thrust optimised parabolic nozzle
 http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
 .................................................................

The radius of the nozzle exit: 
Re = √ε * Rt							[Eqn. 2]
and nozzle length 
LN = 0.8 ((√∈−1) * Rt )/ tan(15)		[Eqn. 3]
.................................................................
For the throat entrant section:
x = 1.5 Rt cosθ
y = 1.5 Rt sinθ + 1.5 Rt + Rt			[Eqn. 4]
where: −135 ≤ θ ≤ −90
(The initial angle isn’t defined and is up to the
combustion chamber designer, -135 degrees is typical.)
.................................................................
For the throat exit section:
x = 0.382 Rt cosθ
y = 0.382 Rt sinθ + 0.382 Rt + Rt		[Eqn. 5]
where: −90 ≤ θ ≤ (θn − 90)
.................................................................
The bell is a quadratic Bézier curve, which has equations:
x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]
.................................................................
Selecting equally spaced divisions between 0 and 1 produces 
the points described earlier in the graphical method, 
for example 0.25, 0.5, and 0.75.
.................................................................
Equations 6 are defined by points N, Q, and E (see the graphical method 
earlier for the locations of these points).

Point N is defined by equations 5 setting the angle to (θn – 90).
Nx = 0.382 Rt cos(θn – 90)
Ny = 0.382 Rt sin(θn – 90) + 0.382 Rt + Rt
.................................................................
Coordinate Ex is defined by equation 3, and coordinate Ey is defined by equation 2.
Ex = 0.8*(((√ε−1)-1)*Rt)/(tan(15)) # degrees in rad
Ey = √ε * Rt
.................................................................
Point Q is the intersection of the lines: ⃗⃗⃗⃗⃗⃗
NQ = m1 x + C1 and: ⃗⃗⃗⃗⃗
QE = m2 x + C2 			[Eqn. 7]

where: gradient 
m1 = tan(θn ) , m2 = tan(θe )	[Eqn. 8]

and: intercept 
C1 = Ny − m1 Nx
C2 = Ey − m2 Ex		[Eqn. 9]
.................................................................
The intersection of these two lines (at point Q) is given by:
Qx = (C2 − C1 ) /(m1 − m2 )
Qy = (m1 C2 − m2 C1 ) / (m1 − m2 ) [Eqn. 10]
.................................................................
"""

# sp.heat, area_ratio, throat_radius, length percentage, 
def bell_nozzle(aratio, Rt, l_percent, cratio, alpha, Lc):
	# upto the nozzle designer, usually -135
	entrant_angle  	= -90 - alpha
	ea_radian 		= math.radians(entrant_angle)
	Rc = Rt * np.sqrt(cratio)

	# nozzle length percntage
	if l_percent == 60:		Lnp = 0.6
	elif l_percent == 80:	Lnp = 0.8
	elif l_percent == 90:	Lnp = 0.9	
	else:					Lnp = 0.8
	# find wall angles (theta_n, theta_e) for given aratio (ar)		
	angles = find_wall_angles(aratio, throat_radius, l_percent)
	# wall angles
	nozzle_length = angles[0]; theta_n = angles[1]; theta_e = angles[2];

	data_intervel  	= 100
	# entrant functions
	ea_start 		= ea_radian
	ea_end 			= -math.pi/2	
	angle_list 		= np.linspace(ea_start, ea_end, data_intervel)
	xe = []; ye = [];
	for i in angle_list:
		xe.append( 1.5 * Rt * math.cos(i) )
		ye.append( 1.5 * Rt * math.sin(i) + 2.5 * Rt )

	# linear convergent section functions
	R2max = (Rc - ye[0]) / (1 - np.cos(math.radians(alpha)))
	R2 = 0.300 * R2max
	diag_x0		= xe[0]
	diag_y0 	= ye[0]	
	diag_yf = Rc - R2 * (1 - np.cos(math.radians(alpha)))
	diag_xf = (diag_yf - diag_y0) / (-np.tan(np.radians(alpha))) + diag_x0
	iters 		= np.linspace(diag_x0, diag_xf, data_intervel)
	xed = []; yed = [];
	for i in iters:
		xed.append( i )
		yed.append( diag_y0 + (abs(i - diag_x0) * np.tan(np.radians(alpha))))

	# Combustion chamber arc functions
	cca_start 		= math.radians(90 - alpha)
	cca_end 		= math.pi/2
	angle_list 		= np.linspace(cca_start, cca_end, data_intervel)
	xeca = []; yeca = [];
	for i in angle_list:
		xeca.append( diag_xf + R2 * math.cos(i) - R2 * math.cos(cca_start))
		yeca.append( diag_yf + R2 * math.sin(i) - R2 * math.sin(cca_start))

	# combustion chamber cylinder functions
	iters 		= np.linspace(xeca[-1], -Lc, data_intervel)
	xecc = []; yecc = [];
	for i in iters:
		xecc.append( i )
		yecc.append( yeca[-1] )

	#exit section
	ea_start 		= -math.pi/2
	ea_end 			= theta_n - math.pi/2
	angle_list 		= np.linspace(ea_start, ea_end, data_intervel)	
	xe2 = []; ye2 = [];
	for i in angle_list:
		xe2.append( 0.382 * Rt * math.cos(i) )
		ye2.append( 0.382 * Rt * math.sin(i) + 1.382 * Rt )

	# bell section
	# Nx, Ny-N is defined by [Eqn. 5] setting the angle to (θn – 90)
	Nx = 0.382 * Rt * math.cos(theta_n - math.pi/2)
	Ny = 0.382 * Rt * math.sin(theta_n - math.pi/2) + 1.382 * Rt 
	# Ex - [Eqn. 3], and coordinate Ey - [Eqn. 2]
	Ex = Lnp * ( (math.sqrt(aratio) - 1) * Rt )/ math.tan(math.radians(15) )
	Ey = math.sqrt(aratio) * Rt 
	# gradient m1,m2 - [Eqn. 8]
	m1 = math.tan(theta_n);  m2 = math.tan(theta_e);
	# intercept - [Eqn. 9]
	C1 = Ny - m1*Nx;  C2 = Ey - m2*Ex;
	# intersection of these two lines (at point Q)-[Eqn.10]
	Qx = (C2 - C1)/(m1 - m2)
	Qy = (m1*C2 - m2*C1)/(m1 - m2)	
	
	# Selecting equally spaced divisions between 0 and 1 produces 
	# the points described earlier in the graphical method
	# The bell is a quadratic Bézier curve, which has equations:
	# x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
	# y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]		
	int_list = np.linspace(0, 1, data_intervel)
	xbell = [(xe2[-1])]; ybell = [(ye2[-1])];
	for t in int_list:		
		xbell.append( ((1-t)**2)*Nx + 2*(1-t)*t*Qx + (t**2)*Ex )
		ybell.append( ((1-t)**2)*Ny + 2*(1-t)*t*Qy + (t**2)*Ey )
	
	# create negative values for the other half of nozzle
	nye 	= [ -y for y in ye]
	nye2  	= [ -y for y in ye2]
	nyed    = [ -y for y in yed]
	nyeca   = [ -y for y in yeca]
	nyecc   = [ -y for y in yecc]
	nybell  = [ -y for y in ybell]

	# return
	return angles, (xe, ye, nye, xe2, ye2, nye2, xed, yed, nyed, xeca, yeca, nyeca, xecc, yecc, nyecc, xbell, ybell, nybell)

# find wall angles (theta_n, theta_e) in radians for given aratio (ar)
def find_wall_angles(ar, Rt, l_percent = 80 ):
	# wall-angle empirical data
	aratio 		= [ 4,    5,    10,   20,   30,   40,   50,   100]
	theta_n_60 	= [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 35.0, 40.0]	
	theta_n_80 	= [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
	theta_n_90 	= [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0]
	theta_e_60 	= [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2]
	theta_e_80 	= [14.0, 13.0, 11.0,  9.0,  8.5,  8.0,  7.5,  7.0]
	theta_e_90 	= [11.5, 10.5,  8.0,  7.0,  6.5,  6.0,  6.0,  6.0]	

	# nozzle length
	f1 = ( (math.sqrt(ar) - 1) * Rt )/ math.tan(math.radians(15) )
	
	if l_percent == 60:
		theta_n = theta_n_60; theta_e = theta_e_60;
		Ln = 0.6 * f1
	elif l_percent == 80:
		theta_n = theta_n_80; theta_e = theta_e_80;
		Ln = 0.8 * f1		
	elif l_percent == 90:
		theta_n = theta_n_90; theta_e = theta_e_90;	
		Ln = 0.9 * f1	
	else:
		theta_n = theta_n_80; theta_e = theta_e_80;		
		Ln = 0.8 * f1

	# find the nearest ar index in the aratio list
	x_index, x_val = find_nearest(aratio, ar)
	# if the value at the index is close to input, return it
	if round(aratio[x_index], 1) == round(ar, 1):
		return Ln, math.radians(theta_n[x_index]), math.radians(theta_e[x_index])

	# check where the index lies, and slice accordingly
	if (x_index>2):
		# slice couple of middle values for interpolation
		ar_slice = aratio[x_index-2:x_index+2]		
		tn_slice = theta_n[x_index-2:x_index+2]
		te_slice = theta_e[x_index-2:x_index+2]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)	
	elif( (len(aratio)-x_index) <= 1):
		# slice couple of values initial for interpolation
		ar_slice = aratio[x_index-2:len(x_index)]		
		tn_slice = theta_n[x_index-2:len(x_index)]
		te_slice = theta_e[x_index-2:len(x_index)]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)	
	else:
		# slice couple of end values for interpolation
		ar_slice = aratio[0:x_index+2]		
		tn_slice = theta_n[0:x_index+2]
		te_slice = theta_e[0:x_index+2]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)						

	return Ln, math.radians(tn_val), math.radians(te_val)

# simple linear interpolation
def interpolate(x_list, y_list, x):
	if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
		raise ValueError("x_list must be in strictly ascending order!")
	intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
	slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals]

	if x <= x_list[0]:
		return y_list[0]
	elif x >= x_list[-1]:
		return y_list[-1]
	else:
		i = bisect_left(x_list, x) - 1
		return y_list[i] + slopes[i] * (x - x_list[i])

# find the nearest index in the list for the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]  

# nozzle contour plot
def plot_nozzle(ax, title, Rt, angles, contour):
	# wall angles
	nozzle_length = angles[0]; theta_n = angles[1]; theta_e = angles[2];

	# contour values
	xe = contour[0];   	ye = contour[1];   	nye = contour[2];
	xe2 = contour[3]; 	ye2 = contour[4];  	nye2 = contour[5];
	xed = contour[6];   yed = contour[7];   nyed = contour[8];
	xeca = contour[9]; yeca = contour[10]; nyeca = contour[11];
	xecc = contour[12]; yecc = contour[13]; nyecc = contour[14];
	xbell = contour[15]; ybell = contour[16]; nybell = contour[17];
	
	# plot

	# set correct aspect ratio
	ax.set_aspect('equal')

	# throat enterant
	ax.plot(xe, ye, linewidth=2.5, color='g')
	ax.plot(xe, nye, linewidth=2.5, color='g')
	
	#convergent diagonal
	ax.plot(xed, yed, linewidth=2.5, color='k')
	ax.plot(xed, nyed, linewidth=2.5, color='k')

	#convergent arc
	ax.plot(xeca, yeca, linewidth=2.5, color='r')
	ax.plot(xeca, nyeca, linewidth=2.5, color='r')	

	#combustion chamber cylinder
	ax.plot(xecc, yecc, linewidth=2.5, color='b')
	ax.plot(xecc, nyecc, linewidth=2.5, color='b')	
	
	# throat inlet line
	x1 = xe[0]; y1 = 0;
	x2 = xe[0]; y2 = nye[0];
	dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	# draw arrow, inlet radial line [x1, y1] to [x2, y2] 
	text = ' Ri = '+ str(round(dist,1))
	ax.plot(xe[0], 0, '+' )
	# draw dimension from [x1, y1] to [x2, y2] 
	ax.annotate( "", [x1, y1], [x2, y2] , arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((x1+x2)/2, (y1+y2)/2, text, fontsize=9 )	

	# nozzle inlet length line [0,0] to [xe[0], 0]
	text = ' Li = ' + str( round( abs(xe[0]), 1) ) 
	ax.plot(0,0, '+' )
	# draw dimension from [0,0] to [xe[0], 0]
	ax.annotate( "", [0,0], [xe[0], 0], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xe[0], 0, text, fontsize=9 )	
		
	# find mid point and draw arc radius
	i = int(len(xe)/2)
	xcenter = 0; 	ycenter = 2.5 * Rt;  
	xarch = xe[i];  yarch = ye[i]
	# draw arrow, enterant radial line [xcenter, ycenter] to [xarch, yarch] 
	text = ' 1.5 * Rt = '+ str( round( 1.5 * Rt, 1 ) ) 
	ax.plot(xcenter, ycenter, '+' )
	# draw dimension from [xcenter, ycenter] to [xarch, yarch]
	ax.annotate( "", [xcenter, ycenter], [xarch, yarch], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((xarch+xcenter)/2, (yarch+ycenter)/2, text, fontsize=9 )	
		
	# throat radius line [0,0] to [xe[-1], ye[-1]]
	text = ' Rt = '+ str(Rt)
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	ax.annotate( "", [0,0], [xe[-1], ye[-1]], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xe[-1]/2, ye[-1]/2, text, fontsize=9 )	

	# throat exit
	ax.plot(xe2, ye2, linewidth=2.5, color='r')
	ax.plot(xe2, nye2, linewidth=2.5, color='r')
	# find mid point and draw arc radius
	i = int(len(xe2)/2)
	xcenter2 = 0; 	ycenter2 = 1.382 * Rt;  
	xarch2 = xe2[i];  yarch2 = ye2[i]
	# draw arrow, exit radial line from [xcenter2,ycenter2] to [xarch2, yarch2]
	text = ' 0.382 * Rt = '+ str( round(0.382 * Rt,1) ) 
	ax.plot(xcenter2, ycenter2, '+' )
	# draw dimension from [xcenter2,ycenter2] to [xarch2, yarch2]
	ax.annotate( "", [xcenter2,ycenter2], [xarch2, yarch2], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((xarch2+xcenter2)/2, (yarch2+ycenter2)/2, text, fontsize=9 )

	# draw theta_n, throat inflexion angle
	adj_text = 2
	origin	= [ xe2[-1], nye2[-1]-adj_text ]
	degree_symbol = r'$\theta$n'	
	draw_angle_arc(ax, np.rad2deg(theta_n), origin, degree_symbol )

	# bell section
	ax.plot(xbell, ybell, linewidth=2.5, color='b')
	ax.plot(xbell, nybell, linewidth=2.5, color='b')

	# throat radius line [0,0] to [xe[-1], ye[-1]]
	text = ' Re = ' + str( round( (math.sqrt(aratio) * Rt), 1) ) 
	ax.plot(xbell[-1],0, '+' )
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	ax.annotate( "", [xbell[-1],0], [xbell[-1], ybell[-1]], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xbell[-1], ybell[-1]/2, text, fontsize=9 )	

	# draw theta_e, throat exit angle
	origin	= [ xbell[-1], nybell[-1] ]
	degree_symbol = r'$\theta$e'	
	draw_angle_arc(ax, (np.rad2deg(theta_e)), origin, degree_symbol )

	# nozzle length line [0,0] to [xe[-1], ye[-1]]
	text = ' Ln = ' + str( round( nozzle_length, 1) ) 
	ax.plot(0,0, '+' )
	# draw dimension from [0,0] to [xbell[-1], 0]
	ax.annotate( "", [0,0], [xbell[-1], 0], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xbell[-1]/2, 0, text, fontsize=9 )	
				
	# axis
	ax.axhline(color='black', lw=0.5, linestyle="dashed")
	ax.axvline(color='black', lw=0.5, linestyle="dashed")		
	
	# grids
	ax.grid()
	ax.minorticks_on()
	ax.grid(which='major', linestyle='-', linewidth='0.5') # , color='red'
	ax.grid(which='minor', linestyle=':', linewidth='0.5') # , color='black'	
	
	# show
	plt.title(title, fontsize=9)
	return

# nozzle contour plot
def plot_nozzle_inches(contour, dia_t, dia_c, dia_e, len_c):

	# contour values
	xe = np.divide(contour[0], 25.4); ye = np.divide(contour[1], 25.4); nye = np.divide(contour[2], 25.4);
	xe2 = np.divide(contour[3], 25.4); ye2 = np.divide(contour[4], 25.4); nye2 = np.divide(contour[5], 25.4);
	xed = np.divide(contour[6], 25.4); yed = np.divide(contour[7], 25.4); nyed = np.divide(contour[8], 25.4);
	xeca = np.divide(contour[9], 25.4); yeca = np.divide(contour[10], 25.4); nyeca = np.divide(contour[11], 25.4);
	xecc = np.divide(contour[12], 25.4); yecc = np.divide(contour[13], 25.4); nyecc = np.divide(contour[14], 25.4);
	xbell = np.divide(contour[15], 25.4); ybell = np.divide(contour[16], 25.4); nybell = np.divide(contour[17], 25.4);
	
	# plot

	fig2 = plt.figure(2, figsize=(12,9))

	# throat enterant
	plt.plot(xe, ye, linewidth=2.5, color='k')
	plt.plot(xe, nye, linewidth=2.5, color='k')
	
	#convergent diagonal
	plt.plot(xed, yed, linewidth=2.5, color='k')
	plt.plot(xed, nyed, linewidth=2.5, color='k')

	#convergent arc
	plt.plot(xeca, yeca, linewidth=2.5, color='k')
	plt.plot(xeca, nyeca, linewidth=2.5, color='k')	

	#combustion chamber cylinder
	plt.plot(xecc, yecc, linewidth=2.5, color='k')
	plt.plot(xecc, nyecc, linewidth=2.5, color='k')	

	# throat exit
	plt.plot(xe2, ye2, linewidth=2.5, color='k')
	plt.plot(xe2, nye2, linewidth=2.5, color='k')

	# bell section
	plt.plot(xbell, ybell, linewidth=2.5, color='k')
	plt.plot(xbell, nybell, linewidth=2.5, color='k')

	# throat diameter line
	text = str(round(dia_t,2)) + 'in'
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	plt.annotate( "", [xe[-1], 0.95 * nye[-1]], [xe[-1], 0.95 * ye[-1]], arrowprops=dict(lw=1, arrowstyle='|-|') )
	plt.text(0.1,0.1, text, fontsize=12 )	

	# chamber diameter line
	text = str(round(dia_c,2)) + 'in'
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	plt.annotate( "", [xecc[-1], 0.95 * nyecc[-1]], [xecc[-1], 0.95 * yecc[-1]], arrowprops=dict(lw=1, arrowstyle='|-|') )
	plt.text(xecc[-1] + 0.1,0.1, text, fontsize=12 )	

	# exit diameter line
	text = str(round(dia_e,2)) + 'in'
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	plt.annotate( "", [xbell[-1], 0.95 * nybell[-1]], [xbell[-1], 0.95 * ybell[-1]], arrowprops=dict(lw=1, arrowstyle='|-|') )
	plt.text(4.5,0.1, text, fontsize=12 )	

	# chamber length line
	text = str(round(len_c,2)) + 'in'
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	plt.annotate( "", [xecc[-1], 1.2 * ybell[-1]], [0, 1.2 * ybell[-1]], arrowprops=dict(lw=1, arrowstyle='|-|') )
	plt.text((-len_c / 2), 1.25* ybell[-1], text, fontsize=12 )

	# divergent section length line
	text = str(round(xbell[-1],2)) + 'in'
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	plt.annotate( "", [0, 1.2 * ybell[-1]], [xbell[-1], 1.2 * ybell[-1]], arrowprops=dict(lw=1, arrowstyle='|-|') )
	plt.text((xbell[-1] / 2), 1.25* ybell[-1], text, fontsize=12 )

	# axis
	plt.axhline(color='black', lw=0.5, linestyle="dashed")
	plt.axvline(color='black', lw=0.5, linestyle="dashed")		
	
	# grids
	plt.grid()
	plt.minorticks_on()
	plt.grid(which='major', linestyle='-', linewidth='0.5') # , color='red'
	plt.grid(which='minor', linestyle=':', linewidth='0.5') # , color='black'	
	
	# show
	plt.xlabel('Inches', fontsize=9)
	plt.ylabel('Inches', fontsize=9)
	plt.axis('equal')
	fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
	
	plt.savefig(r"TCA\Countour Exports\nozzle_contour_inches_plot.png")

	plt.show()
	return

# theta_n in rad,  origin =[startx, starty], degree symbol
def draw_angle_arc(ax, theta_n, origin, degree_symbol=r'$\theta$'):
	length = 50
	# start point
	startx = origin[0]; starty = origin[1];
	# find the end point
	endx = startx + np.cos(-theta_n) * length * 0.5
	endy = starty + np.sin(-theta_n) * length * 0.5
	# draw the angled line
	ax.plot([startx,endx], [starty,endy], linewidth=0.5, color='k')
	# horizontal line
	# ax.hlines(y=starty, xmin=startx, xmax=length, linewidth=0.5, color='k')
	# angle
	arc_obj = Arc([startx, starty], 1, 1, angle=0, theta1=0, theta2=math.degrees(theta_n), color='k' )
	ax.add_patch(arc_obj)
	ax.text(startx+0.5, starty+0.5, degree_symbol + ' = ' + str(round(theta_n,1)) + u"\u00b0")	
	return

# ring of radius r, height h, base point a
def ring(r, h, a=0, n_theta=30, n_height=10):
    theta = np.linspace(0, 2*np.pi, n_theta)
    v = np.linspace(a, a+h, n_height )
    theta, v = np.meshgrid(theta, v)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    z = v
    return x, y, z

# Set 3D plot axes to equal scale. 
# Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
def set_axes_equal_3d(ax: plt.Axes):
    """	
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)
    return

# set axis limits
def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
    return

# 3d plot
def plot3D(ax, contour):
	# unpack the contour values
	xe = contour[0];   	ye = contour[1];   	nye = contour[2];
	xe2 = contour[3]; 	ye2 = contour[4];  	nye2 = contour[5];
	xed = contour[6];   yed = contour[7];   nyed = contour[8];
	xeca = contour[9]; yeca = contour[10]; nyeca = contour[11];
	xecc = contour[12]; yecc = contour[13]; nyecc = contour[14];
	xbell = contour[15]; ybell = contour[16]; nybell = contour[17];

	# collect and append array values
	x = []; y = [];
	x = np.append(x, xecc);  y = np.append(y, yecc)
	x = np.append(x, xeca);  y = np.append(y, yeca)
	x = np.append(x, xed);  y = np.append(y, yed)
	x = np.append(x, xe);  y = np.append(y, ye)	
	x = np.append(x, xe2);  y = np.append(y, ye2)	
	x = np.append(x, xbell);  y = np.append(y, ybell)	
	
	# ring thickness
	thick = 5 * (x[1] - x[0]) # 0.01
	
	# draw multiple rings to create 3d structure
	for i in range(len(y)):
		X, Y, Z = ring(y[i], thick, x[i])
		ax.plot_surface(X, Y, Z, color='g')	
		
	# set correct aspect ratio
	ax.set_box_aspect([1,1,1])
	set_axes_equal_3d(ax)
	# set view
	ax.view_init(-170, -15)
	return

def plot(title, throat_radius, angles, contour):
	# Plot 3d view
	fig1 = plt.figure(1, figsize=(12,9))
	# plot some 2d information
	ax1 = fig1.add_subplot(121)
	plot_nozzle(ax1, title, throat_radius, angles, contour)
	# plot 3d view
	ax2 = fig1.add_subplot(122, projection='3d')
	plot3D(ax2, contour)	
	# show
	fig1.tight_layout(rect=[0, 0.03, 1, 0.95])

	plt.savefig(r"TCA\Countour Exports\nozzle_contour_plot.png")

	plt.show()
	return

######################### EXPORT TO CSV FILE #########################
def export_nozzle_csv(contour, filename=r"TCA\Countour Exports\nozzle_contour.csv"):
    """
    Write ONE CSV that unites all contour points.
    Expected contour structure from bell_nozzle(...):
      [xe, ye, nye,  xe2, ye2, nye2,  xbell, ybell, nybell]
    Only x,y arrays are written.
    CSV columns: segment,x,y,index
    """
    if not isinstance(contour, (list, tuple)) or len(contour) < 9:
        raise ValueError("Unexpected contour structure. Need at least 9 elements as returned by bell_nozzle().")

    xe,   ye   = contour[0], contour[1]
    xe2,  ye2  = contour[3], contour[4]
    xed,  yed  = contour[6], contour[7]
    xeca, yeca = contour[9], contour[10]
    xecc, yecc = contour[12], contour[13]
    xbell,ybell= contour[15], contour[16]

    segments = [
        ("chamber_wall", xecc, yecc),
		("convergent_arc", xeca, yeca),
		("convergent_diagonal", xed, yed),
		("throat_arc", xe,    ye),
        ("inlet_arc",  xe2,   ye2),
        ("bell",       xbell, ybell),
    ]

    with open(filename, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["segment","x","y","index"])
        for name, xs, ys in segments:
            n = min(len(xs), len(ys))
            for i in range(n):
                w.writerow([name, float(xs[i]), float(ys[i]), i])

    return filename


def export_nozzle_dxf(contour):
	
    xe,   ye   = np.divide(contour[0], 1000), np.divide(contour[1], 1000)
    xe2,  ye2  = np.divide(contour[3], 1000), np.divide(contour[4], 1000)
    xed,  yed  = np.divide(contour[6], 1000), np.divide(contour[7], 1000)
    xeca, yeca = np.divide(contour[9], 1000), np.divide(contour[10], 1000)
    xecc, yecc = np.divide(contour[12], 1000), np.divide(contour[13], 1000)
    xbell,ybell= np.divide(contour[15], 1000), np.divide(contour[16], 1000)
    
    xed = xed[::-1]
    yed = yed[::-1]
    xeca = xeca[::-1]
    yeca = yeca[::-1]    
    xecc = xecc[::-1]
    yecc = yecc[::-1]

    xs = list(chain(xecc, xeca, xed, xe, xe2, xbell))
    ys = list(chain(yecc, yeca, yed, ye, ye2, ybell))

    points = list(zip(xs, ys)) 

    doc = ezdxf.new("R2010")
    msp = doc.modelspace()

    msp.add_lwpolyline(points, close=False)
    doc.saveas(r"TCA\Countour Exports\nozzle_contour.dxf")
    return()

# __main method__
if __name__=="__main__":

	with open(r'/Users/dl/Documents/GitHub/Turbopump/TCA/TCA_params.yaml') as file: ##CHANGED FOR DANIEL
		tca_params = yaml.safe_load(file)

	l_percent = tca_params['bell_nozzle_l_percent']	# nozzle length percntage

	# typical upper stage values
	aratio = tca_params['tca_expansion_ratio'] # Ae / At	
	cratio = tca_params['tca_contraction_ratio']  #Ac / At
	cangle = tca_params['tca_convergent_half_angle']   #contraction angle (alpha)
	clength = tca_params['tca_chamber_length']  * 25.4   #chamber length (mm)
	throat_radius = tca_params['tca_throat_diameter'] * 25.4 / 2		#throat radius (mm)

	###CHANGE TO EXPORT WHOLE CONTOUR
	# rao_bell_nozzle_contour
	angles, contour = bell_nozzle(aratio, throat_radius, l_percent, cratio, cangle, clength)
	title = 'Bell Nozzle \n [Area Ratio = ' + str(round(aratio,1)) + ', Throat Radius = ' + str(round(throat_radius,2)) + 'mm]' 
	export_nozzle_csv(contour, filename=r"TCA\Countour Exports\nozzle_contour.csv")
	export_nozzle_dxf(contour)
	plot(title, throat_radius, angles, contour)

	Dt_in = tca_params['tca_throat_diameter']
	Dc_in = tca_params['tca_chamber_diameter']
	De_in = tca_params['tca_exit_diameter']
	Lc_in = tca_params['tca_chamber_length']
	plot_nozzle_inches(contour, Dt_in, Dc_in, De_in, Lc_in)

	
	
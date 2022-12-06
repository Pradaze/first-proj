import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = input("enter date in yy-mm-dd format: ")
G = float(input('enter value of Gravitational constant(conventional value=6.674e-11): ')) * 2.2290714e-24
mass_sun = float(input('Enter mass of the sun(conventional value=1.989e30): '))
dur = int(input('enter duration of simulation(in years): '))
sim_duration = dur * 365
list = ["Uranus's axis of rotation is so tilted, that it appears to be rolling\n on its side while orbiting."
    ,'Jupiter has over 79 moons discovered so far. 4 of them potentially\nhabitable for life.','The tallest mountain in the solar system is Olympus Mons, located\n on Mars which is about 3 times higher than Mt.Everest.',
        'Saturn is gaseous planet whose density is less than water. So in theory, \nSaturn would float in a large enough ocean.','All bodies in the solar system rotate in the same plane in same direction.\nUranus and Venus are the only planets which orbit in the opposite direction.','The Kuiper belt is a disk of icy objects beyond the orbit of Neptune. It is\nmuch larger than the astroid belt even containing dwarf planets like\nPluto and Makemake.',"Neptune has the wildest winds in the solar system, with wind speeds\ncrossing 2000km/h. For comparison, earth's fastest winds only reach 400km/h."]


class Obj:
    def __init__(self, name, rad, color, r, v):
        self.name = name
        self.r = np.array(r, dtype=np.float)
        self.v = np.array(v, dtype=np.float)
        self.xs = []
        self.ys = []
        self.plot = ax.scatter(r[0], r[1], color=color, s=rad**2, zorder=10)
        self.line, = ax.plot([], [], color=color, linewidth=1.4)

class slr_sys:
    def __init__(self, sun):
        self.sun = sun
        self.planets = []
        self.timestamp = ax.text(.05, .97, 'Date: ', color='w', transform=ax.transAxes, fontsize='small')
    def add_planet(self, planet):
        self.planets.append(planet)
    def evolve(self):
        dt = 1.0
        self.time += dt
        digi = str(Time(self.time, format='jd', out_subfmt='str').iso)[5]+ str(Time(self.time, format='jd', out_subfmt='str').iso)[6]
        ft = int(digi)//2
        self.fact = ax.text(-2.15,1.6,'\nFACT:'+list[ft],color='green')
        plots = []
        lines = []
        for p in self.planets:
            p.r += p.v * dt
            acc = -G * mass_sun * p.r / np.sum(p.r**2)**(3./2)
            p.v += acc * dt
            p.xs.append(p.r[0])
            p.ys.append(p.r[1])
            p.plot.set_offsets(p.r[:2])
            p.line.set_xdata(p.xs)
            p.line.set_ydata(p.ys)
            plots.append(p.plot)
            lines.append(p.line)
        self.timestamp.set_text('Date: ' + str(Time(self.time, format='jd', out_subfmt='str').iso)[:9] )
        self.fact.set_text('\nFACT:'+list[ft])
        return plots + lines + [self.timestamp] + [self.fact]

plt.style.use('dark_background')
fig = plt.figure(figsize=[8, 8])
ax = plt.axes([0., 0., 1., 1.], xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
ax.set_aspect('equal')
ax.axis('on')
ss = slr_sys(Obj("Sun", 28/1.2, 'red', [0, 0, 0], [0, 0, 0]))
ss.time = Time(sim_start_date).jd
colors = ['#646464', '#FFE873', '#02ebf2', '#cb5808']
sizes = [0.38/1.2, 0.95/1.2, 1./1.2, 0.53/1.2]
names = ['Mercury', 'Venus', 'Earth', 'Mars']
texty = [.47, .73, 1, 1.5]
for i, nasaid in enumerate([1, 2, 3, 4]):
    obj = Horizons(id=nasaid, location="@sun", epochs=ss.time, id_type='id').vectors()
    ss.add_planet(Obj(nasaid, 20 * sizes[i], colors[i],
                         [np.double(obj[xi]) for xi in ['x', 'y', 'z']],
                         [np.double(obj[vxi]) for vxi in ['vx', 'vy', 'vz']]))
    ax.text(0, - (texty[i] + 0.1), names[i], color=colors[i], zorder=1000, ha='center', fontsize='medium')
def animate(i):
    return ss.evolve()
ani = anm.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=True, interval=20,)
plt.show()









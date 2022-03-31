import matplotlib.animation as anim
import matplotlib.pyplot as plt

from scipy.constants import au, G
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from math import (sqrt, sin, cos, atan2)
from datetime import datetime as dt, timedelta

SCALE = 250 / (149.6e6 * 1000)

class SolarSystem:

    def __init__(self, bodies, date, steps, timestep):
        
        self.bodies = bodies
        self.date = date
        self.steps = steps
        self.timestep = timestep  

    def update_simulation(self):    
        
        self.update_time()

        for body in self.bodies:
            if not hasattr(body, 'plot'):
                body.get_initial_condition(self.date)
                body.plot = ax.scatter(
                    body.position['x'] * SCALE, 
                    body.position['y'] * SCALE, 
                    color=body.color, 
                    s=body.radius**2
                )       
            else:         
                body.calc_body_position(self.bodies, self.timestep)
                body.plot.set_offsets([body.position['x'] * SCALE, body.position['y'] * SCALE]) 

    def update_time(self):
        if not hasattr(self, 'timestamp'):
            self.time = dt.strptime(self.date, '%Y-%m-%d')   
            self.timestamp = ax.text(
                .05, .95, 
                f'Date: {self.time.strftime("%Y-%m-%d")}', 
                color='w', 
                transform=ax.transAxes, 
                fontsize='x-large'
            )
        else:
            self.time += timedelta(seconds=self.timestep)
            self.timestamp.set_text(f'Date: {self.time.strftime("%Y-%m-%d")}')   

class Body:

    def __init__(self, id, name, mass, radius, color):
        
        self.id = id
        self.name = name
        self.mass = mass
        self.radius = radius
        self.color = color

    def get_initial_condition(self, date):
        data = Horizons(id=self.id, epochs=Time(date).jd).vectors() 
        self.position = {p: data[p][0] * au for p in ['x', 'y']}
        self.velocity = {v: data[v][0] * au / (24 * 60 * 60) for v in ['vx', 'vy']}

    def calc_body_position(self, bodies, timestep):
        
        total_fx = total_fy = 0.0
        for other_body in bodies:
            
            if self is other_body: 
                continue
            
            dx = other_body.position['x'] - self.position['x']
            dy = other_body.position['y'] - self.position['y']

            d = sqrt((dx ** 2) + (dy ** 2))
            
            f = G * self.mass * other_body.mass / (d ** 2)

            theta = atan2(dy, dx)
            fx = cos(theta) * f
            fy = sin(theta) * f                

            total_fx += fx
            total_fy += fy

        self.velocity['vx'] += total_fx / self.mass * timestep
        self.velocity['vy'] += total_fy / self.mass * timestep
        
        self.position['x'] += self.velocity['vx'] * timestep
        self.position['y'] += self.velocity['vy'] * timestep

if __name__ == '__main__':

    # create plot
    plt.style.use('dark_background')
    fig = plt.figure(figsize=[10, 10])
    ax = plt.axes([0., 0., 1., 1.], xlim=(-1000, 1000), ylim=(-1000, 1000))
    ax.set_aspect('equal')
    ax.axis('off')
    
    # start date
    start_date = '2021-01-01'

    # steps in days
    steps = 365

    # day in seconds
    timestep = 24 * 60 * 60 
    
    sun = Body(10, 'Sun', 1.989e30, 20, 'yellow')
    mercury = Body(1, 'Mercury', 0.330e24, 5, 'grey')
    venus = Body(2, 'Venus', 4.87e24, 10, 'orange')
    earth = Body(3, 'Earth', 5.97e24, 10, 'blue')
    mars = Body(4, 'Mars', 0.642e24, 5, 'red')
    jupiter = Body(5, 'Jupiter', 1898e24, 15, 'brown')
    saturn = Body(6, 'Saturn', 568e24, 13, 'brown')
    uranus = Body(7, 'Uranus', 86.8e24, 12, 'blue')
    neptune = Body(8, 'Neptune', 102e24, 12, 'blue')

    bodies = [sun, mercury, venus, earth, mars]

    solar_system = SolarSystem(bodies, start_date, steps, timestep)
 
    def animate(i):
        return solar_system.update_simulation()

    a = anim.FuncAnimation(fig, animate, interval=20, frames=steps, repeat=False)
    # writer_gif = anim.PillowWriter(fps=30)
    # a.save('solar_system.gif', writer=writer_gif)
    plt.show()
    

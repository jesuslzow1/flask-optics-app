from flask import render_template, request, url_for, redirect, Response
from app import app
import io
from matplotlib.backends.backend_svg import FigureCanvasSVG
from matplotlib.figure import Figure
import matplotlib.lines as mlines
import numpy as np

@app.route('/')
@app.route('/home', methods = ['GET', 'POST'])
def home():
    n1 = float(request.args.get("rng_n1", 1.0))
    n2 = float(request.args.get("rng_n2", 1.0))
    Ang_i = int(request.args.get("rng_Ang_i", 60))

    return render_template('index.html', activar_home = True, n1 = n1, n2 = n2, Ang_i = Ang_i)

@app.route("/matplot-as-image-<float:n1>-<float:n2>-<int:Ang_i>.svg")
def plot_svg(n1=1, n2=1, Ang_i=50):
    """ renders the plot on the fly.
    """
    svg_figure = SvgCanvas(n1, n2, Ang_i)

    return Response(svg_figure.drawPlot(), mimetype="image/svg+xml")


class SvgCanvas:
    def __init__(self, n1 = 1.0, n2 = 1.52, Ang_i = 60):
        self.fig = Figure()
        self.n1 = n1
        self.n2 = n2
        self.Ang_i = Ang_i
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.Ang_refr = 0
        self.drawPlot()
    def drawPlot(self):
        self.ax.clear()
        self.Ang_ref = self.Ang_i
        #ley de refraccion
        self.Ang_refr=180.0/np.pi*np.arcsin(self.n1*np.sin(np.pi/180.0*self.Ang_i)/self.n2)
        x= np.linspace(-3.0,3.0,1024)
        y= np.linspace(-3.0,3.0,1024)
        y= np.sin(2*np.pi*x/(self.n1/1000)+self.n2)
        bdy= mlines.Line2D([-3,3],[0,0],color='k')
        self.ax.add_line(bdy)
        #dibuja la normal 
        nml=mlines.Line2D([0,0],[-3,3],ls='dashed',color='k',label='Normal')
        self.ax.add_line(nml)
        #Draw the incident ray
        lIncRay= lRefrRay = lReflRay = 2.8
        xIncRay=-lIncRay*np.sin(np.pi/180.0*self.Ang_i)
        yIncRay=lIncRay*np.cos(np.pi/180.0*self.Ang_i)
        incRay= mlines.Line2D([xIncRay,0],[yIncRay,0],color='y',label='Rayo de incidencia')
        self.ax.add_line(incRay)
        y = (yIncRay / xIncRay) * (0.005 - xIncRay / 2) + yIncRay / 2
        self.ax.arrow(xIncRay/2,yIncRay/2,0.005,y,head_width=0.1,color='y',shape='full')
		
		#draw the reflected ray
        xReflRay= -xIncRay
        yReflRay= yIncRay
        reflRay= mlines.Line2D([xReflRay,0],[yReflRay,0],color='g',label='Rayo reflejado')
        self.ax.add_line(reflRay)
        y = (yReflRay / xReflRay) * (0.005 - xReflRay / 2) + yReflRay / 2
        self.ax.arrow(xReflRay/2,yReflRay/2,0.005,y,head_width=0.1,color='g',shape='full')


		#rayo refractado
        if not isinstance (self.Ang_refr, complex):
            xRefrRay= lRefrRay*np.sin(np.pi/180.0*self.Ang_refr)
            yRefrRay= -lRefrRay*np.cos(np.pi/180.0*self.Ang_refr)
            refrRay= mlines.Line2D([0,xRefrRay],[0,yRefrRay],color='r',label='Rayo refractado')
            self.ax.add_line(refrRay)
            y = (yRefrRay / xRefrRay) * (0.005 - xRefrRay / 2) + yRefrRay / 2
            self.ax.arrow(xRefrRay/2,yRefrRay/2,0.005,y,length_includes_head=True,head_width=0.1,color='r',shape='full')
        self.ax.legend()
        self.ax.set_ylim(-3,3)
        self.ax.set_xlim(-3,3)
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax.text(-2.5,2.5,'Medio 1')
        self.ax.text(2.0,-2.5,'Medio 2')
        self.ax.text(0.0,-3.5,'Ángulo de refracción= '+str(self.Ang_refr))
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        return output.getvalue()
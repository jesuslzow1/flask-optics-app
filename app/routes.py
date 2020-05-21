from flask import render_template, request, url_for, redirect, Response
from app import app
import io
from matplotlib.backends.backend_svg import FigureCanvasSVG
from matplotlib.figure import Figure
import matplotlib.lines as mlines
from matplotlib import patches as mpatches
import numpy as np

@app.route('/')
@app.route('/home', methods = ['GET', 'POST'])
def home():
    activar_panel =[ 
        True,
        False,
        False]
    
    n1 = float(request.args.get("rng_n1", 1.0))
    n2 = float(request.args.get("rng_n2", 1.0))
    Ang_i = int(request.args.get("rng_Ang_i", 60))

    return render_template('index.html', activar_panel = activar_panel, activar_home = True, n1 = n1, n2 = n2, Ang_i = Ang_i)

@app.route("/matplot-as-image-<float:n1>-<float:n2>-<int:Ang_i>.svg")
def plot_svg(n1=1, n2=1, Ang_i=50):
    """ renders the plot on the fly.
    """
    svg_figure = SvgCanvas(n1, n2, Ang_i)

    return Response(svg_figure.drawPlot(), mimetype="image/svg+xml")

x = 10
@app.route('/lentes', methods = ['GET', 'POST'])
def lentes():
    activar_panel = [
         False,
         True,
         False
    ]
    n1 = float(request.args.get("rng_n1_lentes", 1.0))
    n2 = float(request.args.get("rng_n2_lentes", 2.22))
    radio_curvatura = float(request.args.get("rng_rc_lentes", 45.0))
    distancia_objeto = float(request.args.get("rng_do_lentes", 90.0))
    altura_objeto = float(request.args.get("rng_ao_lentes", 45.0))
    return render_template('lentes.html', activar_panel = activar_panel, activar_home = False, n1 = n1, n2 = n2, radio_curvatura = radio_curvatura, distancia_objeto = distancia_objeto, altura_objeto = altura_objeto)

@app.route("/lentes-as-image-<float:n1>-<float:n2>-<float:radio_curvatura>-<float:distancia_objeto>-<float:altura_objeto>.svg")
def plot_svg_lentes(n1=1, n2=1, radio_curvatura = 25.0, distancia_objeto = 50.0, altura_objeto = 25.0):
    """ renders the plot on the fly.
    """
    
    svg_figure = MplCanvas(nOne = n1, nTwo = n2, rad = radio_curvatura, objDist=distancia_objeto, objHt=altura_objeto)
    #x = 
    return Response(svg_figure.drawPlot(n1, n2, radio_curvatura, distancia_objeto, altura_objeto), mimetype="image/svg+xml")


@app.route('/about')
def about():
    activar_panel = [
         False,
         False,
         True]
    
    return render_template('about.html', activar_home = False, activar_panel = activar_panel)

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

class MplCanvas:
    def __init__(self, parent=None, width=10, height=10, nOne=1.0, nTwo=1.5, rad=30, objDist=50, objHt=5):
        self.fig = Figure()
        self.ax = self.fig.add_subplot(111, aspect='equal')
        self.imgDist=0
        self.imgHt=0
        self.fOne=0
        self.fTwo=0
        self.mag=0
        self.drawPlot(nOne, nTwo, rad, objDist, objHt)
        self.x_i = 0

    def yCoordLine(self, xOne, xTwo, yOne, yTwo, x):
        slope = (yTwo - yOne)/(xTwo - xOne)
        return slope*(x - xOne) + yOne

    def drawPlot(self, nOne, nTwo, rad, objDist, objHt):
        self.ax.clear()
        power = (nTwo - nOne)/(rad/100.0)
        if not power == 0:
            self.fOne = -nOne/power*100.0
            self.fTwo = nTwo/power*100
        else:
            self.fOne = -1e10
            self.fTwo = 1e10
        objVerg = -nOne/(objDist/100)
        imgVerg = objVerg + power
        if not imgVerg == 0:
            self.imgDist = nTwo/(imgVerg) * 100 # valor en cm
            self.mag = objVerg/imgVerg
            self.imgHt = self.mag * objHt
        else:
            self.imgDist = self.mag = self.imgHt = 1e10
        cm2inch = 1.0/2.54
        rad1 =cm2inch*rad
        l = -cm2inch*objDist
        lp = cm2inch*self.imgDist
        h = cm2inch*objHt
        hp = cm2inch*self.imgHt
        f = cm2inch*self.fOne
        fp =cm2inch*self.fTwo
        x = cm2inch*np.linspace(-200, 400, 1024)
        y = cm2inch*np.linspace(-200, 200, 1024)
		#eje optico
        optAx = mlines.Line2D([x[0],x[-1]], [0, 0], color='k',label='eje optico')
        self.ax.add_line(optAx)
        sfc = mpatches.Arc([rad1, 0], 2*rad1, 2*rad1, angle=0, theta1=85, theta2=265, color='blue')
        self.ax.add_patch(sfc)
        #draw object
        self.ax.arrow(l, 0, 0, h, length_includes_head=True, head_width=0.5, color='k', shape='full')
        self.ax.text(l-7, 0.35*h, 'h')
        #1rayo incidente paralelo al eje optico
        sag = rad1 - (rad1**2-h**2)**0.5
        iR1 = mlines.Line2D([l, sag], [h, h], color='yellow',label='Rayo 1')
        self.ax.add_line(iR1)
        self.ax.arrow(l/2, h, 2, 0, length_includes_head=True, head_width=1, color='k', shape='full',label='Rayo 1')
        #rayo que pasa por el punto focal hasta la imagen
        endPtX = 1.4*lp
        endPtY = self.yCoordLine(sag, fp, h, 0, endPtX)
        rR1 = mlines.Line2D([sag, endPtX], [h, endPtY], color='yellow')
        self.ax.add_line(rR1)
        arBaseX = (endPtX-sag)/2
        arBaseY = self.yCoordLine(sag, fp, h, 0, arBaseX)
        dX = 1
        dY = h/(sag-fp)
        self.ax.arrow(arBaseX, arBaseY, dX, dY, length_includes_head=True, head_width=1, color='k',shape='full')
        #2 rayo incidente que pasa a traves del foco y luego se va paralelo al eje
        iR2Slope = (h-0)/(l-f) 
        iR2Intercept = iR2Slope*(-l) + h
        if abs(iR2Intercept) < rad1:
            iR2EndPtY = iR2Intercept
            rR2 = mlines.Line2D([0, 1.4*lp], [iR2Intercept , iR2Intercept], color='r',label='Rayo 2')
            self.ax.add_line(rR2)
        else:
            iR2EndPtY = iR2Slope*rad1 + iR2Intercept
        iR2 = mlines.Line2D([l, 0], [h, iR2EndPtY], color='r')
        self.ax.add_line(iR2)
        #tercer rayo que pasa por el centro de curvatura del lente
        iR3Slope = (h-0)/(l-rad1)
        iR3Intercept = iR3Slope*(-l) + h
        iR3EndPtX = 1.4*lp
        iR3EndPtY = iR3Slope*iR3EndPtX + iR3Intercept
        iR3 = mlines.Line2D([l, iR3EndPtX], [h, iR3EndPtY], color='green',label='Rayo 3')
        self.ax.add_line(iR3)
        #dibuja la imagen
        self.ax.arrow(lp, 0, 0, hp, length_includes_head=True, head_width=1, color='k',shape='full')
        self.ax.text(1.025*lp, 0.5*hp, 'yi') 
        irV = mlines.Line2D([0, 0], [0, 0.95*y[0]], color='k', ls='dashed')
        self.ax.add_line(irV)
        self.ax.arrow(0, -rad1-5, f, 0, length_includes_head=True, head_width=1, color='k',shape='full')
        self.ax.text(0.5*f, -rad1-10, 'fo')
        self.ax.arrow(0, -rad1-20, l, 0, length_includes_head=True, head_width=1, color='k',shape='full')
        self.ax.text(0.5*l, -rad1-25, 'so')
        self.ax.arrow(0, -rad1+2.5, rad1, 0, length_includes_head= True, head_width=1, color='k',shape='full')
        self.ax.text(0.5*rad1, -rad1+2.5, 'r')
        self.ax.arrow(0, -rad1-12.5, fp, 0, length_includes_head=True, head_width=1, color='k',shape='full')
        self.ax.text(0.5*fp, -rad1-12, 'fi')
        self.ax.arrow(0, -rad1-26.5, lp, 0, length_includes_head=True, head_width=1, color='k',shape='full')
        self.ax.text(0.5*lp, -rad1-26, 'si')
        self.ax.legend()
        self.ax.set_xlim(x[0],x[-1])
        self.ax.set_ylim(y[0],y[-1])
        self.ax.set_xticklabels([])
        self.ax.set_yticklabels([])
        self.ax.text(0.5*x[0], 0.2*y[-1], str('n1 = %.2f' %(nOne)))
        self.ax.text(0.5*x[-1],0.2*y[-1],str('n2 = %.2f' %(nTwo)))
        self.ax.text(-80,-90,'si='+str(round(self.imgDist,2))+'cm')
        self.ax.text(-30,-90,'yi='+str(round(self.imgHt,2))+'cm')
        self.ax.text(20,-90,'fo='+str(round(self.fOne,2))+'cm')
        self.ax.text(-80,-100,'fi='+str(round(self.fTwo,2))+'cm')
        self.ax.text(20,-100,'Mt='+str(round(self.mag,2)))
        #self.draw()
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        return output.getvalue()
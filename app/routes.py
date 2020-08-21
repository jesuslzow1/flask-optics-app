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
        False,
        False,
        False,
        False]
    
    n1 = float(request.args.get("rng_n1", 1.0))
    n2 = float(request.args.get("rng_n2", 1.52))
    Ang_i = float(request.args.get("rng_Ang_i", 60.0))

    return render_template('index.html', activar_panel = activar_panel, activar_home = True, n1 = n1, n2 = n2, Ang_i = Ang_i)

@app.route("/matplot-as-image-<float:n1>-<float:n2>-<float:Ang_i>.svg")
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
         False,
         True,
         False,
         False,
         False
    ]
    n1 = float(request.args.get("rng_n1_lentes", 1.0))
    n2 = float(request.args.get("rng_n2_lentes", 1.5))
    radio_curvatura = float(request.args.get("rng_rc_lentes", 25.0))
    distancia_objeto = float(request.args.get("rng_do_lentes", 68.0))
    altura_objeto = float(request.args.get("rng_ao_lentes", 9.0))
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
         False,
         False,
         True]
    
    return render_template('about.html', activar_home = False, activar_panel = activar_panel)

@app.route('/RLC', methods = ['GET', 'POST'])
def RLC():
    activar_panel = [
         False,
         False,
         False,
         True,
         False,
         False]
    vi = float(request.args.get("rng_vi", 10.0))
    Ii = float(request.args.get("rng_Ii", -1.0))
    R = float(request.args.get("rng_R", 100.0))
    C = float(request.args.get("rng_C",0.002))
    L= float(request.args.get("rng_L",3.0))
    circuito = int(request.args.get("rng_circuito", 1))
    if circuito == 1:
        circuito_b = True 
    else:
        circuito_b = False 
    return render_template('RLC.html', activar_home = False, activar_panel = activar_panel, vi = vi, Ii = Ii, R = R, C=C,L=L,circuito=circuito, circuito_b = circuito_b)

@app.route("/matplot-RLC-<float:vi>-<float:Ii>-<float:R>-<float:C>-<float:L>-<int:circuito>.svg")
def plot_svg_RLC(vi=10, Ii=-1, R=100,C=0.002,L=3, circuito=1):
    """ renders the plot on the fly.
    """
    svg_figure = SvgCanvasRLC(vi, Ii, R,C,L,circuito)

    return Response(svg_figure.drawPlot(), mimetype="image/svg+xml")

@app.route('/espejos', methods = ['GET', 'POST'])
def espejos():
    activar_panel = [
         False,
         True,
         False,
         False,
         False,
         False]
    so = float(request.args.get("rng_so", 10.0))
    ho = float(request.args.get("rng_ho", 5.0))
    Ra = float(request.args.get("rng_Ra" , 8.0))
    espejo = int(request.args.get("rng_espejo",1))
    if espejo == 1:
        espejo_b = True
    else:
        espejo_b =  False
    return render_template('Espejos.html', activar_home = False, activar_panel = activar_panel,so=so,ho=ho,Ra=Ra,espejo=espejo, espejo_b=espejo_b)

@app.route("/matplot-<float:so>-<float:ho>-<float:Ra>-<int:espejo>.svg")
def plot_svg_Espejo(so=10, ho=1, Ra=100,espejo=1):
    """ renders the plot on the fly.
    """
    svg_figure = SvgCanvasEspejo(so, ho, Ra,espejo)
    return Response(svg_figure.drawPlot(), mimetype="image/svg+xml")


@app.route('/CL', methods = ['GET', 'POST'])
def CL():
    activar_panel = [
         False,
         False,
         False,
         False,
         True,
         False]
    angulo = float(request.args.get("rng_angulo", 45.0))
    vo = float(request.args.get("rng_vo", 50.0))
    yi = float(request.args.get("rng_yi" , 48.0))
    masa = float(request.args.get("rng_masa",100.0))
    
    return render_template('Cl.html', activar_home = False, activar_panel = activar_panel,angulo=angulo,vo=vo,yi=yi,masa=masa)

@app.route("/matplotCL-<float:angulo>-<float:vo>-<float:yi>-<float:masa>.svg")
def plot_svg_CL(angulo=45, vo=40, yi=100,masa=100):
    """ renders the plot on the fly.
    """
    svg_figure = SvgCanvasCL(angulo, vo, yi,masa)
    return Response(svg_figure.drawPlot(), mimetype="image/svg+xml")
#***************************CLASES SIMULACIONES************************************************************
#REFLEXION Y REFRACCION
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
        lIncRay= lRefrRay = 2.8
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
        self.ax.set_xlabel('Ángulo de refracción= '+str(round(self.Ang_refr,2)))
        
        #self.figures = Figure()
        #self.axes = self.figures.add_subplot(1, 1, 1)
        #self.axes.plot([1,2,3,4], [1,2,3,4])
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        #FigureCanvasSVG(self.figures).print_svg(output)
        return output.getvalue()
#LENTES
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
#CIRCUITOS
class SvgCanvasRLC:
    def __init__(self, vi = 10, Ii = -1, R = 100,C=0.002,L=3, circuito=1):
        self.fig = Figure()
        self.vi = vi
        self.Ii=Ii-2
        self.I_pi = 0
        self.v_pi = 0
        self.R = R
        self.C=C
        self.L=L
        self.circuito=circuito
        self.ax = self.fig.add_subplot(1,1,1)
        self.drawPlot()
    def drawPlot(self):
        self.ax.clear()
        if(self.circuito==1):
            h=1e-4
            self.I_pi=-(self.Ii*self.R+self.vi)/self.L
            I_pi=self.I_pi
            Ii=self.Ii
            R=self.R
            C=self.C
            L=self.L
            corriente=[Ii]
            ti=0
            tiempo=[ti]

            while(ti<3):
                k1=h*I_pi
                l1=h*(-((R/L)*I_pi+(1/(L*C))*Ii))
                k2=h*(I_pi+0.5*l1)
                l2=h*(-((R/L)*(I_pi+0.5*k1)+(1/(L*C))*(Ii+0.5*l1)))
                k3=h*(I_pi+0.5*l2)
                l3=h*(-((R/L)*(I_pi+0.5*k2)+(1/(L*C))*(Ii+0.5*l2)))
                k4=h*(I_pi+l3)
                l4=h*(-((R/(L))*(I_pi+k3)+(1/(L*C))*(Ii+0.5*l3)))
                       
                If=Ii+(1/6)*(k1+2*k2+2*k3+k4)
                I_pf=I_pi+(1/6)*(l1+2*l2+2*l3+l4)
                tf=ti+h
                tiempo.append(tf)
                ti=tf
                Ii=If
                corriente.append(Ii)
                I_pi=I_pf
            self.ax.plot(tiempo,corriente)
            self.ax.set_xlabel('Tiempo (s)')
            self.ax.set_ylabel('Corriente I(A)')
            self.ax.set_title('Gráfica de corriente')
        elif (self.circuito==2):
            vi=self.vi
            h=1e-4
            self.v_pi=-(1/self.C)*((vi/self.R)+self.Ii)
            v_pi=self.v_pi
            R=self.R
            C=self.C
            L=self.L
            voltaje=[vi]
            ti=0
            tiempo=[ti]
            while(ti<4):
                k1=h*v_pi
                l1=h*(-((1/(R*C))*v_pi+(1/(L*C))*vi))
                k2=h*(v_pi+0.5*l1)
                l2=h*(-((1/(R*C))*(v_pi+0.5*k1)+(1/(L*C))*(vi+0.5*l1)))
                k3=h*(v_pi+0.5*l2)
                l3=h*(-((1/(R*C))*(v_pi+0.5*k2)+(1/(L*C))*(vi+0.5*l2)))
                k4=h*(v_pi+l3)
                l4=h*(-((1/(R*C))*(v_pi+k3)+(1/(L*C))*(vi+0.5*l3)))

                vf=vi+(1/6)*(k1+2*k2+2*k3+k4)
                v_pf=v_pi+(1/6)*(l1+2*l2+2*l3+l4)
                tf=ti+h
                tiempo.append(tf)
                ti=tf
                vi=vf
                voltaje.append(vi)
                v_pi=v_pf
            self.ax.plot(tiempo,voltaje)
            self.ax.set_xlabel('Tiempo (s)')
            self.ax.set_ylabel('Voltaje (v)')
            self.ax.set_title('Gráfica de Voltaje')
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        return output.getvalue()
#ESPEJOS
class SvgCanvasEspejo:
    def __init__(self, so = 10, ho = 1, Ra = 8,espejo=1):
        self.fig = Figure()
        self.espejo= espejo
        self.so = so
        self.ho = ho
        self.hi = 0
        self.Ra  = Ra
        self.si = 0
        self.f  = 0
        self.ax = self.fig.add_subplot(1,1,1)
        self.drawPlot()
    def drawPlot(self):
        self.ax.clear()
        if(self.espejo==1):
            centro= self.Ra
            ho=self.ho
            so=self.so
            f=self.Ra/2
            self.si=1/((1/f)-(1/self.so))
            si=self.si
            Mt=-si/self.so
            hi=ho*Mt
            a=si/f
            if(abs(so)>abs(si)):
                x=np.linspace(-self.so-10,self.so+10)
            else:
                x=np.linspace(-abs(si)-10,abs(si)+10)
            #eje optico
            eje_optico= mlines.Line2D([x[0],x[-1]], [0, 0], color='k')
            self.ax.text(x[-1]-10,-2,'Eje óptico')
            self.ax.add_line(eje_optico)

            if(abs(ho)<abs(hi)):
                Y=np.linspace(-abs(hi)-5,abs(hi)+5)
            else:
                Y= np.linspace(-abs(ho)-5,abs(ho)+5)
        
            vrt=mlines.Line2D([0,0],[Y[0],Y[-1]],ls='dashed',color='k')
            self.ax.add_line(vrt)
            #espejo
            mirr = mpatches.Arc([-centro,0],2*centro, 2*centro, angle=0, theta1=280, theta2=80, color='blue')
            self.ax.add_patch(mirr)

            #rayo 1 paralelo al eje y se refleja a traves de foco
            rayo_1 = mlines.Line2D([-self.so,0], [ho,ho], color='yellow',label='Rayo 1')
            self.ax.arrow(0,ho,-a*f,-a*ho,length_includes_head = True,head_width =0.5,color='yellow',shape= 'full')
            self.ax.add_line(rayo_1) 
        
            #rayo2 
            rayo_2= mlines.Line2D([-self.so,0],[ho,hi],color='red',label='Rayo 2')
            self.ax.add_line(rayo_2)
            if (so>f):
                self.ax.arrow(0,hi,-si,0,length_includes_head = True,head_width =0.5,color='red',shape= 'full')
            elif (so<f):
                r2=mlines.Line2D([-so,-f],[ho,0],ls='dashed',color='red')
                self.ax.add_line(r2)
                r1=mlines.Line2D([-f,0],[0,ho],ls='dashed',color='yellow')
                self.ax.add_line(r1)
                self.ax.arrow(0,hi,-si,0,length_includes_head = True,head_width =0.5,color='red',shape= 'full')
        
            #rayo3
            if(so>centro):
                y=hi* centro/(centro-si)
                rayo_3 = mlines.Line2D([-self.so,0],[ho,y],color='green',label='Rayo 3')
                self.ax.add_line(rayo_3)
                self.ax.arrow(0,y,-si,(-y+hi),length_includes_head = True, head_width = 0.5, color='green',shape='full')
            elif(so<centro):
                y=ho* centro/(centro-so)
                if(so>f):
                    rayo_3 = mlines.Line2D([-self.so,0],[ho,y],ls='dashed',color='green',label='Rayo 3')
                    self.ax.add_line(rayo_3)
                    self.ax.arrow(-so,ho,-si+so,hi-ho,length_includes_head= True,head_width=0.5,color='green',shape='full')

                elif(so<f):
                    rayo_3 = mlines.Line2D([-centro,-self.so],[0,ho],ls='dashed',color='green',label='Rayo 3')        
                    self.ax.add_line(rayo_3)
                    r3= mlines.Line2D([-so,0],[ho,y],color='green')
                    self.ax.add_line(r3)
                    self.ax.arrow(0,y,-si,hi-y,length_includes_head= True,head_width=0.5,color='green',shape='full')
                
            #self.ax.axis('equal')
            #objeto
            self.ax.arrow(-self.so, 0,0, ho, length_includes_head=True, head_width=0.5, color='k', shape='full') #objeto
            self.ax.text(-self.so-2.5,0.5*ho,'ho')
        
            self.ax.scatter(-centro,0,color='blue',label='centro')
            #self.ax.text(-centro,-2.5,'C',label='centro')
            self.ax.scatter(-f,0,color='k',label='foco')
            #self.ax.text(-f-0.5,-2.5,'f')
            #imagen
            self.ax.arrow(-self.si,0,0,hi,length_includes_head = True, head_width=0.5,color='k',shape='full') #imagen
            self.ax.text(-self.si+1,0.5*hi,'hi')

            self.ax.set_xlim(x[0]-3,x[-1]+3)
            self.ax.set_ylim(Y[0]-3,Y[-1]+3)

            self.ax.set_xticklabels([])
            self.ax.set_yticklabels([])
            self.ax.set_xlabel('si= '+ str(round(si,2))+'cm  ,'+'hi = '+str(round(hi,2))+'cm  ,'+'foco= '+str(round(f,2))+'cm  ,'+'Mt= '+ str(round(Mt,2)))

            self.ax.legend(loc='lower right')
   
        elif (self.espejo==2):
            centro= self.Ra
            ho=self.ho
            so=self.so
            f=-self.Ra/2
            self.si=1/((1/f)-(1/so))
            si=self.si
            Mt=-si/so
            hi=ho*Mt
            a=si/f
            x=np.linspace(-self.so-10,centro+5)
            
            #eje optico
            eje_optico= mlines.Line2D([x[0],x[-1]], [0, 0], color='k')
            self.ax.text(x[-1]-10,-2,'Eje óptico')
            self.ax.add_line(eje_optico)

            Y=np.linspace(-abs(ho)-5,abs(ho)+5)
            
            vrt=mlines.Line2D([0,0],[Y[0],Y[-1]],ls='dashed',color='k')
            self.ax.add_line(vrt)
            #espejo
            mirr = mpatches.Arc([centro,0],2*centro, 2*centro, angle=0, theta1=80, theta2=280, color='blue')
            self.ax.add_patch(mirr)

            #rayo 1 paralelo al eje y se refleja a traves de foco
            rayo_1 = mlines.Line2D([-self.so,0], [ho,ho], color='yellow',label='Rayo 1')
            r1 = mlines.Line2D([0,-f], [ho,0],ls='dashed', color='yellow')
            self.ax.arrow(0,ho,f,ho,length_includes_head = True,head_width =0.5,color='yellow',shape= 'full')
            self.ax.add_line(rayo_1) 
            self.ax.add_line(r1)
        
            #rayo2 
            rayo_2= mlines.Line2D([-self.so,0],[ho,hi],color='red',label='Rayo 2')
            self.ax.add_line(rayo_2)
            r2=mlines.Line2D([0,-f],[hi,0],ls='dashed',color='red')
            self.ax.add_line(r2)
            r1=mlines.Line2D([0,-si],[hi,hi],ls='dashed',color='red')
            self.ax.add_line(r1)
            self.ax.arrow(0,hi,-so,0,length_includes_head = True,head_width =0.5,color='red',shape= 'full')
        
            #rayo3
            y=hi* centro/(centro+si)
            rayo_3 = mlines.Line2D([-self.so,0],[ho,y],color='green',label='Rayo 3')
            self.ax.add_line(rayo_3)
            self.ax.arrow(0,y,-0.5*so,0.5*(ho-y),length_includes_head = True, head_width = 0.5, color='green',shape='full')
            r3 = mlines.Line2D([0,centro],[y,0],ls='dashed' ,color='green')
            self.ax.add_line(r3)
            #self.ax.axis('equal')
            #objeto
            self.ax.arrow(-self.so, 0,0, ho, length_includes_head=True, head_width=0.5, color='k', shape='full') #objeto
            self.ax.text(-self.so-2.5,0.5*ho,'ho')
        
            self.ax.scatter(centro,0,color='blue',label='centro')
            #self.ax.text(-centro,-2.5,'C',label='centro')
            self.ax.scatter(-f,0,color='k',label='foco')
            #self.ax.text(-f-0.5,-2.5,'f')
            #imagen
            self.ax.arrow(-self.si,0,0,hi,length_includes_head = True, head_width=0.5,color='k',shape='full') #imagen
            self.ax.text(-self.si+1,0.5*hi,'hi')

            self.ax.set_xlim(x[0]-3,x[-1]+3)
            self.ax.set_ylim(Y[0]-3,Y[-1]+3)

            self.ax.set_xticklabels([])
            self.ax.set_yticklabels([])
            self.ax.set_xlabel('si= '+ str(round(si,2))+'cm  ,'+'hi = '+str(round(hi,2))+'cm  ,'+'foco= '+str(round(f,2))+'cm  ,'+'Mt= '+ str(round(Mt,2)))

            self.ax.legend(loc='lower right')
        
        
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        return output.getvalue()
#CAIDA LIBRE
class SvgCanvasCL:
    def __init__(self,angulo=40,vo=30,yi=25,masa=1):
        self.fig =Figure()
        self.yi=yi
        self.vi=vo
        self.angulo= angulo
        self.masa=masa
        self.ax= self.fig.add_subplot(1,1,1)
        self.drawPlot()
    def drawPlot(self):
        self.ax.clear()
        angulo = (self.angulo*np.pi)/180
        dt = 0.001
        g = -9.8
        m = self.masa
        v_i = self.vi
        y0=self.yi
        v_x = v_i*np.cos(angulo)
        v_y = v_i*np.sin(angulo)
        Vy=[v_y]
        Vx=[v_x]
        aceleracion=[-g]
        x0 = 0
        vg = 0
        yy = [y0]
        xx = [x0]
        aux=10
        
        while(aux>0):
            vg += m*g * dt
            Vy.append(v_y+vg)
            Vx.append(v_x)
            aceleracion.append(-g)
            y0 += (vg + v_y) * dt
            x0 += v_x * dt
            if y0 > 0:
                yy.append(y0)
                xx.append(x0)
            else:
                break
        ymax=yy[0]
        cx=0
        for i in range(len(yy)-1):
            if(yy[i]<yy[i+1]):
                ymax=yy[i+1]
            else:
                ymax=yy[i]
                cx=xx[i]
                break
        tf=xx[-1]/v_x
        tiempo= np.linspace(0,tf,num=len(Vy))
        """
        self.ax.plot(tiempo,Vy,color='green')
        self.ax.plot(tiempo,Vx)
        self.ax.plot(tiempo, aceleracion)
        """
        alcance= mlines.Line2D([xx[0],xx[-1]], [0, 0],ls='dashed', color='k',label='Alcance máximo= '+str(round(xx[-1],2))+'m')
        self.ax.add_line(alcance)
        altura_max= mlines.Line2D([cx,cx],[0,ymax],ls='dashed',color='brown',label='Altura máxima= '+str(round(ymax,2))+'m')
        self.ax.add_line(altura_max)
        self.ax.scatter(xx[-1],yy[-1])
        self.ax.scatter(cx,ymax,color='red')
        self.ax.plot(xx, yy,color='blue')
        self.ax.set_xlabel('X (m)')
        self.ax.set_ylabel('Y (m)')
        self.ax.legend(loc='best')

              
        output = io.BytesIO()
        FigureCanvasSVG(self.fig).print_svg(output)
        return output.getvalue()
"""Concept-led PhotoCraft artwork: hybrid workflow, 3D chamber, energy saddle.

Geometry is authored as editable PhotoCraft paths. The 3D drawings are analytic
projections of explicitly defined spatial coordinates, exported as vectors.
"""
import math
import numpy as np

from gambar_photocraft import CircuitArtwork, INK, GRAY, PALE, WINE, PT, SCALE, knot

BLUE, GREEN, RED, AMBER = "#0288D1", "#059669", "#DC2626", "#F57C00"


class Artwork(CircuitArtwork):
    def __init__(self, width, height):
        super().__init__()
        self.width, self.height = width, height
        self.rect(0, 0, width, height, "#FFFFFF", None, name="Paper")

    def poly(self, points, fill=None, edge=INK, lw=.8, name="Geometry", dashed=False):
        path={"subpaths":[{"closed":True,"op":"combine",
                           "knots":[knot(self.xy(*p)) for p in points]}],"fillRule":"nonzero"}
        self.add("path",name,fill=fill,stroke=edge,lw=lw,dashed=dashed,path=path)

    def curve(self, points, color=INK, lw=1.2, dashed=False, name="Curve"):
        path={"subpaths":[{"closed":False,"op":"combine",
                           "knots":[knot(self.xy(*p)) for p in points]}],"fillRule":"nonzero"}
        self.add("path",name,stroke=color,lw=lw,dashed=dashed,path=path)

    def bezier(self, start, control1, control2, end, color=INK, lw=1.2, name="Bezier"):
        a,b=knot(self.xy(*start)),knot(self.xy(*end))
        a["out"],b["in"]=self.xy(*control1),self.xy(*control2)
        self.add("path",name,stroke=color,lw=lw,
                 path={"subpaths":[{"closed":False,"op":"combine","knots":[a,b]}],"fillRule":"nonzero"})

    def arrow(self, start, end, color=INK, lw=1.2, size=1.5, dashed=False):
        self.line(*start,*end,color,lw,dashed)
        self.tip(start,end,color,size)

    def tip(self, start, end, color=INK, size=1.5):
        a,b=np.asarray(start,dtype=float),np.asarray(end,dtype=float)
        direction=(b-a)/np.linalg.norm(b-a)
        normal=np.array([-direction[1],direction[0]])
        self.poly([b,b-size*direction+.42*size*normal,b-size*direction-.42*size*normal],color,None)

    def dot(self, x, y, r=.75, fill=INK, edge=None, lw=.8):
        self.add("ellipse","Point",fill=fill,stroke=edge,lw=lw,
                 rect=self.xy(x-r,y+r)+[2*r*PT*SCALE]*2)

    def rounded(self, x,y,w,h,fill,edge="none",radius=2):
        self.add("roundedRect","Computing domain",fill=fill,
                 stroke=None if edge=="none" else edge,lw=.8,
                 rect=self.xy(x,y+h)+[w*PT*SCALE,h*PT*SCALE],radii=radius*PT*SCALE)


def workflow():
    a=Artwork(170,50)
    a.rounded(1,4,77,44,"#F1F5F9","#CBD5E1")
    a.rounded(81,4,44,44,"#ECFDF5","#A7D8C3")
    a.rounded(128,4,41,44,"#F1F5F9","#CBD5E1")
    for x,title,color in ((39.5,"CLASSICAL COMPUTER",GRAY),(103,"QUANTUM COMPUTER",GREEN),(148.5,"CLASSICAL COMPUTER",GRAY)):
        a.label(x,44,title,7.3,color,True)
    a.label(5,39,"(a) Molecular model",7.6,bold=True,ha="left")
    a.label(39,39,"(b) Mapping & tapering",7.6,bold=True,ha="left")
    a.line(36.5,9,36.5,40,"#CBD5E1",.8)
    a.arrow((4,27),(34,27),BLUE,.8,1.2)
    for x in (8,15,22,29):
        a.dot(x,27,2.25,"#E1F5FE",BLUE,1.2)
        a.label(x,27,"H",7.7,BLUE)
    a.label(34,30,"$z$",7.7,BLUE)
    for x in (8,15,22):
        a.line(x,30,x,35,BLUE,.6)
        a.line(x,33,x+7,33,BLUE,.8)
        a.tip((x+7,33),(x,33),BLUE,.9)
        a.tip((x,33),(x+7,33),BLUE,.9)
        a.label(x+3.5,35.5,"$R$",7.5)
    a.line(29,30,29,35,BLUE,.6)
    a.label(19,19,"Real RHF orbitals",7.3,GRAY)
    a.label(19,14,r"$\hat H=\hat H^{*}$",9)
    a.label(19,8.5,r"H$_4$ / STO-3G",7.3,GRAY)
    top=np.linspace(41,73,8)
    bottom=np.linspace(45,70,5)
    for q,x in zip((2,3,5,6,7),bottom):
        a.line(top[q],27.4,x,24,"#94A3B8",.8,True)
    for i,x in enumerate(top):
        removed=i in (0,1,4)
        a.dot(x,29,1.6,"#FFEBEE" if removed else "#F0F4C3","#D32F2F" if removed else "#689F38",1)
        a.label(x,33.8,str(i),6.2,GRAY)
        if removed:
            a.line(x-1,28,x+1,30,"#D32F2F",1)
            a.line(x-1,30,x+1,28,"#D32F2F",1)
    for i,x in enumerate(bottom):
        a.dot(x,22.5,1.6,"#F0F4C3","#689F38",1)
        a.label(x,22.5,str(i),5.7)
    a.label(57.5,18,r"Jordan–Wigner + $\mathbb{Z}_2$",7.1)
    a.label(57.5,12.4,r"$H=\sum_k h_kP_k$",9)
    a.label(57.5,7.4,"8 spin orbitals → 5 qubits",6.8,GRAY)
    a.arrow((76,27),(87,27),GRAY,1.2)
    a.label(103,39,"(c) State preparation",7.6,bold=True)
    a.label(103,34.3,r"$|\Phi_{\rm HF}\rangle=|00100\rangle$",8)
    a.rounded(92,21,22,10,"#FFFFFF",GREEN,1.2)
    for y in np.linspace(22.5,29.5,5):
        a.line(86,y,120,y,GREEN,1.2)
    a.rect(98,24,10,4,"#FFFFFF",None)
    a.label(103,26,r"$U(\boldsymbol{\theta})$",9.5,GREEN)
    a.label(103,16.6,r"Measure $\langle P_k\rangle$",8,GREEN)
    a.label(103,10.5,"State preparation + readout",6.8,GRAY)
    a.arrow((121,26),(133,26),GRAY,1.2)
    a.label(148.5,39,"Energy + optimization",7.6,bold=True)
    a.arrow((135,22),(135,35.5),INK,.8,1.1)
    a.arrow((135,22),(165,22),INK,.8,1.1)
    a.line(135,24,165,24,GRAY,.8,True)
    t=np.linspace(0,1,70)
    a.curve([(136+27*x,24.4+9*math.exp(-5*x)) for x in t],AMBER,1.5)
    a.label(165.3,26,"$E_0$",7.2)
    a.label(151,33.5,r"$E(\theta)\geq E_0$",7.7)
    a.label(137,36.5,"Energy",6.1,GRAY)
    a.label(151,20.3,"Iteration",6.1,GRAY)
    a.label(149,15.4,r"$E=\sum_k h_k\langle P_k\rangle$",7.3)
    a.label(148.5,9.6,"L-BFGS-B",8,INK,True)
    a.bezier((147,8),(147,1),(84,1),(84,10),GREEN,1.2,"Classical parameter feedback")
    a.line(84,10,84,18,GREEN,1.2)
    a.bezier((84,18),(84,22.5),(88,22.5),(92,22.5),GREEN,1.2)
    a.tip((89,22.5),(92,22.5),GREEN,1.4)
    a.rect(124,.5,13,3.4,"#FFFFFF",None)
    a.label(130.5,2.2,r"Update $\boldsymbol{\theta}$",7.5,GREEN)
    return a.shapes


def weyl():
    a=Artwork(120,75)
    a.label(60,70.5,r"Folded Weyl chamber: $\pi/4\geq\alpha\geq\beta\geq|\gamma|$",9)
    # Orthographic camera at azimuth -55°, elevation 55°. Coordinates are
    # normalized by pi/4; projection uses orthonormal screen basis vectors.
    az,el=math.radians(-55),math.radians(55)
    right=np.array([-math.sin(az),math.cos(az),0])
    up=np.array([-math.sin(el)*math.cos(az),
                 -math.sin(el)*math.sin(az),math.cos(el)])
    def project(p):
        p=np.asarray(p)
        return np.array([15,35])+32.5*np.array([p@right,p@up])
    I,C,SP,SM,Q=[project(p) for p in ((0,0,0),(1,0,0),(1,1,1),(1,1,-1),(1,1,0))]
    # Distinct lit faces and spatial cross-sections give the volume its depth.
    a.poly([I,SP,SM],"#F4F7FA",None,name="Rear tetrahedron face")
    a.poly([C,SP,SM],"#DDE6EF",None,name="Alpha boundary face")
    a.poly([I,C,SP],[.88,.92,.96,.58],None,name="Upper chamber face")
    a.poly([I,C,SM],[.94,.96,.98,.7],None,name="Lower chamber face")
    for sign in (-1,1):
        for t in (.25,.5,.75):
            triangle=[project(p) for p in ((t,t,sign*t),(1,t,sign*t),(1,1,sign*t))]
            a.curve(triangle+[triangle[0]],"#BDCBD8",.45,False,"Constant gamma section")
    a.line(*I,*SM,"#8A99AA",.8,True)
    a.line(*SP,*SM,INK,1.2)
    a.poly([I,C,Q],[.99,.65,.65,.42],WINE,1.2,"Two-CNOT plane")
    for t in (.25,.5,.75):
        a.line(*project((t,0,0)),*project((t,t,0)),"#D18593",.45)
        a.line(*project((t,t,0)),*project((1,t,0)),"#D18593",.45)
    for x,y in ((I,C),(I,SP),(C,SP),(C,SM)):
        a.line(*x,*y,INK,1.2)
    for p in (I,C,SP,SM): a.dot(*p,.8,INK)
    a.dot(*Q,.9,"#FFFFFF",WINE,1.2)
    generic=project((.85,.68,.43))
    floor=project((.85,.68,0))
    a.line(*generic,*floor,GRAY,.8,True)
    a.dot(*generic,.85,GRAY)
    a.dot(*floor,.6,WINE)
    a.label(I[0]-3,I[1]+4,r"$I\otimes I$",9,ha="right")
    a.label(C[0]-1,C[1]-4.5,"CNOT",8)
    a.label(SP[0]+4,SP[1],"SWAP",8,ha="left")
    a.label(SM[0]+4,SM[1],"SWAP",8,ha="left")
    a.label(Q[0]+4,Q[1],r"$i$SWAP",8,ha="left")
    a.rounded(87,42,31,19,"#F1F5F9","#CBD5E1",1.7)
    a.label(102.5,56.5,r"$\gamma\ne0$",10,GRAY)
    a.label(102.5,50.5,"3 CNOTs",9.2,GRAY,True)
    a.label(102.5,45.5,"Generic CAN",7.5,GRAY)
    a.rounded(87,17,31,20,"#FFF1F2","#D9A9B5",1.7)
    a.label(102.5,32.8,r"$\gamma=0$",10,WINE)
    a.label(102.5,26.7,r"$\leq$ 2 CNOTs",9.2,WINE,True)
    a.label(102.5,21.2,"SO(4), up to local gates",6.6,WINE)
    a.bezier(tuple(generic),(76,52),(81,52),(86,52),GRAY,.8)
    a.bezier(tuple(project((.75,.35,0))),(70,26),(78,27),(86,27),WINE,.8)
    # Small 3D coordinate triad uses the same camera, giving spatial orientation.
    origin=np.array([10.,10.])
    for p,label,offset in (((1,0,0),r"$\alpha$",(1,-2)),((0,1,0),r"$\beta$",(2,1)),((0,0,1),r"$\gamma$",(0,2))):
        direction=np.array([np.dot(p,right),np.dot(p,up)])*9
        end=origin+direction
        a.arrow(origin,end,GRAY,.8,1.1)
        a.label(*(end+offset),label,8,GRAY)
    a.label(60,2.6,"Analytic 3D schematic; no sampled gates",7.1,GRAY)
    return a.shapes


def tangents():
    """Local real/imaginary state derivatives, not an energy landscape.

    A realification of the state-vector space is illustrated with two real
    directions and one representative imaginary direction. The displayed disk
    is a linear tangent subspace, not the full invariant parameter manifold.
    """
    a=Artwork(120,65)
    a.label(60,61.3,"Local tangent directions at a real circuit",9.4,bold=True)
    a.label(31,55.3,r"$\theta_e=0,\quad |\Psi\rangle=|\Psi\rangle^{*}$",8.6,GRAY)
    az,el=math.radians(-50),math.radians(28)
    right=np.array([-math.sin(az),math.cos(az),0.])
    up=np.array([-math.sin(el)*math.cos(az),
                 -math.sin(el)*math.sin(az),math.cos(el)])
    def project(x,y,z=0):
        point=np.array([x,y,z])
        return np.array([31.,32.])+23*np.array([point@right,point@up])
    origin=project(0,0)
    radius=1.18
    # The dashed continuation is behind the real tangent disk.
    a.line(*project(0,0,-.73),*origin,GRAY,.8,True)
    rim=[project(radius*math.cos(t),radius*math.sin(t)) for t in np.linspace(0,2*math.pi,180)]
    a.poly(rim,"#E0F2FE",None,name="Real tangent subspace")
    for t in (-.8,-.4,0,.4,.8):
        extent=math.sqrt(radius*radius-t*t)
        for start,end in (((t,-extent),(t,extent)),((-extent,t),(extent,t))):
            a.line(*project(*start),*project(*end),"#B3D8EA",.55)
    a.curve(rim+[rim[0]],BLUE,1.2,name="Boundary of illustrative tangent patch")
    # Projected arrows use a common orthographic 3D frame.
    imaginary=project(0,0,.88)
    real=project(.95,.55)
    a.arrow(origin,imaginary,"#334155",1.2,1.8)
    a.arrow(origin,real,RED,1.5,1.8)
    a.dot(*origin,1.1,"#FFFFFF",INK,1.2)
    a.label(imaginary[0]+3,imaginary[1]-1,r"$\tau_e$",10,GRAY,ha="left")
    a.label(real[0]-1,real[1]+3.2,r"$\tau_o$",10,RED)
    a.label(18.5,25.6,"Real state variations",7.2,BLUE)
    a.line(origin[0]-.9,origin[1]-.7,19,18.5,INK,.8)
    a.label(6,16.5,r"Real state $|\Psi\rangle$",8.2,ha="left")
    a.label(6,12.2,"Local state-space projection",7.1,GRAY,ha="left")
    # Derivative cards show why the energy responses differ. a and b are real
    # vectors; the signs and components depend on the chosen circuit generator.
    a.rect(69,34,49,21,PALE,GRAY,.8,True,name="Imaginary derivative mechanism")
    a.label(72,51.2,r"EVEN $Y$ · imaginary",8,GRAY,True,ha="left")
    a.label(93.5,45.1,r"$\tau_e=\partial_e|\Psi\rangle=i|b\rangle$",9,GRAY)
    a.label(93.5,38.8,r"$\partial E/\partial\theta_e=0$",10,GRAY)
    a.rect(69,11,49,20,"#FFFFFF",WINE,1.2,name="Real derivative mechanism")
    a.label(72,27.3,r"ODD $Y$ · real",8,WINE,True,ha="left")
    a.label(93.5,21.5,r"$\tau_o=\partial_o|\Psi\rangle=|a\rangle$",9,WINE)
    a.label(93.5,15.6,"Energy gradient may be nonzero",7.5,WINE)
    a.line(3,8,117,8,"#CBD5E1",.8)
    a.label(30,4.1,"CAN: 3 CNOTs at generic angles",7.3,GRAY)
    a.label(91,4.1,r"SO(4): $\leq$ 2 CNOTs",7.7,WINE)
    return a.shapes


def symmetry():
    a=Artwork(170,60)
    a.label(3,56,"(a) Conjugation symmetry",8.7,bold=True,ha="left")
    a.label(66,56,"(b) A saddle in parameter space",8.7,bold=True,ha="left")
    a.label(117,51.3,"Illustrative two-parameter energy landscape",7.1,GRAY)
    a.line(8,48,14,48,GRAY,1.2,True)
    a.label(16,48,r"Even $Y$: $E(v)=E(-v)$",7.6,GRAY,ha="left")
    a.line(8,43.4,14,43.4,WINE,1.2)
    a.label(16,43.4,r"Odd $Y$ (example)",7.6,WINE,ha="left")
    a.arrow((10,19),(10,40),INK,.8,1.1)
    a.arrow((10,19),(58,19),INK,.8,1.1)
    x=np.linspace(-math.pi,math.pi,120)
    sx=lambda v:33+22*v/math.pi
    sy=lambda e:29+8*e
    a.line(33,19,33,39,GRAY,.8,True)
    a.curve([(sx(v),sy(math.cos(v))) for v in x],GRAY,1.2,True)
    a.curve([(sx(v),sy(math.cos(v-.8))) for v in x],WINE,1.5)
    a.line(29,37,37,37,GRAY,.8)
    a.dot(33,37,.6,"#FFFFFF",GRAY,.8)
    for px,label in ((11,r"$-\pi$"),(33,"$0$"),(55,r"$\pi$")):
        a.line(px,18.2,px,19,INK,.8)
        a.label(px,16,label,8)
    a.label(33,11.8,r"Angle $\theta_k$",8)
    a.label(32.5,5.8,r"On $\mathcal{M}_0$: $\partial E/\partial v=0$",8.1,GRAY)
    # E(u,v) = u^2 - 0.85 v^2. u is one odd-Y direction and v one even-Y
    # direction. The green curve is the graph of E restricted to v=0.
    def energy(u,v): return u*u-.85*v*v
    def project(u,v,e=None):
        e=energy(u,v) if e is None else e
        return np.array([122+15*u+18*v,29+3*u-7*v+8*e])
    # A real three-axis frame; z is energy, not another parameter coordinate.
    origin=np.array([78.,17.])
    end_u=origin+np.array([17.,3.4])
    end_v=origin+np.array([18.,-7.])
    end_e=origin+np.array([0.,17.])
    for end in (end_u,end_v,end_e): a.arrow(origin,end,"#8FA0B2",.8,1.2)
    a.label(end_u[0],end_u[1]+2.5,r"$u\;(\theta_o)$",8.3,GRAY)
    a.label(end_v[0],end_v[1]-2.7,r"$v\;(\theta_e)$",8.3,GRAY)
    a.label(end_e[0]-1.8,end_e[1]+2.5,"$E$",10,GRAY)
    # Depth-sort a shaded saddle mesh. Each facet is a native vector shape.
    grid_u=np.linspace(-1.2,1.2,19)
    grid_v=np.linspace(-1,1,19)
    facets=[]
    for u0,u1 in zip(grid_u[:-1],grid_u[1:]):
        for v0,v1 in zip(grid_v[:-1],grid_v[1:]):
            u,v=(u0+u1)/2,(v0+v1)/2
            e=energy(u,v)
            t=max(0,min(1,(e+.85)/2.29))
            low=np.array([.76,.86,.93]); high=np.array([.99,.78,.71])
            color=low*(1-t)+high*t
            normal=np.array([-2*u,1.7*v,1.])
            normal/=np.linalg.norm(normal)
            light=.88+.12*max(0,np.dot(normal,np.array([-.3,-.4,.866])))
            color=np.clip(color*light,0,1)
            fill="#"+"".join(f"{round(c*255):02x}" for c in color)
            depth=-.59*u+.49*v+.65*e
            corners=[project(uu,vv) for uu,vv in ((u0,v0),(u1,v0),(u1,v1),(u0,v1))]
            facets.append((depth,corners,fill))
    for _,corners,fill in sorted(facets,key=lambda f:f[0]):
        a.poly(corners,fill,fill,.12,name="Analytic saddle facet")
    for u in np.linspace(-1.2,1.2,9):
        a.curve([project(u,v) for v in np.linspace(-1,1,60)],"#9DADB9",.35)
    for v in np.linspace(-1,1,9):
        a.curve([project(u,v) for u in np.linspace(-1.2,1.2,60)],"#9DADB9",.35)
    a.curve([project(u,0) for u in np.linspace(-1.2,1.2,90)],GREEN,1.2,name="M0 energy section")
    steps=[1.15,.88,.63,.42,.24,.11,0]
    a.curve([project(u,0) for u in np.linspace(steps[0],0,80)],GREEN,2,name="Constrained optimizer trajectory")
    for u in steps: a.dot(*project(u,0),.7,"#FFFFFF",GREEN,1.1)
    a.tip(project(.2,0),project(.11,0),GREEN,1.6)
    saddle=project(0,0)
    a.dot(*saddle,1,"#FFFFFF",INK,1.2)
    a.arrow(saddle,project(0,.88),RED,1.1,1.5,True)
    a.arrow(saddle,project(0,-.88),RED,1.1,1.5,True)
    a.bezier(tuple(project(.72,0)),(149,44),(155,45),(161,45),GREEN,.8)
    a.label(162,48,r"$\mathcal{M}_0$ slice",8.5,GREEN,ha="right")
    a.label(163,43,"Exact-gradient path",7,GREEN,ha="right")
    a.label(156,12.8,"Negative curvature",7.2,RED)
    a.line(*project(0,.83),149,14.8,RED,.8)
    a.label(119,3.2,r"$v=0$: constrained optimum; a transverse direction can still lower $E$",7.2,GRAY)
    return a.shapes


BUILDERS={"fig1_vqe_pipeline":workflow,"fig5_weyl_chamber":weyl,
          "fig6_tangent_directions":tangents,"fig7_symmetry_manifold":symmetry}

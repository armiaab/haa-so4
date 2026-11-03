"""Reproducible English schematics; PDF for LaTeX, editable SVG for Affinity.

Run with Python, NumPy and Matplotlib. These are mathematical illustrations,
not molecular simulation data. No random samples are used in the Weyl diagram.
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse, Polygon, Rectangle

OUT = Path(__file__).resolve().parents[1] / 'figures' / 'jctc'
OUT.mkdir(exist_ok=True)
INK, BLUE, RED, GRAY, GREEN = '#222222', '#005A8D', '#B23A28', '#606060', '#007766'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'mathtext.fontset': 'stixsans',
                     'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.linewidth': 0.7,
                     'text.color': INK, 'axes.labelcolor': INK, 'savefig.dpi': 600})


def canvas(height=2.5, width=6.6):
    fig, ax = plt.subplots(figsize=(width, height))
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.04, top=0.96)
    ax.set(xlim=(0, 10), ylim=(0, 4)); ax.axis('off')
    return fig, ax


def text(ax, x, y, label, **kw):
    ax.text(x, y, label, ha='center', va='center', **kw)


def arrow(ax, a, b, color=INK, **kw):
    ax.annotate('', xy=b, xytext=a,
                arrowprops=dict(arrowstyle='-|>', color=color, lw=0.9, **kw))


def box(ax, x, y, w, h, label, color=BLUE, fill='#F0F5F8'):
    ax.add_patch(Rectangle((x, y), w, h, ec=color, fc=fill, lw=0.9, zorder=3))
    text(ax, x+w/2, y+h/2, label, color=color, zorder=4)


def save(fig, name):
    for ext in ('pdf', 'svg'):
        fig.savefig(OUT / f'{name}.{ext}', metadata={'Creator': 'src/gambar_skema.py'})
    fig.savefig(OUT / f'{name}.png', dpi=180)
    plt.close(fig)


def workflow():
    fig, ax = canvas(2.35)
    for x, title in [(0.1, '(a) Molecular model'), (3.4, '(b) Qubit Hamiltonian'), (6.8, '(c) Variational optimization')]:
        ax.text(x, 3.7, title, fontsize=9, weight='bold')
    xs = np.linspace(0.4, 2.4, 4)
    ax.plot(xs, [2.3]*4, color=GRAY, lw=1)
    for x in xs:
        ax.plot(x, 2.3, 'o', ms=19, mfc='white', mec=BLUE, mew=1)
        text(ax, x, 2.3, 'H', color=BLUE)
    for x in (xs[:-1]+xs[1:])/2:
        text(ax, x, 2.85, '$R$')
    text(ax, 1.4, 1.15, 'Real RHF orbitals\nSTO-3G')
    arrow(ax, (2.65, 2.3), (3.25, 2.3))
    box(ax, 3.4, 2.0, 2.7, 0.85, 'Jordan-Wigner: 8 qubits')
    arrow(ax, (4.75, 2.0), (4.75, 1.6))
    box(ax, 3.4, 0.7, 2.7, 0.85, r'$Z_2$ tapering: 5 qubits')
    text(ax, 4.75, 0.25, r'Real symmetric $H$')
    arrow(ax, (6.15, 2.3), (6.65, 2.3))
    box(ax, 6.8, 2.0, 2.9, 0.85, r'$|\Phi_{\rm HF}\rangle\ \to\ U(\theta)\ \to\ E(\theta)$')
    arrow(ax, (8.25, 2.0), (8.25, 1.6))
    box(ax, 6.8, 0.7, 2.9, 0.85, 'L-BFGS-B: update angles', GREEN, '#EFF6F3')
    arrow(ax, (9.7, 1.1), (9.7, 2.0), GREEN, connectionstyle='arc3,rad=-0.7')
    text(ax, 8.25, 0.25, r'$E(\theta)\geq E_0$')
    save(fig, 'alur')


def ancilla():
    fig, ax = canvas(2.3)
    for y in (3.2, 2.7, 2.2):
        ax.plot([2, 7.5], [y, y], color=GRAY, lw=0.9)
    ax.plot([2, 7.0], [0.9, 0.9], color=RED, lw=1, ls='--')
    box(ax, 3.4, 0.6, 2.2, 2.9, r'$U(\boldsymbol{\theta})$')
    text(ax, 1.0, 2.75, r'$|\Phi_{\rm HF}\rangle$'+'\nsystem qubits')
    text(ax, 1.0, 0.9, r'$|0\rangle$'+'\nancilla', color=RED)
    text(ax, 8.55, 2.8, r'$\rho_s={\rm Tr}_a|\Psi\rangle\langle\Psi|$')
    text(ax, 8.55, 2.15, 'Generally mixed', fontsize=8)
    box(ax, 6.6, 0.55, 1.1, 0.7, r'${\rm Tr}_a$', RED, '#FBEFEB')
    text(ax, 5.0, 0.15, r'$E={\rm Tr}(H\rho_s)=\langle\Psi|H\otimes I_a|\Psi\rangle\geq E_0$')
    save(fig, 'ancilla')


def parity():
    fig, ax = canvas(3.0)
    names = 'IXYZ'
    for r, p in enumerate(names):
        for c, q in enumerate(names):
            x, y = 1.1+c*0.8, 3.1-r*0.8
            odd = (p+q).count('Y') % 2
            ax.add_patch(Rectangle((x,y-0.4),0.8,0.8,ec=GRAY,lw=0.7,
                                   fc='#FBEFEB' if odd else 'white'))
            text(ax,x+0.4,y,p+q,fontsize=10)
            if p+q in ('XX','YY','ZZ'):
                ax.add_patch(Rectangle((x+0.12,y-0.28),0.56,0.56,fill=False,ec=GRAY,lw=1.3))
            if p+q in ('XY','YX','YZ','ZY'):
                ax.add_patch(Ellipse((x+0.4,y),0.65,0.58,fill=False,ec=RED,lw=1.3))
    for i,p in enumerate(names):
        text(ax,1.5+i*0.8,3.72,p,weight='bold')
        text(ax,0.78,3.1-i*0.8,p,weight='bold')
    ax.text(0.15,2.0,'First qubit',rotation=90,va='center',fontsize=8)
    text(ax,2.7,3.99,'Second qubit',fontsize=8)
    for y,t,c in [(3.35,'Odd number of Y factors',RED),(2.88,r'$R_P(\theta)$ real for all angles',RED),
                  (2.2,'Squares: CAN generators',GRAY),(1.7,'Circles: SO(4) block generators',RED),
                  (1.15,r'$IY,\ YI$: local $R_Y$ rotations',BLUE),
                  (0.5,'Unshaded: complex at generic angles',INK)]:
        ax.text(4.9,y,t,va='center',color=c,fontsize=8.5)
    save(fig, 'parity')


def weyl():
    fig, ax = canvas(3.4)
    # Projection of (0,0,0), (a,0,0), (a,a,+a), (a,a,-a), a=pi/4.
    pts = np.array([[1,1.5],[5.0,0.65],[8.2,3.25],[8.2,0.45]])
    mid = (pts[2]+pts[3])/2
    for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
        ax.plot(*pts[[i,j]].T,color=GRAY,lw=0.9,ls='--' if (i,j)==(0,3) else '-')
    ax.add_patch(Polygon([pts[0],pts[1],mid],fc='#FBEFEB',ec=RED,lw=1.2))
    for p in pts:
        ax.plot(*p,'o',ms=4,color=GRAY)
    ax.plot(*mid,'s',ms=4,color=RED)
    for p,t,dx,dy in [(pts[0],r'$I\otimes I$',-.25,.17),(pts[1],'CNOT',0,-.25),
                     (pts[2],'SWAP',.65,0),(pts[3],'SWAP',.65,0),(mid,r'$i$SWAP',.65,0)]:
        text(ax,p[0]+dx,p[1]+dy,t)
    text(ax,5.2,1.37,r'$\gamma=0$: at most 2 CNOTs',color=RED)
    text(ax,4.8,3.65,r'Folded chamber: $\pi/4\geq\alpha\geq\beta\geq|\gamma|$',fontsize=10)
    text(ax,5.9,2.7,r'Generic $\gamma\ne0$: 3 CNOTs',color=GRAY)
    text(ax,4.9,0.05,'Analytic schematic; no sampled gates',fontsize=8)
    save(fig, 'weyl')


def concept():
    fig, ax = canvas(2.7)
    ax.add_patch(Ellipse((4.7,1.7),4.7,1.3,ec=BLUE,fc='#EDF4F8',lw=1))
    text(ax,4.8,1.25,'Real-amplitude directions',color=BLUE)
    ax.plot(3.0,1.8,'o',mfc='white',mec=INK,ms=6,zorder=5)
    arrow(ax,(3.0,1.8),(3.0,3.15),GRAY)
    arrow(ax,(3.0,1.8),(6.1,1.95),RED)
    text(ax,1.4,1.4,r'Real circuit at $\mathcal{M}_0$',fontsize=8)
    text(ax,3.8,3.65,'Even-Y tangent: imaginary',color=GRAY)
    text(ax,3.8,3.22,r'$\partial E/\partial\theta_e=0$',color=GRAY)
    text(ax,7.7,2.5,'Odd-Y tangent: real',color=RED)
    text(ax,7.7,2.05,'Gradient may be nonzero',color=RED,fontsize=8)
    text(ax,2.5,0.45,'CAN: 3 CNOTs at generic angles',color=GRAY,fontsize=8)
    text(ax,7.1,0.45,'SO(4): at most 2 CNOTs',color=RED,fontsize=8)
    save(fig, 'concept')


def manifold():
    fig, axs = plt.subplots(1,2,figsize=(6.6,2.7),gridspec_kw={'width_ratios':[1,1.15]})
    fig.subplots_adjust(left=.08,right=.98,bottom=.23,top=.86,wspace=.32)
    x=np.linspace(-np.pi,np.pi,160); ax=axs[0]
    ax.plot(x,np.cos(x),color=GRAY,lw=1.3,label='Even Y',ls='--')
    ax.plot(x,np.cos(x-.7),color=RED,lw=1.3,label='Odd Y (example)')
    ax.axvline(0,color=GRAY,lw=.7,ls=':')
    ax.plot([-.65,.65],[1,1],color=GRAY,lw=.8)
    ax.set(xticks=[-np.pi,0,np.pi],xticklabels=[r'$-\pi$','0',r'$\pi$'],yticks=[],
           xlabel=r'Angle $\theta_k$',ylabel='Schematic energy')
    ax.set_title('(a) All other even-Y angles zero',fontsize=9,loc='left')
    ax.legend(fontsize=7,frameon=False,loc='upper left')
    ax=axs[1]; ax.set(xlim=(0,10),ylim=(0,4));ax.axis('off')
    ax.set_title('(b) Invariant set in exact arithmetic',fontsize=9,loc='left')
    ax.add_patch(Polygon([(0.5,.8),(6.9,.8),(9.3,2.5),(2.9,2.5)],fc='#EFF6F3',ec=GREEN,lw=1))
    arrow(ax,(1.9,1.2),(1.9,3.6),GRAY)
    text(ax,2.7,3.85,r'$\theta_e$',color=GRAY)
    text(ax,8.4,.55,r'$\theta_o$',color=GREEN)
    ax.plot([2.5,4,5.6,6.7],[1.1,1.5,1.7,1.85],'-o',color=GREEN,lw=1,ms=3)
    text(ax,5.4,2.15,r'$\mathcal{M}_0:\ \theta_e=0$',fontsize=10)
    text(ax,5.1,0.0,'HAA-CAN: real product states',fontsize=8)
    save(fig,'m0')


def breakeven():
    fig,ax=plt.subplots(figsize=(4.5,3.0))
    fig.subplots_adjust(left=.15,right=.98,bottom=.19,top=.96)
    n=np.linspace(3,10,200); count=n/(n-1); depth=n/2
    ax.fill_between(n,.5,count,color='#FBEFEB')
    ax.fill_between(n,count,depth,color='#F7F3E7')
    ax.fill_between(n,depth,5.5,color='#EDF5F2')
    ax.plot(n,count,color=GREEN,lw=1.3,label=r'CNOT threshold $n_s/(n_s-1)$')
    ax.plot(n,depth,color=RED,lw=1.3,ls='--',label=r'Depth threshold $n_s/2$')
    ax.axvline(5,color=GRAY,lw=.7,ls=':')
    ax.plot(5,1.25,'o',color=GREEN,ms=4);ax.plot(5,2.5,'s',color=RED,ms=4)
    ax.text(5.08,1.4,'1.25',fontsize=8);ax.text(4.72,2.65,'2.5',fontsize=8)
    ax.text(6.1,.72,'More CNOTs and depth',fontsize=8)
    ax.text(6.0,2.0,'Fewer CNOTs; more depth',fontsize=8)
    ax.text(4.3,4.55,'Fewer CNOTs and depth',fontsize=8)
    ax.set(xlim=(3,10),ylim=(.5,5.5),xlabel=r'System qubits $n_s$',ylabel=r'Layer ratio $L_0/L_1$')
    ax.legend(fontsize=7,loc='upper left',frameon=False)
    save(fig,'break-even')


def gradient():
    fig,ax=canvas(2.0)
    for y in (1.4,1.9,2.4): ax.plot([.8,9.6],[y,y],color=GRAY,lw=.8)
    for x,w,l in [(1.3,1.3,'$W$'),(3.3,1.65,r'$R_{P_k}(\theta_k)$'),(6,1.3,'$V$'),(8.1,1.3,r'$\tilde H$')]:
        box(ax,x,1.05,w,1.7,l)
    ax.plot([5.4,5.4],[.8,3.15],color=GREEN,lw=1,ls='--')
    text(ax,5.4,3.55,r'$|\varphi_k\rangle$',color=GREEN)
    text(ax,.45,1.9,r'$|\psi_0\rangle$')
    text(ax,7.7,.4,r'$H_k=V^\dagger\tilde H V$')
    text(ax,4.1,3.5,r'$\theta_k$');arrow(ax,(4.1,3.25),(4.1,2.8))
    for label in ax.texts:
        label.set_fontsize(10)
    save(fig,'gradien')


def toc():
    fig, ax = canvas(1.75, 3.25)
    text(ax,5,3.6,r'Real molecular Hamiltonian: $H=H^*$',fontsize=9)
    box(ax,.5,1.35,3.7,1.35,'CAN\n3 CNOTs',GRAY,'#F3F3F3')
    box(ax,5.8,1.35,3.7,1.35,'SO(4)\nat most 2 CNOTs',RED,'#FBEFEB')
    arrow(ax,(4.35,2.0),(5.65,2.0),RED)
    text(ax,2.35,.96,'Generic angles',fontsize=7)
    text(ax,5,.4,r'$\mathrm{H_4}$   LiH   HF   $\mathrm{BeH_2}$',fontsize=9)
    # JCTC upload: TIF at 300 dpi or EPS, at the actual size of 3.25 in by 1.75 in
    fig.savefig(OUT / 'toc.tif', dpi=300)
    from PIL import Image
    Image.open(OUT / 'toc.tif').convert('RGB').save(OUT / 'toc.tif', dpi=(300, 300), compression='tiff_lzw')
    fig.savefig(OUT / 'toc.eps')
    save(fig,'toc')


if __name__ == '__main__':
    for draw in (workflow,ancilla,parity,weyl,concept,manifold,breakeven,gradient,toc):
        draw()
        print(draw.__name__)

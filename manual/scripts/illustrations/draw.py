"""Run from the manual root. Input files, figure sources and reproduction notes are alongside this script."""
import io,zipfile
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import map_coordinates
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
source=Path('scripts/illustrations/input');out=Path('figures');out.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#27343d','axes.edgecolor':'#71808a','axes.titlelocation':'left','axes.titleweight':'bold'})
ink='#243a49';teal='#00766e';ochre='#a14f1d';gray='#dce5e7'
def save(fig,name):fig.savefig(out/(name+'.png'),dpi=170,bbox_inches='tight',pad_inches=.12,facecolor='white');plt.close(fig)
def csv(z,f):
 with zipfile.ZipFile('archives/'+z) as a:return pd.read_csv(a.open(f))
# Reprojection uses the rotation and pixel convention of the archived tracing script.
G=np.array([[-.0548755604,-.8734370902,-.4838350155],[.4941094279,-.4448296300,.7469822445],[-.8676661490,-.1980763734,.4559837762]])
im=np.asarray(Image.open(source/'diffuse-map.png').convert('RGB'));H,W=im.shape[:2]
grid=np.fromfile(source/'classes.bin',dtype=np.uint8).reshape(481,1440)
colors=['#101722','#879a9e','#233d50','#d3ab6c','#fff0c3']
fig,axs=plt.subplots(4,2,figsize=(9.2,8.2),sharey=True,layout='constrained')
for col,centre in enumerate([0.,180.]):
 lon=np.arange(centre+90,centre-90-.01,-.25);lat=np.arange(-36,36.01,.25);L,B=np.meshgrid(lon,lat)
 l,b=np.deg2rad(L),np.deg2rad(B)
 xyz=np.stack([np.cos(l)*np.cos(b),np.sin(l)*np.cos(b),np.sin(b)],axis=-1)@G
 ra=np.arctan2(xyz[...,1],xyz[...,0]);dc=np.arcsin(xyz[...,2]);px=(.5-ra/(2*np.pi))*W-.5;py=(.5-dc/np.pi)*H-.5
 rgb=np.stack([map_coordinates(np.pad(im[:,:,k].astype(float),((0,0),(1,1)),mode='wrap'),[py,(px%W)+1],order=1,mode='nearest') for k in range(3)],axis=-1).clip(0,255).astype('uint8')
 layers=grid[np.rint((B+60)*4).astype(int),np.rint((L+180)*4).astype(int)%1440]
 for row,ax in enumerate(axs[:,col]):
  extent=[lon[0]+.125,lon[-1]-.125,lat[0]-.125,lat[-1]+.125]
  if row<2:ax.imshow(rgb,extent=extent,origin='lower',aspect='auto')
  elif row==2:ax.imshow(layers,extent=extent,origin='lower',aspect='auto',interpolation='nearest',cmap=ListedColormap(colors),vmin=-.5,vmax=4.5)
  else:ax.imshow(np.take([0.,1.,0.,1.35,2.1],layers),extent=extent,origin='lower',aspect='auto',interpolation='nearest',cmap='gray',vmin=0,vmax=2.1)
  if row==1:ax.contour(L,B,layers>0,levels=[.5],colors=['#61efdd'],linewidths=.95)
  ax.set_yticks([-30,0,30]);ax.set_xlim(lon[0],lon[-1]);ax.set_xticks([centre+90,centre+45,centre,centre-45,centre-90],[f'{v%360:.0f}°' for v in [centre+90,centre+45,centre,centre-45,centre-90]])
  if row<3:ax.tick_params(labelbottom=False)
  if col==0:ax.set_ylabel('Latitude (°)')
  if row==3:ax.set_xlabel('Galactic longitude')
  ax.set_title((['Source','Adopted outer edge','Reconstruction classes','Intrinsic diffuse radiance'][row])+(' · centre' if col==0 else ' · anticentre'),fontsize=10)
fig.legend([Patch(color=c) for c in colors],['Outside','Ordinary band','Dark lane','Bright','Brightest'],loc='outside lower center',ncol=5,frameon=False)
save(fig,'milky-way-reconstruction')
# Same synthetic surface as the DEM panorama example, sampled along due east.
d=np.arange(2.5,500.1,2.5);z=100+15*np.exp(-((d-80)/35)**2-(30/90)**2)+88*np.exp(-((d-270)/70)**2-(40/200)**2);eye=100+15*np.exp(-(80/35)**2-(30/90)**2)+88*np.exp(-(270/70)**2-(40/200)**2)+1.6
a=np.degrees(np.arctan2(z-eye,d));j=a.argmax()
fig,(ax,bx)=plt.subplots(2,1,figsize=(7.6,4.5),sharex=True,layout='constrained')
ax.fill_between(d,95,z,color='#c7dacf');ax.plot(d,z,color=teal);ax.plot([0,d[j]],[eye,z[j]],'--',c=ink,label='Sightline to skyline');ax.scatter([0,d[j]],[eye,z[j]],c=[ink,ochre],zorder=3);ax.set_ylabel('Elevation (m)');ax.legend(frameon=False);ax.set_title('Terrain along one eastward ray')
bx.plot(d,a,c=ink);bx.scatter(d[j],a[j],c=ochre);bx.annotate(f'Maximum {a[j]:.1f}° at {d[j]:.0f} m',(d[j],a[j]),xytext=(d[j]+45,a[j]-3),arrowprops={'arrowstyle':'->','color':ochre});bx.set(xlabel='Distance from observer (m)',ylabel='Geometric altitude (°)');save(fig,'terrain-ray')
# Actual synthetic mask and control from the archived worked example.
zname='conventions-field-panorama.zip';frame=csv(zname,'panorama-profile.csv').sort_values('az')
with zipfile.ZipFile('archives/'+zname) as z:mask=np.asarray(Image.open(io.BytesIO(z.read('skyline.png'))))
fig,(ax,bx)=plt.subplots(2,1,figsize=(7.6,4.3),layout='constrained');ax.set_facecolor('#eaf1f4');ax.imshow(mask,extent=[0,360,100,0],aspect='auto');ax.scatter(90,30,c='#e9a84c',s=36,edgecolor='black');ax.annotate('Angular control (90, 30)',(90,30),xytext=(125,16),arrowprops={'arrowstyle':'->'});ax.set(xlabel='Horizontal pixel coordinate',ylabel='Row from top',title='Input mask and control');bx.plot(frame.az,frame.uncalibrated-frame.truth,c=ochre,label='Nominal limits');bx.plot(frame.az,frame.calibrated-frame.truth,c=teal,label='Fitted limits');bx.axhline(0,c='#687278',lw=.7);bx.legend(ncol=2,frameon=False);bx.set(xlabel='True bearing (°)',ylabel='Altitude error (°)',xlim=(0,360));save(fig,'photo-calibration')
# Actual precomputed test, with the original local-cluster p-values.
s=pd.read_csv(source/'fit.csv');cl=pd.read_csv(source/'clusters.csv');fig,(ax,bx)=plt.subplots(2,1,figsize=(7.8,4.15),sharex=True,height_ratios=[3,1],layout='constrained')
ax.fill_between(s.x,s['CE.lower'],s['CE.upper'],color=gray,label='95% pointwise envelope');ax.plot(s.x,s.y,color=ink,label='Observed summed density');ax.fill_between(s.x,s['CE.upper'],s.y,where=s.y>s['CE.upper'],color='#c77d49',interpolate=True,label='Upper excursion');ax.set(title='Ġgantija · saved random-orientation test',ylabel='Summed density');ax.legend(frameon=False,fontsize=8,loc='upper right')
for j,v in cl.iterrows():
 bx.plot([v.start,v.end],[.4,.4],lw=5,c=teal if j==1 else ochre);bx.text((v.start+v.end)/2,.55,f'p = {v["p.value"]:.3f}',ha='center',fontsize=9)
bx.set(xlim=(-60,-14),xlabel='Declination (°)',ylim=(0,1),yticks=[],ylabel='Local excess');save(fig,'significance-clusters')
# Conditional normal mixture at the geometric, refraction-off settings of the archived doorway example.
A=np.arange(88,92.001,.02);H=np.interp(A,[88,89,89.2,90.8,91,92],[6,6,1,1,6,6]);rad=np.pi/180;phi=51*rad
def dec(h):return np.arcsin(np.sin(phi)*np.sin(h*rad)+np.cos(phi)*np.cos(h*rad)*np.cos(A*rad))/rad
lo,hi=dec(H-.8),dec(H+.8);mu=(hi+lo)/2;sd=abs(hi-lo)/8;D=np.arange(-2,8.01,.025);g=np.exp(-.5*((D[:,None]-mu)/sd)**2)/(sd*np.sqrt(2*np.pi))
door=csv('doorway-and-ridge-examples.zip','doorway.csv');fig,(ax,bx)=plt.subplots(1,2,figsize=(8,3.55),sharey=True,width_ratios=[1.7,1],layout='constrained');ax.pcolormesh(A,D,g,shading='auto',cmap='Blues');ax.plot(A,mu,color=ochre,lw=1.3,label='Conditional mean');ax.legend(frameon=False,loc='upper center');ax.set(xlabel='Bearing through opening (°)',ylabel='Declination (°)',title='Conditional density at each bearing');bx.plot(door.notch,door.dec,c=teal);bx.set(ylim=(D[0],D[-1]),xlabel='Density per degree',title='Integrated doorway density');save(fig,'conditional-density')
# Report the signed distribution with an inset-free, separately scaled tail table in the text.
e=pd.read_csv(source/'h01-bearings.csv').error.to_numpy();assert len(e)==334
fig,ax=plt.subplots(figsize=(7.6,2.8),layout='constrained');edges=np.linspace(-1.5,1.5,31);ax.hist(e,bins=edges,color=teal,edgecolor='white');ax.axvline(0,c=ink,lw=1);ax.axvline(-.5,c=ochre,ls='--');ax.axvline(.5,c=ochre,ls='--');ax.set(xlabel='Predicted minus measured altitude (°)',ylabel='Bearings',title='H01 · centre of the field residual distribution');ax.text(.98,.89,f'{sum(e< -1.5)} below −1.5°  ·  {sum(e>1.5)} above +1.5°',transform=ax.transAxes,ha='right',fontsize=9);save(fig,'h01-residuals')
# Two coordinate frames for the archived field reduction example.
fig,axs=plt.subplots(2,1,figsize=(7.6,3.2),layout='constrained')
for ax,start,end,ticks,title in zip(axs,[-1.25,157.73],[2.25,161.23],[[-3,-2,-1,0,1,2,3,4],[156,157,158,159,160,161,162,163]],['Instrument circle · crossing north','True bearing · referred to the solar centre']):
 ax.axhline(0,color=gray,lw=2)
 ax.annotate('',(end,.10),(start,.10),arrowprops={'arrowstyle':'->','color':ink,'lw':1.5})
 ax.scatter([start,end],[0,0],c=[teal,ochre],s=35)
 for x,label in [(start,'Solar sight'),(end,'Target')]:ax.text(x,-.20,f'{label}  {x%360:.2f}°',ha='center',fontsize=9)
 ax.text((start+end)/2,.2,'+3.50°',ha='center',fontsize=9)
 ax.set_xlim(ticks[0],ticks[-1]);ax.set_ylim(-.38,.45);ax.set_yticks([]);ax.set_xticks(ticks,[f'{t%360}°' for t in ticks]);ax.spines['left'].set_visible(False);ax.set_title(title,fontsize=10)
save(fig,'survey-reference')
# Analytic illustration of simultaneous criteria; these are constructed tracks, not ephemerides.
time=np.linspace(0,60,1201);h=1+.15*time;sun=-16+.20*time;ridge=4.0
sin=np.sin(np.deg2rad(h));ext=.2/(sin+.025*np.exp(-11*sin))
values=[h-ridge,(-8-sun),h-sun-(2.1+10),3-(1+ext)]
labels=['Star above 4° ridge','Sun below −8°','Altitude difference ≥ 12.1°','Extincted magnitude ≤ 3'];mask=np.stack([v>0 for v in values]).all(axis=0)
fig,ax=plt.subplots(figsize=(7.5,3.1),layout='constrained')
for i,(v,name) in enumerate(zip(values+[mask],labels+['All four criteria'])):
 k=4-i;sel=(v>0);ax.plot([0,60],[k,k],lw=5,c=gray,solid_capstyle='butt');ax.plot(time[sel],[k]*sel.sum(),lw=5,c=teal if i<4 else ochre,solid_capstyle='butt')
ax.set(yticks=list(range(4,-1,-1)),yticklabels=labels+['All four criteria'],xlabel='Minutes along the illustrative dawn track',ylim=(-.6,4.6),xlim=(0,60),title='When can the star be visible?');ax.tick_params(axis='y',length=0);ax.spines['left'].set_visible(False);save(fig,'visibility-window')

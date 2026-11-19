"""Reproduce the illustrative diffuse-light classification. Requires numpy/scipy/Pillow.
Run from the package root: python data-raw/milky-way/rebuild.py
Then Rscript data-raw/milky-way/save-grid.R. Source image is bundled for exact reproduction.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter1d,map_coordinates,label,binary_fill_holes
BASE=Path(__file__).resolve().parent
G=np.array([[-.0548755604,-.8734370902,-.4838350155],[.4941094279,-.4448296300,.7469822445],[-.8676661490,-.1980763734,.4559837762]])
# Fixed image-tone cutoffs preserve the reviewed brightness classes during dark-patch tests.
BRIGHT=72.47510740164273
BRIGHTEST=115.30130722988615

def build(step):
 lon=np.arange(-180,180,step);lat=np.linspace(-60,60,round(120/step)+1);L,B=np.deg2rad(np.meshgrid(lon,lat))
 rgb=np.asarray(Image.open(BASE/'source-map.png').convert('RGB'));h,w=rgb.shape[:2]
 v=np.stack([np.cos(B)*np.cos(L),np.cos(B)*np.sin(L),np.sin(B)],axis=-1)@G
 x=((.5-np.arctan2(v[...,1],v[...,0])/(2*np.pi))*w-.5)%w;y=(.5-np.arcsin(v[...,2])/np.pi)*h-.5
 def sample(a):return map_coordinates(np.pad(a.astype(float),((0,0),(1,1)),mode='wrap'),[y,x+1],order=1,mode='nearest')
 tone=sample(rgb@np.array([.2126,.7152,.0722]))
 def blur(a,fwhm):
  s=fwhm/2.35482/step
  a=np.array([gaussian_filter1d(row,s/np.cos(np.deg2rad(b)),mode='wrap') for row,b in zip(a,lat)])
  return gaussian_filter1d(a,s,axis=0,mode='nearest')
 I=blur(tone,1)
 weights=np.broadcast_to((step*np.rad2deg(np.sin(np.deg2rad(lat+step/2))-np.sin(np.deg2rad(lat-step/2))))[:,None],I.shape)
 def components(mask):
  tags,n=label(mask,np.ones((3,3)));parent=np.arange(n+1)
  def root(a):
   while parent[a]!=a:a=parent[a]
   return a
  for yy in range(len(lat)):
   for dy in (-1,0,1):
    if 0<=yy+dy<len(lat):
     a,b=tags[yy,0],tags[yy+dy,-1]
     if a and b:parent[root(a)]=root(b)
  return np.array([root(a) for a in range(n+1)])[tags]
 def area(mask,cut):
  tags=components(mask);ar=np.bincount(tags.ravel(),weights.ravel());ar[0]=0
  return (tags>0)&(ar[tags]>=cut),ar,tags
 tags=components(blur(I,4)>=40);ids,ct=np.unique(tags[abs(lat)<12],return_counts=True);ids0=ids>0;key=ids[ids0][np.argmax(ct[ids0])]
 outer=binary_fill_holes(np.tile(tags==key,(1,3)))[:,len(lon):2*len(lon)]
 dark,ar,tags=area(outer&(I<.92*blur(I,8)),7)
 bright,_,_=area(outer&~dark&(I>BRIGHT),10)
 top,_,_=area(bright&(I>BRIGHTEST),10)
 classes=outer.astype('u1');classes[dark]=2;classes[bright]=3;classes[top]=4
 yy=np.argmin(abs(lat+.75));xx=np.argmin(abs(lon+78.25))
 metrics=dict(step=step,shape=list(classes.shape),source_sha256=hashlib.sha256((BASE/'source-map.png').read_bytes()).hexdigest(),source_dimensions=[w,h],class_areas_sqdeg=[float(weights[classes==i].sum()) for i in range(5)],target_dark_patch_retained=bool(dark[yy,xx]),target_patch_area_sqdeg=float(ar[tags[yy,xx]]),dark_patch_count=int(np.count_nonzero(ar>=7)),brightness_cutoffs=[BRIGHT,BRIGHTEST])
 return classes,lon,lat,metrics,np.stack([sample(rgb[...,i]) for i in range(3)],axis=-1).clip(0,255).astype('u1')

if __name__=='__main__':
 c,lon,lat,m,pic=build(.1)
 c.tofile(BASE/'grid.bin');(BASE/'metrics.json').write_text(json.dumps(m,indent=2)+'\n')
 # Legacy outline paths follow the extremal occupied cell edges at each longitude.
 import csv
 with (BASE.parents[1]/'inst/extdata/milky-way.csv').open('w') as f:
  wr=csv.writer(f,lineterminator="\n");wr.writerow(['edge','longitude','latitude'])
  for edge in ['south','north']:
   for l in np.arange(0,360.05,.1):
    j=round(((l+180)%360)/.1)%len(lon);rows=np.flatnonzero(c[:,j]>0)
    b=lat[rows[0]]-.05 if edge=='south' else lat[rows[-1]]+.05
    wr.writerow([edge,f'{l:.6f}',f'{b:.6f}'])
 np.savez_compressed(BASE/'review.npz',classes=c,lon=lon,lat=lat,pic=pic)
 print(json.dumps(m,indent=2))

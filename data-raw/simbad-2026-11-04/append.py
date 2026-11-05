"""Rebuild the 2026-11-04 SIMBAD supplement. Run from any directory; offline, idempotent.
SIMBAD database attribution/ODbL and Hipparcos attribution: see inst/COPYRIGHTS.
"""
from pathlib import Path
import csv, json, lzma, math, re, io
here=Path(__file__).resolve().parent
pkg=here.parent.parent
archive=pkg/'inst/ephemeris/sefstars.txt.xz'
def read_json(name):
 j=json.loads((here/name).read_text());return [dict(zip([c['name'] for c in j['metadata']],r)) for r in j['data']]
rows=read_json('simbad.cfa.harvard.edu.json');assert len(rows)==60
flux={r['id']:r for r in read_json('flux.json')}
hip={}
for line in (here/'vizier.cfa.harvard.edu-hip.tsv').read_text().splitlines():
 parts=line.split('\t')
 if re.fullmatch(r'\s*\d+\s*',parts[0]):
  hip[int(parts[0])]=dict(zip(['HIP','Vmag','RAICRS','DEICRS','Plx','pmRA','pmDE'],parts))
assert len(hip)==60
pm=next(r for r in read_json('pm29.json') if r['bibcode']=='2007A&A...474..653V')
plx=next(r for r in read_json('plx29.json') if r['bibcode']=='2007A&A...474..653V')
old=lzma.decompress(archive.read_bytes()).decode('utf-8');lines=old.splitlines(keepends=True)
existing={r[15].strip():i for i,line in enumerate(lines) if not line.startswith('#') and len(r:=next(csv.reader([line])))>=16}
# Swiss ICRS records are J2000 positions, sexagesimal angles; pmRA includes cos(dec).
def sex(value,hours=False):
 if hours:value/=15
 sign='-' if value<0 else '+'
 total=round(abs(value)*3600,8);d=int(total//3600);m=int(total%3600//60);s=total%60
 return [str(d) if hours else sign+str(d),str(m),f'{s:.8f}']
def num(v):
 assert v is not None and math.isfinite(float(v))
 return format(float(v),'.15g')
manifest=[];appended=[]
for r in sorted(rows,key=lambda r:int(r['id'].split()[1])):
 n=int(r['id'].split()[1]);notes=[]
 if n==108797:
  r['pmra'],r['pmdec'],r['pm_bibcode']=pm['pmra'],pm['pmde'],pm['bibcode']
  r['plx_value'],r['plx_bibcode']=plx['plx'],plx['bibcode']
  notes.append('PM/parallax from SIMBAD measurements (Hipparcos 2007), matching coordinate reference; uncertain binary astrometry')
 rv_unknown=r['rvz_radvel'] is None
 if rv_unknown:notes.append('Radial velocity unavailable: zero placeholder, not a measured zero')
 v=flux.get(r['id'])
 vmag=v['flux'] if v else float(hip[n]['Vmag'])
 vref=v['bibcode'] if v else '1997HIP...C......0E; I/239/hip_main Vmag'
 if v is None:notes.append('V magnitude from Hipparcos via VizieR; no SIMBAD V value')
 values=['',r['id'],'ICRS']+sex(r['ra'],True)+sex(r['dec'])+[num(r['pmra']),num(r['pmdec']),num(0 if rv_unknown else r['rvz_radvel']),num(r['plx_value']),num(vmag),(r['sp_type'] or '~').replace(',',';'),r['id']]
 assert len(values)==16 and r['plx_value']>0
 action='appended'
 if n==63948:
  # SIMBAD explicitly lists '*  39 Com A' as an alias of this HIP object.
  assert any(a['id']=='HIP 63948' and re.sub(r'\s+',' ',a['alias']).strip()=='* 39 Com A' for a in read_json('aliases.json'))
  i=existing.get('HIP 63948',existing.get('39 Com A'));assert i is not None
  original=next(csv.reader([lines[i]]));assert original[1]=='39 Com A'
  ra=(float(original[3])+float(original[4])/60+float(original[5])/3600)*15
  dec=float(original[6])+float(original[7])/60+float(original[8])/3600
  assert abs(ra-r['ra'])<1e-8 and abs(dec-r['dec'])<1e-8
  # Only the cross-reference field changes. All existing numeric data retained.
  original[15]='HIP 63948';lines[i]=','.join(original)+'\n';values=original
  action='HIP cross-reference added to existing 39 Com A'
  notes.append('Existing astrometry and magnitude preserved; avoids duplicate star')
 elif r['id'] not in existing:
  appended.append(','.join(values)+'\n')
 else:
  assert next(csv.reader([lines[existing[r['id']]]]))==values,'Existing supplement differs'
 manifest.append(dict(HIP=n,main_id=r['main_id'],identifier=values[1],action=action,ra_deg=r['ra'],dec_deg=r['dec'],pmra_mas_yr=r['pmra'],pmdec_mas_yr=r['pmdec'],parallax_mas=r['plx_value'],radial_velocity_km_s='' if rv_unknown else r['rvz_radvel'],stored_Vmag=values[13],coo_bibcode=r['coo_bibcode'],pm_bibcode=r['pm_bibcode'],plx_bibcode=r['plx_bibcode'],rv_bibcode=r['rvz_bibcode'] or '',V_bibcode=vref if n!=63948 else 'Pre-existing catalogue value retained',retrieved='2026-11-04',notes='; '.join(notes)))
text=''.join(lines)
if appended:
 assert len(appended)==59
 if not text.endswith('\n'):text+='\n'
 text+='# SIMBAD supplement retrieved 2026-11-04; see star-supplement-2026-11-04.csv and ../COPYRIGHTS.\n'
 text+='# Contains SIMBAD information (CDS), ODbL 1.0: https://opendatacommons.org/licenses/odbl/1-0/\n'
 text+=''.join(appended)
archive.write_bytes(lzma.compress(text.encode('utf-8'),preset=9))
out=pkg/'inst/ephemeris/star-supplement-2026-11-04.csv'
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=manifest[0].keys());w.writeheader();w.writerows(manifest)
print('59 appended stars; 1 existing cross-reference. Unknown RV:',sum(not r['rv_bibcode'] for r in manifest),'compressed bytes:',archive.stat().st_size)

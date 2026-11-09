"""Build the single-file GitHub Pages edition from the rendered Quarto book.
Usage (at the Quarto project root): python scripts/standalone/build.py OUTPUT.html
"""
from pathlib import Path
import sys,re,json,base64,gzip,mimetypes
from urllib.parse import urlsplit,unquote
from lxml import html,etree
from copy import deepcopy
here=Path(__file__).resolve().parent
book=here.parents[1];docs=book/'docs'
config=(book/'_quarto.yml').read_text()
pages=[p[:-4] for p in re.findall(r'^\s*- ([\w-]+\.qmd)$',config,re.M)]
anchors=json.loads((here/'anchors.json').read_text())
parser=html.HTMLParser(encoding='utf-8',huge_tree=True)
root=html.fromstring(gzip.decompress((here/'template.html.gz').read_bytes()),parser=parser)
css=etree.SubElement(root.xpath('//head')[0],'style');css.text=(book/'manual.css').read_text()
main=root.xpath('//main')[0];toc=root.get_element_by_id('TOC')
for x in root.xpath('//*[contains(@class,"subtitle")]'):
 x.text=re.search(r'^  subtitle: "(.*)"',config,re.M)[1]
def drop(el):
 if el.getparent() is not None:el.getparent().remove(el)
def elt(tag,text=None,**attrs):
 x=etree.Element(tag,attrib=attrs);x.text=text;return x
reserved=set(root.xpath('//*[@id]/@id'));parts={};mapping={}
help_files=sorted((docs/'help').glob('*.html'))
def help_id(f):return 'help-'+re.sub('[^a-z0-9]+','-',f.stem.lower()).strip('-')
reserved.update(help_id(f) for f in help_files)
reserved.add('individual-function-help')
for m in anchors.values():reserved.add(m['id'])
for page in pages:
 cfg=anchors[page];rid=cfg['id']
 s=html.fromstring((docs/(page+'.html')).read_bytes(),parser=parser).xpath('//main')[0]
 if page=='index':sec=s.xpath('./section')[0];title=sec.find('h1');previous=list(sec)
 else:
  sec=elt('section',id=rid,**{'class':'level1'});title=elt('h1');sec.append(title)
  title.text=(book/(page+'.qmd')).read_text().splitlines()[0].split(' {')[0].lstrip('# ').strip()
  for child in list(s):
   if child.tag!='header':sec.append(child)
 parts[page]=sec
 for node in sec.xpath('.//span[@class="header-section-number"]'):
  parent=node.getparent();parent.text=(parent.text or "")+(node.tail or "");parent.remove(node)
 oldref={n['id']:n for n in cfg['sections']};oldtitle={n['title']:n['id'] for n in cfg['sections']}
 ids={rid:rid}
 for node in sec.xpath('.//*[@id]'):
  k=node.get('id');pick=k
  if node.tag=='section' and k not in oldref:
   hh=next((c for c in node if c.tag in ('h2','h3','h4','h5','h6')),None)
   if hh is not None:pick=oldtitle.get(''.join(hh.itertext()).strip(),k)
  if re.fullmatch(r'cb\d+(?:-\d+)?',k) or pick in reserved:pick=page+'-'+k
  n=2
  while pick in reserved:pick=f'{page}-{k}-{n}';n+=1
  reserved.add(pick);node.set('id',pick);ids[k]=pick
 mapping[page]=ids
# Recreate each individual help page so offline links remain live.
appendix=elt('section',id='individual-function-help',**{'class':'level1 unnumbered'});appendix.append(elt('h1','Individual function help'))
for f in help_files:
 sec=elt('section',id=help_id(f),**{'class':'level2'});sec.append(elt('h2',f.stem+'()' if (docs/'function-methods.html').read_text().find('id="fn-'+help_id(f)[5:]+'"')!=-1 else f.stem))
 h=html.fromstring(f.read_bytes(),parser=parser);body=h.xpath('//body')[0]
 for el in list(body):
  if el.tag=='nav':
   el.tag='p'
  sec.append(el)
 for a in sec.xpath('.//a[starts-with(@href,"../")]'):a.set('href',a.get('href')[3:])
 appendix.append(sec)
def embed(f):
 return 'data:'+(mimetypes.guess_type(f.name)[0] or 'application/octet-stream')+';base64,'+base64.b64encode(f.read_bytes()).decode()
def update(sec,page):
 for el in sec.xpath('.//*[@data-anchor-id]'):
  k=el.get('data-anchor-id');el.set('data-anchor-id',mapping.get(page,{}).get(k,k))
 for el in sec.xpath('.//img[@src]'):
  url=el.get('src')
  if not url.startswith(('data:','http:','https:')):el.set('src',embed(docs/url.lstrip('./')))
 for a in sec.xpath('.//a[@href]'):
  raw=a.get('href');u=urlsplit(raw)
  if u.scheme or u.netloc:continue
  path=unquote(u.path).lstrip('./');fid=unquote(u.fragment)
  if path.startswith('help/'):target=help_id(Path(path))
  elif not path:target=mapping.get(page,{}).get(fid,fid) or anchors.get(page,{}).get('id','individual-function-help')
  elif path.endswith(('.html','.qmd')) and Path(path).stem in anchors and '/' not in path:
   dest=Path(path).stem;target=mapping[dest].get(fid,fid) if fid else anchors[dest]['id']
  elif (docs/path).is_file():a.set('download',Path(path).name);a.set('href',embed(docs/path));continue
  else:raise ValueError((page,raw))
  a.set('href','#'+target)
for page in pages:update(parts[page],page);main.append(parts[page])
update(appendix,'appendix');main.append(appendix)
# Rebuild the navigation in chapter order, using chapter/section/subsection headings.
ul=elt('ul');toc.append(ul)
def nav(parent,sec,depth):
 title=next((x for x in sec if x.tag in ['h1','h2','h3']),None)
 if title is None:return
 li=elt('li');parent.append(li);key=sec.get('id');a=elt('a',' '.join(''.join(title.itertext()).split()),href='#'+key,id='toc-'+key,**{'class':'nav-link','data-scroll-target':'#'+key});li.append(a)
 subs=sec.xpath('./section[@id]') if depth<3 else []
 if subs:
  child=elt('ul',**{'class':'collapse'});li.append(child)
  for x in subs:nav(child,x,depth+1)
for sec in main.findall('section'):nav(ul,sec,1)
ids=root.xpath('//*[@id]/@id')
from collections import Counter
duplicates={k:v for k,v in Counter(ids).items() if v>1}
if duplicates:raise ValueError(('duplicate IDs',duplicates))
bad=[a.get('href') for a in root.xpath('//a[starts-with(@href,"#")]') if a.get('href')[1:] not in set(ids) and a.get('href')!='#']
if bad:raise ValueError(('unresolved links',sorted(set(bad))))
Path(sys.argv[1]).write_bytes(etree.tostring(root,method='html',encoding='utf-8',doctype='<!DOCTYPE html>'))

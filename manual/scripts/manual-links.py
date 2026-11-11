"""Keep published manual links usable after a full Quarto book render."""
from html.parser import HTMLParser
import json
from pathlib import Path
import re

project = Path(__file__).resolve().parents[1]
site = project / 'docs'
scripts = project / 'scripts'

class Anchors(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get('id'):
            self.ids.add(attrs['id'])

pages = {}
for file in site.rglob('*.html'):
    parser = Anchors()
    parser.feed(file.read_text())
    pages[file.relative_to(site).with_suffix('').as_posix()] = parser.ids

legacy = json.loads((scripts / 'legacy-anchors.json').read_text())
moved = json.loads((scripts / 'moved-pages.json').read_text())
routes = {'index': {}}
for key, (name, anchor) in legacy.items():
    if name not in pages:
        raise ValueError(f'Missing destination for {key}: {name}')
    fragment = '#' + anchor if anchor in pages[name] else ''
    if name != 'index' or fragment and fragment != '#' + key:
        routes['index'][key] = name + '.html' + fragment
for key, name in moved.items():
    chapter, anchor = key.split('#')
    routes.setdefault(chapter, {})[anchor] = name + '.html'

for name, ids in pages.items():
    if '/' in name:
        continue
    file = site / (name + '.html')
    text = file.read_text()
    text = re.sub(r'href="(?:\./|/)"', 'href="index.html"', text)
    text = re.sub(r'href="/([\w-]+\.html)', r'href="\1', text)
    text = re.sub(r'href="([\w-]+)\.qmd([#"])', r'href="\1.html\2', text)
    text = re.sub(r'<script id="manual-legacy-links">.*?</script>\s*', '', text, flags=re.S)
    if name in routes:
        mapping = json.dumps(routes[name], ensure_ascii=True)
        redirect = ('<script id="manual-legacy-links">\n'
                    '(() => {const routes = ' + mapping + ';\n'
                    'function follow() {const key = location.hash.slice(1); '
                    'if (Object.hasOwn(routes, key)) location.replace(routes[key]);}\n'
                    "follow(); window.addEventListener('hashchange', follow);})();\n</script>")
        text = text.replace('</head>', redirect + '\n</head>', 1)
    file.write_text(text)

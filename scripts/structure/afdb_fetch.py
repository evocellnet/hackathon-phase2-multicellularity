#!/usr/bin/env python3
"""AlphaFold DB per-protein model fetching (EBI), resumable, polite.

Files: https://alphafold.ebi.ac.uk/files/AF-<ACC>-F1-model_v<N>.pdb  (stored gzipped as received).
The AFDB API (https://alphafold.ebi.ac.uk/api/prediction/<ACC>) gives the current pdbUrl and
globalMetricValue (mean pLDDT).  A User-Agent header is REQUIRED (default curl UA gets 403).
Do not use the GCS bucket (401).  Licence: CC BY 4.0.
"""
import gzip, json, os, sys, time, random, threading, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

UA = 'mcell-phase2/1.0 (spascarelli@ethz.ch)'
BASE = 'https://alphafold.ebi.ac.uk'
AFDB_VERSION = 6  # current as of 2026-10 (v4 files return 404)
MAX_PAR = 16  # connections; throughput is governed by the global limiter below


class _Limiter:
    """Global rate limiter shared by all threads. EBI answers HTTP 429 above ~3.3 req/s sustained per IP
    (measured 2026-10-28: 3.0/s sustained for 400 s OK; 4/s -> 429 after ~330 requests; block lasts
    ~1-2 min, no Retry-After).  On 429: pause everybody 90 s and lower the rate by 15% (once per event)."""
    def __init__(self, rate=3.0, rmin=1.0):
        self.rate, self.rmin = rate, rmin
        self.lock = threading.Lock(); self.next = 0.0; self.pause_until = 0.0
    def wait(self):
        with self.lock:
            now = time.time()
            t = max(self.next, self.pause_until, now)
            self.next = t + 1.0 / self.rate
        if t > now: time.sleep(t - now)
    def hit429(self):
        with self.lock:
            now = time.time()
            if now >= self.pause_until:        # first thread to notice this block event
                self.rate = max(self.rmin, self.rate * 0.85)
                print(f'[afdb] 429 -> pause 90s, rate now {self.rate:.2f}/s', file=sys.stderr, flush=True)
            self.pause_until = max(self.pause_until, now + 90)
    def success(self):
        pass
LIM = _Limiter(rate=float(os.environ.get('AFDB_RATE', 3.0)))


def _get(url, tries=6, timeout=60):
    """GET with retries/backoff + global rate limit. 404 returned immediately.
    Returns (status, bytes, content_encoding); status -1 on persistent failure."""
    last = None
    n429 = 0
    i = 0
    while i < tries and n429 < 30:
        LIM.wait()
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA, 'Accept-Encoding': 'gzip'})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                LIM.success()
                return r.status, r.read(), r.headers.get('Content-Encoding', '')
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                return e.code, b'', ''
            last = e
            if e.code == 429:       # rate limited: global pause, does not consume a retry
                LIM.hit429(); n429 += 1; continue
        except Exception as e:  # timeouts, resets
            last = e
        i += 1
        time.sleep(min(60, 2 ** i + random.random()))
    return -1, str(last).encode(), ''


def api_entry(acc):
    st, body, enc = _get(f'{BASE}/api/prediction/{acc}')
    if st != 200:
        return None
    if enc == 'gzip' or body[:2] == b'\x1f\x8b':
        body = gzip.decompress(body)
    j = json.loads(body)
    return j[0] if j else None


def mean_plddt_from_pdb_gz(path):
    """Mean pLDDT = mean B-factor over CA atoms of the (gz) PDB."""
    s = n = 0
    with gzip.open(path, 'rt') as f:
        for l in f:
            if l.startswith('ATOM') and l[12:16] == ' CA ':
                s += float(l[60:66]); n += 1
    return (s / n) if n else float('nan')


def fetch_pdb(acc, outdir, version=AFDB_VERSION):
    """Download AF-<acc>-F1 model to outdir/<name>.pdb.gz. Resumable (skips valid existing file).
    Returns (status, path|None, version_used). status in {'ok','cached','missing','error'}."""
    os.makedirs(outdir, exist_ok=True)
    final = os.path.join(outdir, f'AF-{acc}-F1-model_v{version}.pdb.gz')
    if os.path.exists(final) and os.path.getsize(final) > 100:
        return 'cached', final, version
    st, body, enc = _get(f'{BASE}/files/AF-{acc}-F1-model_v{version}.pdb')
    ver = version
    if st == 404:  # fall back to API-advertised URL (newer/other version)
        e = api_entry(acc)
        if e and e.get('pdbUrl'):
            url = e['pdbUrl']
            st, body, enc = _get(url)
            ver = e.get('latestVersion', version)
            final = os.path.join(outdir, os.path.basename(url) + '.gz')
        else:
            return 'missing', None, None
    if st != 200:
        return ('missing' if st == 404 else 'error'), None, None
    data = body if (body[:2] == b'\x1f\x8b') else gzip.compress(body, 6)
    tmp = final + '.part'
    with open(tmp, 'wb') as f:
        f.write(data)
    try:  # integrity check
        gzip.decompress(data)
    except Exception:
        os.remove(tmp); return 'error', None, None
    os.replace(tmp, final)
    return 'ok', final, ver


def fetch_many(accs, outdir, workers=12, log=None, on_progress=None):
    """Parallel fetch; returns dict acc -> (status, path, version). workers capped at 16."""
    workers = min(workers, MAX_PAR)
    res = {}
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(fetch_pdb, a, outdir): a for a in accs}
        for k, fu in enumerate(as_completed(futs), 1):
            a = futs[fu]
            try:
                res[a] = fu.result()
            except Exception as e:
                res[a] = ('error', None, None)
            if on_progress and k % 1000 == 0:
                on_progress(k, len(accs))
    return res

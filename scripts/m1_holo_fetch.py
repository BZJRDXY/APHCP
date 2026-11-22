"""Fetch bounded RCSB metadata for the APHCP M1 inventory (no coordinates)."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from datetime import datetime, timedelta, timezone
import hashlib
from http.client import IncompleteRead
import json
import os
from pathlib import Path
import sys
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from m1_holo_config import date_filters, read_config, read_json, sha256, window

RUN = Path(sys.argv[1])
SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query'
DATA = 'https://data.rcsb.org/graphql'
HOLDINGS_STATUS = 'https://data.rcsb.org/rest/v1/holdings/status/'


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(payload, ensure_ascii=False, separators=(',', ':'))+'\n', encoding='utf-8')
    tmp.replace(path)


def request(path, url, body, expected, *, reusable=None, refresh=False, keep_partial=None, split_incomplete=False):
    previous = json.loads(path.read_text(encoding='utf8')) if path.is_file() else {}
    if (not refresh and previous.get('request') == body
            and (reusable(previous) if reusable else previous.get('usable'))):
        return previous
    tries = previous.get('attempts', []) if previous.get('request') == body else []
    envelope = {'url':url,'method':'POST','request':body,'attempts':tries}
    for attempt in range(3):
        item = {'started_utc':now()}; started = time.monotonic()
        try:
            payload = json.dumps(body).encode()
            with urlopen(Request(url, payload, {'Content-Type':'application/json','Accept':'application/json'}), timeout=38) as response:
                code, raw = response.status, response.read()
            item.update(status=code, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            data = json.loads(raw) if raw else None
            usable = code in expected and not (isinstance(data,dict) and data.get('errors'))
            envelope.update(response=data, usable=usable)
        except Exception as error:
            item['error'] = repr(error)
            if isinstance(error, IncompleteRead):
                item.update(incomplete_read=True, partial_bytes=len(error.partial))
            if isinstance(error,HTTPError):
                item.update(status=error.code,error_body=error.read(8000).decode(errors='replace'))
            envelope.update(usable=False)
        item.update(finished_utc=now(), elapsed_seconds=round(time.monotonic()-started,3))
        tries.append(item); save(path,envelope)
        if (envelope['usable'] or (keep_partial and keep_partial(envelope))
                or item.get('status') in (400,401,403,404,422)
                or (split_incomplete and item.get('incomplete_read'))):
            break
        time.sleep(1+attempt)
    return envelope


def page_values(result, offset):
    if not result.get('usable'):
        return None, []
    status = result['attempts'][-1].get('status')
    if offset == 0 and status == 204:
        return 0, []
    data = result.get('response')
    if not isinstance(data, dict) or type(data.get('total_count')) is not int or data['total_count'] < 0:
        return None, []
    hits = data.get('result_set', [])
    if not isinstance(hits, list) or any(not isinstance(h, dict) or not isinstance(h.get('identifier'), str) or not h['identifier'] for h in hits):
        return None, []
    return data['total_count'], [h['identifier'] for h in hits]


def search():
    started = now()
    cfg = read_config(RUN)
    start, end, cutoff = window(cfg)
    body = cfg['m1_1_request']; width = cfg['page_size']
    previous = read_json(RUN/'search'/'summary.json') or {}
    refresh = set(previous.get('refresh_offsets', []))
    if not previous.get('concordant'):
        refresh.update(x['start'] for x in previous.get('page_checks', [])
                       if not x['usable'] or x['rows'] != x['expected_rows'])

    def page(offset, baseline=None):
        q=deepcopy(body); q['request_options']['paginate']={'start':offset,'rows':width}
        def valid_cached(result):
            total, ids = page_values(result, offset)
            return (total is not None and (baseline is None or total == baseline)
                    and len(ids) == len(set(ids)) == min(width, max(0, total-offset)))
        return offset, request(RUN/'search'/f'page_{offset:04}.json', SEARCH, q,
                               {200, 204} if offset == 0 else {200},
                               reusable=valid_cached, refresh=offset in refresh)

    def boundary(day):
        q = deepcopy(body)
        filters = date_filters(q)
        filters['greater_or_equal']['value'] = day.isoformat()
        filters['less']['value'] = (day + timedelta(days=1)).isoformat()
        q['request_options']={'results_content_type':['experimental'],'return_counts':True}
        return day.isoformat(),request(RUN/'search'/f'count_{day}.json',SEARCH,q,{200,204},
                                        reusable=lambda result:page_values(result,0)[0] is not None)

    _, first = page(0)
    total, _ = page_values(first, 0)
    if total is None:
        save(RUN/'search'/'summary.json', {'started_utc':started,'concordant':False,'refresh_offsets':[0], 'entry_ids':[]})
        raise SystemExit('Search page 0 is unresolved; see request receipt')
    boundary_days = {start, end}
    if cutoff and start <= cutoff < end:
        boundary_days.add(cutoff)
    with ThreadPoolExecutor(max_workers=3) as pool:
        pages={0:first,**dict(pool.map(lambda offset:page(offset,total), range(width,total,width)))}
        bounds = dict(pool.map(boundary, sorted(boundary_days)))
    page_checks=[]; found={}; refresh=set(); differing_totals=set()
    for offset,r in sorted(pages.items()):
        reported_total, ids = page_values(r, offset)
        expected_rows = min(width,total-offset)
        page_checks.append({'start':offset,'expected_rows':expected_rows,'rows':len(ids),'unique':len(set(ids)),'reported_total':reported_total,'usable':r['usable']})
        if reported_total != total or len(ids) != len(set(ids)) or len(ids) != expected_rows:
            refresh.add(offset)
        if reported_total is not None and reported_total != total:
            differing_totals.add(reported_total)
        for pdb in ids:
            if pdb in found:
                refresh.update((offset, found[pdb]))
            else:
                found[pdb] = offset
    # A live total mismatch can also invalidate the first page's baseline.
    if differing_totals:
        refresh.add(0)
    valid = not refresh and len(found) == total
    bound_checks={day:{'usable':page_values(r,0)[0] is not None,'http_status':r['attempts'][-1].get('status'),
                       'count':page_values(r,0)[0]} for day,r in bounds.items()}
    out={'started_utc':started,'finished_utc':now(),'total_count':total,'rows':sum(p['rows'] for p in page_checks),'unique_ids':len(found),'page_checks':page_checks,'boundary_counts':bound_checks,'concordant':valid,'refresh_offsets':sorted(refresh),'entry_ids':list(found)}
    save(RUN/'search'/'summary.json',out)
    print(json.dumps({k:v for k,v in out.items() if k!='entry_ids'}),flush=True)
    if not valid:
        raise SystemExit('Search has unresolved pages; rerun to request only suspect offsets')

GQL = '''query M1Month($ids: [String!]!) {
 entries(entry_ids: $ids) {
  rcsb_id struct { title } exptl { method }
  rcsb_entry_container_identifiers { assembly_ids model_ids polymer_entity_ids non_polymer_entity_ids branched_entity_ids }
  rcsb_accession_info { initial_release_date deposit_date revision_date major_revision minor_revision status_code }
  pdbx_audit_revision_history { data_content_type major_revision minor_revision revision_date ordinal }
  rcsb_entry_info { structure_determination_methodology polymer_entity_count_protein deposited_atom_count deposited_model_count resolution_combined inter_mol_covalent_bond_count }
  rcsb_primary_citation { year pdbx_database_id_DOI pdbx_database_id_PubMed }
  pdbx_molecule_features { prd_id type name class details }
  polymer_entities {
   rcsb_polymer_entity { pdbx_description pdbx_ec }
   entity_poly { type rcsb_entity_polymer_type rcsb_sample_sequence_length nstd_monomer }
   rcsb_polymer_entity_container_identifiers { entity_id asym_ids auth_asym_ids uniprot_ids prd_id }
   rcsb_polymer_entity_annotation { type annotation_id name provenance_source }
  }
  nonpolymer_entities {
   pdbx_entity_nonpoly { comp_id entity_id name }
   rcsb_nonpolymer_entity_container_identifiers { entity_id asym_ids auth_asym_ids prd_id }
  }
  branched_entities {
   pdbx_entity_branch { type rcsb_branched_component_count }
   rcsb_branched_entity_container_identifiers { entity_id asym_ids auth_asym_ids prd_id chem_comp_monomers }
  }
  assemblies {
   rcsb_assembly_container_identifiers { assembly_id }
   pdbx_struct_assembly { oligomeric_count oligomeric_details }
   pdbx_struct_assembly_gen { assembly_id asym_id_list oper_expression }
   rcsb_assembly_info { atom_count polymer_entity_instance_count_protein nonpolymer_entity_instance_count branched_entity_instance_count }
  }
 }
}'''


def graphql_rows(result, ids):
    """Keep only unambiguous requested rows, even in a partial GraphQL response."""
    if (not result.get('attempts') or result['attempts'][-1].get('status') != 200):
        return {}
    body = result.get('response')
    if not isinstance(body, dict) or not isinstance(body.get('data'), dict):
        return {}
    entries = body['data'].get('entries')
    if not isinstance(entries, list):
        return {}
    errored=set(); found={}; duplicates=set(); wanted=set(ids)
    for error in body.get('errors') or []:
        path = error.get('path') if isinstance(error, dict) else None
        if not isinstance(path, list) or len(path)<2 or path[0]!='entries' or type(path[1]) is not int:
            return {}
        errored.add(path[1])
    for i, row in enumerate(entries):
        pdb = row.get('rcsb_id') if isinstance(row, dict) else None
        if i in errored or not isinstance(pdb,str) or pdb not in wanted:
            continue
        if pdb in found:
            duplicates.add(pdb)
        found[pdb] = row
    return {pdb:row for pdb,row in found.items() if pdb not in duplicates}


def metadata_workers():
    # Late in a batch, a single long window can use capacity freed by its peers.
    workers = int(os.environ.get('APHCP_M1_METADATA_WORKERS', '4'))
    if not 1 <= workers <= 24:
        raise ValueError('APHCP_M1_METADATA_WORKERS must be between 1 and 24')
    return workers


def metadata(probe=False):
    ids = (read_json(RUN/'search'/'summary.json') or {})['entry_ids']
    started=now(); outdir=RUN/'data'; entries=outdir/'entries.json'
    prior=read_json(outdir/'summary.json') or {}; wanted=set(ids)
    evidence=read_json(RUN/'inputs'/'confirmed_unavailable.json') or {}
    if (not isinstance(evidence, dict) or set(evidence)-wanted
            or any(not isinstance(v,dict) or not all(isinstance(v.get(k),str) and v[k].strip()
                       for k in ('source','reference','reason')) for v in evidence.values())):
        raise ValueError('Unavailable evidence requires a Search ID and source/reference/reason per entry')

    def chunk(name, group, *, bulk=False):
        path=outdir/f'{name}.json'; body={'query':GQL,'variables':{'ids':group}}
        old=read_json(path) or {}
        # The bulk response is a snapshot. Repair its unresolved rows via per-ID requests.
        if bulk and old.get('request')==body and old.get('attempts'):
            result=old
        else:
            result=request(path,DATA,body,{200},reusable=lambda r:set(graphql_rows(r,group))==set(group),
                           keep_partial=(lambda r:graphql_rows(r,group)) if bulk else None, split_incomplete=bulk and len(group)>1)
        return result, graphql_rows(result,group)

    if probe:
        _, found=chunk('probe',['5R42','5SAK','1HEW'])
        print(json.dumps({'received':list(found)}))
        if len(found)!=3: raise SystemExit(1)
        return
    valid_snapshot=entries.is_file() and (not prior.get('entries_sha256') or sha256(entries)==prior['entries_sha256'])
    known={r['rcsb_id']:r for r in ((read_json(entries) or []) if valid_snapshot else [])
           if r.get('rcsb_id') in wanted and not r.get('_missing') and not r.get('_confirmed_unavailable')}
    bad=list(prior.get('unusable_batch_offsets', [])); pending=[p for p in ids if p not in known and p not in evidence]
    with ThreadPoolExecutor(max_workers=metadata_workers()) as pool:
        if not valid_snapshot:
            groups=[(i,ids[i:i+20]) for i in range(0,len(ids),20)]
            def bulk(task):
                ix, group=task
                group=[p for p in group if p not in evidence]
                return ix, chunk(f'page_{ix:04}',group,bulk=True) if group else (None,{})
            for ix,(result, found) in pool.map(bulk,groups):
                known.update(found)
                if result is not None and not result['usable']: bad.append(ix)
            pending=[p for p in ids if p not in known and p not in evidence]
        for _,found in pool.map(lambda p:chunk(f'retry_{p}',[p]),pending):
            known.update(found)
    unavailable={p:evidence[p] for p in ids if p in evidence and p not in known}
    missing=[p for p in ids if p not in known and p not in unavailable]
    rows=[known[p] if p in known else {'rcsb_id':p,'_confirmed_unavailable':unavailable[p]}
          if p in unavailable else {'rcsb_id':p,'_missing':True} for p in ids]
    save(entries,rows)
    out={'started_utc':started,'finished_utc':now(),'requested':len(ids),'received':len(known),'unusable_batch_offsets':bad,
         'missing_ids':missing,'confirmed_unavailable':unavailable,'entries_sha256':sha256(entries)}
    save(outdir/'summary.json',out)
    print(json.dumps(out),flush=True)


def confirmed_obsolete(receipt, pdb):
    if receipt.get('url') != HOLDINGS_STATUS+pdb or receipt.get('method') != 'GET':
        return False
    if (receipt.get('attempts') or [{}])[-1].get('status') != 200:
        return False
    body = receipt.get('response')
    if not isinstance(body, dict):
        return False
    keys = body.get('rcsb_repository_holdings_combined_entry_container_identifiers') or {}
    status = body.get('rcsb_repository_holdings_combined') or {}
    return (isinstance(keys,dict) and isinstance(status,dict)
            and body.get('rcsb_id') == keys.get('entry_id') == pdb
            and status.get('status') == 'REMOVED' and status.get('status_code') == 'OBS')


def obsolete_holding(inventory, pdb):
    # A GraphQL omission or core REST 404 alone cannot prove archive withdrawal.
    if not pdb.isascii() or not pdb.isalnum():
        raise ValueError(f'Invalid PLINDER-only PDB identifier: {pdb!r}')
    path = inventory/f'period_holdings_{pdb}.json'; url = HOLDINGS_STATUS+pdb
    old = read_json(path) or {}
    if confirmed_obsolete(old, pdb):
        return path.name
    attempts = old.get('attempts', []) if old.get('url') == url else []
    item = {'started_utc':now()}; started = time.monotonic(); data = None
    try:
        with urlopen(Request(url,headers={'Accept':'application/json'}),timeout=15) as r:
            raw = r.read(); item.update(status=r.status, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            data = json.loads(raw)
    except (OSError, ValueError) as error:
        item['error'] = repr(error)
        if isinstance(error, HTTPError):
            item.update(status=error.code, error_body=error.read(4000).decode(errors='replace'))
    item.update(finished_utc=now(),elapsed_seconds=round(time.monotonic()-started,3))
    receipt = {'url':url,'method':'GET','response':data,'attempts':[*attempts,item]}
    save(path, receipt)
    return path.name if confirmed_obsolete(receipt, pdb) else None


def resolve_period_obsolete(rows, inventory):
    pending = [(row,row['entry_pdb_id'].upper()) for row in rows if row.get('_technical_missing')]
    resolved = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        receipts = pool.map(lambda pair:obsolete_holding(inventory,pair[1]),pending)
        for (row,pdb),evidence in zip(pending,receipts):
            if not evidence:
                continue
            row.pop('_technical_missing')
            row.update(rcsb_protein_eligibility='CONFIRMED_OBSOLETE_RCSB',rcsb_status_evidence=evidence)
            resolved.append(pdb)
    return resolved


def repair_period_holdings():
    """Repair only an already persisted PLINDER-only result, without refetching Search or GraphQL."""
    inventory = RUN/'inventory'; path = inventory/'plinder_period_only_reconcile.json'
    original = read_json(path)
    if not isinstance(original,list):
        raise ValueError(f'No prior PLINDER period result: {path}')
    rows = deepcopy(original)
    resolved = resolve_period_obsolete(rows,inventory)
    if resolved:
        prior = inventory/'plinder_period_only_reconcile_before_holdings.json'
        if not prior.is_file():
            save(prior,original)
        save(path,rows)
    print(json.dumps({'run':str(RUN),'confirmed_obsolete':resolved,
                      'unresolved':[r['entry_pdb_id'].upper() for r in rows if r.get('_technical_missing')]}),flush=True)


def period():
    import csv
    path=RUN/'inventory'/'plinder_period_only.csv'
    with path.open(encoding='utf8') as f: rows=list(csv.DictReader(f))
    ids=[r['entry_pdb_id'].upper() for r in rows]
    if not ids:
        save(RUN/'inventory'/'plinder_period_only_reconcile.json', [])
        return
    start, end, _ = window(read_config(RUN))
    start_day, end_day = start.isoformat(), end.isoformat()
    q='''query PeriodReconcile($ids: [String!]!) {entries(entry_ids:$ids) {rcsb_id exptl {method} rcsb_accession_info {initial_release_date} rcsb_entry_info {structure_determination_methodology polymer_entity_count_protein} }}'''
    inv=RUN/'inventory'; raw=inv/'plinder_period_only_data_api.json'; body={'query':q,'variables':{'ids':ids}}
    cached=read_json(raw) or {}
    r=cached if cached.get('request')==body else request(raw,DATA,body,{200})
    found=graphql_rows(r,ids)
    def retry(pdb):
        one=request(inv/f'period_retry_{pdb}.json',DATA,{'query':q,'variables':{'ids':[pdb]}},{200},
                    reusable=lambda result:pdb in graphql_rows(result,[pdb]))
        return graphql_rows(one,[pdb])
    with ThreadPoolExecutor(max_workers=4) as pool:
        for got in pool.map(retry, [p for p in ids if p not in found]):
            found.update(got)
    for row,pdb in zip(rows,ids):
        e=found.get(pdb); a=(e or {}).get('rcsb_accession_info') or {}; i=(e or {}).get('rcsb_entry_info') or {}; date=a.get('initial_release_date'); protein=i.get('polymer_entity_count_protein'); methodology=i.get('structure_determination_methodology')
        if e is None:
            status = 'REVIEW_DATA_API_MISSING_OR_INCOMPLETE'
            row['_technical_missing']=True
        elif not date or not methodology:
            status = 'REVIEW_DATA_API_MISSING_OR_INCOMPLETE'
        elif not start_day <= date[:10] < end_day:
            status = 'OUTSIDE_WINDOW'
        elif methodology != 'experimental':
            status = 'NOT_EXPERIMENTAL'
        elif protein == 0:
            status = 'NO_PROTEIN_ENTITY'
        elif protein is None:
            status = 'REVIEW_PROTEIN_FIELD_MISSING'
        else:
            status = 'REVIEW_SEARCH_DISAGREEMENT'
        row.update(rcsb_protein_eligibility=status,rcsb_release_date=date,rcsb_protein_entity_count=protein,rcsb_methodology=methodology)
    resolve_period_obsolete(rows,inv)
    save(RUN/'inventory'/'plinder_period_only_reconcile.json', rows)
    from collections import Counter
    print(dict(Counter(r['rcsb_protein_eligibility'] for r in rows)),flush=True)


if __name__ == '__main__':
    {'search':search,'probe':lambda:metadata(True),'metadata':metadata,'period':period,'repair-period-holdings':repair_period_holdings}[sys.argv[2]]()

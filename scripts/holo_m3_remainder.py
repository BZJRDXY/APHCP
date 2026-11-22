"""Freeze the exact M2 shortlist difference without excluding other systems at the same PDB."""
import argparse
from pathlib import Path

from holo_m3_register import jsonl, read, require, sha, write


def key(row):
    return str(row['batch']),row['pdb_id'],row['candidate_key']


def keyed(rows, get_key=key):
    out={}
    for row in rows:
        k=get_key(row)
        require(k not in out,'duplicate_pool_key')
        out[k]=row
    return out


def freeze(ref, first, output, expected, expected_entries):
    shortlist=ref/'m3_program_shortlist.json'; history=ref/'historical_deltas.json'; attribution=ref/'system_attribution.jsonl'
    base=keyed(read(shortlist)['candidates']); examined=keyed(read(first))
    delta=keyed((d for d in read(history)['deltas'] if d.get('delta_intake',0)>0),lambda d:(d['pdb_id'],d['candidate_key']))
    matches=[]
    for row in jsonl(attribution):
        k=key(row); d=delta.get(k[1:]); prior=row.get('historical_retest') or {}
        if d and prior and all(d.get(f)==v for f,v in prior.items()):
            matches.append(row|{'system_sha256':d['system_sha256']})
    hist=keyed(matches)
    require(len(hist)==len(delta) and not base.keys()&hist.keys(),'history_missing_or_already_baseline')
    eligible=base|hist; require(examined.keys()<=eligible.keys(),'first_contains_unknown_key')
    remaining=eligible.keys()-examined.keys(); n_entries=len({k[1] for k in remaining})
    require(len(remaining)==expected and n_entries==expected_entries,'unexpected_remainder_denominator')
    def pack(row):
        return {k:row[k] for k in ('batch','pdb_id','candidate_key','revision','assembly','source_sha256','system_sha256')}|{'source_model':row['model']}
    result={'review_count':len(remaining),'entry_count':n_entries,'pool':{'baseline':len(base),'historical_delta':len(hist),'excluded_all_first':len(examined)},
            'source_files':{str(p):sha(p) for p in (shortlist,history,attribution,first)},
            'review_candidates':[pack(eligible[i]) for i in sorted(remaining)]}
    if output.exists():require(read(output)==result,'frozen_remainder_changed')
    else:write(output,result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for field in ('reconciliation','prior-results','output'):p.add_argument('--'+field,type=Path,required=True)
    p.add_argument('--count',type=int,required=True);p.add_argument('--entries',type=int,required=True)
    a=p.parse_args(); result=freeze(a.reconciliation,a.prior_results,a.output,a.count,a.entries)
    print(result['pool'],result['review_count'],result['entry_count'])

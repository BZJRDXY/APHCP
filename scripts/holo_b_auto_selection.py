"""Freeze a held-out Holo B pilot input without using human outcome labels."""
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import random
import sqlite3

import holo_b_auto_research as engine


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',type=Path,required=True)
    p.add_argument('--sites',type=Path,required=True)
    p.add_argument('--previous-sample',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=int,default=20261126)
    args=p.parse_args()
    db=engine.read_db(args.database)
    db.execute('attach database ? as site',(args.sites.resolve().as_uri()+'?mode=ro',))
    prior=json.loads(args.previous_sample.read_text(encoding='utf-8'))
    prior_pdb={'1A6A','1MBN','5MYQ','8RNP','5MO0'}
    prior_pdb.update(x['pdb_id'] for d in prior['strata'].values() for x in d['selected'])
    metadata_cache={}
    def metadata(sid):
        row=engine.load_record(db,sid);obs=json.loads(row['document'])
        path=obs['files']['source_mmcif']['path']
        if path not in metadata_cache:
            b=engine.read_cif(path)
            metadata_cache[path]=(engine.protein_keys(b),engine.publication_keys(engine.publication(b)))
        return obs,metadata_cache[path]
    used_protein=set();used_paper=set()
    marks=','.join('?' for _ in prior_pdb)
    prior_rows=db.execute('select pdb_id,min(system_id) from result where pdb_id in ('+marks+') group by pdb_id',
                          tuple(sorted(prior_pdb))).fetchall()
    for _,sid in prior_rows:
        _,(proteins,papers)=metadata(sid);used_protein.update(proteins);used_paper.update(papers)
    used_pdb=set(prior_pdb)
    groups=defaultdict(lambda:defaultdict(list))
    why='["target_role_evidence_insufficient_or_conflicting"]'
    sql='''select r.system_id,r.pdb_id,r.disposition,r.reasons,r.primary_class,r.match_level,
                  r.annotation_interpretation,r.conflicts,r.annotations,s.model_evidence,s.pocket_agreement
           from result r left join site.site_check s using(system_id)
           where (r.disposition='REVIEW' and r.reasons=?) or r.disposition='PREPARATION_USE_RECORDED'
           order by r.system_id'''
    for row in db.execute(sql,(why,)):
        if row['pdb_id'] in prior_pdb:continue
        if row['disposition']=='PREPARATION_USE_RECORDED':
            groups['preparation_recorded'][row['pdb_id']].append(row['system_id']);continue
        if row['match_level'] in ('NO_ANNOTATION','NO_CCD_MATCH'):
            label={'SMALL_MOLECULE':'no_plinder_small','ION_OR_ATOMIC_COMPONENT':'no_plinder_ion'}.get(row['primary_class'])
            if label:groups[label][row['pdb_id']].append(row['system_id'])
        if (row['match_level']!='SAME_REVISION_INSTANCE_CHAIN_MATCH' or
            row['annotation_interpretation']!='PLINDER_PROPER_RULE_SUPPORT' or
            row['model_evidence']!='SOURCE_SINGLE_MODEL_1' or row['pocket_agreement']!='EXACT_NONEMPTY' or
            json.loads(row['conflicts'])!=[]):continue
        annotations=json.loads(row['annotations'])
        codes={a.get('ligand_ccd_code') for a in annotations}
        cov={a.get('ligand_is_covalent') for a in annotations}
        if len(codes)!=1 or len(cov)!=1:continue
        code=next(iter(codes));covalent=next(iter(cov))
        label='proper_nag_covalent' if code=='NAG' and covalent is True else \
              'proper_other_covalent' if code!='NAG' and covalent is True else \
              'proper_non_covalent' if covalent is False and row['primary_class']=='SMALL_MOLECULE' else None
        if label:groups[label][row['pdb_id']].append(row['system_id'])
    quota={'proper_nag_covalent':10,'proper_other_covalent':10,'proper_non_covalent':10,
           'no_plinder_small':5,'no_plinder_ion':5,'preparation_recorded':8}
    rng=random.Random(args.seed);selected=[];skipped=Counter();eligible={}
    for group,n in quota.items():
        pdbs=sorted(groups[group]);eligible[group]=len(pdbs);rng.shuffle(pdbs)
        picked=0
        for pdb in pdbs:
            if pdb in used_pdb:skipped[group+':pdb_overlap']+=1;continue
            choices=sorted(groups[group][pdb]);sid=choices[rng.randrange(len(choices))]
            try:obs,(proteins,papers)=metadata(sid)
            except (OSError,ValueError,RuntimeError) as exc:
                skipped[group+':selection_source_error']+=1;continue
            if not proteins:skipped[group+':protein_identity_missing']+=1;continue
            if proteins & used_protein:skipped[group+':protein_overlap']+=1;continue
            if papers & used_paper:skipped[group+':publication_overlap']+=1;continue
            used_pdb.add(pdb);used_protein.update(proteins);used_paper.update(papers)
            selected.append({'stratum':group,'system_id':sid,'pdb_id':pdb,'revision':obs['identity']['revision'],
                'source_sha256':obs['files']['source_mmcif']['sha256'],
                'protein_keys':[list(k) for k in sorted(proteins)],
                'publication_keys':[list(k) for k in sorted(papers)],
                'publication_independence':'checked_keys' if papers else 'unverified_no_keys'})
            picked+=1
            if picked==n:break
        if picked!=n:raise ValueError(f'Insufficient distinct candidates in {group}: {picked}/{n}')
    result={'created_utc':engine.now(),'selection_seed':args.seed,'database':str(args.database.resolve()),
            'site_database':str(args.sites.resolve()),'previous_sample':str(args.previous_sample.resolve()),
            'code_sha256':engine.sha(engine.__file__),'selector_sha256':engine.sha(__file__),
            'exclude_prior_pdb_count':len(prior_pdb),'protein_dedup':'exact canonical sequence or UniProt accession',
            'publication_dedup':'matching DOI, PMID or normalized primary-title text',
            'stratum_quotas':quota,'eligible_pdb_before_dedup':eligible,'selection_exclusions':dict(skipped),
            'systems':selected}
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({'selected':len(selected),'strata':quota,'prior_pdb_excluded':len(prior_pdb),
                      'selection_exclusions':dict(skipped)},ensure_ascii=False),flush=True)


if __name__=='__main__':main()

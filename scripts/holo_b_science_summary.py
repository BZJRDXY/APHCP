"""Produce a compact, reproducible summary of the Holo B research database."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3

from holo_b_science_sampling import eligible_candidate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    c = sqlite3.connect(a.database.resolve().as_uri()+'?mode=ro', uri=True)
    c.row_factory = sqlite3.Row
    if c.execute('select count(*) from result').fetchone()[0] != c.execute('select count(*) from observation').fetchone()[0]:
        raise ValueError('Research derivation is incomplete')
    only = "disposition='REVIEW' and reasons='[\"target_role_evidence_insufficient_or_conflicting\"]'"
    def rows(sql):
        return [dict(x) for x in c.execute(sql)]
    def counts(where):
        return dict(c.execute('select count(*) systems,count(distinct instance_key) instances,'
                              'count(distinct evidence_key) evidence_groups from result where '+where).fetchone())
    out = {'all':counts('1'), 'role_only':counts(only)}
    out['disposition'] = rows('select disposition,count(*) systems from result group by disposition order by systems desc')
    out['decision_state'] = rows('select decision_state,count(*) systems from result group by decision_state order by systems desc')
    out['role_only_prior_current_conflicts'] = rows(f'''select coalesce(prior_role,'NOT_RECORDED') prior_role,
        coalesce(current_role,'NOT_RECORDED') current_role,
        case when conflicts='null' then 'NOT_RECORDED' when json_array_length(conflicts)>0 then 'CONFLICT_RECORDED' else 'NO_RECORDED_CONFLICT' end conflict_state,
        count(*) systems from result where {only} group by 1,2,3 order by systems desc''')
    out['role_only_flags'] = {name:counts(f'({only}) and ({condition})') for name,condition in {
        'old_not_pass':"prior_role is not null and prior_role!='PASS'",
        'old_pass':"prior_role='PASS'",
        'current_not_pass':"current_role is not null and current_role!='PASS'",
        'current_pass':"current_role='PASS'",
        'old_not_pass_current_pass':"prior_role is not null and prior_role!='PASS' and current_role='PASS'",
        'old_pass_current_not_pass':"prior_role='PASS' and current_role is not null and current_role!='PASS'",
        'recorded_conflict':"conflicts!='null' and json_array_length(conflicts)>0",
        'multiple_partner_members':"partner_count>1",
        'multiple_target_residues':"target_residue_count>1",
        'multi_member_and_multi_residue':"partner_count>1 and target_residue_count>1",
        'multi_member_and_competing':"partner_count>1 and json_extract(source_evidence,'$.reason') like '%competing%'",
        'multi_residue_and_competing':"target_residue_count>1 and json_extract(source_evidence,'$.reason') like '%competing%'",
        'recorded_conflict_and_competing':"conflicts!='null' and json_array_length(conflicts)>0 and json_extract(source_evidence,'$.reason') like '%competing%'",
        'old_four_qc_pass':"json_extract(old_qc,'$.coordinate_identity')='PASS' and json_extract(old_qc,'$.experimental_fidelity')='PASS' and json_extract(old_qc,'$.topology')='PASS' and json_extract(old_qc,'$.partner_scope_qc')='PASS'",
        'new_four_qc_pass':"json_extract(new_qc,'$.coordinate_identity')='PASS' and json_extract(new_qc,'$.experimental_fidelity')='PASS' and json_extract(new_qc,'$.topology')='PASS' and json_extract(new_qc,'$.partner_scope')='PASS'",
        'source_reason_no_unambiguous_statement':"json_extract(source_evidence,'$.reason') like 'No unambiguous%'",
        'source_reason_competing_interpretations':"json_extract(source_evidence,'$.reason') like '%competing%'",
        'source_evidence_not_recorded':"source_evidence='null'",
        'ambiguous_scope_artifact_flag':"match_level='SAME_REVISION_SCOPE_OR_SITE_AMBIGUOUS' and exists (select 1 from json_each(annotations) a where json_extract(a.value,'$.ligand_is_artifact')=1)",
        'ambiguous_scope_proper_flag':"match_level='SAME_REVISION_SCOPE_OR_SITE_AMBIGUOUS' and exists (select 1 from json_each(annotations) a where json_extract(a.value,'$.ligand_is_proper')=1)",
        'ambiguous_scope_ion_flag':"match_level='SAME_REVISION_SCOPE_OR_SITE_AMBIGUOUS' and exists (select 1 from json_each(annotations) a where json_extract(a.value,'$.ligand_is_ion')=1)",
    }.items()}
    out['role_only_match'] = rows(f'''select match_level,annotation_interpretation,count(*) systems,
        count(distinct instance_key) instances,count(distinct evidence_key) evidence_groups
        from result where {only} group by 1,2 order by systems desc''')
    out['role_only_class_match'] = rows(f'''select primary_class,match_level,annotation_interpretation,count(*) systems
        from result where {only} group by 1,2,3 order by systems desc''')
    out['all_match'] = rows('''select disposition,match_level,annotation_interpretation,count(*) systems
        from result group by 1,2,3 order by systems desc''')
    candidate_systems = 0
    instances, evidence = set(), set()
    classes, codes, covalent = Counter(), Counter(), Counter()
    competing = 0
    for row in c.execute(f'select * from result where {only} and match_level=\'SAME_REVISION_INSTANCE_CHAIN_MATCH\' and annotation_interpretation=\'PLINDER_PROPER_RULE_SUPPORT\''):
        if eligible_candidate(row):
            candidate_systems += 1; instances.add(row['instance_key']); evidence.add(row['evidence_key'])
            classes[row['primary_class']] += 1
            annotations = json.loads(row['annotations'])
            codes.update({v['ligand_ccd_code'] for v in annotations})
            covalent[str(sorted({str(v['ligand_is_covalent']) for v in annotations}))] += 1
            competing += 'competing' in ((json.loads(row['source_evidence']) or {}).get('reason') or '')
    out['frozen_rule_candidates'] = {'systems':candidate_systems,'instances':len(instances),'evidence_groups':len(evidence)}
    out['candidate_classes'] = dict(classes)
    out['candidate_ccd_counts'] = dict(codes.most_common())
    out['candidate_covalent_flags'] = dict(covalent)
    out['candidate_competing_interpretations'] = competing
    with a.output.open('x',encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(json.dumps({k:out[k] for k in ('all','role_only','disposition','role_only_prior_current_conflicts','role_only_match','frozen_rule_candidates')},ensure_ascii=False))


if __name__ == '__main__':
    main()

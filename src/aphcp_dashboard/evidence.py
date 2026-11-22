"""Read bounded local APHCP receipts; registry data is a separate source."""
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
M1 = 'M1_full_scope_20261121'
PILOT = 'M2_2021-09_Pilot_20261121'
V1 = 'Holo_V1_五例验收_20261121'
COPY = 'Holo_V1_3TCT复制补验_20261121'


def load(root=ROOT):
    sources, issues = {}, []
    def receipt(directory, filename, key):
        p = root/'reports'/'inventory'/directory/filename
        try:
            if p.stat().st_size > 1_000_000:
                raise ValueError('超过本地摘要大小限制')
            value = json.loads(p.read_text(encoding='utf-8-sig'))
            if not isinstance(value, dict):
                raise ValueError('不是 JSON 对象')
            sources[key] = {'path': str(p.relative_to(root)), 'mtime': datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat()}
            return value
        except (OSError, ValueError) as err:
            issues.append(f'{directory}/{filename}：{getattr(err, "strerror", None) or str(err)}')
            return {}

    inv = receipt(M1, 'full_inventory_summary.json', 'm1')
    cross = receipt(M1, 'cross_window_summary.json', 'windows')
    audit = receipt(M1, 'execution_validation.json', 'm1_audit')
    plan = receipt(M1, 'plan.json', 'm1_plan')
    report = root/'reports'/'inventory'/M1/'全库执行回执.md'
    try:
        if report.stat().st_size > 100_000:
            raise ValueError('超过文本回执大小限制')
        match = re.search(r'坐标下载\s+(?:为\s*)?\*{0,2}([0-9,]+)', report.read_text(encoding='utf-8-sig'))
        if match:
            sources['m1_download'] = {'path': str(report.relative_to(root)), 'mtime': datetime.fromtimestamp(report.stat().st_mtime, timezone.utc).isoformat()}
        downloaded = int(match.group(1).replace(',', '')) if match else None
    except (OSError, ValueError):
        downloaded = None
    workflow = {}
    state = root/'PROJECT_STATE.md'
    try:
        for line in state.read_text(encoding='utf-8-sig').splitlines():
            cells = [part.strip() for part in line.strip('|').split('|')]
            if len(cells) == 2 and cells[0] in ('当前待办', '授权范围'):
                workflow[cells[0]] = cells[1]
        sources['workflow'] = {'path': str(state.relative_to(root)),
                               'mtime': datetime.fromtimestamp(state.stat().st_mtime, timezone.utc).isoformat()}
    except (OSError, ValueError):
        issues.append('无法读取 PROJECT_STATE.md 的当前待办和授权范围。')
    planned = {x['run_directory']:x for x in plan.get('windows', [])}
    windows = [dict(w, reuse=planned.get(w.get('run_directory'), {}).get('reuse_status')) for w in cross.get('windows', [])]
    pilot = receipt(PILOT, 'summary.json', 'pilot')
    pilot_run = receipt(PILOT, 'run_receipt.json', 'pilot_run')
    diag = receipt(PILOT, 'diagnostic_audit.json', 'pilot_diag')
    v1 = receipt(V1, 'verification.json', 'v1')
    copy = receipt(COPY, 'verification.json', 'copy')
    diagnostic = {x['pdb_id']:x for x in diag.get('entries', [])}
    cases = []
    for row in pilot.get('entry_results', []):
        pdb = row['pdb_id']; d = diagnostic.get(pdb, {}); maps = d.get('index_mapping_details', [])
        fidelity_issue = any(x.get('omitted_positive_altloc') or x.get('occupancy_changes') for x in maps)
        coords = 'PASS' if maps and all(x.get('state') == 'coordinate_identity_pass' for x in maps) else None
        cases.append(dict(pdb=pdb, group='9 月 M2 Pilot', date=pilot.get('finished_utc'), status=row.get('status'),
                          revision=row.get('revision'), assembly=d.get('assembly_ids'), models=d.get('models'),
                          coordinate=coords, fidelity='DIFFERENCE' if fidelity_issue else None, topology=None, role=None,
                          note='；'.join(filter(None, [row.get('review'), row.get('failure'),
                              '索引代表文件存在占有率或替代构象差异。' if fidelity_issue else ''])),
                          run=pilot_run.get('run_path'), artifacts={'原件':row.get('raw'), '原件坐标见证':row.get('native_structure'),
                                  '逐条 QC':row.get('qc_json')}, detail=d, source='pilot'))
    for pdb, row in v1.get('entries', {}).items():
        cases.append(dict(pdb=pdb, group='Holo Master V1 五例', date=v1.get('start_build_utc'), status=row.get('structure_acceptance'),
                          revision=row.get('revision'), assembly=row.get('assembly'), models=row.get('representative_source_model'),
                          coordinate=row.get('coordinate_identity'), fidelity=row.get('experimental_fidelity'),
                          topology=row.get('topology'), role=row.get('role_qc'),
                          note='SWPW 的序列与共价连接已定位，生物角色仍待核。' if pdb=='5R42' and row.get('role_qc')=='REVIEW' else '',
                          run=v1.get('run'), artifacts={k:v.get('path') for k,v in row.get('files', {}).items()},
                          detail=row.get('checks'), source='v1'))
    if copy.get('artifacts'):
        system = copy['artifacts'].get('system', {})
        prep = copy.get('invocations', {}).get('build', {})
        upstream = copy.get('source', {}).get('path', '')
        revision = re.search(r'-v([0-9]+\.[0-9]+)\.cif\.gz$', upstream)
        models = system.get('model_ids')
        assembly = re.search(r'-a([^-]+)-model', copy['artifacts'].get('complete_assembly', {}).get('path',''))
        cases.append(dict(pdb='3TCT', group='V1 组装复制补验', date=prep.get('finished_utc'), status=copy.get('status'),
                          revision=revision.group(1) if revision else None, assembly=assembly.group(1) if assembly else None,
                          models=models, coordinate=system.get('coordinate_identity'), fidelity=system.get('experimental_fidelity'),
                          topology=system.get('topology'), role=copy.get('role_qc'),
                          note='复制组装若带显式连接或碱基配对，具体端点映射仍待补验。3TCT 原件无显式连接。',
                          run=copy.get('run'), artifacts={k:v.get('path') for k,v in copy['artifacts'].items()},
                          detail=copy.get('checks'), source='copy'))
    return dict(m1=inv, cross=cross, audit=audit, plan=plan, windows=windows, downloaded=downloaded,
                pilot=pilot, pilot_run=pilot_run, v1=v1, copy=copy, cases=cases, workflow=workflow,
                sources=sources, issues=issues)


def accepted(e):
    return e['cross'].get('status') == e['audit'].get('cross_window_status') == 'INVENTORY_COMPLETE'


def newest_cases(e, identifier):
    return sorted((row for row in e['cases'] if row['pdb']==identifier), key=lambda row: row.get('date') or '', reverse=True)

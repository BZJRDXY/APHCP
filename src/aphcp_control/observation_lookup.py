"""Bounded read-only Dashboard view of the unified M3 staging catalog."""
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path


def catalog_view(root, controller, bundle, pdb):
    if not re.fullmatch(r'[A-Z0-9]{4}', pdb):
        raise ValueError('Unified M3 accepts legacy four-character PDB IDs')
    if not all(re.fullmatch(r'holo-m3-[A-Za-z0-9-]+', x) for x in (controller, bundle)):
        raise ValueError('Invalid staging Run')
    runs = (root.resolve()/'runs').resolve(strict=True)
    pointer = runs/controller/'catalog_head.json'
    code = runs/bundle/'code'
    if pointer.resolve()!=pointer or code.resolve()!=code:
        raise ValueError('Catalog path differs from selected staging Run')
    pending = not pointer.is_file()
    if pending:
        plan = pointer.with_name('increment_plan.json')
        if not plan.is_file() or plan.resolve()!=plan:
            return {'state':'pending', 'controller':controller, 'pdb_id':pdb, 'records':[]}
        seed = Path(json.loads(plan.read_text(encoding='utf-8'))['seed_run'])
        if seed.parent != runs:
            raise ValueError('Seed is outside the selected data root')
        label = seed.name
        pattern = r'holo-m3-[a-z0-9-]+-phase-\d{3}'
    else:
        label = json.loads(pointer.read_text(encoding='utf-8'))['catalog_run']
        pattern = re.escape(controller)+r'-phase-\d{3}'
    if (not re.fullmatch(pattern, label) or (runs/label).resolve()!=runs/label or
            not (runs/label/'phase_receipt.json').is_file()):
        raise ValueError('Invalid catalog head')
    sys.path[:0] = [str(code/'scripts'),str(code/'src')]
    from holo_m3_observations import lookup_records

    rows = lookup_records(runs/label,pdb=pdb)
    rows.sort(key=lambda r: (not bool(r.get('files')), not bool(r.get('system_id')), r.get('system_id','')))
    counts = {}
    for row in rows:
        counts[row['status']] = counts.get(row['status'], 0) + 1
    results = []
    for row in rows[:16]:
        ident = row.get('identity',{}); cls = row.get('type') or {}; evidence = row.get('entry_evidence',{})
        results.append({k:row[k] for k in ('system_id','status','role','decision_state','reason') if k in row} |
            {'candidate_key':(row.get('key') or [None,None,None])[2],
             'revision':ident.get('revision'),'model_scope':ident.get('model_scope'),'assembly':ident.get('assembly'),
             'chemical_type':cls.get('primary_class'),'type_annotation':cls.get('annotation_status'),
             'files':{kind:asset['path'] for kind,asset in row.get('files',{}).items()},
             'reasons':row.get('m3_reasons',[])[:4], 'failure_stage':evidence.get('failure_stage'),
             'source_run':row.get('source_run')})
    return {'state':'found' if rows else 'absent','pdb_id':pdb,'catalog_run':label,'records':results,
            'using_previous_catalog':pending,'matching':len(rows),'omitted':max(0,len(rows)-len(results)),'status_counts':counts}


def fetch_observations(settings, identifier):
    import shlex
    from pathlib import PurePosixPath
    from .inventory_lookup import canonical_id
    from .queries import query_args
    from .transport import QueryError, request_json

    _, identifier = query_args('pdb', identifier)
    identifier = canonical_id(identifier)
    catalog = settings.get('dashboard', {}).get('observation_catalog') or {}
    controller, bundle = (catalog.get(k, '') for k in ('controller', 'code_run'))
    python = PurePosixPath(catalog.get('python_relative_to_code_root', ''))
    if not controller:
        raise QueryError('未配置', '尚未指定 M3 暂存目录。')
    if (not all(re.fullmatch(r'holo-m3-[A-Za-z0-9-]+', p) for p in (controller,bundle)) or
            python.is_absolute() or '..' in python.parts or python.name != 'python'):
        raise QueryError('配置异常', 'M3 暂存 Run 或 Python 相对路径不正确。')
    if not re.fullmatch(r'[A-Z0-9]{4}', identifier):
        raise QueryError('暂不支持', '当前暂存目录使用四位 PDB ID。')
    server = settings['server']
    command = shlex.join([str(PurePosixPath(server['code_root'])/python), '-B', '-', server['data_root'],
                          controller, bundle, identifier])
    program = Path(__file__).read_text(encoding='utf-8-sig')
    return request_json(settings, command, stdin=program, limit=48_000)


if __name__=='__main__':
    try:
        result=catalog_view(Path(sys.argv[1]),*sys.argv[2:])
        print(json.dumps({'data':result,'queried_at':datetime.now(timezone.utc).isoformat()},ensure_ascii=True))
    except (OSError, ValueError, KeyError, RuntimeError) as exc:
        sys.exit(f'M3 read-only catalog: {exc}')

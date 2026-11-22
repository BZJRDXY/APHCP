"""Export assemblies with source atom lineage and experimental-value readback."""
from collections import Counter, defaultdict
from math import dist, fsum
from sys import float_info

import gemmi

from holo_candidate_batch import sha256
from holo_source_mapping import MappingError, assembly_instances

ATOM = '_atom_site.'
CONN = '_struct_conn.'
PAIR = '_ndb_struct_na_base_pair.'
MODEL = 'pdbx_PDB_model_num'
LABEL = 'label_asym_id'
COORD = ('Cartn_x', 'Cartn_y', 'Cartn_z')
COVALENT = {'covale', 'disulf', 'modres'}
PROJECTABLE_CONNECTIONS = COVALENT | {'metalc', 'hydrog'}
INSTANCE = '_aphcp_assembly_instance.'
LINEAGE = '_aphcp_atom_source.'
CONN_LINEAGE = '_aphcp_conn_source.'
PAIR_LINEAGE = '_aphcp_base_pair_source.'
TRANSFORM_DECIMALS = 6


def rows(block, category):
    cols = block.get_mmcif_category(category, raw=True)
    return [dict(zip(cols, row)) for row in zip(*cols.values())]


def text(token):
    return gemmi.cif.as_string(token) if token else ''


def put(block, category, records):
    if records:
        keys = dict.fromkeys(k for row in records for k in row)
        block.set_mmcif_category(category, {k: [r.get(k,'.') for r in records] for k in keys}, raw=True)


def describe_instances(structure, assembly, model):
    raw = assembly_instances(structure, assembly, model, include_water=True)
    labels = [label for op, label, _ in raw.values()]
    if not labels:
        raise MappingError('assembly_members_unobserved', 'No observed members for the selected Assembly')
    expand = len(labels) != len(set(labels)) or any(not op.transform.is_identity() for op, _, _ in raw.values())
    instances = {}
    for key, (op, label, rr) in raw.items():
        auth = sorted({author for author, _ in rr})
        item = {'label_asym': label, 'operator': op.name, 'entity': rr[0][1].entity_id, 'auth_asym': auth}
        if expand:
            if len(auth) != 1:
                raise MappingError('assembly_author_chain_ambiguous', 'An expanded member has multiple source author chains', instance=key)
            output = 'i'+key.replace('.', '_', 1)
            item.update(output_label_asym=output, output_auth_asym=output, matrix=op.transform.mat.tolist(),
                        vector=op.transform.vec.tolist())
        instances[key] = item
    return instances, expand


def output_label(member):
    return member.get('output_label_asym', member['label_asym'])


def expanded_atoms(source, instances, *, source_atoms=None):
    """Materialize declared Assembly members exactly once, with stable source/atom keys."""
    by_label = defaultdict(list)
    for key, member in instances.items():
        by_label[member['label_asym']].append((key, member))
    projected, lineage, originals = [], [], {}
    for original in (rows(source, ATOM) if source_atoms is None else source_atoms):
        src_label, model, orig_id = (text(original[k]) for k in (LABEL, MODEL, 'id'))
        if src_label not in by_label:
            continue
        if (model, orig_id) in originals:
            raise MappingError('source_atom_identity_nonunique', 'Source model and atom_site.id must identify one atom')
        originals[model, orig_id] = original
        for key, member in by_label[src_label]:
            row = dict(original)
            row['id'] = str(len(projected)+1)
            row[LABEL] = row['auth_asym_id'] = output_label(member)
            if (member['matrix'], member['vector']) != (IDENTITY, ZERO):
                xyz = [float(text(original[axis])) for axis in COORD]
                for n, axis in enumerate(COORD):
                    value = sum(f*x for f, x in zip(member['matrix'][n], xyz)) + member['vector'][n]
                    row[axis] = f'{value:.{TRANSFORM_DECIMALS}f}'
            projected.append(row)
            lineage.append({'atom_site_id': row['id'], 'source_model_id': model, 'source_atom_site_id': orig_id,
                            'assembly_instance_id': gemmi.cif.quote(key), 'source_label_asym_id': original[LABEL],
                            'source_auth_asym_id': original['auth_asym_id']})
    return projected, lineage, originals


IDENTITY = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
ZERO = [0., 0., 0.]


def verify_transformed(atoms, lineage, originals, instances, model_ids):
    """Verify each copied atom against its original row and the declared operator."""
    mapping = {text(r['atom_site_id']): r for r in lineage}
    if len(mapping) != len(atoms):
        raise ValueError('Disk assembly atom lineage must be one-to-one on output atom IDs')
    residuals = defaultdict(lambda: {'atom_count': 0, 'max_coordinate_error_A': 0., 'max_3d_error_A': 0., 'max_axis_bound_A': 0.})
    for row in atoms:
        link = mapping[text(row['id'])]
        model = text(link['source_model_id'])
        original = originals[model, text(link['source_atom_site_id'])]
        key = text(link['assembly_instance_id'])
        member = instances[key]
        assigned = output_label(member)
        if (model not in model_ids or text(row[MODEL]) != model or text(row[LABEL]) != assigned or
                text(row['auth_asym_id']) != assigned or text(original[LABEL]) != member['label_asym'] or
                original[LABEL] != link['source_label_asym_id'] or original['auth_asym_id'] != link['source_auth_asym_id'] or
                any(row[k] != v for k, v in original.items() if k not in {*COORD, 'id', LABEL, 'auth_asym_id'})):
            raise ValueError('A copied atom lost original identity or experimental fields')
        xyz = [float(text(original[axis])) for axis in COORD]
        actual = [float(text(row[axis])) for axis in COORD]
        expected = []
        metrics = residuals[key]
        metrics['atom_count'] += 1
        for n in range(3):
            terms = [f*x for f, x in zip(member['matrix'][n], xyz)] + [member['vector'][n]]
            value = fsum(terms)
            error = abs(actual[n]-value)
            bound = .5 * 10**(-TRANSFORM_DECIMALS) + 8*float_info.epsilon*max(1., sum(map(abs, terms)))
            if error > bound:
                raise ValueError(f'Disk assembly coordinate failed {key} atom {row["id"]}: {error} > {bound} A')
            expected.append(value)
            metrics['max_coordinate_error_A'] = max(error, metrics['max_coordinate_error_A'])
            metrics['max_axis_bound_A'] = max(bound, metrics['max_axis_bound_A'])
        metrics['max_3d_error_A'] = max(dist(actual, expected), metrics['max_3d_error_A'])
    return dict(residuals)


def materialized_assembly_categories(categories, source, assembly, labels):
    """Coordinates are already expanded. Advertise one identity transform, never re-expand."""
    identity_rows(categories, '_pdbx_struct_assembly.', 'id', {assembly}, source)
    categories['_pdbx_struct_assembly_gen.'] = [{'assembly_id': gemmi.cif.quote(assembly), 'oper_expression': '1',
                                                'asym_id_list': gemmi.cif.quote(','.join(sorted(labels)))}]
    op = {'id': '1', 'type': gemmi.cif.quote('identity operation')}
    for n in range(3):
        for j in range(3):
            op[f'matrix[{n+1}][{j+1}]'] = str(IDENTITY[n][j])
        op[f'vector[{n+1}]'] = '0'
    categories['_pdbx_struct_oper_list.'] = [op]


def local_symmetry(value):
    return text(value) in ('', '1_555')


def choose_connected(block, selected):
    """Close over deposited local covalent/disulfide links; do not guess other components."""
    edges = [r for r in rows(block, CONN) if text(r['conn_type_id']) in COVALENT]
    chosen = set(selected)
    while True:
        expanded = set(chosen)
        for r in edges:
            ends = {text(r[f'ptnr{j}_label_asym_id']) for j in (1, 2)}
            if ends & chosen and all(local_symmetry(r.get(f'ptnr{j}_symmetry')) for j in (1, 2)):
                expanded |= ends
        if expanded == chosen:
            return chosen
        chosen = expanded


def atom_index(atoms):
    out = defaultdict(list)
    for row in atoms:
        out[text(row[MODEL]), text(row[LABEL]), text(row['label_comp_id']), text(row['label_atom_id'])].append(row)
    return out


def endpoint_field(j, key):
    prefix = 'pdbx_' if j == 3 or key in ('label_alt_id','PDB_ins_code') else ''
    return f'{prefix}ptnr{j}_{key}'


def resolve_end(edge, j, index, model):
    name = lambda k: endpoint_field(j, k)
    opts = index[model, *(text(edge.get(name(k))) for k in (LABEL, 'label_comp_id', 'label_atom_id'))]
    fields = {name(k): k for k in ('label_seq_id','auth_asym_id','auth_comp_id','auth_seq_id')}
    fields[name('PDB_ins_code')] = 'pdbx_PDB_ins_code'
    alt = text(edge.get(name('label_alt_id')))
    # An absent altloc is resolved only when the other fields identify a single atom.
    matches = [a for a in opts if (not alt or text(a['label_alt_id']) == alt)
               and all(not text(edge.get(k)) or text(edge[k]) == text(a.get(v)) for k, v in fields.items())]
    if len(matches) != 1:
        raise MappingError('connection_endpoint_unresolved', 'Connection must resolve to one atom per endpoint in each exported model',
                           connection=text(edge['id']), model=model, endpoint=j, matches=len(matches))
    a = matches[0]
    return {k: text(a.get(k)) or None for k in ('id', LABEL, 'auth_asym_id', 'label_comp_id', 'label_seq_id',
                                                'auth_seq_id', 'pdbx_PDB_ins_code', 'label_atom_id', 'label_alt_id')}


def connection_scope(block, labels):
    keep, omitted = [], []
    for row in rows(block, CONN):
        ends = {text(row[f'ptnr{j}_label_asym_id']) for j in (1, 2)}
        third = text(row.get('pdbx_ptnr3_label_asym_id'))
        if third:
            ends.add(third)
        if not ends & labels:
            continue
        local = all(local_symmetry(row.get(f'ptnr{j}_symmetry')) for j in (1, 2))
        if text(row['conn_type_id']) in COVALENT and (not local or not ends <= labels or third):
            raise MappingError('connection_scope_unsupported', 'Cannot represent an included covalent connection in this projection',
                               connection=text(row['id']), labels=sorted(ends))
        if ends <= labels and local and not third:
            keep.append(row)
        else:
            omitted.append({'id': text(row['id']), 'kind': text(row['conn_type_id']), 'external_labels': sorted(ends-labels),
                            'reason': 'partner outside scope' if ends-labels else 'nonlocal symmetry or extra endpoint'})
    return keep, omitted


def read_links(records, atoms, models, *, expanded=False):
    idx = atom_index(atoms)
    links = []
    for r in records:
        # Archive categories like hydrogen bonds are preserved as deposited; mandatory
        # atom-level checks apply to covalent/disulfide and ionic coordination rows.
        if not expanded and text(r['conn_type_id']) not in COVALENT | {'metalc'}:
            continue
        for m in models:
            links.append({'id': text(r['id']), 'kind': text(r['conn_type_id']), 'model': int(m),
                          'ends': [resolve_end(r, j, idx, m) for j in ((1, 2, 3) if text(r.get(endpoint_field(3,LABEL))) else (1, 2))]})
    return links



def deposited_copy_matches(source, expanded, qualifiers, labels, *, equivalent_embeddings=False):
    """Find an Assembly embedding of deposited symmetry partners; never use proximity.

    For a nonlocal relation the directly deposited embedding is A_j = S_j.
    Additional equivalent embeddings require independently verified deposition
    evidence and are intentionally not inferred here.
    """
    qualifiers = tuple(qualifiers)
    if all(sym == '1_555' for sym in qualifiers):
        ops = dict.fromkeys(member['operator'] for member in expanded.values() if member['label_asym'] in labels)
        return [[[(key, m) for key, m in expanded.items() if m['label_asym'] == label and m['operator'] == op]
                 for label in labels] for op in ops]
    if any(not sym for sym in qualifiers):
        return []
    operators = defaultdict(set)
    for row in rows(source, '_pdbx_struct_oper_list.'):
        sym = text(row.get('name')); kind = text(row.get('type'))
        # A source-named point operation is usable only for covalent rows, whose
        # transformed endpoint separation is checked against the deposited distance.
        allowed = {'identity operation', 'crystal symmetry operation'}
        if equivalent_embeddings:
            allowed.add('point symmetry operation')
        if sym and kind in allowed and text(row.get('symmetry_operation')):
            operators[sym].add(text(row['id']))
    if not all(len(operators[sym]) == 1 for sym in qualifiers):
        return []
    direct = [[(key, m) for key, m in expanded.items() if m['label_asym'] == label and m['operator'] in operators[sym]]
              for label, sym in zip(labels, qualifiers)]
    if not equivalent_embeddings or any(len(side)!=1 for side in direct):
        return [direct]
    # A mapped connection can be copied under an Assembly displacement D only
    # when every partner lands on an already declared instance: A_j = D S_j.
    def apply(M,v): return [fsum(x*y for x,y in zip(row,v)) for row in M]
    def mt(M,N): return [apply(list(map(list,zip(*N))), row) for row in M]
    def add(l,r): return [x+y for x,y in zip(l,r)]
    def subtract(l,r): return [x-y for x,y in zip(l,r)]
    def close(l,r): return all(abs(x-y)<1e-6 for x,y in zip(l,r))
    _, pivot = direct[0][0]
    invert = list(map(list,zip(*pivot['matrix'])))  # orthogonal rigid symmetry
    if not all(close(row, want) for row,want in zip(mt(pivot['matrix'],invert),IDENTITY)):
        return [direct]
    existing = defaultdict(list)
    for key,m in expanded.items(): existing[m['label_asym']].append((key,m))
    matches = [direct]
    for key, orig in existing[labels[0]]:
        if key == direct[0][0][0]:
            continue
        D=mt(orig['matrix'],invert); t=subtract(orig['vector'],apply(D,pivot['vector']))
        chosen=[[(key,orig)]]
        for side,label in zip(direct[1:],labels[1:]):
            _,sym=side[0]; want=mt(D,sym['matrix']); v=add(apply(D,sym['vector']),t)
            chosen.append([(k,m) for k,m in existing[label] if close(m['vector'],v) and
                           all(close(row,mr) for row,mr in zip(want,m['matrix']))])
        matches.append(chosen)
    return matches


def expanded_pair_context(source, expanded, prepared, model_ids):
    out, blockers = [], []
    all_atoms = prepared[0]
    observed = {(text(a[MODEL]), text(a[LABEL]), text(a['label_seq_id']), text(a['label_comp_id']),
                 text(a['auth_seq_id']), text(a['pdbx_PDB_ins_code'])) for a in all_atoms}
    seen = set()
    for ordinal, row in enumerate(rows(source, PAIR), 1):
        model = text(row['model_number']); sides = ('i', 'j')
        if model not in model_ids:
            continue
        labels = [text(row.get(f'{side}_{LABEL}')) for side in sides]
        syms = [text(row.get(f'{side}_symmetry')) for side in sides]
        affected = {output_label(m) for m in expanded.values() if m['label_asym'] in labels}
        if not affected:
            continue
        rid = text(row.get('pair_number')) or str(ordinal)
        choices = deposited_copy_matches(source, expanded, syms, labels)
        if not choices:
            blockers.append(('expanded_pair_symmetry_unresolved', rid, affected, {'category':PAIR, 'symmetry':syms}))
        for group in choices:
            active = {output_label(m) for opts in group for _,m in opts} or affected
            if any(len(opts) != 1 for opts in group):
                blockers.append(('expanded_pair_instance_unresolved', rid, active, {'category':PAIR, 'symmetry':syms}))
                continue
            assigned = [opts[0] for opts in group]; new = dict(row); ends=[]
            for side, (key, member) in zip(sides, assigned):
                mapped = output_label(member)
                sig = (model, mapped, *(text(row.get(f'{side}_{k}')) for k in ('label_seq_id','label_comp_id','auth_seq_id','PDB_ins_code')))
                if sig not in observed:
                    break
                ends.append(sig)
                new[f'{side}_{LABEL}'] = mapped
                if f'{side}_auth_asym_id' in new:
                    new[f'{side}_auth_asym_id'] = mapped
                new[f'{side}_symmetry'] = '1_555'
            if len(ends)!=2 or ends[0]==ends[1]:
                blockers.append(('expanded_pair_residue_unresolved', rid, active, {'category':PAIR, 'symmetry':syms}))
                continue
            key = tuple(sorted(ends))
            if key in seen:
                blockers.append(('expanded_pair_relation_nonunique', rid, active, {'category':PAIR, 'symmetry':syms}))
                continue
            seen.add(key)
            output_id = str(len(out)+1)
            new['pair_number'] = output_id
            if 'pair_name' in new:
                new['pair_name'] = gemmi.cif.quote(':'.join(f'{x[1]}_{x[3]}{x[2]}' for x in ends))
            trace = {'output_pair_number':output_id, 'source_pair_number':gemmi.cif.quote(rid), 'source_model_id':model,
                     **{f'{side}_assembly_instance_id':gemmi.cif.quote(pair[0]) for side,pair in zip(sides,assigned)},
                     **{f'{side}_source_symmetry':gemmi.cif.quote(sym) for side,sym in zip(sides,syms)}}
            out.append((frozenset(active), new, trace))
    return out, blockers


def expanded_conn_context(source, expanded, prepared, model_ids):
    """Map ASU-local two-ended covalent rows through the same Assembly operation.

    Return a graph and per-instance blockers so unresolved copies do not suppress
    independent projections of the same Assembly.
    """
    by_label = defaultdict(list)
    for key, member in expanded.items():
        by_label[member['label_asym']].append((key, member))
    graph, blockers = [], []
    def hold(reason, connection, affected, **evidence):
        blockers.append((reason, connection, frozenset(affected), evidence))

    atoms, lineage, originals = prepared
    original_ix = atom_index(originals.values())
    atom_by_id = {a['id']: a for a in atoms}
    projected = {(text(l['source_model_id']), text(l['source_atom_site_id']), text(l['assembly_instance_id'])):
                 atom_by_id[l['atom_site_id']] for l in lineage}
    seen_mapped = set()
    for row in rows(source, CONN):
        src_id = text(row['id'])
        third = any(text(row.get(endpoint_field(3,k))) for k in (LABEL, 'label_comp_id', 'label_atom_id'))
        ends = (1,2,3) if third else (1,2)
        target_labels = [text(row.get(endpoint_field(j,LABEL))) for j in ends]
        affected = {output_label(m) for label in target_labels for _,m in by_label[label]}
        if not affected:
            continue
        syms = tuple(text(row.get(endpoint_field(j,'symmetry'))) for j in ends)
        kind = text(row['conn_type_id'])
        if kind not in PROJECTABLE_CONNECTIONS or len(model_ids) != 1:
            hold('expanded_topology_pending', src_id, affected, category=CONN, symmetry=syms, endpoints=len(ends))
            continue
        if third and (not all(syms) or not all(target_labels)):
            hold('expanded_multi_endpoint_missing_symmetry_or_identity', src_id, affected, category=CONN, symmetry=syms, endpoints=len(ends))
            continue
        model, = model_ids
        try:
            source_ends = [resolve_end(row, j, original_ix, model) for j in ends]
        except MappingError as e:
            hold(e.kind, src_id, affected, **{k: v for k,v in e.evidence.items() if k != 'connection'})
            continue
        choices_all = deposited_copy_matches(source, expanded, syms, target_labels, equivalent_embeddings=kind in COVALENT)
        if not choices_all:
            hold('expanded_conn_symmetry_unresolved', src_id, affected, category=CONN, symmetry=syms)
        for choices in choices_all:
            labels = {output_label(m) for side in choices for _,m in side} or affected
            if any(len(side) != 1 for side in choices):
                hold('expanded_connection_instance_unresolved', src_id, labels, symmetry=syms,
                     matches=[[k for k,_ in side] for side in choices])
                continue
            assigned = [side[0] for side in choices]
            local = all(sym == '1_555' for sym in syms)
            if local and any((m['matrix'],m['vector']) != (assigned[0][1]['matrix'],assigned[0][1]['vector']) for _,m in assigned):
                hold('expanded_connection_instance_unresolved', src_id, labels, symmetry=syms)
                continue
            endpoint_atoms = [projected[model,end['id'],key] for end,(key,_) in zip(source_ends,assigned)]
            if len({a['id'] for a in endpoint_atoms}) != len(ends):
                hold('expanded_connection_endpoints_not_distinct', src_id, labels, symmetry=syms)
                continue
            distance = dist(*[[float(text(a[axis])) for axis in COORD] for a in endpoint_atoms[:2]])
            expected = text(row.get('pdbx_dist_value'))
            if not local and ((expected and abs(distance-float(expected)) > .003) or (not expected and kind in COVALENT)):
                hold('expanded_connection_geometry_unresolved', src_id, labels, symmetry=syms,
                     expected_distance=expected, actual_distance=distance)
                continue
            signature = (src_id, model, frozenset(a['id'] for a in endpoint_atoms))
            if signature in seen_mapped:
                continue # The same source bond mapped in reverse is the same observed relation.
            seen_mapped.add(signature)
            out_id = f'aphcp_conn{len(graph)+1}'
            new = dict(row,id=out_id)
            trace = {'id':out_id,'source_conn_id':row['id'],'source_model_id':model,'endpoint_count':str(len(ends))}
            for j, source_end, (member_key,_), a, sym in zip(ends,source_ends,assigned,endpoint_atoms,syms):
                new[endpoint_field(j,LABEL)] = a[LABEL]
                for k,field in [('auth_asym_id','auth_asym_id'),('label_alt_id','label_alt_id')]:
                    if endpoint_field(j,k) in new:
                        new[endpoint_field(j,k)] = a[field]
                new[endpoint_field(j,'symmetry')] = '1_555'
                trace.update({f'ptnr{j}_source_atom_site_id':source_end['id'], f'ptnr{j}_output_atom_site_id':a['id'],
                              f'ptnr{j}_assembly_instance_id':gemmi.cif.quote(member_key), f'ptnr{j}_source_symmetry':sym})
            for key in ('source_atom_site_id','output_atom_site_id','assembly_instance_id','source_symmetry'):
                trace.setdefault('ptnr3_'+key,'.')
            link = {'id':src_id,'kind':kind,'model':int(model),'ends':source_ends,
                    'source_symmetry':list(syms),'assembly_instances':[key for key,_ in assigned]}
            graph.append((labels,new,trace,link))
    pairs, pair_blockers = expanded_pair_context(source, expanded, prepared, model_ids)
    return graph, blockers+pair_blockers, pairs


def expanded_conn_projection(context, labels, *, complete):
    graph, blockers, pairs = context
    chosen = set(labels)
    closure_edges = [ends for ends, _, _, link in graph if link['kind'] in COVALENT] + [ends for ends, *_ in pairs]
    while True:
        added = set().union(*(ends for ends in closure_edges if ends & chosen)) - chosen
        if not added:
            break
        if complete:
            raise MappingError('topology_crosses_assembly_scope', 'Covalent closure adds members outside the complete Assembly', labels=sorted(added))
        chosen |= added
    for reason, src_id, affected, evidence in blockers:
        if affected & chosen:
            raise MappingError(reason, 'An included Assembly connection cannot be mapped to unique output endpoints',
                               connection=src_id, affected_labels=sorted(affected & chosen), **evidence)
    return (chosen, [edge[1:] for edge in graph if edge[0] <= chosen], [edge[1:] for edge in pairs if edge[0] <= chosen],
            [{'id':edge[3]['id'],'kind':edge[3]['kind'],'external_labels':sorted(edge[0]-chosen),'reason':'partner outside scope'}
             for edge in graph if edge[0] & chosen and not edge[0] <= chosen])

def read_pairs(records, atoms):
    observed = {(text(a[MODEL]), text(a[LABEL]), text(a['label_seq_id']), text(a['label_comp_id'])) for a in atoms}
    pairs = []
    for row in records:
        model = text(row['model_number'])
        ends = [{k: text(row.get(f'{side}_{k}')) for k in (LABEL, 'label_seq_id', 'label_comp_id', 'symmetry')}
                for side in ('i', 'j')]
        if any((model, e[LABEL], e['label_seq_id'], e['label_comp_id']) not in observed for e in ends):
            raise ValueError('Deposited base pair component not found on disk')
        pairs.append({'pair': text(row['pair_number']), 'model': int(model), 'ends': ends})
    return pairs


def identity_rows(scope, category, column, keep, source):
    scope[category] = [row for row in rows(source, category) if text(row.get(column)) in keep]


def export(path, source, entry, assembly, models, labels, *, complete=False, verify=False, expanded=None, prepared=None, prepared_connections=None):
    model_ids = {str(m) for m in models}
    if expanded:
        all_atoms, lineage, originals = expanded_atoms(source, expanded) if prepared is None else prepared
    else:
        all_atoms = rows(source, ATOM) if prepared is None else prepared
    original_labels = set(labels)
    expanded_links, expanded_pairs, expanded_omitted = [], [], []
    if expanded:
        ctx = (expanded_conn_context(source, expanded, (all_atoms, lineage, originals), model_ids)
               if prepared_connections is None else prepared_connections)
        labels, expanded_links, expanded_pairs, expanded_omitted = expanded_conn_projection(ctx, labels, complete=complete)
    atoms = [r for r in all_atoms if text(r[MODEL]) in model_ids and text(r[LABEL]) in labels]
    counts = Counter((text(r[MODEL]), text(r[LABEL])) for r in atoms)
    if set(counts) != {(model, label) for model in model_ids for label in labels}:
        raise MappingError('model_member_coverage', 'An Assembly member is unobserved in one or more selected models')
    atom_ids = [(text(r[MODEL]), text(r['id'])) for r in atoms]
    if len(set(atom_ids)) != len(atom_ids):
        raise MappingError('source_atom_identity_nonunique', 'Source model and atom_site.id must identify an atom uniquely')
    conn, omitted = ([x[0] for x in expanded_links], expanded_omitted) if expanded else connection_scope(source, labels)
    pairs = [row for row, _ in expanded_pairs] if expanded else [r for r in rows(source, PAIR) if text(r['model_number']) in model_ids and
             all(text(r[f'{side}_{LABEL}']) in labels and local_symmetry(r.get(f'{side}_symmetry')) for side in ('i', 'j'))]
    if expanded and conn:
        keys = dict.fromkeys(k for row in conn for k in row)
        conn = [{k:r.get(k,'.') for k in keys} for r in conn]
    source_links = read_links(conn, atoms, sorted(model_ids, key=int), expanded=bool(expanded))
    source_pairs = read_pairs(pairs, atoms)
    entities = {text(r['label_entity_id']) for r in atoms}
    components = {text(r['label_comp_id']) for r in atoms}
    categories = {ATOM: atoms, CONN: conn, PAIR: pairs}
    for category, column, keys in (('_struct_asym.', 'id', labels), ('_entity.', 'id', entities),
                                   ('_entity_poly.', 'entity_id', entities), ('_entity_poly_seq.', 'entity_id', entities),
                                   ('_chem_comp.', 'id', components), ('_struct_conn_type.', 'id', {text(r['conn_type_id']) for r in conn})):
        identity_rows(categories, category, column, keys, source)
    if expanded:
        used_ids = {r['id'] for r in atoms}
        instances = {k: v for k, v in expanded.items() if output_label(v) in labels}
        categories[LINEAGE] = [r for r in lineage if r['atom_site_id'] in used_ids]
        categories[CONN_LINEAGE] = [edge[1] for edge in expanded_links]
        categories[PAIR_LINEAGE] = [trace for _, trace in expanded_pairs]
        categories[INSTANCE] = [{'id': gemmi.cif.quote(k), 'output_label_asym_id': output_label(v),
                                 'output_auth_asym_id': v['output_auth_asym'], 'source_label_asym_id': gemmi.cif.quote(v['label_asym']),
                                 'source_auth_asym_id': gemmi.cif.quote(v['auth_asym'][0]),
                                 'source_assembly_id': gemmi.cif.quote(assembly), 'source_operator_id': gemmi.cif.quote(v['operator'])}
                                for k, v in instances.items()]
        asym = {text(r['id']): r for r in rows(source, '_struct_asym.')}
        categories['_struct_asym.'] = [dict(asym[item['label_asym']], id=output_label(item)) for item in instances.values()]
        for rec in categories['_entity_poly.']:
            strand_list = sorted(output_label(item) for item in instances.values() if item['entity'] == text(rec['entity_id']))
            if strand_list and 'pdbx_strand_id' in rec:
                rec['pdbx_strand_id'] = gemmi.cif.quote(','.join(strand_list))
    if complete and expanded:
        materialized_assembly_categories(categories, source, assembly, labels)
    elif complete:
        identity_rows(categories, '_pdbx_struct_assembly.', 'id', {assembly}, source)
        identity_rows(categories, '_pdbx_struct_assembly_gen.', 'assembly_id', {assembly}, source)
        categories['_pdbx_struct_oper_list.'] = rows(source, '_pdbx_struct_oper_list.')
    work = path if verify else path.with_suffix('.cif.tmp')
    if not verify:
        path.parent.mkdir(parents=True, exist_ok=True)
        doc = gemmi.cif.Document(); block = doc.add_new_block(entry)
        block.set_pair('_entry.id', entry)
        for category, selected in categories.items():
            put(block, category, selected)
        doc.write_file(str(work))
    disk = gemmi.cif.read_file(str(work)).sole_block()
    actual = rows(disk, ATOM)
    structure = gemmi.make_structure_from_block(disk)
    if [m.num for m in structure] != list(models) or sum(m.count_atom_sites() for m in structure) != len(atoms):
        raise ValueError('Gemmi model numbers or atom count disagree with selected source rows')
    if text(disk.find_value('_entry.id')) != entry or any(rows(disk, cat) != selected for cat, selected in categories.items()):
        raise ValueError('CIF readback differs from the selected source rows (atom IDs, geometry, experimental values or topology)')
    disk_links = read_links(rows(disk, CONN), actual, sorted(model_ids, key=int), expanded=bool(expanded))
    disk_pairs = read_pairs(rows(disk, PAIR), actual)
    if disk_links != source_links or disk_pairs != source_pairs:
        raise ValueError('Resolved disk connection/base-pair identities differ from source')
    if expanded:
        by_id = {link['id']: link for link in disk_links}
        by_atom = {text(a['atom_site_id']): a for a in rows(disk, LINEAGE)}
        for tr in rows(disk, CONN_LINEAGE):
            output = by_id[text(tr['id'])]
            if len(output['ends']) != int(text(tr['endpoint_count'])):
                raise ValueError('Disk connection endpoint count differs from source mapping')
            for j, end in enumerate(output['ends'], 1):
                origin = by_atom[end['id']]
                if (text(tr[f'ptnr{j}_output_atom_site_id']) != end['id'] or
                        origin['source_model_id'] != tr['source_model_id'] or
                        origin['source_atom_site_id'] != tr[f'ptnr{j}_source_atom_site_id'] or
                        origin['assembly_instance_id'] != tr[f'ptnr{j}_assembly_instance_id']):
                    raise ValueError('Disk connection endpoint does not match its source atom and Assembly copy')
    if complete:
        for model in structure:
            disk_instances, still_needs_expansion = describe_instances(structure, assembly, model)
            if still_needs_expansion or {output_label(x) for x in disk_instances.values()} != labels:
                raise ValueError('Complete Assembly metadata does not reproduce its on-disk members')
    transform_qc = verify_transformed(actual, rows(disk, LINEAGE), originals, expanded, model_ids) if expanded else None
    if not verify:
        work.replace(path)
    present_alt = Counter(text(r['label_alt_id']) for r in actual if text(r['label_alt_id']))
    artifact = {'path': str(path), 'sha256': sha256(path), 'bytes': path.stat().st_size, 'model_ids': list(models),
            'model_members': {m: {label: counts[str(m), label] for label in sorted(labels)} for m in models},
            'atom_count': len(actual), 'connections': [{'source': link, 'output': other} for link, other in zip(source_links, disk_links)],
            'base_pairs': disk_pairs, 'omitted_touching_connections': omitted,
            'atom_mapping': 'Output model + atom_site.id, label/auth chain, residue/ins, atom and altloc are unchanged from the fixed source; Assembly operator is in the instance map.',
            'omitted_models': sorted({int(text(r[MODEL])) for r in all_atoms}-set(models)),
            'omitted_included_model_labels': sorted({text(r[LABEL]) for r in all_atoms if text(r[MODEL]) in model_ids}-labels),
            'omitted_source_atom_rows': len(all_atoms)-len(atoms), 'altloc_policy': 'all deposited alternatives for selected members/models',
            'named_altloc_atoms': dict(present_alt), 'coordinate_identity': 'PASS', 'experimental_fidelity': 'PASS', 'topology': 'PASS'}
    if expanded:
        artifact.update(atom_mapping='CIF _aphcp_atom_source maps each output atom ID to source model and atom_site.id; _aphcp_assembly_instance maps every chain to its original Assembly operation.',
                        output_instances=instances, coordinate_transform_readback=transform_qc,
                        transform_coordinate_decimals=TRANSFORM_DECIMALS,
                        explicit_source_connection_rows=len({x[2]['id'] for x in expanded_links}),
                        included_by_topology_closure=sorted(labels-original_labels),
                        included_by_covalent_closure=sorted((labels-original_labels) &
                            {end[LABEL] for link in disk_links if link['kind'] in COVALENT for end in link['ends']}),
                        omitted_assembly_atom_rows=artifact.pop('omitted_source_atom_rows'))
        artifact['base_pair_source_trace'] = [{k:text(v) for k,v in trace.items()} for _,trace in expanded_pairs]
        artifact['connections'] = [{'source': edge[2], 'output': disk_link}
                                   for edge, disk_link in zip(expanded_links, disk_links)]
    return artifact, actual



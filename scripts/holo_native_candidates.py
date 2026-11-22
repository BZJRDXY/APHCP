"""Inspect experimentally observed candidates independently of index coverage."""

from collections import defaultdict
from hashlib import sha256
import json

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from holo_source_mapping import assembly_instances, connection_edges, observed_endpoint, subchains, xyz


def key_for(entry, revision, model, assembly, members):
    identity = {'entry': entry, 'revision': revision, 'model': model, 'assembly': assembly,
                'members': sorted(members)}
    return sha256(json.dumps(identity, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:24]


def member_key(instance, label, operator, entity, author, seq, comp):
    return [instance, label, operator, entity, author, seq, comp]


def indexed_member_keys(instances, names):
    return [member_key(name, label, op.name, str(res.entity_id), author, str(res.seqid), res.name)
            for name in names for op, label, rows in [instances[name]] for author, res in rows]


def columns(block, category, key, value):
    c = block.get_mmcif_category(category)
    return dict(zip(c.get(key, []), c.get(value, [])))


def inspect_source(block, structure, entry, revision, rules, *, scope=None, calculate_distances=True, preserve_protein_context=False):
    descriptions = columns(block, '_entity.', 'id', 'pdbx_description')
    polymer_types = columns(block, '_entity_poly.', 'entity_id', 'type')
    chemical_names = columns(block, '_chem_comp.', 'id', 'name')
    explicit = list(connection_edges(block))
    contexts = []
    for model in structure:
        for assembly in [a.name for a in structure.assemblies] or [None]:
            if scope is not None and (model.num, assembly) not in scope:
                continue
            instances = (assembly_instances(structure, assembly, model) if assembly else
                         {f'asu.{k}': (None, k, v) for k, v in subchains(model).items()})
            contexts.append((model, assembly, instances))

    def points(op, residues):
        return np.array([xyz(op.transform.apply(a.pos) if op else a.pos) for _, r in residues
                         for a in r if not a.element.is_hydrogen and a.occ > 0]).reshape((-1, 3))

    output = []
    for model, assembly, instances in contexts:
        entity_by_chain = {name: str(rows[0][1].entity_id) for name, (_, _, rows) in instances.items()}
        proteins = [name for name in instances if str(polymer_types.get(entity_by_chain[name], '')).startswith('polypeptide')]
        eligible_proteins = [p for p in proteins if rules.get(f'entity:{entity_by_chain[p]}', ['', '', ''])[1]
                             not in {'Fab partner', 'antibody partner', 'protein partner', 'polymer ligand'}]
        polymer_entities = {entity_by_chain[name] for name in instances if entity_by_chain[name] in polymer_types}
        items, unresolved = [], []
        chain_items = defaultdict(list)
        for name, (op, label, residues) in sorted(instances.items()):
            entity = entity_by_chain[name]
            operator = op.name if op else None
            if entity in polymer_types:
                if len(polymer_entities) > 1 and rules.get(f'entity:{entity}', ['', '', ''])[1] != 'receptor':
                    items.append((name, 'polymer', entity, op, label, residues))
                    chain_items[name].append(len(items) - 1)
                for author, res in residues:
                    if res.het_flag == 'H' and any(a.occ > 0 and not a.element.is_hydrogen for a in res):
                        items.append((name, 'polymer_modification', res.name, op, label, [(author, res)]))
                        chain_items[name].append(len(items) - 1)
            elif residues[0][1].entity_type == gemmi.EntityType.Branched:
                items.append((name, 'branched', entity, op, label, residues))
                chain_items[name].append(len(items) - 1)
            else:
                for author, res in residues:
                    items.append((name, 'component', res.name, op, label, [(author, res)]))
                    chain_items[name].append(len(items) - 1)

        parent = list(range(len(items)))
        def find(i):
            while parent[i] != i:
                i = parent[i]
            return i

        present_links = []
        by_label = defaultdict(list)
        for name, (op, label, rr) in instances.items():
            by_label[label].append((name, op, rr))
        for edge in explicit:
            matched = []
            for end in edge['ends']:
                matched.append([(name, op.name if op else None,
                                 [i for i in chain_items[name]
                                  if any(observed_endpoint(end, r, author) for author, r in items[i][5])])
                                for name, op, rows in by_label.get(end['label_asym'], ())
                                if any(observed_endpoint(end, r, author) for author, r in rows)])
            conventional_symmetry = all(end['symmetry'] in (None, False, '1_555') for end in edge['ends'])
            for side, options in enumerate(matched):
                for _, operator, indexes in options:
                    partners = [(name, ids) for name, other_op, ids in matched[1 - side] if operator == other_op]
                    if indexes and (len(partners) != 1 or not conventional_symmetry):
                        unresolved.extend((i, edge['id']) for i in indexes)
            if not conventional_symmetry:
                continue
            for left, operator, left_items in matched[0]:
                others = [(right, right_items) for right, op, right_items in matched[1] if op == operator]
                if len(others) != 1:
                    continue
                right, right_items = others[0]
                if left_items or right_items:
                    present_links.append((set(left_items + right_items), edge, [left, right]))
                if (edge['kind'] in {'covale', 'disulf'} and len(left_items) == len(right_items) == 1
                    and all(items[i][1] != 'polymer_modification' for i in left_items + right_items)):
                    if not preserve_protein_context or not any(
                        items[i][1] == 'polymer' and items[i][0] in proteins
                        for i in left_items + right_items
                    ):
                        parent[find(right_items[0])] = find(left_items[0])

        groups = defaultdict(list)
        for i in range(len(items)):
            groups[find(i)].append(i)
        trees = ({p: cKDTree(pxyz) for p in eligible_proteins if len(pxyz := points(instances[p][0], instances[p][2]))}
                 if calculate_distances else {})
        for group in groups.values():
            members, assertions, observed_atoms = [], [], []
            protein_set = set()
            for i in group:
                name, kind, selector, op, label, residues = items[i]
                if kind == 'polymer':
                    protein_set.add(name)
                if calculate_distances:
                    coords = points(op, residues)
                    if len(coords):
                        observed_atoms.append(coords)
                selected = ('entity:' if kind in {'polymer', 'branched'} else 'ccd:') + selector
                if selected in rules:
                    assertions.append(rules[selected])
                elif kind == 'polymer_modification':
                    assertions.append(['pending', 'polymer-incorporated residue',
                                       'Source assigns this hetero residue to a polymer entity; role as an independent ligand is not established.'])
                members.append({'instance': name, 'label_asym': label, 'operator': op.name if op else None,
                                'entity': entity_by_chain[name], 'description': descriptions.get(entity_by_chain[name]),
                                'kind': kind, 'polymer_type': polymer_types.get(entity_by_chain[name]),
                                'residues': [{'auth_asym': author, 'auth_seq': str(r.seqid), 'label_seq': r.label_seq,
                                              'ccd': r.name, 'ccd_name': chemical_names.get(r.name),
                                              'positive_heavy_atoms': sum(a.occ > 0 and not a.element.is_hydrogen for a in r)}
                                             for author, r in residues]})
            if not any(r['positive_heavy_atoms'] for m in members for r in m['residues']):
                continue
            key_members = [member_key(m['instance'], m['label_asym'], m['operator'], m['entity'],
                                      r['auth_asym'], r['auth_seq'], r['ccd']) for m in members for r in m['residues']]
            links = [{'source': edge, 'instances': names} for ids, edge, names in present_links if ids.intersection(group)]
            unresolved_links = sorted({edge for i, edge in unresolved if i in group})
            positions = np.concatenate(observed_atoms) if trees else None
            distances = {p: round(float(tree.query(positions, k=1)[0].min()), 4)
                         for p, tree in trees.items() if p not in protein_set}
            contacts = sorted(p for p, distance in distances.items() if distance <= 4.)
            if assertions and len({row[0] for row in assertions}) == 1:
                role = assertions[0][0]
                subtype = '; '.join(dict.fromkeys(row[1] for row in assertions))
                rationale = ' '.join(dict.fromkeys(row[2] for row in assertions))
            else:
                role, subtype, rationale = ('pending', 'unassigned',
                    'CCD/entity and coordinates establish observation; no consistent contextual role evidence was supplied.')
            output.append({'candidate_key': key_for(entry, revision, model.num, assembly, key_members),
                           'entry': entry, 'revision': revision, 'model': model.num, 'assembly': assembly,
                           'members': members, 'source_member_identity': key_members,
                           'protein_instance_candidates': distances, 'contact_instances_4A': contacts,
                           'identity_status': 'ambiguous_or_incomplete' if (not assembly or unresolved_links or len(contacts) != 1)
                           else 'one_contact_candidate',
                           'identity_reason': ('No deposited assembly; ASU used for observations.' if not assembly else
                                              '4 A is contact evidence only; all chain candidates retained and L1 is not confirmed.'),
                           'connections': links, 'unresolved_connection_ids': unresolved_links,
                           'role': role, 'role_detail': subtype, 'role_evidence': rationale,
                           'formal_qc': 'not_run'})
    return output


def assemble_candidate_ledger(block, structure, entry, revision, systems, rules, *, direct=False):
    ledger, unmatched = {}, defaultdict(list)
    for system_id, assembly, annotations, mapped in systems:
        if not mapped or mapped['status'] != 'coordinate_identity_pass':
            unmatched[(structure[0].num, assembly)].append(system_id)
        if not mapped:
            continue
        instances = assembly_instances(structure, assembly, structure[0])
        for tag in annotations:
            instance = tag['ligand_id'].split('__')[-1]
            members = indexed_member_keys(instances, [instance])
            key = key_for(entry, revision, mapped['source_model'], assembly, members)
            proteins = mapped['protein_instances']
            l1 = key_for(entry, revision, mapped['source_model'], assembly, members +
                         [[name, instances[name][1], instances[name][0].name] for name in proteins])
            row = ledger.setdefault(key, {'candidate_key': key, 'entry': entry, 'revision': revision,
                'model': mapped['source_model'], 'assembly': assembly, 'source_member_identity': members,
                'role': 'pending', 'role_evidence': 'Index proper flag describes coverage, not biological role.',
                'formal_qc': 'not_run', 'aliases': [], 'routes': ['plinder']})
            row['aliases'].append({'system_id': system_id, 'ligand_id': tag['ligand_id'],
                                   'receptor_instances': proteins, 'status': mapped['status'],
                                   'verified_l1_id': l1 if mapped['status'] == 'coordinate_identity_pass' else None,
                                   'proper': tag['ligand_is_proper']})
    if direct or unmatched:
        native = inspect_source(block, structure, entry, revision, rules, scope=None if direct else set(unmatched))
        for row in native:
            key = row['candidate_key']
            old_row = ledger.get(key, {})
            members = {tuple(m) for m in row['source_member_identity']}
            constituent_keys = [k for k, item in ledger.items() if k != key and item['assembly'] == row['assembly']
                                and item['model'] == row['model'] and item['routes'] == ['plinder']
                                and {tuple(m) for m in item['source_member_identity']} < members]
            for k in constituent_keys:
                ledger[k]['contained_in_source_observation'] = key
            ledger[key] = old_row | row | {'routes': [*old_row.get('routes', []), 'source_mmcif'],
                'aliases': old_row.get('aliases', []), 'indexed_constituent_keys': constituent_keys,
                'fallback_systems_in_context': unmatched[(row['model'], row['assembly'])]}
    return list(ledger.values())



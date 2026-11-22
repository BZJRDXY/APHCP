"""Match observed mmCIF atoms to a PLINDER representative view."""

from collections import Counter, defaultdict
from math import dist

import gemmi


class MappingError(ValueError):
    def __init__(self, kind, message, **evidence):
        super().__init__(message)
        self.kind = kind
        self.evidence = evidence


def altloc(atom):
    return None if atom.altloc in ('\0', ' ', '.', '?') else atom.altloc


def xyz(pos):
    return (pos.x, pos.y, pos.z)


def subchains(model, *, include_water=False):
    found = defaultdict(list)
    for chain in model:
        for residue in chain:
            if include_water or not residue.is_water():
                found[residue.subchain].append((chain.name, residue))
    return found


def assembly_instances(structure, assembly_id, model, *, include_water=False):
    chains = subchains(model, include_water=include_water)
    assembly = next((a for a in structure.assemblies if a.name == assembly_id), None)
    if assembly is None:
        raise MappingError('assembly_unsupported', f'No explicit assembly {assembly_id}')
    offsets, instances = Counter(), {}
    for generator in assembly.generators:
        for operator in generator.operators:
            for label in generator.subchains:
                offsets[label] += 1
                if label in chains:
                    instances[f'{offsets[label]}.{label}'] = (operator, label, chains[label])
    return instances


def connection_edges(block):
    cat = block.get_mmcif_category('_struct_conn.')
    for i, kind in enumerate(cat.get('conn_type_id', [])):
        if kind not in {'covale', 'disulf', 'metalc'}:
            continue
        def value(key):
            return cat[key][i] if key in cat else None

        yield {'id': value('id'), 'kind': kind, 'ends': [
            {'label_asym': value(f'ptnr{j}_label_asym_id'),
             'ccd': value(f'ptnr{j}_label_comp_id'), 'atom': value(f'ptnr{j}_label_atom_id'),
             'label_seq': value(f'ptnr{j}_label_seq_id'), 'auth_asym': value(f'ptnr{j}_auth_asym_id'),
             'auth_seq': value(f'ptnr{j}_auth_seq_id'), 'insertion': value(f'pdbx_ptnr{j}_PDB_ins_code'),
             'altloc': value(f'pdbx_ptnr{j}_label_alt_id'), 'symmetry': value(f'ptnr{j}_symmetry')}
            for j in (1, 2)
        ]}


def observed_endpoint(end, residue, author=None):
    return (residue.subchain == end['label_asym'] and residue.name == end['ccd']
            and (not end['auth_asym'] or author == end['auth_asym'])
            and (not end['label_seq'] or str(residue.label_seq) == end['label_seq'])
            and (not end['auth_seq'] or str(residue.seqid.num) == end['auth_seq'])
            and (not end['insertion'] or residue.seqid.icode == end['insertion'])
            and any(a.name == end['atom'] and a.occ > 0 and (not end['altloc'] or altloc(a) == end['altloc'])
                    for a in residue))


def numeric_precision(source_block, output_block):
    """Use export precision and a float32 transform error bound, in Angstroms."""
    def decimals(token):
        mantissa, _, exp = token.lower().partition('e')
        return len(mantissa.partition('.')[2]) - int(exp or '0')

    values = (str(v) for axis in 'xyz' for v in output_block.find_values(f'_atom_site.Cartn_{axis}'))
    digits = max((decimals(value) for value in values if value not in {'.', '?'}), default=3)
    # Shorter decimal strings can be exact values with trailing zeros stripped.
    source_numbers = (abs(float(v)) for axis in 'xyz' for v in source_block.find_values(f'_atom_site.Cartn_{axis}')
                      if str(v) not in {'.', '?'})
    output_numbers = (abs(float(v)) for axis in 'xyz' for v in output_block.find_values(f'_atom_site.Cartn_{axis}')
                      if str(v) not in {'.', '?'})
    extent = max(1., max(source_numbers, default=0.), max(output_numbers, default=0.))
    epsilon = 3 ** .5 * (.5 * 10 ** -max(digits, 3) + 8 * 2 ** -24 * extent)
    return {'output_decimal_digits': digits, 'extent_A': extent, 'distance_tolerance_A': epsilon}


def map_atoms(original, output, operator, tolerance):
    sources, used, missing = {}, {}, []
    for atom in original:
        if atom.element.is_hydrogen:
            continue
        key = (atom.name, altloc(atom))
        if key in sources:
            raise MappingError('ambiguous_atom', 'Duplicate source atom identifier', atom=key)
        sources[key] = (xyz(operator.transform.apply(atom.pos)), round(atom.occ, 3), atom.element.name)
    for atom in output:
        if atom.element.is_hydrogen:
            continue
        alt = altloc(atom)
        trials = [(dist(coords, xyz(atom.pos)), key) for key, (coords, _, element) in sources.items()
                  if key[0] == atom.name and (alt is None or key[1] == alt) and element == atom.element.name]
        matches = [(distance, key) for distance, key in trials if distance <= tolerance]
        if len(matches) != 1 or matches[0][1] in used:
            raise MappingError('ambiguous_atom' if len(matches) > 1 or matches else 'atom_coordinate_mismatch',
                               f'Cannot uniquely map {output.subchain} {output.seqid} {atom.name}',
                               nearest=[(delta, key) for delta, key in sorted(trials, key=lambda t: t[0])[:4]],
                               atom=atom.name, tolerance_A=tolerance)
        distance, key = matches[0]
        used[key] = (round(atom.occ, 3), alt, distance)
    if not used:
        raise MappingError('no_observed_atom', 'No mapped heavy atom in output residue')
    for key, (_, occ, _) in sources.items():
        if occ <= 0 or (key in used and used[key][0] > 0):
            continue
        alternate_present = key[1] is not None and any(
            other[0] == key[0] and other[1] not in (None, key[1]) and sources[other][1] > 0 and occ_out > 0
            for other, (occ_out, _, _) in used.items())
        if not alternate_present:
            missing.append({'atom': key[0], 'altloc': key[1], 'occupancy': occ})
    return sources, used, missing


def map_residues(original, output, operator, tolerance):
    by_component = defaultdict(list)
    for i, (_, res) in enumerate(original):
        by_component[res.name].append(i)
    used = set()
    for _, target in output:
        indexes = by_component[target.name]
        # Branch residue numbers in reconstructed CIFs are synthesized by the
        # writer; only polymer label_seq_id represents a source sequence key.
        is_polymer = bool(indexes) and original[indexes[0]][1].entity_type == gemmi.EntityType.Polymer
        if is_polymer and target.label_seq is not None:
            indexes = [i for i in indexes if original[i][1].label_seq == target.label_seq]
        else:
            listed = [i for i in indexes if str(target.seqid) in
                      (str(original[i][1].seqid), str(original[i][1].label_seq))]
            # Preserve insertion/numbering evidence but let the atoms decide;
            # PLINDER may serialize branched numbering into a new sequence.
            indexes = listed + [i for i in indexes if i not in listed]

        candidates, failures = [], []
        anchor = next((a for a in target if not a.element.is_hydrogen), None)
        for i in indexes:
            res = original[i][1]
            if anchor is None or not any(a.name == anchor.name and a.element == anchor.element
                                         and dist(xyz(operator.transform.apply(a.pos)), xyz(anchor.pos)) <= tolerance
                                         for a in res):
                continue
            try:
                atoms, matched, missing = map_atoms(res, target, operator, tolerance)
            except MappingError as exc:
                failures.append({'source_auth_seq': str(res.seqid), 'source_label_seq': res.label_seq,
                                 'kind': exc.kind, 'detail': exc.evidence})
                continue
            candidates.append((i, atoms, matched, missing))
        if len(candidates) != 1 or candidates[0][0] in used:
            raise MappingError('residue_identity_ambiguous' if candidates or failures else 'residue_identity_unmatched',
                               f'Cannot uniquely map {target.subchain} {target.seqid} {target.name}',
                               output_label_seq=target.label_seq,
                               sources=[(str(original[i][1].seqid), original[i][1].label_seq) for i in indexes[:20]],
                               atom_failures=failures, matches=[row[0] for row in candidates])
        i, atoms, matched, missing = candidates[0]
        used.add(i)
        yield original[i], target, atoms, matched, missing
    if len(used) != len(original):
        raise MappingError('source_residues_omitted', 'Source instance has unmapped residues', omitted=[
            {'auth_seq': str(res.seqid), 'label_seq': res.label_seq, 'ccd': res.name,
             'observed_heavy_atoms': sum(not a.element.is_hydrogen and a.occ > 0 for a in res)}
            for i, (_, res) in enumerate(original) if i not in used])


def verify_system(structure, block, derived_block, assembly, proteins, ligands, annotations):
    out = gemmi.make_structure_from_block(derived_block)
    if len(out) != 1 or not out[0].count_atom_sites():
        raise MappingError('model_unsupported', 'Expected one nonempty derivative model', model_count=len(out))
    # The pinned PLINDER reconstruction calls Biotite with model=1 (first
    # model); an NMR representative is not a verification of its ensemble.
    model = structure[0]
    instances = assembly_instances(structure, assembly, model)
    observed = subchains(out[0])
    expected = set(proteins) | set(ligands)
    missing = sorted(expected - instances.keys())
    if set(observed) != expected or missing:
        raise MappingError('assembly_instance_mismatch', 'Unexpected derivative or source chain instances',
                           expected=sorted(expected), observed=sorted(observed), source_missing=missing)
    precision = numeric_precision(block, derived_block)
    tolerance = precision['distance_tolerance_A']
    edges = list(connection_edges(block))
    result = {'source_model': model.num, 'source_model_count': len(structure), 'output_model': out[0].num,
              'assembly': assembly, 'protein_instances': list(proteins), 'ligand_instances': list(ligands),
              'source_operators': {k: instances[k][0].name for k in sorted(expected)}, 'precision': precision,
              'matched_heavy_atoms': 0, 'max_distance_A': 0., 'omitted_source_altloc_atoms': 0,
              'changed_occupancies': 0, 'ligand_atom_map': {}, 'differences': [], 'target_annotations': []}
    for name in sorted(expected):
        operator, label, original = instances[name]
        try:
            for (author, res), target, atoms, matched, omitted in map_residues(original, observed[name], operator, tolerance):
                result['matched_heavy_atoms'] += len(matched)
                result['max_distance_A'] = max(result['max_distance_A'], *(row[2] for row in matched.values()))
                result['omitted_source_altloc_atoms'] += sum(key not in matched and key[1] is not None and val[1] > 0
                                                             for key, val in atoms.items())
                result['changed_occupancies'] += sum(out_occ != atoms[key][1] for key, (out_occ, _, _) in matched.items())
                if omitted:
                    result['differences'].append({'kind': 'positive_source_atoms_omitted', 'instance': name,
                                                  'auth_seq': str(res.seqid), 'component': res.name, 'atoms': omitted})
                if name in ligands:
                    result['ligand_atom_map'].setdefault(name, []).append({
                        'label_asym': label, 'auth_asym': author, 'auth_seq': str(res.seqid), 'label_seq': res.label_seq,
                        'output_seq': str(target.seqid), 'ccd': res.name,
                        'atoms': [{'atom': key[0], 'source_altloc': key[1], 'source_occupancy': val[1],
                                   'output': matched.get(key)} for key, val in atoms.items()]})
        except MappingError as exc:
            exc.evidence.update(instance=name, operator=operator.name, source_model=model.num)
            raise
    for annotation in annotations:
        name = annotation['ligand_id'].split('__')[-1]
        _, label, residues = instances[name]
        actual = Counter(res.name for _, res in residues)
        expected_components = Counter(annotation['ligand_ccd_code'].split('-'))
        links = [edge for edge in edges if any(observed_endpoint(end, res, author)
                  for author, res in residues for end in edge['ends'])]
        absent = expected_components - actual
        crosslinks = []
        for edge in links:
            if edge['kind'] != 'covale':
                continue
            for end in edge['ends']:
                if end['label_asym'] != label and end['ccd'] in absent:
                    partner = [key for key, (op, sub, members) in instances.items()
                               if sub == end['label_asym'] and op.name == instances[name][0].name
                               and any(observed_endpoint(end, r, author) for author, r in members)]
                    crosslinks.append({'edge': edge, 'observed_partner_instances': partner,
                                       'present_in_derivative': any(k in observed for k in partner)})
        if actual != expected_components:
            traceable = {end['ccd'] for link in crosslinks if link['observed_partner_instances']
                         for end in link['edge']['ends'] if end['label_asym'] != label}
            proven = not (actual - expected_components) and set(absent) <= traceable
            omitted = proven and any(not link['present_in_derivative'] for link in crosslinks)
            if not omitted:
                raise MappingError('component_identity_unresolved', f'Unexplained component identity for {name}',
                                   expected=dict(expected_components), within_instance=dict(actual), crosslinks=crosslinks)
            result['differences'].append({'kind': 'linked_source_component_omitted', 'instance': name,
                                          'missing_components': dict(absent), 'crosslinks': crosslinks})
        result['target_annotations'].append({'ligand_id': annotation['ligand_id'], 'ccd': annotation['ligand_ccd_code'],
                                             'source_components': dict(actual), 'source_explicit_links': links})
    if len(structure) > 1:
        alternatives = []
        for other in list(structure)[1:]:
            other_instances = assembly_instances(structure, assembly, other)
            try:
                for name in ligands:
                    op, _, residues = other_instances[name]
                    list(map_residues(residues, observed[name], op, tolerance))
            except (KeyError, MappingError):
                continue
            alternatives.append(other.num)
        result['model_evidence'] = {'selected_by_writer': 'PLINDER/Biotite first deposited model',
                                    'other_models_tested': len(structure) - 1,
                                    'other_models_with_same_ligand_coordinates': alternatives}
    result['status'] = 'confirmed_derivative_difference' if result['differences'] else 'coordinate_identity_pass'
    return result

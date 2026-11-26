"""Provisional §05 decisions for a single experimental receptor chain."""
from collections import defaultdict
import re

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from apo_c_pair import align_coordinates
from apo_c_triplet_regression import positive_heavy, residue_index
from holo_source_mapping import assembly_instances, connection_edges


RULE = 'single-chain-observed-site-v2'
STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def text_value(block, tag):
    return gemmi.cif.as_string(block.find_value(tag) or '')


def deposited_target_role(holo, block, structure):
    historic = holo.get('role')
    if historic == 'PASS':
        return {'status': 'HISTORICAL_EVIDENCE', 'historic_role': historic,
                'evidence': holo.get('role_evidence')}
    members = holo['target_members']
    if historic != 'REVIEW' or len(members) != 1 or holo['external_covalent_connections']:
        return {'status': 'UNRESOLVED', 'historic_role': historic}
    member = members[0]
    ccd = member['component'].upper()
    title = text_value(block, '_struct.title')
    instances = assembly_instances(structure, str(holo['identity']['assembly']), structure[0])
    actual = [(key, author, str(r.seqid)) for key, (_, _, rows) in instances.items()
              for author, r in rows if r.name.upper() == ccd and r.entity_type == gemmi.EntityType.NonPolymer]
    if actual != [(member['instance'], member['author_chain'], member['author_seq'])]:
        return {'status': 'UNRESOLVED', 'historic_role': historic, 'candidate_target_instances': actual}
    if (len(member['observed_heavy_atoms']) < 3 or len(holo['site_residues']) < 2
            or not re.search(r'\b(?:bound|complex(?:ed)?|inhibited)\b', title, re.I)
            or re.search(r'\b(?:unbound|unliganded|ligand.free|apo|without)\b', title, re.I)):
        return {'status': 'UNRESOLVED', 'historic_role': historic}
    # Treat a coordination complex as an explicit target token only if the
    # fixed deposition connects that metal to this particular ligand instance.
    linked = {}
    for edge in holo['target_connections']:
        flags = edge['target_endpoint_membership']
        if edge['kind'] != 'metalc' or sum(flags) != 1:
            continue
        ligand, metal = (edge['ends'][i] for i in (flags.index(True), flags.index(False)))
        if ligand['atom'] not in member['observed_heavy_atom_names']:
            continue
        metal_id = metal['ccd'].upper()
        if any(r.subchain == metal['label_asym'] and r.name == metal_id
               and author == metal['auth_asym']
               and str(r.seqid.num) == str(metal['auth_seq'])
               and any(a.element.is_metal and a.occ > 0 and a.name == metal['atom'] for a in r)
               for _, _, rows in instances.values() for author, r in rows):
            linked[metal_id + ccd] = edge['id']
    terms = [ccd, *linked]
    for row in block.find('_chem_comp.', ['id', 'name']):
        if gemmi.cif.as_string(row[0]).upper() == ccd:
            terms.append(gemmi.cif.as_string(row[1]).strip())
    term = next((t for t in terms if t and re.search(r'(?<!\w)' + re.escape(t) + r'(?!\w)', title, re.I)), None)
    if not term:
        return {'status': 'UNRESOLVED', 'historic_role': historic}
    return {'status': 'FIXED_DEPOSITION_INSTANCE', 'historic_role': historic,
            'field': '_struct.title', 'text': title, 'matched_term': term,
            'observed_instance': actual[0], 'coordination_record': linked.get(term),
            'source_sha256': holo['source_sha256']}


def experiment_record(block, holo_members, holo_block):
    title = text_value(block, '_struct.title')
    growth = block.get_mmcif_category('_exptl_crystal_grow.').get('pdbx_details', [])
    target_codes = {v['component'] for v in holo_members}
    components = holo_block.get_mmcif_category('_chem_comp.')
    names = set()
    for i, code in enumerate(components.get('id', [])):
        if code in target_codes:
            for key in ('name', 'pdbx_synonyms'):
                name = str(components.get(key, [''] * len(components['id']))[i] or '')
                names.update(v.strip().lower() for v in re.split(r'[;]', name) if len(v.strip()) >= 5)
    conflicts = []
    if re.search(r'\b(?:bound|complex(?:ed)? with)\b', title, re.I) and any(
            re.search(r'(?<!\w)' + re.escape(name) + r'(?!\w)', title.lower()) for name in names):
        conflicts.append({'field': '_struct.title', 'text': title, 'reason': 'depositor title names the target'})
    for detail in filter(None, growth):
        lower = detail.lower()
        if any(re.search(r'(?<!\w)' + re.escape(name) + r'(?!\w)', lower) for name in names):
            conflicts.append({'field': '_exptl_crystal_grow.pdbx_details', 'text': detail})
        elif re.search(r'\b(?:co[- ]?crystallized with|soaked (?:with|in)|in the presence of)\b', lower):
            conflicts.append({'field': '_exptl_crystal_grow.pdbx_details', 'text': detail,
                              'reason': 'added component requires context'})
    return {'title': title, 'method': text_value(block, '_exptl.method'),
            'growth_details': growth, 'citation': block.get_mmcif_category('_citation.'),
            'native_or_unliganded_title_signal': bool(re.search(r'\b(?:unliganded|ligand[- ]free|apo[- ]?form|apo|native)\b',
                                                               title, re.I)),
            'preparation_conflicts_or_unresolved_additions': conflicts}


def reference_segments(block, accession, author_chain):
    ranges = []
    rows = block.find('_struct_ref_seq.', ['pdbx_db_accession', 'pdbx_strand_id', 'seq_align_beg', 'seq_align_end',
                                         'db_align_beg', 'db_align_end'])
    for row in rows:
        acc, chains, *positions = [gemmi.cif.as_string(v) for v in row]
        if acc == accession and author_chain in (x.strip() for x in chains.split(',')) and all(p.isdecimal() for p in positions):
            ranges.append(tuple(map(int, positions)))
    return ranges


def domain_conflict(holo, hblock, block, entity, observation):
    common = entity.get('common_uniprot', [])
    if len(common) != 1 or entity.get('site_residues_covered_by_construct_sequence') != 0:
        return None
    author = observation.get('author_chain')
    h_authors = {x['author_chain'] for x in holo['site_residues']}
    if not author or len(h_authors) != 1:
        return None
    h_segments = reference_segments(hblock, common[0], next(iter(h_authors)))
    c_segments = reference_segments(block, common[0], author)
    site = []
    for row in holo['site_residues']:
        n = row['label_seq']
        mapped = {start + n - first for first, last, start, end in h_segments
                  if first <= n <= last and last - first == end - start}
        if len(mapped) != 1:
            return None
        site.extend(mapped)
    if site and c_segments and all(not start <= n <= end for n in site for _, _, start, end in c_segments):
        return {'accession': common[0], 'holo_site_uniprot_positions': sorted(site), 'candidate_reference_segments': c_segments}
    return None


def target_geometry(holo, holo_structure):
    inst = assembly_instances(holo_structure, str(holo['identity']['assembly']), holo_structure[0])
    op, _, rows = inst[holo['binding_instance']]
    residues = residue_index(rows)
    site = {row['label_seq'] for row in holo['site_residues']}
    members = []
    altlocs = set()
    for n in site:
        if any(a['altloc'] for a in positive_heavy(residues[n][1], op)):
            altlocs.add(n)
    for member in holo['target_members']:
        blank_target = [a['xyz'] for a in member['observed_heavy_atoms'] if not a['altloc']]
        if len(blank_target) < 3:
            continue
        tree = cKDTree(blank_target)
        footprint = set()
        for n in site:
            coords = [a['xyz'] for a in positive_heavy(residues[n][1], op) if not a['altloc']]
            if coords and min(tree.query(coords)[0]) <= 4:
                footprint.add(n)
        members.append((member, tree, footprint))
    return op, residues, members, altlocs


def inspect_observation(obs, holo_geometry, structure, block):
    hop, hres, members, holo_altlocs = holo_geometry
    assembly = assembly_instances(structure, str(obs['assembly']), structure[0], include_water=True)
    op, label, rows = assembly[obs['instance']]
    cres = residue_index(rows)
    coordinates, alpha_h, alpha_c = defaultdict(list), [], []
    uncertain_altlocs = []
    for row in obs['mapped_site']:
        q, n = row['holo_label_seq'], row['candidate_label_seq']
        if n not in cres:
            continue
        residue = cres[n][1]
        atoms = positive_heavy(residue, op)
        if any(a['altloc'] for a in atoms):
            uncertain_altlocs.append(q)
        for a in atoms:
            coordinates[q].append(a)
        hca = next((a for a in positive_heavy(hres[q][1], hop) if a['name'] == 'CA' and not a['altloc']), None)
        cca = next((a for a in atoms if a['name'] == 'CA' and not a['altloc']), None)
        hp, cp = (hca['xyz'] if hca else None), (cca['xyz'] if cca else None)
        if hp is not None and cp is not None:
            alpha_h.append(hp)
            alpha_c.append(cp)
    fit = None
    if len(alpha_h) >= 4 and np.linalg.matrix_rank(np.array(alpha_h) - np.mean(alpha_h, axis=0)) >= 2:
        fit = align_coordinates(np.array(alpha_h), np.array(alpha_c))
    linked_site = []
    positions = {str(r['author_seq']) for r in obs['mapped_site'] if 'author_seq' in r}
    label_positions = {str(r['candidate_label_seq']) for r in obs['mapped_site'] if r['candidate_label_seq'] is not None}
    for edge in connection_edges(block):
        if edge['kind'] == 'covale' and any(end['label_asym'] == label and
                (str(end['label_seq']) in label_positions or str(end['auth_seq']) +
                 str(end.get('insertion') or '').replace('?', '').replace('.', '') in positions)
                for end in edge['ends']):
            linked_site.append(edge)
    contacts = obs['actual_component_contacts'].get('nonwater', [])
    near = [c for c in contacts if c['atoms_within_4A_of_mapped_site'] > 0]
    ambiguous_ligands = []
    site_tree = cKDTree([a['xyz'] for values in coordinates.values() for a in values])
    for instance, (comp_op, _, residues) in assembly.items():
        if instance == obs['instance']:
            continue
        for author, residue in residues:
            if residue.entity_type == gemmi.EntityType.Polymer or residue.is_water():
                continue
            zero = []
            for atom in residue:
                if atom.occ == 0 and not atom.element.is_hydrogen:
                    p = comp_op.transform.apply(atom.pos)
                    zero.append([p.x, p.y, p.z])
            if zero and min(site_tree.query(zero)[0]) <= 4:
                ambiguous_ligands.append({'instance': instance, 'component': residue.name, 'author_chain': author,
                                          'author_seq': str(residue.seqid)})
    assignments, unresolved, companions = [], [], []
    # Collect neighboring protein residues by assembly copy for a footprint
    # check against the real atoms of the complete partner.
    protein_contacts = defaultdict(list)
    for contact in near:
        if contact['kind'] == 'polymer':
            protein_contacts[contact['instance']].append(contact)
    projected = [(member, cKDTree(tree.data @ fit[0] + fit[1]), footprint)
                 for member, tree, footprint in members] if fit and fit[2] <= 2 else []
    for instance, contacts_in_copy in protein_contacts.items():
        c_op, _, c_rows = assembly[instance]
        keys = {(x['author_chain'], x['author_seq']) for x in contacts_in_copy}
        atoms = [a['xyz'] for author, r in c_rows if (author, str(r.seqid)) in keys
                 for a in positive_heavy(r, c_op) if not a['altloc']]
        if not atoms:
            continue
        real_tree = cKDTree(atoms)
        touched = {q for q, aa in coordinates.items() if
                   (blank := [a['xyz'] for a in aa if not a['altloc']]) and min(real_tree.query(blank)[0]) <= 4}
        for member, probe, footprint in projected:
            shared = sorted(touched & footprint)
            overlap = int(np.count_nonzero(probe.query(atoms)[0] <= 2))
            if len(shared) >= 2 and overlap >= 3:
                assignments.append({'kind': 'POLYMER_PARTNER', 'component_instance': instance,
                                    'component': 'protein_partner', 'holo_member': member['component'],
                                    'residues': sorted(keys), 'shared_actual_contact_positions': shared,
                                    'candidate_heavy_atoms_within_2A_of_probe': overlap})
    for contact in near:
        if contact['kind'] == 'polymer':
            if not any(x['component_instance'] == contact['instance'] for x in assignments):
                unresolved.append(contact)
            continue
        c_op, _, c_rows = assembly[contact['instance']]
        partners = [res for author, res in c_rows if author == contact['author_chain'] and
                    str(res.seqid) == contact['author_seq'] and res.name == contact['component']]
        if len(partners) != 1:
            unresolved.append(contact)
            continue
        atoms = [a for a in positive_heavy(partners[0], c_op) if not a['altloc']]
        from apo_c_triplet_regression import ccd_heavy_atoms
        expected = ccd_heavy_atoms(block, contact['component'])
        if (len(atoms) == 1 and expected == [atoms[0]['name']]
                and fit and fit[2] <= 1 and projected
                and all(m['missing_ccd_heavy_atom_names'] == [] for m, _, _ in projected)):
            distance = min(float(probe.query(atoms[0]['xyz'])[0]) for _, probe, _ in projected)
            if distance > 4:
                companions.append(contact | {'basis': 'DEPOSITED_MONATOMIC_COMPONENT_OUTSIDE_OBSERVED_TARGET_FOOTPRINT',
                                            'distance_to_locally_aligned_target_A': round(distance, 4)})
                continue
        if len(atoms) < 3:
            unresolved.append(contact)
            continue
        xyz = np.array([a['xyz'] for a in atoms])
        tree = cKDTree(xyz)
        touched = {q for q, aa in coordinates.items() if
                   (blank := [a['xyz'] for a in aa if not a['altloc']]) and min(tree.query(blank)[0]) <= 4}
        before = len(assignments)
        for member, target_tree, footprint in members:
            shared = sorted(touched & footprint)
            if len(shared) < 2:
                continue
            overlap = 0
            if fit and fit[2] <= 2:
                rot, shift, _ = fit
                predicted = target_tree.data @ rot + shift
                overlap = int(np.count_nonzero(cKDTree(predicted).query(xyz)[0] <= 2))
            if contact['component'] == member['component'] or overlap >= 3:
                assignments.append({'kind': 'TARGET' if contact['component'] == member['component'] else 'OTHER',
                                    'component': contact['component'], 'component_instance': contact['instance'],
                                    'component_author_chain': contact['author_chain'],
                                    'component_author_seq': contact['author_seq'], 'holo_member': member['component'],
                                    'shared_actual_contact_positions': shared,
                                    'candidate_heavy_atoms_within_2A_of_probe': overlap,
                                    'local_alignment_ca_count': len(alpha_h),
                                    'compatible_atom_scope': 'positive occupancy, blank altloc on measured atoms',
                                    'local_ca_rmsd_A': round(fit[2], 4) if fit else None})
        if len(assignments) == before:
            unresolved.append(contact)
    return {'operator': str(op.name), 'author_chain': next((row['author_chain'] for row in obs['mapped_site']
                                                                      if row.get('author_chain')), None),
            'local_alignment_ca_count': len(alpha_h), 'local_ca_rmsd_A': round(fit[2], 4) if fit else None,
            'near_nonwater_instances': near, 'occupied_subsites': assignments, 'covalent_site_links': linked_site,
            'noncompetitive_companions': companions, 'unresolved_nonwater_contacts': unresolved,
            'zero_occupancy_nearby_components': ambiguous_ligands,
            'holo_site_altloc_positions': sorted(holo_altlocs),
            'candidate_site_altloc_positions': uncertain_altlocs}


def decide_candidate(holo, prepared, evidence):
    if evidence['status'] != 'SOURCE_AND_ENTITY_READ':
        technical = evidence['status'] not in {'SOURCE_ENTITY_NOT_RESOLVED', 'SOURCE_REVISION_CONFLICT_REQUIRES_REVIEW'}
        return {'rule': RULE, 'status': 'TECHNICAL_GAP' if technical else 'REVIEW', 'reasons': [evidence['status']]}
    hblock, hs, _, _, _ = prepared
    block = gemmi.cif.read(evidence['candidate_source']['path']).sole_block()
    structure = gemmi.make_structure_from_block(block)
    experimental = experiment_record(block, holo['target_members'], hblock)
    hgeom = target_geometry(holo, hs)
    binding_role = deposited_target_role(holo, hblock, hs)
    members = holo['target_members']
    ordinary_target = (len(members) == 1 and members[0]['entity_type'] == 'EntityType.NonPolymer'
                       and len(members[0]['observed_heavy_atoms']) >= 3
                       and not holo['external_covalent_connections'])
    complete_target = ordinary_target and members[0]['missing_ccd_heavy_atom_names'] == []
    target_alt = any(a['altloc'] for m in members for a in m['observed_heavy_atoms'])
    observations = []
    for entity in evidence.get('entities', []):
        obs_rows = entity.get('assembly_observations', [])
        for obs in obs_rows:
            if obs['status'] != 'ACTUAL_CONTACTS_MEASURED':
                conflict = domain_conflict(holo, hblock, block, entity, obs)
                if conflict:
                    obs['decision'] = {'status': 'SCIENTIFIC_EXCLUSION', 'reasons': ['NO_CANONICAL_TARGET_DOMAIN_COVERAGE'],
                                       'domain_conflict': conflict}
                else:
                    obs['decision'] = {'status': 'TECHNICAL_GAP' if obs['status'] == 'TECHNICAL_FAILURE' else 'REVIEW',
                                       'reasons': [obs['status']],
                                       'unknown_holo_positions': [x['holo_label_seq'] for x in obs.get('mapped_site', [])]}
                observations.append(obs)
                continue
            identity = len(entity.get('common_uniprot', [])) == 1
            if 'sifts_chain_mappings' in entity:
                identity = identity and any(m['struct_asym_id'] == obs['label_asym'] for rows in
                                           entity['sifts_chain_mappings'].values() for m in rows)
            facts = inspect_observation(obs, hgeom, structure, block)
            missing = [row['holo_label_seq'] for row in obs['mapped_site']
                       if not row.get('observed_positive_heavy_atoms') or
                       not row.get('same_residue_component') or row.get('missing_relative_holo_heavy_atoms')]
            reasons = []
            if not identity:
                reasons.append('IDENTITY_NOT_RESOLVED_FOR_CHAIN')
            if not complete_target:
                reasons.append('HOLO_TARGET_PARTIAL_OR_COMPLEX')
            if binding_role['status'] == 'UNRESOLVED' or not all(holo.get('qc', {}).get(k) == 'PASS' for k in
                    ('coordinate_identity', 'experimental_fidelity', 'topology', 'partner_scope_qc')):
                reasons.append('HOLO_ROLE_OR_SCOPE_QC_PENDING')
            if missing:
                reasons.append('MAPPED_SITE_NOT_FULLY_OBSERVED')
            al = entity['sequence_alignment']
            if (al['query_length'] != al['candidate_length'] or al['query_coverage'] != 1
                    or al['identity_on_aligned_positions'] != 1):
                reasons.append('CONSTRUCT_DIFF_REQUIRES_REVIEW')
            # When the Holo target itself has a single shared coordinate set, testing
            # the union of its receptor contacts only enlarges the required observed site.
            if target_alt or facts['candidate_site_altloc_positions']:
                reasons.append('ALTERNATE_CONFORMATIONS_REQUIRE_REVIEW')
            if facts['covalent_site_links'] or any(r.get('component') not in STANDARD for r in obs['mapped_site']):
                reasons.append('SITE_MODIFICATION_OR_COVALENT_CONTEXT')
            if facts['unresolved_nonwater_contacts'] or (facts['noncompetitive_companions'] and not
                    experimental['native_or_unliganded_title_signal']):
                reasons.append('NONWATER_CONTACT_IN_MAPPED_SITE')
            if facts['zero_occupancy_nearby_components']:
                reasons.append('DEPOSITED_NEARBY_COMPONENT_HAS_ZERO_OCCUPANCY')
            if experimental['preparation_conflicts_or_unresolved_additions']:
                reasons.append('EXPERIMENT_PREPARATION_CONTEXT_REQUIRES_REVIEW')
            if experimental['method'] != 'X-RAY DIFFRACTION' or len(structure) != 1 or len(hs) != 1:
                reasons.append('METHOD_OR_MODEL_OUTSIDE_ORDINARY_SCOPE')
            if ordinary_target and identity and facts['occupied_subsites']:
                status = 'SCIENTIFIC_EXCLUSION'
                reasons.insert(0, 'OCCUPIED_MAPPED_SUBSITE')
            else:
                status = 'REVIEW' if reasons else 'SUPPORTED_RESEARCH_PAIR'
            obs['decision'] = {'status': status, 'reasons': reasons or ['OBSERVED_MAPPED_SITE_WITHOUT_RELATED_MODELLED_OCCUPANCY'],
                               'unknown_holo_positions': missing, **facts}
            observations.append(obs)
    viable = [o for o in observations if o['decision']['status'] == 'SUPPORTED_RESEARCH_PAIR']
    if viable and len(observations) != 1:
        for obs in viable:
            obs['decision'].update(status='REVIEW', reasons=['CANDIDATE_ASSEMBLY_OR_CHAIN_NOT_UNIQUE'])
    states = [o['decision']['status'] for o in observations]
    if 'SUPPORTED_RESEARCH_PAIR' in states:
        status = 'SUPPORTED_RESEARCH_PAIR'
    elif states and all(s == 'SCIENTIFIC_EXCLUSION' for s in states):
        status = 'SCIENTIFIC_EXCLUSION'
    elif states and all(s == 'TECHNICAL_GAP' for s in states):
        status = 'TECHNICAL_GAP'
    else:
        status = 'REVIEW'
    return {'rule': RULE, 'status': status, 'candidate_experiment': experimental, 'holo_target_role': binding_role,
            'observation_states': states,
            'scope': 'Specified positive-occupancy mapped site and deposited experimental assembly/model. '
                     'Supported research pairs require complete ordinary Holo target and a unique candidate chain observation.'}

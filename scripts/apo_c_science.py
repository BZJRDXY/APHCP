"""Provisional §05 decisions for a single experimental receptor chain."""
from collections import defaultdict
import re

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from apo_c_literature import deposited_apo
from apo_c_pair import align_coordinates, reference_segments
from apo_c_triplet_regression import positive_heavy, residue_index
from holo_source_mapping import assembly_instances, connection_edges, observed_endpoint


RULE = 'single-chain-observed-site-v5'
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


def primary_experiment(block):
    cat = block.get_mmcif_category('_citation.')
    for i, key in enumerate(cat.get('id', [])):
        if key == 'primary':
            return {k: cat[k][i] for k in ('title', 'pdbx_database_id_DOI', 'pdbx_database_id_PubMed') if k in cat}
    return {}


def experiment_record(block, holo_members, holo_block, binding_role):
    title = text_value(block, '_struct.title')
    growth = block.get_mmcif_category('_exptl_crystal_grow.').get('pdbx_details', [])
    target_codes = {v['component'] for v in holo_members}
    components = holo_block.get_mmcif_category('_chem_comp.')
    names = {code.lower() for code in target_codes}
    for i, code in enumerate(components.get('id', [])):
        if code in target_codes:
            for key in ('name', 'pdbx_synonyms'):
                name = str(components.get(key, [''] * len(components['id']))[i] or '')
                names.update(v.strip().lower() for v in re.split(r'[;]', name) if len(v.strip()) >= 5)
    terms = '|'.join(re.escape(v) for v in sorted(names, key=lambda name: (-len(name), name)))
    target = r'(?<!\w)(?:' + terms + r')(?!\w)'
    conflicts, explicit = [], []
    if re.search(r'\b(?:bound|complex(?:ed)? with)\b', title, re.I) and re.search(target, title, re.I):
        conflicts.append({'field': '_struct.title', 'text': title, 'reason': 'depositor title names the target'})
    for detail in filter(None, growth):
        lower = detail.lower()
        if re.search(r'\b(?:without|absence of|no added)\s+' + target, lower):
            explicit.append({'field': '_exptl_crystal_grow.pdbx_details', 'text': detail,
                             'basis': 'TARGET_EXPLICITLY_ABSENT_IN_DEPOSITION'})
        elif re.search(target, lower):
            conflicts.append({'field': '_exptl_crystal_grow.pdbx_details', 'text': detail})
        elif re.search(r'\b(?:co[- ]?crystallized with|soaked (?:with|in)|in the presence of)\b', lower):
            conflicts.append({'field': '_exptl_crystal_grow.pdbx_details', 'text': detail,
                              'reason': 'added component requires context'})
    apo_title = deposited_apo(title)
    cite, holo_cite = primary_experiment(block), primary_experiment(holo_block)
    common_citation = any(cite.get(key) and holo_cite.get(key) and
                          str(cite[key]).casefold() == str(holo_cite[key]).casefold()
                          for key in ('pdbx_database_id_DOI', 'pdbx_database_id_PubMed'))
    study = cite.get('title') or ''
    aliases = [binding_role.get('matched_term'), *sorted(names, key=lambda name: (-len(name), name))]
    bound_term = next((term for term in aliases if term
        and re.search(r'(?<!\w)' + re.escape(term) + r'[- ]bound\b', study, re.I)), None)
    paired_study = (apo_title and common_citation and bound_term
                    and study.casefold() == (holo_cite.get('title') or '').casefold()
                    and re.search(r'\b(?:apo|unliganded|ligand[- ]free)\b', study, re.I))
    if paired_study:
        explicit.append({'field': '_struct.title and _citation.title', 'basis': 'APO_ARM_OF_MATCHED_PRIMARY_EXPERIMENT',
                         'deposition_title': title, 'paired_primary_citation': cite, 'bound_target_term': bound_term})
    status = ('CONFLICT' if conflicts else 'EXPLICIT_UNBOUND' if explicit else
              'DEPOSITION_CONTEXT' if apo_title or re.search(r'\bnative\b', title, re.I) else 'NOT_ESTABLISHED')
    return {'status': status, 'evidence': explicit, 'title': title, 'method': text_value(block, '_exptl.method'),
            'growth_details': growth, 'citation': block.get_mmcif_category('_citation.'),
            'primary_citation': cite, 'holo_primary_citation': holo_cite, 'source_entry': text_value(block, '_entry.id'),
            'preparation_conflicts_or_unresolved_additions': conflicts}


def site_metal_coordination(edges, assembly, instance, residues):
    """Resolve deposited coordination on the selected model and assembly copy."""
    operator, label, _ = assembly[instance]

    def residue_matches(end, author, residue):
        return (residue.subchain == end['label_asym'] and residue.name == end['ccd']
                and (not end['auth_asym'] or author == end['auth_asym'])
                and (not end['label_seq'] or str(residue.label_seq) == str(end['label_seq']))
                and (not end['auth_seq'] or str(residue.seqid.num) == str(end['auth_seq']))
                and residue.seqid.icode.strip() == str(end.get('insertion') or '').strip())

    def endpoint(end, key, op, author, residue):
        if not observed_endpoint(end, residue, author):
            return []
        return [{'instance': key, 'label_asym': residue.subchain, 'author_chain': author,
                 'author_seq': str(residue.seqid), 'label_seq': residue.label_seq, 'component': residue.name,
                 'atom': a.name, 'element': a.element.name, 'is_metal': a.element.is_metal,
                 'altloc': a.altloc if a.altloc != '\0' else '', 'occupancy': a.occ,
                 'xyz': [p.x, p.y, p.z]} for a in residue
                if a.name == end['atom'] and a.occ > 0 and
                (not end['altloc'] or a.altloc == end['altloc'])
                for p in [op.transform.apply(a.pos)]]

    records = []
    for edge in edges:
        if edge['kind'] != 'metalc':
            continue
        for i, site_end in enumerate(edge['ends']):
            if site_end['label_asym'] != label:
                continue
            matched = [(author, r) for author, r in residues if residue_matches(site_end, author, r)]
            if not matched:
                continue
            partner = edge['ends'][1-i]
            site_atoms = [a for author, r in matched for a in endpoint(site_end, instance, operator, author, r)]
            partner_atoms = [a for key, (op, subchain, rows) in assembly.items()
                             if subchain == partner['label_asym'] and operator.transform.approx(op.transform, 1e-6)
                             for author, r in rows if residue_matches(partner, author, r)
                             for a in endpoint(partner, key, op, author, r)]
            resolved = (len(site_atoms) == len(partner_atoms) == 1
                        and all(e['symmetry'] in (None, '1_555') for e in edge['ends'])
                        and (site_atoms[0]['is_metal'] or partner_atoms[0]['is_metal']))
            distance = float(np.linalg.norm(np.subtract(site_atoms[0]['xyz'], partner_atoms[0]['xyz']))) if resolved else None
            if resolved and (any(a['altloc'] for a in site_atoms + partner_atoms) or not 0 < distance <= 4):
                resolved = False
            records.append({'connection': edge, 'site_endpoint': i,
                            'site_atoms': site_atoms, 'partner_atoms': partner_atoms,
                            'status': 'COORDINATION_OBSERVED' if resolved else 'DECLARED_COORDINATION_UNRESOLVED',
                            'observed_distance_A': round(distance, 4) if distance is not None else None})
    return records


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


def complete_altlocs(atoms, required):
    """Check each measured local conformation together with its shared atoms."""
    locations = {a['altloc'] for a in atoms if a['altloc']}
    return all(required <= {a['name'] for a in atoms if a['altloc'] in ('', alt)}
               for alt in locations)


def end_segments_only(positions, length):
    if not positions:
        return True
    shared = set(range(1, length + 1)).difference(positions)
    if not shared:
        return False
    first, last = min(shared), max(shared)
    return all(p < first or p > last for p in positions)


def compare_construct(alignment, hres, hop, holo_site, cres, op, site_coords):
    left = alignment['unaligned_holo_positions']
    right = alignment['unaligned_candidate_positions']
    if not left and not right and alignment['identity_on_aligned_positions'] == 1:
        return {'status': 'IDENTICAL'}
    record = {'status': 'UNRESOLVED', 'unaligned_holo_positions': left,
              'unaligned_candidate_positions': right}
    if (alignment['identity_on_aligned_positions'] != 1
            or not end_segments_only(left, alignment['query_length'])
            or not end_segments_only(right, alignment['candidate_length'])
            or set(left).intersection(holo_site)):
        return record
    holo_coords = [a['xyz'] for p in holo_site for a in positive_heavy(hres[p][1], hop)]
    for positions, index, operator, site in ((left, hres, hop, holo_coords), (right, cres, op, site_coords)):
        extension = [a['xyz'] for p in positions if p in index for a in positive_heavy(index[p][1], operator)]
        if extension and min(cKDTree(site).query(extension)[0]) <= 4:
            return record | {'basis': 'TERMINAL_RESIDUES_CONTACT_THE_OBSERVED_SITE'}
    return record | {'status': 'TERMINAL_ONLY_OUTSIDE_OBSERVED_SITE'}


def inspect_observation(obs, holo_geometry, structure, block, alignment):
    hop, hres, members, holo_altlocs = holo_geometry
    assembly = assembly_instances(structure, str(obs['assembly']), structure[0], include_water=True)
    op, label, rows = assembly[obs['instance']]
    cres = residue_index(rows)
    coordinates, alpha_h, alpha_c = defaultdict(list), [], []
    candidate_altlocs, incomplete_altlocs = [], []
    for row in obs['mapped_site']:
        q, n = row['holo_label_seq'], row['candidate_label_seq']
        if n not in cres:
            continue
        residue = cres[n][1]
        atoms = positive_heavy(residue, op)
        holo_atoms = positive_heavy(hres[q][1], hop)
        required = {a['name'] for a in holo_atoms}
        if not complete_altlocs(holo_atoms, required) or not complete_altlocs(atoms, required):
            incomplete_altlocs.append(q)
        if any(a['altloc'] for a in atoms):
            candidate_altlocs.append(q)
        for a in atoms:
            coordinates[q].append(a)
        hca = next((a for a in holo_atoms if a['name'] == 'CA' and not a['altloc']), None)
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
    metal_links = site_metal_coordination(
        connection_edges(block), assembly, obs['instance'],
        [cres[n] for n in {row['candidate_label_seq'] for row in obs['mapped_site']} if n in cres])
    contacts = obs['actual_component_contacts'].get('nonwater', [])
    near = [c for c in contacts if c['atoms_within_4A_of_mapped_site'] > 0]
    ambiguous_ligands = []
    site_coords = [a['xyz'] for atoms in coordinates.values() for a in atoms]
    site_tree = cKDTree(site_coords)
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
        if any(c['connection']['ends'][1-c['site_endpoint']]['label_asym'] == contact['label_asym']
               and c['connection']['ends'][1-c['site_endpoint']]['ccd'] == contact['component']
               and str(c['connection']['ends'][1-c['site_endpoint']]['auth_seq']) == str(partners[0].seqid.num)
               for c in metal_links):
            unresolved.append(contact | {'basis': 'DEPOSITED_SITE_COORDINATION'})
            continue
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
            'site_metal_coordination': metal_links,
            'noncompetitive_companions': companions, 'unresolved_nonwater_contacts': unresolved,
            'zero_occupancy_nearby_components': ambiguous_ligands,
            'holo_site_altloc_positions': sorted(holo_altlocs),
            'candidate_site_altloc_positions': candidate_altlocs,
            'incomplete_site_altloc_positions': incomplete_altlocs,
            'construct_comparison': compare_construct(alignment, hres, hop,
                [r['holo_label_seq'] for r in obs['mapped_site']], cres, op, site_coords)}


def decide_candidate(holo, prepared, evidence, literature=None):
    if evidence['status'] != 'SOURCE_AND_ENTITY_READ':
        technical = evidence['status'] not in {'SOURCE_ENTITY_NOT_RESOLVED', 'SOURCE_REVISION_CONFLICT_REQUIRES_REVIEW'}
        return {'rule': RULE, 'status': 'TECHNICAL_GAP' if technical else 'REVIEW', 'reasons': [evidence['status']],
                'candidate_experiment': {'status': 'NOT_ASSESSED'}, 'delivery_eligible': False}
    hblock, hs, _, _, _ = prepared
    block = gemmi.cif.read(evidence['candidate_source']['path']).sole_block()
    structure = gemmi.make_structure_from_block(block)
    binding_role = deposited_target_role(holo, hblock, hs)
    experimental = experiment_record(block, holo['target_members'], hblock, binding_role)
    experimental.update(fixed_revision=evidence['candidate_source']['revision'],
                        source_sha256=evidence['candidate_source']['sha256'])
    hgeom = target_geometry(holo, hs)
    models = [m.num for m in structure]
    members = holo['target_members']
    ordinary_target = (len(members) == 1 and members[0]['entity_type'] == 'EntityType.NonPolymer'
                       and len(members[0]['observed_heavy_atoms']) >= 3
                       and not holo['external_covalent_connections'])
    complete_target = ordinary_target and members[0]['missing_ccd_heavy_atom_names'] == []
    target_alt = any(a['altloc'] for m in members for a in m['observed_heavy_atoms'])
    target_alt_complete = all(complete_altlocs(m['observed_heavy_atoms'],
                                              set(m['observed_heavy_atom_names'])) for m in members)
    observations = []
    pending_experiment = []
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
                obs['decision']['site_observation'] = {'status': 'UNAVAILABLE',
                    'observed_model': obs.get('model'), 'assembly': obs.get('assembly'), 'instance': obs.get('instance')}
                obs['decision']['delivery_eligible'] = False
                observations.append(obs)
                continue
            identity = len(entity.get('common_uniprot', [])) == 1
            if 'sifts_chain_mappings' in entity:
                identity = identity and any(m['struct_asym_id'] == obs['label_asym'] for rows in
                                           entity['sifts_chain_mappings'].values() for m in rows)
            facts = inspect_observation(obs, hgeom, structure, block, entity['sequence_alignment'])
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
            if facts['construct_comparison']['status'] == 'UNRESOLVED':
                reasons.append('CONSTRUCT_DIFF_REQUIRES_REVIEW')
            has_altlocs = target_alt or facts['candidate_site_altloc_positions'] or facts['holo_site_altloc_positions']
            # The union of measured conformations can establish absence of nearby components.
            # A real contact with alternate positions needs assessment in compatible conformations.
            altloc_pending = bool(has_altlocs and (not target_alt_complete
                or facts['incomplete_site_altloc_positions'] or facts['near_nonwater_instances']))
            if altloc_pending:
                reasons.append('ALTERNATE_CONFORMATIONS_REQUIRE_REVIEW')
            if facts['covalent_site_links'] or any(r.get('component') not in STANDARD for r in obs['mapped_site']):
                reasons.append('SITE_MODIFICATION_OR_COVALENT_CONTEXT')
            if facts['site_metal_coordination']:
                reasons.append('SITE_METAL_COORDINATION_REQUIRES_REVIEW')
            if facts['unresolved_nonwater_contacts']:
                reasons.append('NONWATER_CONTACT_IN_MAPPED_SITE')
            if facts['zero_occupancy_nearby_components']:
                reasons.append('DEPOSITED_NEARBY_COMPONENT_HAS_ZERO_OCCUPANCY')
            if experimental['preparation_conflicts_or_unresolved_additions']:
                reasons.append('EXPERIMENT_PREPARATION_CONTEXT_REQUIRES_REVIEW')
            if experimental['method'] != 'X-RAY DIFFRACTION' or len(structure) != 1 or len(hs) != 1:
                reasons.append('METHOD_OR_MODEL_OUTSIDE_ORDINARY_SCOPE')
            site_status = ('OCCUPIED' if facts['occupied_subsites'] else
                           'INCOMPLETE' if missing or not complete_target else
                            'UNRESOLVED' if (not identity or altloc_pending
                               or facts['unresolved_nonwater_contacts'] or facts['site_metal_coordination']
                               or facts['covalent_site_links'] or facts['zero_occupancy_nearby_components']
                               or any(r.get('component') not in STANDARD for r in obs['mapped_site'])) else
                           'OBSERVED_NO_COMPETING_OCCUPANCY')
            if ordinary_target and identity and facts['occupied_subsites']:
                if len(models) == 1:
                    status = 'SCIENTIFIC_EXCLUSION'
                    reasons.insert(0, 'OCCUPIED_MAPPED_SUBSITE')
                else:
                    status = 'REVIEW'
                    reasons.insert(0, 'OCCUPIED_IN_INSPECTED_MODEL_OTHER_MODELS_NOT_CHECKED')
            else:
                status = 'REVIEW' if reasons else 'SUPPORTED_RESEARCH_PAIR'
            obs['decision'] = {'status': status, 'reasons': reasons,
                               'unknown_holo_positions': missing, 'delivery_eligible': status == 'SUPPORTED_RESEARCH_PAIR',
                               'experimental_state': experimental['status'],
                               'site_observation': {'status': site_status, 'observed_model': obs['model'],
                                                    'assembly': obs['assembly'], 'instance': obs['instance'],
                                                    'unexamined_models': [m for m in models if m != obs['model']]}, **facts}
            if site_status != 'OCCUPIED':
                pending_experiment.append(obs['decision'])
            observations.append(obs)
    if (literature is not None and experimental['status'] in {'DEPOSITION_CONTEXT', 'NOT_ESTABLISHED'}
            and any(d['delivery_eligible'] for d in pending_experiment)):
        article = literature(experimental)
        experimental['literature'] = article
        if article.get('conflict'):
            experimental['status'] = 'CONFLICT'
        elif article.get('support'):
            experimental['status'] = 'EXPLICIT_UNBOUND'
            experimental['evidence'].append(article['support'])
    for decision in pending_experiment:
        decision['experimental_state'] = experimental['status']
        if experimental.get('literature', {}).get('conflict'):
            decision['reasons'].append('EXPERIMENT_DOCUMENTED_ACTIVE_SITE_INHIBITION')
            decision.update(status='REVIEW', delivery_eligible=False)
        elif experimental['status'] not in {'EXPLICIT_UNBOUND', 'CONFLICT'}:
            decision['reasons'].append('EXPERIMENT_UNBOUND_NOT_ESTABLISHED')
            decision.update(status='REVIEW', delivery_eligible=False)
        elif not decision['reasons']:
            decision['reasons'].append('OBSERVED_MAPPED_SITE_WITHOUT_RELATED_MODELLED_OCCUPANCY')
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
            'delivery_eligible': status == 'SUPPORTED_RESEARCH_PAIR',
            'candidate_model_scope': {'deposited': models, 'inspected': models[:1], 'unexamined': models[1:]},
            'observation_states': states,
            'scope': 'Specified positive-occupancy mapped site and deposited experimental assembly/model. '
                     'Each supported observation has its own assembly, model and receptor chain copy.'}

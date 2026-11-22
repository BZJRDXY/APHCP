"""Rebuild the four fixed controls and assert source recovery using the shared mapper."""

import argparse
from pathlib import Path

import gemmi
from plinder.core.release import PlinderRelease

import holo_pilot as pilot
from holo_candidate_batch import atomic_json, parse_system_id, revision_from_cif, sha256
from holo_candidate_regression import json_from
from holo_source_mapping import MappingError, verify_system


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline_fixed', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('plinder_release', type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False, parents=True)
    (args.output / 'raw').symlink_to((args.baseline_fixed / 'raw').resolve(), target_is_directory=True)
    prior = {r['pdb_id']: r for r in json_from(args.baseline_fixed / 'holo_systems.json')}
    release = PlinderRelease(data_dir=args.plinder_release)
    records = []
    for case in pilot.CASES:
        pdb, revision, _, ccd, system_id, _, ligand, _ = case
        old = prior[pdb.upper()]
        source = args.output / old['source']
        assert sha256(source) == old['source_sha256']
        block, parsed_revision = revision_from_cif(source, pdb)
        assert parsed_revision == revision
        structure = gemmi.make_structure_from_block(block)
        asm, proteins, ligands = parse_system_id(system_id)
        tag = {'ligand_id': f'{pdb}__{asm}__{ligand}', 'ligand_ccd_code': ccd}
        # The fixed pilot also records retained/contact chains, all source
        # poses and the exact branched-link endpoint-to-output mapping.
        regenerated = pilot.sample(args.output, release, case)
        derivative = args.output / regenerated['derivative']
        checks = verify_system(structure, block, gemmi.cif.read(str(derivative)).sole_block(),
                               asm, proteins, ligands, [tag])
        assert checks['status'] == 'coordinate_identity_pass'
        prior_readback = verify_system(structure, block,
            gemmi.cif.read(str(args.baseline_fixed / old['derivative'])).sole_block(), asm, proteins, ligands, [tag])
        assert prior_readback['status'] == 'coordinate_identity_pass'
        details = regenerated['coordinate_checks']
        evidence = {'pdb_id': pdb.upper(), 'revision': revision, 'status': checks['status'],
                    'source_sha256': sha256(source), 'derivative_sha256': sha256(derivative),
                    'derivative': regenerated['derivative'], 'old_derivative_readback': 'pass',
                    'retained_proteins': details['retained_protein_chains'], 'ligand_instance': ligand,
                    'ligand_heavy_atoms': sum(a['output'] is not None for m in checks['ligand_atom_map'][ligand]
                                              for a in m['atoms'])}
        if pdb == '1m17':
            assert ccd == 'AQ4' and evidence['ligand_heavy_atoms'] == 29
        elif pdb == '1b86':
            assert details['retained_protein_chains'] == details['ligand_contact_chains'] == ['1.B', '1.D']
            assert len(details['assembly_member_chains']) == len(regenerated['source_hem_instances']) == 4
            wrong_tag = dict(tag, ligand_ccd_code='HEM')
            try:
                verify_system(structure, block, gemmi.cif.read(str(derivative)).sole_block(),
                              asm, proteins, ligands, [wrong_tag])
            except MappingError as exc:
                assert exc.kind == 'component_identity_unresolved'
            else:
                raise ValueError('Incorrect alternative component was accepted')
            heme_labels = regenerated['source_hem_instances']
            assert ligand.split('.', 1)[1] not in heme_labels
            wrong_ligand = f'1.{heme_labels[0]}'
            try:
                verify_system(structure, block, gemmi.cif.read(str(derivative)).sole_block(), asm,
                              proteins, [wrong_ligand], [{'ligand_id': f'{pdb}__{asm}__{wrong_ligand}',
                                                        'ligand_ccd_code': 'HEM'}])
            except MappingError as exc:
                assert exc.kind == 'assembly_instance_mismatch'
            else:
                raise ValueError('Incorrect HEM chain instance was accepted')
            evidence.update(contact_chains=details['ligand_contact_chains'], observed_hem_labels=heme_labels,
                            wrong_target_rejected=True, ligand_heavy_atoms=evidence['ligand_heavy_atoms'])
        elif pdb == '1hew':
            assert evidence['ligand_heavy_atoms'] == 43
            assert [m['ccd'] for m in checks['ligand_atom_map'][ligand]] == ['NAG'] * 3
            links = details['ligand_links_from_source']
            edges = {(l['endpoints'][0]['source_author_seq_id'], l['endpoints'][0]['atom'],
                      l['endpoints'][1]['source_author_seq_id'], l['endpoints'][1]['atom']) for l in links}
            assert edges == {('1', 'O4', '2', 'C1'), ('2', 'O4', '3', 'C1')}
            evidence['links'] = links
        elif pdb == '3tct':
            assert details['assembly_member_chains'] == ['1.A', '1.B', '2.A', '2.B']
            assert details['ligand_contact_chains'] == ['1.A', '2.A']
            atoms = checks['ligand_atom_map'][ligand][0]['atoms']
            assert all(sum(a['source_altloc'] == alt and a['source_occupancy'] == .25 for a in atoms) == 20
                       for alt in ('A', 'B'))
            assert {a['source_altloc'] for a in atoms if a['output'] is not None} == {'A'}
            assert {a['output'][0] for a in atoms if a['output'] is not None} == {1.}
            evidence.update(assembly_members=details['assembly_member_chains'], contacts=details['ligand_contact_chains'],
                            source_pose_atoms={'A': 20, 'B': 20}, source_pose_occupancy=.25,
                            output_pose='A', output_occupancy=1.)
        atomic_json(args.output / 'cases' / f'{pdb}.json', {'fixed_checks': evidence,
                    'readback': checks, 'contact_and_topology_source_proof': details})
        records.append(evidence)
    atomic_json(args.output / 'summary.json', {'status': 'PASS', 'cases': records, 'fixed_cases': len(records),
                'new_downloads': 0, 'fixed_script_sha256': sha256(Path(__file__))})


if __name__ == '__main__':
    main()

"""Conservative, instance-tied interpretation of explicit depositor role statements."""
from collections import defaultdict
import re

from holo_v1_export import rows, text

VERSION = 'm2-role-source-evidence-1'
CRYSTAL = ('_exptl_crystal_grow.pdbx_details', '_exptl_crystal_grow.details')
EXPERIMENT = (*CRYSTAL, '_exptl.details', '_diffrn.details', '_struct.pdbx_model_details')
CUE = r'(?:in complex with|complexed with|complex with|bound to|bound with|soaked with|in the presence of)'
ROLE = r'(?:substrate(?: analogue| analog| mimic)?|inhibitor|agonist|antagonist)'
# These words describe a role on the exact entity row, rather than a term elsewhere in an entry.
BOUND_SITE = re.compile(r'\bbound (?:at|to|in) (?:the |a |an )?(?:(?:allosteric|regulatory|active|catalytic|substrate|binding) site)\b')
ENV_CONFLICT = re.compile(r'\b(?:not|instead|without|rather than|apo|unbound|unliganded|control|replaced|no)\b')
ROLE_CONFLICT = re.compile(r'\b(?:not|without|replaced|absent|putative|probable|possible|potential|only|unbound|unliganded|cofactor|cryo(?: |)protect\w*|detergent|buffer|precipitant|solvent)\b')


def normalize(value):
    return ' '.join(re.findall('[a-z0-9]+', str(value or '').lower()))


def role_source(block):
    return {'entity': {text(x['id']): {k:text(x.get(k)) for k in ('type','pdbx_description','details')}
                       for x in rows(block,'_entity.')},
            'experiment': {k:[text(v) for v in block.find_values(k) if text(v)] for k in EXPERIMENT}}


def token(name):
    return r'\b'+re.escape(name)+r'\b'


def chemical_names(options):
    return sorted({normalize(n) for _,n in options if len(normalize(n))>=4}, key=lambda n:(-len(n), n))


def unique_name(name, titles, own):
    """An alias cannot point to two different deposited entities/CCD assignments."""
    return sum(name in choices for key,choices in titles.items() if key != own)==0


def build_role_context(source, names):
    aliases={cc:chemical_names(a) for cc,a in names.items()}
    entities=defaultdict(set)
    for eid,row in source['entity'].items():
        name=normalize(row.get('pdbx_description'))
        if name:
            entities[name].add(eid)
    return aliases,entities


def assess(row, prior, title, names, polylen, source, context=None):
    """Return only new source-evidence assertions; leave old positive rules separate.

    A primary-paper title is deliberately not used here: a publication may cover several
    different entries. Contact ties an assertion to an observation; it cannot create one.
    """
    if prior=='PASS' or not row['contact_instances_4A'] or len(row['members'])!=1:
        return None
    member=row['members'][0]; eid=member['entity']; kind=member['kind']; residues=member['residues']
    ent=source['entity'].get(eid,{}); details=ent.get('details',''); d=normalize(details); t=normalize(title)
    aliases,entities=context if context is not None else build_role_context(source,names)
    if kind=='component' and len(residues)==1:
        cc=residues[0]['ccd']; code=normalize(cc); exact_names=aliases.get(cc,[])
        if ent.get('type')!='non-polymer':
            return None
        # Depositor's own instance entity explicitly assigns a binding site to this chemical.
        for name in exact_names:
            if not ROLE_CONFLICT.search(d) and unique_name(name,aliases,cc) and re.match(token(name)+r' (?:is |was )?',d):
                if BOUND_SITE.search(d):
                    return 'PASS', {'rule':'entity_chemical_bound_site','field':'_entity.details','entity':eid,
                                    'ccd':cc,'matched':name,'text':details}
        # A short literal CCD is useful only if the immediately following role and context disambiguate it.
        if len(code)>=3 and not ROLE_CONFLICT.search(t) and any(re.search(r'\b'+CUE+r'\b',t[:m.start()]) for m in
                re.finditer(token(code)+r' (?:as (?:an? |the )?)?'+ROLE+r'\b',t)):
            return 'PASS', {'rule':'literal_ccd_with_explicit_title_role','field':'_struct.title','entity':eid,'ccd':cc,
                            'matched':cc,'text':title}
        # Explicit designation as cryoprotectant. A reservoir or soak mention alone is insufficient.
        if not any(re.search(token(n),t) for n in exact_names+[code]) and not d:
            for field in CRYSTAL:
                for raw in source['experiment'].get(field,[]):
                    for clause in re.split(r'(?<=[.!?])\s+|[;\n]',raw):
                        c=normalize(clause)
                        if ENV_CONFLICT.search(c):
                            continue
                        for name in exact_names:
                            if not unique_name(name,aliases,cc) or re.search(token(name)+r' \d',c):
                                continue
                            not_target = re.search(
                                token(name)+r'\s+as (?:an? |the )?cryo ?protectant\b|'
                                r'\bcryo ?protected (?:with|using) (?:a )?(?:\d+(?: \d+)? )?'+token(name)+r'(?:$| and | for | before | prior )|'
                                r'\bcryo ?protectant\s*:? (?:\d+(?: \d+)? )?'+token(name)+r'$',c)
                            if not_target and not any(re.search(token(name)+r' (?:as (?:an? |the )?)?'+ROLE+r'\b',z) for z in (d,t)):
                                return 'NON_TARGET_EVIDENCED', {'rule':'explicit_crystal_cryoprotectant','field':field,
                                        'entity':eid,'ccd':cc,'matched':name,'text':clause.strip()}
    if kind=='polymer' and str(member.get('polymer_type') or '').startswith('polypeptide'):
        if not 2<=polylen.get(eid,0)<=35 or ROLE_CONFLICT.search(t):
            return None
        # A non-generic sequence alias is extracted from the deposited entity's own prose.
        # It must name a single entity and occur as the named peptide in the complex statement.
        m=re.match(r'^(?:cyclic|macrocyclic|bicyclic) peptide ([a-z]{2,8} [a-z]{0,3}\d{1,4})(?: |$)',d)
        if m:
            alias=m[1]
            if sum(bool(re.search(token(alias),normalize(v.get('details')))) for v in source['entity'].values())==1:
                if re.search(r'\b'+CUE+r'\b.{0,90}\b(?:cyclic|macrocyclic|bicyclic) peptide '+re.escape(alias)+r'\b',t):
                    return 'PASS', {'rule':'entity_named_cyclic_peptide','field':'_entity.details + _struct.title',
                                    'entity':eid,'matched':alias,'text':details,'deposition_title':title}
        desc=normalize(ent.get('pdbx_description'))
        if entities.get(desc)=={eid} and re.fullmatch(r'(?:bisubstrate |bicyclic |cyclic |macrocyclic )?peptide inhibitor',desc):
            match=desc if 'bisubstrate ' not in desc else 'bisubstrate inhibitor'
            if not any(v['type']=='non-polymer' and 'inhibitor' in normalize(v.get('pdbx_description')) for v in source['entity'].values()) and re.search(
                    r'\b'+CUE+r'\b(?: (?:an?|the))? '+re.escape(match)+r'\b',t):
                return 'PASS', {'rule':'entity_explicit_peptide_inhibitor','field':'_entity.pdbx_description + _struct.title',
                                'entity':eid,'matched':match,'text':ent['pdbx_description'],'deposition_title':title}
    return None

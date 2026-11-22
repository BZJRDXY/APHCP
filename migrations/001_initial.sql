CREATE TABLE schema_migrations (
    version INTEGER PRIMARY KEY, filename TEXT NOT NULL, applied_at TEXT NOT NULL
);
CREATE TABLE project (
    id TEXT PRIMARY KEY, name TEXT NOT NULL, milestone TEXT NOT NULL, status TEXT NOT NULL
);
CREATE TABLE module (
    id TEXT PRIMARY KEY, name TEXT NOT NULL, status TEXT NOT NULL
);
CREATE TABLE batch (
    id TEXT PRIMARY KEY, module_id TEXT NOT NULL REFERENCES module(id), status TEXT NOT NULL,
    source_release TEXT, window_start TEXT, window_end TEXT, input_manifest TEXT,
    config_version TEXT, git_commit TEXT, candidate_count INTEGER, next_step TEXT
);
CREATE TABLE task (
    id TEXT PRIMARY KEY, batch_id TEXT NOT NULL REFERENCES batch(id), name TEXT NOT NULL,
    task_key TEXT NOT NULL UNIQUE, status TEXT NOT NULL
);
CREATE TABLE run (
    id TEXT PRIMARY KEY, batch_id TEXT NOT NULL REFERENCES batch(id),
    task_id TEXT NOT NULL REFERENCES task(id), status TEXT NOT NULL, pid INTEGER,
    git_commit TEXT NOT NULL, code_version TEXT NOT NULL,
    config_version TEXT NOT NULL, config_hash TEXT NOT NULL, source_release TEXT NOT NULL,
    input_manifest_hash TEXT NOT NULL, output_manifest_hash TEXT, environment_lock_hash TEXT NOT NULL,
    host TEXT NOT NULL, start_time_utc TEXT NOT NULL, end_time_utc TEXT, elapsed_seconds REAL,
    input_records INTEGER NOT NULL, success_records INTEGER, failed_records INTEGER,
    uncertain_records INTEGER, bytes_read INTEGER, bytes_written INTEGER, downloaded_bytes INTEGER,
    cpu_seconds REAL, peak_disk_usage INTEGER, peak_memory INTEGER, error_summary TEXT, report_path TEXT
);
CREATE TABLE pdb_entry (id TEXT PRIMARY KEY COLLATE NOCASE, title TEXT, released_on TEXT);
CREATE TABLE entry_version (
    id TEXT PRIMARY KEY, pdb_id TEXT NOT NULL COLLATE NOCASE REFERENCES pdb_entry(id),
    source_name TEXT NOT NULL, source_release TEXT NOT NULL, source_record_id TEXT NOT NULL,
    revision TEXT NOT NULL, UNIQUE (source_name, source_release, source_record_id, revision)
);
CREATE INDEX entry_version_pdb ON entry_version(pdb_id);
CREATE TABLE asset_file (
    id TEXT PRIMARY KEY, entry_version_id TEXT NOT NULL REFERENCES entry_version(id),
    uri TEXT NOT NULL, sha256 TEXT, bytes INTEGER, status TEXT NOT NULL, fetched_at TEXT,
    fetch_run_id TEXT REFERENCES run(id)
);
CREATE TABLE work_item (
    id TEXT PRIMARY KEY, entry_version_id TEXT NOT NULL REFERENCES entry_version(id),
    parent_id TEXT REFERENCES work_item(id), kind TEXT NOT NULL CHECK
        (kind IN ('assembly', 'chain', 'ligand', 'system')),
    label TEXT NOT NULL, component_id TEXT, next_step TEXT
);
CREATE TABLE system_member (
    system_id TEXT NOT NULL REFERENCES work_item(id),
    member_id TEXT NOT NULL REFERENCES work_item(id),
    role TEXT NOT NULL CHECK (role IN ('protein', 'ligand')),
    PRIMARY KEY (system_id, member_id, role)
);
CREATE TABLE run_item (
    run_id TEXT NOT NULL REFERENCES run(id), object_type TEXT NOT NULL, object_id TEXT NOT NULL,
    pdb_id TEXT COLLATE NOCASE REFERENCES pdb_entry(id), stage TEXT NOT NULL, status TEXT NOT NULL,
    error_code TEXT, error_message TEXT, output_uri TEXT, updated_at TEXT NOT NULL,
    PRIMARY KEY (run_id, object_type, object_id, stage)
);
CREATE INDEX run_item_pdb ON run_item(pdb_id, updated_at DESC);
CREATE INDEX run_item_object ON run_item(object_id, updated_at DESC);
CREATE TABLE config_version (id TEXT PRIMARY KEY, sha256 TEXT NOT NULL, path TEXT);
CREATE TABLE decision (id TEXT PRIMARY KEY, opened_on TEXT, topic TEXT, details TEXT);
CREATE TABLE threshold (
    id TEXT PRIMARY KEY, name TEXT, prior_value TEXT, value TEXT, changed_at TEXT,
    rationale TEXT, impact TEXT, config_version TEXT REFERENCES config_version(id), status TEXT, evidence TEXT
);
CREATE TABLE issue (id TEXT PRIMARY KEY, stage TEXT, status TEXT, description TEXT, next_evidence TEXT);
CREATE TABLE storage_snapshot (sampled_at TEXT PRIMARY KEY, report_json TEXT NOT NULL);
CREATE TABLE release (id TEXT PRIMARY KEY, created_at TEXT, kind TEXT, manifest TEXT);

INSERT INTO project VALUES ('APHCP', 'APHCP 蛋白–配体资源平台', 'M0', 'IN_PROGRESS');
INSERT INTO module VALUES
    ('holo', 'Holo–Ligand', 'NOT_STARTED'),
    ('apo', 'Apo 结构', 'NOT_STARTED'),
    ('affinity', '实验亲和力', 'NOT_STARTED'),
    ('plip', '相互作用注释', 'NOT_STARTED'),
    ('similarity', '相似性', 'NOT_STARTED'),
    ('aphcp', 'APHCP 配对与发布', 'NOT_STARTED');

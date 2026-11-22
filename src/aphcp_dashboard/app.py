"""Local APHCP research console."""
import csv
import html
import io
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import yaml

from aphcp_control import transport
from aphcp_control.inventory_lookup import canonical_id
from aphcp_control.queries import query_args
from aphcp_dashboard import evidence

ROOT = evidence.ROOT
ZONE = timezone(timedelta(hours=8))
PAGES = {'home':'流程管理', 'm1':'全库 M1', 'm2':'结构与 QC', 'pdb':'PDB 搜索', 'batches':'批次与任务',
         'runs':'运行记录', 'server':'服务器与代码', 'govern':'模块与治理'}
NAMES = {'PASS':'通过', 'REVIEW':'待核', 'DIFFERENCE':'发现差异', 'FAIL':'失败', 'NO_TARGET':'无目标',
         'INVENTORY_COMPLETE':'已验收', 'PILOT_COMPLETE':'试跑完成', 'NOT_STARTED':'未启动', 'SUCCEEDED':'成功', 'FAILED':'失败', 'RUNNING':'运行中',
         'PENDING':'待执行', 'PLANNED':'已计划', 'BLOCKED':'受阻', 'READY':'就绪', 'OPEN':'待处理',
         'COMPLETED':'已完成', 'ACTIVE':'进行中', 'APPROVED':'已批准', 'CLOSED':'已关闭', 'PAUSED':'已暂停',
         'CANCELLED':'已取消', 'DRAFT':'草案', 'UNSET':'未设', 'PROVISIONAL':'暂定', 'QUEUED':'排队中',
         'PLANNING':'规划中', 'ACCEPTED':'已验收', 'SKIPPED':'已跳过', 'not_run':'未运行',
         'MATCH':'一致', 'DIFFERENT':'不同', 'UNKNOWN':'未知', 'REUSE_ACCEPTED':'复用', 'NEW_WINDOW':'新取'}
ROLES = {'A_NONPOLYMER':'非聚合物', 'B_BRANCHED':'分支实体', 'C_PARTNER_CLUE':'伙伴线索', 'INSUFFICIENT':'线索不足'}
LABELS = {'id':'标识', 'status':'状态', 'name':'名称', 'title':'标题', 'batch_id':'批次', 'task_id':'任务',
          'source':'来源', 'run_id':'运行', 'start_time_utc':'开始（北京时间）', 'end_time_utc':'结束（北京时间）',
          'updated_at':'更新（北京时间）','run_count':'运行数', 'pdb_id':'PDB ID', 'stage':'阶段',
          'error_message':'失败原因', 'error_code':'错误代码', 'report_path':'报告路径', 'output_uri':'产物路径'}


def esc(v):
    return html.escape(str(v))


def n(v):
    return f'{v:,}' if isinstance(v, int) and not isinstance(v, bool) else '—'


def stamp(value):
    if not value:
        return '未知'
    try:
        t = datetime.fromisoformat(value.replace('Z', '+00:00'))
        return t.astimezone(ZONE).strftime('%Y-%m-%d %H:%M:%S') if t.tzinfo else str(value)
    except (AttributeError, ValueError):
        return str(value)


def status(v):
    return NAMES.get(v, v or '未知')


def badge(label, tone='neutral'):
    return f'<span class="ap-badge ap-{tone}">● {esc(label)}</span>'


def metric_grid(values):
    cards = ''.join(f'<div class="ap-metric"><label>{esc(k)}</label><span class="ap-value">{esc(v)}</span><small>{esc(h)}</small></div>'
                    for k,v,h in values)
    st.html(f'<div class="ap-grid {"three" if len(values)==3 else ""}">{cards}</div>')


def heading(title, desc):
    st.html(f'<div class="ap-compact"><h1>{esc(title)}</h1><p>{esc(desc)}</p></div>')


def go_to(page, pdb=None):
    st.session_state.pending_nav = page
    st.query_params['page'] = page
    if pdb:
        st.session_state.pending_pdb_input = pdb
        st.query_params['pdb'] = pdb
    else:
        st.query_params.pop('pdb', None)
    st.rerun()


def records(data, empty='暂无匹配记录。', *, key=None, height=None, registry=False):
    if not data:
        st.caption(empty)
        return
    rows = [data] if isinstance(data, dict) else data
    if registry:
        rows = [{LABELS.get(k,k): (status(v) if k=='status' else stamp(v) if k in ('start_time_utc','end_time_utc','updated_at') and v else v)
                 for k,v in row.items()} for row in rows]
    rows = [{k: ', '.join(map(str,v)) if isinstance(v,list) else str(v) if isinstance(v,dict) else v for k,v in row.items()} for row in rows]
    st.dataframe(pd.DataFrame(rows), hide_index=True, width='stretch', key=key, height=height or 'auto')


def lozenge(text, tone='neutral', note=''):
    st.html(f'<div class="ap-strip">{badge(text,tone)}{esc(note)}</div>')


def source_panel(e, *keys):
    with st.expander('来源与更新时间'):
        rows = [{'本地证据':e['sources'][k]['path'], '本地文件更新时间（北京时间）':stamp(e['sources'][k]['mtime'])} for k in keys if k in e['sources']]
        records(rows, '当前没有可读取的本地证据。')
        if e['issues']:
            st.caption('未读取的回执：' + '；'.join(e['issues']))


def chart(fig, key, height=210):
    fig.update_layout(height=height, margin=dict(l=5,r=8,t=13,b=5), paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)', font=dict(size=12),
                      legend=dict(orientation='h',y=-.11))
    st.plotly_chart(fig, key=key, width='stretch', config={'displayModeBar':False,'responsive':True})


def server_key(settings):
    ssh, server = settings['ssh'], settings['server']
    return ssh['user'], ssh['host'], ssh['identity'], server['data_root'], server['code_root']


def cached_payload(key, retrieval, settings, *, force=False, ttl=45):
    now = datetime.now(timezone.utc)
    endpoint = server_key(settings) if settings else None
    key = (endpoint,*key)
    store = st.session_state.setdefault('registry_cache', {})
    item = store.get(key) or (next((v for k,v in reversed(list(store.items())) if k[1:]==key[1:]), None) if settings is None else None)
    if settings is None:
        return dict(payload=(item or {}).get('payload'), state='配置异常', failure=st.session_state.get('settings_issue'), attempted=None, stale=bool(item))
    if not force and item and (now-item['attempted']).total_seconds() < ttl:
        return item
    old = (item or {}).get('payload')
    try:
        new = dict(payload=retrieval(), state='在线', attempted=now, stale=False, failure=None)
    except transport.QueryError as exc:
        new = dict(payload=old, state=exc.state, attempted=now, stale=bool(old), failure=str(exc))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        new = dict(payload=old, state='配置或查询异常', attempted=now, stale=bool(old), failure=str(exc))
    store[key] = new
    if len(store) > 40:
        store.pop(next(iter(store)))
    return new


def registry_fetch(settings, view, value=None, *, force=False):
    return cached_payload(('registry',view,value), lambda:transport.fetch(settings,view,value), settings, force=force)


def query_status(item, name='Pro5000 实时查询'):
    q = item.get('payload') or {}
    attempted = item.get('attempted')
    label = item['state']; last = f"上次成功：{stamp(q.get('queried_at'))}（北京时间）" if q else '没有成功查询缓存'
    if attempted:
        last += f' · 最近检查：{stamp(attempted.isoformat())}'
    if q.get('detail_limit'):
        last += f" · 登记明细上限 {q['detail_limit']} 条"
    lozenge(f'{name}：{label}', 'good' if label=='在线' else 'warn', last)
    if item.get('failure'):
        st.caption(('以下实时区块显示上次缓存。' if item.get('stale') else '') + item['failure'])
    return q.get('data')


def overview(e, settings, force):
    counts = e['m1'].get('counts', {})
    m1_done = evidence.accepted(e)
    v1 = e['v1'].get('execution', {}).get('disk_verify', {})
    v1_pass = v1.get('structure_states', {}).get('PASS')
    v1_total = len(e['v1'].get('entries', {}))
    copy_done = e['copy'].get('status') == 'PASS'
    workflow = e.get('workflow', {})

    heading('流程管理', '看当前关口、决定下一步，并追到支持决定的批次和验收证据。')
    with st.container(border=True):
        st.markdown('### 当前要推进的工作')
        note = workflow.get('当前待办')
        if note:
            st.markdown(note)
        else:
            st.caption('当前待办尚无可读的项目状态记录。')
        if scope := workflow.get('授权范围'):
            with st.expander('当前执行范围'):
                st.markdown(scope)
        current, batches = st.columns(2)
        if current.button('查看结构批次的验收依据' if m1_done else '查看候选盘点进度',
                          key='flow-current'):
            go_to('m2' if m1_done else 'm1')
        if batches.button('查看批次与任务', key='flow-batches'):
            go_to('batches')

    with st.container(border=True):
        st.markdown('### 各阶段做到哪里')
        m1_summary = (f"{n(e['m1'].get('window_count'))} 窗已对账 · {n(counts.get('entries'))} 个唯一候选 Entry"
                      if m1_done else '等待可核对的窗口验收回执')
        m2_checks = []
        if v1_pass is not None:
            m2_checks.append(f'V1 文件回读 {v1_pass}/{v1_total}')
        if copy_done:
            m2_checks.append('3TCT 组装复制通过')
        steps = [
            ('M1', '全库候选盘点', m1_summary, '已验收' if m1_done else '待确认', 'good' if m1_done else 'neutral'),
            ('M2', '结构试跑和扩量', ' · '.join(m2_checks) or '等待结构验收回执',
             '部分验收' if m2_checks else '待确认', 'warn' if m2_checks else 'neutral'),
            ('M3', 'Holo Master 发布', '发布依据来自正式资源发布记录', '查看登记', 'neutral'),
        ]
        st.html(''.join(f'<div class="ap-step"><span class="ap-num">{key}</span><div><b>{esc(title)}</b>'
                        f'<small>{esc(detail)}</small></div>{badge(state, tone)}</div>'
                        for key, title, detail, state, tone in steps))

        waiting = []
        newest = {case['pdb']: evidence.newest_cases(e, case['pdb'])[0] for case in e['cases']}
        for case in newest.values():
            if case.get('role') == 'REVIEW' and case.get('note'):
                waiting.append(f"{case['pdb']} · {case['note']}")
        if waiting:
            st.markdown('**需要判断的结构项**')
            for item in waiting:
                st.markdown(f'- {item}')

    links = st.columns(3)
    for col, (target, label) in zip(links, [('m1', 'M1 验收与覆盖'), ('runs', '运行与产物'), ('pdb', '按 PDB 查证据')]):
        if col.button(label, width='stretch', key='flow-' + target):
            go_to(target)

    with st.expander('正式登记与服务器查询'):
        live = st.button('读取当前登记快照', key='home-live')
        if live:
            st.session_state.home_live = True
        if st.session_state.get('home_live'):
            d = query_status(registry_fetch(settings, 'project', force=live or force))
            if d:
                st.markdown('**在途批次**')
                records(d.get('active_batches'), '当前登记没有在途批次。', registry=True)
                st.markdown('**运行中的任务**')
                records(d.get('active_tasks'), '当前登记没有运行中的任务。', registry=True)
                st.markdown('**登记中的阻塞或问题**')
                records(d.get('blockers'), '当前登记没有阻塞或开放问题。', registry=True)
    source_panel(e, 'workflow', 'm1_audit', 'v1', 'copy')


def csv_download(name, rows, label='导出当前摘要 CSV'):
    if rows:
        stream=io.StringIO(); writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader();writer.writerows(rows)
        st.download_button(label, data=stream.getvalue().encode('utf-8-sig'), file_name=name, mime='text/csv', on_click='ignore',key=name)


def budget(e):
    scenarios=e['m1'].get('coordinate_budget',{}); verified=e['audit'].get('file_cache_status',{})
    st.markdown('**原始文件与坐标预算**')
    st.caption(f"唯一计划文件 {n(e['m1'].get('unique_coordinate_files'))} 份 · 已验证复用 {n(verified.get('VERIFIED_REUSE'))} 份 · 本轮实测坐标下载 {n(e['downloaded'])} 份。预算是估算。")
    fields={'download_bytes':'预计网络下载','compressed_original_bytes':'压缩原件','derived_cif_scenario_bytes':'派生 CIF 情景',
            'largest_window_serial_scratch_bytes':'最大单窗串行临时空间'}
    records([{'情景':k.upper(), **{label: round(s[field]/1e9,2) if field in s else None for field,label in fields.items()}}
             for k,s in scenarios.items()])
    st.caption('单位：GB（10⁹ B）。临时空间取最大单窗串行估算，并发运行峰值需另估。')


def window_rows(e):
    rows=[]
    for w in e['windows']:
        begin,_,end=w.get('bounds','').partition('/')
        rows.append({'窗口':w['window'],'开始（含）':begin,'结束（不含）':end,'计划 Entry':w.get('planned_count'),
                     '实收 Entry':w.get('search_total_count'),'状态':status(w.get('receipt_status')),
                     '来源':status(w.get('reuse')), '元数据待补':len(w.get('metadata_missing_ids',[])),
                     '范围待核':len(w.get('period_missing_ids',[])), '查询时间（北京时间）':stamp(w.get('search_utc'))})
    return rows


def m1_page(e):
    heading('全库 M1 · 候选盘点','RCSB 实验蛋白 Entry · 窗口按首次公开日期左闭右开 · Pro5000 保留完整母表')
    inv=e['m1']; c=inv.get('counts',{}); audit=e['audit']; cross=e['cross']; rs=audit.get('plinder_revision_status_annotated_only',{})
    lozenge('盘点已验收' if evidence.accepted(e) else '盘点状态未知', 'good' if evidence.accepted(e) else 'neutral',
            f"范围：{inv.get('scope') or '未知'} · 复验：{stamp(audit.get('verification_finished_utc'))}（北京时间）")
    metric_grid([('唯一 Entry',n(c.get('entries')),'RCSB 母表'),('有 PLINDER 注释',n(c.get('plinder_annotation_entries')),'候选索引线索'),
                 ('无 PLINDER 注释',n(c.get('plinder_no_annotation_entries')),'保留原件识别分支'),('候选 System',n(c.get('plinder_candidate_systems')),'尚未做全库结构 QC')])
    with st.container(border=True):
        st.markdown('**分母与版本对账**')
        st.write(f"冻结分母 {n(inv.get('denominator',{}).get('total_count'))} · {n(inv.get('window_count'))} 窗 · 分母差额 {n(cross.get('denominator_minus_unique_ids'))}。")
        st.caption('首次分母查询：'+stamp(inv.get('denominator',{}).get('queried_utc'))+'（北京时间）')
        st.write(f"注释命中内：修订一致 {n(rs.get('MATCH'))}、不同 {n(rs.get('DIFFERENT'))}、未知 {n((c.get('plinder_annotation_entries')-rs['MATCH']-rs['DIFFERENT']) if isinstance(c.get('plinder_annotation_entries'),int) and {'MATCH','DIFFERENT'} <= rs.keys() else None)}。")
        st.caption(f"Assembly 对 {n(c.get('rcsb_assembly_pairs'))} · 配体注释行 {n(c.get('plinder_annotation_rows'))} · 全量修订比较含无注释 Entry，见覆盖详情。")
    with st.container(border=True): budget(e)
    st.subheader('窗口台账')
    rows=window_rows(e)
    c1,c2,c3=st.columns([2,1,1]); word=c1.text_input('搜索窗口／起止日期',placeholder='例如 2021-09'); states=c2.selectbox('状态',['全部',*sorted({x['状态'] for x in rows})]); reuse=c3.selectbox('取得方式',['全部',*sorted({x['来源'] for x in rows})]);
    rows=[w for w in rows if (not word.strip() or word.strip().casefold() in ' '.join(map(str,w.values())).casefold()) and (states=='全部' or w['状态']==states) and (reuse=='全部' or w['来源']==reuse)]
    records(rows,'没有窗口符合筛选。',height=300,key='m1-windows'); st.caption(f'{len(rows)} 个显示窗口'); csv_download('aphcp-m1-windows.csv',rows)
    if rows:
        chosen=st.selectbox('窗口详情',[None,*[w['窗口'] for w in rows]],format_func=lambda value:value or '选择窗口')
        if chosen:
            data=next(w for w in e['windows'] if w['window']==chosen)
            st.caption(f"{data['bounds']} · 验收状态：{status(data.get('receipt_status'))}")
            st.code(data['run_directory'],language=None)
            with st.expander('此窗的回执字段'): st.json(data,expanded=False)
    with st.expander('覆盖分布与审核口径'):
        coverage=inv.get('coverage',{})
        for k,title in [('methods','实验方法（Entry，可重叠）'),('discovery_margins','候选发现路径（Entry，可重叠）'),('family_entry_evidence','蛋白家族证据（Entry，可重叠）'),('indexed_ligand_row_classes','配体类别（注释行）'),('pli_binding_evidence','结合线索（注释行，可重叠）'),('all_source_revision_comparisons','全母表修订对账（Entry）')]:
            st.markdown(f'**{title}**'); records([{'类别':ROLES.get(x,status(x)),'数目':y} for x,y in coverage.get(k,{}).items()])
    source_panel(e,'m1','windows','m1_audit','m1_plan','m1_download')


def case_row(c):
    def join(v): return ', '.join(map(str,v)) if isinstance(v,list) else str(v) if v is not None else None
    return {'PDB':c['pdb'],'验收批次':c['group'],'修订':c['revision'],'Assembly':join(c['assembly']),'模型':join(c['models']),
            '阶段结论':status(c['status']),'坐标身份':status(c['coordinate']),'实验保真':status(c['fidelity']),
            '拓扑':status(c['topology']),'配体角色':status(c['role'])}


def case_detail(c, key):
    r=case_row(c)
    st.markdown(f"**{r['PDB']} · {r['验收批次']}**")
    st.caption(f"固定修订 {r['修订'] or '未知'} · Assembly {r['Assembly'] or '未知'} · 模型 {r['模型'] or '未知'} · {stamp(c.get('date'))}（北京时间）")
    metric_grid([(k,r[k], '阶段结论：'+r['阶段结论'] if k=='坐标身份' else '') for k in ('坐标身份','实验保真','拓扑','配体角色')])
    if c.get('note'):
        reasons={'TARGET_ROLE_PENDING':'目标角色待核','REPRESENTATIVE_MODEL_IS_NOT_COMPLETE_ENSEMBLE':'代表模型未覆盖全部实验模型',
                 'DNA_PARTNER_STRANDS_NOT_GROUPED_AS_DUPLEX':'早期试跑未归组双链', 'INDEX_GEOMETRY_OR_COMPONENT_REVIEW':'索引几何或组分待核',
                 'SOURCE_REVISION_DIFFERS':'来源修订不同，走原件核查','TARGET_CONTACT_OR_CONNECTION_UNRESOLVED':'目标接触或连接未定位',
                 'INDEX_LOOKUP_TECHNICAL_FALLBACK':'索引查询异常，走原件核查'}
        note=c['note']
        for word, label in reasons.items(): note=note.replace(word,label)
        st.write(note)
    with st.expander('产物路径和逐项证据'):
        if c.get('run'): st.code(c['run'],language=None)
        labels={'system':'体系 CIF','complete_assembly':'完整 Assembly CIF'}
        records([{'产物':labels.get(k,k),'服务器路径':(c['run']+'/'+v if c.get('run') and not v.startswith('/') else v)} for k,v in c['artifacts'].items() if v])
        if c.get('detail'): st.json(c['detail'],expanded=False)


def m2_page(e):
    heading('结构与 QC','9 月 M2 Pilot、V1 五例与组装复制补验 · 数字按各批分母，案例可重叠')
    p=e['pilot']; v=e['v1'].get('execution',{}).get('disk_verify',{}); cp=e['copy'];
    v_total=len(e['v1'].get('entries',{})); vs=v.get('structure_states',{}).get('PASS')
    metric_grid([('M2 Pilot 处理',n(p.get('processed')),'固定版候选 Entry'),('索引几何通过',n(p.get('index_mapping',{}).get('coordinate_identity_pass')),'候选 System · 不等于保真'),
                 ('V1 五例结构',f'{vs} / {v_total}' if vs is not None else '—','实物回读'),('3TCT 复制',status(cp.get('status')),'模型、坐标与来源')])
    states=p.get('states',{}); roles=v.get('role_states',{})
    if states or roles:
        st.caption('9 月试跑：'+' / '.join(f'{status(k)} {n(x)}' for k,x in states.items())+'。V1 五例角色：'+' / '.join(f'{status(k)} {n(x)}' for k,x in roles.items())+'。')
    if cp:
        st.info('复制组装带显式连接或碱基配对时，端点映射仍待实物补验。3TCT 原件没有显式连接。')
    rows=e['cases']; a,b,c,d=st.columns([1.3,1.4,1,1]); term=a.text_input('筛选 PDB',placeholder='例如 5R42').strip().upper();
    cohort=b.selectbox('证据批次',['全部',*dict.fromkeys(x['group'] for x in rows)]); qc=c.selectbox('阶段结论',['全部',*sorted({status(x['status']) for x in rows})]); role=d.selectbox('角色',['全部',*sorted({status(x['role']) for x in rows})]);
    rows=[x for x in rows if (not term or term in x['pdb']) and (cohort=='全部' or cohort==x['group']) and (qc=='全部' or qc==status(x['status'])) and (role=='全部' or role==status(x['role']))]
    records([case_row(x) for x in rows],'没有案例符合筛选。',height=300,key='cases');
    if rows:
        selection=st.selectbox('打开案例详情',range(len(rows)),format_func=lambda i:f"{rows[i]['pdb']} · {rows[i]['group']}",key='m2-case');
        with st.container(border=True): case_detail(rows[selection], 'm2-case')
    source_panel(e,'pilot','pilot_run','pilot_diag','v1','copy')


def pdb_page(e, settings, force):
    heading('PDB 搜索','精确查询：已验收 M1 母表、M2 结构证据和正式 SQLite 分开呈现')
    if 'pdb_input' not in st.session_state:
        st.session_state.pdb_input=st.query_params.get('pdb','')
    with st.form('pdb_query'):
        typed=st.text_input('PDB ID',key='pdb_input',placeholder='例如 3TCT 或 1MBN',max_chars=16)
        submitted=st.form_submit_button('查询 PDB',type='primary')
    if submitted:
        if not typed.strip():
            st.query_params.pop('pdb',None)
            st.warning('请输入 PDB ID。'); return
        st.query_params['pdb']=typed.strip()
    typed=st.query_params.get('pdb','')
    if not typed:
        st.caption('可试查'); a,b=st.columns(2)
        if a.button('3TCT · 结构及组装回读',width='stretch'): go_to('pdb','3TCT')
        if b.button('1MBN · M1 母表候选',width='stretch'): go_to('pdb','1MBN')
        return
    try:
        _, identifier=query_args('pdb',typed)
    except ValueError:
        st.error('PDB ID 格式不正确。请输入 4 位字母数字，或 pdb_ 加 8 位字母数字。'); return
    identifier=canonical_id(identifier)
    cases=evidence.newest_cases(e,identifier)
    st.subheader(identifier)
    force=submitted or force
    if settings:
        # Start the independent reads together, but mutate Streamlit state on its own thread.
        with ThreadPoolExecutor(max_workers=2) as pool:
            preload={}
            for key, fetch in [(('inventory',identifier),lambda:transport.fetch_inventory(settings,identifier)),
                               (('registry','pdb',identifier),lambda:transport.fetch(settings,'pdb',identifier))]:
                endpoint=server_key(settings)
                item=st.session_state.get('registry_cache',{}).get((endpoint,*key)); now=datetime.now(timezone.utc)
                if force or not item or (now-item['attempted']).total_seconds()>=45:
                    preload[key]=pool.submit(fetch)
            inv=cached_payload(('inventory',identifier),lambda:preload['inventory',identifier].result() if ('inventory',identifier) in preload else transport.fetch_inventory(settings,identifier),settings,force=force)
            reg=cached_payload(('registry','pdb',identifier),lambda:preload['registry','pdb',identifier].result() if ('registry','pdb',identifier) in preload else transport.fetch(settings,'pdb',identifier),settings,force=force)
    else:
        inv=cached_payload(('inventory',identifier),None,settings); reg=registry_fetch(settings,'pdb',identifier)
    inventory=(inv.get('payload') or {}).get('data') or {}; official=(reg.get('payload') or {}).get('data') or {}
    entry=inventory.get('candidate')
    no_entry = '未找到' if inventory.get('state')=='absent' else '未能确认'
    metric_grid([('M1 候选盘点','已盘点' if entry else no_entry,'冻结母表 · 非结构验收'),
                 ('M2 结构证据','已有案例' if cases else '未见本地回执',f'{len(cases)} 份本地案例回执'),
                 ('正式 SQLite','有登记条目' if official.get('entry') else ('无登记条目' if reg.get('payload') else '未能确认'),'登记条目不等于正式结构入库或发布')])
    if not entry and not cases and not official.get('entry'):
        if inventory.get('state')=='absent' and reg.get('payload') and not inv.get('stale') and not reg.get('stale'):
            st.info('三个来源均未找到该 PDB ID。')
        else: st.info('可用来源未找到；至少一个远端来源不可用或不完整。')
    st.markdown('**M1 · 已验收母表**'); query_status(inv,'M1 只读查询')
    if entry:
        st.write(entry.get('title') or identifier)
        hit=entry['plinder_index_presence']
        records({'首次公开日':(entry.get('first_release_date') or '')[:10], '最新修订':entry.get('rcsb_latest_revision'),
                 '拟取修订':entry.get('selected_revision'),'Assembly':entry.get('rcsb_assembly_ids'),'实验模型':entry.get('model_count'),
                 '实验方法':entry.get('methods'),'PLINDER 注释':'有' if hit is True else '无' if hit is False else '未知',
                 '候选 System':entry.get('plinder_candidate_system_count'),'索引修订':entry.get('plinder_source_revision_if_available'),
                 '修订比较':status(entry.get('plinder_revision_status'))})
        st.caption('候选发现线索：'+('、'.join(ROLES.get(x,x) for x in entry.get('discoveries') or []) or '未知')+
                   ' · 非聚合物 CCD：'+('、'.join(entry.get('nonpolymer_ccds') or []) or '未见'))
        with st.expander('分窗来源与候选证据'):
            st.code(inventory['source'],language=None); st.caption('左闭右开窗口：'+' / '.join(inventory['bounds']))
            st.json(entry,expanded=False)
    elif inventory.get('state')=='absent': st.caption('冻结母表未匹配此 ID。')
    elif inventory.get('state')=='incomplete': st.warning('分窗清单或条目文件不完整，暂不能确认。')
    if cases:
        st.markdown('**M2 · 本地结构证据**')
        for i,item in enumerate(cases):
            with st.container(border=True): case_detail(item, f'pdb-case-{i}')
    st.markdown('**正式 SQLite 明细**'); query_status(reg,'登记查询')
    if official.get('entry'):
        for name,label in [('entry','登记条目'),('versions','结构版本'),('files','文件及取数'),('history','阶段记录')]:
            with st.expander(label,expanded=name=='entry'): records(official.get(name),registry=True)
        for name,label in [('assembly','组装体'),('chain','链实例'),('ligand','配体实例'),('system','体系')]:
            with st.expander(label): records(official.get('objects',{}).get(name),registry=True)
    elif reg.get('payload'):
        st.caption('正式 SQLite 无此条目。M1 和 M2 的回执结果仍可独立查看。')


def evidence_item(e):
    batch=e['plan'].get('batch_id'); wins=e['windows']
    if not batch or not e['cross']: return None
    return {'id':batch,'status':e['cross'].get('status'),'run_count':len(wins),'name':'全库 M1 · 验收回执','source':'独立验收回执'}


def receipt_runs(e):
    out=[]
    if evidence_item(e):
        out += [{'id':w['window'],'batch_id':e['plan']['batch_id'],'status':w.get('receipt_status'),'start_time_utc':w.get('search_utc'),
                 'receipt':w,'source':'独立验收回执'} for w in e['windows']]
    for label, item, date, identifier in [('M2 Pilot',e['pilot_run'],e['pilot_run'].get('finished_utc'),e['pilot'].get('batch_id')),
            ('V1 五例',e['v1'],e['v1'].get('start_build_utc'),None),
            ('3TCT 补验',e['copy'],e['copy'].get('invocations',{}).get('build',{}).get('finished_utc'),None)]:
        path=item.get('run_path') or item.get('run')
        if path:
            out.append({'id':Path(path).name,'batch_id':identifier,'status':item.get('status') or ('PASS' if e['v1'].get('execution',{}).get('disk_verify',{}).get('structure_states',{}).get('PASS') else None),
                        'start_time_utc':None,'end_time_utc':date,'name':label,'receipt':{'产物目录':path},'source':'独立验收回执'})
    return out


def choose_registry_rows(rows, mode):
    a,b=st.columns([2,1]); word=a.text_input('搜索 ID、名称或批次',key=mode+'-search').strip().casefold(); filt=b.selectbox('状态',['全部',*sorted({status(x.get('status')) for x in rows})],key=mode+'-state')
    view=[x for x in rows if (not word or word in ' '.join(str(x.get(k,'')) for k in ('id','name','batch_id','task_id')).casefold()) and (filt=='全部' or status(x.get('status'))==filt)]
    keys=('id','name','batch_id','status','source','run_count','start_time_utc','end_time_utc'); rows=[{k:x.get(k) for k in keys if x.get(k) is not None} for x in view]
    records(rows,'没有符合筛选的记录。',height=270,key=mode+'-table',registry=True)
    return st.selectbox('选择一条查看详情',[None,*view],format_func=lambda x: f"{x['id']} · {x.get('source','SQLite')}" if x else '请选择', key=mode+'-selected') if view else None


def batch_page(e, settings, force):
    heading('批次与任务','筛选后打开批次；独立验收数据与正式 SQLite 登记分别呈现')
    d=query_status(registry_fetch(settings,'batches',force=force),'批次登记') or {}; official=d.get('batches',[])
    local=evidence_item(e); ids={r['id'] for r in official}; rows=[local] if local and local['id'] not in ids else []
    rows+=official
    item=choose_registry_rows(rows,'batch')
    if not item: return
    st.subheader(item['id']);
    if local and item['id']==local['id']:
        lozenge('独立验收回执', 'good' if evidence.accepted(e) else 'warn','本次全库运行未写入正式登记。')
        st.caption(f"{n(len(e['windows']))} 窗 · 母表 Entry {n(e['m1'].get('counts',{}).get('entries'))}")
        records(window_rows(e),height=250,key='batch-windows')
        if st.button('打开全库 M1'): go_to('m1')
    if item['id'] in ids:
        data=query_status(registry_fetch(settings,'batch',item['id'],force=force),'批次详情查询')
        if data:
            if not data.get('batch'): st.info('未找到批次详情。')
            else:
                records(data['batch'],registry=True)
                for name,title in [('tasks','任务'),('runs','相关运行'),('latest_stage_outcomes','最新阶段结果'),('usage','资源消耗')]:
                    with st.expander(title,expanded=name=='tasks'):
                        items=data.get(name)
                        records(items,registry=True)
                        if name=='tasks' and items:
                            task=st.selectbox('任务详情',[None,*[x['id'] for x in items]],format_func=lambda v:v or '选择任务')
                            if task:
                                records(next(x for x in items if x['id']==task),registry=True)


def run_page(e, settings, force):
    heading('运行记录','状态与时间统一显示；只展示服务器返回的有界明细和已知产物路径')
    d=query_status(registry_fetch(settings,'runs',force=force),'运行登记') or {}; official=d.get('runs',[]); official_ids={x['id'] for x in official}
    local={x['id']:x for x in receipt_runs(e)}
    item=choose_registry_rows([*official,*[r for k,r in local.items() if k not in official_ids]],'run')
    if not item:return
    st.subheader(item['id']); rid=item['id']
    if rid in local:
        with st.expander('独立验收与产物',expanded=True):
            val=local[rid]['receipt']
            if 'run_directory' in val: st.code(val['run_directory'],language=None)
            elif '产物目录' in val: st.code(val['产物目录'],language=None)
            st.json(val,expanded=False)
    if rid in official_ids:
        data=query_status(registry_fetch(settings,'run',rid,force=force),'运行详情查询')
        if data:
            record=data.get('run')
            if not record: st.info('未找到运行详情。');return
            if record.get('error_message'): st.warning(record['error_message'])
            records(record,registry=True)
            for k in ('report_path','output_uri'):
                if record.get(k): st.code(record[k],language=None)
            with st.expander(f"单项结果 · {data.get('item_count',0)} 条登记",expanded=True): records(data.get('items'),'此运行没有单项登记。',registry=True)
            if data.get('report_preview'):
                with st.expander('运行报告预览'): st.code(data['report_preview'],language=None)


def server_page(e, settings, force):
    heading('服务器与代码','Pro5000 实时资源、登记快照与本地 Git 状态')
    mode=st.radio('服务器面板',['资源与存储','代码版本'],horizontal=True)
    if mode=='代码版本':
        data=query_status(registry_fetch(settings,'deployment',force=force)) or {}
        try: local=transport.git_state(ROOT)
        except OSError: local={}
        remote=data.get('current')
        if not local or not remote: st.info('同步状态未知：需要本地和服务器两个版本。')
        else: lozenge('代码已同步' if transport.in_sync(local,remote) else '代码存在差异或本地有在手修改','good' if transport.in_sync(local,remote) else 'warn')
        st.markdown('**本地 Git**'); records(local)
        st.markdown('**服务器部署**'); records(remote)
        with st.expander('历史部署与运行版本'): records(data.get('history'),registry=True);records(data.get('run_commits'),registry=True)
        return
    data=query_status(registry_fetch(settings,'storage',force=force)) or {}
    if not data:
        st.caption('资源实时值未知。验收回执仍可从总览、M1 和 M2 查看。'); return
    host=data.get('host',{}); total=host.get('MemTotal_bytes'); remain=host.get('MemAvailable_bytes'); cpu=host.get('cpu_used_percent')
    metric_grid([('CPU 使用率',f'{cpu:.1f}%' if cpu is not None else '—',f"逻辑核心 {n(host.get('cpu_count'))}"),
                 ('可用内存',f'{remain/2**30:.1f} GiB' if remain is not None else '—',f"总量 {total/2**30:.1f} GiB" if total else '总量未知'),
                 ('活跃运行',n(len(data.get('active_runs',[]))),'正式登记状态')])
    disks=data.get('disks',{}); keys=[k for k in ('data','code') if k in disks]
    if keys:
        names=['数据盘' if k=='data' else '代码盘' for k in keys]; fig=go.Figure()
        for field,label,color in [('used','已用','#14899d'),('free','可用','#bfd2e3')]:
            fig.add_bar(x=[disks[k][field]/2**30 for k in keys],y=names,orientation='h',name=label,marker_color=color,hovertemplate='%{x:,.1f} GiB<extra></extra>')
        fig.update_layout(barmode='stack',xaxis_title='GiB'); chart(fig,'disk',190)
    st.caption('资源采样：'+stamp(data.get('sampled_at'))+'（北京时间）')
    with st.expander('运行、分层磁盘与登记资源'):
        records(data.get('active_runs'),registry=True); records(data.get('layer_snapshot')); records(data.get('transfers')); records(data.get('resources'))


def governance_page(e, settings, force):
    heading('模块与治理','模块状态、决策参数与本地历史台账')
    mode=st.radio('治理面板',['决策与参数','模块状态'],horizontal=True)
    if mode=='模块状态':
        data=query_status(registry_fetch(settings,'modules',force=force)) or {}; records(data.get('modules'),registry=True); return
    data=query_status(registry_fetch(settings,'governance',force=force)) or {}
    for name,title,filename in [('decision','决策','decision_log.csv'),('threshold','阈值','threshold_registry.csv'),('issue','问题','issue_register.csv'),('config_version','配置版本',None)]:
        with st.expander(title+' · SQLite / 本地台账'):
            records(data.get(name),'没有可用的正式记录。',registry=True)
            if filename:
                st.caption('本地台账最近 25 条')
                try:
                    with (ROOT/'management'/filename).open(encoding='utf-8-sig',newline='') as f: records(list(csv.DictReader(f))[-25:])
                except OSError: st.caption('本地台账不可用。')
    with st.expander('本地存储预算台账'):
        try:
            with (ROOT/'management'/'storage_ledger.csv').open(encoding='utf-8-sig',newline='') as f:
                records(list(csv.DictReader(f))[-25:])
        except OSError:
            st.caption('本地台账不可用。')


def main():
    st.set_page_config(page_title='APHCP 科研控制台',page_icon='🔬',layout='wide',initial_sidebar_state='auto')
    st.html((Path(__file__).with_name('dashboard.html')).read_text(encoding='utf-8-sig'))
    try:
        settings=transport.load_settings(ROOT)
        refresh=max(10,min(300,int(settings['dashboard']['refresh_seconds'])))
        st.session_state.pop('settings_issue',None)
    except (OSError, ValueError, KeyError, TypeError, yaml.YAMLError) as error:
        settings=None; refresh=30; st.session_state.settings_issue=str(error)
    initial=st.query_params.get('page','home')
    if 'nav' not in st.session_state:
        st.session_state.nav=initial if initial in PAGES else 'home'
    if pending := st.session_state.pop('pending_nav',None):
        st.session_state.nav=pending
    if pending := st.session_state.pop('pending_pdb_input',None):
        st.session_state.pdb_input=pending
    def change_nav():
        st.query_params['page']=st.session_state.nav
        if st.session_state.nav!='pdb': st.query_params.pop('pdb',None)
    with st.sidebar:
        st.html('<div class="ap-brand"><div class="ap-mark">A</div><div><strong>APHCP</strong><small>实验蛋白–配体资源平台</small></div></div>')
        st.caption('科研进度与证据')
        page=st.radio('导航',list(PAGES),format_func=PAGES.__getitem__,label_visibility='collapsed',key='nav',on_change=change_nav)
        st.divider(); st.caption('PRO5000 · 本机只读控制台')
        if settings is None: st.warning('服务器配置异常，本地验收回执仍可用。')
        else: st.caption('实时状态与缓存时间显示在各查询页。')
        online=page in ('batches','runs','server','govern')
        auto=st.toggle('定时刷新当前页',value=False,disabled=not online or settings is None)
        if auto and online: st.caption(f'每 {refresh} 秒更新当前页。')
    @st.fragment(run_every=refresh if auto and online and settings else None)
    def render():
        with st.container(horizontal=True,horizontal_alignment='right'):
            force=st.button('刷新本页',key='refresh-'+page)
        e=evidence.load(ROOT)
        if page=='home': overview(e,settings,force)
        elif page=='m1': m1_page(e)
        elif page=='m2': m2_page(e)
        elif page=='pdb': pdb_page(e,settings,force)
        elif page=='batches': batch_page(e,settings,force or auto)
        elif page=='runs': run_page(e,settings,force or auto)
        elif page=='server': server_page(e,settings,force or auto)
        else: governance_page(e,settings,force or auto)
    render()


main()

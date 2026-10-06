"""Build a Chinese PDF report directly from the frozen numerical comparison."""
from pathlib import Path
import json, csv, itertools, html, os, subprocess, tempfile
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'report';OUT.mkdir(exist_ok=True)
render_dir=tempfile.TemporaryDirectory(prefix='unified_report_figures_')
for name in ['baseline_effect','source_comparison','primary_session_effects','geometry_comparison']:
 subprocess.run(['pdftoppm','-r','160','-singlefile','-png',str(ROOT/'figures'/(name+'.pdf')),str(Path(render_dir.name)/name)],check=True,capture_output=True)
def read(n):return json.loads((ROOT/n).read_text())
rows=read('summary.json');look={(x['scenario'],x['metric'],x['group']):x for x in rows}
base=read('fig2_baseline_summary.json');blook={(x['metric'],x['baseline_rule'],x['group']):x for x in base}
raw=read('raw_rebuild_checks.json');cache=read('fig2_cache_checks.json');counts=read('count_matching.json')
primary='union_common240__all__thirds'
# Equal-animal summaries are deterministic descriptive statistics, not a new inference model.
balanced=[]
for x in rows:
 if x['group']!='pooled':continue
 a=look.get((x['scenario'],x['metric'],'M1'));b=look.get((x['scenario'],x['metric'],'M2'))
 if not a or not b:continue
 v=[a['mean'],b['mean']];mu=sum(v)/2;null=[(i*v[0]+j*v[1])/2 for i,j in itertools.product([-1,1],repeat=2)]
 balanced.append(dict(scenario=x['scenario'],metric=x['metric'],M1=v[0],M2=v[1],animal_balanced_mean=mu,n_animals=2,exact_signflip_p=sum(abs(t)>=abs(mu)-1e-12 for t in null)/4))
(ROOT/'animal_balanced_summary.json').write_text(json.dumps(balanced,indent=2)+'\n')
alook={(x['scenario'],x['metric']):x for x in balanced}

def fmt(v,d=3):return 'NA' if v is None else f'{v:.{d}f}'
metrics=[('var_within_sdi_DiD','类别内 SDI DiD (%)'),('var_within_svm_loo_DiD','类别内 SVM DiD (%; LOO)'),('var_between_sdi','variation 类别间 SDI'),('var_between_svm','variation decoding (pp)'),('global_between_sdi','global 类别间 SDI'),('global_between_svm','global decoding (pp)'),('global_orthogonality_Polar_minus_Grating','Polar-Grating 角度变化 (deg)'),('global_RSA_Polar_net','Polar net RSA')]
main_table=[['指标','M1','M2','两动物等权']]
for key,title in metrics:
 x=alook[(primary,key)];main_table.append([title,fmt(x['M1']),fmt(x['M2']),fmt(x['animal_balanced_mean'])])

font_path=Path(os.environ.get('UNIFIED_REPORT_FONT','/System/Library/Fonts/Supplemental/Arial Unicode.ttf'))
if font_path.is_file():pdfmetrics.registerFont(TTFont('CJKReport',str(font_path)))
else:pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
FONT='CJKReport' if font_path.is_file() else 'STSong-Light'
styles=getSampleStyleSheet()
styles.add(ParagraphStyle(name='CBody',fontName=FONT,fontSize=10.2,leading=16,spaceAfter=8,wordWrap='CJK'))
styles.add(ParagraphStyle(name='CTitle',fontName=FONT,fontSize=21,leading=29,spaceAfter=15,textColor=colors.HexColor('#16344b')))
styles.add(ParagraphStyle(name='CHeading',fontName=FONT,fontSize=15,leading=21,spaceAfter=10,textColor=colors.HexColor('#16344b')))
styles.add(ParagraphStyle(name='CSmall',fontName=FONT,fontSize=8.6,leading=13,spaceAfter=6,wordWrap='CJK'))
styles.add(ParagraphStyle(name='CCell',fontName=FONT,fontSize=9,leading=13,wordWrap='CJK'))
story=[];md=[]
def heading(t):story.append(Paragraph(html.escape(t),styles['CHeading']));md.append('\n## '+t+'\n')
def para(t,small=False):story.append(Paragraph(html.escape(t),styles['CSmall' if small else 'CBody']));md.append(t+'\n')
def table(data,widths=None):
 cells=[[Paragraph(html.escape(str(v)),styles['CCell']) for v in r] for r in data]
 t=Table(cells,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e8eff4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#a7bac8')),('LINEBELOW',(0,1),(-1,-1),.25,colors.HexColor('#d7e0e6')),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6)]));story.append(t);story.append(Spacer(1,10))
 md.extend(['| '+' | '.join(map(str,data[0]))+' |','|'+'---|'*len(data[0]),*['| '+' | '.join(map(str,r))+' |' for r in data[1:]],''])
def page():story.append(PageBreak())
def figure(name,width=500):
 p=Path(render_dir.name)/(name+'.png')
 from PIL import Image as PILImage
 with PILImage.open(p) as im:w,h=im.size
 story.append(Image(str(p),width=width,height=width*h/w));story.append(Spacer(1,8));md.append(f'Companion figure: `unified/figures/{name}.pdf`\n')

story.append(Paragraph('统一分析与结果一致性报告',styles['CTitle']))
para('2026-10-06 | 分支 unified-analysis | 原 neuron 集合，不做 SNR 或 V2 重筛选',True)
heading('1. 结论与建议')
para('数据来源已经查清，并建立了可独立重跑的统一候选流程。但本次测试没有找到一套同时保持统一输入、统一比较范围和原图全部数值／显著性的方案。最关键的反例是 Fig2 的 M2 改善：仅把 baseline 改成最早可用 FR recording，同样的 post 日期就会使 SDI 与 SVM 的改善方向反转。')
para('建议把 union_common240 作为严格统一的描述性主流程：每个 recording 使用历史版本中已经出现过的 unit 集合并集；所有 recording 使用 630-870 ms 反应窗；统一类别定义、日期清单及指标计算。union_native 保留 630 ms 到记录的刺激结束时间，作为时间窗敏感性对照。方案选择基于输入一致性，而不是显著性或与旧图的相似程度。')
para('共同窗控制跨动物的积分时长，native 窗保留到各记录的刺激结束；两者都有明确解释。共同窗在本报告中是统一比较的锚点，不意味着删去 M2 后 130 ms 的响应在生理上必然更优。最终论文主窗应由实验设计或独立生理依据决定。',True)
para('当前 early/late thirds 仅是 recording 时间顺序分组，不等于已核实的训练前后。五个 session 来自两个动物；报告提供 session 效应及 M1/M2 分别汇总，再给两个动物等权的描述性平均，不把类别对和重复天数当作独立动物。')
table(main_table,[245,80,80,95])
para('表中为 late-early 或训练相对对照的 DiD。SVM DiD 是相对百分比差；类别间 decoding 是准确率百分点变化，两者不能混用。类别间下降与 Polar 正向变化较稳定；类别内改善、M2 扩张和几何显著性没有整体保留。',True)

page();heading('2. 不一致来自哪些输入差异')
para(f'核对了 52 个主 recording：22 个 variation、30 个 global。另一组 variation/prototype 输入有 50 个文件。所有 52 个 canonical FR 都从原始 CDT 重建；最大绝对差为 {max(x["max_abs_fr_difference"] for x in raw):.2e} Hz。原始文件 SHA-256 与 channel/unit 清单已保留。')
table([['同一 recording 的版本','canonical units','另一版本 units','差异'],['s1 variation 041116','28','34','unit 集合不同'],['s3 variation 071217','52','61','unit 集合不同'],['s2 global 031217','28','42','630 ms vs 500 ms 起始'],['s5 variation 082418','62','40','只有 38 个 unit 相同；并集 64']],[210,75,75,140])
para('因此不能把另一版本简单当作 canonical 的超集，也不能靠重新编号解决分析差异。union 是同一原始 recording 内的历史集合并集；它没有根据结果选择 unit，也没有跨 recording 假定 channel/unit 标签代表同一个生物学 neuron。unit-0 标签被保留，其是否代表独立单细胞尚需记录定义确认；不能将整套混合数据直接称为 V2-only single-unit data。')
table([['候选族','neuron / unit','反应窗与类别'],['stored_cross','原 variation/prototype 集合','原窗；含额外 U 标签'],['stored_cross_common_labels','与上行相同','只保留共同类别；用于隔离标签效应'],['canonical','原 var/global 集合','630 ms 至记录的刺激结束'],['union_native','原历史集合并集','630 ms 至刺激结束；共同类别'],['union_common240','同一并集','所有 recording 均 630-870 ms']],[135,170,195])
para('所有重建保留 350-500 ms 的 recording 级 pooled baseline 和原有 balanced-first-trial 规则。前三组窗口／类别并不完全相同，源族比较不能直接解释为纯 neuron 效应；canonical→union_native 才主要隔离 unit 集合，union_native→union_common240 隔离窗口长度。',True)

page();heading('3. Fig2：baseline 选择本身足以改变结论')
table([['session','缓存中的 baseline','最早可用 FR'],['s3 (M2)','071617','071217'],['s5 (M2)','081518','081318']],[150,175,175])
para('在同一套 fresh 重算结果中，只改变 baseline，保留完全相同的 post recording、neuron、SDI 公式和 LOO SVM 计算。下表沿用原图 post-day 平均权重，仅用于定位原因，不用作独立样本推断。')
data=[['M2 指标','fresh + 旧 baseline','fresh + 最早 FR','原缓存参考']]
for metric in ['sdi','svm']:
 old=blook[(metric,'fresh_cache_baseline','M2')];new=blook[(metric,'fresh_earliest_fr_baseline','M2')];ref=blook[(metric,'cache_baseline','M2')];data.append([metric.upper(),fmt(old['day_mean'],2)+'%',fmt(new['day_mean'],2)+'%',fmt(ref['day_mean'],2)+'%'])
table(data,[100,135,135,130]);figure('baseline_effect',500)
para('这并不证明最早可用 FR 是真正训练前 baseline；它证明当前结论依赖 baseline 的定义。不能为了保留正向改善而只选择更晚的 baseline。训练起点和累计 exposure 必须来自独立日志。')
para(f'旧 SDI 缓存与对应 fresh 类别均值的最大差约 {max(x["max_abs_difference"] for x in cache if x["metric"]=="sdi"):.2e}；SVM 类别均值可有小幅重算差异，因此 baseline 因果定位使用 fresh 对 fresh，而不是混用旧 post 与新 baseline。',True)

page();heading('4. 输入版本对方向的影响')
para('下图固定两族都存在的 recording 日期、相同 early/late thirds 和相同估计器。类别间 SVM 为五折，类别内原 LOO 另行检验。误差条为 session SEM，不是动物层面的置信区间。')
figure('source_comparison',500)
para('variation 类别间 SDI／decoding 在不同版本与分期中总体向下；Polar net RSA 和 Polar 相对 Grating 的角度变化总体向上。global 类别间变化较小，canonical 下部分 stage 规则可反向。native-union 与共同窗-union 的方向更一致，但均不能因此升级为已证实的显著学习效应。',True)
para('原 prototype 有额外 U 类别。统一到 G/H/R/S/T 改变了 global 类别间比较的 estimand 和 pair 数；不能把原 70 对与新五类别结果直接当作同一组数据。中间候选 stored_cross_common_labels 单独记录了这个类别范围影响。',True)

page();heading('5. 主候选的五个 session 不是同一个故事')
figure('primary_session_effects',500)
para('M1 的类别内 SDI 平均提高，而 M2 平均为负；LOO SVM 的 session 效应也有正有负。因此一个正的 pooled bar 不能代替两个动物分别成立的结论。原图 post-day pooling 还会让 recording 较多的 session 获得更大权重。')
para('主表采用两个动物等权的描述性平均；完整 summary.csv 同时给出按 session 汇总、各动物汇总和两侧精确 sign-flip 检验。这些都是探索性比较，未预注册、未作模型选择确认，不应仅从很多 sensitivity 组合中挑一个 p<0.05 的结果发表。',True)

page();heading('6. RMS／扣噪几何与 MFT 必须分开')
figure('geometry_comparison',500)
para('直接 RMS 半径测量全部 trial cloud 的离散程度；扣噪 signal radius 估计 condition mean 的离散程度。半径除以 sqrt(unit 数)，避免总范数随 neuron 数机械增长。扣噪维度为 signal covariance 的 participation ratio，负 eigenvalue 做 PSD 截断。它们不等于 MFT 的归一化半径、维度或容量。')
para('variation 的 DiD 为 T 类别相对 U 类别。global 直接几何在 s1/s3/s5 用 R/S/T 相对 G，在 s2/s4 用 S 相对 G/H/R/T 的平均；Fig5 的表征指标则在所有 session 用 Polar R/S/T 相对 Grating G。两类问题的对照定义不同，已明确保存，不能互换解释。',True)
para('共同窗主流程中，variation signal radius 的平均 DiD 在 M1 为正、M2 为负；global signal radius 在两动物都为负。trained-condition 距离增加也主要来自 M1。不能概括为两个动物中所有 manifold 都一致扩张。',True)

page();heading('7. 分期、统计与 neuron 数控制')
para('比较了 first/last thirds、first/last halves、first/last recording 和同一 session 内跨两种刺激的 calendar thirds。calendar thirds 下，s2 和 s4 的 variation 没有足够早期记录，variation 只剩三个 session；强行让两种刺激共享同一时间段会改变纳入范围，不能默默补齐。')
table([['规则','解释','限制'],['按 modality 的 thirds','相同 rank 分法，保留五个 session','不保证相同训练 exposure'],['按 session 的 calendar thirds','两种刺激用相同绝对日期边界','variation 只剩三个 session'],['真实训练／exposure 分期','由独立训练记录决定','缺少已验证的逐 recording 对应表']],[135,195,170])
para('找到了手工 exposure 计数函数，但它没有提供经过核对的完整日期映射。例如一个主 session 的 exposure 函数有四行，而当前 variation FR 只有三个 recording。报告没有把它自动当作完整训练日志。training_timeline_template.csv 列出需要补的原始证据字段。')
para('原 Fig3E/F 使用类别对作样本，Fig4E/F 的 SDI 与 SVM 还用了不同分期。只保持旧数据而改为 session 汇总：Fig3E 的双侧 t p≈0.159，Fig3F≈0.125，Fig4E≈0.215，Fig4F≈0.0025。最后一个仍依赖小样本参数检验假设；四者不能统一宣称保持原星号。')
para('对于五个 session，两侧精确 sign-flip 的最小 p 为 2/32=0.0625；按两个动物检验时最小 p 为 0.5。这是样本数与检验分辨率的限制，不代表真实效应一定为零。session 还嵌套在动物内，因此不能据这些 p 推断动物总体。原 Fig4B-D 的 stale-p 标星问题也不能通过选数据“复现”成正确统计。')
para(f'另做了 neuron 数敏感性：每个 recording 随机保留 {counts["units_per_recording"]} 个历史 unit，固定种子 {counts["seeds"]}，三次抽样。主分析仍保留全部 unit，未进行 SNR 筛选。三个抽样只用于数量控制，不是完整 Monte Carlo 置信区间。')
fixed='fixed_count_common240__all__thirds';ft=[['指标','全部 unit (session mean)','固定 unit 数 (session mean)']]
for key,title in [('var_within_sdi_DiD','类别内 SDI DiD (%)'),('var_between_svm','variation decoding (pp)'),('global_between_svm','global decoding (pp)'),('global_RSA_Polar_net','Polar net RSA')]:
 a=look[(primary,key,'pooled')];b=look[(fixed,key,'pooled')];ft.append([title,fmt(a['mean']),fmt(b['mean'])])
table(ft,[220,140,140])

page();heading('8. manifold 容量的验证范围与交付内容')
para('本次已完成全 recording 的 SDI、类别间五折 SVM、直接几何、orthogonality 和 RSA；canonical 与共同窗候选还完成了全 variation recording 的原 LOO SVM。没有把直接几何冒充 MFT。历史日级 MFT 缓存另作分期诊断，源 FR 的 hash lineage 未证明，不能把它算作统一数据下的容量重建。')
benchmark=ROOT/'manifold_benchmark.json'
if benchmark.exists():
 b=read('manifold_benchmark.json');para(f'另外完成 s1 variation 041116 的 1,000-vector MFT benchmark，耗时 {b["elapsed_seconds"]/60:.1f} 分钟。该单 recording 证明估计器可运行；它不能验证 52 个 recording 的容量趋势或显著性。')
else:para('全量 1,000-vector MFT 尚未验证；本报告不据此宣称容量结论得到复现。')
para('类别间五折 decoding 的训练与测试可包含同一 stimulus condition 的不同 trial，因此它测量 trial 层面的类别区分，不能直接解释为未见 stimulus 的泛化。Hyperbolic 只有两个相关 session，其相对 RSA 分母也可能不稳定，应单独标为探索性结果。',True)
para('所以，最可复现的一套输入与代码已经建立，若目标是严格统一数据与计算，应使用共同窗 union 候选并把 native 窗作为敏感性；若目标是复现旧图全部星号，当前证据不支持这一承诺。若真实训练日志最终指定的 baseline 与 earliest-FR 不同，需按该独立证据重新分期，再判断学习结论。')
table([['文件／目录','用途'],['unified/inputs/union_common240','统一共同窗 FR，保留原历史 unit 并集'],['unified/inputs/union_native','同样 unit 的 native 窗对照'],['unified/daily','各候选的 recording 级结果'],['summary.csv / session_effects.json','完整统计对照与 session 效应'],['stage_manifest.json','每个比较实际使用的日期'],['raw_sources.json / raw_rebuild_checks.json','原始 CDT 来源与重建误差'],['fig2_baseline_*.json','baseline-only 对照'],['unified/figures/*.pdf','比较图和几何／session 图'],['unified/run_unified.m','独立运行入口']],[230,270])
para('重跑：在 MATLAB 中 addpath(\'unified\'); run_unified()。若需要从 raw 重建，传入本地 CDT 根目录。默认复用已保存的日级结果；要强制重算，应先移走对应 daily 目录。原始 CDT 未上传，处理后的统一 FR 和结果均随分支保存。',True)
para('main 未改动。报告、统一 FR、代码和完整对照保存在 unified-analysis 分支。部分原图为示意／PCA 展示，未在此逐像素重画；本报告检验的是分析输入、效应与统计范围，不能视为所有 Fig1-5 panel 的完整版式复刻。',True)

def footer(c,doc):
 c.setFont(FONT,8);c.setFillColor(colors.HexColor('#637789'));c.drawString(42,24,'统一分析与结果一致性 | 2026-10-06');c.drawRightString(A4[0]-42,24,str(doc.page))
doc=SimpleDocTemplate(str(OUT/'unified_analysis_report.pdf'),pagesize=A4,rightMargin=42,leftMargin=42,topMargin=40,bottomMargin=42,title='Unified analysis consistency report',author='Manifold transformation analysis')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
(OUT/'unified_analysis_report.md').write_text('# 统一分析与结果一致性报告\n\n2026-10-06\n'+'\n'.join(md))
print(OUT/'unified_analysis_report.pdf')

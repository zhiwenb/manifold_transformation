from pathlib import Path
import csv,json,shutil,statistics
from PIL import Image,ImageOps,ImageDraw
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
r=Path(__file__).resolve().parents[1];a=r/'reference_check';a.mkdir(exist_ok=True)
import subprocess
render_root=r/'outputs_reference'/'audit_reference_images';render_root.mkdir(parents=True,exist_ok=True)
for i in range(1,6):
    subprocess.run(['pdftoppm','-scale-to','1600','-singlefile','-png',str(r/'references'/f'fig{i}.pdf'),str(render_root/f'fig{i}')],check=True)
# Boxes are normalized coordinates in the supplied reference pages. They are
# comparison crops only, never substituted for independently generated results.
refs={i:render_root/f'fig{i}.png' for i in range(1,6)}
boxes={
'fig2D':(2,(.04,.245,.96,.665)), 'fig2E':(2,(.025,.675,.292,.95)), 'fig2F':(2,(.30,.675,.574,.95)),
'fig3A':(3,(.025,.045,.54,.51)), 'fig3BCD':(3,(0,0,1,1)), 'fig3E':(3,(.36,.51,.647,.975)), 'fig3F':(3,(.67,.51,.971,.975)),
'fig4BCD':(4,(.70,.00,.995,.97)), 'fig4E':(4,(.005,.52,.343,.93)), 'fig4F':(4,(.385,.52,.682,.947)),
'fig5B':(5,(.02,.555,.465,.98)), 'fig5C':(5,(.495,.55,.955,.98))}
notes={
'fig2D':'Cloud geometry and selected recording agree visually. Standalone layout differs from the assembled page.',
'fig2E':'SDI/SVM bars and stars agree. Corrected panel label C to E.',
'fig2F':'Distance/variance bars and stars agree. Corrected panel label D to F.',
'fig3A':'Cloud geometry agrees. Stimulus thumbnails and full-page montage are not rebuilt.',
'fig3BCD':'Original variation tables: means and star pattern agree. Combined alpha/R/D output differs from page placement.',
'fig3E':'Original variation decoding: paired distribution and *** agree.',
'fig3F':'Restored Early/Late violin. Original script skips mismatched s2 pairs: 46 observations; *** agrees.',
'fig4BCD':'Bar means reproduce mixed inputs. WARNING: original stars share a stale p value; reference significance is not validated.',
'fig4E':'Restored original prototype decoding cache; ** and distribution agree.',
'fig4F':'Restored original global SDI cache; *** and distribution agree.',
'fig5B':'Means and Polar/Grating star agree. Errorbar/point structure matches visually.',
'fig5C':'Means and *, ns, ***, ns agree. Individual dot jitter is not pixel-identical.'}
manifest=[];doc=canvas.Canvas(str(a/'comparison.pdf'),pagesize=(1000,740))
for panel,(num,box) in boxes.items():
 orig=Image.open(refs[num]).convert('RGB');w,h=orig.size;crop=orig.crop(tuple(round(v*(w if k%2==0 else h)) for k,v in enumerate(box)))
 out=r/'outputs_reference'/panel/'panel_1.png';assert out.exists(),out
 dst=r/'results'/f'{panel}.png';shutil.copy2(out,dst);shutil.copy2(out.with_suffix('.pdf'),r/'results'/f'{panel}.pdf')
 crop.save(a/f'{panel}_reference.png');new=Image.open(out).convert('RGB')
 board=Image.new('RGB',(1200,720),'white');draw=ImageDraw.Draw(board);draw.text((20,12),panel+' | supplied reference (left) / regenerated (right)',fill='black')
 for x,im in [(20,crop),(620,new)]:
  im.thumbnail((560,640));board.paste(im,(x+(560-im.width)//2,55+(640-im.height)//2))
 board.save(a/f'{panel}_comparison.png')
 doc.setFont('Helvetica-Bold',16);doc.drawString(30,710,panel+' - Reference comparison')
 doc.setFont('Helvetica',10)
 # Wrap findings onto two lines as needed.
 import textwrap
 for j,line in enumerate(textwrap.wrap(notes[panel],145)):doc.drawString(30,690-j*13,line)
 doc.drawImage(ImageReader(board),20,50,width=960,height=576)
 doc.setFont('Helvetica',9);doc.drawString(30,25,'Visual/source audit. Not a claim of exact pixel or full-page assembly identity.');doc.showPage()
 manifest.append({'panel':panel,'status':('visual_match_with_invalid_star_logic' if panel=='fig4BCD' else 'visually_consistent_after_correction'),'scope':'analytical_content; not pixel-identical full-page layout','finding':notes[panel]})
for num,letters in [(1,'ABCD'),(2,'ABCG'),(4,'A'),(5,'A')]:
 for letter in letters:manifest.append({'panel':f'fig{num}{letter}','status':'not_regenerated','scope':'reference PDF only','finding':'Complete generation source has not been established; no equality claim.'})
doc.setFont('Helvetica-Bold',16);doc.drawString(30,710,'Panels not independently regenerated')
doc.setFont('Helvetica',12);y=665
for row in manifest:
 if row['status']=='not_regenerated':doc.drawString(35,y,row['panel']+' - reference PDF only; complete source not found.');y-=28
doc.setFont('Helvetica',11);doc.drawString(30,110,'The supplied figure pages combine standalone analyses, schematic assets and manual layout.')
doc.drawString(30,88,'The analytical panels above were checked; complete Fig1-5 page identity is not established.')
doc.save();(a/'coverage.json').write_text(json.dumps(manifest,indent=2)+'\n')
with (a/'coverage.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(manifest[0]));wr.writeheader();wr.writerows(manifest)
html='<!doctype html><meta charset="utf-8"><title>Reference figure audit</title><style>body{font:16px sans-serif;max-width:1200px;margin:30px auto}img{max-width:100%}article{margin:35px 0;border-top:1px solid #ddd}p{line-height:1.5}</style><h1>Reference figure audit</h1><p>Left: original PDF crop. Right: regenerated analytical panel. Matching scientific content does not establish pixel-identical full-page reproduction.</p>'
for row in manifest:
 html+=f'<article><h2>{row["panel"]}</h2><p>{row["finding"]}</p>'
 if row['status']!='not_regenerated':html+=f'<img src="{row["panel"]}_comparison.png">'
 html+='</article>'
(a/'index.html').write_text(html)
# Record means from actual regenerated output tables, rather than reading bar
# heights from the rendered PDF. Reference raw plotting vectors are not supplied.
summary={}
for panel in ['fig3E','fig3F','fig4E','fig4F']:
 name='SDI_paired_per_pair_data.csv' if panel.endswith('F') else 'paired_per_pair_data.csv'
 rows=list(csv.DictReader((r/'outputs_reference'/panel/name).open()))
 early='early_SDI' if panel.endswith('F') else 'early_acc';late='late_SDI' if panel.endswith('F') else 'late_acc'
 if early not in rows[0]:early='early';late='late'
 summary[panel]={'n_pairs':len(rows),'early_mean':statistics.mean(float(x[early]) for x in rows),'late_mean':statistics.mean(float(x[late]) for x in rows),'mean_delta':statistics.mean(float(x['delta']) for x in rows)}
(a/'regenerated_means.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Audit PDF, comparison images, HTML, coverage table and means saved to',a)

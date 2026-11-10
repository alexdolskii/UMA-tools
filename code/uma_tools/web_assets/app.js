/* Offline views of saved measurements. No filtering, fitting or tests. */
(() => {
  'use strict';
  const index = window.UMA_INDEX;
  const $ = id => document.getElementById(id);
  let current, selection = 0;
  const el = (tag, text, cls) => {
    const node = document.createElement(tag);
    if (text !== undefined && text !== null) node.textContent = String(text);
    if (cls) node.className = cls;
    return node;
  };
  const append = (parent, ...nodes) => nodes.forEach(n => parent.append(n));
  const button = (text, action, cls) => {
    const node = el('button', text, cls); node.type = 'button';
    node.addEventListener('click', action); return node;
  };
  const section = (id, title, note) => {
    const node = $(id); node.replaceChildren(el('h2', title));
    if (note) node.append(el('p', note)); return node;
  };
  const download = (file, label) => {
    const link = el('a', label || file.name); link.href = file.path;
    link.download = file.name.split('/').pop(); return link;
  };
  const num = x => x === '' || x === null ? NaN : Number(x);
  const yes = x => String(x).toLowerCase() === 'true' || x === '1';
  const dataRows = (run, table) => run.tables[table]?.rows || [];
  const fmt = value => Number.isFinite(Number(value)) ? Number(value).toLocaleString('en-US',{maximumSignificantDigits:4}) : String(value);
  const svgNS = 'http://www.w3.org/2000/svg';
  const svg = (tag, attrs = {}, text) => {
    const n = document.createElementNS(svgNS, tag);
    Object.entries(attrs).forEach(([k,v]) => n.setAttribute(k, String(v)));
    if (text !== undefined) n.textContent = String(text); return n;
  };
  const loadRun = id => new Promise((resolve, reject) => {
    if (window.UMA_RUNS[id]) return resolve(window.UMA_RUNS[id]);
    const entry = index.runs.find(r => r.id === id);
    const script = el('script'); script.src = entry.script;
    script.onload = () => resolve(window.UMA_RUNS[id]);
    script.onerror = () => reject(new Error('Missing run data. Copy the entire report folder.'));
    document.head.append(script);
  });

  function overview() {
    const root = $('overview'); root.replaceChildren(el('div', 'MATRIX ANALYSIS / SAVED EXPERIMENT', 'eyebrow'));
    append(root, el('h1', index.experiment), el('p', 'Explore completed reports, compare saved masks, and trace every point back to its image. Each run keeps its own parameters and population.'));
    const metrics = el('div', null, 'metrics');
    for (const [value,label] of [[index.runs.length,'Recognized saved runs'], [index.runs.filter(r=>r.kind==='Report').length,'UMA report runs'], [index.runs.filter(r=>r.kind==='Area').length,'Area measurement runs'], [index.state.warning_count,'Export warnings']]) {
      const item = el('div', null, 'metric'); append(item, el('strong', value), el('span', label)); metrics.append(item);
    }
    append(root, metrics, el('span', index.state.status, 'pill'), el('span', 'UMA-tools '+index.state.package_version, 'pill'));
    const all = el('details'); all.append(el('summary', 'All runs and export diagnostics'));
    for (const run of index.runs) all.append(button(`${run.kind} · ${run.status} · ${run.name}`, () => {$('run').value=run.id; selectRun(run.id);}));
    if (index.state.warnings.length) {
      const warnings = el('pre', index.state.warnings.join('\n')); all.append(warnings);
    }
    root.append(all);
  }

  function marker(shape, x, y, color, edge) {
    const a = {fill:color, stroke:edge, 'stroke-width':edge==='#d62728'?2:1};
    if (shape==='s') return svg('rect', {...a,x:x-4,y:y-4,width:8,height:8});
    if (shape==='D' || shape==='d') return svg('path',{...a,d:`M${x},${y-5}l5,5 -5,5 -5,-5Z`});
    if (['^','v','<','>'].includes(shape)) {
      const rotation = {'^':0,'>':90,'v':180,'<':270}[shape];
      return svg('path',{...a,d:`M${x},${y-5}l5,9 -10,0Z`,transform:`rotate(${rotation} ${x} ${y})`});
    }
    if (shape==='P' || shape==='+') return svg('path',{...a,d:`M${x-2},${y-5}h4v3h3v4h-3v3h-4v-3h-3v-4h3Z`});
    if (shape==='X' || shape==='x') return svg('path',{...a,d:`M${x-2},${y-5}h4v3h3v4h-3v3h-4v-3h-3v-4h3Z`,transform:`rotate(45 ${x} ${y})`});
    if (shape==='*') {
      const points=Array.from({length:10},(_,i)=>{const a=i*Math.PI/5-Math.PI/2,r=i%2?2.5:6;return `${x+r*Math.cos(a)},${y+r*Math.sin(a)}`;});
      return svg('polygon',{...a,points:points.join(' ')});
    }
    if (['p','h','8'].includes(shape)) {
      const sides={'p':5,'h':6,'8':8}[shape];
      const points=Array.from({length:sides},(_,i)=>{const angle=i*2*Math.PI/sides-Math.PI/2;return `${x+5*Math.cos(angle)},${y+5*Math.sin(angle)}`;});
      return svg('polygon',{...a,points:points.join(' ')});
    }
    if (/^\$\d+\$$/.test(shape)) return svg('text',{...a,x,y:y+4,'text-anchor':'middle','font-size':11,'font-weight':'bold','stroke-width':.4,'paint-order':'stroke'},shape.slice(1,-1));
    return svg('circle',{...a,cx:x,cy:y,r:4.2});
  }

  function drawPanel(plot, panel, run) {
    const n = panel.groups.length, width = Math.max(420, n*145+100), height = 500;
    const labels=Object.fromEntries(panel.groups.map(group=>{
      const lines=[];
      for(const line of String(panel.labels?.[group]||group).split('\n')) {
        let buffer='';
        for(const word of line.split(/\s+/)) {
          if(buffer.length+word.length>21 && buffer){lines.push(buffer);buffer='';}
          let rest=word;while(rest.length>21){if(buffer){lines.push(buffer);buffer='';}lines.push(rest.slice(0,21));rest=rest.slice(21);}
          buffer+=(buffer?' ':'')+rest;
        }
        if(buffer)lines.push(buffer);
      }
      return [group,lines];
    }));
    const comparisons = (plot.statistical_comparisons || []).filter(c => c.Panel===panel.id && c.Bracket_Level!==null);
    const levels = Math.max(0,...comparisons.map(c=>Number(c.Bracket_Level)+1));
    const top = 48+levels*23, bottom = 355+levels*23;
    const fullHeight=Math.max(height+levels*23,bottom+80+Math.max(1,...Object.values(labels).map(a=>a.length))*14);
    const chart = svg('svg',{viewBox:`0 0 ${width} ${fullHeight}`,width,height:fullHeight,role:'img','aria-label':`${plot.title}: ${panel.title}`});
    const yMin=Number(plot.y_min), yMax=Number(plot.y_max);
    const x = v => 65+(v-.4)/(n+.2)*(width-95);
    const y = v => bottom-(v-yMin)/(yMax-yMin||1)*(bottom-top);
    chart.append(svg('text',{x:65,y:22,'font-size':14,'font-weight':'bold',fill:'#152a2a'},panel.title||panel.id));
    for(let i=0;i<=4;i++) {
      const value=yMin+(yMax-yMin)*i/4;
      chart.append(svg('line',{x1:65,x2:width-30,y1:y(value),y2:y(value),stroke:'#dce4df'}));
      chart.append(svg('text',{x:57,y:y(value)+4,'text-anchor':'end','font-size':11,fill:'#586966'},fmt(value)));
    }
    chart.append(svg('line',{x1:65,x2:width-30,y1:bottom,y2:bottom,stroke:'#82958a'}));
    const cutoff = num(run.parameters.fn_area_threshold_percent);
    if(plot.name==='Fibronectin' && Number.isFinite(cutoff) && cutoff>=yMin && cutoff<=yMax) chart.append(svg('line',{x1:65,x2:width-30,y1:y(cutoff),y2:y(cutoff),stroke:'#d62728','stroke-dasharray':'5 5'}));
    panel.groups.forEach((group,i)=>{
      const pos=x(i+1), b=plot.boxes[group], style=plot.condition_styles[group];
      if(b) {
        chart.append(svg('line',{x1:pos,x2:pos,y1:y(b.whislo),y2:y(b.whishi),stroke:'#64748b'}));
        for(const val of [b.whislo,b.whishi]) chart.append(svg('line',{x1:pos-12,x2:pos+12,y1:y(val),y2:y(val),stroke:'#64748b'}));
        chart.append(svg('rect',{x:pos-27,y:y(b.q3),width:54,height:Math.max(.5,y(b.q1)-y(b.q3)),fill:style.color,stroke:'#334155'}));
        chart.append(svg('line',{x1:pos-27,x2:pos+27,y1:y(b.med),y2:y(b.med),stroke:style.median_color,'stroke-width':2}));
      }
      const words=labels[group];
      words.forEach((text,line)=>chart.append(svg('text',{x:pos,y:bottom+22+line*14,'text-anchor':'middle','font-size':11,fill:'#152a2a'},text)));
      chart.append(svg('text',{x:pos,y:bottom+30+words.length*14,'text-anchor':'middle','font-size':11,fill:'#586966'},`n=${plot.group_counts[group]} images · ${plot.group_well_counts[group]} wells`));
    });
    for(const p of plot.points.filter(p=>p.Panel===panel.id)) {
      const point=marker(p.Point_Marker,x(Number(p.Panel_X)),y(Number(p.Value)),p.Point_Color,yes(p.Red_Outline)?'#d62728':plot.point_edge_color||'#334155');
      const description=`${p.Group}, ${p.Well}: ${p.Value} ${p.Unit}. ${p.Image_ID}`;
      point.setAttribute('class','point');point.setAttribute('tabindex','0');point.setAttribute('role','button');point.setAttribute('aria-label',description);
      point.append(svg('title',{},description));
      const open=()=>showImage(run,run.observations.findIndex(o=>o.image===p.Image_ID));
      point.addEventListener('click',open);point.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();open();}});chart.append(point);
    }
    for(const c of comparisons) {
      const yy=top-14-Number(c.Bracket_Level)*23,left=x(Number(c.Bracket_Left)),right=x(Number(c.Bracket_Right));
      chart.append(svg('path',{d:`M${left},${yy+5}V${yy}H${right}V${yy+5}`,fill:'none',stroke:'#334155'}));
      const label=svg('text',{x:(left+right)/2,y:yy-4,'text-anchor':'middle','font-size':12},c.Annotation);
      label.append(svg('title',{},`${c.Treatment} vs ${c.Control}; adjusted p=${c.P_Holm}; ${c.Stats_Unit}`));chart.append(label);
    }
    return chart;
  }

  function results(run) {
    const root=section('results','Results','One point = one image. Colors, well shapes, populations, boxes and statistics are taken from the saved report. Click or focus a point to inspect its image.');
    if(!run.charts.length) {root.append(el('p','No validated interactive plots in this run. Available measurements, previews and original files are listed below.','warning'));return;}
    for(const plot of run.charts) {
      const card=el('article',null,'chart');append(card,el('h3',plot.title),el('p',`${plot.view} · ${plot.unit} · ${plot.point_count} images`));
      const panels=el('div',null,'panels');
      for(const panel of plot.panels){const block=el('div',null,'panel');block.append(drawPanel(plot,panel,run));panels.append(block);}
      append(card,panels,el('p',plot.caption,'caption'),el('p',plot.statistics_note,'caption'));
      const links=el('div',null,'downloads');
      for(const file of run.downloads.filter(f=>[plot.pdf_file,plot.png_file].some(p=>p && p.replaceAll('\\','/').split('/').pop()===f.name.split('/').pop()))) links.append(download(file));
      card.append(links);
      const legend=el('details');legend.append(el('summary','Well shapes and display rules'));
      legend.append(el('p','Gray points represent individual images; a red outline identifies low Fiber coverage in the unfiltered view. Each condition keeps the recorded box color. Well shapes are preserved after filtering.','caption'));
      const symbols=el('div',null,'downloads');
      Object.entries(plot.well_markers||{}).forEach(([well,shape])=>{
        const item=el('span'),icon=svg('svg',{width:22,height:22,viewBox:'0 0 22 22','aria-label':shape});icon.append(marker(shape,11,10,plot.point_color,plot.point_edge_color));append(item,icon,el('span',well));symbols.append(item);
      });legend.append(symbols);card.append(legend);root.append(card);
    }
  }

  function showTable(parent, table, run) {
    if(!table) return;
    let page=0,query='',sortKey='',reverse=false;
    const toolbar=el('div',null,'toolbar'),search=el('input');search.type='search';search.placeholder='Search this table';search.setAttribute('aria-label','Search this table');
    const count=el('span'),previous=button('Previous',()=>{page--;render();}),next=button('Next',()=>{page++;render();});
    append(toolbar,search,count,previous,next);const scroll=el('div',null,'scroll');append(parent,toolbar,scroll);
    function render() {
      let values=table.rows.filter(row=>!query||Object.values(row).some(v=>String(v).toLowerCase().includes(query)));
      if(sortKey) values=[...values].sort((a,b)=>{const va=num(a[sortKey]),vb=num(b[sortKey]);const difference=Number.isFinite(va)&&Number.isFinite(vb)?va-vb:String(a[sortKey]).localeCompare(String(b[sortKey]));return reverse?-difference:difference;});
      page=Math.max(0,Math.min(page,Math.ceil(values.length/50)-1));count.textContent=`${values.length} rows · page ${page+1}/${Math.max(1,Math.ceil(values.length/50))}`;
      previous.disabled=page===0;next.disabled=(page+1)*50>=values.length;
      const tableEl=el('table'),thead=el('thead'),head=el('tr');
      for(const column of table.columns){const th=el('th',column);th.tabIndex=0;const sort=()=>{reverse=sortKey===column?!reverse:false;sortKey=column;render();};th.addEventListener('click',sort);th.addEventListener('keydown',e=>{if(e.key==='Enter')sort();});head.append(th);}
      thead.append(head);tableEl.append(thead);const body=el('tbody');
      for(const row of values.slice(page*50,(page+1)*50)) {
        const tr=el('tr');table.columns.forEach(c=>tr.append(el('td',row[c]??'')));
        const image=row.Image_ID||row.Source_File_Name||row.File_Name;
        const index=run.observations.findIndex(o=>o.image===image);
        if(index>=0||row.Well){tr.className='data-row';tr.tabIndex=0;const open=()=>index>=0?showImage(run,index):showWell(run,row.Well);tr.addEventListener('click',open);tr.addEventListener('keydown',e=>{if(e.key==='Enter')open();});}
        body.append(tr);
      }
      tableEl.append(body);scroll.replaceChildren(tableEl);
    }
    search.addEventListener('input',()=>{query=search.value.trim().toLowerCase();page=0;render();});render();
  }

  function tables(run) {
    const root=section('tables','Measurements & statistics','Search and sorting change the table view only. Click an image row, or a well-mean row to inspect its constituent images. Original files retain their full precision.');
    const chooser=el('select'),area=el('div');chooser.setAttribute('aria-label','Select measurement table');
    Object.keys(run.tables).sort().forEach(name=>{const option=el('option',`${name} (${run.tables[name].row_count} rows)`);option.value=name;chooser.append(option);});
    const initial=['merged_data','plot_data','statistics'].find(n=>run.tables[n]);if(initial)chooser.value=initial;
    const render=()=>{area.replaceChildren();showTable(area,run.tables[chooser.value],run);};chooser.addEventListener('change',render);append(root,chooser,area);render();
  }

  function plate(run) {
    const root=section('plate','Plate & channels','Saved scheme for this run. The current file beside the raw images is not reread. Click a well to see its images.');
    const design=run.documents['plate_design.json'];
    if(!design?.Wells){root.append(el('p','No structured plate snapshot in this run. Download the saved workbook if available.'));return;}
    const grid=el('div',null,'plate-grid');grid.append(el('span',''));for(let c=1;c<=12;c++)grid.append(el('span',c,'plate-head'));
    for(const letter of 'ABCDEFGH') {
      grid.append(el('span',letter,'plate-head'));
      for(let c=1;c<=12;c++) {
        const well=letter+String(c).padStart(2,'0'),record=design.Wells[well];
        const dot=button(well,()=>showWell(run,well),'well');dot.title=record?JSON.stringify(record):'Unassigned well';
        const condition=record?.condition||record?.group||record?.Group;
        const style=run.charts[0]?.condition_styles[condition];
        if(style)dot.style.background=style.color;
        if(!record)dot.disabled=true;grid.append(dot);
      }
    }
    const wrapper=el('div',null,'scroll');wrapper.append(grid);root.append(wrapper);
    const detail=el('details');append(detail,el('summary','Saved well assignments and plot layout'),el('pre',JSON.stringify(design,null,2)));root.append(detail);
  }

  function quality(run) {
    const root=section('qc','Quality control','Measurement status and web-export status are separate. Missing files are unavailable, never zero measurements.');
    append(root,el('span',`Source: ${run.status}`,'pill'),el('span',`Export: ${index.state.status}`,'pill'));
    for(const warning of run.warnings)root.append(el('p',warning,'warning'));
    const tables=['qc','group_filter_counts','processing_exclusions','fia_image_check','image_check','errors'];
    for(const name of tables) if(run.tables[name]){const details=el('details');details.append(el('summary',`${name} (${run.tables[name].rows.length} rows)`));showTable(details,run.tables[name],run);root.append(details);}
  }

  function previewFigure(file) {
    const asset=index.assets[file.asset],fig=el('figure');fig.append(el('h3',file.label));
    if(asset?.status==='AVAILABLE') {
      const image=el('img');image.src=asset.preview;image.loading='lazy';image.alt=file.label;fig.append(image);
      const range=asset.display_range?` · shared display range ${asset.display_range.join('–')}`:'';
      fig.append(el('p',`${asset.width} × ${asset.height} · ${asset.dtype}${range}`,'caption'));
    } else fig.append(el('p',asset?.error||'Not available','warning'));
    fig.append(el('p',file.source,'paths'));return fig;
  }

  function gallery(run) {
    const root=section('gallery','Image examples','Up to five reproducible examples per condition and saved population (seed 0). Every linked image remains accessible through plots and tables.');
    const chooser=el('select');chooser.setAttribute('aria-label','Gallery condition');const all=el('option','All conditions');all.value='';chooser.append(all);
    [...new Set(run.gallery.map(g=>g.group))].sort().forEach(group=>{const option=el('option',group||'Unannotated');option.value=group;chooser.append(option);});
    const area=el('div');append(root,chooser,area);
    const render=()=>{
      area.replaceChildren();
      for(const sample of run.gallery.filter(g=>!chooser.value||g.group===chooser.value)) {
        area.append(el('h3',`${sample.group||'Unannotated'} · ${sample.population} · ${sample.indices.length} of ${sample.total}`));const cards=el('div',null,'cards');
        for(const id of sample.indices) {
          const observation=run.observations[id],card=el('article',null,'card');
          const file=observation.files.find(f=>f.label==='Area SUM projection' && index.assets[f.asset]?.preview)||observation.files.find(f=>index.assets[f.asset]?.preview);
          if(file){const image=el('img');image.src=index.assets[file.asset].preview;image.alt=file.label;image.loading='lazy';card.append(image);}else card.append(el('p','Preview unavailable','warning'));
          append(card,el('p',observation.image,'caption'),button('Inspect image',()=>showImage(run,id)));cards.append(card);
        }
        area.append(cards);
      }
    };chooser.addEventListener('change',render);render();
  }

  function showWell(run,well) {
    const content=$('image-content');content.replaceChildren(el('h2',`Well ${well}`));
    const records=run.observations.map((o,i)=>({o,i})).filter(r=>r.o.well===well);
    content.append(el('p',`${records.length} available images in this saved run. For well statistics, the saved well mean and constituent images are available in Measurements & statistics.`));
    records.forEach(({o,i})=>content.append(button(o.image,()=>showImage(run,i))));
    const assignment=run.documents['plate_design.json']?.Wells?.[well];if(assignment)content.append(el('pre',JSON.stringify(assignment,null,2)));
    if(!$('image-dialog').open)$('image-dialog').showModal();
  }

  function showImage(run,id) {
    const observation=run.observations[id];if(!observation)return;
    const content=$('image-content');content.replaceChildren(el('h2',observation.image));
    append(content,el('p',`${observation.group} · ${observation.well} · ${run.name}`),el('p',observation.source||'Original source path not recorded','paths'));
    if(observation.reason)content.append(el('p',observation.reason,'warning'));
    observation.warnings.forEach(w=>content.append(el('p',w,'warning')));
    const images=el('div',null,'image-grid');observation.files.forEach(f=>images.append(previewFigure(f)));content.append(images);
    const matches=index.comparisons.filter(c=>c.key===observation.comparison_key && !(c.run===run.id && c.index===id) && c.files.length);
    if(matches.length) {
      const compare=el('section');compare.append(el('h3','Compare saved Area masks'));
      compare.append(el('p','Same original image SHA256, Fiber marker and source channel. Each mask keeps its own intensity thresholds.'));
      const chooser=el('select');chooser.setAttribute('aria-label','Comparison Area run');
      matches.forEach((match,i)=>{const item=index.runs.find(r=>r.id===match.run);const option=el('option',`${match.threshold.Threshold_Lower} to ${match.threshold.Threshold_Upper||match.threshold.Effective_Threshold_Upper} · ${item.name}`);option.value=i;chooser.append(option);});
      const pair=el('div',null,'compare');
      const render=()=>{pair.replaceChildren();const original=observation.files.find(f=>f.label==='Area mask');if(original){const fig=previewFigure(original);fig.prepend(el('p',JSON.stringify(observation.threshold),'caption'));pair.append(fig);}const alternative=matches[Number(chooser.value)];if(alternative.files[0]){const fig=previewFigure(alternative.files[0]);fig.prepend(el('p',JSON.stringify(alternative.threshold),'caption'));pair.append(fig);}};
      chooser.addEventListener('change',render);append(compare,chooser,pair);render();content.append(compare);
    }
    const detail=el('details');detail.open=true;append(detail,el('summary','Saved measurements & metadata'),el('pre',JSON.stringify(observation.row,null,2)));content.append(detail);
    if(!$('image-dialog').open)$('image-dialog').showModal();
  }

  function provenance(run) {
    const root=section('provenance','Methods & provenance','This page reads recorded results. It does not reopen ND2 stacks, rerun segmentation, change thresholds, recompute statistics or pool repeated runs. Scientific units remain those recorded by UMA.');
    append(root,el('p','Signal previews use a shared linear display range within the same measurement run, marker/channel, projection type, Z count (when recorded) and data type. Binary masks use black/white. Previews are reduced to at most 1024 pixels on their longest side; measurements and TIFFs are unchanged.'),el('pre',JSON.stringify({source:run.source,route:run.route,parameters:run.parameters},null,2)));
    const files=el('details');files.append(el('summary',`Download original documents (${run.downloads.length})`));
    const list=el('div',null,'downloads');run.downloads.forEach(file=>list.append(download(file)));files.append(list);root.append(files);
    root.append(download({path:'manifest.json',name:'manifest.json'},'Export manifest & source paths'));
  }

  async function selectRun(id) {
    const request=++selection;$('loading').hidden=false;$('report-content').hidden=true;
    try {
      const run=await loadRun(id);if(request!==selection)return;current=run;
      $('run-context').textContent=`${run.kind} · ${run.status} · ${run.observations.length} images · ${run.route.marker_name||'Marker not recorded'}${run.route.source_channel?' / C'+run.route.source_channel:''}`;
      results(run);plate(run);quality(run);gallery(run);tables(run);provenance(run);
      $('loading').hidden=true;$('report-content').hidden=false;
    } catch(error){$('loading').textContent=error.message;}
  }

  if(!index){$('loading').textContent='Report data missing. Keep the data directory beside index.html.';return;}
  overview();
  index.runs.forEach(run=>{const option=el('option',`${run.kind} · ${run.status} · ${run.name}`);option.value=run.id;$('run').append(option);});
  const initial=[...index.runs].reverse().find(r=>r.kind==='Report' && ['SUCCESS','PARTIAL'].includes(r.status))||index.runs[0];
  $('run').value=initial.id;$('run').addEventListener('change',()=>selectRun($('run').value));selectRun(initial.id);
  $('image-dialog').querySelector('.close').addEventListener('click',()=>$('image-dialog').close());
  try {$('annotation').value=JSON.parse($('annotation-data').textContent);}catch(_){$('annotation').value='';}
  $('save-annotation').addEventListener('click',()=>{
    const clone=document.documentElement.cloneNode(true);
    clone.querySelector('#annotation-data').textContent=JSON.stringify($('annotation').value).replaceAll('<','\\u003c');
    clone.querySelector('#image-dialog').removeAttribute('open');
    clone.querySelectorAll('head script[src]').forEach(n=>n.remove());
    const blob=new Blob(['<!doctype html>\n'+clone.outerHTML],{type:'text/html;charset=utf-8'});
    const url=URL.createObjectURL(blob),link=el('a');link.href=url;link.download='UMA_annotated_report.html';link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
    $('save-note').textContent='Save the downloaded HTML beside index.html. Keep the entire folder together.';
  });
})();

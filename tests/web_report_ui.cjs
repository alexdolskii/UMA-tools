/* DOM regressions; these do not replace a browser rendering check. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const {JSDOM, VirtualConsole, ResourceLoader} = require('jsdom');

class LocalOnly extends ResourceLoader {
  fetch(url, options) {
    assert.ok(url.startsWith('file:'), `Unexpected network request: ${url}`);
    return super.fetch(url, options);
  }
}

(async () => {
  const root=process.argv[2], errors=[];
  const virtualConsole=new VirtualConsole();
  virtualConsole.on('jsdomError', error=>errors.push(String(error)));
  let blob;
  const dom=await JSDOM.fromFile(path.join(root,'index.html'),{
    resources:new LocalOnly(),runScripts:'dangerously',virtualConsole,
    beforeParse(window) {
      window.HTMLDialogElement.prototype.showModal=function(){this.open=true;};
      window.HTMLDialogElement.prototype.close=function(){this.open=false;};
      window.URL.createObjectURL=value=>{blob=value;return 'blob:annotation';};
      window.URL.revokeObjectURL=()=>{};
      window.HTMLAnchorElement.prototype.click=function(){};
    }
  });
  const {window}=dom, doc=window.document;
  async function ready() {
    for(let i=0;i<100;i++) {
      if(doc.querySelectorAll('.point').length && !doc.querySelector('#report-content').hidden)return;
      await new Promise(resolve=>setTimeout(resolve,30));
    }
    throw new Error('Run did not render: '+errors.join('\n'));
  }
  await ready();
  const runId=doc.querySelector('#run').value, run=window.UMA_RUNS[runId];
  assert.equal(doc.querySelectorAll('.chart').length,run.charts.length);
  assert.equal(doc.querySelectorAll('.point').length,run.charts.reduce((n,p)=>n+p.points.length,0));
  assert.match(doc.querySelector('#results').textContent,/n=0/);
  assert.match(doc.querySelector('#results').textContent,/ns|Not tested|\*/);
  const first=doc.querySelector('.point');
  first.dispatchEvent(new window.MouseEvent('click',{bubbles:true}));
  assert.equal(doc.querySelector('#image-dialog').open,true);
  assert.match(doc.querySelector('#image-content').textContent,/Compare saved Area masks/);
  assert.equal(doc.querySelectorAll('.compare img').length,2);
  assert.ok(doc.querySelectorAll('#image-content img').length>=7);
  const compare=doc.querySelector('select[aria-label="Comparison Area run"]');
  compare.value=compare.options.length-1;compare.dispatchEvent(new window.Event('change'));
  assert.equal(doc.querySelectorAll('.compare img').length,2);
  doc.querySelector('.close').click();assert.equal(doc.querySelector('#image-dialog').open,false);
  const search=doc.querySelector('#tables input');search.value='nonexistent observation';search.dispatchEvent(new window.Event('input'));
  assert.equal(doc.querySelectorAll('#tables tbody tr').length,0);
  search.value='';search.dispatchEvent(new window.Event('input'));assert.ok(doc.querySelectorAll('#tables tbody tr').length);
  doc.querySelector('.well:not(:disabled)').click();assert.equal(doc.querySelector('#image-dialog').open,true);
  const annotation='Experiment note </script><script>window.injected=true</script> & <b>plain</b>';
  doc.querySelector('#annotation').value=annotation;doc.querySelector('#save-annotation').click();
  const text=await new Promise(resolve=>{const reader=new window.FileReader();reader.onload=()=>resolve(reader.result);reader.readAsText(blob);});
  assert.ok(text.includes('\\u003c/script>'));
  const parsed=new JSDOM(text);
  assert.equal(JSON.parse(parsed.window.document.querySelector('#annotation-data').textContent),annotation);
  assert.equal(window.injected,undefined);
  for(const file of run.downloads)assert.ok(fs.existsSync(path.join(root,file.path)));
  assert.deepEqual(errors,[]);
  dom.window.close();parsed.window.close();
  console.log('PASS: offline DOM, saved points, n=0, mask comparison, tables, wells, annotations, downloads');
})().catch(error=>{console.error(error);process.exit(1);});

ObjC.import('PDFKit'); ObjC.import('Foundation');
function run(argv){ var outdir=argv[0]; var n=0;
 for (var p of argv.slice(1)){ var d=$.PDFDocument.alloc.initWithURL($.NSURL.fileURLWithPath(p)); if(!d||d.isNil()) continue;
  var t=''; for (var i=0;i<d.pageCount;i++){ var s=d.pageAtIndex(i).string; t+='\n=== page '+(i+1)+' ===\n'+(s.isNil()?'[no text layer]':ObjC.unwrap(s)); }
  var name=p.split('/').pop().replace(/\.pdf$/,'.txt');
  $(t).writeToFileAtomicallyEncodingError(outdir+'/'+name,true,$.NSUTF8StringEncoding,null); n++; }
 return 'wrote '+n; }

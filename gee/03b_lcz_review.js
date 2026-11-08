// =====================================================================
// Step 4b. Review the LCZ candidates: accept or reject each 100 m square.
// Work through one class at a time. When a class is finished, click
// "Save this class": it exports the accepted squares as an asset (LCZ_ok_<code>).
// Then run the merge lines at the bottom (uncomment them) to build LCZ_training.
// =====================================================================
var PROJECT = 'projects/lunar-tine-512302-q0/assets/';
// Which set to review:
//   'LCZ_candidates'        GEE candidates (first batch)
//   'LCZ_candidates_extra'  GEE candidates (second batch)
//   'LCZ_samples_clean'     researcher-located site points from the CSV (each checked one by one)
var BATCH = 'LCZ_samples_clean';
var SUFFIX = BATCH === 'LCZ_candidates' ? '' : BATCH === 'LCZ_candidates_extra' ? '_extra' : '_csv';
var UTM = ee.Projection('EPSG:32630');
var metricsImg = ee.Image(PROJECT + 'LCZ_metrics_100m');
function num(v) { return ee.Number(ee.Algorithms.If(ee.Algorithms.ObjectType(v).equals('String'), ee.Number.parse(v), v)); }
var cands = BATCH === 'LCZ_samples_clean'
  ? ee.FeatureCollection(PROJECT + BATCH).map(function (f) {
      var pt = ee.Geometry(ee.Algorithms.If(f.get('longitude'),
        ee.Geometry.Point([num(f.get('longitude')), num(f.get('latitude'))]), f.geometry()));
      return ee.Feature(pt.buffer(50).bounds(1, UTM), {lcz: num(f.get('lcz_code')).int(),
        cand_id: f.get('sample_id'), site: f.get('site'), split: f.get('split')});
    })
  : ee.FeatureCollection(PROJECT + BATCH);
var NAMES = {'2': 'LCZ 2 Compact mid-rise', '3': 'LCZ 3 Compact low-rise', '5': 'LCZ 5 Open mid-rise',
  '6': 'LCZ 6 Open low-rise', '8': 'LCZ 8 Large low-rise', '9': 'LCZ 9 Sparsely built', '10': 'LCZ 10 Heavy industry', '11': 'LCZ A Dense trees', '12': 'LCZ B Scattered trees',
  '14': 'LCZ D Low plants', '16': 'LCZ F Bare soil', '17': 'LCZ G Water'};
var dry = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(cands.geometry())
  .filterDate('2020-01-01', '2023-01-01').filter(dry).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)).median();

// Convert to latitude/longitude first: the CSV squares are built in UTM (metres)
var all = cands.map(function (f) { return f.setGeometry(f.geometry().transform('EPSG:4326', 1)); })
               .toList(2000).getInfo();          // candidates as client-side GeoJSON
var cls = null, items = [], i = 0, decisions = {};

Map.setOptions('HYBRID');
var info = ui.Label('Choose a class to start.');
var metrics = ui.Label('');
var tally = ui.Label('');
var sv = ui.Label('Open Street View here', {color: 'blue'});
var earth = ui.Label('Open Google Earth 3D here', {color: 'blue'});
var showS2 = ui.Checkbox('Show Sentinel-2 2020-2022 (the image the classifier uses)', false, function (on) {
  Map.layers().forEach(function (l) { if (l.getName() === 'Sentinel-2 2020-2022') l.setShown(on); });
});

function centre(f) {
  var ring = f.geometry.coordinates[0], x = 0, y = 0;
  for (var j = 0; j < ring.length - 1; j++) { x += ring[j][0]; y += ring[j][1]; }
  return [y / (ring.length - 1), x / (ring.length - 1)];
}
function show() {
  if (!items.length) { info.setValue('No candidates for this class.'); return; }
  var f = items[i], p = f.properties, ll = centre(f);
  Map.layers().reset();
  Map.addLayer(s2, {bands: ['B4', 'B3', 'B2'], min: 0, max: 2500}, 'Sentinel-2 2020-2022', showS2.getValue());
  Map.addLayer(ee.FeatureCollection([ee.Feature(f)]).style({color: 'ff0000', fillColor: '00000000', width: 3}), {}, 'Candidate');
  Map.setCenter(ll[1], ll[0], 17);
  var d = decisions[p.cand_id];
  info.setValue(NAMES[cls] + ': candidate ' + (i + 1) + ' of ' + items.length +
                (d === true ? '  [ACCEPTED]' : d === false ? '  [REJECTED]' : ''));
  metrics.setValue('Measuring...');
  metricsImg.reduceRegion({reducer: ee.Reducer.mean(), geometry: ee.Feature(f).geometry(), scale: 100})
    .evaluate(function (v) {
      if (!v || v.bf === null || v.bf === undefined) { metrics.setValue('No measurements here (outside the study area?)'); return; }
      metrics.setValue('Building cover ' + Math.round(v.bf * 100) + '%, height ' + v.height.toFixed(1) +
        ' m, vegetated ' + Math.round(v.veg * 100) + '%, water ' + Math.round(v.water * 100) +
        '%, bare ' + Math.round(v.bare * 100) + '%' + (p.site ? '  |  site: ' + p.site : ''));
    });
  sv.setUrl('https://www.google.com/maps/@?api=1&map_action=pano&viewpoint=' + ll[0] + ',' + ll[1]);
  earth.setUrl('https://earth.google.com/web/@' + ll[0] + ',' + ll[1] + ',30a,400d,35y,0h,45t,0r');
  var acc = items.filter(function (f) { return decisions[f.properties.cand_id] === true; }).length;
  tally.setValue('Accepted so far in this class: ' + acc + ' (target 50)');
}
function decide(v) { decisions[items[i].properties.cand_id] = v; if (i < items.length - 1) i++; show(); }

var select = ui.Select({items: Object.keys(NAMES).map(function (k) { return {label: NAMES[k], value: k}; }),
  placeholder: 'Choose a class', onChange: function (v) {
    cls = v; i = 0;
    items = all.filter(function (f) { return String(f.properties.lcz) === v; });
    show();
  }});
var saveBtn = ui.Button('Save this class', function () {
  var ids = items.filter(function (f) { return decisions[f.properties.cand_id] === true; })
                 .map(function (f) { return f.properties.cand_id; });
  Export.table.toAsset({collection: cands.filter(ee.Filter.inList('cand_id', ids)),
    description: 'LCZ_ok_' + cls + SUFFIX, assetId: PROJECT + 'LCZ_ok_' + cls + SUFFIX});
  info.setValue('Export task created for ' + NAMES[cls] + ' (' + ids.length + ' squares). Run it in the Tasks tab.');
});
// Hand-drawn squares: layers named lcz2, lcz3, lcz5, lcz6, lcz8, lcz9, lcz10, lczA, lczB, lczD, lczF, lczG
// (Import as FeatureCollection). The button saves each non-empty layer as LCZ_drawn_<code>.
var drawnBtn = ui.Button('Save drawn squares', function () {
  // Layer names: letter classes may have a number after them (lczG, lczG2, lczG3 ...);
  // numbered classes may have an underscore suffix (lcz3, lcz3_b ...). All layers of one class are combined.
  var CODES = {'2': 2, '3': 3, '5': 5, '6': 6, '8': 8, '9': 9, '10': 10, A: 11, B: 12, D: 14, F: 16, G: 17};
  var byCode = {}, counts = {};
  Map.drawingTools().layers().forEach(function (layer) {
    var name = layer.getName();
    var m = name.match(/^lcz([ABDFG])\d*$/) || name.match(/^lcz(10|2|3|5|6|8|9)(_\w*)?$/);
    if (!m) return;
    var code = CODES[m[1]];
    layer.geometries().forEach(function (g) {
      (byCode[code] = byCode[code] || []).push(ee.Feature(g, {lcz: code}));
    });
    counts[code] = (counts[code] || 0) + layer.geometries().length();
  });
  var made = [];
  Object.keys(byCode).forEach(function (code) {
    Export.table.toAsset({collection: ee.FeatureCollection(byCode[code]), description: 'LCZ_drawn_' + code,
                          assetId: PROJECT + 'LCZ_drawn_' + code});
    made.push('class ' + code + ': ' + counts[code] + ' squares');
  });
  info.setValue(made.length ? 'Export tasks created (' + made.join('; ') + '). Run them in the Tasks tab.'
                            : 'No drawn layers found. Name layers like lczG, lczG2 or lcz3, lcz3_b and draw squares first.');
});

ui.root.insert(0, ui.Panel([
  ui.Label('LCZ candidate review', {fontWeight: 'bold', fontSize: '16px'}),
  select, info, metrics, sv, earth, showS2,
  ui.Panel([ui.Button('< Back', function () { if (i > 0) { i--; show(); } }),
            ui.Button('Accept', function () { decide(true); }),
            ui.Button('Reject', function () { decide(false); }),
            ui.Button('Skip >', function () { if (i < items.length - 1) { i++; show(); } })],
           ui.Panel.Layout.flow('horizontal')),
  tally, saveBtn, drawnBtn,
  ui.Label('Check: uniform for 200 m around? Storeys and roof cover match the class? Same in the Sentinel-2 layer?',
           {fontSize: '11px'})
], null, {width: '330px'}));

// ---- When every class is saved, set MERGE to true and run: builds LCZ_training from every saved
// review asset (LCZ_ok_<code> and LCZ_ok_<code>_extra), every hand-placed asset (LCZ_manual_<code>)
// and any hand-drawn layers listed in DRAWN (for example lcz2, lcz8, lcz10). ----
var MERGE = false;
if (MERGE) {
  var DRAWN = [];                       // e.g. [lcz2, lcz8, lcz10] if these layers were drawn in this script
  var parts = [];
  ['2', '3', '5', '6', '8', '9', '10', '11', '12', '14', '16', '17'].forEach(function (k) {
    ['LCZ_ok_' + k, 'LCZ_ok_' + k + '_extra', 'LCZ_ok_' + k + '_csv', 'LCZ_manual_' + k].forEach(function (name) {
      try { ee.data.getAsset(PROJECT + name); parts.push(ee.FeatureCollection(PROJECT + name).select(['lcz'])); } catch (e) {}
    });
  });
  DRAWN.forEach(function (fc) { parts.push(fc.select(['lcz'])); });
  var training = ee.FeatureCollection(parts).flatten();
  print('Training areas per class', training.aggregate_histogram('lcz'));
  Export.table.toAsset({collection: training, description: 'LCZ_training', assetId: PROJECT + 'LCZ_training'});
}

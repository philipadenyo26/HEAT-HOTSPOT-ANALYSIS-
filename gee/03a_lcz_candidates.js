// =====================================================================
// Step 4a (automated). Candidate LCZ training areas from measured criteria
// Measures each 100 m cell: building cover and height (Google Open Buildings 2.5D, 2021),
// mean footprint (Open Buildings v3 outlines), vegetation, water and bare ground
// (Sentinel-2, 2020-2022 dry season). Applies the Stewart and Oke (2012) ranges,
// keeps cells whose 3 x 3 block is uniform, and samples candidates per class.
//
// HOW TO RUN (set MODE below):
//   'TEST'   quick check on a 10 x 10 km box over central Accra (no exports, cheap)
//   'STAGE1' computes the 100 m measurements for the whole region and saves them
//            as the image asset LCZ_metrics_100m (run the export task, wait for it)
//   'STAGE2' reads LCZ_metrics_100m, applies the rules and saves LCZ_candidates
//            (fast and cheap; can be re-run, e.g. after changing UNIFORM)
// =====================================================================
var MODE = 'STAGE2';
var UNIFORM = 7;                 // cells out of 9 in the 3 x 3 block that must share the class
var UNIFORM_RARE = 5;            // lower threshold for classes that are rare in Greater Accra
var RARE = [2, 8, 12];           // LCZ 2, LCZ 8 and LCZ B
var N_PER_CLASS = 80;            // candidates per class (about 50 should survive review)
var EXTRA = false;               // true: a second batch of candidates that skips every cell in LCZ_candidates
                                 //       (saved as LCZ_candidates_extra; set N_PER_CLASS to e.g. 60)

var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var testBox   = ee.Geometry.Rectangle([-0.25, 5.53, -0.15, 5.63]);
var region    = MODE === 'TEST' ? testBox : districts.geometry();
var UTM  = ee.Projection('EPSG:32630');
var P100 = UTM.atScale(100);
var dry  = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));

// ---------------------------------------------------------------------
// Measurements (used by TEST and STAGE1)
// ---------------------------------------------------------------------
function measurements(region) {
  function to100(img, reducer) { return img.reduceResolution({reducer: reducer, maxPixels: 1024}).reproject(P100); }
  // Open Buildings 2.5D is stored at 0.5 m: each 100 m cell reads it on a 5 m sampling grid
  function ob100(img, reducer) {
    return img.setDefaultProjection(UTM.atScale(5)).reduceResolution({reducer: reducer, maxPixels: 1024}).reproject(P100);
  }
  var obCol = ee.ImageCollection('GOOGLE/Research/open-buildings-temporal/v1')
    .filterBounds(region).filter(ee.Filter.calendarRange(2021, 2021, 'year'));
  var ob = obCol.mosaic().setDefaultProjection(obCol.first().select('building_presence').projection());
  var isBldg = ob.select('building_presence').gt(0.5);
  var bf = ob100(isBldg.unmask(0).toFloat(), ee.Reducer.mean()).rename('bf');                       // building cover (0-1)
  var h  = ob100(ob.select('building_height').updateMask(isBldg), ee.Reducer.mean()).unmask(0).rename('height');

  // Mean footprint from Open Buildings v3 outlines (area of each building, assigned to the cell of its centre)
  var outlines = ee.FeatureCollection('GOOGLE/Research/open-buildings/v3/polygons')
    .filterBounds(region).filter(ee.Filter.gte('confidence', 0.75))
    .map(function (f) { return ee.Feature(f.geometry().centroid(1), {area: f.get('area_in_meters')}); });
  var foot = outlines.reduceToImage(['area'], ee.Reducer.mean()).reproject(P100).unmask(0).rename('footprint');

  // Sentinel-2 dry-season composite
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(region)
    .filterDate('2020-01-01', '2023-01-01').filter(dry).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(function (img) {
      var scl = img.select('SCL');
      var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
      return img.updateMask(ok).divide(10000).select(['B3', 'B4', 'B8']);
    }).median().setDefaultProjection(UTM.atScale(10));
  var ndvi = s2.normalizedDifference(['B8', 'B4']);
  var ndwi = s2.normalizedDifference(['B3', 'B8']);
  return ee.Image.cat([
    bf, h, foot,
    to100(ndvi.gt(0.3).toFloat(), ee.Reducer.mean()).rename('veg'),          // vegetated share
    to100(ndvi.gt(0.55).toFloat(), ee.Reducer.mean()).rename('dense'),       // dense canopy share
    to100(ndwi.gt(0.1).toFloat(), ee.Reducer.mean()).rename('water'),
    to100(ndvi.lt(0.15).and(ndwi.lte(0)).toFloat(), ee.Reducer.mean()).rename('bare'),
    to100(ndvi, ee.Reducer.mean()).rename('ndvi_mean'),
    to100(ndvi, ee.Reducer.stdDev()).rename('ndvi_sd')
  ]).toFloat();
}

// ---------------------------------------------------------------------
// Class rules and uniformity (used by TEST and STAGE2)
// ---------------------------------------------------------------------
function classify(m) {
  var bf = m.select('bf'), h = m.select('height'), foot = m.select('footprint'), veg = m.select('veg');
  var dense = m.select('dense'), water = m.select('water'), bare = m.select('bare');
  var nMean = m.select('ndvi_mean'), nSd = m.select('ndvi_sd');
  var none = bf.lt(0.02);
  var c = ee.Image(0).reproject(P100)
    .where(water.gte(0.95), 17)                                                                   // G water
    .where(none.and(bare.gte(0.9)), 16)                                                           // F bare soil or sand
    .where(none.and(nMean.gte(0.2)).and(nMean.lte(0.45)).and(nSd.lt(0.06)).and(dense.lt(0.1)).and(water.lt(0.05)), 14) // D low plants
    .where(none.and(veg.gte(0.4)).and(dense.gte(0.05)).and(dense.lt(0.5)).and(nSd.gte(0.07)), 12) // B scattered trees (dry-season canopy)
    .where(none.and(dense.gte(0.9)), 11)                                                          // A dense trees
    .where(bf.gte(0.05).and(bf.lt(0.2)).and(veg.gte(0.6)).and(h.lte(10)), 9)                     // 9 sparsely built
    .where(bf.gte(0.2).and(bf.lt(0.4)).and(h.gte(3)).and(h.lte(10)).and(foot.gt(0)).and(foot.lt(500)).and(veg.gte(0.1)), 6) // 6 open low-rise
    .where(bf.gte(0.2).and(bf.lte(0.5)).and(h.lte(10)).and(veg.lt(0.2)).and(foot.gte(500)), 8)  // 8 large low-rise
    .where(bf.gte(0.4).and(h.gte(3)).and(h.lte(10)).and(veg.lt(0.3)).and(foot.gt(0)).and(foot.lt(300)), 3) // 3 compact low-rise
    .where(bf.gte(0.4).and(h.gt(10)).and(h.lte(25)).and(veg.lt(0.2)), 2);                       // 2 compact mid-rise
  var lcz = ee.Image(0).reproject(P100);
  [2, 3, 6, 8, 9, 11, 12, 14, 16, 17].forEach(function (k) {
    var same = c.eq(k).reduceNeighborhood({reducer: ee.Reducer.sum(),
      kernel: ee.Kernel.square({radius: 1, units: 'pixels', normalize: false})}).reproject(P100);
    var need = RARE.indexOf(k) >= 0 ? UNIFORM_RARE : UNIFORM;
    lcz = lcz.where(c.eq(k).and(same.gte(need)), k);
  });
  return {raw: c, lcz: lcz.selfMask().rename('lcz')};
}

var pal = ['d10000', 'ff0000', 'ff9955', 'bcbcbc', 'ffccaa', '006a00', '00aa00', 'b9db79', 'fbf7ae', '6a6aff'];
var hist = ee.Reducer.frequencyHistogram();
Map.setOptions('HYBRID');
Map.addLayer(districts.style({color: 'yellow', fillColor: '00000000'}), {}, 'Districts');

if (MODE === 'TEST') {
  var mT = measurements(testBox);
  var rT = classify(mT);
  Map.centerObject(testBox, 13);
  print('Average measurements', mT.select(['bf', 'height', 'veg', 'dense', 'water', 'bare', 'ndvi_mean'])
    .addBands(mT.select('footprint').updateMask(mT.select('footprint').gt(0)))
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: testBox, scale: 100, maxPixels: 1e9}));
  print('Cells per class before uniformity', rT.raw.reduceRegion({reducer: hist, geometry: testBox, scale: 100, maxPixels: 1e9}));
  print('Cells per class after uniformity', rT.lcz.reduceRegion({reducer: hist, geometry: testBox, scale: 100, maxPixels: 1e9}));
  Map.addLayer(rT.lcz.clip(testBox), {min: 2, max: 17, palette: pal}, 'Candidate cells (test box)');
}

if (MODE === 'STAGE1') {
  Map.centerObject(districts, 10);
  Export.image.toAsset({image: measurements(districts.geometry().buffer(500)).clip(districts.geometry()),
    description: 'LCZ_metrics_100m', assetId: PROJECT + 'LCZ_metrics_100m',
    region: districts.geometry().bounds(), crs: 'EPSG:32630', scale: 100, maxPixels: 1e10});
  print('Stage 1: run the LCZ_metrics_100m task in the Tasks tab. When it has finished, set MODE to STAGE2.');
}

if (MODE === 'STAGE2') {
  var m = ee.Image(PROJECT + 'LCZ_metrics_100m');
  var r = classify(m);
  var lcz = r.lcz.clip(districts.geometry());
  print('Uniform cells per class (whole region)', lcz.reduceRegion({reducer: hist, geometry: districts.geometry(),
    scale: 100, maxPixels: 1e10}));
  var OUT = EXTRA ? 'LCZ_candidates_extra' : 'LCZ_candidates';
  var pool = lcz;
  if (EXTRA) {                                    // leave out cells already proposed in the first batch
    var used = ee.Image(0).paint(ee.FeatureCollection(PROJECT + 'LCZ_candidates'), 1).reproject(P100);
    pool = lcz.updateMask(used.not());
  }
  var cands = pool.addBands(m).stratifiedSample({numPoints: N_PER_CLASS, classBand: 'lcz',
    region: districts.geometry(), projection: P100, seed: EXTRA ? 7 : 42, geometries: true, tileScale: 4})
    .map(function (f) {
      return f.setGeometry(f.geometry().buffer(50).bounds(1, UTM))
              .set('cand_id', ee.String(EXTRA ? 'x' : '').cat(ee.String(f.get('system:index'))));
    });
  print('Candidates per class', cands.aggregate_histogram('lcz'));
  Map.centerObject(districts, 10);
  Map.addLayer(lcz, {min: 2, max: 17, palette: pal}, 'Uniform candidate cells', false);
  Map.addLayer(cands.style({color: 'ff0000', fillColor: '00000000', width: 1}), {}, 'Candidates (100 m squares)');
  Export.table.toAsset({collection: cands, description: OUT, assetId: PROJECT + OUT});
  Export.table.toDrive({collection: cands, description: OUT + '_metrics', folder: 'GEE_thesis',
    fileFormat: 'CSV', selectors: ['cand_id', 'lcz', 'bf', 'height', 'footprint', 'veg', 'dense', 'water', 'bare', 'ndvi_mean', 'ndvi_sd']});
}

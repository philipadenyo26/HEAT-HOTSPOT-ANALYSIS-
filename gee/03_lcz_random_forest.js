// =====================================================================
// LCZ classification with Random Forest - FINAL RUN (Chapter 3, Section 3.8)
// Stages run automatically: first Run saves the sampled pixels; the next Run creates the accuracy and map tasks.
//
// Training data (100 m squares):
//   (a) LCZ_samples_clean - CSV site points (split train/valid BY SITE); used only for classes
//       without any reviewed asset
//   (b) every reviewed asset: LCZ_ok_<code>[_extra|_csv], LCZ_manual_<code>, LCZ_drawn_<code><letters>
//       (squares without a split are split 75/25 at random)
// Training: up to 150 squares per class, spatially balanced - up to 38 from each quadrant
//   (NE, NW, SE, SW about the centre of the 29 districts), then filled to 150 from the rest.
// Features: eight-variable Sentinel-2 stack (2020-2022 dry season), from the saved asset, plus building
//   cover, height and footprint from Open Buildings (LCZ_metrics_100m): eleven inputs.
// Random Forest: 100 trees, minimum leaf population 5, bag fraction 0.5, seed 42.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var region    = districts.geometry().buffer(2000);
var UTM    = ee.Projection('EPSG:32630');
var FOLDER = 'GEE_thesis';
var CLASSES = [2, 3, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17];
var MAX_TRAIN_PER_CLASS = 150;
var PER_QUADRANT = Math.ceil(MAX_TRAIN_PER_CLASS / 4);   // 38
var SAMPLE_SCALE = 20;            // pixels sampled every 20 m: about 25 per square
var USE_SITE_POINTS = true;       // CSV points are used only for classes that have no reviewed asset
var MORPHOLOGY = true;            // true: add building cover, height, footprint (Open Buildings) to the 8 Sentinel-2 inputs
var TAG = MORPHOLOGY ? 'morph' : 'base';
var SUMMARY = false;              // one Console request with all accuracy figures (set false if it times out;
                                  // the same figures are in the lcz_accuracy export)

// ---- (a) CSV site points ----
function num(v) { return ee.Number(ee.Algorithms.If(ee.Algorithms.ObjectType(v).equals('String'), ee.Number.parse(v), v)); }
var HAS_CSV = true;
try { ee.data.getAsset(PROJECT + 'LCZ_samples_clean'); } catch (e) { HAS_CSV = false; }
if (!HAS_CSV) { USE_SITE_POINTS = false; print('LCZ_samples_clean not found: using reviewed squares only (expected once every class is reviewed).'); }
var csv = !HAS_CSV ? ee.FeatureCollection([]) : ee.FeatureCollection(PROJECT + 'LCZ_samples_clean').map(function (f) {
  var pt = ee.Geometry(ee.Algorithms.If(f.get('longitude'),
    ee.Geometry.Point([num(f.get('longitude')), num(f.get('latitude'))]), f.geometry()));
  return ee.Feature(pt, {lcz: num(f.get('lcz_code')).int(), split: f.get('split'), source: 'site_points'});
});

// ---- (b) Reviewed assets (any ending) ----
var extra = [], reviewedCodes = [], reviewedAssets = [];
ee.data.listAssets(PROJECT.slice(0, -1)).assets.forEach(function (a) {
  var name = a.name.split('/').pop();
  var m = name.match(/^LCZ_(ok|manual|drawn)_(\d+)[a-z]*(_extra|_csv)?$/);
  if (!m || a.type !== 'TABLE') return;
  var k = parseInt(m[2], 10);
  if (CLASSES.indexOf(k) < 0) return;
  var origin = /^LCZ_ok_\d+(_extra)?$/.test(name) ? 'gee_candidate' : 'user';
  extra.push(ee.FeatureCollection(a.name).select(['lcz', 'split']).map(function (f) { return f.set('source', origin); }));
  reviewedAssets.push(name);
  if (reviewedCodes.indexOf(k) < 0) reviewedCodes.push(k);
});
print('Reviewed assets found (client-side, no compute)', reviewedAssets);
var missing = CLASSES.filter(function (k) { return reviewedCodes.indexOf(k) < 0; });
if (missing.length) print('WARNING - classes with NO reviewed squares' + (USE_SITE_POINTS ? ' (CSV points used for these)' : ' (they will be MISSING from the map)'), missing);
else print('All 12 classes have reviewed squares.');

csv = csv.filter(ee.Filter.inList('lcz', reviewedCodes).not());
var reviewed = ee.FeatureCollection(extra).flatten().randomColumn('r', 42).map(function (f) {
  var rnd = ee.Algorithms.If(ee.Number(f.get('r')).lt(0.75), 'train', 'valid');
  return f.set('split', ee.Algorithms.If(f.get('split'), f.get('split'), rnd));
});

// ---- Standard 100 m squares, whole-number class codes, quadrant about the regional centre ----
var centre = districts.geometry().centroid(1).transform('EPSG:4326', 1).coordinates();
var cLon = ee.Number(centre.get(0)), cLat = ee.Number(centre.get(1));
var samples = (USE_SITE_POINTS ? csv.merge(reviewed) : reviewed).map(function (f) {
  var c = f.geometry().centroid(1).transform('EPSG:4326', 1).coordinates();
  var lon = ee.Number(c.get(0)), lat = ee.Number(c.get(1));
  var quad = ee.String(ee.Algorithms.If(lat.gte(cLat), 'N', 'S')).cat(ee.Algorithms.If(lon.gte(cLon), 'E', 'W'));
  return ee.Feature(f.geometry().centroid(1).buffer(50).bounds(1, UTM), {
    lcz: ee.Number(f.get('lcz')).int(), split: f.get('split'), source: f.get('source'),
    quadrant: quad, lon: lon, lat: lat, uid: f.get('system:index')});
}).filter(ee.Filter.inList('lcz', CLASSES));

// ---- Spatially balanced training selection: up to 38 per quadrant, then fill to 150 ----
var trainAll = samples.filter(ee.Filter.eq('split', 'train'));
var trainSq = ee.FeatureCollection(CLASSES.map(function (k) {
  var cls = trainAll.filter(ee.Filter.eq('lcz', k)).randomColumn('p', k);
  var first = ee.FeatureCollection(['NE', 'NW', 'SE', 'SW'].map(function (q) {
    return cls.filter(ee.Filter.eq('quadrant', q)).sort('p').limit(PER_QUADRANT);
  })).flatten();
  var rest = cls.filter(ee.Filter.inList('uid', first.aggregate_array('uid')).not()).sort('p');
  var room = ee.Number(MAX_TRAIN_PER_CLASS).subtract(first.size()).max(0);
  return first.merge(rest.limit(room));
})).flatten();
var validSq = samples.filter(ee.Filter.eq('split', 'valid'));

// ---- Sentinel-2 stack (saved asset; built and exported only if missing) ----
var STACK_ASSET = PROJECT + 'S2_stack_2020_2022';
var stack;
try {
  ee.data.getAsset(STACK_ASSET);
  stack = ee.Image(STACK_ASSET);
  print('Using saved Sentinel-2 stack');
} catch (e) {
  var dry = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(region).filterDate('2020-01-01', '2023-01-01').filter(dry)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(function (img) {
      var scl = img.select('SCL');
      var ok = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
      return img.updateMask(ok).divide(10000).select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).resample('bilinear');
    }).median();
  stack = s2.addBands([s2.normalizedDifference(['B8', 'B4']).rename('NDVI'),
                       s2.normalizedDifference(['B11', 'B8']).rename('NDBI')]).clip(region).toFloat();
  Export.image.toAsset({image: stack, description: 'S2_stack_2020_2022', assetId: STACK_ASSET,
    region: region, scale: 10, crs: 'EPSG:32630', maxPixels: 1e11});
  print('Saved stack NOT found: the composite is computed on the fly (heavy). Run the S2_stack_2020_2022 task if possible.');
}
var BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12', 'NDVI', 'NDBI'];
if (MORPHOLOGY) {
  stack = stack.addBands(ee.Image(PROJECT + 'LCZ_metrics_100m').select(['bf', 'height', 'footprint']).unmask(0).toFloat());
  BANDS = BANDS.concat(['bf', 'height', 'footprint']);
}

// ---- Stage PIXELS: sample once and save as assets (tiled batch exports, no time-out) ----
var PIX_TRAIN = PROJECT + 'LCZ_pix_train_' + TAG, PIX_VALID = PROJECT + 'LCZ_pix_valid_' + TAG;
// The stage is chosen automatically: no editing between runs.
//   saved pixels missing  -> PIXELS: creates the two pixel tasks
//   saved pixels present  -> RESULTS: creates the accuracy task and the map task
function exists(id) { try { ee.data.getAsset(id); return true; } catch (e) { return false; } }
var STAGE = (exists(PIX_TRAIN) && exists(PIX_VALID)) ? 'RESULTS' : 'PIXELS';
print('STAGE (automatic):', STAGE);
if (STAGE === 'PIXELS') {
  var trainPixS = stack.select(BANDS).sampleRegions({collection: trainSq, properties: ['lcz'],
    scale: SAMPLE_SCALE, projection: UTM.atScale(SAMPLE_SCALE), tileScale: 8, geometries: true});   // assets need a location
  var validPixS = stack.select(BANDS).sampleRegions({collection: validSq, properties: ['lcz', 'source'],
    scale: SAMPLE_SCALE, projection: UTM.atScale(SAMPLE_SCALE), tileScale: 8, geometries: true});
  Export.table.toAsset({collection: trainPixS, description: 'LCZ_pix_train_' + TAG, assetId: PIX_TRAIN});
  Export.table.toAsset({collection: validPixS, description: 'LCZ_pix_valid_' + TAG, assetId: PIX_VALID});
  Export.table.toDrive({collection: trainSq.map(function (f) { return f.set('set', 'train'); })
      .merge(validSq.map(function (f) { return f.set('set', 'valid'); })),
    description: 'lcz_training_validation_squares_' + TAG, folder: FOLDER, fileFormat: 'CSV',
    selectors: ['lcz', 'set', 'quadrant', 'source', 'lat', 'lon']});
  print('NEXT: run LCZ_pix_train_' + TAG + ' and LCZ_pix_valid_' + TAG +
        ' in Tasks. When both are Completed, just click Run again.');
}

// ---- Stages ACCURACY and MAP: train on the saved pixels ----
if (STAGE === 'RESULTS') {
  var trainPix = ee.FeatureCollection(PIX_TRAIN);
  var validPix = ee.FeatureCollection(PIX_VALID);
  var rf = ee.Classifier.smileRandomForest({numberOfTrees: 100, minLeafPopulation: 5, bagFraction: 0.5, seed: 42})
    .train({features: trainPix, classProperty: 'lcz', inputProperties: BANDS});
  print('Classifier inputs', BANDS);
}

if (STAGE === 'RESULTS') {
  var validClassified = validPix.classify(rf);
  var cm = validClassified.errorMatrix('lcz', 'classification', CLASSES);
  var cmUser = validClassified.filter(ee.Filter.neq('source', 'gee_candidate')).errorMatrix('lcz', 'classification', CLASSES);
  var importance = ee.Dictionary(rf.explain().get('importance'));
  var LABELS = {2: 'LCZ 2', 3: 'LCZ 3', 5: 'LCZ 5', 6: 'LCZ 6', 8: 'LCZ 8', 9: 'LCZ 9', 10: 'LCZ 10',
                11: 'LCZ A', 12: 'LCZ B', 14: 'LCZ D', 16: 'LCZ F', 17: 'LCZ G'};
  var names = CLASSES.map(function (k) { return LABELS[k]; });
  var producers = ee.Dictionary.fromLists(names, cm.producersAccuracy().project([0]).toList());
  var users = ee.Dictionary.fromLists(names, cm.consumersAccuracy().project([1]).toList());
  if (SUMMARY) print('SUMMARY', ee.Dictionary({overall_accuracy: cm.accuracy(), kappa: cm.kappa(),
    OA_user_placed_only: cmUser.accuracy(), producer_accuracy: producers, user_accuracy: users}));
  Export.table.toDrive({collection: ee.FeatureCollection([ee.Feature(null, {
      inputs: TAG, class_order: ee.List(CLASSES).join(','), confusion_matrix: ee.String.encodeJSON(cm.array().toList()),
      OA: cm.accuracy(), kappa: cm.kappa(), OA_user_placed: cmUser.accuracy(), kappa_user_placed: cmUser.kappa(),
      producer_accuracy: ee.String.encodeJSON(producers), user_accuracy: ee.String.encodeJSON(users),
      feature_importance: ee.String.encodeJSON(importance)})]),
    description: 'lcz_accuracy_' + TAG, folder: FOLDER, fileFormat: 'CSV'});
  // Accuracy of the 100 m LCZ map itself: each validation square takes the majority class of its
  // 10 m predictions (the same rule that builds LCZ_100m), compared with its reference label.
  var sqPred = stack.select(BANDS).classify(rf).rename('pred')
    .reduceRegions({collection: validSq, reducer: ee.Reducer.mode(), scale: 10, crs: 'EPSG:32630', tileScale: 4})
    .filter(ee.Filter.notNull(['mode']))
    .map(function (f) { return f.set('mapped', ee.Number(f.get('mode')).int()); });
  var cmSq = sqPred.errorMatrix('lcz', 'mapped', CLASSES);
  var cmSqUser = sqPred.filter(ee.Filter.neq('source', 'gee_candidate')).errorMatrix('lcz', 'mapped', CLASSES);
  Export.table.toDrive({collection: ee.FeatureCollection([ee.Feature(null, {
      inputs: TAG, level: '100 m map (majority per validation square)', squares: sqPred.size(),
      class_order: ee.List(CLASSES).join(','), confusion_matrix: ee.String.encodeJSON(cmSq.array().toList()),
      OA: cmSq.accuracy(), kappa: cmSq.kappa(), OA_user_placed: cmSqUser.accuracy(), kappa_user_placed: cmSqUser.kappa(),
      producer_accuracy: ee.String.encodeJSON(ee.Dictionary.fromLists(names, cmSq.producersAccuracy().project([0]).toList())),
      user_accuracy: ee.String.encodeJSON(ee.Dictionary.fromLists(names, cmSq.consumersAccuracy().project([1]).toList()))})]),
    description: 'lcz_accuracy_100m_' + TAG, folder: FOLDER, fileFormat: 'CSV'});
  print('NEXT: run lcz_accuracy_100m_' + TAG + ' (accuracy of the 100 m map).');
}

if (STAGE === 'RESULTS') {
  var lcz10 = stack.select(BANDS).classify(rf).rename('LCZ');
  var lcz100 = lcz10.reproject({crs: UTM, scale: 10})
    .reduceResolution({reducer: ee.Reducer.mode(), maxPixels: 1024}).reproject({crs: UTM, scale: 100});
  Export.image.toDrive({image: lcz100.toByte(), description: 'LCZ_100m_' + TAG, folder: FOLDER,
    fileNamePrefix: 'LCZ_100m', region: region, scale: 100, crs: 'EPSG:32630', maxPixels: 1e10});
  print('THEN: run LCZ_100m_' + TAG + ' (downloads as LCZ_100m.tif).');
}

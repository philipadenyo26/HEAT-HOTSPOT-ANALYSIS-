// =====================================================================
// Step 4e. Audit the training samples against measured values
// Each 100 m square is checked against LCZ_metrics_100m (building cover, water,
// vegetation, dense canopy, bare ground). Squares whose measurements cannot match
// their label (e.g. a "water" square with no water) are flagged as SUSPECT.
// Outputs: counts per class (Console), red squares on the map, and a CSV
// (lcz_sample_audit) listing every sample with its measurements and verdict.
// =====================================================================
var PROJECT = 'projects/lunar-tine-512302-q0/assets/';
var UTM = ee.Projection('EPSG:32630');
var m = ee.Image(PROJECT + 'LCZ_metrics_100m');

function num(v) { return ee.Number(ee.Algorithms.If(ee.Algorithms.ObjectType(v).equals('String'), ee.Number.parse(v), v)); }
var csv = ee.FeatureCollection(PROJECT + 'LCZ_samples_clean').map(function (f) {
  var pt = ee.Geometry(ee.Algorithms.If(f.get('longitude'),
    ee.Geometry.Point([num(f.get('longitude')), num(f.get('latitude'))]), f.geometry()));
  return ee.Feature(pt.buffer(50).bounds(1, UTM), {lcz: num(f.get('lcz_code')).int(),
    sample_id: f.get('sample_id'), site: f.get('site'), source: 'site_points'});
});

// Mean measurements inside each square
var measured = m.select(['bf', 'height', 'veg', 'dense', 'water', 'bare', 'ndvi_mean'])
  .reduceRegions({collection: csv, reducer: ee.Reducer.mean(), scale: 100, crs: 'EPSG:32630'});

// Minimum conditions a square must meet to be plausible for its label (deliberately lenient)
function verdict(f) {
  var k = ee.Number(f.get('lcz'));
  var bf = ee.Number(f.get('bf')), water = ee.Number(f.get('water')), dense = ee.Number(f.get('dense'));
  var veg = ee.Number(f.get('veg')), bare = ee.Number(f.get('bare'));
  var built = bf.gte(0.10);
  var ok = ee.Algorithms.If(k.eq(17), water.gte(0.5),
           ee.Algorithms.If(k.eq(11), dense.gte(0.4).and(bf.lt(0.05)),
           ee.Algorithms.If(k.eq(12), veg.gte(0.2).and(bf.lt(0.05)).and(water.lt(0.2)),
           ee.Algorithms.If(k.eq(14), bf.lt(0.05).and(dense.lt(0.4)).and(water.lt(0.2)),
           ee.Algorithms.If(k.eq(16), bare.gte(0.4).and(bf.lt(0.05)).and(water.lt(0.2)),
           ee.Algorithms.If(k.eq(9), bf.gte(0.02).and(bf.lt(0.3)),
           built.and(water.lt(0.3))))))));                       // classes 2, 3, 5, 6, 8, 10: must be built up
  return f.set('ok', ee.Number(ok));
}
var audited = measured.filter(ee.Filter.notNull(['bf'])).map(verdict);

print('Samples per class', audited.aggregate_histogram('lcz'));
print('SUSPECT samples per class', audited.filter(ee.Filter.eq('ok', 0)).aggregate_histogram('lcz'));
print('Plausible samples per class', audited.filter(ee.Filter.eq('ok', 1)).aggregate_histogram('lcz'));

// Share of plausible samples at each site (1 = every point fits the 2020-2022 dry-season data, 0 = none do)
var bySite = audited.reduceColumns({reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'site'}),
                                    selectors: ['ok', 'site']});
print('Plausible share by site (look for sites near 0)', bySite.get('groups'));

// Dry-season Sentinel-2 composite the classifier learns from, for visual checks of red squares
var dry = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(audited.geometry())
  .filterDate('2020-01-01', '2023-01-01').filter(dry).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)).median();

Map.setOptions('HYBRID');
Map.centerObject(audited, 10);
Map.addLayer(s2, {bands: ['B4', 'B3', 'B2'], min: 0, max: 2500}, 'Sentinel-2 dry season 2020-2022', false);
Map.addLayer(audited.filter(ee.Filter.eq('ok', 1)).style({color: '00ff00', fillColor: '00000000', width: 2}), {}, 'Plausible');
Map.addLayer(audited.filter(ee.Filter.eq('ok', 0)).style({color: 'ff0000', fillColor: '00000000', width: 2}), {}, 'SUSPECT');

Export.table.toDrive({collection: audited, description: 'lcz_sample_audit', folder: 'GEE_thesis', fileFormat: 'CSV',
  selectors: ['sample_id', 'site', 'lcz', 'ok', 'bf', 'height', 'veg', 'dense', 'water', 'bare', 'ndvi_mean']});

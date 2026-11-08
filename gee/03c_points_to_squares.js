// =====================================================================
// Step 4c. Turn locations identified by the researcher (e.g. 100 factories) into
// 100 m training squares for one LCZ class.
// 1. Make a CSV with three columns: name, latitude, longitude
// 2. Upload it to GEE (Assets > New > CSV file), e.g. as 'LCZ10_points'
// 3. Set CLASS and TABLE below, Run, check the squares on the map, run the export task
// Points closer than 150 m to an earlier point are dropped so squares do not overlap.
// =====================================================================
var PROJECT = 'projects/lunar-tine-512302-q0/assets/';
var CLASS = 10;                        // LCZ code: 10 heavy industry, 8 large low-rise, 2 compact mid-rise
var TABLE = 'LCZ10_points';            // name of the uploaded CSV asset
var MIN_GAP = 150;                     // metres between points

var UTM = ee.Projection('EPSG:32630');
var pts = ee.FeatureCollection(PROJECT + TABLE);
var list = pts.toList(5000).getInfo();

// Keep points that are at least MIN_GAP metres from every point already kept
function dist(a, b) {
  var R = 6371000, toRad = Math.PI / 180;
  var dLat = (b[1] - a[1]) * toRad, dLon = (b[0] - a[0]) * toRad;
  var x = Math.sin(dLat / 2) * Math.sin(dLat / 2) +
          Math.cos(a[1] * toRad) * Math.cos(b[1] * toRad) * Math.sin(dLon / 2) * Math.sin(dLon / 2);
  return 2 * R * Math.asin(Math.sqrt(x));
}
var kept = [];
list.forEach(function (f) {
  var c = f.geometry.coordinates;
  var ok = kept.every(function (k) { return dist(k.geometry.coordinates, c) >= MIN_GAP; });
  if (ok) kept.push(f);
});
print('Points uploaded', list.length, 'Points kept (no overlap)', kept.length);

var squares = ee.FeatureCollection(kept.map(function (f) {
  return ee.Feature(ee.Geometry.Point(f.geometry.coordinates).buffer(50).bounds(1, UTM),
                    {lcz: CLASS, name: f.properties.name || ''});
}));

var dry = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'), ee.Filter.calendarRange(1, 3, 'month'));
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(squares.geometry())
  .filterDate('2020-01-01', '2023-01-01').filter(dry).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20)).median();
Map.setOptions('HYBRID');
Map.centerObject(squares, 11);
Map.addLayer(s2, {bands: ['B4', 'B3', 'B2'], min: 0, max: 2500}, 'Sentinel-2 2020-2022', false, 0.7);
Map.addLayer(squares.style({color: 'ff00ff', fillColor: '00000000', width: 2}), {}, 'Training squares');

Export.table.toAsset({collection: squares, description: 'LCZ_manual_' + CLASS, assetId: PROJECT + 'LCZ_manual_' + CLASS});

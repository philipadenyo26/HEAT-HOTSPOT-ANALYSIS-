// =====================================================================
// Step 7. CMIP6 delta-change factors (NASA NEX-GDDP-CMIP6, daily tas)
// Chapter 3, Phase 7 (Objective 4)
// Delta = dry-season mean tas (future window) - dry-season mean tas (1981-2010), per model,
// then the ensemble mean and standard deviation across models.
// =====================================================================
var PROJECT   = 'projects/lunar-tine-512302-q0/assets/';
var districts = ee.FeatureCollection(PROJECT + 'GAR_districts');
var box = districts.geometry().bounds().buffer(60000).bounds();   // extra margin for bilinear interpolation
var FOLDER = 'GEE_thesis';
var MODELS = ['ACCESS-CM2', 'CNRM-ESM2-1', 'GFDL-ESM4', 'INM-CM5-0', 'MPI-ESM1-2-HR', 'MRI-ESM2-0'];
var drySeason = ee.Filter.or(ee.Filter.calendarRange(11, 12, 'month'),
                             ee.Filter.calendarRange(1, 3, 'month'));
var col = ee.ImageCollection('NASA/GDDP-CMIP6').filterBounds(box).select('tas');

function dryMean(model, scenario, start, end) {
  return col.filter(ee.Filter.eq('model', model)).filter(ee.Filter.eq('scenario', scenario))
            .filterDate(start, end).filter(drySeason).mean();
}

// Check that every model has data for each experiment
['historical', 'ssp245', 'ssp585'].forEach(function (s) {
  print(s, col.filter(ee.Filter.eq('scenario', s)).filter(ee.Filter.inList('model', MODELS))
              .aggregate_array('model').distinct());
});

var WINDOWS = {'2041_2060': ['2041-01-01', '2061-01-01'], '2081_2100': ['2081-01-01', '2101-01-01']};
var rows = [];
['ssp245', 'ssp585'].forEach(function (scen) {
  Object.keys(WINDOWS).forEach(function (w) {
    var deltas = ee.ImageCollection(MODELS.map(function (m) {
      var base = dryMean(m, 'historical', '1981-01-01', '2011-01-01');
      var fut  = dryMean(m, scen, WINDOWS[w][0], WINDOWS[w][1]);
      return fut.subtract(base).rename('dT').set('model', m);
    }));
    var mean = deltas.mean().rename('dT_mean');
    var sd   = deltas.reduce(ee.Reducer.stdDev()).rename('dT_sd');
    var name = 'delta_' + scen + '_' + w;
    Export.image.toDrive({image: mean.addBands(sd).toFloat(), description: name, folder: FOLDER,
      region: box, scale: 27830, crs: 'EPSG:4326', maxPixels: 1e9});
    rows.push(deltas.map(function (img) {
      var v = img.reduceRegion({reducer: ee.Reducer.mean(), geometry: districts.geometry(), scale: 27830});
      return ee.Feature(null, {scenario: scen, window: w, model: img.get('model'), dT: v.get('dT')});
    }));
  });
});
Export.table.toDrive({collection: ee.FeatureCollection(rows).flatten(), description: 'delta_by_model',
  folder: FOLDER, fileFormat: 'CSV', selectors: ['scenario', 'window', 'model', 'dT']});

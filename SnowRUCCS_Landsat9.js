/**
 * SnowRUCCS - Snow Recognition Under Cloud Cover System
 * 基于 Landsat 9 的云检测与积雪识别算法
 *
 * 功能说明：
 *   1. 对 Landsat 9 TOA 影像进行预处理（波段重命名、添加辅助数据）
 *   2. Pass One: 潜在云像元识别 (PCP) —— 包含基础测试、白度测试、HOT 测试、波段比值测试、卷云测试
 *   3. 水体掩膜与陆地区域提取
 *   4. Pass Two: 云概率计算 —— 陆地温度概率、陆地变化概率、综合云概率
 *   5. 最终云层与云阴影层生成
 *   6. 积雪识别与云下积雪估计 (AllSnow)
 *   7. 结果导出至 Google Drive
 *
 * 数据源：
 *   - Landsat 9 TOA (LANDSAT/LC09/C02/T1_TOA)
 *   - JRC Global Surface Water (JRC/GSW1_2/GlobalSurfaceWater)
 *   - SRTM DEM (USGS/SRTMGL1_003)
 *
 * 使用说明：
 *   在 Google Earth Engine Code Editor 中运行本脚本。
 *   需要提前定义研究区变量 roi（可使用 ee.Geometry.Polygon 或 ee.FeatureCollection）。
 */

// ============================================================
// 0. 研究区定义与全局参数
// ============================================================

// 研究区 (ROI) - 请根据实际需求修改
// 方式1: 手动定义多边形
// var roi = ee.Geometry.Polygon(
//   [[[99.55872144025032, 38.53071351009249],
//     [99.06982983868782, 36.80901308821066],
//     [101.12976636212532, 36.41216852020269],
//     [101.66260327618782, 38.129968875136],
//     [99.55872144025032, 38.53071351009249]]]);

// 方式2: 使用 FeatureCollection 矢量数据
// var roi = ee.FeatureCollection("users/ycqysl/BBH");

Map.centerObject(roi, 10);
Map.addLayer(roi, {}, 'BBH');

// 全局数据集
var GSWO = ee.Image("JRC/GSW1_2/GlobalSurfaceWater");
var TOAL9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA");
var SRTM = ee.Image("USGS/SRTMGL1_003");
var worldWideGeom = ee.Geometry.Rectangle([-180, -90, 180, 90], 'EPSG:4326', false);

// 可视化参数
var VisParam = {
  opacity: 1,
  bands: ["NIR", "Red", "Green"],
  min: 0.067354,
  max: 0.818746,
  gamma: 1.364
};
var CloudVis = {
  opacity: 1,
  bands: ["Cloud"],
  palette: ["fdfff5", "172dc8"]
};

// 时间范围与轨道参数
var s_date = ee.Date.fromYMD(2022, 02, 28);
var e_date = ee.Date.fromYMD(2022, 03, 01);

// ============================================================
// 1. 影像筛选与预处理
// ============================================================

// 筛选 Landsat 9 影像集
var imgCol8 = TOAL9.filterDate(s_date, e_date)
                    .filterBounds(roi)
                    .filter(ee.Filter.eq('TARGET_WRS_PATH', 133))
                    .filter(ee.Filter.eq('TARGET_WRS_ROW', 34));

/**
 * 波段重命名与辅助数据添加
 * @param {ee.Image} img - Landsat 9 TOA 影像
 * @return {ee.Image} 重命名波段并添加 BT、水体、地形等辅助波段后的影像
 */
function RenameAddBandsL8(img) {
  // 亮温 (BT)：将热红外波段转换为摄氏度
  var bt = img.select('B10').subtract(273.15).toFloat().rename('BT');

  // 辅助数据
  var water = GSWO.select('occurrence');
  var slope = ee.Terrain.slope(SRTM).rename('slope');
  var aspect = ee.Terrain.aspect(SRTM).rename('aspect');
  var elevation = SRTM.select('elevation').rename('elevation');

  // 波段重命名：B2-B7, B9 → Blue, Green, Red, NIR, SWIR1, SWIR2, Cirrus
  var renamed = img.select(
    ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B9'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Cirrus']
  );

  // 合并所有波段
  var addBands = renamed.addBands(bt)
                        .addBands(water)
                        .addBands(slope)
                        .addBands(aspect)
                        .addBands(elevation);

  // 保留时间戳和索引属性
  var time_start = img.get("system:time_start");
  var index = img.get('system:index');
  var out = addBands.set("system:time_start", time_start)
                    .set('system:index', index);
  return out.copyProperties(img);
}

var L8 = imgCol8.map(RenameAddBandsL8).sort("system:time_start");
print('L8', L8);
Map.addLayer(L8, VisParam, 'JCSl8');

// 导出预处理影像
Export.image.toDrive({
  image: L8,
  description: '20220401_L9',
  folder: 'Landsat9',
  region: roi,
  scale: 30,
  maxPixels: 1.0E13,
});

// ============================================================
// 2. 光谱指数计算
// ============================================================

/** 归一化差异雪盖指数 (NDSI) = (Green - SWIR1) / (Green + SWIR1) */
function NDSI(img) {
  return img.normalizedDifference(['Green', 'SWIR1']);
}

/** 归一化植被指数 (NDVI) = (NIR - Red) / (NIR + Red) */
function NDVI(img) {
  return img.normalizedDifference(['NIR', 'Red']);
}

/** 归一化建筑指数 (NDBI) = (SWIR1 - NIR) / (SWIR1 + NIR) */
function NDBI(img) {
  return img.normalizedDifference(['SWIR1', 'NIR']);
}

// ============================================================
// 3. Pass One - 潜在云像元 (PCP) 识别
// ============================================================

/**
 * 3.1 基础测试
 * 条件: SWIR2 > 0.03 且 NDSI < 0.8 且 NDVI < 0.8 且 BT < 27°C
 */
function Basictest(img) {
  var basic = (img.select('SWIR2').gt(0.03))
              .and(NDSI(img).lt(0.8))
              .and(NDVI(img).lt(0.8))
              .and(img.select('BT').lt(27))
              .rename('basic');
  return basic;
}

/**
 * 3.2 白度测试
 * 计算可见光波段与均值的偏差程度，衡量像元"白色"特征
 */
function Whitetest(img) {
  var ADD = img.select('Blue').add(img.select('Green')).add(img.select('Red'));
  var mean = ADD.reduce(ee.Reducer.mean());
  var a = img.expression(
    '(abs(B1 - B5) + abs(B2 - B5) + abs(B3 - B5)) / B5', {
      'B1': img.select('Blue'),
      'B2': img.select('Green'),
      'B3': img.select('Red'),
      'B5': mean
    });
  var white = a.lt(0.7).multiply(a);
  return white.rename('whiteness');
}

/**
 * 3.3 HOT (Haze Optimized Transformation) 测试
 * HOT = Blue - 0.5 * Red - 0.08，云像元通常 HOT > 0
 */
function HOT(img) {
  var hot = img.expression(
    'B1 - 0.5 * B3 - 0.08', {
      'B1': img.select('Blue'),
      'B3': img.select('Red')
    });
  return hot.gt(0).rename('HOT');
}

/**
 * 3.4 NIR/SWIR1 波段比值测试
 * 条件: NIR / SWIR1 > 0.75
 */
function B4B5(img) {
  var ratio = img.select('NIR').divide(img.select('SWIR1'));
  return ratio.gt(0.75).rename('Ratio4_5');
}

/**
 * 3.5 卷云概率
 * 利用 Cirrus 波段计算卷云概率
 */
function Cirrus_p(img) {
  var p = img.select('Cirrus').divide(4);
  return p.where(p.lt(0), 0);
}

/**
 * 3.6 PCP（潜在云像元）综合判别
 * 综合以上所有测试结果
 */
function PCPL8(img) {
  var pcp = (Basictest(img).eq(1))
            .and(Whitetest(img).lt(0.8))
            .and(HOT(img).eq(1))
            .and(B4B5(img).eq(1))
            .and(Cirrus_p(img).gt(0.01));
  return pcp.rename('PCP');
}

var pcpcover = L8.map(PCPL8).first();
Map.addLayer(pcpcover, {}, 'PCPL8');

// ============================================================
// 4. 水体与陆地掩膜
// ============================================================

/**
 * 水体识别
 * 结合 NDVI、NIR 阈值与 JRC 全球地表水数据识别水体
 */
function Water(img) {
  var condi1 = (NDVI(img).lt(0.01)).and(img.select('NIR').lt(0.11));
  var condi2 = (NDVI(img).lt(0.1)).and(NDVI(img).gt(0)).and(img.select('NIR').lt(0.05));
  var water0 = (condi1.eq(1)).or(condi2.eq(1));
  var wat = water0.eq(1);
  var occur_gs = img.select('occurrence').unmask().neq(0);
  var water = wat.and(occur_gs.eq(1));
  return water.rename('water');
}

/**
 * 陆地区域提取
 * 非水体且非 PCP 的区域
 */
function LandL8(img) {
  var a = (Water(img).eq(0)).and(PCPL8(img).eq(0));
  return a.rename('land');
}

// ============================================================
// 5. Pass Two - 云概率计算
// ============================================================

/**
 * 5.1 陆地温度概率
 * 基于亮温的百分位数计算温度概率
 */
function lTempL8_p(img) {
  var reducer1 = ee.Reducer.percentile([17.5]);
  var reducer2 = ee.Reducer.percentile([82.5]);

  var land_BT = img.select('BT').updateMask(LandL8(img));

  var Tlow = ee.Number(land_BT.reduceRegion({
    reducer: reducer1,
    maxPixels: 10e13
  }).get('BT'));

  var Thigh = ee.Number(land_BT.reduceRegion({
    reducer: reducer2,
    maxPixels: 10e13
  }).get('BT'));

  var lTemp = img.expression(
    '(Thigh + 4 - BT) / (Thigh + 4 - (Tlow - 4))', {
      Thigh: Thigh,
      BT: land_BT,
      Tlow: Tlow
    });

  return lTemp.where(lTemp.lt(0), 0);
}

/**
 * 5.2 陆地变化概率
 * 综合 NDSI、NDVI、NDBI 和白度指标
 */
function Vari_pL8(img) {
  // NDSI 修正：处理绿波段饱和情况
  var gre_satu = img.select('Green').updateMask(LandL8(img)).gte(1);
  var ndsi = NDSI(img).updateMask(LandL8(img));
  var condition1 = gre_satu.and(ndsi.lt(0)).eq(1);
  var modify_ndsi = ndsi.where(condition1, 0).abs();

  // NDVI 修正：处理红波段饱和情况
  var red_satu = img.select('Red').updateMask(LandL8(img)).gte(1);
  var ndvi = NDVI(img).updateMask(LandL8(img));
  var condition2 = red_satu.and(ndvi.gt(0)).eq(1);
  var modify_ndvi = ndvi.where(condition2, 0).abs();

  // NDBI
  var ndbi = NDBI(img).updateMask(LandL8(img));

  // 取各指标最大值
  var maxi = modify_ndsi.max(modify_ndvi).max(ndbi).max(Whitetest(img));

  // 变化概率 = 1 - max
  var lvari_P = maxi.expression('1 - max', { max: maxi });

  return lvari_P;
}

/**
 * 5.3 综合云概率
 * Cloud_p = 温度概率 × 变化概率 + 卷云概率 × 0.45 + 0.3
 */
function L8Cloud_p(img) {
  var a = lTempL8_p(img).multiply(Vari_pL8(img))
          .add(Cirrus_p(img).multiply(0.45))
          .add(0.3);
  return a.rename('Cloud_p');
}

var L8_lcp = L8.map(L8Cloud_p);
Map.addLayer(L8_lcp.first(), {}, 'L8_lcp');

// ============================================================
// 6. 最终云层生成
// ============================================================

/**
 * 潜在云层提取
 * 综合 PCP、水体掩膜与云概率阈值，形态学开运算优化
 */
function PotentialCloud_L8(img) {
  var cloud1 = PCPL8(img).eq(1).updateMask(Water(img).not()).unmask();
  var cloud2 = cloud1.eq(1).or(L8Cloud_p(img).gt(2.5)).unmask(1);

  // 形态学开运算（先腐蚀后膨胀）去除噪点
  var open = cloud2.focal_min({ radius: 1, kernelType: 'circle', units: 'pixels', iterations: 1 })
                   .focal_max({ radius: 1, kernelType: 'circle', units: 'pixels', iterations: 1 });

  return cloud2.rename('Cloud');
}

var TCL8 = L8.map(PotentialCloud_L8).first();
print('TCL8', TCL8);
Map.addLayer(TCL8, {}, 'PotentialCloud_L8');

// 导出云层结果
Export.image.toDrive({
  image: TCL8,
  description: '20220401cloud_L9_O',
  folder: 'Landsat9',
  region: roi,
  scale: 30,
  maxPixels: 1.0E13,
});

// ============================================================
// 7. 云阴影检测
// ============================================================

/**
 * 云阴影检测
 * 利用 FMask 填充最小值法检测 NIR 和 SWIR1 波段的阴影区域
 */
function CloudShadows_L8(img) {
  var nir = img.select('NIR').toInt();
  var swir1 = img.select('SWIR1').toInt();

  // FMask 填充最小值
  var finir = ee.Algorithms.FMask.fillMinima(nir, nir.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: worldWideGeom,
    maxPixels: 10e13
  }));
  var fiswir1 = ee.Algorithms.FMask.fillMinima(swir1, swir1.reduceRegion({
    reducer: ee.Reducer.max(),
    geometry: worldWideGeom,
    maxPixels: 10e13
  }));

  // 计算残差
  var res1 = finir.subtract(nir);
  var res2 = fiswir1.subtract(swir1);

  var result1 = res1.gt(0.02).mask(Water(img).not())
                     .mask(PotentialCloud_L8(img).not())
                     .unmask(0);
  var result2 = res2.gt(0.02).mask(Water(img).not())
                     .mask(PotentialCloud_L8(img).not())
                     .unmask(0);

  // 潜在阴影判定
  var pshad = img.select('NIR').lt(0.2).and(img.select('SWIR1').lt(0.15));

  var cloud1 = result1.add(result2).gt(1);
  var cloud2 = cloud1.eq(1).or(pshad.eq(1));
  var cloud3 = cloud2.eq(1).updateMask(Water(img).eq(0));
  var shadow = pshad.updateMask(PCPL8(img).not()).unmask().add(cloud3);

  return pshad.rename('cloudshadow');
}

var shadL8 = L8.map(CloudShadows_L8).first();
Map.addLayer(shadL8, {}, 'shadowL8');

// 导出云阴影结果
Export.image.toDrive({
  image: shadL8,
  description: '20220401cloudshadow_L9_O',
  folder: 'Landsat9',
  region: roi,
  scale: 30,
  maxPixels: 1.0E13,
});

// ============================================================
// 8. 积雪识别
// ============================================================

/**
 * 积雪识别
 * 基于 NDSI、NIR、Blue 和 Blue/Red 比值的多条件阈值法
 */
function SnowL8(img) {
  img = ee.Image(img);
  var s = (NDSI(img).gt(0.2))
          .and(img.select('NIR').gt(0.15))
          .and(img.select('Blue').gt(0.28))
          .and(img.select('Blue').divide(img.select('Red')).gt(0.85));
  return s.rename('snow');
}

var snowL8 = L8.map(SnowL8);
print('snowL8', snowL8);
Map.addLayer(snowL8, {}, 'snowL8');

// ============================================================
// 9. 云下积雪估计 (AllSnow)
// ============================================================

/**
 * 云下积雪估计
 * 利用积雪区域的平均高程，将云和云阴影覆盖区域中高于该高程的像元判定为积雪
 */
function AllsnowALL(img) {
  img = ee.Image(img);

  // 计算积雪区域的平均高程
  var snoweleva = img.select('elevation').updateMask(SnowL8(img));
  var a = snoweleva.where(snoweleva.lt(0), 0).selfMask().rename('elevation');
  var meanheight = ee.Number(a.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: a.geometry(),
    maxPixels: 10e13
  }).get('elevation'));

  // 云与云阴影联合掩膜
  var cloudandcloudshadow = PotentialCloud_L8(img)
                            .add(CloudShadows_L8(img)).neq(0).unmask();

  // 云覆盖区域中高于平均积雪高程的像元判定为新增积雪
  var newsnow = img.select('elevation')
                   .gt(ee.Number(meanheight))
                   .updateMask(cloudandcloudshadow)
                   .unmask().gt(0);

  // 合并直接识别的积雪与云下估计的积雪
  var allsnow = newsnow.add(SnowL8(img).unmask()).gt(0).rename('Allsnow');

  return allsnow;
}

// 取消注释以下行以运行 AllSnow 计算
// var allsnowall = L8.map(AllsnowALL);
// print('allsnowall', allsnowall);
// Map.addLayer(allsnowall, {}, 'allsnow');

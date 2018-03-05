const DEBUG_MODE = false;
const VIEW_ALL_OFF = false;
const VIEW2D_OFF = VIEW_ALL_OFF || false;
const VIEW3D_OFF = VIEW_ALL_OFF || false;
const IMAGE_OFF = VIEW_ALL_OFF || false;
const VOLREN_TIME = .3; // second
const PresetColor = {
  red: '#b82e2e',
  green: '#109618',
  blue: '#3366cc',
  white: '#ffffff',
  black: '#000000',
  gray: '#333333',
  redHex: 0xb82e2e,
  greenHex: 0x109618,
  blueHex: 0x3366cc,
  whiteHex: 0xffffff,
  blackHex: 0x000000
};

var data = {
  // tree view 
  treeData: {},
  d3treeRoot: undefined,
  treeNodeCount: undefined,
  treeMinMax: {},
  treeDepth: 0,
  
  // tf view
  tfArray: undefined,

  // 2d view
  range: undefined, // 2d 
  labelRange: undefined, // 2d

  // angle view
  isoValue: .2,
  angleData: undefined,
  coords_centroid_y: undefined,
  coords_centroid_x: undefined,
  thetaMin: Infinity,
  thetaMax: -Infinity,
  angleValueMin: Infinity,
  angleValueMax: -Infinity,

  // volren image
  enableAngle: false,
  startAngle: 0,
  endAngle: Math.PI * 2,
  enableShading: true,
  Ks: .2,
  Kd: .3,
  Ka: .04,
  lightingDirectionX: -1.0,
  lightingDirectionY: 0,
  lightingDirectionZ: 0,
  psiStart: undefined,
  psiEnd: undefined
};

var doneRendering = true;
var volrenTimer = undefined;
var globalStatus = {
  tfArray: undefined,
  tfControlPoints: undefined,
  doneRendering: true,
  volrenTimer: undefined,
  updateEnableSameTFEditor: false,
  lastUrl: undefined,
};

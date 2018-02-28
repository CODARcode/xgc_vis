const DEBUG_MODE = false;
const VIEW2D_OFF = false;
const VIEW3D_OFF = false;
const IMAGE_OFF = false;
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
  treeData: {},
  d3treeRoot: undefined,
  treeNodeCount: undefined,
  treeMinMax: {},
  treeDepth: 0,
  tfArray: undefined
};

var doneRendering = true;
var volrenTimer = undefined;
var globalStatus = {
  tfArray: undefined,
  doneRendering: true,
  volrenTimer: undefined
};

var doneRendering = true;
var data = {
  treeData: {},
  d3treeRoot: undefined,
  treeNodeCount: undefined,
  treeMinMax: {},
  treeDepth: 0,
  tfArray: undefined
};
const DEBUG_MODE = false;
const VIEW2D_OFF = false;
const VIEW3D_OFF = true;
const PresetColor = {
  red: '#b82e2e',
  green: '#109618',
  blue: '#3366cc',
  white: '#ffffff',
  black: '#000000',
  redHex: 0xb82e2e,
  greenHex: 0x109618,
  blueHex: 0x3366cc,
  whiteHex: 0xffffff,
  blackHex: 0x000000
};
var volrenTimer = undefined;
const VOLREN_TIME = .3; // second
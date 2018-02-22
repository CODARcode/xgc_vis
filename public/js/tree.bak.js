var clock = new THREE.Clock();
var scene = new THREE.Scene();
var targetList = [];
var mouse = { x: 0, y: 0 };

var camera = new THREE.PerspectiveCamera(60, window.innerWidth/window.innerHeight, 0.1, 100); 
camera.position.z = 3;

var renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0xffffff, 1);
document.body.appendChild(renderer.domElement);
renderer.domElement.id = 'webGL-canvas';

var pointLight = new THREE.PointLight(0xffffff);
pointLight.position.x = 100;
pointLight.position.y = 100;
pointLight.position.z = 100;
// pointLight.castShadow = true;
// pointLight.shadowDarkness = 0.5;
scene.add(pointLight);

var directionalLight = new THREE.DirectionalLight(0xffffff);
scene.add(directionalLight);

var controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.mouseButtons = {
  ZOOM: THREE.MOUSE.MIDDLE,
  ORBIT: THREE.MOUSE.LEFT,
  PAN: THREE.MOUSE.RIGHT,
}
var raycaster = new THREE.Raycaster();
var mousePos = new THREE.Vector2();

function render() {
  stats.begin();

  raycaster.setFromCamera(mousePos, camera);

  // scene
  var delta = clock.getDelta();
  requestAnimationFrame(render);
  // cameraControls.update(delta);
  directionalLight.position.copy(camera.position);
  renderer.render(scene, camera);

  stats.end();
}

function onResize() {
  camera.aspect = window.innerWidth/window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  // cameraControls.handleResize();
}

function getCameraMatrix() {
  return camera.matrix.toArray();
}

function setCameraMatrix(str) {
  var cameraState = JSON.parse(str);
  // ... read cameraState somehow ...
  camera.matrix.fromArray(cameraState);
  // Get back position/rotation/scale attributes
  camera.matrix.decompose(camera.position, camera.quaternion, camera.scale); 
}

function renderRadialTree() {
  var PARTICLE_SIZE = 20;
  var pointPositions = [];
  var pointColors = [];
  var linePositions = [];
  var lineColors = [];

  for (var id in treeData) {
    var i = Number(id);
    
    pointPositions.push(treeData[id].x, treeData[id].y, 0);
    pointColors.push( 0.01 + 0.1 * ( i / nodeCount ), 1.0, 0.5 );
    var parent = treeData[id].parent;
    if (parent != undefined) {
      var parentId = parent.id;
      linePositions.push(treeData[id].x, treeData[id].y, 0);
      linePositions.push(treeData[parentId].x, treeData[parentId].y, 0);
      lineColors.push(0.01 + 0.1 * ( i / nodeCount ), 1.0, 0.5);
      lineColors.push(0.01 + 0.1 * ( i / nodeCount ), 1.0, 0.5);
    }
  }
  console.log(pointPositions, pointColors);
  var pointGeometry = new THREE.BufferGeometry();
  pointGeometry.addAttribute('position', new THREE.Float32BufferAttribute(pointPositions, 3));
  pointGeometry.addAttribute('color', new THREE.Float32BufferAttribute(pointColors, 3));
  pointGeometry.computeBoundingSphere();
  var pointSprite = new THREE.TextureLoader().load("img/disc.png");
  var pointMaterial = new THREE.PointsMaterial({
    size: .01, 
    vertexColors: THREE.VertexColors, 
    map: pointSprite,
    alphaTest: 0.5,
    transparent: true
  });
  var points = new THREE.Points( pointGeometry, pointMaterial );
  console.log(points);
  scene.add( points );

  var lineGeometry = new THREE.BufferGeometry();
  var lineGeometry = new THREE.BufferGeometry();
  lineGeometry.addAttribute( 'position', new THREE.Float32BufferAttribute( linePositions, 3 ) );
  lineGeometry.addAttribute( 'color', new THREE.Float32BufferAttribute( lineColors, 3 ) );
  var lineMaterial = new THREE.LineBasicMaterial({color: 0x0000ff});
  var lines = new THREE.LineSegments(lineGeometry, lineMaterial);
  console.log(lines);
  scene.add( lines );
  
}

render();

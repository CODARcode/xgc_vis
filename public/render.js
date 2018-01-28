var stats = new Stats();
stats.showPanel(0);
stats.dom.id="stats";
document.body.appendChild(stats.dom);
$("#stats").css({visibility: "hidden"});

var clock = new THREE.Clock();
var scene = new THREE.Scene();

var camera = new THREE.PerspectiveCamera(30, window.innerWidth/window.innerHeight, 0.1, 100); 
camera.position.z = 3;

var renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0xffffff, 1);
document.body.appendChild(renderer.domElement);

var pointLight = new THREE.PointLight(0xffffff);
pointLight.position.x = 100;
pointLight.position.y = 100;
pointLight.position.z = 100;
// pointLight.castShadow = true;
// pointLight.shadowDarkness = 0.5;
scene.add(pointLight);

var directionalLight = new THREE.DirectionalLight(0xffffff);
scene.add(directionalLight);

cameraControls = new THREE.TrackballControls(camera, renderer.domElement);
cameraControls.target.set(0, 0, 0);
cameraControls.zoomSpeed = 0.04;
cameraControls.panSpeed = 0.04;
// cameraControls.addEventListener("change", render); // not working.. sigh

var raycaster = new THREE.Raycaster();
var mousePos = new THREE.Vector2();

window.addEventListener("mousedown", onMouseDown, false );
window.addEventListener("mousemove", onMouseMove, false);
window.addEventListener("resize", onResize, false);

function render() {
  stats.begin();

  raycaster.setFromCamera(mousePos, camera);
  var intersects = raycaster.intersectObjects(scene.children);
  // for (i=0; i<intersects.length; i++)
  //   intersects[i].object.material.color.set(0xff0000);

  // scene
  var delta = clock.getDelta();
  requestAnimationFrame(render);
  cameraControls.update(delta);
  directionalLight.position.copy(camera.position);
  renderer.render(scene, camera);

  stats.end();
}

function onMouseDown(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onMouseMove(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onResize() {
  camera.aspect = window.innerWidth/window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  cameraControls.handleResize();
}

function updateMesh(mesh) {
  console.log("updating mesh...");

  var geom = new THREE.Geometry();
  for (var i = 0; i < mesh.nNodes; i ++) {
    geom.vertices.push(new THREE.Vector3(mesh.coords[i*2] - 1.7, mesh.coords[i*2+1], 0));
  }
  for (var i = 0; i < mesh.nTriangles; i ++) {
    var face = new THREE.Face3(mesh.conn[i*3], mesh.conn[i*3+1], mesh.conn[i*3+2]);
    face.color.setRGB( Math.random(), Math.random(), Math.random() );
    geom.faces.push(face);
  }

  var material = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: false, vertexColors: THREE.FaceColors } );
  var obj = new THREE.Mesh( geom, material );
  scene.add(obj);

  requestData();
}

initializeControlPanel();
render();
